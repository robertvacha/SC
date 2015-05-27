#include "updater.h"

#include <iomanip>

void Updater::openFilesClusterStatistics(FILE** cl_stat, FILE** cl, FILE** cl_list, FILE** ef, FILE** statf) {

    // Opening files for cluster statistics
    *cl_stat = *cl = *cl_list = *ef = *statf = NULL;
    if(sim->write_cluster){
        // Empty file
        *cl_stat = fopen(files->clusterstatfile, "w");
        fclose(*cl_stat);
        *cl_stat = fopen(files->clusterstatfile, "a");
        // Empty file
        *cl = fopen(files->clusterfile, "w");
        fclose(*cl);
        *cl = fopen(files->clusterfile, "a");
    }
    // write energy
    if (report <= nsweeps){
        // Empty file
        *ef = fopen(files->energyfile, "w");
        fclose(*ef);
        *ef = fopen(files->energyfile, "a");
        fprintf (*ef, "# sweep    energy\n");
        *statf = fopen(files->statfile, "w");
        fclose(*statf);
        *statf = fopen(files->statfile, "a");
        fprintf (*statf, "# sweep    volume\n");
    }
}



void Updater::initValues(long& next_adjust, long& next_calc, long& next_dump, long& next_frame) {

    // Initialize some values at begining
    conf->partVecInit();
    next_adjust = adjust;
    next_calc = paramfrq;
    next_dump = report;
    next_frame = sim->movie;
    //nem = vol = shapex = shapey = shapez = nullstat;
    //for (i=0; i<MAXF; i++) smec[i] = nullstat;

    sim->wl.wl_meshsize = 0;
    sim->wl.radiushole = NULL;
    sim->wl.radiusholeold = NULL;
    sim->wl.radiusholemax = 0;
    sim->wl.partincontactold = 0;
    sim->wl.partincontact = 0;
    sim->wl.wlmdim = 0;
    sim->wl.wlmdim = 0;
    sim->wl.length[0]=0;
    sim->wl.length[1]=0;
    sim->wl.currorder[0]=0;
    sim->wl.currorder[1]=0;
    sim->wl.neworder[0]=0;
    sim->wl.neworder[1]=0;
    sim->wl.weights = NULL;
    sim->wl.hist = NULL;

    conf->massCenter();
}




void Updater::simulate(long nsweeps, long adjust, long paramfrq, long report) {
    long i,j;

    this->nsweeps = nsweeps;
    this->adjust = adjust;
    this->paramfrq = paramfrq;
    this->report = report;

    long next_adjust;  // Next sweep number for step size adjustment
    long next_calc;    // Next sweep number for order parameter calculation
    long next_dump;    // Next sweep number for reporting statistics
    long next_frame;   // Next sweep number for dumping a movie fram
    long step;         // Step number within a given sweep
    long sweep;        // Current sweep number

    //struct stat nem;   // Nematic order parameter
    //struct stat vol;   // Volume statistic
    //struct stat shapex, shapey, shapez;   // geo.box shape statistics
    //struct stat smec[MAXF];       // Smectic order parameters (Fourier coeeficients)

    FILE *mf=NULL;                                  // Handle for movie file
    FILE *cl_stat=NULL, *cl=NULL, *cl_list=NULL;    // Handle for cluster statistics
    FILE *ef=NULL, *statf=NULL;                     // Handle for energy file and statistical file

    double edriftstart;     // Energy drift calculation - start
    double edriftchanges;   // Energy drift calculation - accumulate all changes through moves
    double edriftend;       // Energy drift calculation - end
    double pvdriftstart;    // PV drift calculation - start
    double pvdriftend;      // PV drift calculation - end
    double volume;          // volume of geo.box
    double moveprobab;      // random number selecting the move

    openFilesClusterStatistics(&cl_stat, &cl, &cl_list, &ef, &statf);
    //=== Initialise counters etc. ===//
    sim->shprob = sim->shave/(double)conf->pvec.size();

    initValues(next_adjust, next_calc, next_dump, next_frame);
    if (sim->movie > 0) {
        mf = fopen(files->moviefile, "a");
    } else {
        mf = NULL;
    }

    sim->wl.init(files->wlinfile);
    //do moves - START OF REAL MC
    if(sim->pairlist_update){
        genPairList(); // Does that solve the problem?
    }

    //do energy drift check - start calculation
    volume = conf->geo.volume();
    edriftstart = calcEnergy(0, 0, 0);
    pvdriftstart = sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper;
    //printf("starting energy: %.15f \n",calc_energy(0, intfce, 0, topo, conf, sim,0));
    //printf("press: %.15f\n",sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper);
    edriftchanges = 0.0;

    /********************************************************/
    /*                 Simulation Loop                      */
    /********************************************************/
    for (sweep=1; sweep <= nsweeps; sweep++) {
//        if(nsweeps>=10 && sweep%(nsweeps/10) == 0) {
//            volume = conf->geo.volume();
//            edriftend = calcEnergy.allToAll();
//            pvdriftend =  sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper;
//            cout << "sweep: " << sweep << " particles: " << conf->pvec.size()
//                 << " drift:" << edriftend - edriftstart - edriftchanges +pvdriftend -pvdriftstart << endl;
//        }

        // Try replica exchange
        if((sim->nrepchange) && (sweep % sim->nrepchange == 0)){
            edriftchanges += move.replicaExchangeMove(sweep);

            if(sim->pairlist_update)
                genPairList();
        }
        /*____________GrandCanonical Move____________*/
        if(sim->nGrandCanon != 0 && sweep%sim->nGrandCanon == 0) {
            edriftchanges += move.muVTMove();

            if(sim->pairlist_update)
                genPairList();
        }
        /*____________Cluster Move____________*/
        if(sim->nClustMove != 0 && sweep%sim->nClustMove == 0) {
            edriftchanges += move.clusterMove();

            if(sim->pairlist_update)
                genPairList();
        }
        if( (sim->pairlist_update) && // pair_list allowed
                (
                    (sweep % sim->pairlist_update == 0) && // on scheduled sweep
                    !(sim->nGrandCanon != 0 && sweep%sim->nGrandCanon == 0) && // not on grandCanon sweep
                    !(sim->nClustMove != 0 && sweep%sim->nClustMove == 0) &&
                    !((sim->nrepchange) && (sweep % sim->nrepchange == 0))  // not on replica exchange sweep
                )
                ) {
            genPairList();
        }
        //normal moves
        for (step=1; step <= (long)conf->pvec.size(); step++) {
            moveprobab = ran2();

            if ( moveprobab < sim->shprob) {
                edriftchanges += move.pressureMove();
                continue;
            }
            if (moveprobab < sim->shprob + sim->chainprob) {
                edriftchanges += move.chainMove();
                continue;
            }
            if (moveprobab < sim->shprob + sim->chainprob + sim->switchprob){
                //=== This is an attempt to switch a type ===
                edriftchanges += move.switchTypeMove();
            } else {
                // single particle moves
                edriftchanges += move.particleMove();
            }
            //TEST OVERLAPS
            if(conf->checkall(topo.ia_params)){
                cout<<"OVERLAP DETECTED"<<endl;
            }
        } // End of step loop for this sweep

        //=== Start of end-of-sweep housekeeping ===
        // Adjustment of maximum step sizes during equilibration
        if (sweep == next_adjust) {
            for (i = 0; i < MAXT ;i++) {
                if ((sim->trans[i].acc > 0)||(sim->trans[i].rej >0))
                    optimizeStep (sim->trans + i, 1.5, 0.0);
                if ((sim->rot[i].acc > 0)||(sim->rot[i].rej >0))
                    optimizeRot (sim->rot + i, 5.0, 0.01);
            }
            for (i = 0; i < MAXMT; i++) {
                if ((sim->chainm[i].acc > 0)||(sim->chainm[i].rej > 0))
                    optimizeStep (sim->chainm + i, 1.5, 0.0);
                if ((sim->chainr[i].acc > 0)||(sim->chainr[i].rej > 0))
                    optimizeRot (sim->chainr + i, 5.0, 0.01);
            }
            optimizeStep (&(sim->edge), 1.0, 0.0);
            next_adjust += adjust;
        }

        if ( (sim->wl.wlm[0] > 0) && (sim->wl.alpha > WL_ZERO) && !(sweep % 1000) ) {
            // recalculate system CM to be sure there is no accumulation of errors by +- rejection moves
            /* BUG - not used any longer: caused problems with PBC normal moves systemCM movement
              can be calculated from CM movements of individual particles
              present center of mass calculation use pbc and thus particles that moved across the geo.box
              is in this calculation used in pripary geo.box but in other moves in in in the particles position
             if ( (sim->wlm[0] == 1) || (sim->wlm[1] == 1) )
                  masscenter(topo.npart,topo.ia_params, conf);
            */
            sim->wl.min = sim->wl.hist[0];
            sim->wl.max = sim->wl.hist[0];
            for (i=0;i < sim->wl.length[0];i++) {
                j=0;
                if ( sim->wl.hist[i+j*sim->wl.length[0]] > sim->wl.max ) sim->wl.max = sim->wl.hist[i+j*sim->wl.length[0]];
                if ( sim->wl.hist[i+j*sim->wl.length[0]] < sim->wl.min ) sim->wl.min = sim->wl.hist[i+j*sim->wl.length[0]];
                for (j=1;j < sim->wl.length[1];j++) {
                    if ( sim->wl.hist[i+j*sim->wl.length[0]] > sim->wl.max ) sim->wl.max = sim->wl.hist[i+j*sim->wl.length[0]];
                    if ( sim->wl.hist[i+j*sim->wl.length[0]] < sim->wl.min ) sim->wl.min = sim->wl.hist[i+j*sim->wl.length[0]];
                }
            }
            if ( sim->wl.min > WL_MINHIST ) {
                if ( sim->temper * log(sim->wl.max/sim->wl.min) < WL_GERR ) {
                    /*DEBUG
                      for (i=1;i<wl.length;i++) {
                      printf (" %15.8e %15ld %15.8f\n",sim->wl.weights[i],sim->wl.hist[i],pvec[0].pos.z);
                      fflush(stdout);
                      }
                     */
                    if ( sim->wl.alpha < WL_ALPHATOL) break;
                    sim->wl.alpha/=2;
                    printf("%f \n", sim->wl.alpha);
                    fflush (stdout);
                    sim->wl.wmin = sim->wl.weights[0];

                    for (i=0;i < sim->wl.length[0];i++) {
                        j=0;
                        sim->wl.hist[i+j*sim->wl.length[0]] = 0;
                        sim->wl.weights[i+j*sim->wl.length[0]] -= sim->wl.wmin;
                        for (j=1;j < sim->wl.length[1];j++) {
                            sim->wl.hist[i+j*sim->wl.length[0]] = 0;
                            sim->wl.weights[i+j*sim->wl.length[0]] -= sim->wl.wmin;
                        }
                    }

                }
            }
        }

        if (!(sweep % 100000)) {
            //reinitialize patch vectors to avoid cummulation of errors
            conf->partVecInit();
        }

        /// Sampling of statistics
        if (sweep == next_calc)
        {
            /*s2 = nematic(npart, particle);
              accumulate (&nem, s2);
              for (i=0; i<terms; i++) {
              ci = smectic(npart, particle, i+1);
              accumulate (&smec[i], ci);
              }
              accumulate (&shapex, (*geo.box).x);
              accumulate (&shapey, (*geo.box).y);
              accumulate (&shapez, (*geo.box).z);
              volume = (*geo.box).x * (*geo.box).y * (*geo.box).z;
              accumulate (&vol, volume);
              next_calc += paramfrq;
             */
        }
        /// Writing of statistics
        if (sweep == next_dump) {
            /*printf ("Statistics after %ld sweeps:\n", sweep);
              printf ("   Mean and RMS fluctuation of S2:  %13.8f %13.8f\n",
              nem.mean, nem.rms);
              for (i=0; i<terms; i++) {
              printf ("   Mean & fluc. Fourier coeff. %3ld: %13.8f %13.8f\n",
              i+1, smec[i].mean, smec[i].rms);
              }
              printf ("   Mean & fluc geo.box dimensions:  x   %13.8f %13.8f\n",
              shapex.mean, shapex.rms);
              printf ("                                y   %13.8f %13.8f\n",
              shapey.mean, shapey.rms);
              printf ("                                z   %13.8f %13.8f\n",
              shapez.mean, shapez.rms);
              printf ("   Mean & fluctuation volume:      %13.8f %13.8f\n",
              vol.mean, vol.rms);
              printf ("   Mean & fluc. volume over volume of particles:    %13.8f %13.8f\n",
              vol.mean/pvolume, vol.rms/pvolume);
              printf ("\n");
              fflush (stdout);
             */

            fprintf (statf, " %ld; %.10f\n", sweep, conf->geo.box.x * conf->geo.box.y * conf->geo.box.z);
            fprintf (ef, " %ld; %.10f  %f \n", sweep, calcEnergy(0, 0, 0), alignmentOrder());
            if (sim->wl.wlm[0] > 0) {
                sim->wl.write(files->wloutfile);
            }
            //print mesh distribution
            //mesh_findholesdistrib(&sim->wl.mesh);
            next_dump += report;
        }

        // Writing of movie frame
        if (sweep == next_frame) {
            //fprintf (mf, "> box %.10f %.10f %.10f ; num_part %ld ; sweep %ld <\n", conf->geo.box.x, conf->geo.box.y, conf->geo.box.z, (long)conf->pvec.size(), sweep);
            fprintf (mf, "%ld\nsweep %ld; box %.10f %.10f %.10f\n",(long)conf->pvec.size(), sweep, conf->geo.box.x, conf->geo.box.y, conf->geo.box.z);
            conf->draw(mf);
            fflush (mf);
            next_frame += sim->movie;
        }

        // Writing out cluster statistics
        if(sim->write_cluster && (sweep % sim->write_cluster == 0)){
            writeCluster(cl_stat, cl, cl_list, false, sweep);
        }
        //=== End of housekeeping ===

    }  // End of sweeps loop

    //do energy drift check - at the end calculation
    volume = conf->geo.volume();
    edriftend = calcEnergy.allToAll();
    pvdriftend =  sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper;
    printf("Energy drift: %.5e \n",edriftend - edriftstart - edriftchanges +pvdriftend -pvdriftstart);
//    printf("EdriftChanges: %.5e\n", edriftchanges);
//    printf("PVChanges: %.5e\n", pvdriftend -pvdriftstart);
    printf("Starting energy: %.8f \n",edriftstart);
    printf("Ending energy: %.8f \n",edriftend);
    printf("Starting energy+pv: %.8f \n",edriftstart+pvdriftstart);
    printf("System:\n");

    for(i=0; i < conf->pvec.molTypeCount; i++)
        printf("%s %d\n", topo.moleculeParam[i].name, conf->pvec.molCountOfType(i));

    if(sim->nGrandCanon != 0) {
        cout << "Acceptance:\n";
        cout << "  Type   insAcc insRej delAcc delRej <num of particles>\n";
        cout << std::setprecision(3) <<  std::fixed << std::left;
        for(int i=0; i<conf->pvec.molTypeCount; i++) {
            if(topo.moleculeParam[i].activity != -1.0) {
                cout << "  " << std::setw(6) <<topo.moleculeParam[i].name << " "
                     << std::setw(6)
                     << (double)topo.moleculeParam[i].insAcc/(topo.moleculeParam[i].insAcc+topo.moleculeParam[i].insRej) << " "
                     << std::setw(6)
                     << (double)topo.moleculeParam[i].insRej/(topo.moleculeParam[i].insAcc+topo.moleculeParam[i].insRej) << " "
                     << std::setw(6)
                     << (double)topo.moleculeParam[i].delAcc/(topo.moleculeParam[i].delAcc+topo.moleculeParam[i].delRej) << " "
                     << std::setw(6)
                     << (double)topo.moleculeParam[i].delRej/(topo.moleculeParam[i].delAcc+topo.moleculeParam[i].delRej) << " "
                     << std::setw(6)
                     << (double)topo.moleculeParam[i].muVtAverageParticles / topo.moleculeParam[i].muVtSteps << endl;
            }
        }
        cout << std::setprecision(6);
    }
    fflush(stdout);

    sim->wl.endWangLandau(files->wloutfile);

    //end movie
    if (sim->movie > 0)
        fclose (mf);
    //end cluster
    if(sim->write_cluster){
        fclose(cl_stat);
        fclose(cl);
    }
    if (report < nsweeps) {
        fclose(ef);
        fclose(statf);
    }
}


void Updater::optimizeStep(Disp *x, double hi, double lo) {
    double newrmsd;

    newrmsd = (*x).mx * RATIO(*x);
    if ((*x).oldrmsd > 0) {
        if ( newrmsd  < (*x).oldrmsd ) {
            if ( (*x).oldmx > 1 ) {
                (*x).mx /= 1.05;
                (*x).oldmx = 0.95;
            } else {
                (*x).mx *= 1.05;
                (*x).oldmx = 1.05;
            }
        } else {
            if ( (*x).oldmx > 1 ) {
                (*x).mx *= 1.05;
                (*x).oldmx = 1.05;
            } else {
                (*x).mx /= 1.05;
                (*x).oldmx = 0.95;
            }
        }
    }
    if (newrmsd > 0 ) (*x).oldrmsd = newrmsd;
    else {
        (*x).oldrmsd = 0.0;
        (*x).mx /= 1.05;
        (*x).oldmx = 0.95;
    }

    if ( (*x).mx > hi ) (*x).mx = hi;
    if ( (*x).mx < lo ) (*x).mx = lo;

    (*x).acc = (*x).rej = 0;
}

void Updater::optimizeRot(Disp *x, double hi, double lo) {
    double newrmsd;

    newrmsd = (*x).mx * RATIO((*x)) ;
    if ((*x).oldrmsd > 0) {
        if ( newrmsd  > (*x).oldrmsd ) {
            if ( (*x).oldmx > 1) {
                (*x).mx *= 0.99;
                (*x).oldmx *= 0.99;
            } else {
                (*x).mx *= 1.01;
                (*x).oldmx *= 1.01;
            }
        } else {
            if ( (*x).oldmx > 1) {
                (*x).mx *= 1.01;
                (*x).oldmx *= 1.01;
            } else {
                (*x).mx *= 0.99;
                (*x).oldmx *= 0.99;
            }
        }
    }
    if (newrmsd > 0 ) (*x).oldrmsd = newrmsd;
    else {
        (*x).oldrmsd = 0.0;
        (*x).mx *= 1.01;
        (*x).oldmx = 1.01;
    }

    if ( (*x).mx > hi ) (*x).mx = hi;
    if ( (*x).mx < lo ) (*x).mx = lo;

    (*x).acc = (*x).rej = 0;
}

double Updater::alignmentOrder() {
    double sumdot=0;
    long i,j;
    Vector r_cm;

    for (i = 0; i < (long)(long)conf->pvec.size() - 1; i++) {
        for (j = i + 1; j < (long)(long)conf->pvec.size(); j++) {
            r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[j].pos);
            if ( DOT(r_cm,r_cm) < 1.5*1.5 ) {
                sumdot+= DOT(conf->pvec[i].dir,conf->pvec[j].dir);
            }
        }
    }

    return sumdot;
}

void Updater::genPairList() {
    genSimplePairList();
}

void Updater::genSimplePairList() {
    Vector r_cm;
    double r_cm2;
    double max_dist;
    // Set the pairlist to zero
    //DEBUG_INIT("Gen Pairlist")
    for(unsigned int i = 0; i < conf->neighborList.size(); i++){
        //DEBUG_INIT("%ld", i);
        conf->neighborList[i].neighborCount = 0;
    }
    if(conf->pvec.size() <= 0)
        return ;
    for(unsigned int i = 0; i < conf->pvec.size()-1; i++){
        for(unsigned int j = i + 1; j < conf->pvec.size(); j++){
            assert(conf->pvec.size() == conf->neighborList.size());
            assert(&conf->pvec[i] != NULL);
            assert(&conf->pvec[j] != NULL);

            r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[j].pos);

            /*r_cm.x = conf->pvec[i].pos.x - conf->pvec[j].pos.x;
            r_cm.y = conf->pvec[i].pos.y - conf->pvec[j].pos.y;
            r_cm.z = conf->pvec[i].pos.z - conf->pvec[j].pos.z;
            if ( r_cm.x < 0  )
                r_cm.x = conf->geo.box.x * (r_cm.x - (double)( (long)(r_cm.x-0.5) ) );
            else
                r_cm.x = conf->geo.box.x * (r_cm.x - (double)( (long)(r_cm.x+0.5) ) );
            if ( r_cm.y < 0  )
                r_cm.y = conf->geo.box.y * (r_cm.y - (double)( (long)(r_cm.y-0.5) ) );
            else
                r_cm.y = conf->geo.box.y * (r_cm.y - (double)( (long)(r_cm.y+0.5) ) );
            if ( r_cm.z < 0  )
                r_cm.z = conf->geo.box.z * (r_cm.z - (double)( (long)(r_cm.z-0.5) ) );
            else
                r_cm.z = conf->geo.box.z * (r_cm.z - (double)( (long)(r_cm.z+0.5) ) );*/

            r_cm2 = DOT(r_cm,r_cm);
            max_dist = AVER(sim->trans[conf->pvec[i].type].mx, \
                    sim->trans[conf->pvec[j].type].mx);
            max_dist *= (1 + sim->pairlist_update) * 2;
            max_dist += topo.maxcut;
            max_dist *= max_dist; /* squared */

            if (r_cm2 <= max_dist){
                conf->neighborList[i].neighborID[conf->neighborList[i].neighborCount] = j;
                conf->neighborList[j].neighborID[conf->neighborList[j].neighborCount] = i;

                conf->neighborList[i].neighborCount++;
                conf->neighborList[j].neighborCount++;

                /*sim->pairlist[i].pairs[sim->pairlist[i].num_pairs++] = j; // DEL AFTER
                sim->pairlist[j].pairs[sim->pairlist[j].num_pairs++] = i;*/
            }
        }
    }
    ////Check for too many pairs
    //for(i = 0; i < (long)conf->pvec.size(); i++){
    //    //if (sim->pairlist.list[i].num_pairs >= (long)conf->pvec.size())
    //    if (sim->pairlist[i].num_pairs >= (long)conf->pvec.size()){
    //        fprintf(stderr, "ERROR: Too many pairs for particle %ld!!!\n", i);
    //        exit(1);
    //    }
    //}
}


int Updater::writeCluster(FILE *cl_stat, FILE *cl, FILE *cl_list, bool decor, long sweep) {
    genClusterList();
    sortClusterList();
    calcClusterEnergies();
    if(cl_stat){
        if(decor == false){
            // if no decor, this means usually into a file. Hence print info
            // about number of line per frame
            fprintf(cl_stat, "Sweep: %ld | Maximal size: %ld\n",
                    sweep, sim->max_clust);
        }
        printStat::printClusterStat(cl_stat, decor, sim);
        /*
           print_clstat_oneline(cl_stat, sweep, sim);
         */
    }
    if(cl){
        if(decor == false){
            fprintf(cl, "Sweep: %ld | Number of clusters: %ld\n",
                    sweep, sim->num_cluster);
        }
        printStat::printClusters(cl, decor, sim);
    }
    if(cl_list){
        if(decor == false){
            fprintf(cl_list, "Sweep: %ld | Number of particles: %ld\n",
                    sweep, (long)(long)conf->pvec.size());
        }
        printStat::printClusterList(cl, decor, sim, conf);
    }
    return 0;
}


int Updater::genClusterList() {
    bool change = true; /* does it still change? */
    //long neighbour;
    long i, j, fst, snd, tmp, minnumber, maxnumber;

    // Set clusterindex to the corresponding index
    for( i = 0; i < (long)conf->pvec.size(); i++){
        sim->clusterlist[i] = i;
    }

    // Start determining the cluster
    while(change){
        change = false;
        for(i = 0; i < (long)conf->pvec.size(); i++){
            /*If nore pairlist go over all pairs*/
            maxnumber = (long)conf->pvec.size();
            minnumber = i ;
            if (sim->pairlist_update) {
                maxnumber = conf->neighborList[i].neighborCount;
                //maxnumber = sim->pairlist[i].num_pairs; // del after
                minnumber=0;
            }
            /* Go over pairs to see if they are in the cluster */
            for(j = minnumber; j < maxnumber; j++){
                fst = i;
                snd = j;
                if (sim->pairlist_update) {
                    snd = conf->neighborList[i].neighborID[j];
                    //snd = sim->pairlist[i].pairs[j];
                }
                /*do cluster analysis only for spherocylinders*/
                if ( (topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].geotype[0] < SP) && \
                    (topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].geotype[1] < SP) ) {

                    /* if they are close to each other */
                    if(sameCluster(fst, snd)){
                        if(fst > snd){
                            tmp = snd;
                            snd = fst;
                            fst = tmp;
                        }

                        if(sim->clusterlist[fst] < sim->clusterlist[snd]){
                            sim->clusterlist[snd] = sim->clusterlist[fst];
                            change = true;
                            break;
                            /* => will eventually start the i loop from new */
                        }
                        if(sim->clusterlist[snd] < sim->clusterlist[fst]){
                            sim->clusterlist[fst] = sim->clusterlist[snd];
                            change = true;
                            break;
                            /* => will eventually start the i loop from new */
                        }
                    }
                }
            }
            if(change){
                break;
            }
        }
    }

    return 0;
}

int Updater::sortClusterList() {
    long cluster_indices[(long)conf->pvec.size()];   /* holds the different cluster indices.
                        (currently too much memory) */
    long num_cluster = 0;                /* number of clusters, temporary needed */

    /* how many clusters are there? */
    long max_index = -1;
    for(int i = 0; i < (long)conf->pvec.size(); i++){
        if(max_index < sim->clusterlist[i]){
            max_index = sim->clusterlist[i];
            cluster_indices[num_cluster++] = max_index;
        }
    }

    /* free the memory from the old clusters */
    if(sim->clusters){
        for(int i = 0; i < sim->num_cluster; i++){
            if(sim->clusters[i].particles){
                free(sim->clusters[i].particles);
            }
        }
        free(sim->clusters);
    }

    /* Allocate memory for the clusters */
    sim->clusters = (Cluster*) malloc(sizeof(Cluster) * num_cluster);

    if (!sim->clusters){
        fprintf(stderr, "Couldn't allocate any memory!\n");
        exit(1);
    }

    for(int i = 0; i < num_cluster; i++){
        /* allocate maximal space for all the clusters */
        sim->clusters[i].particles = (long int*) malloc(sizeof(long) * (long)conf->pvec.size());

        if (!sim->clusters[i].particles){
            fprintf(stderr, "Couldn't allocate any memory!\n");
            exit(1);
        }

        sim->clusters[i].npart = 0;
    }

    /* fill in the particles belonging to one cluster */
    for(int i = 0; i < num_cluster; i++){
        for(int j = 0; j < (long)conf->pvec.size(); j++){
            if(sim->clusterlist[j] == cluster_indices[i]){
                sim->clusters[i].particles[sim->clusters[i].npart++] = j;
            }
        }
    }
    sim->num_cluster = num_cluster;

    /* Find the biggest size */
    sim->max_clust = 0;
    for(int i = 0; i < num_cluster; i++){
        if(sim->clusters[i].npart > sim->max_clust){
            sim->max_clust = sim->clusters[i].npart;
        }
    }
    /* Set the statistics to zero */
    sim->clusterstat = (long int*) malloc(sizeof(long) * sim->max_clust);

    if (!sim->clusterstat){
        fprintf(stderr, "Couldn't allocate any memory!\n");
        exit(1);
    }

    for(int i = 0; i < sim->max_clust; i++){
        sim->clusterstat[i] = 0;
    }
    /* Do the statistics */
    for(int i = 0; i < num_cluster; i++){
        sim->clusterstat[sim->clusters[i].npart - 1]++;
    }

    return 0;
}

int Updater::calcClusterEnergies() {
    for(int i = 0; i < sim->num_cluster; i++) {
        sim->clustersenergy[i]=0.0;
        for(int j = 0; j < sim->clusters[i].npart; j++) {
            for(int k = j+1; k < sim->clusters[i].npart; k++) {
                sim->clustersenergy[i]+= calcEnergy.p2p(sim->clusters[i].particles[j], // particle 1
                                                        sim->clusters[i].particles[k]);        // particle 2
            }
        }
    }
    return 0;
}



int Updater::sameCluster(long fst, long snd) {

    ConList conFst = conf->pvec.getConlist(fst);
    ConList conSnd = conf->pvec.getConlist(snd);

    /*if two particles are bonded they belong to the same cluster*/
    if ( ((topo.moleculeParam[conf->pvec[fst].molType]).bond1c >= 0) ||
        ((topo.moleculeParam[conf->pvec[fst].molType]).bonddc >= 0) ){
        if ( (&conf->pvec[snd] == conFst.conlist[1]) || (&conf->pvec[snd] == conFst.conlist[0]) ) {
          return true;
        }
    }
    if ( ((topo.moleculeParam[conf->pvec[snd].molType]).bond1c >= 0) ||
        ((topo.moleculeParam[conf->pvec[snd].molType]).bonddc >= 0) ){
        if ( (&conf->pvec[fst] == conSnd.conlist[1]) || (&conf->pvec[fst] == conSnd.conlist[0]) ) {
          return false;
        }
    }

    /*cluster is made of particles closer tna some distance*/
/*	struct vector2 image(struct vector2 r1, struct vector2 r2, struct vector2 geo.box);
    struct vector2 r_cm = image(conf->pvec[fst].pos,
            conf->pvec[snd].pos,
            conf->geo.box);
    double dist2 = DOT(r_cm, r_cm);
    * TODO: Make it much more efficient => define cluster_dist!!! *
    if(dist2 > topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].sigma * topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].sigma*4.0){
        return false;
    }
    else {
        return true;
    }*/

    /*cluster is made of attractively interacting particles*/
    /*double paire(long, long, double (* intfce[MAXT][MAXT])(struct interacts *),
            struct topo * topo, struct conf * conf); Redeclaration*/

    if(calcEnergy.p2p(fst, snd) > -0.10 ){
        return false;
    }
    else {
        return true;
    }
}

