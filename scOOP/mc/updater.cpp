#include "updater.h"

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
    if (report < nsweeps){
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
    conf->partVecInit(topo);
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

    conf->massCenter(topo);
}

void Updater::initWangLandau() {
    if ( sim->wlm[0] >0 ) {
        if (sim->wl.initCalc(files->wlinfile) != 0)
            return;
        sim->wl.wlmdim = 1 ;
        if ( sim->wlm[1] > 0 )
            sim->wl.wlmdim = 2 ;
        for (long wli=0; wli < sim->wl.wlmdim; wli++) {
            switch (sim->wlm[wli]) {
                case 1:
                    conf->massCenter(topo);
                    sim->wl.currorder[wli] = sim->wl.zOrder(wli);
                    break;
                case 2:
                    sim->wl.wl_meshsize = (topo->ia_params[sim->wl.wlmtype][sim->wl.wlmtype].sigma) / 3.0; // TODO
                    sim->wl.mesh.data = NULL;
                    sim->wl.mesh.tmp = NULL;
                    sim->wl.origmesh.data = NULL;
                    sim->wl.origmesh.tmp = NULL;
                    sim->wl.currorder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize, (long)conf->particleStore.size(), sim->wl.wlmtype) - sim->wl.minorder[wli]);
                    break;
                case 3:
                    sim->wl.currorder[wli] = (long) floor( (conf->particleStore[0].dir.z - sim->wl.minorder[wli])/ sim->wl.dorder[wli] );
                    break;
                case 4:
                    sim->wl.currorder[wli] = sim->wl.twoPartDist(wli);
                    break;
                case 5:
                    conf->massCenter(topo);
                    sim->wl.radiusholemax = 0;
                    sim->wl.radiushole = NULL;
                    sim->wl.radiusholeold = NULL;
                    sim->wl.currorder[wli] = move.radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    sim->wl.radiusholemax = 0;
                    sim->wl.radiushole = NULL;
                    sim->wl.radiusholeold = NULL;
                    sim->wl.currorder[wli] = move.radiusholeAll(wli,&(conf->particleStore[0].pos));
                    break;
                case 7:
                    sim->wl.currorder[wli] = move.contParticlesAll(wli);
                    break;
                default:
                    sim->wl.currorder[wli] = 0;
                    break;
            }
            if ( (sim->wl.currorder[wli] >= sim->wl.length[wli] ) || (sim->wl.currorder[wli] < 0) ) {
                printf("Error: starting Wang-Landau method with order parameter %f out of range(%f - %f)\n\n", sim->wl.dorder[wli]*sim->wl.currorder[wli] + \
                   sim->wl.minorder[wli], sim->wl.minorder[wli], sim->wl.minorder[wli]+sim->wl.dorder[wli]*sim->wl.length[wli]  );
                sim->wl.end();
                return;
            }
        }
        if (sim->wl.alpha < WL_ALPHATOL/100) sim->wl.alpha = WL_ZERO;
        fflush (stdout);
    }
}


void Updater::endWangLandau() {
    long i=0,j=0;
    if (sim->wlm[0] > 0) {
        sim->wl.min = sim->wl.hist[0];
        for (i=0;i < sim->wl.length[0];i++) {
            j=0;
            if ( sim->wl.hist[i+j*sim->wl.length[0]] < sim->wl.min ) sim->wl.min = sim->wl.hist[i+j*sim->wl.length[0]];
            for (j=1;j < sim->wl.length[1];j++) {
                if ( sim->wl.hist[i+j*sim->wl.length[0]] < sim->wl.min ) sim->wl.min = sim->wl.hist[i+j*sim->wl.length[0]];
            }
        }
        sim->wl.wmin = sim->wl.weights[0];
        for (i=0;i < sim->wl.length[0];i++) {
            j=0;
            sim->wl.weights[i+j*sim->wl.length[0]] -= sim->wl.wmin;
            for (j=1;j < sim->wl.length[1];j++) {
                sim->wl.weights[i+j*sim->wl.length[0]] -= sim->wl.wmin;
            }
        }
        sim->wl.write(files->wloutfile);
        sim->wl.end();
        if ( (sim->wlm[0] == 2)||(sim->wlm[1] == 2) ) {
            sim->wl.mesh.end();
            sim->wl.origmesh.end();
        }
        if ( (sim->wlm[0] == 5)||(sim->wlm[1] == 5)||(sim->wlm[0] == 6)||(sim->wlm[1] == 6)  ) {
            if ( sim->wl.radiushole != NULL ) free(sim->wl.radiushole);
            if ( sim->wl.radiusholeold != NULL ) free(sim->wl.radiusholeold);
        }
    }
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
    //struct stat shapex, shapey, shapez;   // Box shape statistics
    //struct stat smec[MAXF];       // Smectic order parameters (Fourier coeeficients)

    FILE *mf=NULL;                                  // Handle for movie file
    FILE *cl_stat=NULL, *cl=NULL, *cl_list=NULL;    // Handle for cluster statistics
    FILE *ef=NULL, *statf=NULL;                     // Handle for energy file and statistical file

    double edriftstart;     // Energy drift calculation - start
    double edriftchanges;   // Energy drift calculation - accumulate all changes through moves
    double edriftend;       // Energy drift calculation - end
    double pvdriftstart;    // PV drift calculation - start
    double pvdriftend;      // PV drift calculation - end
    double volume;          // volume of box
    double moveprobab;      // random number selecting the move

    openFilesClusterStatistics(&cl_stat, &cl, &cl_list, &ef, &statf);

    //=== Initialise counters etc. ===//
    sim->shprob = sim->shave/(double)conf->particleStore.size();

    initValues(next_adjust, next_calc, next_dump, next_frame);

    if (sim->movie > 0) {
        mf = fopen(files->moviefile, "a");
    } else {
        mf = NULL;
    }

    initWangLandau();

    //do moves - START OF REAL MC
    if(sim->pairlist_update){
        genPairList(); // Does that solve the problem?
    }

    //do energy drift check - start calculation
    volume = conf->box.x * conf->box.y * conf->box.z;
    edriftstart = calcEnergy(0, 0, 0);
    pvdriftstart = sim->press * volume - (double)conf->particleStore.size() * log(volume) / sim->temper;
    //printf("starting energy: %.15f \n",calc_energy(0, intfce, 0, topo, conf, sim,0));
    //printf("press: %.15f\n",sim->press * volume - (double)conf->particleStore.size() * log(volume) / sim->temper);
    edriftchanges = 0.0;

    for (sweep=1; sweep <= nsweeps; sweep++) {

        // Try replica exchange
        if((sim->nrepchange) && (sweep % sim->nrepchange == 0)){
            edriftchanges += move.replicaExchangeMove(sweep);

            if(sim->pairlist_update)
                genPairList();
        }

        // Try muVT insert delete moves
        if(sim->nGrandCanon != 0 && sweep%sim->nGrandCanon == 0) {
            edriftchanges += move.muVTMove();

            if(sim->pairlist_update)
                genPairList();
        }

        // Generate the pairlist, also generate after each muVT move
        if( (sim->pairlist_update) && // pair_list allowed
                (
                    (sweep % sim->pairlist_update == 0) && // on scheduled sweep
                    !(sim->nGrandCanon != 0 && sweep%sim->nGrandCanon == 0) && // not on grandCanon sweep
                    !((sim->nrepchange) && (sweep % sim->nrepchange == 0))  // not on replica exchange sweep
                ) ) {
            genPairList();
        }

        //normal moves
        for (step=1; step <= (long)(long)conf->particleStore.size(); step++) {
            moveprobab = ran2();
            if ( moveprobab < sim->shprob) {
                // pressure moves
                edriftchanges += move.pressureMove();
            } else {
                if (moveprobab < sim->shprob + sim->chainprob) {
                    // single particle moves
                    edriftchanges += move.chainMove();
                }
                else if (moveprobab < sim->shprob + sim->chainprob + sim->switchprob){
                    //=== This is an attempt to switch a type ===
                    edriftchanges += move.switchTypeMove();

                } else {
                    // single particle moves
                    edriftchanges += move.particleMove();

                } // end of else next to chain moves
            } // end of else next to volume moves        
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

        if ( (sim->wlm[0] > 0) && (sim->wl.alpha > WL_ZERO) && !(sweep % 1000) ) {
            // recalculate system CM to be sure there is no accumulation of errors by +- rejection moves
            if ( (sim->wlm[0] == 1) || (sim->wlm[1] == 1) )
                  conf->massCenter(topo);
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
                      printf (" %15.8e %15ld %15.8f\n",sim->wl.weights[i],sim->wl.hist[i],particleStore[0].pos.z);
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
            conf->partVecInit(topo);
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
              accumulate (&shapex, (*box).x);
              accumulate (&shapey, (*box).y);
              accumulate (&shapez, (*box).z);
              volume = (*box).x * (*box).y * (*box).z;
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
              printf ("   Mean & fluc box dimensions:  x   %13.8f %13.8f\n",
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

            fprintf (statf, " %ld; %.10f\n", sweep, conf->box.x * conf->box.y * conf->box.z);
            fprintf (ef, " %ld; %.10f  %f \n", sweep, calcEnergy(0, 0, 0), alignmentOrder());
            if (sim->wlm[0] > 0) {
                sim->wl.write(files->wloutfile);
            }
            //print mesh distribution
            //mesh_findholesdistrib(&sim->wl.mesh);
            next_dump += report;
        }

        // Writing of movie frame
        if (sweep == next_frame) {
            fprintf (mf, "%ld\n", (long)conf->particleStore.size());
            fprintf (mf, "sweep %ld;  box %.10f %.10f %.10f\n", sweep, conf->box.x, conf->box.y, conf->box.z);
            printStat::draw(mf, conf);
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
    volume = conf->box.x * conf->box.y * conf->box.z;
    edriftend = calcEnergy(0, 0, 0);
    pvdriftend =  sim->press * volume - (double)conf->particleStore.size() * log(volume) / sim->temper;
    printf("Energy drift: %.5e \n",edriftend - edriftstart - edriftchanges +pvdriftend -pvdriftstart);
    printf("Starting energy: %.8f \n",edriftstart);
    printf("Starting energy+pv: %.8f \n",edriftstart+pvdriftstart);
    printf("System:\n");
    for(i=0; i < conf->molTypeCount; i++)
        printf("%s %d\n", topo->chainparam[i].name, conf->molCountOfType(i));
    fflush(stdout);

    endWangLandau();

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

    for (i = 0; i < (long)(long)conf->particleStore.size() - 1; i++) {
        for (j = i + 1; j < (long)(long)conf->particleStore.size(); j++) {
            r_cm = image(&conf->particleStore[i].pos, &conf->particleStore[j].pos, &conf->box);
            if ( DOT(r_cm,r_cm) < 1.5*1.5 ) {
                sumdot+= DOT(conf->particleStore[i].dir,conf->particleStore[j].dir);
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

    for(unsigned int i = 0; i < (long)(long)conf->neighborList.size(); i++){
        //DEBUG_INIT("%ld", i);
        conf->neighborList[i].neighborCount = 0;
    }
    for(unsigned int i = 0; i < conf->particleStore.size()-1; i++){
        for(unsigned int j = i + 1; j < conf->particleStore.size(); j++){
            assert(conf->particleStore.size() == conf->neighborList.size());
            assert(conf->particleStore[i] != NULL);
            assert(conf->particleStore[j] != NULL);

            r_cm.x = conf->particleStore[i].pos.x - conf->particleStore[j].pos.x;
            r_cm.y = conf->particleStore[i].pos.y - conf->particleStore[j].pos.y;
            r_cm.z = conf->particleStore[i].pos.z - conf->particleStore[j].pos.z;
            if ( r_cm.x < 0  )
                r_cm.x = conf->box.x * (r_cm.x - (double)( (long)(r_cm.x-0.5) ) );
            else
                r_cm.x = conf->box.x * (r_cm.x - (double)( (long)(r_cm.x+0.5) ) );
            if ( r_cm.y < 0  )
                r_cm.y = conf->box.y * (r_cm.y - (double)( (long)(r_cm.y-0.5) ) );
            else
                r_cm.y = conf->box.y * (r_cm.y - (double)( (long)(r_cm.y+0.5) ) );
            if ( r_cm.z < 0  )
                r_cm.z = conf->box.z * (r_cm.z - (double)( (long)(r_cm.z-0.5) ) );
            else
                r_cm.z = conf->box.z * (r_cm.z - (double)( (long)(r_cm.z+0.5) ) );

            r_cm2 = DOT(r_cm,r_cm);
            max_dist = AVER(sim->trans[conf->particleStore[i].type].mx, \
                    sim->trans[conf->particleStore[j].type].mx);
            max_dist *= (1 + sim->pairlist_update) * 2;
            max_dist += topo->maxcut;
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
    //for(i = 0; i < (long)conf->particleStore.size(); i++){
    //    //if (sim->pairlist.list[i].num_pairs >= (long)conf->particleStore.size())
    //    if (sim->pairlist[i].num_pairs >= (long)conf->particleStore.size()){
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
                    sweep, (long)(long)conf->particleStore.size());
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
    for( i = 0; i < (long)conf->particleStore.size(); i++){
        sim->clusterlist[i] = i;
    }

    // Start determining the cluster
    while(change){
        change = false;
        assert(conf->particleStore.size() == conf->neighborList.size());
        for(i = 0; i < (long)conf->particleStore.size(); i++){
            /*If nore pairlist go over all pairs*/
            maxnumber = (long)conf->particleStore.size();
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
                if ( (topo->ia_params[conf->particleStore[fst].type][conf->particleStore[snd].type].geotype[0] < SP) && \
                    (topo->ia_params[conf->particleStore[fst].type][conf->particleStore[snd].type].geotype[1] < SP) ) {

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
    long cluster_indices[(long)conf->particleStore.size()];   /* holds the different cluster indices.
                        (currently too much memory) */
    long num_cluster = 0;                /* number of clusters, temporary needed */

    /* how many clusters are there? */
    long max_index = -1;
    for(int i = 0; i < (long)conf->particleStore.size(); i++){
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
        sim->clusters[i].particles = (long int*) malloc(sizeof(long) * (long)conf->particleStore.size());

        if (!sim->clusters[i].particles){
            fprintf(stderr, "Couldn't allocate any memory!\n");
            exit(1);
        }

        sim->clusters[i].npart = 0;
    }

    /* fill in the particles belonging to one cluster */
    for(int i = 0; i < num_cluster; i++){
        for(int j = 0; j < (long)conf->particleStore.size(); j++){
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
                sim->clustersenergy[i]+= (calcEnergy.pairE)(&conf->particleStore[sim->clusters[i].particles[j]]
                        , sim->clusters[i].particles[j]
                        , &conf->particleStore[sim->clusters[i].particles[k]]
                        , sim->clusters[i].particles[k]);
            }
        }
    }
    return 0;
}



int Updater::sameCluster(long fst, long snd) {

    /*if two particles are bonded they belong to the same cluster*/
    if ( ((topo->chainparam[conf->particleStore[fst].molType]).bond1c >= 0) ||
        ((topo->chainparam[conf->particleStore[fst].molType]).bonddc >= 0) ){
        if ( (snd == conf->neighborList[fst].conlist[1]) || (snd == conf->neighborList[fst].conlist[0]) ) {
          return true;
        }
    }
    if ( ((topo->chainparam[conf->particleStore[snd].molType]).bond1c >= 0) ||
        ((topo->chainparam[conf->particleStore[snd].molType]).bonddc >= 0) ){
        if ( (fst == conf->neighborList[snd].conlist[1]) || (fst == conf->neighborList[snd].conlist[0]) ) {
          return false;
        }
    }

    /*cluster is made of particles closer tna some distance*/
/*	struct vector2 image(struct vector2 r1, struct vector2 r2, struct vector2 box);
    struct vector2 r_cm = image(conf->particleStore[fst].pos,
            conf->particleStore[snd].pos,
            conf->box);
    double dist2 = DOT(r_cm, r_cm);
    * TODO: Make it much more efficient => define cluster_dist!!! *
    if(dist2 > topo->ia_params[conf->particleStore[fst].type][conf->particleStore[snd].type].sigma * topo->ia_params[conf->particleStore[fst].type][conf->particleStore[snd].type].sigma*4.0){
        return false;
    }
    else {
        return true;
    }*/

    /*cluster is made of attractively interacting particles*/
    /*double paire(long, long, double (* intfce[MAXT][MAXT])(struct interacts *),
            struct topo * topo, struct conf * conf); Redeclaration*/

    if((calcEnergy.pairE)(&conf->particleStore[fst], fst, &conf->particleStore[snd], snd) > -0.10 ){
        return false;
    }
    else {
        return true;
    }
}

