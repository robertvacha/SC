#include "updater.h"

#include <iomanip>
#include <algorithm>

void Updater::emptyFiles() {

    FILE* file = NULL;
    if(sim->write_cluster){
        // Empty file
        file = fopen(files->clusterstatfile, "w");
        fclose(file);
        // Empty file
        file = fopen(files->clusterfile, "w");
        fclose(file);
    }
    if (report <= nsweeps){
        file = fopen(files->energyfile, "w");
        fclose(file);

        file = fopen(files->statfile, "w");
        fclose(file);
    }
}



void Updater::initValues() {

    // Initialize some values at begining
    conf->partVecInit();
    next_adjust = adjust;
    next_calc = paramfrq;
    next_dump = report;
    next_frame = sim->movie;
    //nem = vol = shapex = shapey = shapez = nullstat;
    //for (i=0; i<MAXF; i++) smec[i] = nullstat;

    wl.wl_meshsize = 0;
    wl.radiushole = NULL;
    wl.radiusholeold = NULL;
    wl.radiusholemax = 0;
    wl.partincontactold = 0;
    wl.partincontact = 0;
    wl.wlmdim = 0;
    wl.wlmdim = 0;
    wl.length[0]=0;
    wl.length[1]=0;
    wl.currorder[0]=0;
    wl.currorder[1]=0;
    wl.neworder[0]=0;
    wl.neworder[1]=0;
    wl.weights = NULL;
    wl.hist = NULL;

    conf->massCenter();
}




void Updater::simulate(long nsweeps, long adjust, long paramfrq, long report) {

    //cout << (calcEnergy.allToAll()/sim->temper) / (4.0 / 1.38064852 / 6.02214086 *1e3) / 7.5 *300  << " kJ/mol" << endl;
    //cout << calcEnergy.allToAll()/(4*sim->temper) * 1.38064852 * 6.02214086 /1000 / 7.5   << " kJ/mol" << endl;
    //cout << calcEnergy.allToAll()/(4*sim->temper) << " kT" << endl;
    //cout << calcEnergy.allToAll()/sim->temper << endl;
    //cout << calcEnergy.allToAll() << endl;
    bool mpi = false;
#ifdef ENABLE_MPI
    mpi = true;
#endif
    long i,j;

    size_t time = 0;
    size_t temp = 0;

    this->nsweeps = nsweeps;
    this->adjust = adjust;
    this->paramfrq = paramfrq;
    this->report = report;

    long step;         // Step number within a given sweep
    long sweep;        // Current sweep number

    //struct stat nem;   // Nematic order parameter
    //struct stat vol;   // Volume statistic
    //struct stat shapex, shapey, shapez;   // geo.box shape statistics
    //struct stat smec[MAXF];       // Smectic order parameters (Fourier coeeficients)

    emptyFiles();

    FILE* ef;
    FILE* statf;
    ef = fopen(files->energyfile, "a");
    fprintf (ef, "# sweep    energy\n");
    fclose(ef);

    statf = fopen(files->statfile, "a");
    fprintf (statf, "# sweep    volume\n");
    fclose(statf);


    //=== Initialise counters etc. ===//
    if(!conf->pvec.empty()) {
        sim->shprob = sim->shave/(double)conf->pvec.size();
    } else {
        sim->shprob = 0.0;
    }

    initValues();

    wl.init(files->wlinfile);
    //do moves - START OF REAL MC
    if(sim->pairlist_update){
        temp = clock();
        genPairList();
        sim->pairList += clock() - temp;
    }

    double edriftchanges = 0.0;   // Energy drift calculation - accumulate all changes through move
    double edriftend;       // Energy drift calculation - end
    double pvdriftend;      // PV drift calculation - end
    double moveprobab;      // random number selecting the move
    double edriftstart = 0;

    edriftstart = calcEnergy.allToAll(conf->energyMatrix);     // Energy drift calculation - start

    double volume = conf->geo.volume();          // volume of geo.box
    const double pvdriftstart = sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper;    // PV drift calculation - start

    //printf("starting energy: %.15f \n",calc_energy(0, intfce, 0, topo, conf, sim,0));
    //printf("press: %.15f\n",sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper);

    /********************************************************/
    /*                 Simulation Loop                      */
    /********************************************************/
    time = clock();
    for (sweep=1; sweep <= nsweeps; sweep++) {
        if(nsweeps>=10 && sweep%(nsweeps/10) == 0 && !mpi) {
            volume = conf->geo.volume();
            edriftend = calcEnergy.allToAll();
            pvdriftend =  sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper;
            time = clock()-time;
            cout << "sweep: " << sweep << " particles: " << conf->pvec.size()
                 << " drift: " << edriftend - edriftstart - edriftchanges +pvdriftend -pvdriftstart;
            cout <<"\nInteractionEnergy: " << edriftend;
            cout << ", sweeps per hour: " << (3600.0)/((double)time/CLOCKS_PER_SEC)*nsweeps/10<< "\n" << endl;
            time = clock();
        }

#ifdef ENABLE_MPI
    // receive MPI data
#endif

        //____________Replica Exchange Move____________
        if((sim->nrepchange) && (sweep % sim->nrepchange == 0)){

            edriftchanges += move.replicaExchangeMove(sweep); // sending reference to pvdriftstart
            files->initMPIRank(sim->pseudoRank);

            if(sim->pairlist_update) {
                temp = clock();
                genPairList();
                sim->pairList += clock() - temp;
            }
        }
        //____________GrandCanonical Move____________
        if(sim->nGrandCanon != 0 && sweep%sim->nGrandCanon == 0) {
            unsigned int size = conf->pvec.size();
            edriftchanges += move.muVTMove();

            if(size != conf->pvec.size())
                calcEnergy.allToAll(conf->energyMatrix);

            if(sim->pairlist_update) {
                temp = clock();
                genPairList();
                sim->pairList += clock() - temp;
            }
        }
        //____________Cluster Move____________
        if(sim->nClustMove != 0 && sweep%sim->nClustMove == 0) {
            edriftchanges += move.clusterMove();

            if(sim->pairlist_update) {
                temp = clock();
                genPairList();
                sim->pairList += clock() - temp;
            }
        }
        if( (sim->pairlist_update) && // pair_list allowed
                (
                    (sweep % sim->pairlist_update == 0) && // on scheduled sweep
                    !(sim->nGrandCanon != 0 && sweep%sim->nGrandCanon == 0) && // not on grandCanon sweep
                    !(sim->nClustMove != 0 && sweep%sim->nClustMove == 0) && // not on Cluster move sweep
                    !((sim->nrepchange) && (sweep % sim->nrepchange == 0))  // not on replica exchange sweep
                )
                ) {
            temp = clock();
            genPairList();
            sim->pairList += clock() - temp;
        }
        //normal moves
        for (step=1; step <= (long)conf->pvec.size(); step++) {
            moveprobab = ran2();

            //assert(testEnergyMatrix()); // Super expensive

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
            //if(conf->checkall(topo.ia_params)){
            //    cout<<"OVERLAP DETECTED"<<endl;
            //}
        } // End of step loop for this sweep

        //=== Start of end-of-sweep housekeeping ===
        // Adjustment of maximum step sizes during equilibration
        if (sweep == next_adjust) {
            for (i = 0; i < MAXT ;i++) {
                if ((sim->stat.trans[i].acc > 0)||(sim->stat.trans[i].rej >0))
                    optimizeStep (sim->stat.trans + i, 1.5, 0.0);
                if ((sim->stat.rot[i].acc > 0)||(sim->stat.rot[i].rej >0))
                    optimizeRot (sim->stat.rot + i, 5.0, 0.01);
            }
            for (i = 0; i < MAXMT; i++) {
                if ((sim->stat.chainm[i].acc > 0)||(sim->stat.chainm[i].rej > 0))
                    optimizeStep (sim->stat.chainm + i, 1.5, 0.0);
                if ((sim->stat.chainr[i].acc > 0)||(sim->stat.chainr[i].rej > 0))
                    optimizeRot (sim->stat.chainr + i, 5.0, 0.01);
            }
            optimizeStep (&(sim->stat.edge), 1.0, 0.0);
            next_adjust += adjust;
        }

        if ( (wl.wlm[0] > 0) && (wl.alpha > WL_ZERO) && !(sweep % 1000) ) {
            // recalculate system CM to be sure there is no accumulation of errors by +- rejection moves
            /* BUG - not used any longer: caused problems with PBC normal moves systemCM movement
              can be calculated from CM movements of individual particles
              present center of mass calculation use pbc and thus particles that moved across the geo.box
              is in this calculation used in pripary geo.box but in other moves in in in the particles position
             if ( (sim->wlm[0] == 1) || (sim->wlm[1] == 1) )
                  masscenter(topo.npart,topo.ia_params, conf);
            */
            wl.min = wl.hist[0];
            wl.max = wl.hist[0];
            for (i=0;i < wl.length[0];i++) {
                j=0;
                if ( wl.hist[i+j*wl.length[0]] > wl.max ) wl.max = wl.hist[i+j*wl.length[0]];
                if ( wl.hist[i+j*wl.length[0]] < wl.min ) wl.min = wl.hist[i+j*wl.length[0]];
                for (j=1;j < wl.length[1];j++) {
                    if ( wl.hist[i+j*wl.length[0]] > wl.max ) wl.max = wl.hist[i+j*wl.length[0]];
                    if ( wl.hist[i+j*wl.length[0]] < wl.min ) wl.min = wl.hist[i+j*wl.length[0]];
                }
            }
            if ( wl.min > WL_MINHIST ) {
                if ( sim->temper * log(wl.max/wl.min) < WL_GERR ) {
                    /*DEBUG
                      for (i=1;i<wl.length;i++) {
                      printf (" %15.8e %15ld %15.8f\n",wl.weights[i],wl.hist[i],pvec[0].pos.z);
                      fflush(stdout);
                      }
                     */
                    if ( wl.alpha < WL_ALPHATOL) break;
                    wl.alpha/=2;
                    printf("%f \n", wl.alpha);
                    fflush (stdout);
                    wl.wmin = wl.weights[0];

                    for (i=0;i < wl.length[0];i++) {
                        j=0;
                        wl.hist[i+j*wl.length[0]] = 0;
                        wl.weights[i+j*wl.length[0]] -= wl.wmin;
                        for (j=1;j < wl.length[1];j++) {
                            wl.hist[i+j*wl.length[0]] = 0;
                            wl.weights[i+j*wl.length[0]] -= wl.wmin;
                        }
                    }

                }
            }
        }

        if (!(sweep % 100000)) { //reinitialize patch vectors to avoid cummulation of errors
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

            ef = statf = NULL;
            ef = fopen(files->energyfile, "a");
            statf = fopen(files->statfile, "a");

            fprintf (statf, " %ld; %.10f\n", sweep, conf->geo.box.x * conf->geo.box.y * conf->geo.box.z);
            fprintf (ef, " %ld; %.10f  %f \n", sweep, calcEnergy(0, 0, 0), alignmentOrder());
            if (wl.wlm[0] > 0) {
                wl.write(files->wloutfile);
            }
            //print mesh distribution
            //mesh_findholesdistrib(&wl.mesh);
            next_dump += report;

            fclose(ef);
            fclose(statf);
            ef = statf = NULL;
        }

        // Writing of movie frame
        if (sweep == next_frame)
            dumpMovie(sweep);


        // Writing out cluster statistics
        if(sim->write_cluster && (sweep % sim->write_cluster == 0)){
            writeCluster(false, sweep);
        }
        //=== End of housekeeping ===

    }  // End of sweeps loop

    //do energy drift check - at the end calculation
    volume = conf->geo.volume();
    edriftend = calcEnergy.allToAll();
    pvdriftend =  sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper;
    if(sim->mpinprocs > 1) {
        printf("%d Energy drift: %.5e \n", sim->pseudoRank, edriftend - edriftstart - edriftchanges +pvdriftend -pvdriftstart);
        /*printf("%d Starting energy: %.8f \n", sim->pseudoRank, edriftstart);
        printf("%d Ending energy: %.8f \n", sim->pseudoRank, edriftend);
        printf("%d EdriftChanges: %.5e\n", sim->pseudoRank, edriftchanges);
        printf("%d Starting pV: %.8f \n", sim->pseudoRank, pvdriftstart);
        printf("%d Ending pV: %.8f \n", sim->pseudoRank, pvdriftend);*/
    }
    else {
        printf("Energy drift: %.5e \n",edriftend - edriftstart - edriftchanges +pvdriftend -pvdriftstart);
        printf("Starting energy+pv: %.8f \n",edriftstart+pvdriftstart);
        printf("Starting energy: %.8f \n",edriftstart);
        printf("Ending energy: %.8f \n",edriftend);

        cout << std::setprecision(2) << std::left;
        cout << "\n\n******************************************************************************" << endl;
        cout << "*                               Moves Statistics                             *" << endl;
        cout << "******************************************************************************" << endl;
        cout << setw(30) << "Move" << setw(10) << "Acc (%)" << setw(10) << "Rej (%)" << setw(10) << "Steps" << endl;

        if(sim->stat.stepsSTrans() > 0) {
            cout << setw(30) << "Single particle translation: "
                 << setw(10) <<  (double)sim->stat.accSTrans()/sim->stat.stepsSTrans()*100.0
                 << setw(10) << (double)sim->stat.rejSTrans()/sim->stat.stepsSTrans()*100.0
                 << setw(10) << sim->stat.stepsSTrans() << endl;
        }

        if(sim->stat.stepsSRot() > 0) {
            cout << setw(30) << "Single particle rotation: "
                 << setw(10) << 100.0*sim->stat.accSRot()/sim->stat.stepsSRot()
                 << setw(10) << (double)sim->stat.rejSRot()/sim->stat.stepsSRot()*100.0
                 << setw(10) << sim->stat.stepsSRot() << endl;
        }

        if(sim->stat.stepsCTrans() > 0) {
            cout << setw(30) << "Chain translation: "
                 << setw(10) <<  (double)sim->stat.accCTrans()/sim->stat.stepsCTrans()*100.0
                 << setw(10) << (double)sim->stat.rejCTrans()/sim->stat.stepsCTrans()*100.0
                 << setw(10) << sim->stat.stepsCTrans() << endl;
        }

        if(sim->stat.stepsCRot() > 0) {
            cout << setw(30) << "Chain rotation: "
                 << setw(10) <<  (double)sim->stat.accCRot()/sim->stat.stepsCRot()*100.0
                 << setw(10) << (double)sim->stat.rejCRot()/sim->stat.stepsCRot()*100.0
                 << setw(10) << sim->stat.stepsCRot() << endl;
        }

        if(sim->stat.edge.acc + sim->stat.edge.rej > 0) {
            cout << setw(30) << "Pressure move: "
                 << setw(10) << (double)sim->stat.edge.acc / (sim->stat.edge.acc + sim->stat.edge.rej)*100.0
                 << setw(10) << (double)sim->stat.edge.rej / (sim->stat.edge.acc + sim->stat.edge.rej)*100.0
                 << setw(10) << (sim->stat.edge.acc + sim->stat.edge.rej) << endl;
        }

        if(sim->stat.stepsSwitch() > 0) {
            cout << setw(30) << "Switch type move: "
                 << setw(10) <<  (double)sim->stat.accSwitch()/sim->stat.stepsSwitch()*100.0
                 << setw(10) << (double)sim->stat.rejSwitch()/sim->stat.stepsSwitch()*100.0
                 << setw(10) << sim->stat.stepsSwitch() << endl;
        }

        if(sim->nGrandCanon != 0) {
            for(int i=0; i<conf->pvec.molTypeCount; i++) {
                cout << setw(20) << "Insert move of type " << setw(10) << topo.moleculeParam[i].name
                     << setw(10) << (double) sim->stat.grand[i].insAcc / (sim->stat.grand[i].insAcc + sim->stat.grand[i].insRej)*100.0
                     << setw(10) << (double) sim->stat.grand[i].insRej / (sim->stat.grand[i].insAcc + sim->stat.grand[i].insRej)*100.0
                     << setw(10) << (sim->stat.grand[i].insAcc + sim->stat.grand[i].insRej) << endl;

                cout << setw(20) << "Remove move of type " << setw(10) << topo.moleculeParam[i].name
                     << setw(10) << (double) sim->stat.grand[i].delAcc / (sim->stat.grand[i].delAcc + sim->stat.grand[i].delRej)*100.0
                     << setw(10) << (double) sim->stat.grand[i].delRej / (sim->stat.grand[i].delAcc + sim->stat.grand[i].delRej)*100.0
                     << setw(10) << (sim->stat.grand[i].delAcc + sim->stat.grand[i].delRej) << endl;

                cout << setw(20) << "Average particles of type " << setw(10) << topo.moleculeParam[i].name
                     << setw(10) << (double) sim->stat.grand[i].muVtAverageParticles / sim->stat.grand[i].muVtSteps << endl;
            }
        }

        //printf("EdriftChanges: %.5e\n", edriftchanges);
    }

//    printf("PVChanges: %.5e\n", pvdriftend -pvdriftstart);

    printf("System:\n");

    for(i=0; i < conf->pvec.molTypeCount; i++)
        printf("%s %d\n", topo.moleculeParam[i].name, conf->pvec.molCountOfType(i));

    fflush(stdout);

    wl.endWangLandau(files->wloutfile);
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

    newrmsd = (*x).mx * x->ratio() ;
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

void Updater::genSimplePairList() {
    Vector r_cm;
    double r_cm2;
    // Set the pairlist to zero
    //DEBUG_INIT("Gen Pairlist")
    for(unsigned int i = 0; i < conf->neighborList.size(); i++){
        //DEBUG_INIT("%ld", i);
        conf->neighborList[i].neighborCount = 0;
    }

    if(conf->pvec.size() <= 0)
        return;

    for(unsigned int i = 0; i < conf->pvec.size()-1; i++){
        for(unsigned int j = i + 1; j < conf->pvec.size(); j++){
            assert(conf->pvec.size() == conf->neighborList.size());
            assert(&conf->pvec[i] != NULL);
            assert(&conf->pvec[j] != NULL);

            r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[j].pos);

            r_cm2 = DOT(r_cm,r_cm);
#ifndef NDEBUG
            double max_dist;
            max_dist = AVER(sim->stat.trans[conf->pvec[i].type].mx, \
                    sim->stat.trans[conf->pvec[j].type].mx);
            max_dist *= (1 + sim->pairlist_update) * 2;
            max_dist += topo.maxcut;
            max_dist *= max_dist; /* squared */
#endif
            assert(max_dist == sim->max_dist_squared[conf->pvec[i].type][conf->pvec[j].type]);

            if (r_cm2 <= sim->max_dist_squared[conf->pvec[i].type][conf->pvec[j].type]){
                conf->neighborList[i].neighborID[conf->neighborList[i].neighborCount] = j;
                conf->neighborList[j].neighborID[conf->neighborList[j].neighborCount] = i;

                /*conf->neighborList[i].energyID[conf->neighborList[i].neighborCount] = calcEnergy.p2p(i,j);
                conf->neighborList[j].energyID[conf->neighborList[j].neighborCount] = conf->neighborList[i].energyID[conf->neighborList[i].neighborCount];*/

                conf->neighborList[i].neighborCount++;
                conf->neighborList[j].neighborCount++;
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


int Updater::writeCluster(bool decor, long sweep) {

    FILE *cl_stat = NULL;
    FILE *cl = NULL;

    cl_stat = fopen(files->clusterstatfile, "a");
    cl = fopen(files->clusterfile, "a");

    genClusterList();
    sortClusterList();
    calcClusterEnergies();
    if(cl_stat){
        if(decor == false){
            // if no decor, this means usually into a file. Hence print info
            // about number of line per frame
            fprintf(cl_stat, "Sweep: %ld | Maximal size: %ld\n", sweep, sim->max_clust);
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
    FILE *cl_list = NULL; // CL_LIST NEVER USED
    if(cl_list){
        if(decor == false){
            fprintf(cl_list, "Sweep: %ld | Number of particles: %ld\n", sweep, (long)conf->pvec.size());
        }
        printStat::printClusterList(cl, decor, sim, conf);
    }

    fclose(cl_stat);
    fclose(cl);
    //fclose(cl_list);

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
                //if ( (topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].geotype[0] < SP) && \
                 //   (topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].geotype[1] < SP) ) {

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
                //}
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
    sim->clusterstat = (long int*) realloc( sim->clusterstat, sizeof(long) * sim->max_clust); // OLD, MISTAKE? memmory dont have to be 0, no free
    memset(sim->clusterstat, 0, sizeof(long) * sim->max_clust);

    if (!sim->clusterstat){
        fprintf(stderr, "Couldn't allocate any memory!\n");
        exit(1);
    }

    for(int i = 0; i < sim->max_clust; i++) {
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

    /*Here we add particle into cluster if particles fst and snd are from same molecule*/
    Molecule fstMol = conf->pvec.getMolOfPart(fst);
    if( std::find(fstMol.begin(), fstMol.end(), snd) != fstMol.end() ){
        return true;
    }

    if(calcEnergy.p2p(fst, snd) > -0.10 ) {
        return false;
    }
    else {
        return true;
    }
}

