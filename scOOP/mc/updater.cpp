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

    conf->massCenter();
}




void Updater::simulate(long nsweeps, long adjust, long paramfrq, long report) {

    if(sim->nGrandCanon != 0) {
        cout << "\n!!!!!!!!!!!! WARNING !!!!!!!!!!!!\n" << endl;
        cout << "Energy matrix for GrandCanonical move is not implemented efficiently" << endl;
        cout << "Go to end of file totalenergycalculator.h" << endl;
        cout << "Comment out typedef TotalEMatrix<PairE> TotalEnergyCalculator" << endl;
        cout << "Uncomment typedef TotalEFull<PairE> TotalEnergyCalculator for Full calculation\n" << endl;
    }

    //cout << (calcEnergy.allToAll()/sim->temper) / (4.0 / 1.38064852 / 6.02214086 *1e3) / 7.5 *300  << " kJ/mol" << endl;
    //cout << calcEnergy.allToAll()/(4*sim->temper) * 1.38064852 * 6.02214086 /1000 / 7.5   << " kJ/mol" << endl;
    //cout << calcEnergy.allToAll()/(4*sim->temper) << " kT" << endl;
    //cout << calcEnergy.allToAll()/sim->temper << endl;
    //cout << calcEnergy.allToAll() << endl;

    size_t time = 0;
    size_t temp = 0;

    this->nsweeps = nsweeps;
    this->adjust = adjust;
    this->paramfrq = paramfrq;
    this->report = report;

    long step;         // Step number within a given sweep
    long sweep;        // Current sweep number

    //
    // Erase content of energy.dat, stat.dat, clusterstat and clusterfile
    //
    emptyFiles();

    FILE* ef;
    ef = fopen(files->energyfile, "a");
    fprintf (ef, "# sweep    energy\n");
    fclose(ef);

    FILE* statf;
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

    if(showPairInteractions) {
        for(int i=0; i< conf->pvec.size(); ++i) {
            conf->pvec[i].testInit(PSC, i);
        }
        double e;
        // chain
        int len = 2;
        for(int i=0; i< conf->pvec.size()-1; i += len) {
            e = calcEnergy.p2p(i,i+1);
            /*if(e < 1000.0)
                printf("%.5lf\n", e);
            else printf("%lf\n", 1000.0);*/
        }
        exit(0);
        //single part
        for(int i=0; i< 1/* conf->pvec.size()-1*/; ++i) {
            for(int j=i+1; j< conf->pvec.size(); ++j) {
                e = calcEnergy.p2p(i,j);
                if(e < 1000.0)
                    printf("%.5lf\n", e);
                else printf("%lf\n", 1000.0);
            }
            printf("\n");
        }
        exit(0);
    }

    move.wl.init(files->wlinfile);
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

    calcEnergy.initEM();
    edriftstart = calcEnergy.allToAll();     // Energy drift calculation - start

    double volume = conf->geo.volume();          // volume of geo.box
    const double pvdriftstart = sim->press * volume - (double)conf->pvec.size() * log(volume) * sim->temper;    // PV drift calculation - start

    /********************************************************/
    /*                 Simulation Loop                      */
    /********************************************************/
    time = clock();
    mcout.get() << std::left;
    mcout.get() << setw(10)  << "sweeps" << " "
                << setw(10) << "particles" << " "
                << setw(14) << "Drift" << " "
                << setw(14) << "PE" << " "
                << setw(8) << "Sweeps per hour" << endl;
    time = clock();
    for (sweep=1; sweep <= nsweeps; sweep++) {

        // SIMULATION LOG
        if(sim->thermo>0 && sweep%(sim->thermo) == 0) {
            volume = conf->geo.volume();
            edriftend = calcEnergy.allToAll();
            pvdriftend =  sim->press * volume - (double)conf->pvec.size() * log(volume) * sim->temper;
            time = clock()-time;
            mcout.get() << std::left;
            mcout.get() << setw(10)  << sweep << " "
                        << setw(10) << conf->pvec.size() << " "
                        << setw(14) << edriftend - edriftstart - edriftchanges +pvdriftend -pvdriftstart << " "
                        << setw(14) << edriftend << " "
                        << setw(8) << (3600.0)/((double)time/CLOCKS_PER_SEC)*nsweeps/10 << endl;
            time = clock();
        }

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

            //
            // Quick and dirty recalc of energy matrix, TODO: make it more efficient
            //
            if(size != conf->pvec.size())
                calcEnergy.initEM();

            assert(testEnergyMatrix() && "muVT move is messing up energy matrix");

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

        //if(sweep%100 == 0) {
        //    cout << sweep << " NEIGHBORS: " << calcEnergy.neighAvg/calcEnergy.div << " " << calcEnergy.interAvg/calcEnergy.div << " " << calcEnergy.interAvg2/calcEnergy.div << endl;
        //}

        //=== Start of end-of-sweep housekeeping ===
        // Adjustment of maximum step sizes during equilibration
        if (sweep == next_adjust) {
            for (int i = 0; i < MAXT ;i++) {
                if ((sim->stat.trans[i].acc > 0)||(sim->stat.trans[i].rej >0))
                    optimizeStep (sim->stat.trans[i], 1.5, 0.0);
                if ((sim->stat.rot[i].acc > 0)||(sim->stat.rot[i].rej >0))
                    optimizeRot (sim->stat.rot[i], 5.0, 0.01);
            }
            for (int i = 0; i < MAXMT; i++) {
                if ((sim->stat.chainm[i].acc > 0)||(sim->stat.chainm[i].rej > 0))
                    optimizeStep (sim->stat.chainm[i], 1.5, 0.0);
                if ((sim->stat.chainr[i].acc > 0)||(sim->stat.chainr[i].rej > 0))
                    optimizeRot (sim->stat.chainr[i], 5.0, 0.01);
            }
            optimizeStep ( sim->stat.edge, 1.0, 0.0);
            next_adjust += adjust;
        }

        if ( (move.wl.wlm[0] > 0) && (move.wl.alpha > WL_ZERO) && !(sweep % 1000) ) {
            // recalculate system CM to be sure there is no accumulation of errors by +- rejection moves
            /* BUG - not used any longer: caused problems with PBC normal moves systemCM movement
              can be calculated from CM movements of individual particles
              present center of mass calculation use pbc and thus particles that moved across the geo.box
              is in this calculation used in pripary geo.box but in other moves in in in the particles position
             if ( (sim->wlm[0] == 1) || (sim->wlm[1] == 1) )
                  masscenter(topo.npart,topo.ia_params, conf);
            */
            if(move.wl.update(sim->temper))
                break;
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
            fprintf (ef, " %ld; %.10f  %f \n", sweep, calcEnergy.allToAll(), alignmentOrder());

            if (move.wl.wlm[0] > 0 && mcout.rank == 0) { // Write WL file only for rank 0 file
                move.wl.write(files->wloutfile);
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
            clust.writeCluster(false, sweep);
        }
        //=== End of housekeeping ===

    }  // End of sweeps loop

    //do energy drift check - at the end calculation
    volume = conf->geo.volume();
    edriftend = calcEnergy.allToAll();
    pvdriftend =  sim->press * volume - (double)conf->pvec.size() * log(volume) * sim->temper;
    mcout.get() << "\nEnergy drift: " << edriftend - edriftstart - edriftchanges +pvdriftend -pvdriftstart << endl;
    mcout.get() << "Starting energy+pv: " << edriftstart+pvdriftstart << endl;
    mcout.get() << "Starting energy: " << edriftstart << endl;
    mcout.get() << "Ending energy: " << edriftend << endl;

    sim->stat.print();

    if(sim->nGrandCanon != 0) {
        for(int i=0; i<conf->pvec.molTypeCount; i++) {
            if(sim->stat.grand[i].insAcc + sim->stat.grand[i].insRej > 0) {
                mcout.get() << setw(20) << "Insert move of type " << setw(10) << topo.moleculeParam[i].name
                            << setw(10) << (double) sim->stat.grand[i].insAcc / (sim->stat.grand[i].insAcc + sim->stat.grand[i].insRej)*100.0
                            << setw(10) << (double) sim->stat.grand[i].insRej / (sim->stat.grand[i].insAcc + sim->stat.grand[i].insRej)*100.0
                            << setw(10) << (sim->stat.grand[i].insAcc + sim->stat.grand[i].insRej) << endl;

                mcout.get() << setw(20) << "Remove move of type " << setw(10) << topo.moleculeParam[i].name
                            << setw(10) << (double) sim->stat.grand[i].delAcc / (sim->stat.grand[i].delAcc + sim->stat.grand[i].delRej)*100.0
                            << setw(10) << (double) sim->stat.grand[i].delRej / (sim->stat.grand[i].delAcc + sim->stat.grand[i].delRej)*100.0
                            << setw(10) << (sim->stat.grand[i].delAcc + sim->stat.grand[i].delRej) << endl;

                mcout.get() << setw(20) << "Average particles of type " << setw(10) << topo.moleculeParam[i].name
                            << setw(10) << (double) sim->stat.grand[i].muVtAverageParticles / sim->stat.grand[i].muVtSteps << endl;
            }
        }
    }
    //printf("EdriftChanges: %.5e\n", edriftchanges);


    //    printf("PVChanges: %.5e\n", pvdriftend -pvdriftstart);

    mcout.get() << "\nSystem:" << endl;

    for(int i=0; i < conf->pvec.molTypeCount; i++)
        mcout.get() << topo.moleculeParam[i].name << " " << conf->pvec.molCountOfType(i) << endl;

    mcout.get() << endl;

    move.wl.endWangLandau(files->wloutfile);
}





void Updater::optimizeStep(Disp& x, double hi, double lo) {
    double newrmsd;

    newrmsd = x.mx * x.ratio();
    if (x.oldrmsd > 0) {
        if ( newrmsd  < x.oldrmsd ) {
            if ( x.oldmx > 1 ) {
                x.mx /= 1.05;
                x.oldmx = 0.95;
            } else {
                x.mx *= 1.05;
                x.oldmx = 1.05;
            }
        } else {
            if ( x.oldmx > 1 ) {
                x.mx *= 1.05;
                x.oldmx = 1.05;
            } else {
                x.mx /= 1.05;
                x.oldmx = 0.95;
            }
        }
    }
    if (newrmsd > 0 ) x.oldrmsd = newrmsd;
    else {
        x.oldrmsd = 0.0;
        x.mx /= 1.05;
        x.oldmx = 0.95;
    }

    if ( x.mx > hi ) x.mx = hi;
    if ( x.mx < lo ) x.mx = lo;

    x.acc = x.rej = 0;
}

void Updater::optimizeRot(Disp& x, double hi, double lo) {
    double newrmsd;

    newrmsd = x.mx * x.ratio() ;
    if (x.oldrmsd > 0) {
        if ( newrmsd  > x.oldrmsd ) {
            if ( x.oldmx > 1) {
                x.mx *= 0.99;
                x.oldmx *= 0.99;
            } else {
                x.mx *= 1.01;
                x.oldmx *= 1.01;
            }
        } else {
            if ( x.oldmx > 1) {
                x.mx *= 1.01;
                x.oldmx *= 1.01;
            } else {
                x.mx *= 0.99;
                x.oldmx *= 0.99;
            }
        }
    }
    if (newrmsd > 0 ) x.oldrmsd = newrmsd;
    else {
        x.oldrmsd = 0.0;
        x.mx *= 1.01;
        x.oldmx = 1.01;
    }

    if ( x.mx > hi ) x.mx = hi;
    if ( x.mx < lo ) x.mx = lo;

    x.acc = x.rej = 0;
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
    int pos=-1;
    double r_cm2;
    bool bonded = false;
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

            // Are bonds, Angles defined?
            if(i+2 >= j) {
                pos = (i - conf->pvec.first[conf->pvec[i].molType]) % topo.moleculeParam[ conf->pvec[i].molType ].molSize();

                if (topo.moleculeParam[conf->pvec[i].molType].bond1c >= 0.0 || topo.moleculeParam[conf->pvec[i].molType].bonddc >= 0.0 || topo.moleculeParam[conf->pvec[i].molType].bondhc >= 0.0) {
                    if(pos+1 < topo.moleculeParam[conf->pvec[i].molType].molSize() && i+1 == j) {
                        bonded = true;
                    }
                }
                if (topo.moleculeParam[conf->pvec[i].molType].bond2c >= 0.0) {
                    if(pos+2 < topo.moleculeParam[conf->pvec[i].molType].molSize() && i+2 == j ) {
                        bonded = true;
                    }
                }
            }

            if (r_cm2 <= sim->max_dist_squared[conf->pvec[i].type][conf->pvec[j].type] || bonded){
                conf->neighborList[i].neighborID[conf->neighborList[i].neighborCount] = j;
                conf->neighborList[j].neighborID[conf->neighborList[j].neighborCount] = i;

                conf->neighborList[i].neighborCount++;
                conf->neighborList[j].neighborCount++;
                bonded = false;
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
