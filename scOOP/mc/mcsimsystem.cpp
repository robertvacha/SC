/** @file mcsimsystem.cpp*/

#include "mcsimsystem.h"


void MCSimSystem::init(int argc, char** argv) {

    cout << "\nPatchy Spherocylinders version 3.6\n"
         << "-------------------------------------" << endl;

    Inicializer init(&topo, &sim, &conf, &files);

    init.initWriteFiles();

    init.initMPI(argc,argv);

    cout << "Reading options..." << endl;
    init.readOptions();

#ifdef EXTRA_HYDROPHOBIC_ALL_BODY_ATTRACTION
    printf("\n!!! Extra hydrophobic interaction in e_cpsc_cpsc added\n\n");
#endif

    init.initTop(); // here particleStore filled in setParticleParams
    init.testChains(); // if no chains -> move probability of chains 0

    cout << "\nReading configuration...\n";
    init.initConfig();

    cout << "Equilibration of maximum step sizes: " << sim.nequil/2 << " sweeps" << endl;

    clearOutfiles();

    if (sim.pairlist_update) {
        cout << "\nAllocating memory for pairlist..." << endl;
        init.initPairlist();
    }

    updater = new Updater(&topo, &sim, &conf, &files);

    sim.wl.setTopoConf(&topo, &conf);
    if(sim.pairlist_update)
        conf.pairlist_update = true;
    else conf.pairlist_update = false;
}




void MCSimSystem::equilibrate() {
    if (sim.nequil) {

        printf("\nStart equilibration...\n");

        updater->simulate(sim.nequil/2, sim.adjust, 0, 0);
        updater->simulate(sim.nequil/2, 0,          0, 0);

        printf ("   Equilibrated maximum displacement / acceptance ratio:            \n");
        printStat::printEqStat(sim.trans,2.0,MAXT);
        printf ("   Equilibrated maximum rotation / acceptance ratio:                       \n");
        printStat::printEqStat(sim.rot,1.0,MAXT);
        printf ("   Equilibrated maximum box length change / acceptance ratio:              \n");
        printf ("                     %.6e  /  %.6e\n", sim.edge.mx/2.0,RATIO(sim.edge));
        printf ("   Equilibrated maximum displacement of chain / acceptance ratio:   \n");
        printStat::printEqStat(sim.chainm,2.0,MAXMT);
        printf ("   Equilibrated maximum rotation of chain / acceptance ratio:              \n");
        printStat::printEqStat(sim.chainr,1.0,MAXMT);
        printf ("\n");

        printf ("Further equilibration of configuration:  %ld sweeps\n", sim.nequil/2);
        fflush (stdout);

        outfile = fopen("config.eq", "w");
        fprintf (outfile, "%15.8e %15.8e %15.8e\n", conf.box.x, conf.box.y, conf.box.z);
        printStat::draw(outfile, &conf);

        fclose (outfile);

        printf ("   Equilibrated configuration written to config.eq\n");
        printf ("   Box dimensions: %.10f, %.10f, %.10f\n\n", conf.box.x, conf.box.y, conf.box.z);
    }
}




void MCSimSystem::productionRun() {
    printf ("Production run:  %ld sweeps\n\n", sim.nsweeps);
    fflush (stdout);

    updater->simulate(sim.nsweeps, 0, sim.paramfrq, sim.report);

#ifdef ENABLE_MPI
        printf ("   MPI replica changeT / changeP / acceptance ratio: \t %.6f   /   %.6f  /  %.6f\n\n", sim.mpiexch.mx,sim.mpiexch.angle,RATIO(sim.mpiexch));
#endif
    outfile = fopen(files.configurationoutfile, "w");
    fprintf (outfile, "%15.8e %15.8e %15.8e\n", conf.box.x, conf.box.y, conf.box.z);
    printStat::draw (outfile, &conf);
    fclose (outfile);


    if(sim.nGrandCanon != 0) {
        FILE* inFile = fopen(files.topologyInFile, "r");
        outfile = fopen(files.topologyOutFile, "w");

        char line[128];

        while(strncmp(line, "[System]", 8) != 0) {
            if(fgets(line,127, inFile) == NULL ) {
                printf("Error writing Topology [System] not found\n");
                break;
            }
            fputs(line, outfile);
        }
        for(int i=0; i < conf.pvecGroupList.molTypeCount; i++)
            fprintf(outfile, "%s %d\n", topo.chainparam[i].name, conf.molCountOfType(i));

        fclose (outfile);
        fclose (inFile);
    }

    /// For testing the pairlist
    //gen_pairlist(&topo, &sim, &conf);
    //FILE * fpairlist;
    //fpairlist = fopen("pairlist.dat", "w");
    //print_pairlist(fpairlist, &sim, &topo);
    //fclose(fpairlist);
    //printf("sqmaxcut = %f\n", topo.sqmaxcut);

    /// For testing the cluster algorithm
    //gen_clusterlist(&topo, &sim, &conf);
    //print_clusterlist(stdout, TRUE, &topo, &sim, &conf);
    //sort_clusterlist(&topo, &sim);
    //print_clusters(stdout, TRUE, &sim);
    //print_clusterstat(stdout, TRUE, &sim);
}




void MCSimSystem::dealloc(){

    if (memoryDealloc())
        exit(1);
#ifdef ENABLE_MPI
    MPI_Finalize();
#endif
    printf ("\nDone\n\n");
}




/************************************************************************************************
 *                                      PRIVATE METHODS                                         *
 ************************************************************************************************/




void MCSimSystem::clearOutfiles() {
    if ( sim.wlm[0] > 0 ) {
        outfile = fopen(files.wlinfile, "r");
        if (outfile == NULL) {
            printf ("ERROR: Cannot open file for Wang-Landau method (%s).\n",files.wlinfile);
            memoryDealloc();
            exit(1);
        }
        fclose (outfile);
    }

    // Empty movie file
    mov = fopen("movie", "w");
    fclose (mov);
}




int MCSimSystem::memoryDealloc() {
    printf ("Deallocating memory...\n");

    for(unsigned int i=0; i < conf.neighborList.size(); i++){
        delete conf.neighborList[i].neighborID;
        conf.neighborList[i].neighborID = NULL;
    }


    /*if (conf.particle != NULL) // left just for peace of mind, before Particles* --> std::vector<Particles* >
        free(conf.particle);
    conf.particle = NULL;*/

    if (sim.clusterlist != NULL)
        free(sim.clusterlist);

    if (sim.clustersenergy != NULL)
        free(sim.clustersenergy);

    if(topo.switchlist)
        free(topo.switchlist);

    if (sim.pairlist_update) {
        if(deallocPairlist()) {
            return 1;
        }
    }

    for(int i=0; i<MAXMT; i++) {
        free(topo.chainparam[i].name);
    }
    return 0;
}




int MCSimSystem::deallocPairlist() {
    /*long i;
    if(sim.pairlist != NULL){
        for(i = 0; i < topo.npart; i++){
            if(sim.pairlist[i].pairs != NULL){
                free(sim.pairlist[i].pairs);
            }
        }
        free(sim.pairlist);
    }*/
    return 0;
}




