/** @file main.cpp*/

#include "mc/inicializer.h"
#include "mc/updater.h"

using namespace std;

Topo topo; // Global instance of topology

int main(int argc, char** argv) {

    FILE *outfile,*mov;       // Handle for writing configuration

    Sim sim;                  // Should contain the simulation options.
    Conf conf;                // Should contain fast changing particle and box(?) information
    FileNames files;

    Updater* updater; // need to get an instance of updater after initialization, because of initFCE

    sim.wl.setConf(&conf);

    cout << "\nPatchy Spherocylinders version 3.6\n"
         << "-------------------------------------" << endl;

#ifdef EXTRA_HYDROPHOBIC_ALL_BODY_ATTRACTION
    cout << "\n!!! Extra hydrophobic interaction in e_cpsc_cpsc added\n" << endl;
#endif


    /********************************************************/
    /*                  INITIALIZATION                      */
    /********************************************************/

    Inicializer init(&sim, &conf, &files);

    init.initWriteFiles();
    init.initMPI(argc,argv);
    init.readOptions();
    init.initTop(); // here particleStore filled in setParticleParams
    init.testChains(); // if no chains -> move probability of chains 0

    cout << "\nReading configuration...\n";
    if(init.poolConfig)
        init.initConfig(files.configurationPool, conf.pool);
    init.initConfig(files.configurationInFile, conf.pvec);

    cout << "Equilibration of maximum step sizes: " << sim.nequil/2 << " sweeps" << endl;

    if ( sim.wl.wlm[0] > 0 ) {
        outfile = fopen(files.wlinfile, "r");
        if (outfile == NULL) {
            printf ("ERROR: Cannot open file for Wang-Landau method (%s).\n",files.wlinfile);
            sim.~Sim(); topo.~Topo(); conf.~Conf();
            exit(1);
        }
        fclose (outfile);
    }

    // Empty movie file
    mov = fopen("movie", "w");
    fclose (mov);

    if (sim.pairlist_update) {
        init.initNeighborList();
        conf.pairlist_update = true;
    }

    updater = new Updater(&sim, &conf, &files);


    /********************************************************/
    /*                  EQUILIBRATION                       */
    /********************************************************/

    if (sim.nequil) {
        printf("\nStart equilibration...\n");

        updater->simulate(sim.nequil/2, sim.adjust, 0, 0);
        updater->simulate(sim.nequil/2, 0,          0, 0);

        sim.printEqStat();

        cout << "Further equilibration of configuration:  " << sim.nequil/2 <<  " sweeps" << endl;

        outfile = fopen("config.eq", "w");
        fprintf (outfile, "%15.8e %15.8e %15.8e\n", conf.box.x, conf.box.y, conf.box.z);
        conf.draw(outfile);
        fclose (outfile);

        printf ("   Equilibrated configuration written to config.eq\n");
        printf ("   Box dimensions: %.10f, %.10f, %.10f\n\n", conf.box.x, conf.box.y, conf.box.z);
    }

    /********************************************************/
    /*                  PRODUCTION RUN                      */
    /********************************************************/

    cout << "Production run:  "<< sim.nsweeps << " sweeps\n" << endl;

    updater->simulate(sim.nsweeps, 0, sim.paramfrq, sim.report);

#ifdef ENABLE_MPI
        printf ("   MPI replica changeT / changeP / acceptance ratio: \t %.6f   /   %.6f  /  %.6f\n\n", sim.mpiexch.mx,sim.mpiexch.angle,RATIO(sim.mpiexch));
#endif

    outfile = fopen(files.configurationoutfile, "w");
#ifdef TESTING
    fprintf (outfile, "%15.7e %15.7e %15.7e\n", conf.box.x, conf.box.y, conf.box.z);
#else
    fprintf (outfile, "%15.8e %15.8e %15.8e\n", conf.box.x, conf.box.y, conf.box.z);
#endif
    conf.draw (outfile);
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
            fprintf(outfile, "%s %d\n", topo.moleculeParam[i].name, conf.pvecGroupList.molCountOfType(i));

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

    /********************************************************/
    /*                   MEMORY DEALLOC                     */
    /********************************************************/

    //sim.~Sim(); topo.~Topo(); conf.~Conf(); // note: called automatically at end of main

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif

    printf ("\nDone\n\n");

    return 0;
}
