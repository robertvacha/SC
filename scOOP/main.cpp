/** @file main.cpp*/

#include <fstream>
#include <ctime>
#include <iomanip>

#include "mc/inicializer.h"
#include "mc/updater.h"
#include "mc/mygetline.h"
#include "mc/randomGenerator.h"
#include "unitTests/pvectester.h"

#include "mc/analysis.h"

using namespace std;

Topo topo; // Global instance of topology

#ifdef RAN2
    Ran2 ran2;
#else
  #ifdef DSFMT
    Dsfmt ran2;
  #else
    MersenneTwister ran2;
  #endif
#endif

int main(int argc, char** argv) {
    int rank=0, procs=1;

#ifdef ENABLE_MPI
    cout << "MPI SIMULATION" << endl;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procs );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank );

#endif

#ifdef EXTRA_HYDROPHOBIC_ALL_BODY_ATTRACTION
    cout << "\n!!! Extra hydrophobic interaction in e_cpsc_cpsc added, strenght: " << E_ISO << endl;
#endif

    cout << "\nPatchy Spherocylinders version 3.6\n-------------------------------------" << endl;

    FILE *infile,*outfile,*mov;       // Handle for writing configuration

    FileNames files(rank);
    Conf conf;                // Should contain fast changing particle and box(?) information
    Sim sim(&conf, &files, rank, procs);                  // Should contain the simulation options.

    /********************************************************/
    /*                  INITIALIZATION                      */
    /********************************************************/

    Inicializer init(&sim, &conf, &files);
    init.initTop(); // here particleStore filled in setParticleParams
    init.testChains(); // if no chains -> move probability of chains 0

    cout << "\nReading configuration...\n";
    if(init.poolConfig) {
        infile = fopen(files.configurationPool, "r");
        if (infile == NULL) {
            fprintf (stderr, "\nERROR: Could not open %s file.\n\n", files.configurationPool);
            exit (1);
        }
        if(!init.initConfig(&infile, conf.pool))
            exit(1);
        fclose (infile);
    }
    infile = fopen(files.configurationInFile, "r");
    if (infile == NULL) {
        fprintf (stderr, "\nERROR: Could not open %s file.\n\n", files.configurationInFile);
        exit (1);
    }
    if(!init.initConfig(&infile, conf.pvec))
        exit(1);
    conf.geo.info();
    fclose (infile);

    cout << "Equilibration of maximum step sizes: " << sim.nequil/2 << " sweeps" << endl;

    // Empty movie file
    mov = fopen("movie", "w");
    fclose (mov);

    if (sim.pairlist_update) {
        init.initNeighborList();
        conf.pairlist_update = true;
    }

    // count grand canonically active species
    for(int i=0; i<conf.pvec.molTypeCount; i++) {
        if(topo.moleculeParam[i].activity != -1.0)
            topo.gcSpecies++;
    }

    conf.initEMatrix();

    //
    // THESE REQUIRES SIM AND TOP DATA
    //

    //
    //  PRE-COMPUTE PAIRLIST CUTOFF, After topo init
    //
    for(int i=0; i< MAXT; i++) {
        for(int j=0; j<MAXT; j++) {
            sim.max_dist_squared[i][j] = AVER(sim.stat.trans[i].mx, sim.stat.trans[j].mx);
            sim.max_dist_squared[i][j] *= (1 + sim.pairlist_update) * 2;
            sim.max_dist_squared[i][j] += topo.maxcut;
            sim.max_dist_squared[i][j] *= sim.max_dist_squared[i][j]; /* squared */
            if(sim.cell < sim.max_dist_squared[i][j])
                sim.cell = sim.max_dist_squared[i][j];
        }
    }

    if(sim.nGrandCanon == 0) {
        bool test = false;
        for(int i=0; i < MAXMT; ++i) {
            if( topo.moleculeParam[i].name != NULL && topo.moleculeParam[i].activity != -1 ) {
                test = true; // we have a gc active species
                break;
            }
        }
        if(test) {
            cout << "Activity stated in top.init, But nGrandCanon=0 in options" << endl;
            exit(1);
        }
    }

    Updater updater(&sim, &conf, &files); // need to get an instance of updater after initialization, because of initFCE

    /********************************************************/
    /*                      ANALYZE                         */
    /********************************************************/

    /*ofstream user;
    user.open ("user.dat");

    ofstream myfile;
    myfile.open ("curvatureFull");

    ofstream myfileS;
    myfileS.open ("curvatureSelect");

    ofstream myfileF;
    myfileF.open ("curvature");

    //
    // r1, r2 - principal curvatures
    // fi - angle between principal planes
    //
    double r1, r2, fi;
    int R1i, R1j, R2i, R2j, center;
    int frames = 0;

    infile = fopen("movieX", "r");
    if (infile == NULL) {
        fprintf (stderr, "\nERROR: Could not open %s file.\n\n", "movie");
        exit (1);
    }
    fseek ( infile , 0 , SEEK_SET );

    char * line;
    long pos;
    size_t line_size = (STRLEN + 1) * sizeof(char);
    int i=0, ii=0, prev = 0, offset=0;

    cout << "Start:" << endl;

    pos = ftell(infile);
    if(myGetLine(&line, &line_size, infile) == -1){ // number of part
        fprintf (stderr, "ERROR: Could not read box size (Inicializer::initConfig)\n\n");
        return false;
    }
    if(myGetLine(&line, &line_size, infile) == -1){ // frame info
        fprintf (stderr, "ERROR: Could not read box size (Inicializer::initConfig)\n\n");
        return false;
    }
    strip_comment(line);
    trim(line);

    sscanf(line, "%*s %d", &i);
    fseek(infile, pos, SEEK_SET);
    double e[1027];
    while (init.initConfig(&infile, conf.pvec) ) {

        for(unsigned int j=0; j< conf.pvec.size(); j++) {
            e[j] += updater->calcEnergy.oneToAllBasic(j);
            //user << e[j] << " "<< e[j] << " "<< e[j] << " "<< e[j] << " "<< e[j] << " "<< e[j] << " ";
        }
        //user << "\n";
        frames++;
        cout << updater->calcEnergy.allToAllBasic() << " " << i << endl;

        /*analyzeCur(r1, r2, fi, R1i, R1j, R2i, R2j, center, &conf, &sim);
        cout << std::setw(8) << std::setprecision(2) << std::fixed;
        cout << i << " "<< std::setw(7) <<  r1 << " " << std::setw(7)<< r2 << " " << std::setw(5) <<fi
             << " " << std::setw(4)<< R1i << " " << std::setw(4)<< R1j << " " << std::setw(4)<< R2i << " " << std::setw(4)<< R2j  << " " << center<< endl;

        myfileF << std::setw(8) << std::setprecision(2) << std::fixed;
        myfileF << i << " "<< std::setw(7) <<  r1 << " " << std::setw(7)<< r2 << " " << std::setw(5) <<fi
             << " " << std::setw(4)<< R1i << " " << std::setw(4)<< R1j << " " << std::setw(4)<< R2i << " " << std::setw(4)<< R2j  << " " << center<< endl;

        if(R1i != -1 && R1j != -1 && R2i != -1 && R2j != -1) {
            myfile << i << " " << r1 << " " << r2 << " " << fi << endl;
            if(fi > 88) {
                myfileS << i << " " << r1 << " " << r2 << " " << fi << endl;
            }
        }*//*

        pos = ftell(infile);
        if(myGetLine(&line, &line_size, infile) == -1){ // number of part
            fprintf (stderr, "ERROR: Could not read box size (Inicializer::initConfig)\n\n");
            break;
        }
        if(myGetLine(&line, &line_size, infile) == -1){ // frame info
            fprintf (stderr, "ERROR: Could not read box size (Inicializer::initConfig)\n\n");
            break;
        }
        strip_comment(line);
        trim(line);

        prev = ii;
        sscanf(line, "%*s %d", &ii);
        i=ii;
        if(prev > ii) {
            offset += prev;
        }
        i += offset;
        fseek(infile, pos, SEEK_SET);


    }
    cout << frames << endl;
    for(unsigned int j=0; j< 1027; j++) {
        e[j] /= frames;
    }
    for(int k=0; k<frames; k++) {
        for(unsigned int j=0; j< 1027; j++) {
            user << e[j] << " "<< e[j] << " "<< e[j] << " "<< e[j] << " "<< e[j] << " "<< e[j] << " ";
        }
        user << "\n";
    }

    user.close();
    myfile.close();
    fclose (infile);
    cout << "Exit" << endl;
    exit(1);*/

    /********************************************************/
    /*                  EQUILIBRATION                       */
    /********************************************************/

    assert((int)conf.pvec.size() == conf.pvec.first[conf.pvec.molTypeCount]);

    if (sim.nequil) {
        printf("\nStart equilibration...\n");

        updater.simulate(sim.nequil/2, sim.adjust, 0, 0);
        updater.simulate(sim.nequil/2, 0,          0, 0);

        sim.printEqStat();

        cout << "Further equilibration of configuration:  " << sim.nequil/2 <<  " sweeps" << endl;

        outfile = fopen("config.eq", "w");
        fprintf (outfile, "%15.8e %15.8e %15.8e\n", conf.geo.box.x, conf.geo.box.y, conf.geo.box.z);
        conf.draw(outfile);
        fclose (outfile);

        printf ("   Equilibrated configuration written to config.eq\n");
        printf ("   Box dimensions: %.10f, %.10f, %.10f\n\n", conf.geo.box.x, conf.geo.box.y, conf.geo.box.z);
    }

    /********************************************************/
    /*                  PRODUCTION RUN                      */
    /********************************************************/

    cout << "Production run:  "<< sim.nsweeps << " sweeps\n" << endl;

    sim.all = clock();
    updater.simulate(sim.nsweeps, 0, sim.paramfrq, sim.report);
    sim.all = clock() - sim.all;

#ifdef ENABLE_MPI
    cout << sim.pseudoRank << "p, MPI replica changeT / changeP / acceptance ratio: " << sim.stat.mpiexch.mx << ", " << sim.stat.mpiexch.angle << ", " << sim.stat.mpiexch.ratio() << endl;
    //cout << sim.pseudoRank << "p, rej: " << sim.mpiexch.rej << ", acc: " << sim.mpiexch.acc << endl;
#endif

    outfile = fopen(files.configurationoutfile, "w");
#ifdef TESTING
    fprintf (outfile, "%15.6e %15.6e %15.6e\n", conf.geo.box.x, conf.geo.box.y, conf.geo.box.z);
#else
    fprintf (outfile, "%15.8e %15.8e %15.8e\n", conf.geo.box.x, conf.geo.box.y, conf.geo.box.z);
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
        for(int i=0; i < conf.pvec.molTypeCount; i++)
            fprintf(outfile, "%s %d\n", topo.moleculeParam[i].name, conf.pvec.molCountOfType(i));

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

    cout << "\nSimulation time:" << (double)sim.all / CLOCKS_PER_SEC << " s" << endl;
    cout << "PairList generation time:" << 100.0*(double)sim.pairList / sim.all << " %" << endl;

    cout << "\nDONE" << endl;

    return 0;
}


