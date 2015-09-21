/** @file main.cpp*/

#include <fstream>

#include "mc/inicializer.h"
#include "mc/updater.h"
#include "mc/mygetline.h"
#include "mc/randomGenerator.h"
#include "unitTests/pvectester.h"

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

    void analyzeCur(double& r1, double& r2, double& fi, Conf* conf, int &mid, int &mid2);

int main(int argc, char** argv) {
#ifdef OMP1
    cout << "OPENMP SIMULATION" << endl;
#endif

#ifdef ENABLE_MPI
    cout << "MPI SIMULATION" << endl;
#endif

    FILE *infile,*outfile,*mov;       // Handle for writing configuration

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

    // count grand canonically active species
    for(int i=0; i<conf.pvec.molTypeCount; i++) {
        if(topo.moleculeParam[i].activity != -1.0)
            topo.gcSpecies++;
    }

    updater = new Updater(&sim, &conf, &files);

    /********************************************************/
    /*                  SOME TESTS                          */
    /********************************************************/

#ifndef NDEBUG
    //pVecTester pTest;
    //assert(pTest.test());
#endif

    /********************************************************/
    /*                      ANALYZE                         */
    /********************************************************/

    /*ofstream myfile;
    myfile.open ("curvature");

    double r1, r2, fi, a=0.0, b=0.0, c=0.0;
    double aver1=0.0, aver2=0.0;
    vector<double> array1;
    vector<double> array2;
    double s1=0.0, s2=0.0;
    int N = 0, mid=0, mid2=0;

    infile = fopen("movie", "r");
    if (infile == NULL) {
        fprintf (stderr, "\nERROR: Could not open %s file.\n\n", "movie");
        exit (1);
    }
    fseek ( infile , 0 , SEEK_SET );

    char * line;
    size_t line_size = (STRLEN + 1) * sizeof(char);
    for(int i=0; i<200; i++) {
        if(!init.initConfig(&infile, conf.pvec))
            break;
        analyzeCur(r1, r2, fi, &conf, mid, mid2);
        //if(r1 > r2-1.0 && r1 < r2+1.0) {
            cout << r1 << " " << r2 << " " << fi << endl;
            aver1 += r1;
            aver2 += r2;
            N++;
            array1.push_back(r1);
            array2.push_back(r2);
            cout << mid <<" "<< mid2 << endl;
        //}

        a += r1;
        b += r2;
        c += fi;
    }
    aver1 /= N;
    aver2 /= N;
    for(unsigned int w=0; w<array1.size(); w++) {
        s1 += (aver1 - array1[w])*(aver1 - array1[w]);
        s2 += (aver2 - array2[w])*(aver2 - array2[w]);

    }
    s1 /= N; s1 = sqrt(s1);
    s2 /= N; s2 = sqrt(s2);

    cout << "N= " << N << endl;
    cout << "r1="<< aver1 << ", s1= " << s1 << ", r2= " << aver2 << ", s2= " << s2 << endl;
    cout << "H= " << 1/aver1 - 1/aver2 << ", s= " << s1*s1/(aver1*aver1) + s2*s2/(aver2*aver2) << endl;
    cout << "K= " << 1/aver1 * 1/aver2 * (-1.0) << ", s= " << 1/(aver1*aver2) * (s1*s1/aver1 + s2*s2/aver2) << endl;

    myfile << (aver1+aver2)/2;
    myfile.close();
    fclose (infile);
    exit(1);*/

    /********************************************************/
    /*                  EQUILIBRATION                       */
    /********************************************************/

    assert((int)conf.pvec.size() == conf.pvec.first[conf.pvec.molTypeCount]);

    if (sim.nequil) {
        printf("\nStart equilibration...\n");

        updater->simulate(sim.nequil/2, sim.adjust, 0, 0);
        updater->simulate(sim.nequil/2, 0,          0, 0);

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

    updater->simulate(sim.nsweeps, 0, sim.paramfrq, sim.report);

#ifdef ENABLE_MPI
        printf ("   MPI replica changeT / changeP / acceptance ratio: \t %.6f   /   %.6f  /  %.6f\n\n", sim.mpiexch.mx,sim.mpiexch.angle,RATIO(sim.mpiexch));
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

    printf ("\nDone\n\n");

    return 0;
}


void analyzeCur(double &r1, double &r2, double &fi, Conf* conf, int& mid, int& mid2) {
    int result = 0;
    int size = 19;
    Vector dist, dist1, dist2;
    Vector offset;
    Vector min;
    Vector a,b;
    double minimum=999.9;
    int index_fi_max,index_fi_max2;
    double fi_max=0.0, r_min=999, fi_max2=0.0, r_min2=999;
    int index;
    double sinFi;
    double r;
    int indexBase = 1026, indexBase2=18;
    for(int i=size; i< size+size-1; i++) { // 0 - 495 with 1026
        result += i;

        dist = conf->pvec[result].pos;
        dist-= conf->pvec[indexBase].pos;

        offset = conf->pvec[indexBase].pos;
        dist.scale(0.5);
        offset += dist;
        dist.scale(2.0);

        minimum=999.9;
        for(unsigned int q=0; q<conf->pvec.size(); q++) {
            min = offset;
            min -= conf->pvec[q].pos;
            if(sqrt(min.dot(min)) < minimum) {
                index = q;
                minimum = sqrt(min.dot(min));
            }
        }
        a = conf->pvec[result].pos - conf->pvec[index].pos;
        b = conf->pvec[indexBase].pos - conf->pvec[index].pos;



        sinFi = (a.cross(b)).size() / (a.size() * b.size());

        dist *= 30.0;
        r = dist.size() / (2*sinFi);

        if(r < r_min) {
            dist1 = dist;
            r_min = r;
            fi_max=sinFi;
            index_fi_max = result;
            mid = index;
        }

    }

    //cout << index_fi_max << " " << indexBase << endl;
    //cout << "dist=" <<dist1.size() << " r=" <<r_min << " Fi: " << 180- asin(fi_max)*57.2957795 << "\n" << endl;

    for(int i=size+size-1; i> size; i--) { // 495 - 1008 with 18
        result += i;

        dist = conf->pvec[result].pos;
        dist-= conf->pvec[indexBase2].pos;

        offset = conf->pvec[indexBase2].pos;
        dist.scale(0.5);
        offset += dist;
        dist.scale(2.0);

        minimum=999.9;
        for(unsigned int q=0; q<conf->pvec.size(); q++) {
            min = offset;
            min -= conf->pvec[q].pos;
            if(sqrt(min.dot(min)) < minimum) {
                index = q;
                minimum = sqrt(min.dot(min));
            }
        }
        a = conf->pvec[result].pos - conf->pvec[index].pos;
        b = conf->pvec[indexBase2].pos - conf->pvec[index].pos;



        sinFi = (a.cross(b)).size() / (a.size() * b.size());

        dist *= 30.0;
        r = dist.size() / (2*sinFi);

        if(r < r_min2) {
            dist2 = dist;
            r_min2 = r;
            fi_max2=sinFi;
            index_fi_max2 = result;
            mid2 = index;
        }
    }



    //cout << index_fi_max2 << " " << indexBase2 << endl;
    //cout << "dist=" << dist2.size() << " r=" <<r_min2 << " Fi: " << 180- asin(fi_max2)*57.2957795 << "\n" << endl;

    //cout << "angle dist1, dist2 = " << asin((dist1.cross(dist2)).size() / (dist1.size() * dist2.size()))*57.2957795 << endl;

    r1 = r_min;
    r2 = r_min2;
    fi = asin((dist1.cross(dist2)).size() / (dist1.size() * dist2.size()))*57.2957795;
}


