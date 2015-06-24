/** @file main.cpp*/

#include "mc/inicializer.h"
#include "mc/updater.h"
#include "mc/mygetline.h"

using namespace std;

Topo topo; // Global instance of topology

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
        init.initConfig(&infile, conf.pool);
        fclose (infile);
    }
    infile = fopen(files.configurationInFile, "r");
    if (infile == NULL) {
        fprintf (stderr, "\nERROR: Could not open %s file.\n\n", files.configurationInFile);
        exit (1);
    }
    init.initConfig(&infile, conf.pvec);
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

    updater = new Updater(&sim, &conf, &files);

    /********************************************************/
    /*                  ANALYZE                             */
    /********************************************************/

    /*double r1, r2, fi, a=0.0, b=0.0, c=0.0;
    double aver1=0.0, aver2=0.0;
    vector<double> array1;
    vector<double> array2;
    double s1=0.0, s2=0.0;
    int N = 0, mid=0, mid2=0;

    infile = fopen("movieAll4", "r");
    if (infile == NULL) {
        fprintf (stderr, "\nERROR: Could not open %s file.\n\n", "movieAll");
        exit (1);
    }
    char * line;
    size_t line_size = (STRLEN + 1) * sizeof(char);
    for(int i=0; i<200; i++) {
        init.initConfig(&infile, conf.pvec);
        analyzeCur(r1, r2, fi, &conf, mid, mid2);
        if(r1 > r2-0.1 && r1 < r2+0.1) {
            cout << r1 << " " << r2 << " " << fi << endl;
            aver1 += r1;
            aver2 += r2;
            N++;
            array1.push_back(r1);
            array2.push_back(r2);
            cout << mid <<" "<< mid2 << endl;
        }

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
    fclose (infile);
    exit(1);*/

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

/*! \mainpage My Personal Index Page
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 *
 \subsection todo TODO in future s
  - Non-equilibrium candidate moves
  - check scaling of particles of different sizes - should scale with contact area!
  - cell list - divide simulation box in cells where
  particles interact with each other and outside is definitely 0 - safe time
  better scaling with system size, possibly long spherovylinders could be in several
  celles to keep good scaling
  - better cluster algorithm - put in wang-landau
  - cluster list work for spherocylinders only now

\subsection verOOP scOOP 1.0
- Grandcanonical ensemble added MoveCreator::muVTMove()
  - insertion and deletion moves
  - example topologies:

\subsection ver35 Version 3.5
   - linear bond at spherocylinders, where second spherocilinder is harmonicaly
     attached to a point that is in distance of bondlength from the first spherocylinder
     and it follows the direction of spherocylinder
   - bonded particles belong to the same cluster
   - print energy at statistical reports
   - have particles of different lengths
   - interaction scaling back to v1+v2 (no addition of 1.0) - more physical

\subsection ver34 Version 3.4
    - New handling of the option file
    - reaction coordinate radius around z axis for a pore calculations
    - reaction coordinate as number of particles in contact (defined by distance of CMs)
    - 2D Wang-Landau method
    - New Wang-Landau coordinate - radius pore in vesicle around begining of xy plane
    - New models TPSC, TCPSC, TCHPSC, TCHCPSC- models with two patches
      note that switch function on sides of patch are linear in cos angle not in angle
      as a results two patches with overlaping sides do not compensate easily to a flat profile
    - FIX chirality was doubled (angle twice as large)
    - Added posibility of exluded interactions [EXCLUDE] in topology file
    - MPI replica exchange with different temperatures and pressure (paraltemp paralpress)
      input configuration is #{number of process}config.init, if it does not exist config.init is used
      each replica is with different random seed = seed+mpirank
    - config.init can look like movie snapshot
    - MPI exchange with Wang-Landau
    - added angular interaction between neighboring spherocylinders (in chain)
      angle1 is angle between sc directions and angle2 ins angle between the patches

\subsection ver33  Version 3.3
  - external potantial can be added as a part of topology (it can be hard or attractive wall)


\subsection ver32 Changes made by Noah S. Bieler and Robert Vacha: New version 3.2

 * - The length has now to be specified in the topology file, but they are not
 * allowed to differ from each other. The option file shall no longer contain
 * a length option.
 * - The particles can now switch their type based on the chemical potential
 * delta_mu (= energy difference from state 2 to state 1).
 * - For that a new option was introduced: Average attempts per sweep to switch
 * a type.
 * - A lot of variables are now combined in either topo, sim or conf. The rule
 *   should be:
 *   - topo: Everything that belongs to the topology and that should not change
 *           during the game.
 *   - sim:  Options and stuff, that has to do with the simulation. (Maybe the
 *           current target and so should be saved in there as well)
 *   - conf: What changes every step concerning the particles and the box or
 *           in other words: what has been read from conf.init
 * - added a cluster determing routine => sim->clusterlist + sim->clusters
 * - added macros for TRUE and FALSE
 * - Added Option for the random seed
 * - Basic Neighbour list implemented
 * - New types: chiral CPSC (CHCPSC) and chiral PSC (CHPSC) and their interactions


 \subsection ver31 Patchy Spherocylinder Version 3.1

  - Wang-Landau method of free energy calculations
    - It is set in options file as:
      - O = none,
      - 1 = z-distance of 1st paticle from system CM
      - 2 = hole in xyplane of SCA = membrane hole
  - It reads a file wl.dat and write wl-new at the end. There is value of alpha at the first line and then
  there are three columns:
    - 1- order parameter,
    - 2- weights,
    - 3- histogram
  - Interaction of spherocylinders is scaled based on the volume of attractive patch, the unit of one
  is that two spheres of diameter sigma =1.0 are attracting each other by 1.0. Using this in interaction
  among lipids and spherocylinders should be consistent.
  - Start up configuration "config.init" file has a box size at the first line now.
  (I tested performance: compilation with optimization -O2 speed up 10%
  rest has negligible effect including usage of static arrays instead of dynamic
  most of the time consumes paire function.
  6,519,638,177  :simulate
  6,492,411,300  :energyone
  5,705,685,593  :paire
  542,561,887  :bondenergy
  489,463,361  :eattractive11
  450,443,970  :image
  115,126,519  :erepulsive



 \subsection ver30  Patchy Spherocylinder Version 3.0

 - Beads were added to the particle list.
 - bead(10) - repulsive
 - bead(11) - isotropocally attractive
 - It is necessary to provide also a topology file (top.init)
 - Particles are placed in chains according to the topology order including connections
 - Particle arryas are allocated dynamicly on heap now
 - dispacement and rotation are optimized for highest RMSD performace
 - NPT ensemble with isotropic and anisotropic couplings, in pressure moves all
   particles are rescaled with their center (chains are not rescaled with CM)
    -0 - anisotropic coupling,
    -1 - isotropic coupling,
    -2 - isotropic in xy z=const
 - bead types and their interactions
 - repulsive(10) purely repulsive shpere with WCA potential on closest distance
 - parameters: Patch repulsion sigma - defined where repulsion reaches zero
   isotropic(11) - isotropic cos^2 potential is acting isotropicaly dependent only on
   closest distance between obejcts.
 - Parameters: distance of attractivity (should be at least
   sigma*2^(1/6)) defines how far is attraction constant -e. After this distance
  follows switch length on which attraction goes to zero as cos^2.
 - Rest as repulsive model.

\subsection ver2 Patchy Spherocylinder Version 2.0

 - It is possible to make chains of spherocylinders that are connected through
hemispherical caps by harmonic bond.
 - There are two parameters eq distance and
strength of harmonic spring, note that units are in 1 kT/e, the MC strength of bond
is changing with parameter temperature..

\subsection ver1 Patchy Spherocylinder Version 1

 - Includes diffferent types of possible interactions:
   - repulsive(0) - purely repulsive spherocylinder with WCA potential on closest distance.
     - parameters: Patch repulsion sigma - defined where repulsion reaches zero.
   - isotropic(1) - isotropic cos^2 potential is acting isotropicaly dependent only on
     closest distance between spherocylinders.
     - Parameters: distance of patch, Interaction distance of patch (should be at least
       sigma*2^(1/6)) defines how far is attraction constant -e. After this distance
       follows Switch length on which attraction goes to zero as cos^2.
   - Rest as repulsive model.

 - patchy(2)
   - Attractive potential in limited to an angular wedge on spherocylinder.
   - Patch goes all the way through, making also hemispherical caps on end attractive.
     - Parameters:Anglular part has a parameter defining it size "Angular size of patch
       (degrees)" and witdh of switch function "Angular switch off of patch (degrees)" on which
       attraction reaches zero
     - it is a linear function
     - Rest as isotropic model.

 - cylindrical(3)
   - Attractive potential in limited to an angular wedge on cylindrical part
     of spherocylinders. The hemispherical caps on ends are repulsive.
   - Rest as patchy model.

Note particles are inside numbered from 0, there is prealocated size of particles MAXN
because in future there can be grand canonical ensamble and number of particles may vary

Follows mc of hard wall spherocylinder version 7 by Mark Miller -description below


\subsection ver0 Version 1

Performs basic constant volume MC simulation of hard spherocylinders with rigid
cuboidal boundary conditions.

Run parameters are read in from the file "options".  The template for this file
appears at the end of the code.  The values must be inserted before the colons.

The initial configuration is read from the file "config.init".  The first line contain size
of box The format for the file is nine columns: three for the positions and three for the
direction vector and three for direction of pathc.  The direction vectors are normalised
after being read in.  The configuration is checked for particle overlaps.

The unit of length is taken as the spherocylinder diameter.  Hence the ratio
L/D is equal to the length of the cylinder.

Order parameters for nematic and smectic order are evaluated.  The nematic order
parameter is related to the coefficient of the quadratic term in the Legendre
expansion of the orientational distribution function.  Any smectic order is
assumed to be directed along the z axis, and is detected by the coefficients
of the Fourier expansion of the position distribution function.

MM 12.vii.01

..................................................................................

Version 2

The aspect ratio of the box may now fluctuate, keeping the volume constant.
Two new parameters are required in the options file to specify the average number
of attempted shape changes per sweep, and the initial maximum trial change in
a box dimension.

Shape changes are made by picking one of the three box lengths at random,
making a random change, evenly distributed between plus and minus a finite
interval, choosing a second direction and doing the same, then determining
the new length in the remaining direction from the condition of constant
volume.

The step-size equilibration period is now split into three parts: displacement,
rotation, and shape change.

The most important change to the code is that the particle coordinates are
now stored as fractions of the box dimensions.  However, input and output
configurations are still communicated in units of the cylinder diameter, D=1.

Note that the displacement maximum step size is now specified as a fraction of
the box length, not as an absolute distance.

MM 18.vii.01

..................................................................................

Version 3

Constant pressure MC.  The volume may fluctuate.  Volume changes are attempted
by altering just one box length at a time, chosen at random.  The running
average of the density is calculated and reported.

MM 24.vii.01

..................................................................................

Version 7

The composite translation-plus-rotation moves have been split into separate
move types, each of which is attempted with equal probability.  This enables
acceptance ratios to be accumulated separately for these degrees of freedom, so
that maximum step sizes can be adjusted more sensibly.

A few other things have been tidied up, such as defining structures for the
book-keeping of statistics and acceptance ratios.

MM 9.v.02

--------------------------------------------------------------------------------
 */
