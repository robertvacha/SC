/*MODIFICATIONS:
when wl=3 .. rotation of 0th particle it is allowed to move

pokousime se odstranit bug ze castice pri wl=1 se obcas zasekne
je to kvuli pohybu center of mass

zrejme je chybou tim ze kvuli akumulaci numerickych chyb se cas od casu prepocita centrum hmotnosti systemu
u toho se pouziji periodicke okrajove podminky, ktere to cele mohou zkazit protoze castice co odletela do druheho boxu
pri navratu od primarniho boxu (-z/2 do z/2) posune stred hmostnosti

testujeme - odstranili jsme prepocet stredu hmotnosti z obav z numerickych chyb


vypis "F I N I S H E D" do stdout po splneni wl podminek

*/


/*TODO in future s
  - Non-equilibrium candidate moves
  - check scaling of particles of different sizes - should scale with contact area!
  - cell list - divide simulation box in cells where
  particles interact with each other and outside is definitely 0 - safe time
  better scaling with system size, possibly long spherovylinders could be in several 
  celles to keep good scaling
  - better cluster algorithm - put in wang-landau 
  - cluster list work for spherocylinders only now
 */

/*------------------------------------------------------------------------------
 Version 3.5
   - linear bond at spherocylinders, where second spherocilinder is harmonicaly 
     attached to a point that is in distance of bondlength from the first spherocylinder 
     and it follows the direction of spherocylinder
   - bonded particles belong to the same cluster
   - print energy at statistical reports
   - have particles of different lengths
   - interaction scaling back to v1+v2 (no addition of 1.0) - more physical
 
*/ 
/*------------------------------------------------------------------------------
    Version 3.4
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
*/
/*-------------------------------------------------------------------------------
  Version 3.3
  -external potantial can be added as a part of topology - it can be hard or attractive wall
 */

/** 
 * Changes made by Noah S. Bieler and Robert Vacha:
 *
 * New version 3.2
 *
 * - The length has now to be specified in the topology file, but they are not 
 * allowed to differ from each other. The option file shall no longer contain
 * a length option.
 * - The particles can now switch their type based on the chemical potential
 * delta_mu (= energy difference from state 2 to state 1).
 * - For that a new option was introduced: Average attempts per sweep to switch
 * a type.
 * - A lot of variables are now combined in either topo, sim or conf. The rule 
 *   should be:
 *   > topo: Everything that belongs to the topology and that should not change
 *           during the game.
 *   > sim:  Options and stuff, that has to do with the simulation. (Maybe the 
 *           current target and so should be saved in there as well)
 *   > conf: What changes every step concerning the particles and the box or
 *           in other words: what has been read from conf.init
 * - added a cluster determing routine => sim->clusterlist + sim->clusters
 * - added macros for TRUE and FALSE
 * - Added Option for the random seed
 * - Basic Neighbour list implemented
 * - New types: chiral CPSC (CHCPSC) and chiral PSC (CHPSC) and their interactions
 */

/*--------------------------------------------------------------------------------
  sc31.c
  Patchy Spherocylinder Version 3.1
  Wang-Landau method of free energy calculations
  It is set in options file as:
  O = none, 1 = z-distance of 1st paticle from system CM, 2 = hole in xyplane of SCA = membrane hole 
  It reads a file wl.dat and write wl-new at the end. There is value of alpha at the first line and then
  there are three columns:
  1- order parameter, 2- weights, 3- histogram
  Interaction of spherocylinders is scaled based on the volume of attractive patch, the unit of one
  is that two spheres of diameter sigma =1.0 are attracting each other by 1.0. Using this in interaction
  among lipids and spherocylinders should be consistent.
  Start up configuration "config.init" file has a box size at the first line now.
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
 */
/*  --------------------------------------------------------------------------------
    sc3.c
    Patchy Spherocylinder Version 3.0

Beads were added to the particle list.
bead(10) - repulsive
bead(11) - isotropocally attractive
-It is necessary to provide also a topology file (top.init)
-Particles are placed in chains according to the topology order including connections
-Particle arryas are allocated dynamicly on heap now
-dispacement and rotation are optimized for highest RMSD performace
-NPT ensemble with isotropic and anisotropic couplings, in pressure moves all 
particles are rescaled with their center (chains are not rescaled with CM)
0 - anisotropic coupling, 1 - isotropic coupling, 2 - isotropic in xy z=const
bead types and their interactions 
repulsive(10) purely repulsive shpere with WCA potential on closest distance
parameters: Patch repulsion sigma - defined where repulsion reaches zero
isotropic(11) - isotropic cos^2 potential is acting isotropicaly dependent only on 
closest distance between obejcts. 
Parameters: distance of attractivity (should be at least 
sigma*2^(1/6)) defines how far is attraction constant -e. After this distance
follows switch length on which attraction goes to zero as cos^2. 
Rest as repulsive model.

sc2.c
Patchy Spherocylinder Version 2.0

It is possible to make chains of spherocylinders that are connected through
hemispherical caps by harmonic bond. There are two parameters eq distance and 
strength of harmonic spring, note that units are in 1 kT/e, the MC strength of bond
is changing with parameter temperature..

Patchy Spherocylinder Version 1.0

Includes diffferent types of possible interactions:
repulsive(0) - purely repulsive spherocylinder with WCA potential on closest distance. 
parameters: Patch repulsion sigma - defined where repulsion reaches zero.
isotropic(1) - isotropic cos^2 potential is acting isotropicaly dependent only on 
closest distance between spherocylinders. 
Parameters: distance of patch, Interaction distance of patch (should be at least 
sigma*2^(1/6)) defines how far is attraction constant -e. After this distance
follows Switch length on which attraction goes to zero as cos^2. Rest as repulsive model.
patchy(2) - Attractive potential in limited to an angular wedge on spherocylinder. Patch
goes all the way through, making also hemispherical caps on end attractive.
Parameters:Anglular part has a parameter defining it size "Angular size of patch 
(degrees)" and witdh of switch function "Angular switch off of patch (degrees)" on which
attraction reaches zero - it is a linear function. Rest as isotropic model.
cylindrical(3) - Attractive potential in limited to an angular wedge on cylindrical part
of spherocylinders. The hemispherical caps on ends are repulsive. Rest as 
patchy model.

Note particles are inside numbered from 0, there is prealocated size of particles MAXN
because in future there can be grand canonical ensamble and number of particles may vary

Follows mc of hard wall spherocylinder version 7 by Mark Miller -description below

sc.c

Version 1

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

--------------------------------------------------------------------------------*/

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef MACOS
# include "getline.h"
#endif

#ifdef MPI
# include <mpi.h>
#endif

      /* Macros for DEBUG messages */
#ifdef DEBUGGING_INIT
#define DEBUG_INIT(...) fprintf(stderr, "DB in INIT: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); fflush(stderr);
#else 
#define DEBUG_INIT(...)
#endif

#ifdef DEBUGGING_SIM
#define DEBUG_SIM(...) fprintf(stderr, "DB in SIM: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); fflush(stderr);
#else 
#define DEBUG_SIM(...)
#endif

#ifdef DEBUGGING
#define DEBUG(...) fprintf(stderr, "DB: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); fflush(stderr);
#else 
#define DEBUG(...)
#endif
      /* End of DEBUG macros */

      /* With  pairlist ? */
#define WITH_PAIRLIST

      /* Boolean Macros */
#define BOOL int
#define TRUE 1
#define FALSE 0
      /* End of Boolean Macros */

#define MAXF 20             /* Maximum number of Fourier terms */
#define MAXN 14000           /* Maximum number of particles */
#define MAXCHL 10          /* Maximum length of chain */
#define ZEROTOL 1.0e-12     /* Dot products below ZEROTOL are deemed zero */
#define ZEROTOL2 1.0e-8     /* numbers below ZEROTOL are deemed zero */
#define PI 3.141592653589793238462643383279  /* pi */
#define PIH 1.57079632679489661923132169163975  /* pi half*/

      /*Particle types*/
#define SC 10             /*spherocylinder*/
#define SCN SC+0          /*spherocylinder non-attractive*/
#define SCA SC+1          /*spherocylinder isotropicaly attractive*/
#define PSC SC+2          /*spherocylinder with patchy attraction*/
#define CPSC SC+3         /*spherocylinder with cylindrical patchy attraction*/
#define CHPSC SC+4        /* chiral psc */
#define CHCPSC SC+5       /* chiral cpsc */
#define TPSC SC+6          /*spherocylinder with two patches*/
#define TCPSC SC+7         /*spherocylinder with two cylindrical patches*/
#define TCHPSC SC+8        /* chiral 2psc */
#define TCHCPSC SC+9       /* chiral 2cpsc */

#define SP 30             /*sphere - should be over all apherocylinders*/
#define SPN SP+0          /* sphere non-attractive*/
#define SPA SP+1          /* spherocylinder isotropicaly attractive*/

#define MAXT 30           /* Maximum number of types we have*/
#define MAXMT 100         /* Maximum number of molecular types */

      /*Reading topology*/
#define SMSTR 64           /* Small string length */
#define STRLEN 400         /* maximum length of line*/
#define CONTINUE    '\\'   /* symbol for line continue*/
#define COMMENTSIGN '#'    /* symbol for comment*/
#define OPENKEY  '['       /* starting sign for keyword*/ 
#define CLOSEKEY ']'       /* ending sign for keyword*/ 
#define SEPARATOR ':'      /* sign for separator*/ 
#define OPENMOL '{'        /* starting sign for molecules*/
#define CLOSEMOL '}'       /* ending sign for molecules*/
#define BOXSEP 'x'       /* extraction of box*/

      /* Wang Landau method */
#define WL_GERR 0.0001           /* Max roughnes in histogram */
#define WL_ALPHATOL 0.0000001     /* Covergence crietria for detailed balance */
#define WL_MINHIST 1000         /* Minimum histogram sampling for considering roughness */
#define WL_ZERO 0.000000000000  /* Zero for histogram with current weights*/  
#define WL_CONTACTS 36.0          /* Square distance under which are particles in contact */ 

      /* Math */
#define DOT(a,b) ((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)   /* Dot product */
#define AVER(a,b) ((a+b)*0.5)                                      /* Arithmetic average*/
#define ROUND(a) (a > 0.0) ? floor(a + 0.5) : ceil(a - 0.5);       /* Round double*/
#define PMONE(a) (1 - 2 * a)                                       /* Takes 1 or 0, return +-1 */

      /* Acceptance ratio */
#define RATIO(a) ( ((a).acc+(a).rej) > 0 ? 1.0*(a).acc/((a).acc+(a).rej) : 0.0 )

#define INBOX(a,b) ( a > 0 ? modf(a,&b) : modf(a,&b)+1 )


/*................................................................
 Structure definitions
*/

struct vector {             /* Define a 3D vector structure */
        double x;
        double y;
        double z;
};

struct quat {             /* Define a quaternion structure */
	double w;
	double x;
	double y;
	double z;
};

struct particles {     /* Define a particle */
	struct vector pos;       /* Position vector */
	struct vector dir;       /* Unit direction vector of axis */
	struct vector patchdir[2];  /* Vector defining orientation of patch */
	struct vector patchsides[4];  /* Vector defining sides of patch */
	struct vector chdir[2];     /* Direction for chirality - keep in memory to increase speed */
	long chaint;             /* Chain type*/
	long chainn;             /* Chain number*/
	int type;                /* Type of the particle */ 
	int switchtype;          /* With which kind of particle do you want to switch?*/
	double delta_mu;         /* Chemical potential for the switch */
	int switched;            /* 0: in initial stat; 1: in the switched stat */
};

struct ia_param{      /* Contatins properties and parameters of particle types */
	char name[SMSTR];           /* The name of the particle type */
	char other_name[SMSTR];     /* The name of the particle type */
	int geotype[2];             /* The geometrical type: spherocylinder (0-repulsive, 1-isotropic, 2-patchy, 3-cylindrical) 
				       or sphere (0-repulsive, 1-isotropic) */
	double sigma;               /* Repulsion wca*/
	double epsilon;             /* Repulsion strength*/
	double pdis;                /* Interaction distance of patch */
	double pswitch;             /* Switch of distance of patch */
	double pangl[4];            /* angular size of patch as was specifid in input */
	double panglsw[4];          /* angular size of patchswitch as was specifid in input */
	double pcangl[4];           /* cosine of half size angle - rotation from patch direction to side */
	double pcanglsw[4];         /* cosine of half size angle plus switch - rotation from patch direction to side */
	double rcut;                /* Cutoff for attraction */
	double rcutwca;             /* Cutoff for repulsion*/
	double pcoshalfi[4];        /* Cosine of half angle going to side of interaction */
	double psinhalfi[4];        /* Sine of half angle going to side of interaction -useful for quaterion rotation */
	double csecpatchrot[2];      /* Cosine of Rotation of second patches in 2psc models*/
	double ssecpatchrot[2];      /* Sine of Rotation of second patches in 2psc models*/
	double volume;              /* Volume of particle for geometrical center calculations*/
	double pvolscale;           /* Scale of patch volume size*/
	double len[2];              /* Length of the PSC */
	double half_len[2];         /* Half length of the PSC */
	double chiral_cos[2];       /* Coctains the cosinus for the chiral rotation of the patch */
	double chiral_sin[2];       /* Contains the sinus for the chiral rotation of the patch */
};

struct interacts {         /* Parameters pased to functions of interaction calculation */
	double dist;                    /* closest distance */
	struct vector distvec;          /* vector of closes distance */
	struct particles * part1;       /* particle 1 */
	struct particles * part2;       /* particle 2 */
	struct vector box;              /* box size */
	struct ia_param * param;        /* interaction parameters */
	struct vector r_cm;             /* vector connecting center of masses */
	double distcm;                  /* distance between center of masses */
	double dotrcm;                  /* square size of r_cm*/
	double contt;                   /* closest point on spherocylinder to sphere */
};

struct chainparams {      /*Parameters for inner interaction in chains*/
	double bond1eq;          /* Equilibrium distance of harmonic bond between nearest neighbours*/
	double bond1c;           /* Spring constant for harmonic bond between nearest neighbours*/
	double bond2eq;          /* Equilibrium distance of harmonic bond between second nearest neighbours*/
	double bond2c;           /* Spring constant for harmonic bond between second nearest neighbours*/
	double bonddeq;          /* Equilibrium distance of directional harmonic bond between the nearest neighbours*/
	double bonddc;           /* Spring constant for directional harmonic bond between the nearest neighbours*/
	double angle1eq;          /* Equilibrium angle between two spherocylinders -neerest neighbours*/
	double angle1c;           /* Spring constant angle between two spherocylinders -nearest neighbours*/
	double angle2eq;          /* Equilibrium angle between two spherocylinder patches -nearest neighbours*/
	double angle2c;           /* Spring constant for angle between two spherocylinder patches -nearest neighbours*/
};

struct molecule {                 /* This structure is for io only */
	char * name;                    /* The name of the molecule */
	long * type;                    /* The type of the particle */
	long * switchtype;              /* The switchtype of the particle */
	double * delta_mu;              /* The chemical potential for the switch */
};

struct disp {               /* Define step size and acceptance ratio statistics */
	double mx;               /* Maximum value displacement, cos(angle), etc.  */
	double angle;            /* Maximum angle, since in .mx cos(angle) is saved */
	long acc;                /* Number of accepted steps */
	long rej;                /* Number of rejected steps */
	double oldrmsd;          /* Averaged mx value in previous equilibration round */
	double oldmx;            /* Change in mx in last equlibrium step */
};

struct stat {               /* Define statistics counters */
	double sum;
	double sum2;
	long samples;
	double mean;
	double rms;
};

struct meshs {    /* Mesh for hole order parameter  */
	int    dim[2];   /* Mesh dimensions */
	int    *data;    /* Mesh data */
	int    *tmp;     /* tmpporary list for hole search */
};

struct wls {          /* Wang landau method (wl) */
	double *weights;        /* Array of weights for wl method */
	long   *hist;           /* Array of histogram for wl method */
	long   length[2];       /* Length of above arrays */
	double dorder[2];       /* Increments of order parameter */
	double minorder[2];     /* Minimum order parameter */
	double alpha;           /* Current modifier of weights */
	long   currorder[2];    /* Walue of current order parameter*/
	long   neworder[2];     /* wl order parameter in new step */
	long   max;             /* wl maximum of histogram */
	long   min;             /* wl minimum of histogram */
	double wmin;            /* weights minimum */
	int    wlmdim;          /* Dimwnsionality of wang landau */
	int    wlmtype;         /* Atom type for the Wang landau method (wl) */
	double wl_meshsize;         /* Size of mesh bin for hole order paremeter*/
	struct meshs mesh;          /* Mesh for hole order */
	struct meshs origmesh;      /* Mesh store for rejected moves */
	long * radiushole;          /* Array for hole radius around origin */
	long * radiusholeold;       /* Array for hole radius around origin-bigmove */
	long   radiusholemax;         /* Size of array for hole radius*/
	long   partincontact;         /* Number of particles in contact */
	long   partincontactold;      /* Number of particles in contact - old for move*/
};


struct pairs{               /* The structure holding the particle numbers of the pairs and the number of pairs */
	long num_pairs;         /* The number of pairs */
	long * pairs;           /* The paritcle numbers of the paris */
};

struct pairlist{    /* I think, this is done too complicated: just sim->pairs[npart] should be enough */
	struct pairs * list;    /* contains the pairlist of all paritcles */
};

struct cluster{             /* contains all the particles of one cluster */
	long npart;
	long * particles;
};

struct exters{
	BOOL exist;                           /* existence of external potential*/ 
	double thickness;                     /* external wall thicnkess*/
	double epsilon;                       /* depth of attraction*/
	double attraction;                    /* distance of attraction*/ 
	double sqmaxcut;                      /* distance when nothing can interact*/
	struct ia_param interactions[MAXT];   /* Interaction parameters with particle types generated from above params*/
};

struct topo{                /* It would be nice, if this struct would contain all the topo stuff in the end*/
	long * switchlist;                          /* List containing the number of all the particles with switchtypes */
	long n_switch_part;                         /* number of particles with switchtype */
	double sqmaxcut;                            /* square of distance over which even spherocylinders cannot interact (distance between CM) */
	double maxcut;                              /* distance over which even spherocylinders cannot interact (distance between CM) */
	long conlist[MAXN][4];                      /* Connectivity list, we have connection to tail and head and secon neighbours so far*/
	long chainlist[MAXN][MAXCHL];               /* List of chains*/
	long chainnum;                              /* Number of chains */
	struct chainparams chainparam[MAXMT];       /* parameters for chains */
	struct ia_param ia_params[MAXT][MAXT];      /* parametrization of particles for all interations*/
	long npart;                                 /* Number of particles */
	struct exters exter;                        /* external potential - wall */
};

struct sim{                 /* Should contain mostly all the simulation options and variables, that can change in every step. */
	double press;               /* Pressure */
	double paralpress;          /* Parallel pressure for replica exachnge*/
	double dpress;		        /* Pressure change for replica exchange*/
	double shave;               /* Average number of volume changes to attempt per sweep */
	double shprob;              /* Probability of attempting a volume change */
	double chainprob;           /* Average number of chain move attempt per sweep */
	double switchprob;          /* Average number of type switch attempt per sweep */
	int pairlist_update;        /* Number of sweep per upedating the pairlist */
	double temper;              /* Temperature*/
	double paraltemper;         /* Temperature for parallel tempering */
	double dtemp;               /* Temprature step */
	int ptype;                  /* Type of pressure coupling*/
	long adjust;                /* Number of sweeps between step size adjustments */
	long movie;                 /* Number of sweeps between movie frames */
	long nequil;                /* Number of equilibration sweeps */
	long nsweeps;               /* Number of production sweeps */
	long paramfrq;              /* Number of sweeps between order parameter samples */
	long report;                /* Number of sweeps between statistics reports */
	//    long terms;                 /* Number of Fourier terms as smectic order parameters */
	long nrepchange;            /* Number of sweeps between replica exchanges */
	int wlm[2];                 /* Wang landau method (wl) */
	struct disp edge;           /* Maximum box length change and statistics */
	struct disp rot[MAXT];      /* Maximum rotation and statistics */
	struct disp trans[MAXT];    /* Maximum translation  and statistics*/
	struct disp chainm[MAXMT];  /* Maximum translation for chain  and statistics*/
	struct disp chainr[MAXMT];  /* Maximum rotation for chain and statistics */
	struct disp mpiexch;        /* MPI statistics*/
	struct pairs * pairlist;    /* The pairlist */
	long write_cluster;         /* Number of sweeps per writing out cluster info */
	long * clusterlist;         /* clusterlist[i] = cluster index of particle i */
	struct cluster * clusters;  /* informations about the single clusters */
	double *clustersenergy;     /* list of energies of clusters*/
	long num_cluster;           /* number of single clusters */
	long * clusterstat;         /* Statistics about the seize of cluster */
	long max_clust;             /* maximal clustersize */
	struct wls wl;              /* Wang landau data */
	int mpirank;                /* MPI number for given process*/
	int mpinprocs;              /* MPI number of processes */  
};

typedef enum { /* Holds the type of a variable in struct option */
    Int,
    Int2,
    Long,
    Double
} Type;

typedef struct { /* for reading in the options */
    char *id;       /* The name of the value in the option file*/
    Type type;      /* The type (int, double or long) */
    BOOL set;       /* Wheter the variable has been set */
    void *var;      /* The variable */
} Option;

struct conf{                      /* Configuration of the system*/
	struct particles * particle;  /* All particles*/
	struct vector box;            /* Box size*/
	double sysvolume;	          /* Something like total mass*/
	struct vector syscm;          /* System center of mass*/
};

struct filenames {
	/* input files */
	char configurationinfile[30];
	char topologyfile[30];
	char optionsfile[30];
	char wlinfile[30];
	/* output files */
	char configurationoutfile[30];
	char moviefile[30];
	char wloutfile[30];
	char statfile[30];
	char clusterfile[30];
	char clusterstatfile[30];
	char energyfile[30];
};

struct mpiexchangedata{ /* extra type for mpi communication*/
	struct vector box;    /* box of configuration */
	double energy;	      /* energy of configuration */
	double volume;        /* volume of configuration */
	int accepted;         /* bool if accepted */
	struct vector syscm;  /* system CM of configuration */
	long radiusholemax;   /* size of array for WL*/
	long wl_order[2];     /* wang-landau order parameter*/
};

#ifdef MPI
MPI_Datatype MPI_vector, MPI_Particle, MPI_exchange;
#endif

const struct stat nullstat = {0.0, 0.0, 0, 0.0, 0.0};

long seed = 6;             /* Seed for random number generator */


/*..............................................................................*/

int main(int argc, char **argv)
{
	DEBUG("start");
	FILE *outfile,*mov;                                 /* Handle for writing configuration */
	double (* intfce[MAXT][MAXT])(struct interacts *);  /*array of interaction functions*/
	struct topo topo;                                   /* will maybe contain all the topo stuff in future */
	struct sim sim;                                     /* Should contain the simulation options.  */
	struct conf conf;                                   /* Should contain fast changing particle and box(?) information */
	struct filenames files;
	
	int memoryalloc(struct conf * conf);
	int memorydealloc(struct conf * conf, struct topo * topo, struct sim * sim);

	void read_options(struct sim* sim, char filename[30]);
	void init_top(struct topo *, struct conf * conf, struct sim * sim, char filename[30]);
	void init_config(struct topo * topo, struct conf * conf, struct sim * sim, char filename[30]);
	void init_intfce(double (* intfce[MAXT][MAXT])(struct interacts *), struct topo * topo);

	void draw(FILE *, struct conf * conf, struct topo * topo);
	void printeqstat(struct disp *, double, int);
	void simulate(long nsweeps, long adjust, long paramfrq, long report, 
			double (* intfce[MAXT][MAXT])(struct interacts *),
			struct topo * topo, struct sim * sim, struct conf * conf, struct filenames *files);

	void init_pairlist(struct topo * topo, struct sim * sim);
	void gen_pairlist(struct topo * topo, struct sim * sim, struct conf * conf);
	void print_pairlist(FILE * stream, struct sim * sim, struct topo * topo);

	int gen_clusterlist(struct topo * topo, struct sim * sim, struct conf * conf, double (* intfce[MAXT][MAXT])(struct interacts *) );
	int print_clusterlist(FILE * stream, BOOL decor, struct topo * topo, struct sim * sim, struct conf * conf);
	int print_clusters(FILE * stream, BOOL decor, struct sim * sim);
	int print_clusterstat(FILE * stream, BOOL decor, struct sim * sim);
	int sort_clusterlist(struct topo * topo, struct sim * sim);

	printf ("\nPatchy Spherocylinders version 3.5 ");
	sprintf(files.configurationinfile, "config.init");
	sprintf(files.configurationoutfile, "config.last");
	sprintf(files.optionsfile, "options");
	sprintf(files.topologyfile, "top.init");
	sprintf(files.moviefile, "movie");
	sprintf(files.wlinfile, "wl.dat");
	sprintf(files.wloutfile, "wl-new.dat");
	sprintf(files.statfile, "stat.dat");
	sprintf(files.clusterfile, "cluster.dat");
	sprintf(files.clusterstatfile, "cluster_stat.dat");
	sprintf(files.energyfile, "energy.dat");
#ifdef MPI
    FILE *infile;
	printf(" MPI version");
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &(sim.mpinprocs) );                                                                                         
	MPI_Comm_rank(MPI_COMM_WORLD, &(sim.mpirank) ); 

	sprintf(files.configurationoutfile, "%dconfig.last", sim.mpirank);
	sprintf(files.moviefile, "%dmovie", sim.mpirank);
	sprintf(files.wloutfile, "%dwl-new.dat", sim.mpirank);
	sprintf(files.clusterfile, "%dcluster.dat", sim.mpirank);
	sprintf(files.clusterstatfile, "%dcluster_stat.dat", sim.mpirank);
	sprintf(files.energyfile, "%denergy.dat", sim.mpirank);
	sprintf(files.statfile, "%dstat.dat", sim.mpirank);
	
	/*test if there is a specific input configuration for mpi run*/
	sprintf(files.configurationinfile, "%dconfig.init", sim.mpirank);
	infile = fopen(files.configurationinfile, "r");
	if (infile != NULL) 
		fclose (infile);
	else  sprintf(files.configurationinfile, "config.init");

	/*test if there is a specific input wang-landau for mpi run*/
	sprintf(files.wlinfile, "%dwl.dat", sim.mpirank);
	infile = fopen(files.wlinfile, "r");
	if (infile != NULL) 
		fclose (infile);
	else  sprintf(files.wlinfile, "wl.dat");
#endif 
	
	printf ("\n-------------------------------------\n");
	printf ("Reading options...\n");
	read_options(&sim,files.optionsfile);
	init_top(&topo, &conf, &sim,files.topologyfile);
	if (topo.chainnum ==0) {
	    /*no chain make the probability of moving them 0*/
	    if (sim.chainprob > 0)
		printf ("No chains... chain move probability set to 0.\n");
	    sim.chainprob = 0;  
	}
	printf ("\nReading configuration...\n");
	init_config(&topo, &conf, &sim, files.configurationinfile);
	printf ("Equilibration of maximum step sizes:  %ld sweeps\n", sim.nequil/2);
	fflush (stdout);
	if ( sim.wlm[0] > 0 ) {
		outfile = fopen(files.wlinfile, "r");
		if (outfile == NULL) {
			printf ("ERROR: Cannot open file for Wang-Landau method (%s).\n",files.wlinfile);
			memorydealloc(&conf, &topo, &sim);
			exit(1);
		}
		fclose (outfile);
	}
	/* Empty movie file */
	mov = fopen("movie", "w");
	fclose (mov);
	printf ("\nInitializing energy functions...\n");
	init_intfce(intfce, &topo);
	if (sim.pairlist_update) {
		init_pairlist(&topo, &sim);
	}
	if (sim.nequil) {
		printf("\nStart equilibration...\n");
		simulate(sim.nequil/2, sim.adjust, 0, 0, intfce, &topo, &sim, &conf,&files);
		simulate(sim.nequil/2, 0,          0, 0, intfce, &topo, &sim, &conf,&files);
		printf ("   Equilibrated maximum displacement / acceptance ratio:            \n");
		printeqstat(sim.trans,2.0,MAXT);
		printf ("   Equilibrated maximum rotation / acceptance ratio:                       \n");
		printeqstat(sim.rot,1.0,MAXT);
		printf ("   Equilibrated maximum box length change / acceptance ratio:              \n");
		printf ("                     %.6le  /  %.6le\n", sim.edge.mx/2.0,RATIO(sim.edge));
		printf ("   Equilibrated maximum displacement of chain / acceptance ratio:   \n");
		printeqstat(sim.chainm,2.0,MAXMT);
		printf ("   Equilibrated maximum rotation of chain / acceptance ratio:              \n");
		printeqstat(sim.chainr,1.0,MAXMT);
		printf ("\n");

		printf ("Further equilibration of configuration:  %ld sweeps\n", sim.nequil/2);
		fflush (stdout);

		outfile = fopen("config.eq", "w");
		fprintf (outfile, "%15.8le %15.8le %15.8le\n", conf.box.x, conf.box.y, conf.box.z);
		draw (outfile, &conf, &topo);
		fclose (outfile);
		printf ("   Equilibrated configuration written to config.eq\n");
		printf ("   Box dimensions: %.10lf, %.10lf, %.10lf\n\n", conf.box.x, conf.box.y, conf.box.z);
	}
	
	printf ("Production run:  %ld sweeps\n\n", sim.nsweeps);
	fflush (stdout);

	simulate(sim.nsweeps, 0, sim.paramfrq, sim.report, intfce, &topo, &sim, &conf,&files);
#ifdef MPI
	printf ("   MPI replica changeT / changeP / acceptance ratio: \t %.6lf   /   %.6lf  /  %.6lf\n\n", sim.mpiexch.mx,sim.mpiexch.angle,RATIO(sim.mpiexch));
#endif
	outfile = fopen(files.configurationoutfile, "w");
	fprintf (outfile, "%15.8le %15.8le %15.8le\n", conf.box.x, conf.box.y, conf.box.z);
	draw (outfile, &conf, &topo);
	fclose (outfile);

	// For testing the pairlist
	//gen_pairlist(&topo, &sim, &conf);
	//FILE * fpairlist;
	//fpairlist = fopen("pairlist.dat", "w");
	//print_pairlist(fpairlist, &sim, &topo);
	//fclose(fpairlist);
	//printf("sqmaxcut = %lf\n", topo.sqmaxcut);

	//// For testing the cluster algorithm
	//gen_clusterlist(&topo, &sim, &conf);
	//print_clusterlist(stdout, TRUE, &topo, &sim, &conf);
	//sort_clusterlist(&topo, &sim);
	//print_clusters(stdout, TRUE, &sim);
	//print_clusterstat(stdout, TRUE, &sim);
	
	if (memorydealloc(&conf, &topo, &sim))  
		exit(1);
#ifdef MPI
	MPI_Finalize();
#endif
	printf ("\nDone\n\n");

	return 0;
}

/*..............................................................................*/
/*.........................SIMULATION RUN.......................................*/
/*..............................................................................*/

void simulate(long nsweeps, long adjust, long paramfrq, long report, 
		double (* intfce[MAXT][MAXT])(struct interacts *),
		struct topo * topo, struct sim * sim, struct conf * conf, struct filenames *files)
{

	long i,j,wli;
	long next_adjust;  /* Next sweep number for step size adjustment */
	long next_calc;    /* Next sweep number for order parameter calculation */
	long next_dump;    /* Next sweep number for reporting statistics */
	long next_frame;   /* Next sweep number for dumping a movie frame */
	long step;         /* Step number within a given sweep */
	long sweep;        /* Current sweep number */
	//struct stat nem;   /* Nematic order parameter */
	//struct stat vol;   /* Volume statistics */
	//struct stat shapex, shapey, shapez;   /* Box shape statistics */
	//struct stat smec[MAXF];       /* Smectic order parameters (Fourier coeeficients) */
	FILE *mf;              /* Handle for movie file */
	FILE *cl_stat, *cl, *cl_list; /* Handle for cluster statistics */
	FILE *ef, *statf;              /* Handle for energy file and statistical file*/
	double edriftstart; /* Energy drift calculation - start */
	double edriftchanges;  /* Energy drift calculation - accumulate all changes through moves */
	double edriftend;   /* Energy drift calculation - end */
	double pvdriftstart; /* PV drift calculation - start */
	double pvdriftend;   /* PV drift calculation - end */
	double volume;		 /* volume of box*/
	double moveprobab;   /* random number selecting the move*/
	
	/* function declarations */
	//double nematic(long, struct particles *);
	double ran2(long *);
	//double smectic(long, struct particles *, long);
	double calc_energy(long target, double (* intfce[MAXT][MAXT])(struct interacts *), 
			int mode, struct topo * topo, struct conf * conf, struct sim * sim,int chainn);
	void accumulate(struct stat *, double);
	void draw(FILE *, struct conf * conf, struct topo * topo);
	void optimizestep(struct disp *, double, double);
	void optimizerot(struct disp *, double, double);
	void partvecinit(struct topo * topo, struct sim * sim, struct conf * conf );
	int wlinit(struct wls *, char filename[30]);
	int wlwrite(struct wls *, char filename[30]);
	int wlend(struct wls *);
	int mesh_init(struct meshs *, double, long, struct conf * conf, struct sim * sim);
	int mesh_end(struct meshs *);
	long z_order(struct wls *, struct conf * conf,int);
	long twopartdist(struct wls *, struct conf * conf,int);
	void mesh_print (struct meshs *);
	void masscenter(long, struct ia_param [MAXT][MAXT], struct conf * conf);
	void gen_pairlist(struct topo * topo, struct sim * sim, struct conf * conf);
	int write_cluster(FILE * cl_stat, FILE * cl, FILE * cl_list, 
			BOOL decor, long sweep,	struct sim * sim, struct topo * topo, 
			struct conf * conf, double (* intfce[MAXT][MAXT])(struct interacts *));
	double particlemove(struct topo * topo, struct sim * sim, struct conf * conf,
						double (* intfce[MAXT][MAXT])(struct interacts *));
	double chainmove(struct topo * topo, struct sim * sim, struct conf * conf,
					 double (* intfce[MAXT][MAXT])(struct interacts *));
	double switchtypemove(struct topo * topo, struct sim * sim, struct conf * conf,
						  double (* intfce[MAXT][MAXT])(struct interacts *));
	double pressuremove(struct topo * topo, struct sim * sim, struct conf * conf,
						double (* intfce[MAXT][MAXT])(struct interacts *));
	double replicaexchangemove(struct topo * topo, struct sim * sim, struct conf * conf,
						double (* intfce[MAXT][MAXT])(struct interacts *), long sweep);
	long radiushole_all(struct topo *topo, struct conf *conf, struct sim * sim,int, struct vector *);
	long radiushole_position(double, struct sim *,int);
	long contparticles_all(struct topo *topo, struct conf *conf, struct sim * sim,int wli);
	double alignment_order(struct conf * conf, struct topo * topo);
    int memorydealloc(struct conf * conf, struct topo * topo, struct sim * sim);

    double paire(long, long, double (* intfce[MAXT][MAXT])(struct interacts *),
                    struct topo * topo, struct conf * conf);

	/* Opening files for cluster statistics */
	cl_stat = cl = cl_list = ef = statf = NULL;
	if(sim->write_cluster){
		// Empty file
		cl_stat = fopen(files->clusterstatfile, "w");
		fclose(cl_stat);
		cl_stat = fopen(files->clusterstatfile, "a");
		// Empty file
		cl = fopen(files->clusterfile, "w");
		fclose(cl);
		cl = fopen(files->clusterfile, "a");
	}
	/* write energy*/
	if (report <= nsweeps){
		// Empty file
		ef = fopen(files->energyfile, "w");
		fclose(ef);
		ef = fopen(files->energyfile, "a");
		fprintf (ef, "# sweep    energy\n");
		statf = fopen(files->statfile, "w");
		fclose(statf);
		statf = fopen(files->statfile, "a");
		fprintf (statf, "# sweep    volume\n");
	}

	/*=== Initialise counters etc. ===*/
//	double pvolume;    /* Volume of all particles*/
/*	pvolume =0.0;
	for (i=0;i < topo->npart;i++) {
		if  (conf->particle[i].type>=0 )  
			pvolume += topo->ia_params[conf->particle[i].type][conf->particle[i].type].volume;
	}*/
	sim->shprob = sim->shave/(double)topo->npart;
	for (i=0;i<MAXT;i++){
		sim->rot[i].acc = 0;
		sim->rot[i].rej = 0;
		sim->rot[i].oldrmsd = 0;
		sim->rot[i].oldmx = 0;
		sim->trans[i].acc = 0;
		sim->trans[i].rej = 0;
		sim->trans[i].oldrmsd = 0;
		sim->trans[i].oldmx = 0;
	}
	for (i=0;i<MAXMT;i++){
		sim->chainm[i].acc = 0;
		sim->chainm[i].rej = 0;
		sim->chainm[i].oldrmsd = 0;
		sim->chainm[i].oldmx = 0;
		sim->chainr[i].acc = 0;
		sim->chainr[i].rej = 0;
		sim->chainr[i].oldrmsd = 0;
		sim->chainr[i].oldmx = 0;
	}
	//(*edge).acc = (*edge).rej = (*edge).oldrmsd = (*edge).oldmx = 0;
	sim->edge.acc = sim->edge.rej = sim->edge.oldrmsd = sim->edge.oldmx = 0;
	sim->mpiexch.acc = sim->mpiexch.rej = sim->mpiexch.oldrmsd = sim->mpiexch.oldmx = 0;

	/*Initialize some values at begining*/
	partvecinit(topo,sim,conf);
	next_adjust = adjust;
	next_calc = paramfrq;
	next_dump = report;
	next_frame = sim->movie;
	//nem = vol = shapex = shapey = shapez = nullstat;
	//for (i=0; i<MAXF; i++) smec[i] = nullstat;
	if (sim->movie > 0) {
		mf = fopen(files->moviefile, "a");
	} else {
		mf = NULL;
	}
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
	masscenter(topo->npart,topo->ia_params, conf);
	/* Initialization of wang-landaou method*/
	if ( sim->wlm[0] >0 ) {  
		if (wlinit(&sim->wl,files->wlinfile) != 0) 
			return;
		sim->wl.wlmdim = 1 ;
		if ( sim->wlm[1] > 0 )
			sim->wl.wlmdim = 2 ;
		for (wli=0;wli<sim->wl.wlmdim;wli++) {
			switch (sim->wlm[wli]) {
				case 1: 
					masscenter(topo->npart,topo->ia_params, conf);
					sim->wl.currorder[wli] = z_order(&sim->wl,conf,wli);
					break;
				case 2:
					sim->wl.wl_meshsize = (topo->ia_params[sim->wl.wlmtype][sim->wl.wlmtype].sigma) / 3.0; // TODO
					sim->wl.mesh.data = NULL;
					sim->wl.mesh.tmp = NULL;
					sim->wl.origmesh.data = NULL;
					sim->wl.origmesh.tmp = NULL;
					sim->wl.currorder[wli] = (long) (mesh_init(&sim->wl.mesh,sim->wl.wl_meshsize, topo->npart, conf, sim) - sim->wl.minorder[wli]);
					break;
				case 3:
					sim->wl.currorder[wli] = (long) floor( (conf->particle[0].dir.z - sim->wl.minorder[wli])/ sim->wl.dorder[wli] );
					break;
				case 4:
					sim->wl.currorder[wli] = twopartdist(&sim->wl,conf,wli);
					break;
				case 5:
					masscenter(topo->npart,topo->ia_params, conf);
					sim->wl.radiusholemax = 0;
					sim->wl.radiushole = NULL;
					sim->wl.radiusholeold = NULL;
					sim->wl.currorder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->syscm));
					break;
				case 6:
					sim->wl.radiusholemax = 0;
					sim->wl.radiushole = NULL;
					sim->wl.radiusholeold = NULL;
					sim->wl.currorder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
					break;
				case 7:	
					sim->wl.currorder[wli] = contparticles_all(topo,conf,sim,wli);
					break;
				default:
					sim->wl.currorder[wli] = 0;
					break;
			}
			if ( (sim->wl.currorder[wli] >= sim->wl.length[wli] ) || (sim->wl.currorder[wli] < 0) ) {
				printf("Error: starting Wang-Landau method with order parameter %f out of range(%f - %f)\n\n", sim->wl.dorder[wli]*sim->wl.currorder[wli] + \
				   sim->wl.minorder[wli], sim->wl.minorder[wli], sim->wl.minorder[wli]+sim->wl.dorder[wli]*sim->wl.length[wli]  );
				wlend(&sim->wl);
				if (memorydealloc(conf, topo, sim))  
		            exit(1);
				exit(2);
				return;
			}
		}
		if (sim->wl.alpha < WL_ALPHATOL/100) sim->wl.alpha = WL_ZERO;
		fflush (stdout);
	}

double e;
        for(int i=0; i< 1/*topo->npart-1*/; ++i) {
            for(int j=i+1; j< topo->npart; ++j) {
                e = paire(i, j, intfce, topo, conf);
                if(e < 1000.0)
                    printf("%.5lf\n", e);
                else printf("%lf\n", 1000.0);
            }
            printf("\n");
        }
        exit(0);


	/*do moves - START OF REAL MC*/
	if(sim->pairlist_update){
		gen_pairlist(topo, sim, conf); // Does that solve the problem?
	}

	/*do energy drift check - start calculation*/
	volume = conf->box.x * conf->box.y * conf->box.z;
	edriftstart = calc_energy(0, intfce, 0, topo, conf, sim,0);
	pvdriftstart = sim->press * volume - (double)topo->npart * log(volume) / sim->temper;
	//printf("starting energy: %.15f \n",calc_energy(0, intfce, 0, topo, conf, sim,0)); 
	//printf("press: %.15f\n",sim->press * volume - (double)topo->npart * log(volume) / sim->temper);
	edriftchanges = 0.0;


	for (sweep=1; sweep <= nsweeps; sweep++) {
		// Try replica exchange
		if((sim->nrepchange) && (sweep % sim->nrepchange == 0)){
			edriftchanges += replicaexchangemove(topo,sim,conf,intfce,sweep);
		}
		// Generate the pairlist
		if((sim->pairlist_update) && (sweep % sim->pairlist_update == 0)){
			gen_pairlist(topo, sim, conf);
		}
		//normal moves
		for (step=1; step <= topo->npart; step++) {
			moveprobab = ran2(&seed);
			if ( moveprobab < sim->shprob) {
				/* pressure moves*/
				edriftchanges += pressuremove(topo,sim,conf,intfce);
			} else {
				if (moveprobab < sim->shprob + sim->chainprob) {
					/* single particle moves*/
					edriftchanges += chainmove(topo,sim,conf,intfce);
				}
				else if (moveprobab < sim->shprob + sim->chainprob + sim->switchprob){
					/*=== This is an attempt to switch a type ===*/
					edriftchanges += switchtypemove(topo,sim,conf,intfce);
					
				} else {
					/* single particle moves*/
					edriftchanges += particlemove(topo,sim,conf,intfce);
					
				} /* end of else next to chain moves */ 
			} /* end of else next to volume moves */
		}
		/**** End of step loop for this sweep ****/

		/*=== Start of end-of-sweep housekeeping ===*/

		/* Adjustment of maximum step sizes during equilibration */
		if (sweep == next_adjust) {
			for (i = 0; i < MAXT ;i++) {
				if ((sim->trans[i].acc > 0)||(sim->trans[i].rej >0)) 
					optimizestep (sim->trans + i, 1.5, 0.0);
				if ((sim->rot[i].acc > 0)||(sim->rot[i].rej >0)) 
					optimizerot (sim->rot + i, 5.0, 0.01);
			}
			for (i = 0; i < MAXMT; i++) {
				if ((sim->chainm[i].acc > 0)||(sim->chainm[i].rej > 0)) 
					optimizestep (sim->chainm + i, 1.5, 0.0);
				if ((sim->chainr[i].acc > 0)||(sim->chainr[i].rej > 0)) 
					optimizerot (sim->chainr + i, 5.0, 0.01);
			}
			optimizestep (&(sim->edge), 1.0, 0.0);
			next_adjust += adjust;
		}

		if ( (sim->wlm[0] > 0) && (sim->wl.alpha > WL_ZERO) && !(sweep % 1000) ) {
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
					  printf (" %15.8le %15ld %15.8f\n",sim->wl.weights[i],sim->wl.hist[i],particle[0].pos.z);
					  fflush(stdout);
					  }
					 */
					if ( sim->wl.alpha < WL_ALPHATOL) {
						printf("\nF I N I S H E D\n\n");
						fflush (stdout);
						break;
					}
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

		if (!(sweep % 10000)) {
		    /*reinitialize pach vectors to avoid cummulation of errors*/
		    partvecinit(topo,sim,conf);
		}
		
		/* Sampling of statistics */
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
		/* Writing of statistics */
		if (sweep == next_dump) {
			/*printf ("Statistics after %ld sweeps:\n", sweep);
			  printf ("   Mean and RMS fluctuation of S2:  %13.8lf %13.8lf\n",
			  nem.mean, nem.rms);
			  for (i=0; i<terms; i++) {
			  printf ("   Mean & fluc. Fourier coeff. %3ld: %13.8lf %13.8lf\n",
			  i+1, smec[i].mean, smec[i].rms);
			  }
			  printf ("   Mean & fluc box dimensions:  x   %13.8lf %13.8lf\n",
			  shapex.mean, shapex.rms);
			  printf ("                                y   %13.8lf %13.8lf\n",
			  shapey.mean, shapey.rms);
			  printf ("                                z   %13.8lf %13.8lf\n",
			  shapez.mean, shapez.rms);
			  printf ("   Mean & fluctuation volume:      %13.8lf %13.8lf\n",
			  vol.mean, vol.rms);
			  printf ("   Mean & fluc. volume over volume of particles:    %13.8lf %13.8lf\n",
			  vol.mean/pvolume, vol.rms/pvolume);
			  printf ("\n");
			  fflush (stdout);
			 */
			fprintf (statf, " %ld; %.10lf\n", sweep, conf->box.x * conf->box.y * conf->box.z);
			fprintf (ef, " %ld; %.10lf  %f \n", sweep, calc_energy(0, intfce, 0, topo, conf, sim,0), alignment_order(conf,topo));
			if (sim->wlm[0] > 0) {
				wlwrite(&sim->wl,files->wloutfile);
			}
			next_dump += report;
		}

		/* Writing of movie frame */
		if (sweep == next_frame) {
			fprintf (mf, "%ld\n", topo->npart);
			fprintf (mf, "sweep %ld;  box %.10lf %.10lf %.10lf\n", sweep, conf->box.x, conf->box.y, conf->box.z);
			draw (mf, conf, topo);
			fflush (mf);
			next_frame += sim->movie;
		}

		/* Writing out cluster statistics */
		if(sim->write_cluster && (sweep % sim->write_cluster == 0)){
			write_cluster(cl_stat, cl, cl_list, FALSE, sweep, sim, topo, conf, intfce);
		}
		/*=== End of housekeeping ===*/

	}
	/**** End of sweeps loop ****/

	/*do energy drift check - at the end calculation*/
	volume = conf->box.x * conf->box.y * conf->box.z;
	edriftend = calc_energy(0, intfce, 0, topo, conf, sim,0); 
	pvdriftend =  sim->press * volume - (double)topo->npart * log(volume) / sim->temper;
	printf("Energy drift: %.15lf \n",edriftend - edriftstart - edriftchanges +pvdriftend -pvdriftstart);
        printf("Starting energy+pv: %.8lf \n",edriftstart+pvdriftstart);
	printf("Starting energy: %.8lf \n",edriftstart);
	fflush(stdout);

	/* End wang-landau*/
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
		wlwrite(&sim->wl,files->wloutfile);
		wlend(&sim->wl);
		if ( (sim->wlm[0] == 2)||(sim->wlm[1] == 2) ) {
			mesh_end(&sim->wl.mesh);
			mesh_end(&sim->wl.origmesh);
		}
		if ( (sim->wlm[0] == 5)||(sim->wlm[1] == 5)||(sim->wlm[0] == 6)||(sim->wlm[1] == 6)  ) {
			if ( sim->wl.radiushole != NULL ) free(sim->wl.radiushole);
			if ( sim->wl.radiusholeold != NULL ) free(sim->wl.radiusholeold);
		}
	}
	/*end movie*/
	if (sim->movie > 0) 
		fclose (mf);
	/*end cluster*/
	if(sim->write_cluster){
		fclose(cl_stat);
		fclose(cl);
	}
	if (report < nsweeps) {
		fclose(ef);
		fclose(statf);
	}

}

/*..................................MOVES.........................................*/
/*................................................................................*/

/*..............................PARTICLE MOVES....................................*/
double particlemove(struct topo * topo, struct sim * sim, struct conf * conf,
				  double (* intfce[MAXT][MAXT])(struct interacts *))
{
	double edriftchanges =0.0;
	long target;
	
	double ran2(long *);
	double partdisplace(struct topo * topo, struct sim * sim, struct conf * conf,
						double (* intfce[MAXT][MAXT])(struct interacts *),long target);
	double partrotate(struct topo * topo, struct sim * sim, struct conf * conf,
					  double (* intfce[MAXT][MAXT])(struct interacts *),long target);
	
	/*=== This is a particle move step ===*/
	target = ran2(&seed) * topo->npart;
	if ( ((ran2(&seed) < 0.5) || (topo->ia_params[conf->particle[target].type][conf->particle[target].type].geotype[0] >= SP)) ) { /* no rotation for spheres */
		//target = 1;
		//printf ("displacement\n\n");
		edriftchanges = partdisplace(topo,sim,conf,intfce,target);
		
	} else {
		/*=== Rotation step ===*/
		edriftchanges = partrotate(topo,sim,conf,intfce,target);
	}
	/*=== End particle move step ===*/
	return edriftchanges;
}

/*................................................................................*/
double partdisplace(struct topo * topo, struct sim * sim, struct conf * conf,
				  double (* intfce[MAXT][MAXT])(struct interacts *),long target)
{
	double edriftchanges,energy,enermove,wlener;
	struct vector orig, dr, origsyscm;
	int reject=0,wli;
	double radiusholemax_orig=0;
	
	double calc_energy(long target, double (* intfce[MAXT][MAXT])(struct interacts *), 
					   int mode, struct topo * topo, struct conf * conf, struct sim * sim,int chainn);
	int movetry(double, double, double);
	void wlreject(struct sim *, long);
	void wlaccept(int, struct wls *);
	long meshorder_moveone(struct vector, struct vector, struct meshs *, long, long, \
						   struct conf * conf, struct sim * sim, int wli);
	int mesh_cpy(struct meshs *, struct meshs *);
	//void mesh_print (struct meshs *);
	long z_order(struct wls *, struct conf * conf, int wli);
	long twopartdist(struct wls *, struct conf *conf, int wli);
	struct vector ranvec(void);
	int longarray_cpy (long **, long **, long, long);
	long radiusholeorder_moveone(struct vector *oldpos, struct conf *conf, struct sim * sim, long target, int wli,struct vector *);
	long radiushole_all(struct topo *topo, struct conf *conf, struct sim * sim,int, struct vector *);
	long contparticles_moveone(struct vector *oldpos, struct conf *conf, struct sim * sim, long target,int wli);
	long contparticles_all(struct topo *topo, struct conf *conf, struct sim * sim,int wli);
	
	/*=== Displacement step ===*/
	edriftchanges =0.0;
	origsyscm.x = 0;
	origsyscm.y = 0;
	origsyscm.z = 0;
	energy = calc_energy(target, intfce, 1, topo, conf, sim,0);
	orig = conf->particle[target].pos;
	dr = ranvec();
	//ran = sqrt(ran2(&seed));
	dr.x *= sim->trans[conf->particle[target].type].mx/conf->box.x; 
	dr.y *= sim->trans[conf->particle[target].type].mx/conf->box.y;
	dr.z *= sim->trans[conf->particle[target].type].mx/conf->box.z;
	conf->particle[target].pos.x += dr.x;
	conf->particle[target].pos.y += dr.y;
	conf->particle[target].pos.z += dr.z;
	//} while (conf->particle[target].pos.x < 0.25 || conf->particle[target].pos.x > 0.50);
	
	reject = 0;
	wlener = 0.0;
	if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
		for (wli=0;wli<sim->wl.wlmdim;wli++) {
			switch (sim->wlm[wli]) {
				case 1: origsyscm = conf->syscm;
					conf->syscm.x += dr.x * topo->ia_params[conf->particle[target].type][conf->particle[target].type].volume / conf->sysvolume;
					conf->syscm.y += dr.y * topo->ia_params[conf->particle[target].type][conf->particle[target].type].volume / conf->sysvolume;
					conf->syscm.z += dr.z * topo->ia_params[conf->particle[target].type][conf->particle[target].type].volume / conf->sysvolume;
					sim->wl.neworder[wli] = z_order(&sim->wl, conf,wli);
					break;
				case 2: mesh_cpy(&sim->wl.origmesh,&sim->wl.mesh);
					sim->wl.neworder[wli] = meshorder_moveone(orig, conf->particle[target].pos, &sim->wl.mesh, topo->npart, target, conf, sim,wli);
					break;
				case 4:
					sim->wl.neworder[wli] = twopartdist(&sim->wl,conf,wli);
					break;
				case 5:
					radiusholemax_orig = sim->wl.radiusholemax;
					origsyscm = conf->syscm;
					conf->syscm.x += dr.x * topo->ia_params[conf->particle[target].type][conf->particle[target].type].volume / conf->sysvolume;
					conf->syscm.y += dr.y * topo->ia_params[conf->particle[target].type][conf->particle[target].type].volume / conf->sysvolume;
					conf->syscm.z += dr.z * topo->ia_params[conf->particle[target].type][conf->particle[target].type].volume / conf->sysvolume;
					longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
					sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->syscm));
					break;
				case 6:
					radiusholemax_orig = sim->wl.radiusholemax;
					longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
					if ( target == 0 )
					  sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
					else 
					  sim->wl.neworder[wli] = radiusholeorder_moveone(&orig, conf, sim,target,wli,&(conf->particle[0].pos));
					break;
				case 7:
					sim->wl.partincontactold = sim->wl.partincontact;
					if ( target == 0 )
						sim->wl.neworder[wli] = contparticles_all(topo,conf,sim,wli);
					else
						sim->wl.neworder[wli] = contparticles_moveone(&orig,conf,sim,target,wli);
					break;
				default: 
					sim->wl.neworder[wli] = sim->wl.currorder[wli];
					break;
			}
			if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
		}
		if (!reject) { 
			wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
			energy += wlener;
		}
		
	}
	
	if (!reject) {  /* wang-landaou ok, try move - calcualte energy */
		enermove =  calc_energy(target, intfce, 1, topo, conf, sim,0);
	}
	if ( reject || movetry(energy, enermove, sim->temper) ) {  /* probability acceptance */
		conf->particle[target].pos = orig;
		sim->trans[conf->particle[target].type].rej++;
		if ( (sim->wlm[0] == 1) || (sim->wlm[0] == 5) || (sim->wlm[1] == 1) || (sim->wlm[1] == 5) ) 
			conf->syscm = origsyscm;
		wlreject(sim,radiusholemax_orig);
	} else { /* move was accepted */
		sim->trans[conf->particle[target].type].acc++;
		wlaccept(sim->wlm[0],&sim->wl);
		edriftchanges = enermove - energy + wlener;
		//printf("%lf\t%lf\n", conf->particle[0].pos.z * conf->box.z , enermove);
		//printf("%.12f\t%.12f\t%.12f\n", energy , enermove,edriftchanges);
	}
	
	return edriftchanges;
}
/*................................................................................*/
double partrotate(struct topo * topo, struct sim * sim, struct conf * conf,
				double (* intfce[MAXT][MAXT])(struct interacts *),long target)
{
	double edriftchanges,energy,enermove,wlener;
	struct particles origpart;
	int reject=0,wli;
	
	double calc_energy(long target, double (* intfce[MAXT][MAXT])(struct interacts *), 
					   int mode, struct topo * topo, struct conf * conf, struct sim * sim,int chainn);
	int movetry(double, double, double);
	void wlreject(struct sim *,long);
	void wlaccept(int, struct wls *);
	void normalise(struct vector *);
	void ortogonalise(struct vector *,struct vector);
	void psc_rotate(struct particles *,double,int);
	
	/*=== Rotation step ===*/
	//printf ("rotation %ld npart %ld\n\n",target,npart);
	energy = calc_energy(target, intfce, 1, topo, conf, sim,0);
	origpart = conf->particle[target];
	psc_rotate(&conf->particle[target],sim->rot[conf->particle[target].type].angle, topo->ia_params[conf->particle[target].type][conf->particle[target].type].geotype[0]);
	/*should be normalised and ortogonal but we do for safety*/
	normalise (&conf->particle[target].dir);
	ortogonalise(&conf->particle[target].patchdir[0],conf->particle[target].dir);

	reject = 0;
	edriftchanges =0.0;
	wlener = 0.0;
	if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
		for (wli=0;wli<sim->wl.wlmdim;wli++) {
			switch (sim->wlm[wli]) {
				case 3: 
					if (target == 0)  sim->wl.neworder[wli] = (long) floor( (conf->particle[0].dir.z - sim->wl.minorder[wli])/ sim->wl.dorder[wli] );
					else sim->wl.neworder[wli] = sim->wl.currorder[wli];
					/* only rotation change direction */
					break;
				default: 
					sim->wl.neworder[wli] = sim->wl.currorder[wli];
					break;
			}
			if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
		}
		if (!reject) { 
			wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
			energy += wlener;
		}
	}
	
	if (!reject) {  /* wang-landaou ok, try move - calcualte energy */
		enermove =  calc_energy(target, intfce, 1, topo, conf, sim,0);
	}
	if ( reject || movetry(energy,enermove,sim->temper) ) {  /* probability acceptance */
		conf->particle[target] = origpart;
		sim->rot[conf->particle[target].type].rej++;
		wlreject(sim,sim->wl.radiusholemax);
	} else { /* move was accepted */
		// DEBUG
		//fprintf(fenergy, "%lf\t%lf\n", conf->particle[1].pos.x * conf->box.x , enermove);
		sim->rot[conf->particle[target].type].acc++;
		wlaccept(sim->wlm[0],&sim->wl);
		edriftchanges = enermove - energy + wlener;
		//printf("%lf\t%lf\n", conf->particle[0].patchdir[0].z, enermove);
	}
	
	return edriftchanges;
}

/*..................... This is an attempt to switch a type.................................*/
double switchtypemove(struct topo * topo, struct sim * sim, struct conf * conf,
					  double (* intfce[MAXT][MAXT])(struct interacts *) )
{
	double edriftchanges,energy,enermove,wlener;
	int reject=0,wli;
	long target;
	double radiusholemax_orig=0;
	
	double ran2(long *);
	double calc_energy(long target, double (* intfce[MAXT][MAXT])(struct interacts *), 
					   int mode, struct topo * topo, struct conf * conf, struct sim * sim,int chainn);
	int movetry(double, double, double);
	void wlreject(struct sim *,long);
	void wlaccept(int, struct wls *);
	void int_partvec(long, struct ia_param *, struct conf *);
	int mesh_init(struct meshs *, double, long, struct conf * conf, struct sim * sim);
	int mesh_cpy(struct meshs *, struct meshs *);
	long z_order(struct wls *, struct conf * conf,int);
	long twopartdist(struct wls *, struct conf *conf,int);
	long radiushole_all(struct topo *topo, struct conf *conf, struct sim * sim,int, struct vector *);
	int longarray_cpy (long **target, long **source,long,long);
	long contparticles_all(struct topo *topo, struct conf *conf, struct sim * sim,int wli);
	
	/*=== This is an attempt to switch a type ===*/
	edriftchanges =0.0;
	wlener = 0.0;
	target = ran2(&seed) * topo->n_switch_part;
	target = topo->switchlist[target];
	DEBUG_SIM("Switching the particle type");
	DEBUG_SIM("PARTICLE: %ld", target);
	energy = calc_energy(target, intfce, 1, topo, conf, sim,0);
	// Start switching the type
	int switched = conf->particle[target].switched;
	int pmone = PMONE(switched);
	DEBUG_SIM("switched = %d", switched);
	DEBUG_SIM("pmone = %d", pmone);
	int tmp_type = conf->particle[target].type;
	conf->particle[target].type = conf->particle[target].switchtype;
	conf->particle[target].switchtype = tmp_type; 
	conf->particle[target].switched += pmone;
	int_partvec(target,&(topo->ia_params[conf->particle[target].type][conf->particle[target].type]),conf);
	DEBUG_SIM("Particle %ld is %d switched", target, switched);
	//DEBUG
#ifdef DEBUGGING_SIM
	if ((abs(pmone) != 1) || (conf->particle[target].type == conf->particle[target].switchtype)){
		fprintf(stderr, "ERROR: Something went wrong, when switching the type of particle %ld\n", target);
		exit(1);
	}
#endif
	if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
		for (wli=0;wli<sim->wl.wlmdim;wli++) {
			switch (sim->wlm[wli]) {
				/*case 1: sim->wl.neworder = z_order(&sim->wl, conf,wli);
					 break;*/
				case 2: mesh_cpy(&sim->wl.origmesh,&sim->wl.mesh);
					sim->wl.neworder[wli] = (long) (mesh_init(&sim->wl.mesh,sim->wl.wl_meshsize,topo->npart, conf, sim) - sim->wl.minorder[wli]);
					break;
				/*case 4:
					sim->wl.neworder = twopartdist(&sim->wl,conf,wli);
					break;*/
				case 5:
					radiusholemax_orig = sim->wl.radiusholemax;
					longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
					sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->syscm));
					break;
				case 6:
					radiusholemax_orig = sim->wl.radiusholemax;
					longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
					sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
					break;
				case 7:
					sim->wl.partincontactold = sim->wl.partincontact;
					sim->wl.neworder[wli] = contparticles_all(topo,conf,sim,wli);
					break;
				default: 
					sim->wl.neworder[wli] = sim->wl.currorder[wli];
					break;
			}
		if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
		}
		if (!reject) { 
			wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
			energy += wlener;
		}
	}
	
	if (!reject) {
		enermove = conf->particle[target].delta_mu * pmone;
		// DEBUG
		//double dmu = enermove;
		//particle[target].switched += pmone;
		enermove += calc_energy( target, intfce, 1, topo, conf, sim,0);
		//printf("energy: %lf \t %lf\t%lf\n",particle[target].delta_mu, dmu, enermove);
	}
	
	// If not accepted: switch back
	if ( reject || movetry(energy,enermove,sim->temper) ) {  /* probability acceptance */
		DEBUG_SIM("Did NOT switch it\n");
		conf->particle[target].switchtype = conf->particle[target].type; 
		conf->particle[target].type = tmp_type;
		conf->particle[target].switched -= pmone;
		int_partvec(target,&(topo->ia_params[conf->particle[target].type][conf->particle[target].type]),conf);
		wlreject(sim,radiusholemax_orig);
	} else { /* move was accepted */
		wlaccept(sim->wlm[0],&sim->wl);
		edriftchanges = enermove - energy + wlener;
	}
	
	return edriftchanges;
}

/*.................................CHAIN MOVES....................................*/
/*................................................................................*/
double chainmove(struct topo * topo, struct sim * sim, struct conf * conf,
			   double (* intfce[MAXT][MAXT])(struct interacts *))
{
	double edriftchanges =0.0;
	long target;
	
	double ran2(long *);
	double chaindisplace(struct topo * topo, struct sim * sim, struct conf * conf,
						 double (* intfce[MAXT][MAXT])(struct interacts *), long target);
	double chainrotate(struct topo * topo, struct sim * sim, struct conf * conf,
						 double (* intfce[MAXT][MAXT])(struct interacts *), long target);
	
	/*=== This is a chain move step ===*/
	target = ran2(&seed) * topo->chainnum;
	if (ran2(&seed) < 0.5) {
		/*=== Displacement step of cluster/chain ===*/
		edriftchanges = chaindisplace(topo,sim,conf,intfce,target);
	} else {
		/*=== Rotation step of cluster/chain ===*/
		edriftchanges = chainrotate(topo,sim,conf,intfce,target);
		
	} /* ==== END OF CHAIN MOVES ===== */
	return edriftchanges;
}
/*................................................................................*/
double chaindisplace(struct topo * topo, struct sim * sim, struct conf * conf,
					 double (* intfce[MAXT][MAXT])(struct interacts *), long target)
{
	double edriftchanges,energy,enermove,wlener;
	struct vector dr, origsyscm;
	int reject=0,wli;
	struct vector cluscm;
	long current,i;
	struct particles chorig[MAXCHL];
	double radiusholemax_orig=0;
	
	double calc_energy(long target, double (* intfce[MAXT][MAXT])(struct interacts *), 
					   int mode, struct topo * topo, struct conf * conf, struct sim * sim,int chainn);
	int movetry(double, double, double);
	void wlreject(struct sim *,long);
	void wlaccept(int, struct wls *);
	long meshorder_movechain(long [MAXN], struct meshs *, long, struct conf * conf, \
							 struct sim * sim, struct particles chorig[MAXCHL],int);
	int mesh_cpy(struct meshs *, struct meshs *);
	//void mesh_print (struct meshs *);
	long z_order(struct wls *, struct conf * conf,int);
	long twopartdist(struct wls *, struct conf *conf,int);
	struct vector ranvec(void);
	int longarray_cpy (long **target, long **source,long,long);
	long radiusholeorder_movechain(long chain[MAXN], struct conf * conf, \
					struct sim * sim,struct particles chorig[MAXCHL],int,struct vector *);
	long radiushole_all(struct topo *topo, struct conf *conf, struct sim * sim,int, struct vector *);
	long contparticles_movechain(long chain[MAXN], struct conf * conf, struct sim * sim,struct particles chorig[MAXCHL],int wli);
	long contparticles_all(struct topo *topo, struct conf *conf, struct sim * sim,int wli);
	
	/*=== Displacement step of cluster/chain ===*/
	//printf ("move chain\n\n");
	energy =0.0;
	wlener = 0.0;
	edriftchanges=0.0;
	i=0;
	current = topo->chainlist[target][0];
	cluscm.x = 0;
	cluscm.y = 0;
	cluscm.z = 0;
	origsyscm.x = 0;
	origsyscm.y = 0;
	origsyscm.z = 0;
	while (current >=0 ) {   /* store old configuration calculate energy*/
		chorig[i].pos = conf->particle[current].pos;
		energy += calc_energy(current, intfce, 2, topo, conf, sim, target);
		i++;
		current = topo->chainlist[target][i];
	}
	dr = ranvec();
	dr.x *= sim->chainm[conf->particle[target].chaint].mx/conf->box.x;
	dr.y *= sim->chainm[conf->particle[target].chaint].mx/conf->box.y;
	dr.z *= sim->chainm[conf->particle[target].chaint].mx/conf->box.z;
	i=0;
	current = topo->chainlist[target][0];
	while (current >=0 ) { /* move chaine to new position  */
		if ( (sim->wlm[0] == 1) || (sim->wlm[0] == 5) || (sim->wlm[1] == 1) || (sim->wlm[1] == 5) ) { /* calculate move of center of mass  */
			cluscm.x += dr.x*topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
			cluscm.y += dr.y*topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
			cluscm.z += dr.z*topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
		}
		conf->particle[current].pos.x += dr.x;
		conf->particle[current].pos.y += dr.y;
		conf->particle[current].pos.z += dr.z;
		i++;
		current = topo->chainlist[target][i];
	}
	enermove = 0.0;
	
	reject = 0;
	if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
		for (wli=0;wli<sim->wl.wlmdim;wli++) {
			switch (sim->wlm[wli]) {
				case 1: origsyscm = conf->syscm;
					conf->syscm.x += cluscm.x / conf->sysvolume;
					conf->syscm.y += cluscm.y / conf->sysvolume;
					conf->syscm.z += cluscm.z / conf->sysvolume;
					sim->wl.neworder[wli] = z_order(&sim->wl, conf,wli);
					break;
				case 2: mesh_cpy(&sim->wl.origmesh,&sim->wl.mesh);
					sim->wl.neworder[wli] = meshorder_movechain(topo->chainlist[target], &sim->wl.mesh, topo->npart, conf, sim, chorig,wli);
					break;
				case 4:
					sim->wl.neworder[wli] = twopartdist(&sim->wl,conf,wli);
					break;
				case 5: 
					radiusholemax_orig = sim->wl.radiusholemax;
					origsyscm = conf->syscm;
					conf->syscm.x += cluscm.x / conf->sysvolume;
					conf->syscm.y += cluscm.y / conf->sysvolume;
					conf->syscm.z += cluscm.z / conf->sysvolume;
					longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
					sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->syscm));
					break;
				case 6: 
					radiusholemax_orig = sim->wl.radiusholemax;
					longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
					if ( target == 0 )
					  sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
					else 
					  sim->wl.neworder[wli] = radiusholeorder_movechain(topo->chainlist[target], conf, sim, chorig,wli,&(conf->particle[0].pos));
					break;
				case 7:
					sim->wl.partincontactold = sim->wl.partincontact;
					if ( target == 0 )
						sim->wl.neworder[wli] = contparticles_all(topo,conf,sim,wli);
					else 
						sim->wl.neworder[wli] = contparticles_movechain(topo->chainlist[target],conf,sim,chorig,wli);
					break;
				default: 
					sim->wl.neworder[wli] = sim->wl.currorder[wli];
					break;
			}
			if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
		}
		if (!reject) { 
			wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
			energy += wlener;
		}
	}
	if (!reject) { /* wang-landaou ok, try move - calcualte energy */
		i=0;
		current = topo->chainlist[target][0];
		while (current >=0 ) {
			enermove += calc_energy(current, intfce, 2, topo, conf, sim,target);
			i++;
			current = topo->chainlist[target][i];
		}
	}
	if ( reject || movetry(energy, enermove, sim->temper) ) {  /* probability acceptance */
		i=0;
		current = topo->chainlist[target][0];
		while (current >=0 ) {
			conf->particle[current].pos = chorig[i].pos;
			i++;
			current = topo->chainlist[target][i];
		}
		sim->chainm[conf->particle[target].chaint].rej++;
		if ( (sim->wlm[0] == 1) || (sim->wlm[0] == 5) || (sim->wlm[1] == 1) || (sim->wlm[1] == 5) ) 
			conf->syscm = origsyscm;
		wlreject(sim,radiusholemax_orig);
	} else { /* move was accepted */
		sim->chainm[conf->particle[target].chaint].acc++;
		wlaccept(sim->wlm[0],&sim->wl);
		edriftchanges = enermove - energy + wlener;
	}
	
	return edriftchanges;
}
/*................................................................................*/
double chainrotate(struct topo * topo, struct sim * sim, struct conf * conf,
				   double (* intfce[MAXT][MAXT])(struct interacts *), long target)
{
	double edriftchanges,energy,enermove,wlener;
	int reject=0,wli;
	struct vector cluscm;
	double chainvolume;
	long current, i;
	struct particles chorig[MAXCHL];
	double radiusholemax_orig=0;
	
	double calc_energy(long target, double (* intfce[MAXT][MAXT])(struct interacts *), 
					   int mode, struct topo * topo, struct conf * conf, struct sim * sim,int chainn);
	int movetry(double, double, double);
	void wlreject(struct sim *,long);
	void wlaccept(int, struct wls *);
	long meshorder_movechain(long [MAXN], struct meshs *, long, struct conf * conf, \
							 struct sim * sim, struct particles chorig[MAXCHL],int);
	int mesh_cpy(struct meshs *, struct meshs *);
	void cluster_rotate(long, struct vector, double, struct topo * topo, struct conf * conf);
	long z_order(struct wls *, struct conf * conf,int);
	long twopartdist(struct wls *, struct conf *conf,int);
	int longarray_cpy (long **target, long **source,long,long);
	long radiusholeorder_movechain(long chain[MAXN], struct conf * conf, struct sim * sim,\
					struct particles chorig[MAXCHL],int,struct vector *);
	long radiushole_all(struct topo *topo, struct conf *conf, struct sim * sim,int, struct vector *);
	long contparticles_movechain(long chain[MAXN], struct conf * conf, struct sim * sim,struct particles chorig[MAXCHL],int wli);
	long contparticles_all(struct topo *topo, struct conf *conf, struct sim * sim,int wli);
	
	/*=== Rotation step of cluster/chain ===*/
	//printf ("rotation of chain\n\n");
	energy=0.0; /* set values to zero*/
	edriftchanges=0.0;
	wlener = 0.0;
	current = topo->chainlist[target][0];
	cluscm.x = conf->particle[current].pos.x*topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
	cluscm.y = conf->particle[current].pos.y*topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
	cluscm.z = conf->particle[current].pos.z*topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
	chorig[0] = conf->particle[current];
	chainvolume = topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
	energy += calc_energy(current, intfce, 2, topo, conf, sim,target);
	i=1;
	current = topo->chainlist[target][i];
	while (current >=0 ) {   /* store old configuration calculate energy*/
		chorig[i] = conf->particle[current];
		/*We have chains whole! don't have to do PBC*/
		/*r_cm.x = conf->particle[current].pos.x - conf->particle[first].pos.x;
		 r_cm.y = conf->particle[current].pos.y - conf->particle[first].pos.y;
		 r_cm.z = conf->particle[current].pos.z - conf->particle[first].pos.z;
		 if ( r_cm.x < 0  ) 
		 r_cm.x -= (double)( (long)(r_cm.x-0.5) );
		 else 
		 r_cm.x -= (double)( (long)(r_cm.x+0.5) );
		 if ( r_cm.y < 0  ) 
		 r_cm.y -= (double)( (long)(r_cm.y-0.5) );
		 else 
		 r_cm.y -= (double)( (long)(r_cm.y+0.5) );
		 if ( r_cm.z < 0  ) 
		 r_cm.z -= (double)( (long)(r_cm.z-0.5) );
		 else 
		 r_cm.z -= (double)( (long)(r_cm.z+0.5) );
		 */
		cluscm.x += conf->particle[current].pos.x*topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
		cluscm.y += conf->particle[current].pos.y*topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
		cluscm.z += conf->particle[current].pos.z*topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
		chainvolume += topo->ia_params[conf->particle[current].type][conf->particle[current].type].volume;
		energy += calc_energy(current, intfce, 2, topo, conf, sim,target);
		i++;
		current = topo->chainlist[target][i];
	}
	cluscm.x = cluscm.x/chainvolume;
	cluscm.y = cluscm.y/chainvolume;
	cluscm.z = cluscm.z/chainvolume;
	/*do actual rotations around geometrical center*/
	cluster_rotate(target, cluscm, sim->chainr[conf->particle[target].chaint].angle, topo, conf);
	enermove=0.0;
	
	reject = 0;
	if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
		for (wli=0;wli<sim->wl.wlmdim;wli++) {
			switch (sim->wlm[wli]) {
				case 1: 
					if (target == 0) sim->wl.neworder[wli] = z_order(&sim->wl, conf,wli);
					else sim->wl.neworder[wli] = sim->wl.currorder[wli];
					/* if we rotated cluster it is around its CM so no change*/
					break;
				case 2: 
					mesh_cpy(&sim->wl.origmesh,&sim->wl.mesh);
					sim->wl.neworder[wli] = meshorder_movechain(topo->chainlist[target], &sim->wl.mesh, topo->npart, conf, sim, chorig,wli);
					break;
				case 3: 
					if (target == 0)  sim->wl.neworder[wli] = (long) floor( (conf->particle[0].dir.z - sim->wl.minorder[wli])/ sim->wl.dorder[wli] );
					else sim->wl.neworder[wli] = sim->wl.currorder[wli];
					/* only rotation change direction */
					break;
				case 4:
					sim->wl.neworder[wli] = twopartdist(&sim->wl,conf,wli);
					break;	
				case 5:
					radiusholemax_orig = sim->wl.radiusholemax;
					longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
					sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->syscm));
					break;
				case 6:
					radiusholemax_orig = sim->wl.radiusholemax;
					longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
					if ( target == 0 )
					    sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
					else 
					    sim->wl.neworder[wli] = radiusholeorder_movechain(topo->chainlist[target], conf, sim, chorig,wli,&(conf->particle[0].pos));
					break;
				case 7:
					sim->wl.partincontactold = sim->wl.partincontact;
					if ( target == 0 )
						sim->wl.neworder[wli] = contparticles_all(topo,conf,sim,wli);
					else 
						sim->wl.neworder[wli] = contparticles_movechain(topo->chainlist[target],conf,sim,chorig,wli);
					break;
				default: 
					sim->wl.neworder[wli] = sim->wl.currorder[wli];
					break;
			}
			if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
		}
		if (!reject) { 
			wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
			energy += wlener;
		}
	}
	if (!reject) { /* wang-landaou ok, try move - calcualte energy */
		i=0;
		current = topo->chainlist[target][0];
		while (current >=0 ) {
			enermove +=  calc_energy(current, intfce, 2, topo, conf, sim,target);
			i++;
			current = topo->chainlist[target][i];
		}
	}
	if ( reject || movetry(energy, enermove, sim->temper) ) { /* probability acceptance */
		i=0;
		current = topo->chainlist[target][0];
		while (current >=0 ) {
			conf->particle[current] = chorig[i];
			i++;
			current = topo->chainlist[target][i];
		}
		sim->chainr[conf->particle[target].chaint].rej++;
		wlreject(sim,radiusholemax_orig);
	} else { /* move was accepted */
		sim->chainr[conf->particle[target].chaint].acc++;
		wlaccept(sim->wlm[0],&sim->wl);
		edriftchanges = enermove - energy + wlener;
	}
	
	return edriftchanges;
}

/*..............................PRESSURE MOVES....................................*/
/*................................................................................*/
double pressuremove(struct topo * topo, struct sim * sim, struct conf * conf,
					double (* intfce[MAXT][MAXT])(struct interacts *))
{
	double edriftchanges,energy,enermove,wlener;
	int reject=0,wli;
	double old_side;   /* Box length before attempted change */
	double *side;      /* Box dimension to try changing */
	double psch;       /* Size of a box change during pressure */
	double pvol;       /* Size of a volume during pressure */
	double pvoln;      /* Size of a new volume during pressure */
	double rsave;      /* Saved random number */
	double area;  
	double radiusholemax_orig=0;

	
	double ran2(long *);
	double calc_energy(long target, double (* intfce[MAXT][MAXT])(struct interacts *), 
					   int mode, struct topo * topo, struct conf * conf, struct sim * sim,int chainn);
	int movetry(double, double, double);
	void wlreject(struct sim *,long);
	void wlaccept(int, struct wls *);	
	int mesh_init(struct meshs *, double, long, struct conf * conf, struct sim * sim);
	int mesh_cpy(struct meshs *, struct meshs *);
	long z_order(struct wls *, struct conf * conf,int);
	long twopartdist(struct wls *, struct conf *conf,int);
	long radiushole_all(struct topo *topo, struct conf *conf, struct sim * sim,int, struct vector *);
	int longarray_cpy (long **target, long **source,long,long);
	long contparticles_all(struct topo *topo, struct conf *conf, struct sim * sim,int wli);
	
	/*=== This is a volume change step ===*/
	/*calculate energy*/
	edriftchanges=0.0;
	wlener = 0.0;
	energy = calc_energy(0, intfce, 0, topo, conf, sim,0);
	/* Choose an edge */
	switch (sim->ptype) {
		case 0:
			/* Anisotropic pressure coupling */
			rsave = ran2(&seed);   
			if (rsave < 1.0/3.0) {
				side = &(conf->box.x);
				area = conf->box.y * conf->box.z;
			} else if (rsave < 2.0/3.0) {
				side = &(conf->box.y);
				area = conf->box.x * conf->box.z;
			} else {
				side = &(conf->box.z);
				area = conf->box.x * conf->box.y;
			}
			old_side = *side;
			*side += sim->edge.mx * (ran2(&seed) - 0.5);
			
			reject = 0;
			if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
				for (wli=0;wli<sim->wl.wlmdim;wli++) {
					switch (sim->wlm[wli]) {
						case 1: 
							sim->wl.neworder[wli] = z_order(&sim->wl, conf,wli);
							break;
						case 2: 
							mesh_cpy(&sim->wl.origmesh,&sim->wl.mesh);
							sim->wl.neworder[wli] = (long) (mesh_init(&sim->wl.mesh,sim->wl.wl_meshsize, topo->npart, conf, sim) - sim->wl.minorder[wli]);
							break;
						case 4:
							sim->wl.neworder[wli] = twopartdist(&sim->wl,conf,wli);
							break;
						case 5:
							radiusholemax_orig = sim->wl.radiusholemax;
							longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
							sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->syscm));
							break;
						case 6:
							radiusholemax_orig = sim->wl.radiusholemax;
							longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
							sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
							break;
						case 7:
							sim->wl.partincontactold = sim->wl.partincontact;
							sim->wl.neworder[wli] = contparticles_all(topo,conf,sim,wli);
							break;
						default: 
							sim->wl.neworder[wli] = sim->wl.currorder[wli];
							break;
					}
					if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
				}
				if (!reject) { 
					wlener = sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
					energy += wlener;
				}
			}
			if (!reject) { /* wang-landaou ok, try move - calculate energy */
				enermove = sim->press * area * (*side - old_side) - (double)topo->npart * log(*side/old_side) / sim->temper;
				enermove += calc_energy(0, intfce, 0, topo, conf, sim,0);
			}
			if ( reject || *side <= 0.0 || ( movetry(energy,enermove,sim->temper) ) ) { /* probability acceptance */
				*side = old_side;
				sim->edge.rej++;
				wlreject(sim,radiusholemax_orig);
			} else {  /* move was accepted */
				sim->edge.acc++;
				wlaccept(sim->wlm[0],&sim->wl);
				edriftchanges = enermove - energy + wlener;
			}
			break;
		case 1:
			/* Isotropic pressure coupling */
			psch = sim->edge.mx * (ran2(&seed) - 0.5);
			pvol = conf->box.x * conf->box.y * conf->box.z;
			conf->box.x += psch;
			conf->box.y += psch;
			conf->box.z += psch;
			pvoln = conf->box.x * conf->box.y * conf->box.z;
			
			reject = 0;
			if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
				for (wli=0;wli<sim->wl.wlmdim;wli++) {
					switch (sim->wlm[wli]) {
						case 1: sim->wl.neworder[wli] = z_order(&sim->wl,conf,wli);
							break;
						case 2: mesh_cpy(&sim->wl.origmesh,&sim->wl.mesh);
							sim->wl.neworder[wli] = (long) (mesh_init(&sim->wl.mesh,sim->wl.wl_meshsize, topo->npart, conf, sim) - sim->wl.minorder[wli]);
							break;
						case 4:
							sim->wl.neworder[wli] = twopartdist(&sim->wl,conf,wli);
							break;
						case 5:
							radiusholemax_orig = sim->wl.radiusholemax;
							longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
							sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->syscm));
							break;
						case 6:
							radiusholemax_orig = sim->wl.radiusholemax;
							longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
							sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
							break;
						case 7:
							sim->wl.partincontactold = sim->wl.partincontact;
							sim->wl.neworder[wli] = contparticles_all(topo,conf,sim,wli);
							break;
						default: 
							sim->wl.neworder[wli] = sim->wl.currorder[wli];
							break;
					}
					if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
				}
				if (!reject) { 
					wlener = sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
					energy += wlener;
				}
			}
			if (!reject) { /* wang-landaou ok, try move - calcualte energy */
				enermove = sim->press * (pvoln - pvol) - (double)topo->npart * log(pvoln/pvol) / sim->temper;
				enermove += calc_energy(0, intfce, 0, topo, conf, sim,0);
			}
			if ( reject || movetry(energy,enermove,sim->temper) )  { /* probability acceptance */
				conf->box.x -= psch;
				conf->box.y -= psch;
				conf->box.z -= psch;
				sim->edge.rej++;
				wlreject(sim,radiusholemax_orig);
			} else { /* move was accepted */
				sim->edge.acc++;
				wlaccept(sim->wlm[0],&sim->wl);
				edriftchanges = enermove - energy + wlener;
			}
			break;
		case 2:
			/* Isotropic pressure coupling in xy, z constant */
			psch = sim->edge.mx * (ran2(&seed) - 0.5);
			pvol = conf->box.x * conf->box.y;
			conf->box.x += psch;
			conf->box.y += psch;
			pvoln = conf->box.x * conf->box.y;
			
			reject = 0;
			if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
				for (wli=0;wli<sim->wl.wlmdim;wli++) {
					switch (sim->wlm[wli]) {
						/*no change in case 1, it does not change box.z*/
						case 2: mesh_cpy(&sim->wl.origmesh,&sim->wl.mesh);
							sim->wl.neworder[wli] = (long) (mesh_init(&sim->wl.mesh,sim->wl.wl_meshsize,topo->npart, conf, sim) - sim->wl.minorder[wli]);
							break;
						case 4:
							sim->wl.neworder[wli] = twopartdist(&sim->wl,conf,wli);
							break;
						case 5:
							radiusholemax_orig = sim->wl.radiusholemax;
							longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
							sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->syscm));
							break;
						case 6:
							radiusholemax_orig = sim->wl.radiusholemax;
							longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
							sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
							break;
						case 7:
							sim->wl.partincontactold = sim->wl.partincontact;
							sim->wl.neworder[wli] = contparticles_all(topo,conf,sim,wli);
							break;
						default: 
							sim->wl.neworder[wli] = sim->wl.currorder[wli];
							break;
					}
					if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
				}
				if (!reject) { 
					wlener = sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
					energy += wlener;
				}
			}
			if (!reject) { /* wang-landaou ok, try move - calculate energy */
				enermove = sim->press * conf->box.z * (pvoln - pvol) - (double)topo->npart * log(pvoln/pvol) / sim->temper;
				enermove += calc_energy(0, intfce, 0, topo, conf, sim,0);
			}
			if ( reject || movetry(energy,enermove,sim->temper) )  { /* probability acceptance */
				conf->box.x -= psch;
				conf->box.y -= psch;
				sim->edge.rej++;
				wlreject(sim,radiusholemax_orig);
			} else { /* move was accepted */
				sim->edge.acc++;
				wlaccept(sim->wlm[0],&sim->wl);
				edriftchanges = enermove - energy + wlener;
			}
			break;
		case 3:
			/* Isotropic pressure coupling in xy, z coupled to have fixed volume */
			psch = sim->edge.mx * (ran2(&seed) - 0.5);
			pvol = conf->box.x * conf->box.y * conf->box.z;
			conf->box.x += psch;
			conf->box.y += psch;
			conf->box.z = pvol / conf->box.x / conf->box.y;
			
			reject = 0;
			if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
				for (wli=0;wli<sim->wl.wlmdim;wli++) {
					switch (sim->wlm[wli]) {
						case 1: sim->wl.neworder[wli] = z_order(&sim->wl, conf,wli);
							break;
						case 2: mesh_cpy(&sim->wl.origmesh,&sim->wl.mesh);
							sim->wl.neworder[wli] = (long) (mesh_init(&sim->wl.mesh,sim->wl.wl_meshsize,topo->npart, conf, sim) - sim->wl.minorder[wli]);
							break;
						case 4:
							sim->wl.neworder[wli] = twopartdist(&sim->wl,conf,wli);
							break;
						case 5:
							radiusholemax_orig = sim->wl.radiusholemax;
							longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
							sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->syscm));
							break;
						case 6:
							radiusholemax_orig = sim->wl.radiusholemax;
							longarray_cpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
							sim->wl.neworder[wli] = radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
							break;
						case 7:
							sim->wl.partincontactold = sim->wl.partincontact;
							sim->wl.neworder[wli] = contparticles_all(topo,conf,sim,wli);
							break;
						default: 
							sim->wl.neworder[wli] = sim->wl.currorder[wli];
							break;
					}
					if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
				}
				if (!reject) { 
					wlener = sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
					energy += wlener;
				}
			}
			if (!reject) { /* wang-landaou ok, try move - calculate energy */
				enermove = calc_energy(0, intfce, 0, topo, conf, sim,0);
			}
			if ( reject || movetry(energy,enermove,sim->temper) )  { /* probability acceptance */
				conf->box.x -= psch;
				conf->box.y -= psch;
				conf->box.z = pvol / conf->box.x / conf->box.y;
				sim->edge.rej++;
				wlreject(sim,radiusholemax_orig);
			} else { /* move was accepted */
				sim->edge.acc++;
				wlaccept(sim->wlm[0],&sim->wl);
				edriftchanges = enermove - energy + wlener;
			}
			break;
			
		default:
			fprintf (stderr, "ERROR: unknown type of pressure coupling %d",sim->ptype);
			exit(1);
	}
	
	/*=== End volume change step ===*/
	return edriftchanges;
}




/*..................... Switch replicas move in  MPI ..............................*/
/*.................................................................................*/
double replicaexchangemove(struct topo * topo, struct sim * sim, struct conf * conf,
					  double (* intfce[MAXT][MAXT])(struct interacts *), long sweep )
{
	double edriftchanges=0.0;
#ifdef MPI	
	double change, *recwlweights;
	MPI_Status status;
	int oddoreven,count,wli,sizewl = 0;
	struct mpiexchangedata localmpi,receivedmpi;
	BOOL reject;
	long localwl,receivedwl;
	
	double ran2(long *);
	void gen_pairlist(struct topo * topo, struct sim * sim, struct conf * conf);
	int longarray_cpy (long **target, long **source,long,long);
	int mesh_init(struct meshs *, double, long, struct conf * conf, struct sim * sim);
	double calc_energy(long target, double (* intfce[MAXT][MAXT])(struct interacts *), 
					   int mode, struct topo * topo, struct conf * conf, struct sim * sim,int chainn);
	void wlaccept(int, struct wls *);
	//int mpi_newdatatypes();
	
	//mpi_newdatatypes();
	int i;
	struct vector vec;
	struct particles part;
	struct mpiexchangedata exch;
	MPI_Aint     dispstart;
					
	MPI_Datatype MPI_vector; 
	MPI_Datatype type[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; 
	int          blocklen[3] = {1, 1, 1}; 
	MPI_Aint     disp[3]; 
	MPI_Address( &vec, &dispstart);
	MPI_Address( &(vec.x), &disp[0]); 
	MPI_Address( &(vec.y), &disp[1]); 
	MPI_Address( &(vec.z), &disp[2]);  
	for (i=0; i <3; i++) disp[i] -= dispstart;	
	MPI_Type_struct( 3, blocklen, disp, type, &MPI_vector);
	MPI_Type_commit( &MPI_vector); 
	
	MPI_Datatype MPI_Particle; 
	MPI_Datatype type2[11] = {MPI_vector,MPI_vector,MPI_vector,MPI_vector,MPI_vector, MPI_LONG, MPI_LONG, MPI_INT,MPI_INT,MPI_DOUBLE, MPI_INT}; 
	int          blocklen2[11] = {1, 1, 2,4,2,1,1,1,1,1,1,}; 
	MPI_Aint     disp2[11]; 
	MPI_Address( &part, &dispstart);
	MPI_Address( &(part.pos), &disp2[0]); 
	MPI_Address( &(part.dir), &disp2[1]); 
	MPI_Address( &(part.patchdir), &disp2[2]); 
	MPI_Address( &(part.patchsides), &disp2[3]);
	MPI_Address( &(part.chdir), &disp2[4]);
	MPI_Address( &(part.chaint), &disp2[5]);
	MPI_Address( &(part.chainn), &disp2[6]);
	MPI_Address( &(part.type), &disp2[7]);
	MPI_Address( &(part.switchtype), &disp2[8]);
	MPI_Address( &(part.delta_mu), &disp2[9]);
	MPI_Address( &(part.switched), &disp2[10]);
	for (i=0; i <11; i++) disp2[i] -= dispstart;
	MPI_Type_struct( 11, blocklen2, disp2, type2, &MPI_Particle);
	MPI_Type_commit( &MPI_Particle); 
	
	if (sim->wl.length[1] > 0) {
	    sizewl = sim->wl.length[1] * sim->wl.length[0];
	} else {
	    sizewl = sim->wl.length[0];
	}
	MPI_Datatype MPI_exchange; 
	MPI_Datatype type3[7] = {MPI_vector, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_vector, MPI_LONG, MPI_LONG}; 
	int          blocklen3[7] = {1, 1, 1, 1, 1, 1, 2}; 
	MPI_Aint     disp3[7]; 
	MPI_Address( &exch, &dispstart);
	MPI_Address( &(exch.box), &disp3[0]); 
	MPI_Address( &(exch.energy), &disp3[1]); 
	MPI_Address( &(exch.volume), &disp3[2]); 
	MPI_Address( &(exch.accepted), &disp3[3]); 
	MPI_Address( &(exch.syscm), &disp3[4]); 
	MPI_Address( &(exch.radiusholemax), &disp3[5]); 
	MPI_Address( &(exch.wl_order), &disp3[6]); 
	for (i=0; i <7; i++) disp3[i] -= dispstart;	
	MPI_Type_struct(7, blocklen3, disp3, type3, &MPI_exchange);
	MPI_Type_commit( &MPI_exchange);
	/*=== This is an attempt to switch replicas ===*/

	localmpi.box = conf->box;
	localmpi.energy = calc_energy(0, intfce, 0, topo, conf, sim,0);
	localmpi.volume = conf->box.x * conf->box.y * conf->box.z;
	localmpi.accepted = 0;
	localmpi.syscm = conf->syscm;
	localmpi.radiusholemax = sim->wl.radiusholemax;
	recwlweights =  malloc( sizeof(double) * sizewl  );
	for (wli=0;wli<2;wli++) {
		localmpi.wl_order[wli] = 0;
		receivedmpi.wl_order[wli] = 0;
	}
	for (wli=0;wli<sim->wl.wlmdim;wli++) {
		localmpi.wl_order[wli] = sim->wl.currorder[wli];
		//fprintf(stdout,"wli %d %ld  %ld\n\n", wli, localmpi.wl_order[wli], sim->wl.currorder[wli] );
	}

	if ( (sweep % (2*sim->nrepchange)) == 0) 
		/* exchange odd ones with even ones*/
		oddoreven=1;
	else 
		/* exchange even ones with odd ones*/
		oddoreven=0;
	if (sim->mpinprocs == 2) 
		oddoreven=1;
	count = 1;
		
	if (sim->mpirank % 2 == oddoreven) {
		if (sim->mpirank > 0) {
			MPI_Send(&localmpi, 1, MPI_exchange, sim->mpirank-1, count, MPI_COMM_WORLD);
			MPI_Send(sim->wl.weights, sizewl, MPI_DOUBLE, sim->mpirank-1, count, MPI_COMM_WORLD);
			//printf("send data: rank: %d energy: %f volume: %f pressure: %f \n",sim->mpirank,localmpi.energy,localmpi.volume,localmpi.pressure);
			
			MPI_Recv(&receivedmpi, 1, MPI_exchange, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			/*decision of accepting or rejecting the exchange was done on other process
			here we took received configuration (if move was accepted))*/
			//printf("received data: rank: %d energy: %f volume: %f pressure: %f \n",sim->mpirank,receivedmpi.energy,receivedmpi.volume,receivedmpi.pressure);
			
			if (receivedmpi.accepted == 1) {
				sim->mpiexch.acc++;
				struct particles *temppart;
				temppart = malloc(topo->npart*sizeof(struct particles));
				MPI_Recv(temppart, topo->npart, MPI_Particle, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
/*				printf("received data: rank: %d\n", sim->mpirank);
				printf("part0  x %f y %f z %f\n",temppart[0].pos.x, temppart[0].pos.y, temppart[0].pos.z);
				printf("part1  x %f y %f z %f\n",temppart[1].pos.x, temppart[1].pos.y, temppart[1].pos.z);
				printf("part0  chaint %ld chainn %ld type %d\n",temppart[0].chaint,temppart[0].chainn,temppart[0].type);
*/
				MPI_Send(conf->particle, topo->npart, MPI_Particle, sim->mpirank-1, count, MPI_COMM_WORLD);
/*				printf("send data: rank: %d\n",sim->mpirank);
				printf("part0  x %f y %f z %f\n",conf->particle[0].pos.x,conf->particle[0].pos.y,conf->particle[0].pos.z);
				printf("part1  x %f y %f z %f\n",conf->particle[1].pos.x,conf->particle[1].pos.y,conf->particle[1].pos.z);
				printf("part0  chaint %ld chainn %ld type %d\n",conf->particle[0].chaint,conf->particle[0].chainn,conf->particle[0].type);
*/
				localmpi.accepted = receivedmpi.accepted;
				conf->box = receivedmpi.box;
				conf->syscm = receivedmpi.syscm;
				memcpy(conf->particle,temppart,topo->npart*sizeof(struct particles));
				edriftchanges = receivedmpi.energy - localmpi.energy;
				edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)topo->npart * log(receivedmpi.volume / localmpi.volume) / sim->temper;
				if ( sim->wlm[0] >0 ) {
					for (wli=0;wli<sim->wl.wlmdim;wli++) {
						sim->wl.neworder[wli] = receivedmpi.wl_order[wli];
					}
					wlaccept(sim->wlm[0],&sim->wl);
					//exchange wl data mesh size and radius hole s
					for (wli=0;wli<sim->wl.wlmdim;wli++) {
						switch (sim->wlm[wli]) {
							case 2:
								/*it is complicated to send because of different sizes 
								 we would have to send sizes first and realocate corrrect mesh size and then send data
								 it is better to recalculate (a bit slower though)*/
								mesh_init(&sim->wl.mesh,sim->wl.wl_meshsize, topo->npart, conf, sim);
								break;
							case 5:
								//radiushole_all(topo,conf,sim,wli,&(conf->syscm));
								sim->wl.radiusholeold = (long*) realloc(sim->wl.radiusholeold,sizeof(long)*receivedmpi.radiusholemax);
								MPI_Recv(sim->wl.radiusholeold,receivedmpi.radiusholemax, MPI_LONG, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
								MPI_Send(sim->wl.radiushole,sim->wl.radiusholemax, MPI_LONG, sim->mpirank-1, count, MPI_COMM_WORLD);
								longarray_cpy(&sim->wl.radiushole,&sim->wl.radiusholeold,sim->wl.radiusholemax,receivedmpi.radiusholemax);
								sim->wl.radiusholemax=receivedmpi.radiusholemax;
								break;
							case 6:
								//radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
								sim->wl.radiusholeold = (long*) realloc(sim->wl.radiusholeold,sizeof(long)*receivedmpi.radiusholemax);
								MPI_Recv(sim->wl.radiusholeold,receivedmpi.radiusholemax, MPI_LONG, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
								MPI_Send(sim->wl.radiushole,sim->wl.radiusholemax, MPI_LONG, sim->mpirank-1, count, MPI_COMM_WORLD);
								longarray_cpy(&sim->wl.radiushole,&sim->wl.radiusholeold,sim->wl.radiusholemax,receivedmpi.radiusholemax);
								sim->wl.radiusholemax=receivedmpi.radiusholemax;
								break;
							case 7:	
								//contparticles_all(topo,conf,sim,wli);
								MPI_Recv(&(sim->wl.partincontactold),1, MPI_LONG, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
								MPI_Send(&(sim->wl.partincontact),1, MPI_LONG, sim->mpirank-1, count, MPI_COMM_WORLD);
								sim->wl.partincontact=sim->wl.partincontactold;
								break;
						}
					}
				}
				
				free(temppart);
			} else {
				sim->mpiexch.rej++;
				if ( sim->wlm[0] > 0 ) {
					sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]] -= sim->wl.alpha;
					sim->wl.hist[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]]++;
				}
			}
			
		}
	} else {
		if (sim->mpirank+1 < sim->mpinprocs) {
			/*there is above process*/
			MPI_Recv(&receivedmpi, 1, MPI_exchange, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(recwlweights, sizewl, MPI_DOUBLE, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			/*we got new configuration*/
			//printf("received data: rank: %d energy: %f volume: %f \n",sim->mpirank,receivedmpi.energy,receivedmpi.volume);
					
			/*evaluate if accepte or reject the configuration*/
			/*acc = exp( (1/sim->temper - 1/(sim->temper + sim.dtemp)) * (E_here - E_received) + 
			(sim->press /sim->temper - pressure_received /(sim.temper + sim->dtemp)) * (V_here - V_received)
			if pressure the same it it simplier*/
			reject = FALSE;
			change = (1/sim->temper - 1/(sim->temper + sim->dtemp)) * (localmpi.energy - receivedmpi.energy);
			//printf("acceptance decision: change: %f localE: %f receivedE: %f tempf: %f \n",change,localmpi.energy,receivedmpi.energy,(1/sim->temper - 1/(sim->temper + sim->dtemp)));
			change += (sim->press/sim->temper - (sim->press + sim->dpress)/(sim->temper + sim->dtemp)) * (localmpi.volume - receivedmpi.volume);
			//printf("pressf: %f  \n",(sim->press/sim->temper - (sim->press + sim->dpress)/(sim->temper + sim->dtemp)));
			if (sim->wlm[0] > 0) {
				localwl = sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0];
				receivedwl = receivedmpi.wl_order[0] + receivedmpi.wl_order[1]*sim->wl.length[0];
				//fprintf(stdout,"decide wl   %ld %ld %ld energychange: %f \n", receivedmpi.wl_order[0],  receivedmpi.wl_order[1], receivedwl, change );
				//fprintf(stdout,"local weights %ld %f %ld %f \n",localwl,sim->wl.weights[localwl],receivedwl,sim->wl.weights[receivedwl]);
				change += (-sim->wl.weights[localwl] + sim->wl.weights[receivedwl] )/sim->temper + ( -recwlweights[receivedwl] + recwlweights[localwl])/(sim->temper + sim->dtemp) ;
				//fprintf(stdout,"wlchange %f \n\n",change);
			}
			if (  (!(reject)) && ( (change > 0) || (ran2(&seed) < exp(change))  )  ) {
				/* Exchange ACCEPTED send local stuff*/
				//printf("exchange accepted \n");
				sim->mpiexch.acc++;
				localmpi.accepted = 1;
				conf->box = receivedmpi.box;
				conf->syscm = receivedmpi.syscm;
				edriftchanges = receivedmpi.energy - localmpi.energy;
				edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)topo->npart * log(receivedmpi.volume / localmpi.volume) / sim->temper;
				//printf("edrift %f\n",edriftchanges);
				if ( sim->wlm[0] > 0 ) {
					for (wli=0;wli<sim->wl.wlmdim;wli++) {
						sim->wl.neworder[wli] = receivedmpi.wl_order[wli];
					}
					wlaccept(sim->wlm[0],&sim->wl);
				}
				MPI_Send(&localmpi, 1, MPI_exchange, sim->mpirank+1, count, MPI_COMM_WORLD);
				//printf("send data: rank: %d energy: %f volume: %f pressure: %f \n",sim->mpirank,localmpi.energy,localmpi.volume,localmpi.pressure);
				/*send and receive configuration*/
				MPI_Send(conf->particle, topo->npart, MPI_Particle, sim->mpirank+1, count, MPI_COMM_WORLD);
/*				printf("send data: rank: %d\n",sim->mpirank);
				printf("part0  x %f y %f z %f\n",conf->particle[0].pos.x,conf->particle[0].pos.y,conf->particle[0].pos.z);
				printf("part1  x %f y %f z %f\n",conf->particle[1].pos.x,conf->particle[1].pos.y,conf->particle[1].pos.z);
				printf("part0  chaint %ld chainn %ld type %d\n",conf->particle[0].chaint,conf->particle[0].chainn,conf->particle[0].type);
*/
				MPI_Recv(conf->particle, topo->npart, MPI_Particle, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
/*				printf("recieved data: rank: %d\n",sim->mpirank);
				printf("part0  x %f y %f z %f\n",conf->particle[0].pos.x,conf->particle[0].pos.y,conf->particle[0].pos.z);
				printf("part1  x %f y %f z %f\n",conf->particle[1].pos.x,conf->particle[1].pos.y,conf->particle[1].pos.z);
				printf("part0  chaint %ld chainn %ld type %d\n",conf->particle[0].chaint,conf->particle[0].chainn,conf->particle[0].type);
*/				
				if ( sim->wlm[0] > 0 ) {
					//exchange wl data mesh size and radius hole s
					for (wli=0;wli<sim->wl.wlmdim;wli++) {
						switch (sim->wlm[wli]) {
							case 2:
								/*it is complicated to send because of different sizes 
								  we would have to send sizes first and realocate corrrect mesh size and then send data
								  it is better to recalculate (a bit slower though)*/
								mesh_init(&sim->wl.mesh,sim->wl.wl_meshsize, topo->npart, conf, sim);
								break;
							case 5:
								//radiushole_all(topo,conf,sim,wli,&(conf->syscm));
								sim->wl.radiusholeold = (long*) realloc(sim->wl.radiusholeold,sizeof(long)*receivedmpi.radiusholemax);
								MPI_Send(sim->wl.radiushole,sim->wl.radiusholemax, MPI_LONG, sim->mpirank+1, count, MPI_COMM_WORLD);
								MPI_Recv(sim->wl.radiusholeold,receivedmpi.radiusholemax, MPI_LONG, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
								longarray_cpy(&sim->wl.radiushole,&sim->wl.radiusholeold,sim->wl.radiusholemax,receivedmpi.radiusholemax);
								sim->wl.radiusholemax=receivedmpi.radiusholemax;
								break;
							case 6:
								//radiushole_all(topo,conf,sim,wli,&(conf->particle[0].pos));
								sim->wl.radiusholeold = (long*) realloc(sim->wl.radiusholeold,sizeof(long)*receivedmpi.radiusholemax);
								MPI_Send(sim->wl.radiushole,sim->wl.radiusholemax, MPI_LONG, sim->mpirank+1, count, MPI_COMM_WORLD);
								MPI_Recv(sim->wl.radiusholeold,receivedmpi.radiusholemax, MPI_LONG, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
								longarray_cpy(&sim->wl.radiushole,&sim->wl.radiusholeold,sim->wl.radiusholemax,receivedmpi.radiusholemax);
								sim->wl.radiusholemax=receivedmpi.radiusholemax;
								break;
							case 7:	
								//contparticles_all(topo,conf,sim,wli);
								MPI_Send(&(sim->wl.partincontact),1, MPI_LONG, sim->mpirank+1, count, MPI_COMM_WORLD);
								MPI_Recv(&(sim->wl.partincontact),1, MPI_LONG, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
								break;
						}
					}
				}
			} else {
				/*if exchange rejected send back info */
				//printf("exchange rejected\n");
				sim->mpiexch.rej++;
				MPI_Send(&localmpi, 1, MPI_exchange, sim->mpirank+1, count, MPI_COMM_WORLD);
				if ( sim->wlm[0] > 0 ) {
					sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]] -= sim->wl.alpha;
					sim->wl.hist[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]]++;
				}
			}
		}
	}
	if ( (localmpi.accepted) && (sim->pairlist_update) )
		gen_pairlist(topo, sim, conf);
	MPI_Type_free(&MPI_exchange);
	MPI_Type_free(&MPI_Particle);
	MPI_Type_free(&MPI_vector);	
	free(recwlweights);
#endif	
	return edriftchanges;
}

/*int mpi_newdatatypes()
{	
	int i;
	struct vector vec;
	struct particles part;
	struct mpiexchangedata exch;
	MPI_Aint     dispstart;
	
	MPI_Datatype MPI_vector; 
	MPI_Datatype type[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; 
	int          blocklen[3] = {1, 1, 1}; 
	MPI_Aint     disp[3]; 
	MPI_Address( &vec, &dispstart);
	MPI_Address( &(vec.x), &disp[0]); 
	MPI_Address( &(vec.y), &disp[1]); 
	MPI_Address( &(vec.z), &disp[2]);  
	for (i=0; i <3; i++) disp[i] -= dispstart;	
	MPI_Type_struct( 3, blocklen, disp, type, &MPI_vector);
	MPI_Type_commit( &MPI_vector); 
	
	MPI_Datatype MPI_Particle; 
	MPI_Datatype type2[11] = {MPI_vector,MPI_vector,MPI_vector,MPI_vector,MPI_vector, MPI_LONG, MPI_LONG, MPI_INT,MPI_INT,MPI_DOUBLE, MPI_INT}; 
	int          blocklen2[11] = {1, 1, 2,4,2,1,1,1,1,1,1,}; 
	MPI_Aint     disp2[11]; 
	MPI_Address( &part, &dispstart);
	MPI_Address( &(part.pos), &disp2[0]); 
	MPI_Address( &(part.dir), &disp2[1]); 
	MPI_Address( &(part.patchdir), &disp2[2]); 
	MPI_Address( &(part.patchsides), &disp2[3]);
	MPI_Address( &(part.chdir), &disp2[4]);
	MPI_Address( &(part.chaint), &disp2[5]);
	MPI_Address( &(part.chainn), &disp2[6]);
	MPI_Address( &(part.type), &disp2[7]);
	MPI_Address( &(part.switchtype), &disp2[8]);
	MPI_Address( &(part.delta_mu), &disp2[9]);
	MPI_Address( &(part.switched), &disp2[10]);
	for (i=0; i <11; i++) disp2[i] -= dispstart;
	MPI_Type_struct( 11, blocklen2, disp2, type2, &MPI_Particle);
	MPI_Type_commit( &MPI_Particle); 
	
	MPI_Datatype MPI_exchange; 
	MPI_Datatype type3[5] = {MPI_vector, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT}; 
	int          blocklen3[5] = {1, 1, 1, 1, 1}; 
	MPI_Aint     disp3[5]; 
	MPI_Address( &exch, &dispstart);
	MPI_Address( &(exch.box), &disp3[0]); 
	MPI_Address( &(exch.energy), &disp3[1]); 
	MPI_Address( &(exch.volume), &disp3[2]); 
	MPI_Address( &(exch.pressure), &disp3[3]); 
	MPI_Address( &(exch.accepted), &disp3[4]); 
	for (i=0; i <5; i++) disp3[i] -= dispstart;	
	MPI_Type_struct( 5, blocklen3, disp3, type3, &MPI_exchange);
	MPI_Type_commit( &MPI_exchange);
	
	return 0;
}*/

/*................................................................................*/
/*................................................................................*/
/*....................END OF MOVES, INTERACTION FUNCTIONS FOLLOW..................*/
/*................................................................................*/


/*..............................................................................*/
/* 
   Determines total energy of two spherocylinders type PSC PSC 
 */
double e_psc_psc(struct interacts * interact)
{
	double atrenergy, repenergy;
	
	void closestdist(struct interacts *);
	double erepulsive(struct interacts *);
	double eattractive_psc_psc(struct interacts *,int,int);
	
	closestdist(interact);
	repenergy = erepulsive(interact);
	if ( ( interact->dist > interact->param->rcut ) || ( interact->param->epsilon == 0.0 ) ) 
		atrenergy = 0.0;
	else {
		BOOL firstCH=FALSE, secondCH=FALSE;
		struct vector olddir1 = interact->part1->dir;
		struct vector olddir2 = interact->part2->dir;
		if ( (interact->param->geotype[0] == CHPSC)||(interact->param->geotype[0] == TCHPSC) )
			firstCH = TRUE;
		if ( (interact->param->geotype[1] == CHPSC)||(interact->param->geotype[1] == TCHPSC) )
			secondCH = TRUE;
		if (firstCH)
			interact->part1->dir = interact->part1->chdir[0];
		if (secondCH)
			interact->part2->dir = interact->part2->chdir[0];
		
		if ((firstCH) || (secondCH) ) { 
		    closestdist(interact);
		}
		atrenergy = eattractive_psc_psc(interact,0,0);

		/*addition of interaction of second patches*/
		if ( (interact->param->geotype[0] == TPSC) || (interact->param->geotype[0] == TCHPSC) || 
		  (interact->param->geotype[1] == TPSC) ||(interact->param->geotype[1] == TCHPSC) ) {
		    BOOL firstT=FALSE, secondT=FALSE;
		    if ( (interact->param->geotype[0] == TPSC) || (interact->param->geotype[0] == TCHPSC) )
			firstT = TRUE;
		    if ( (interact->param->geotype[1] == TPSC) ||(interact->param->geotype[1] == TCHPSC)  )
			secondT = TRUE;
		    
		    if (firstT) {
                        if (firstCH && secondCH) {
                            interact->part1->dir = interact->part1->chdir[1];
                            interact->part2->dir = interact->part2->chdir[0];
                            closestdist(interact);
                        }
                        if (firstCH && !secondCH) {
                            interact->part1->dir = interact->part1->chdir[1];
                            closestdist(interact);
                        }
                        if (!firstCH && secondCH) {
                            interact->part2->dir = interact->part2->chdir[0];
                            closestdist(interact);
                        }
			atrenergy += eattractive_psc_psc(interact,1,0);
                    }

		    if ( (firstT) && (secondT) ) {
                        if (firstCH && secondCH) {
                            interact->part1->dir = interact->part1->chdir[1];
                            interact->part2->dir = interact->part2->chdir[1];
                            closestdist(interact);
                        }
                        if (firstCH && !secondCH) {
                            interact->part1->dir = interact->part1->chdir[1];
                            closestdist(interact);
                        }
                        if (!firstCH && secondCH) {
			    interact->part2->dir = interact->part2->chdir[1];
			    closestdist(interact);
			}
			atrenergy += eattractive_psc_psc(interact,1,1);
		    }

		    if (secondT) {
                        if (firstCH && secondCH) {
                            interact->part1->dir = interact->part1->chdir[0];
                            interact->part2->dir = interact->part2->chdir[1];
                            closestdist(interact);
                        }
                        if (firstCH && !secondCH) {
                            interact->part1->dir = interact->part1->chdir[0];
                            closestdist(interact);
                        }
                        if (!firstCH && secondCH) {
                            interact->part2->dir = interact->part2->chdir[0];
                            closestdist(interact);
                        }
			atrenergy += eattractive_psc_psc(interact,0,1);
		    }
		}
		
		if (firstCH) 
		    interact->part1->dir = olddir1;
		if (secondCH)
		    interact->part2->dir = olddir2;
		
	}
	return repenergy+atrenergy;
	
}

/* 
   Determines attractive energy of two spherocylinders type PSC PSC 
 */
double eattractive_psc_psc(struct interacts * interact,int patchnum1,int patchnum2)
{
	int i, intrs;
	double rcut, atrenergy, ndist;
	double v1, v2, f0, f1, f2, T1, T2, S1, S2, a;
	double intersections[5];
	struct vector vec1, vec2, vec_intrs, vec_mindist;

	struct vector vec_sub(struct vector, struct vector);
	struct vector vec_sum(struct vector, struct vector);
	struct vector vec_create(double, double, double);
	struct vector vec_scale(struct vector, double);
	struct vector vec_perpproject(struct vector *, struct vector*);
	struct vector mindist_segments(struct vector, double, struct vector, double, struct vector);
	void normalise(struct vector *);
	int psc_intersect(struct particles *, struct particles *,
			double, double, struct vector, double *,double, struct ia_param *, int which, int patchnum);
	double fanglscale(double, struct ia_param *, int which);

	rcut = interact->param->rcut;

	//interact->halfl = interact->param->half_len[0];
	//DEBUG_SIM("halfl = %lf", interact->halfl);
	for(i=0;i<5;i++) 
		intersections[i]=0;
	//cospatch = param.pcanglsw;
	//cospatchinr = param.pcangl;
	/*1- do intersections of spherocylinder2 with patch of spherocylinder1 at. 
	  cut distance C*/
	//DEBUG_SIM("first intersection");
	intrs=psc_intersect(interact->part1,interact->part2,interact->param->half_len[0],interact->param->half_len[1],interact->r_cm, intersections, rcut, interact->param,0, patchnum1);
	if (intrs <2){ 
		//DEBUG_SIM("No intersection :(");
		return 0.0; /*sc is all outside patch, attractive energy is 0*/
	}

	T1=intersections[0]; /*points on sc2*/
	T2=intersections[1];
	/*2- now do the same oposite way psc1 in patch of psc2*/
	for(i=0;i<5;i++) 
		intersections[i]=0;
	//DEBUG_SIM("get vector");
	vec1=vec_scale(interact->r_cm,-1.0);
	//DEBUG_SIM("second intersection");
	intrs=psc_intersect(interact->part2,interact->part1,interact->param->half_len[1],interact->param->half_len[0],vec1, intersections, rcut, interact->param,1, patchnum2);
	if (intrs <2) 
		return 0.0; /*sc is all outside patch, attractive energy is 0*/

	S1=intersections[0]; /*points on sc1*/
	S2=intersections[1];

	/*3- scaling function1: dependence on the length of intersetions*/
	v1=fabs(S1-S2);
	v2=fabs(T1-T2);
	f0=0.5*(v1+v2);
	
	/*4a- with two intersection pices calculate vector between their CM 
	  -this is for angular orientation*/
	vec1=vec_scale(interact->part1->dir,(S1+S2)*0.5);
	vec2=vec_scale(interact->part2->dir,(T1+T2)*0.5);
	vec_intrs.x=vec2.x-vec1.x-interact->r_cm.x;
	vec_intrs.y=vec2.y-vec1.y-interact->r_cm.y;
	vec_intrs.z=vec2.z-vec1.z-interact->r_cm.z;
	/*vec_intrs should be from sc1 to sc2*/
	//fprintf (stderr, "segments_CM: %.8f %.8f %.8f \n",vec_intrs.x,vec_intrs.y,vec_intrs.z);

	/*4b - calculate closest distance attractive energy from it*/
	vec_mindist = mindist_segments(interact->part1->dir,v1,interact->part2->dir,v2,vec_intrs);
	//fprintf (stderr, "segments closest dist: %.8f %.8f %.8f \n",vec_mindist.x,vec_mindist.y,vec_mindist.z);
	ndist=sqrt(DOT(vec_mindist,vec_mindist));
	//dist=DOT(vec_intrs,vec_intrs);
	if (ndist < interact->param->pdis) 
		atrenergy = -interact->param->epsilon;
	//atrenergy = -1.0;
	else {
		atrenergy = cos(PIH*(ndist-interact->param->pdis)/interact->param->pswitch);
		atrenergy *= -atrenergy*interact->param->epsilon ;
	}
	/*5- scaling function2: angular dependence of patch1*/
	vec1=vec_scale(vec_intrs,1.0);
	//vec1=vec_scale(vec_mindist,-1.0);
	vec1=vec_perpproject(&vec1, &interact->part1->dir);
	normalise(&vec1);
	a = DOT(vec1,interact->part1->patchdir[patchnum1]);
	f1 = fanglscale(a,interact->param, 0+2*patchnum1);

	/*6- scaling function3: angular dependence of patch2*/
	vec1=vec_scale(vec_intrs,-1.0);
	//vec1=vec_scale(vec_mindist,1.0);
	vec1=vec_perpproject(&vec1, &interact->part2->dir);
	normalise(&vec1);
	a = DOT(vec1,interact->part2->patchdir[patchnum2]);
	f2 = fanglscale(a,interact->param, 1+2*patchnum2);
	//printf("v1: %f v2: %f f0: %f f1: %f f2: %f ener: %f\n",v1,v2,f0,f1,f2,atrenergy);

	/*7- put it all together*/
	atrenergy *=f0*f1*f2;
	//if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
	//	    fprintf (stderr, "attraction  %.8f \n",atrenergy);
	//	    exit(1);

	return atrenergy;
}

/* a = r_ij * n_i */
double fanglscale(double a, struct ia_param * param, int which)
{
	double f;

	// TODO for different types
	if (a <= param->pcanglsw[which]) 
		f=0.0;
	else {
		if (a >= param->pcangl[which]) 
			f=1.0;
		else {
			f = 0.5 - ((param->pcanglsw[which] + param->pcangl[which])*0.5 - a )/(param->pcangl[which] - param->pcanglsw[which]);
		}
	}

	return f;
}


/*CPSC..............................................................................*/
/* 
   Determines total energy of two spherocylinders of type 3 -cylindrical psc -CPSC
 */
double e_cpsc_cpsc(struct interacts * interact)
{
	double atrenergy, repenergy;

	void closestdist(struct interacts *);
	double erepulsive(struct interacts *);
	double eattractive_cpsc_cpsc(struct interacts *,int,int);

	//DEBUG_SIM("do energy 33") ;
	closestdist(interact);
	repenergy = erepulsive(interact);
	//DEBUG_SIM("got the rep. energy");
	if ( ( interact->dist > interact->param->rcut ) || ( interact->param->epsilon == 0.0 ) ) 
		atrenergy = 0.0;
	else {
		BOOL firstCH=FALSE, secondCH=FALSE;
		struct vector olddir1 = interact->part1->dir;
		struct vector olddir2 = interact->part2->dir;
		if ( (interact->param->geotype[0] == CHCPSC)||(interact->param->geotype[0] == TCHCPSC) )
			firstCH = TRUE;
		if ( (interact->param->geotype[1] == CHCPSC)||(interact->param->geotype[1] == TCHCPSC) )
			secondCH = TRUE;
		if(firstCH)
			interact->part1->dir = interact->part1->chdir[0];
		if(secondCH)
			interact->part2->dir = interact->part2->chdir[0];
		
		if ((firstCH) || (secondCH) ) {
		    closestdist(interact);
		}
		atrenergy = eattractive_cpsc_cpsc(interact,0,0);
		
		/*addition of interaction of second patches*/
		if ( (interact->param->geotype[0] == TCPSC) || (interact->param->geotype[0] == TCHCPSC) || 
		  (interact->param->geotype[1] == TCPSC) ||(interact->param->geotype[1] == TCHCPSC) ) {
		    BOOL firstT=FALSE, secondT=FALSE;
		    if ( (interact->param->geotype[0] == TCPSC) || (interact->param->geotype[0] == TCHCPSC) )
			firstT = TRUE;
		    if ( (interact->param->geotype[1] == TCPSC) ||(interact->param->geotype[1] == TCHCPSC)  )
			secondT = TRUE;
		    
		    if (firstT) {
			if (firstCH) {
			    interact->part1->dir = interact->part1->chdir[1];
			    closestdist(interact);
			}
			atrenergy += eattractive_cpsc_cpsc(interact,1,0);
		    }
		    if ( (firstT) && (secondT) ) {
			if (secondCH) {
			    interact->part2->dir = interact->part2->chdir[1];
			    closestdist(interact);
			}
			atrenergy += eattractive_cpsc_cpsc(interact,1,1);
		    }
		    if (secondT) {
			if (firstT && firstCH ) {
			    interact->part1->dir = interact->part1->chdir[0];
			    closestdist(interact);
			}
			atrenergy += eattractive_cpsc_cpsc(interact,0,1);
		    }
		}
		
		if (firstCH) 
		    interact->part1->dir = olddir1;
		if (secondCH)
		    interact->part2->dir = olddir2;
	}
	return repenergy+atrenergy;
}

/* 
   Determines attractive energy of two spherocylinders of type 3 -cylindrical psc -CPSC
 */
double eattractive_cpsc_cpsc(struct interacts * interact, int patchnum1, int patchnum2)
{
	int i, intrs;
	double rcut, atrenergy, v1, v2, f0, f1, f2, T1, T2, S1, S2, a, ndist;
	double intersections[5];
	struct vector vec1, vec2, vec_intrs, vec_mindist;

	struct vector vec_sub(struct vector, struct vector);
	struct vector vec_sum(struct vector, struct vector);
	struct vector vec_create(double, double, double);
	struct vector vec_scale(struct vector, double);
	struct vector vec_perpproject(struct vector*, struct vector*);
	struct vector mindist_segments(struct vector, double, struct vector, double, struct vector);
	void normalise(struct vector *);
	int cpsc_intersect(struct particles *, struct particles *,
			double, double, struct vector, double *,double, struct ia_param *, int which, int patchnum);
	double fanglscale(double, struct ia_param *, int which);

	rcut = interact->param->rcut;

//	interact->halfl = interact->param->half_len[0];
	for(i=0;i<5;i++) 
		intersections[i]=0;
	/*1- do intersections of spherocylinder2 with patch of spherocylinder1 at. 
	  cut distance C*/
	intrs=cpsc_intersect(interact->part1,interact->part2,interact->param->half_len[0],interact->param->half_len[1],interact->r_cm, intersections, rcut, interact->param,0, patchnum1);
	if (intrs <2) 
		return 0.0; /*sc is all outside patch, attractive energy is 0*/
	T1=intersections[0]; /*points on sc2*/
	T2=intersections[1];
	/*2- now do the same oposite way psc1 in patch of psc2*/
	for(i=0;i<5;i++) 
		intersections[i]=0;
	vec1=vec_scale(interact->r_cm,-1.0);
	intrs=cpsc_intersect(interact->part2,interact->part1,interact->param->half_len[1],interact->param->half_len[0],vec1, intersections, rcut, interact->param,1, patchnum2);
	if (intrs <2) 
		return 0.0; /*sc is all outside patch, attractive energy is 0*/
	S1=intersections[0]; /*points on sc1*/
	S2=intersections[1];

	/*3- scaling function1: dependence on the length of intersetions*/
	v1=fabs(S1-S2);
	v2=fabs(T1-T2);
	f0=0.5*(v1+v2);
	
	/*4a- with two intersection pices calculate vector between their CM 
	  -this is for angular orientation*/
	vec1=vec_scale(interact->part1->dir,(S1+S2)*0.5);
	vec2=vec_scale(interact->part2->dir,(T1+T2)*0.5);
	vec_intrs.x=vec2.x-vec1.x-interact->r_cm.x;
	vec_intrs.y=vec2.y-vec1.y-interact->r_cm.y;
	vec_intrs.z=vec2.z-vec1.z-interact->r_cm.z;
	/*vec_intrs should be from sc1 to sc2*/
	//	    fprintf (stderr, "segments_CM: %.8f %.8f %.8f \n",vec_intrs.x,vec_intrs.y,vec_intrs.z);

	/*4b - calculate closest distance attractive energy from it*/
	vec_mindist = mindist_segments(interact->part1->dir,v1,interact->part2->dir,v2,vec_intrs);
	//	    fprintf (stderr, "segments closest dist: %.8f %.8f %.8f \n",vec_mindist.x,vec_mindist.y,vec_mindist.z);
	ndist=sqrt(DOT(vec_mindist,vec_mindist));
	//dist=DOT(vec_intrs,vec_intrs);
	if (ndist < interact->param->pdis) 
		atrenergy = -interact->param->epsilon;
	else {
		atrenergy = cos(PIH*(ndist-interact->param->pdis)/interact->param->pswitch);
		atrenergy *= -atrenergy*interact->param->epsilon ;
	}

	/*5- scaling function2: angular dependence of patch1*/
	vec1=vec_scale(vec_intrs,1.0);
	//vec1=vec_scale(vec_mindist,-1.0);
	vec1=vec_perpproject(&vec1, &interact->part1->dir);
	normalise(&vec1);
	a = DOT(vec1,interact->part1->patchdir[patchnum1]);
	f1 = fanglscale(a,interact->param, 0+2*patchnum1);

	/*6- scaling function3: angular dependence of patch2*/
	vec1=vec_scale(vec_intrs,-1.0);
	//vec1=vec_scale(vec_mindist,1.0);
	vec1=vec_perpproject(&vec1, &interact->part2->dir);
	normalise(&vec1);
	a = DOT(vec1,interact->part2->patchdir[patchnum2]);
	f2 = fanglscale(a,interact->param, 1+2*patchnum2);

	/*7- put it all together*/
	atrenergy *=f0*f1*f2;
	//if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
	//	    fprintf (stderr, "attraction  %.8f \n",atrenergy);
	//	    exit(1);

	return atrenergy;
}


/*..............................................................................*/
/* 
   Determines total energy of spherocylinders type PSC and CPSC 
 */

double e_psc_cpsc(struct interacts * interact)
{
	double atrenergy, repenergy;

	void closestdist(struct interacts *);
	double erepulsive(struct interacts *);
	double eattractive_psc_cpsc(struct interacts *,int,int);

	//DEBUG_SIM("do energy 23") ;
	closestdist(interact);
	repenergy = erepulsive(interact);
	//DEBUG_SIM("got the rep. energy");
	if ( ( interact->dist > interact->param->rcut ) || ( interact->param->epsilon == 0.0 ) ) 
		atrenergy = 0.0;
	else {
	    BOOL firstCH=FALSE, secondCH=FALSE;
		struct vector olddir1 = interact->part1->dir;
		struct vector olddir2 = interact->part2->dir;
		if ((interact->param->geotype[0] == CHPSC) || (interact->param->geotype[0] == CHCPSC)||
		  (interact->param->geotype[0] == TCHPSC) || (interact->param->geotype[0] == TCHCPSC) )
			firstCH = TRUE;
		if ((interact->param->geotype[1] == CHPSC) || (interact->param->geotype[1] == CHCPSC)||
		  (interact->param->geotype[1] == TCHPSC) || (interact->param->geotype[1] == TCHCPSC) )
			secondCH = TRUE;
		if(firstCH)
			interact->part1->dir = interact->part1->chdir[0];
		if(secondCH)
			interact->part2->dir = interact->part2->chdir[0];
		
		if ((firstCH) || (secondCH) ) {
		    closestdist(interact);
		}
		atrenergy = eattractive_psc_cpsc(interact,0,0);
		
		
		/*addition of interaction of second patches*/
		if ( (interact->param->geotype[0] == TCPSC) || (interact->param->geotype[0] == TCHCPSC) || 
		      (interact->param->geotype[0] == TPSC) || (interact->param->geotype[0] == TCHPSC) ||
		      (interact->param->geotype[1] == TCPSC) ||(interact->param->geotype[1] == TCHCPSC) || 
		      (interact->param->geotype[1] == TPSC) || (interact->param->geotype[1] == TCHPSC) ) {
		    BOOL firstT=FALSE, secondT=FALSE;
		    if ( (interact->param->geotype[0] == TCPSC) || (interact->param->geotype[0] == TCHCPSC) || 
			    (interact->param->geotype[0] == TPSC) || (interact->param->geotype[0] == TCHPSC) )
			firstT = TRUE;
		    if ( (interact->param->geotype[1] == TCPSC) || (interact->param->geotype[1] == TCHCPSC) || 
			    (interact->param->geotype[1] == TPSC) || (interact->param->geotype[1] == TCHPSC) )
			secondT = TRUE;
		    
		    if (firstT) {
			if (firstCH) {
			    interact->part1->dir = interact->part1->chdir[1];
			    closestdist(interact);
			}
			atrenergy += eattractive_psc_cpsc(interact,1,0);
		    }
		    if ( (firstT) && (secondT) ) {
			if (secondCH) {
			    interact->part2->dir = interact->part2->chdir[1];
			    closestdist(interact);
			}
			atrenergy += eattractive_psc_cpsc(interact,1,1);
		    }
		    if (secondT) {
			if (firstT && firstCH ) {
			    interact->part1->dir = interact->part1->chdir[0];
			    closestdist(interact);
			}
			atrenergy += eattractive_psc_cpsc(interact,0,1);
		    }
		}
		
		if (firstCH) 
		    interact->part1->dir = olddir1;
		if (secondCH)
		    interact->part2->dir = olddir2;
		
	}
	return repenergy+atrenergy;
}

/* 
   Determines attractive energy of spherocylinders type PSC and CPSC 
 */
double eattractive_psc_cpsc(struct interacts * interact,int patchnum1,int patchnum2)
{
	int i, intrs;
	double rcut, atrenergy, ndist;
	double v1, v2, f0, f1, f2, T1, T2, S1, S2, a;
	double intersections[5];
	struct vector vec1, vec2, vec_intrs, vec_mindist;

	struct vector vec_sub(struct vector, struct vector);
	struct vector vec_sum(struct vector, struct vector);
	struct vector vec_create(double, double, double);
	struct vector vec_scale(struct vector, double);
	struct vector vec_perpproject(struct vector*, struct vector*);
	struct vector mindist_segments(struct vector, double, struct vector, double, struct vector);
	void normalise(struct vector *);
	int psc_intersect(struct particles *, struct particles *,
			double, double, struct vector, double *,double, struct ia_param *, int which,int patchnum);
	int cpsc_intersect(struct particles *, struct particles *,
			double, double, struct vector, double *,double, struct ia_param *, int which,int patchnum);
	double fanglscale(double, struct ia_param *, int which);

	rcut = interact->param->rcut;

	//interact->halfl = interact->param->half_len[0];
	//DEBUG_SIM("halfl = %lf", interact->halfl);
	for(i=0;i<5;i++) 
		intersections[i]=0;
	BOOL first;
	if ( (interact->param->geotype[0] == PSC)||(interact->param->geotype[0] == CHPSC)||(interact->param->geotype[0] == TPSC)||(interact->param->geotype[0] == TCHPSC) ){
		first = TRUE;
	} else {
		first = FALSE;
	}
	//cospatch = param.pcanglsw;
	//cospatchinr = param.pcangl;
	/*1- do intersections of spherocylinder2 with patch of spherocylinder1 at. 
	  cut distance C*/
	//DEBUG_SIM("first intersection");
	if (first) { 
		intrs=psc_intersect(interact->part1,interact->part2,interact->param->half_len[0],interact->param->half_len[1],interact->r_cm, intersections, rcut, interact->param,0, patchnum1);
	} else {
		intrs=cpsc_intersect(interact->part1,interact->part2,interact->param->half_len[0],interact->param->half_len[1],interact->r_cm, intersections, rcut, interact->param,0, patchnum1);
	}
	//DEBUG_SIM("first intersection: done");
	if (intrs <2){ 
		//DEBUG_SIM("No intersection :(");
		return 0.0; /*sc is all outside patch, attractive energy is 0*/
	}

	T1=intersections[0]; /*points on sc2*/
	T2=intersections[1];
	/*2- now do the same oposite way psc1 in patch of psc2*/
	for(i=0;i<5;i++) 
		intersections[i]=0;
	//DEBUG_SIM("get vector");
	vec1=vec_scale(interact->r_cm,-1.0);
	//DEBUG_SIM("second intersection");
	if (first) {
		intrs=cpsc_intersect(interact->part2,interact->part1,interact->param->half_len[1],interact->param->half_len[0],vec1, intersections, rcut, interact->param,1, patchnum2);
	} else {
		intrs=psc_intersect(interact->part2,interact->part1,interact->param->half_len[1],interact->param->half_len[0],vec1, intersections, rcut, interact->param,1, patchnum2);
	}
	if (intrs <2) 
		return 0.0; /*sc is all outside patch, attractive energy is 0*/

	S1=intersections[0]; /*points on sc1*/
	S2=intersections[1];

	/*3- scaling function1: dependence on the length of intersetions*/
	v1=fabs(S1-S2);
	v2=fabs(T1-T2);
	f0=0.5*(v1+v2);
	
	/*4a- with two intersection pices calculate vector between their CM 
	  -this is for angular orientation*/
	vec1=vec_scale(interact->part1->dir,(S1+S2)*0.5);
	vec2=vec_scale(interact->part2->dir,(T1+T2)*0.5);
	vec_intrs.x=vec2.x-vec1.x-interact->r_cm.x;
	vec_intrs.y=vec2.y-vec1.y-interact->r_cm.y;
	vec_intrs.z=vec2.z-vec1.z-interact->r_cm.z;
	/*vec_intrs should be from sc1 to sc2*/
	//	    fprintf (stderr, "segments_CM: %.8f %.8f %.8f \n",vec_intrs.x,vec_intrs.y,vec_intrs.z);

	/*4b - calculate closest distance attractive energy from it*/
	vec_mindist = mindist_segments(interact->part1->dir,v1,interact->part2->dir,v2,vec_intrs);
	//	    fprintf (stderr, "segments closest dist: %.8f %.8f %.8f \n",vec_mindist.x,vec_mindist.y,vec_mindist.z);
	ndist=sqrt(DOT(vec_mindist,vec_mindist));
	//dist=DOT(vec_intrs,vec_intrs);
	if (ndist < interact->param->pdis) 
		atrenergy = -interact->param->epsilon;
	//atrenergy = -1.0;
	else {
		atrenergy = cos(PIH*(ndist-interact->param->pdis)/interact->param->pswitch);
		atrenergy *= -atrenergy*interact->param->epsilon ;
	}

	/*5- scaling function2: angular dependence of patch1*/
	vec1=vec_scale(vec_intrs,1.0);
	//vec1=vec_scale(vec_mindist,-1.0);
	vec1=vec_perpproject(&vec1, &interact->part1->dir);
	normalise(&vec1);
	a = DOT(vec1,interact->part1->patchdir[patchnum1]);
	f1 = fanglscale(a,interact->param, 0+2*patchnum1);

	/*6- scaling function3: angular dependence of patch2*/
	vec1=vec_scale(vec_intrs,-1.0);
	//vec1=vec_scale(vec_mindist,1.0);
	vec1=vec_perpproject(&vec1, &interact->part2->dir);
	normalise(&vec1);
	a = DOT(vec1,interact->part2->patchdir[patchnum2]);
	f2 = fanglscale(a,interact->param, 1+2*patchnum2);

	/*7- put it all together*/
	atrenergy *=f0*f1*f2;
	//if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
	//	    fprintf (stderr, "attraction  %.8f \n",atrenergy);
	//	    exit(1);

	return atrenergy;
}

/*..............................................................................*/
/* 
 * Determines total energy of spherocylinder type 1 and sphere type 11 
 */

double e_spa_sca(struct interacts * interact)
{
	double atrenergy, repenergy, b, f0, halfl;

	struct vector vec_perpproject(struct vector *, struct vector *);
	void normalise(struct vector *);
	void closestdist(struct interacts *);
	double erepulsive(struct interacts *);
	double fanglscale(double, struct ia_param *, int which);

	//DEBUG    printf ("do energy 111 \n\n");
	closestdist(interact);
	repenergy = erepulsive(interact);

	if ( ( interact->dist > interact->param->rcut ) || ( interact->param->epsilon == 0.0 ) ) 
		atrenergy = 0.0;
	else {
		/*calculate closest distance attractive energy*/
		if (interact->dist < interact->param->pdis) 
			atrenergy = -interact->param->epsilon;
		else {
			atrenergy = cos(PIH*(interact->dist-interact->param->pdis)/interact->param->pswitch);
			atrenergy *= -atrenergy*interact->param->epsilon ;
		}

		/*scaling function for the length of spherocylinder within cutoff*/

		if (interact->param->geotype [0] < SP)
			halfl = interact->param->half_len[0];
		else
		halfl = interact->param->half_len[1];		
		b = sqrt(interact->param->rcut*interact->param->rcut-interact->dist*interact->dist);
		if ( interact->contt + b > halfl ) 
			f0 = halfl;
		else 
			f0 = interact->contt + b;
		if ( interact->contt - b < -halfl ) 
			f0 -= -halfl;
		else 
			f0 -= interact->contt - b;
		atrenergy *= f0;
		//if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
		//fprintf (stderr, "attraction211  %.8f x: %.8f y: %.8f z: %.8f \n",atrenergy,vec1.x,vec1.y,vec1.z);
		//exit(1);
	}

	return repenergy+atrenergy;
}


/*..............................................................................*/
/* 
 * Determines total energy of spherocylinder type 2 and sphere type 11 
 */
double e_psc_spa(struct interacts * interact)
{
	double atrenergy, repenergy;

	void closestdist(struct interacts *);
	double erepulsive(struct interacts *);
	double eattractive_psc_spa(struct interacts *, int);

	//DEBUG_SIM("do energy 211") ;
	closestdist(interact);
	repenergy = erepulsive(interact);
	//DEBUG_SIM("got the rep. energy");
	if ( ( interact->dist > interact->param->rcut ) || ( interact->param->epsilon == 0.0 ) ) 
		atrenergy = 0.0;
	else {
		BOOL firstCH=FALSE, secondCH=FALSE;
		struct vector olddir1 = interact->part1->dir;
		struct vector olddir2 = interact->part2->dir;
		if ( (interact->param->geotype[0] == CHPSC) || (interact->param->geotype[0] == TCHPSC) )
			firstCH = TRUE;
		if ( (interact->param->geotype[1] == CHPSC) || (interact->param->geotype[1] == TCHPSC) )
			secondCH = TRUE;
		if(firstCH)
			interact->part1->dir = interact->part1->chdir[0];
		if(secondCH)
			interact->part2->dir = interact->part2->chdir[0];
		
		if ((firstCH) || (secondCH) ) {
		    closestdist(interact);
		}
		atrenergy = eattractive_psc_spa(interact,0);
		
		/*addition of interaction of second patches*/
		if ( (interact->param->geotype[0] == TPSC) || (interact->param->geotype[0] == TCHPSC) || 
		  (interact->param->geotype[1] == TPSC) ||(interact->param->geotype[1] == TCHPSC) ) {
		    BOOL firstT=FALSE, secondT=FALSE;
		    if ( (interact->param->geotype[0] == TPSC) || (interact->param->geotype[0] == TCHPSC) )
			firstT = TRUE;
		    if ( (interact->param->geotype[1] == TPSC) ||(interact->param->geotype[1] == TCHPSC)  )
			secondT = TRUE;
		    
		    if (firstT) {
			if (firstCH) {
			    interact->part1->dir = interact->part1->chdir[1];
			    closestdist(interact);
			}
			atrenergy += eattractive_psc_spa(interact,1);
		    }
		    if (secondT) {
			if(secondCH) {
			    interact->part2->dir = interact->part2->chdir[1];
			    closestdist(interact);
			}
			atrenergy += eattractive_psc_spa(interact,1);
		    }
		    if ( (firstT) && (secondT) ) {
		      fprintf (stderr, "ERROR PSC should interact s SPA but got two PSC \n");
		      exit(1);
		    }
		}
		

		if (firstCH) 
		    interact->part1->dir = olddir1;
		if (secondCH)
		    interact->part2->dir = olddir2;
		
	}
	return repenergy+atrenergy;
}
/* 
 * Determines attractive energy of spherocylinder type 2 and sphere type 11 
 */

double eattractive_psc_spa(struct interacts * interact, int patchnum1)
{
	double atrenergy, a, b, f0, halfl;
	struct vector vec1;

	struct vector vec_perpproject(struct vector *, struct vector*);
	void normalise(struct vector *);
	double fanglscale(double, struct ia_param *, int which);
	int which;

	/*calculate closest distance attractive energy*/
	if (interact->dist < interact->param->pdis)
		atrenergy = -interact->param->epsilon;
	else {
		atrenergy = cos(PIH*(interact->dist-interact->param->pdis)/interact->param->pswitch);
		atrenergy *= -atrenergy*interact->param->epsilon ;
	}
	/*scaling function: angular dependence of patch1*/
	if (interact->param->geotype[0] < SP) {
		which = 0;
		vec1=vec_perpproject(&interact->distvec, &interact->part1->dir);
		normalise(&vec1);
		a = DOT(vec1,interact->part1->patchdir[patchnum1]);
        halfl=interact->param->half_len[0];
	} else {
		which = 1;
		vec1=vec_perpproject(&interact->distvec, &interact->part2->dir);
		normalise(&vec1);
		a = DOT(vec1,interact->part2->patchdir[patchnum1]);
        halfl=interact->param->half_len[1];
	}
	/*scaling function for the length of spherocylinder within cutoff*/

	b = sqrt(interact->param->rcut*interact->param->rcut-interact->dist*interact->dist);
	if ( interact->contt + b > halfl ) 
		f0 = halfl;
	else 
		f0 = interact->contt + b;
	if ( interact->contt - b < -halfl ) 
		f0 -= -halfl;
	else 
		f0 -= interact->contt - b;
	atrenergy *= fanglscale(a,interact->param, which)*f0;
	//if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
	//fprintf (stderr, "attraction211  %.8f x: %.8f y: %.8f z: %.8f \n",atrenergy,vec1.x,vec1.y,vec1.z);
	//exit(1);

	return atrenergy;
}

/*..............................................................................*/
/* 
   Determines total energy of spherocylinder type 3 and sphere type 11 
 */

double e_cpsc_spa(struct interacts * interact)
{
	double atrenergy, repenergy, halfl;

	void closestdist(struct interacts *);
	double erepulsive(struct interacts *);
	double eattractive_cpsc_spa(struct interacts *,int);

	//DEBUG_SIM("do energy 311") ;
	closestdist(interact);
	repenergy = erepulsive(interact);
	//DEBUG_SIM("got the rep. energy");
    
    if (interact->param->geotype[0] < SP) {
        halfl=interact->param->half_len[0];
	} else {
        halfl=interact->param->half_len[1];
	}
	if ( ( interact->dist > interact->param->rcut ) || ( interact->param->epsilon == 0.0 ) || ( interact->dist > interact->param->rcut )
	  || (interact->contt > halfl) || (interact->contt < -halfl) ) 
		atrenergy = 0.0;
	else {
		BOOL firstCH=FALSE, secondCH=FALSE;
		struct vector olddir1 = interact->part1->dir;
		struct vector olddir2 = interact->part2->dir;
		if ( (interact->param->geotype[0] == CHCPSC) || (interact->param->geotype[0] == TCHCPSC) )
			firstCH = TRUE;
		if ( (interact->param->geotype[1] == CHCPSC) || (interact->param->geotype[1] == TCHCPSC) )
			secondCH = TRUE;
		if(firstCH)
			interact->part1->dir = interact->part1->chdir[0];
		if(secondCH)
			interact->part2->dir = interact->part2->chdir[0];
		
		if ((firstCH) || (secondCH) ) {
		    closestdist(interact);
		}
		atrenergy = eattractive_cpsc_spa(interact,0);
		
		/*addition of interaction of second patches*/
		if ( (interact->param->geotype[0] == TCPSC) || (interact->param->geotype[0] == TCHCPSC) || 
		  (interact->param->geotype[1] == TCPSC) ||(interact->param->geotype[1] == TCHCPSC) ) {
		    BOOL firstT=FALSE, secondT=FALSE;
		    if ( (interact->param->geotype[0] == TCPSC) || (interact->param->geotype[0] == TCHCPSC) )
			firstT = TRUE;
		    if ( (interact->param->geotype[1] == TCPSC) ||(interact->param->geotype[1] == TCHCPSC)  )
			secondT = TRUE;
		    
		    if (firstT) {
			if (firstCH) {
			    interact->part1->dir = interact->part1->chdir[1];
			    closestdist(interact);
			}
			atrenergy += eattractive_cpsc_cpsc(interact,1,0);
		    }
		    if (secondT) {
			if(secondCH) {
			    interact->part2->dir = interact->part2->chdir[1];
			    closestdist(interact);
			}
			atrenergy += eattractive_cpsc_cpsc(interact,0,1);
		    }
		    if ( (firstT) && (secondT) ) {
			fprintf (stderr, "ERROR PSC should interact s SPA but got two PSC \n");
			exit(1);
		    }
		}

		if (firstCH) 
		    interact->part1->dir = olddir1;
		if (secondCH)
		    interact->part2->dir = olddir2;
		
	}
	return repenergy+atrenergy;
}

/* 
   Determines attractive energy of spherocylinder type 3 and sphere type 11 
 */

double eattractive_cpsc_spa(struct interacts * interact,int patchnum1)
{
	double atrenergy, a, b, f0, halfl;
	struct vector vec1;
	int which;

	struct vector vec_perpproject(struct vector *, struct vector*);
	void normalise(struct vector *);
	double fanglscale(double, struct ia_param *, int which);


	/*if it is in cylindrical part c>-halfl and c<halfl*/
	/*calculate closest distance attractive energy*/
	if (interact->dist < interact->param->pdis) 
		atrenergy = -interact->param->epsilon;
	else {
		atrenergy = cos(PIH*(interact->dist-interact->param->pdis)/interact->param->pswitch);
		atrenergy *= -atrenergy*interact->param->epsilon ;
	}
	/*scaling function: angular dependence of patch1*/
	if (interact->param->geotype[0] < SP) {
		which = 0;
		vec1=vec_perpproject(&interact->distvec, &interact->part1->dir);
		normalise(&vec1);
		a = DOT(vec1,interact->part1->patchdir[patchnum1]);
        halfl = interact->param->half_len[0];
	} else {
		which = 1;
		vec1=vec_perpproject(&interact->distvec, &interact->part2->dir);
		normalise(&vec1);
		a = DOT(vec1,interact->part2->patchdir[patchnum1]);
        halfl = interact->param->half_len[1];
	}
	/*scaling function for the length of spherocylinder within cutoff*/
	b = sqrt(interact->param->rcut*interact->param->rcut-interact->dist*interact->dist);
	if ( interact->contt + b > halfl ) 
		f0 = halfl;
	else 
		f0 = interact->contt + b;
	if ( interact->contt - b < -halfl ) 
		f0 -= -halfl;
	else 
		f0 -= interact->contt - b;

	atrenergy *= fanglscale(a,interact->param, which)*f0;

	//if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
	//fprintf (stderr, "attraction311  %.8f a: %.8f\n",atrenergy,a);
	//exit(1);

	return atrenergy;
}



/*..............................................................................*/
/* 
   Determines total energy of two spherocylinders type 11 
 */

double e_2sca_or_2spa(struct interacts * interact)
{
	double repenergy, atrenergy;

	double erepulsive(struct interacts *);
	void closestdist(struct interacts *);

	closestdist(interact);
	repenergy = erepulsive(interact);
	if ( ( interact->dist > interact->param->rcut ) || ( interact->param->epsilon == 0.0 ) ) 
		atrenergy = 0.0;
	else {
		if (interact->dist < interact->param->pdis) 
			atrenergy = -interact->param->epsilon;
		else  {
			atrenergy = cos(PIH*(interact->dist-interact->param->pdis)/interact->param->pswitch);
			atrenergy *= -atrenergy*interact->param->epsilon ;
		}
	}

	return repenergy+atrenergy;
}


/*..............................................................................*/
/* 
   Determines total energy with purely repulsive types 
 */

double e_spn_or_scn(struct interacts * interact)
{
	double repenergy;

	double erepulsive(struct interacts *);
	void closestdist(struct interacts *);

	closestdist(interact);
	repenergy = erepulsive(interact);

	return repenergy;
}


/*..............................................................................*/
/*
   Determines repulsive energy of two spherocylinders 
 */

double erepulsive(struct interacts * interact)
{
	double repenergy, en6;
    
	/* WCA repulsion */
	if (interact->dist > interact->param->rcutwca) repenergy = 0.0;
	else {
            en6 = pow((interact->param->sigma/interact->dist),6);
            repenergy = interact->param->epsilon * ( 4*en6*(en6-1) + 1.0);
	}
        //int Digs = 20;
        //printf("dist: %.*e, repenergy: %.*e\n",Digs, interact->dist, Digs, repenergy);

	return repenergy;
}

/*..............................................................................*/
/*
   Indicates not yet programmed interaction
 */

double enoexist(struct interacts * interact)
{
	double energy=0.0;

	fprintf (stderr, "ERROR: We have not programed interaction of types %d and %d\n",
			interact->part1->type,interact->part2->type);
	exit (1);

	return energy;
}

/* function for calculation of harmonic potential*/
double harmonic(double aktualvalue, double eqvalue, double springconst)
{
	return springconst*(aktualvalue-eqvalue)*(aktualvalue-eqvalue)*0.5;
}

/*..............................................................................*/
/*
   Determines bond energy 
 */
double bondenergy(long num1, long num2, struct interacts * interact, struct topo * topo, struct conf * conf)
{
	double energy=0.0, bondlength, halfl;
	struct vector vec1, vec2, vecbond;
	int * geotype = interact->param->geotype;

	struct vector image(struct vector, struct vector, struct vector);
	double harmonic(double, double, double);


	/*interaction with nearest neighbours -harmonic*/
	if ((topo->chainparam[conf->particle[num1].chaint]).bond1c >= 0) {
		 
		if (num2 == topo->conlist[num1][1]) {
			/*num1 is connected to num2 by tail*/
			if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
				energy = harmonic(interact->distcm,topo->chainparam[conf->particle[num1].chaint].bond1eq,topo->chainparam[conf->particle[num1].chaint].bond1c);
			else {
				if (geotype[0] < SP) 
					halfl=interact->param->half_len[0]; 
				else 
					halfl = 0.0;
				vec1.x=conf->particle[num1].pos.x - conf->particle[num1].dir.x * halfl /conf->box.x;
				vec1.y=conf->particle[num1].pos.y - conf->particle[num1].dir.y * halfl /conf->box.y;
				vec1.z=conf->particle[num1].pos.z - conf->particle[num1].dir.z * halfl /conf->box.z;
				if (geotype[1] < SP) 
					halfl=interact->param->half_len[1]; 
				else 
					halfl = 0.0;
				vec2.x=conf->particle[num2].pos.x + conf->particle[num2].dir.x * halfl /conf->box.x;
				vec2.y=conf->particle[num2].pos.y + conf->particle[num2].dir.y * halfl /conf->box.y;
				vec2.z=conf->particle[num2].pos.z + conf->particle[num2].dir.z * halfl /conf->box.z;
				vecbond = image(vec1, vec2, conf->box);
				bondlength = sqrt(DOT(vecbond,vecbond));
				energy = harmonic(bondlength,topo->chainparam[conf->particle[num1].chaint].bond1eq,topo->chainparam[conf->particle[num1].chaint].bond1c);
			}
		} else {
			if (num2 == topo->conlist[num1][0]) {
				/*num1 is connected to num2 by head*/
				if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
					energy = harmonic(interact->distcm,topo->chainparam[conf->particle[num1].chaint].bond1eq,topo->chainparam[conf->particle[num1].chaint].bond1c);
				else {
					if (geotype[0] < SP) 
						halfl=interact->param->half_len[0]; 
					else 
						halfl = 0.0;
					vec1.x=conf->particle[num1].pos.x + conf->particle[num1].dir.x * halfl /conf->box.x;
					vec1.y=conf->particle[num1].pos.y + conf->particle[num1].dir.y * halfl /conf->box.y;
					vec1.z=conf->particle[num1].pos.z + conf->particle[num1].dir.z * halfl /conf->box.z;
					if (geotype[0] < SP) 
						halfl=interact->param->half_len[0]; 
					else 
						halfl = 0.0;
					vec2.x=conf->particle[num2].pos.x - conf->particle[num2].dir.x * halfl /conf->box.x;
					vec2.y=conf->particle[num2].pos.y - conf->particle[num2].dir.y * halfl /conf->box.y;
					vec2.z=conf->particle[num2].pos.z - conf->particle[num2].dir.z * halfl /conf->box.z;
					vecbond = image(vec1, vec2, conf->box);
					bondlength = sqrt(DOT(vecbond,vecbond));
					energy = harmonic(bondlength,topo->chainparam[conf->particle[num1].chaint].bond1eq,topo->chainparam[conf->particle[num1].chaint].bond1c);
				}
			}
		} 
	}
	/*interaction with second nearest neighbours -harmonic*/
	if (topo->chainparam[conf->particle[num1].chaint].bond2c >= 0) {
		if (num2 == topo->conlist[num1][2]) {
			/*num1 is connected to num2 by tail*/
			if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
				energy = harmonic(interact->distcm,topo->chainparam[conf->particle[num1].chaint].bond2eq,topo->chainparam[conf->particle[num1].chaint].bond2c);
			else {
				vecbond = image(conf->particle[num1].pos, conf->particle[num2].pos, conf->box);
				bondlength = sqrt(DOT(vecbond,vecbond));
				energy = harmonic(bondlength,topo->chainparam[conf->particle[num1].chaint].bond2eq,topo->chainparam[conf->particle[num1].chaint].bond2c);
			}
		} else {
			if (num2 == topo->conlist[num1][3]) {
				/*num1 is connected to num2 by head*/
				if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
					energy = harmonic(interact->distcm,topo->chainparam[conf->particle[num1].chaint].bond2eq,topo->chainparam[conf->particle[num1].chaint].bond2c);
				else {
					vecbond = image(conf->particle[num1].pos, conf->particle[num2].pos, conf->box);
					bondlength = sqrt(DOT(vecbond,vecbond));
					energy = harmonic(bondlength,topo->chainparam[conf->particle[num1].chaint].bond2eq,topo->chainparam[conf->particle[num1].chaint].bond2c);
				}
			}
		}
	}
    /*interaction with nearest neighbours - direct harmonic bond*/
	if ((topo->chainparam[conf->particle[num1].chaint]).bonddc > 0) {
        
		if (num2 == topo->conlist[num1][1]) {
			/*num1 is connected to num2 by tail*/
			if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
				energy = harmonic(interact->distcm,topo->chainparam[conf->particle[num1].chaint].bonddeq,topo->chainparam[conf->particle[num1].chaint].bonddc);
			else {
				if (geotype[0] < SP) 
					halfl=interact->param->half_len[0]; 
				else 
					halfl = 0.0;
				vec1.x=conf->particle[num1].pos.x - conf->particle[num1].dir.x * halfl /conf->box.x;
				vec1.y=conf->particle[num1].pos.y - conf->particle[num1].dir.y * halfl /conf->box.y;
				vec1.z=conf->particle[num1].pos.z - conf->particle[num1].dir.z * halfl /conf->box.z;
				if (geotype[1] < SP) 
					halfl=interact->param->half_len[1]; 
				else 
					halfl = 0.0;
				vec2.x=conf->particle[num2].pos.x + conf->particle[num2].dir.x * (halfl + topo->chainparam[conf->particle[num1].chaint].bonddeq) /conf->box.x ;
				vec2.y=conf->particle[num2].pos.y + conf->particle[num2].dir.y * (halfl + topo->chainparam[conf->particle[num1].chaint].bonddeq) /conf->box.y ;
				vec2.z=conf->particle[num2].pos.z + conf->particle[num2].dir.z * (halfl + topo->chainparam[conf->particle[num1].chaint].bonddeq) /conf->box.z ;
				vecbond = image(vec1, vec2, conf->box);
				bondlength = sqrt(DOT(vecbond,vecbond));
				energy = harmonic(bondlength,0.0,topo->chainparam[conf->particle[num1].chaint].bonddc);
			}
		} else {
			if (num2 == topo->conlist[num1][0]) {
				/*num1 is connected to num2 by head*/
				if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
					energy = harmonic(interact->distcm,topo->chainparam[conf->particle[num1].chaint].bond1eq,topo->chainparam[conf->particle[num1].chaint].bond1c);
				else {
					if (geotype[0] < SP) 
						halfl=interact->param->half_len[0]; 
					else 
						halfl = 0.0;
					vec1.x=conf->particle[num1].pos.x + conf->particle[num1].dir.x * (halfl + topo->chainparam[conf->particle[num1].chaint].bonddeq) /conf->box.x ;
					vec1.y=conf->particle[num1].pos.y + conf->particle[num1].dir.y * (halfl + topo->chainparam[conf->particle[num1].chaint].bonddeq) /conf->box.y ;
					vec1.z=conf->particle[num1].pos.z + conf->particle[num1].dir.z * (halfl + topo->chainparam[conf->particle[num1].chaint].bonddeq) /conf->box.z ;
					if (geotype[0] < SP) 
						halfl=interact->param->half_len[0]; 
					else 
						halfl = 0.0;
					vec2.x=conf->particle[num2].pos.x - conf->particle[num2].dir.x * halfl /conf->box.x;
					vec2.y=conf->particle[num2].pos.y - conf->particle[num2].dir.y * halfl /conf->box.y;
					vec2.z=conf->particle[num2].pos.z - conf->particle[num2].dir.z * halfl /conf->box.z;
					vecbond = image(vec1, vec2, conf->box);
					bondlength = sqrt(DOT(vecbond,vecbond));
					energy = harmonic(bondlength,0.0,topo->chainparam[conf->particle[num1].chaint].bonddc);
				}
			}
		} 
	}
	//printf("bondlength: %f\n",bondlength);
	//    printf("bondener: %f\n",energy);
	return energy;
}


/*..............................................................................*/
/*
 Determines angle energy between spherocylinders 
 */
double angleenergy(long num1, long num2, struct interacts * interact, struct topo * topo, struct conf * conf)
{
	double energy=0.0, currangle, halfl;
	struct vector vec1, vec2;
	int * geotype = interact->param->geotype;
    
	struct vector image(struct vector, struct vector, struct vector);
    void normalise(struct vector *);
	double harmonic(double, double, double);
    
    
	/*angle interaction with nearest neighbours -harmonic*/
	if ((topo->chainparam[conf->particle[num1].chaint]).angle1c >= 0) {
		if (num2 == topo->conlist[num1][0]) {
			/*num1 is connected to num2 by tail*/
			if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                /*spheres do not have this interaction*/
				energy += 0.0;
			else {
				if (geotype[0] < SP) 
                    vec1 = conf->particle[num1].dir;
				else {
					halfl=interact->param->half_len[1];
                    //sphere angle is defined versus the end of spherocylinder
                    vec1.x=conf->particle[num2].pos.x - conf->particle[num2].dir.x * halfl /conf->box.x;
                    vec1.y=conf->particle[num2].pos.y - conf->particle[num2].dir.y * halfl /conf->box.y;
                    vec1.z=conf->particle[num2].pos.z - conf->particle[num2].dir.z * halfl /conf->box.z;
                    vec1 = image(vec1, conf->particle[num1].pos, conf->box);
                }
                if (geotype[1] < SP) 
					vec2 = conf->particle[num2].dir;
				else {
                    halfl=interact->param->half_len[0];
                    vec2.x=conf->particle[num1].pos.x + conf->particle[num1].dir.x * halfl /conf->box.x;
                    vec2.y=conf->particle[num1].pos.y + conf->particle[num1].dir.y * halfl /conf->box.y;
                    vec2.z=conf->particle[num1].pos.z + conf->particle[num1].dir.z * halfl /conf->box.z;
                    vec2 = image(vec2, conf->particle[num2].pos, conf->box);
				}
                normalise(&vec1);
                normalise(&vec2);
				currangle = acos(DOT(vec1,vec2));
				energy += harmonic(currangle,topo->chainparam[conf->particle[num1].chaint].angle1eq,topo->chainparam[conf->particle[num1].chaint].angle1c);
			}
		} else {
			if (num2 == topo->conlist[num1][1]) {
				/*num1 is connected to num2 by head*/
				if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                /*spheres do not have this interaction*/
                    energy += 0.0;
                else {
                    if (geotype[0] < SP) 
                        vec1 = conf->particle[num1].dir;
                    else {
                        halfl=interact->param->half_len[1];
                        //sphere angle is defined versus the end of spherocylinder
                        vec1.x=conf->particle[num2].pos.x + conf->particle[num2].dir.x * halfl /conf->box.x;
                        vec1.y=conf->particle[num2].pos.y + conf->particle[num2].dir.y * halfl /conf->box.y;
                        vec1.z=conf->particle[num2].pos.z + conf->particle[num2].dir.z * halfl /conf->box.z;
                        vec1 = image(vec1, conf->particle[num1].pos, conf->box);
                    }
                    if (geotype[1] < SP) 
                        vec2 = conf->particle[num2].dir;
                    else {
                        halfl=interact->param->half_len[0];
                        vec2.x=conf->particle[num1].pos.x - conf->particle[num1].dir.x * halfl /conf->box.x;
                        vec2.y=conf->particle[num1].pos.y - conf->particle[num1].dir.y * halfl /conf->box.y;
                        vec2.z=conf->particle[num1].pos.z - conf->particle[num1].dir.z * halfl /conf->box.z;
                        vec2 = image(vec2, conf->particle[num2].pos, conf->box);
                    }
                    normalise(&vec1);
                    normalise(&vec2);
                    currangle = acos(DOT(vec1,vec2));
                    energy += harmonic(currangle,topo->chainparam[conf->particle[num2].chaint].angle1eq,topo->chainparam[conf->particle[num2].chaint].angle1c);
                }
			}
		} 
	}
    
	/*interaction between the orientation of  spherocylinders patches -harmonic*/
	if (topo->chainparam[conf->particle[num1].chaint].angle2c >= 0) {
        if (num2 == topo->conlist[num1][0]) {
			/*num1 is connected to num2 by tail*/
			if ( (geotype[0] < SP) && (geotype[1] < SP) ) {
                currangle = acos(DOT(conf->particle[num1].patchdir[0],conf->particle[num2].patchdir[0]) - DOT(conf->particle[num1].dir,conf->particle[num2].patchdir[0])  );
                energy += harmonic(currangle,topo->chainparam[conf->particle[num1].chaint].angle2eq,topo->chainparam[conf->particle[num1].chaint].angle2c);
            } else {
				energy += 0.0;
            }
		} else {
			if (num2 == topo->conlist[num1][1]) {
				/*num1 is connected to num2 by head*/
				if ( (geotype[0] < SP) && (geotype[1] < SP) ) {
                    currangle = acos(DOT(conf->particle[num2].patchdir[0],conf->particle[num1].patchdir[0]) - DOT(conf->particle[num2].dir,conf->particle[num1].patchdir[0])  );
                    energy += harmonic(currangle,topo->chainparam[conf->particle[num2].chaint].angle2eq,topo->chainparam[conf->particle[num2].chaint].angle2c);
                } else {
                    energy += 0.0;
                }
			}
		}
	}
	//    printf("angleener: %f\n",energy);
	return energy;
}




/*   cluses distance calculation*/

void closestdist(struct interacts * interact)
{
	double c, d, halfl;

    struct vector mindist_segments(struct vector dir1, double halfl1,
                                   struct vector dir2, double halfl2, struct vector r_cm);
	double linemin(double, double);

	
	//printf("we have %d %d  ",interact->param->geotype[0],interact->param->geotype[1] );
	if ((interact->param->geotype[0] >= SP) && (interact->param->geotype[1] >= SP)) { /*we have two spheres - most common, do nothing*/
		//printf("we have two spheres ");
		interact->distvec = interact->r_cm;
		interact->dist = sqrt(interact->dotrcm);
		interact->distcm = interact->dist;
	} else {
		if ((interact->param->geotype[0] < SP) && (interact->param->geotype[1] < SP)) { /*we have two spherocylinders*/
			interact->distvec = mindist_segments(interact->part1->dir,interact->param->half_len[0],
                                                 interact->part2->dir, interact->param->half_len[1], interact->r_cm);
			interact->dist=sqrt(DOT(interact->distvec,interact->distvec));
		} else {
			if (interact->param->geotype[0] < SP) { /*We have one spherocylinder -it is first one*/
				halfl=interact->param->half_len[0];/*finding closest vector from sphyrocylinder to sphere*/
				c = DOT(interact->part1->dir,interact->r_cm);
				if (c >= halfl) d = halfl;
				else {
					if (c > -halfl) d = c;
					else d = -halfl;
				}
				interact->contt = c;
				interact->distvec.x = - interact->r_cm.x + interact->part1->dir.x * d;
				interact->distvec.y = - interact->r_cm.y + interact->part1->dir.y * d;
				interact->distvec.z = - interact->r_cm.z + interact->part1->dir.z * d;
				interact->dist=sqrt(DOT(interact->distvec,interact->distvec));
			} else { /*lst option first one is sphere second one spherocylinder*/
				halfl=interact->param->half_len[1]; /*finding closest vector from sphyrocylinder to sphere*/
				c = DOT(interact->part2->dir,interact->r_cm); 
				if (c >= halfl) d = halfl;
				else {
					if (c > -halfl) d = c;
					else d = -halfl;
				}
				interact->contt = -c;
				interact->distvec.x = interact->r_cm.x - interact->part2->dir.x * d;
				interact->distvec.y = interact->r_cm.y - interact->part2->dir.y * d;
				interact->distvec.z = interact->r_cm.z - interact->part2->dir.z * d;
				interact->dist=sqrt(DOT(interact->distvec,interact->distvec));
			}
		}
	}

}


/*..............................................................................*/
/*
   Determines energy of two particles
 */

double paire(long num1, long num2, double (* intfce[MAXT][MAXT])(struct interacts *), 
		struct topo * topo, struct conf * conf)
{
	double energy=0.0;  /* energy*/
	struct vector r_cm;  /* Vector between centres of mass from part2 to part1*/
	struct interacts interact; /*interaction parameters*/

	double bondenergy(long, long, struct interacts *, struct topo * topo, struct conf * conf);
	double angleenergy(long, long, struct interacts *, struct topo * topo, struct conf * conf);

	/*Placing interactin particle in unit box and finding vector connecting CM*/
	/*r_cm = image(part1.pos, part2.pos, box); explicit statement below for performance optimization*/
	r_cm.x = conf->particle[num1].pos.x - conf->particle[num2].pos.x;
	r_cm.y = conf->particle[num1].pos.y - conf->particle[num2].pos.y;
	r_cm.z = conf->particle[num1].pos.z - conf->particle[num2].pos.z;
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

	interact.dotrcm = DOT(r_cm,r_cm);
	if ( interact.dotrcm > topo->sqmaxcut) return 0.0;  /* distance so far that even spherocylinders cannot be within cutoff  */
	interact.r_cm=r_cm;
	interact.contt = 0;
	interact.distvec.x = 0;
	interact.distvec.y = 0;
	interact.distvec.z = 0;
	
	interact.box = conf->box;
	interact.part1 = &conf->particle[num1];
	interact.part2 = &conf->particle[num2];
	interact.param = topo->ia_params[conf->particle[num1].type] + conf->particle[num2].type;
	
	if(intfce[conf->particle[num1].type][conf->particle[num2].type] == NULL){
		fprintf(stderr, "interaction function for type %d and %d not defined!\n",
				conf->particle[num1].type, conf->particle[num2].type);
	}
	energy = (*intfce[conf->particle[num1].type][conf->particle[num2].type])( &interact);

	//printf("num: %ld  %ld   e: %f dist: %f",num1,num2,energy,interact.dist);
	energy += bondenergy ( num1, num2, &interact, topo, conf);
	energy += angleenergy ( num1, num2, &interact, topo, conf);

	//printf("  e: %f\n",energy);
	return energy;
}


/*...........................................................................*/

/*Calculates interaction of target particle and external field version 2
  calculate projection of spherocylinder in direction of patch and calculate 
  interacting line segment within cutoff
 */

double extere2 (long target, struct topo * topo, struct conf * conf)
{
	double repenergy=0.0,atrenergy=0.0;  /* energy*/
	double rcmz;                         /* z distance between*/
	double ndist;                        /* distance for CM of interacting line segment*/
	double interendz;                    /* z coordinate of interaction end*/
	struct interacts interact;           /* interaction parameters*/
	double orient; 
    double halfl;
	BOOL positive, orientin;  
	struct vector olddir;  
	struct vector project; /*vector for projection down to plane */
	
	double erepulsive(struct interacts *);
//	struct vector vec_perpproject(struct vector*, struct vector*);
//	void normalise(struct vector *);
	double fanglscale(double, struct ia_param *, int which);
	void exter2_closestdist(struct interacts * interact, BOOL *positive, BOOL *orientin, double *orient, 
                            double *rcmz,double *interendz, struct vector *project);
	double exter2_atre(struct interacts * interact,int *orientin, double *rcmz, double *interendz, BOOL *positive, 
                       double orient,struct vector *project, double *ndist,int, double );
	
	/* calcualte distance to center of mass*/
	if (  conf->particle[target].pos.z < 0  ) {
		rcmz = conf->box.z * (conf->particle[target].pos.z - (double)( (long)(conf->particle[target].pos.z - 0.5) ) );
		
	} else {
		rcmz = conf->box.z * (conf->particle[target].pos.z - (double)( (long)(conf->particle[target].pos.z + 0.5) ) );
		
	}
	project.x=0;
	project.y=0;
	if (rcmz < 0) {
		interact.dist = -rcmz;
		positive = FALSE;
		interendz = -1.0;
		project.z = 1.0;
	} else {
		interact.dist = rcmz;
		positive = TRUE;
		interendz = 1.0;
		project.z = -1.0;
	}
	
	interact.dotrcm = rcmz * rcmz;
	if ( interact.dotrcm > topo->exter.sqmaxcut) return 0.0;  /* distance so far that even spherocylinders cannot be within cutoff  */
	interact.distvec.z = interact.r_cm.z;
	interact.distcm = interact.dist;
	interact.box = conf->box;
	interact.part1 = &conf->particle[target];
	interact.param = &topo->exter.interactions[conf->particle[target].type];
	halfl = 0.5* topo->exter.interactions[conf->particle[target].type].len[0];
	ndist = interact.dist;
	orientin = TRUE;
	orient = 0.0;
	
	
	exter2_closestdist(&interact,&positive,&orientin,&orient,&rcmz,&interendz,&project);
	/* now we have closest distance so we can calculate repulsion*/
	repenergy = erepulsive(&interact);
	//printf("dist: %f",interact.dist);
	
	
	/*save chiral stuff*/
	olddir = interact.part1->dir;
	if ((interact.param->geotype[0] == CHCPSC)||(interact.param->geotype[0] == CHPSC)) {
		interact.part1->dir = interact.part1->chdir[0];
		exter2_closestdist(&interact,&positive,&orientin,&orient,&rcmz,&interendz,&project);
	}
	if (( interact.dist > interact.param->rcut ) || (interact.param->epsilon == 0.0 ) || 
		( (interact.part1->patchdir[0].z >0)&&(positive) ) || ( (interact.part1->patchdir[0].z <0)&&(!(positive)) )  )
		atrenergy = 0.0;
	else {
		atrenergy = exter2_atre(&interact,&orientin,&rcmz,&interendz,&positive,orient,&project,&ndist,0,halfl);
	}
	if ((interact.param->geotype[0] == TCPSC)||(interact.param->geotype[0] == TPSC)||
	  (interact.param->geotype[0] == TCHCPSC)||(interact.param->geotype[0] == TCHPSC)) {
		if ((interact.param->geotype[0] == TCHCPSC)||(interact.param->geotype[0] == TCHPSC)) {
		    interact.part1->dir = interact.part1->chdir[1];
		    exter2_closestdist(&interact,&positive,&orientin,&orient,&rcmz,&interendz,&project);
		}
		exter2_closestdist(&interact,&positive,&orientin,&orient,&rcmz,&interendz,&project);
		if (( interact.dist > interact.param->rcut ) || (interact.param->epsilon == 0.0 ) || 
		      ( (interact.part1->patchdir[1].z >0)&&(positive) ) || ( (interact.part1->patchdir[1].z <0)&&(!(positive)) )  )
		    atrenergy += 0.0;
		else {
		    atrenergy += exter2_atre(&interact,&orientin,&rcmz,&interendz,&positive,orient,&project,&ndist,1,halfl);
		}
		
	}
	
	if ((interact.param->geotype[0] == CHCPSC)||(interact.param->geotype[0] == CHPSC)||
	  (interact.param->geotype[0] == TCHCPSC)||(interact.param->geotype[0] == TCHPSC) ) {
		interact.part1->dir = olddir;
	}

	//printf("%f   %f \n",conf->particle[target].pos.z*conf->box.z,repenergy+atrenergy);
	return repenergy+atrenergy;
}

double exter2_atre(struct interacts * interact,int *orientin, double *rcmz, double *interendz, BOOL *positive, double orient,struct vector *project, double *ndist,int numpatch, double halfl)
{
	struct vector pbeg,pend;             /* projected spherocylinder begining and end*/
	double a,length1,length2, f0,f1;
	struct vector cm1,cm2;               /* centrar of interacting segments */	
	int line;
	struct vector partbeg,partend;        /*closest and furthest point of particle*/
	struct vector inters;
	double atrenergy=0.0;
	
	int cpsc_wall(struct vector* pbeg, struct vector* pend, struct vector* projectdir,  
		struct vector* partdir, double* halfl,BOOL* orientin,BOOL* positive, double* rcmz,
		 double * cut, struct vector* partbeg, struct vector* partend);
	int psc_wall(struct vector* pbeg, struct vector* pend, struct vector* projectdir,  
		struct vector* partdir, BOOL* positive, double * cut, struct vector* partbeg, struct vector* partend);
	
	/*interaction with PATCHY SPHEROCYLINDERS*/
	if ((interact->param->geotype[0] < SP)&&(interact->param->geotype[0] > SCA)) {
		//printf("partdir: %f %f %f \n ",interact->part1->dir.x,interact->part1->dir.y,interact->part1->dir.z);
		//printf("patchdir: %f %f %f \n ",interact->part1->patchdir[0].x,interact->part1->patchdir[0].y,interact->part1->patchdir[0].z);
		/* calculate position of closest and furthest point (begining and end of spherocylinder)*/
		a = (*orientin-0.5)*2;  /*if orientin a =1 else a=-1 */
		partbeg.x =  a * interact->part1->dir.x * halfl;
		partbeg.y =  a * interact->part1->dir.y * halfl;
		partbeg.z = *rcmz + a * interact->part1->dir.z *halfl;
		partend.x =  - a * interact->part1->dir.x * halfl;
		partend.y =  - a * interact->part1->dir.y * halfl;
		partend.z = *rcmz - a * interact->part1->dir.z * halfl;
			
			
		//printf("partbeg %f %f %f partend %f %f %f  \n",partbeg.x,partbeg.y,partbeg.z,partend.x,partend.y,partend.z);			
		/*calculate interacting line segment and its cm of spherocylinder*/
		/*calculate end point z*/
		if ( (interact->param->rcut - interact->dist)/fabs(interact->part1->dir.z) < 2.0*halfl ){
			/*if cutoff goes through spherocylinder the end point is at cutoff*/
			*interendz *= interact->param->rcut;
		} else {
			/*endpoint is at the end of spherocylinders*/
			*interendz = partend.z;
		}
		/*calculate CM of interacting line segment of spherocylinder*/
		if (*positive) {
			cm1.z = AVER(*interendz,interact->dist);
		} else {
			cm1.z = AVER(*interendz,-interact->dist);
		}
		if (interact->part1->dir.z != 0.0 ) {
			a = (*interendz - cm1.z ) / interact->part1->dir.z;
			length1= -orient*2.0*a;
			a = a + orient*halfl;
		} else {
			a = 0.0;
			length1 = 2.0*halfl;
		}
		//printf("len1: %f rcm %f interz %f cutoff %f  \n",length1,rcmz, interendz,interact.dist);
		cm1.x = interact->part1->dir.x * a;
		cm1.y = interact->part1->dir.y * a;
			
		/* we have interacting segment*/
		if ((interact->param->geotype[0] == CPSC)||(interact->param->geotype[0] == CHCPSC)) {
			/*CPSC type*/
			if ( ((*interendz >= interact->dist)&&(*positive)) || ((*interendz <= -interact->dist)&&(!(*positive))) ){ 
				/*test if projection is not all out of interaction*/
				line = cpsc_wall(&pbeg,&pend,project,&interact->part1->dir, \
					 &interact->param->half_len[0],orientin,positive,rcmz,&interact->param->rcut,&partbeg,&partend);
				//printf("line: %d beg %f %f end %f %f  \n",line,pbeg.x,pbeg.y,pend.x,pend.y);	
			} else {
				line = 0;
			}
		} else {
			/*PSC and CHPSC interaction with wall */
			line = psc_wall(&pbeg,&pend,project,&interact->part1->dir, \
				positive,&interact->param->rcut,&partbeg,&partend);	
			//printf("line: %d beg %f %f end %f %f  \n",line,pbeg.x,pbeg.y,pend.x,pend.y);
		}		
			
		if (line > 0) {
			/*cm2 by average begining and end*/
			cm2.x = AVER(pbeg.x,pend.x);
			cm2.y = AVER(pbeg.y,pend.y);
			cm2.z = 0.0;
				
			/*length by size of end-benining*/
			length2 = sqrt( (pend.x-pbeg.x)*(pend.x-pbeg.x)+(pend.y-pbeg.y)*(pend.y-pbeg.y) );
				
			inters.x = cm2.x - cm1.x;
			inters.y = cm2.y - cm1.y;
			inters.z = cm2.z - cm1.z;
			//printf("cm2 %f %f %f inters %f %f %f \n",cm2.x,cm2.y,cm2.z,inters.x,inters.y,inters.z);
			*ndist = sqrt(DOT(inters,inters));
			if (*ndist < interact->param->pdis) { 
				atrenergy = -interact->param->epsilon;
			}
			else {
				atrenergy= cos(PIH*(*ndist-interact->param->pdis)/interact->param->pswitch);
				atrenergy *= -atrenergy*interact->param->epsilon;
			}
			/*  scaling function1: dependence on the length of intersetions plus*/
			f0=(length1 + length2)*0.5;
			/*scaling with angle*/
			f1 = fabs(interact->part1->patchdir[numpatch].z);
			atrenergy *= f0*f1;
				
			//printf(" %f   %f    %f %f %f %f %f \n",conf->particle[target].pos.z*conf->box.z,atrenergy, area, length1, length2,f0,ndist);
			//printf("%f %f %f %f\n",pbeg.x,pbeg.y,pend.x,pend.y);
		} else {
			atrenergy = 0.0;
		}
			
	} else {
		if (*ndist < interact->param->pdis) 
			atrenergy = -interact->param->epsilon;
		else  {
			atrenergy= cos(PIH*(*ndist-interact->param->pdis)/interact->param->pswitch);
			atrenergy *= -atrenergy*interact->param->epsilon;
		}
		/*add wall scaling wall area/ particle arear.. to reflect that we have a wall not sphere */
		atrenergy *= (interact->param->rcut*interact->param->rcut - (*ndist)*(*ndist))/(interact->param->sigma*interact->param->sigma) ;
	}
	
	return atrenergy;
}

void exter2_closestdist(struct interacts * interact, BOOL *positive, BOOL *orientin, double *orient, double *rcmz,double *interendz, struct vector *project)
{
	if (*rcmz < 0) {
		interact->dist = -(*rcmz);
		*positive = FALSE;
		*interendz = -1.0;
		project->z = 1.0;
	} else {
		interact->dist = (*rcmz);
		*positive = TRUE;
		*interendz = 1.0;
		project->z = -1.0;
	}
	/*psc closest is allways end closer to wall*/
	if (interact->param->geotype[0] < SP ){
		/*calculate closest point distance*/
		if (interact->part1->dir.z > 0) {
			if (*positive) {
				*orientin = FALSE;
				*orient = -1.0;
				interact->dist = *rcmz -interact->part1->dir.z * interact->param->half_len[0];
			} else {
				*orientin = TRUE;
				*orient = 1.0;
				interact->dist = -( *rcmz + interact->part1->dir.z * interact->param->half_len[0]);
			}
		} else {
			if (*positive) {
				*orientin = TRUE;
				*orient = 1.0;
				interact->dist = *rcmz + interact->part1->dir.z * interact->param->half_len[0];
			} else {
				*orientin = FALSE;
				*orient = -1.0;
				interact->dist = -( *rcmz -interact->part1->dir.z * interact->param->half_len[0]);
			}
		}
	}
}


/*...........................................................................*/

/*Calculates interaction of target particle and external field
  calculate projection of patch of spherocylinder on wall 
  evaluate intersection area and calculate interaction from that
 */

double exter_atre(struct interacts * interact,int *orientin, double *rcmz, double *interendz, BOOL *positive, double orient,struct vector *project, double *ndist,int numpatch,double halfl)
{
	double area,a,b,c,r2;
	double atrenergy=0.0;  /* energy*/
	BOOL countend;
	struct vector cm1,cm2;               /* centrar of interacting segments */
	struct vector pbeg,pend;             /* projected spherocylinder begining and end*/
	struct vector inters,newdir;  
	struct vector pbeg1,pend1,pbeg2,pend2,pextr1,pextr2,pextr3,pextr4; /*additinal point of projected patch for calculation of area */
	double length1, cuttoproject, f0;
	int line, line1, line2,extra;
	struct vector partbeg,partend;        /*closest and furthest point of particle*/

	double erepulsive(struct interacts *);
	struct vector vec_perpproject(struct vector*, struct vector*);
	void normalise(struct vector *);
	double fanglscale(double, struct ia_param *, int which);
	struct vector vec_create(double, double, double);
	double areaeightpoints(struct vector*,struct vector*,struct vector*,struct vector*,struct vector*,struct vector*,struct vector*,struct vector*);
	int cpsc_wall(struct vector* pbeg, struct vector* pend, struct vector* projectdir,  
			struct vector* partdir, double* halfl,BOOL* orientin,BOOL* positive, double* rcmz,
			double * cut, struct vector* partbeg, struct vector* partend);
	int psc_wall(struct vector* pbeg, struct vector* pend, struct vector* projectdir,  
			struct vector* partdir, BOOL* positive, double * cut, struct vector* partbeg, struct vector* partend);
	int cutprojectatwall(struct vector* pextr1, struct vector* pextr2, struct vector* pextr3, struct vector* pextr4, 
						 struct vector* projectdir, struct vector* partdir, double * cutdist, struct vector *partbeg, 
						 struct vector *partend, struct vector *pend, double *cuttoproject, BOOL* orientin);
	void exter2_closestdist(struct interacts * interact, BOOL *positive, BOOL *orientin, double *orient, double *rcmz,double *interendz, struct vector *project);
	

	/*interaction with PATCHY SPHEROCYLINDERS*/
	if ((interact->param->geotype[0] < SP)&&(interact->param->geotype[0] > SCA)) {
		//printf("partdir: %f %f %f \n ",interact->part1->dir.x,interact->part1->dir.y,interact->part1->dir.z);
		//printf("patchdir: %f %f %f \n ",interact->part1->patchdir[numpatch].x,interact->part1->patchdir[numpatch].y,interact->part1->patchdir[numpatch].z);
		/* calculate position of closest and furthest point (begining and end of spherocylinder)*/
		a = (*orientin-0.5)*2;  /*if orientin a =1 else a=-1 */
		partbeg.x =  a * interact->part1->dir.x * halfl;
		partbeg.y =  a * interact->part1->dir.y * halfl;
		partbeg.z = *rcmz + a * interact->part1->dir.z * halfl;
		partend.x =  - a * interact->part1->dir.x * halfl;
		partend.y =  - a * interact->part1->dir.y * halfl;
		partend.z = *rcmz - a * interact->part1->dir.z * halfl;
			
		//printf("partbeg %f %f %f partend %f %f %f  \n",partbeg.x,partbeg.y,partbeg.z,partend.x,partend.y,partend.z);			
		/*calculate interacting line segment and its cm of spherocylinder*/
		/*calculate end point z*/
		if ( (interact->param->rcut - interact->dist)/fabs(interact->part1->dir.z) < halfl*2.0 ){
			/*if cutoff goes through spherocylinder the end point is at cutoff*/
			*interendz *= interact->param->rcut;
		} else {
			/*endpoint is at the end of spherocylinders*/
			*interendz = partend.z;
		}
		/*calculate CM of interacting line segment of spherocylinder*/
		if (*positive) {
			cm1.z = AVER(*interendz,interact->dist);
		} else {
			cm1.z = AVER(*interendz,-interact->dist);
		}
		if (interact->part1->dir.z != 0.0 ) {
		    a = (*interendz - cm1.z ) / interact->part1->dir.z;
			length1= -orient*2.0*a;
			a = a + orient*halfl;
		} else {
		    a = 0.0;
			length1 = 2.0*halfl;
		}
		//printf("len1: %f rcm %f interz %f cutoff %f  \n",length1,rcmz, interendz,interact->dist);
		cm1.x = interact->part1->dir.x * a;
		cm1.y = interact->part1->dir.y * a;
		/*calculate projection on wall as infinite line and make it interacting segment*/
		if (interact->part1->patchdir[numpatch].z != 0) {
			cuttoproject = -interact->param->rcut*interact->part1->patchdir[numpatch].z; /*z coordinate of point where projection is in cut distance*/
			if ( ((partend.z < cuttoproject)&&(*positive)) || ((cuttoproject < partend.z)&&(!(*positive))) ){ 
				cuttoproject = partend.z;
			}
		} else {
			cuttoproject = partbeg.z;
		}	
		//printf("cutproject %f \n",cuttoproject);
		//printf("cm1 %f %f %f \n",cm1.x, cm1.y,cm1.z );		
			
		/* we have interacting segment*/
		if ((interact->param->geotype[0] == CPSC)||(interact->param->geotype[0] == CHCPSC)) {
			/*CPSC type*/
			if ( ((cuttoproject >= interact->dist)&&(*positive)) || ((cuttoproject <= -interact->dist)&&(!(*positive))) ){ 
				/*test if projection is not all out of interaction*/
				line = cpsc_wall(&pbeg,&pend,&interact->part1->patchdir[numpatch],&interact->part1->dir, \
						&interact->param->half_len[0],orientin,positive,rcmz,&interact->param->rcut,&partbeg,&partend);
				//printf("line: %d beg %f %f end %f %f  \n",line,pbeg.x,pbeg.y,pend.x,pend.y);	
			} else {
				line = 0;
			}
		} else {
			/*PSC and CHPSC interaction with wall */
			line = psc_wall(&pbeg,&pend,&interact->part1->patchdir[numpatch],&interact->part1->dir, \
					positive,&interact->param->rcut,&partbeg,&partend);	
			//printf("line: %d beg %f %f end %f %f  \n",line,pbeg.x,pbeg.y,pend.x,pend.y);
		}		
		
		if (line > 0) {
			area = 0.0;	    
			/*project cutoff boudaries*/
			if (line == 2 ) {
				/*if projection end is on sphere of begining don't care about cylinder cutoff*/
				extra = 0;
			} else {
				extra = cutprojectatwall(&pextr1, &pextr2, &pextr3, &pextr4, &interact->part1->patchdir[numpatch],  \
					 &interact->part1->dir, &interact->param->rcut, &partbeg, &partend,&pend,&cuttoproject,orientin);
			}
			//printf("extr1: %d %f %f extr2 %f %f extr3 %f %f extr4 %f %f \n",extra,pextr1.x,pextr1.y,pextr2.x,pextr2.y,pextr3.x,pextr3.y,pextr4.x,pextr4.y);
					
			/*project patch boundaries on the first side*/
			newdir=interact->part1->patchsides[0+2*numpatch];
			line1 = cpsc_wall(&pbeg1,&pend1,&newdir,&interact->part1->dir, \
				&interact->param->half_len[0],orientin,positive,rcmz,&interact->param->rcut,&partbeg,&partend);
			if ( ((interact->param->geotype[0] == PSC)||(interact->param->geotype[0] == CHPSC)) ) {
				line1 = psc_wall(&pbeg1,&pend1,&newdir,&interact->part1->dir, \
					  positive,&interact->param->rcut,&partbeg,&partend);
			}
			//printf("line1: %d beg1 %f %f end1 %f %f  \n",line1,pbeg1.x,pbeg1.y,pend1.x,pend1.y);
			
			/*project patch boundaries on the second side*/
			newdir=interact->part1->patchsides[1+2*numpatch];
			line2 = cpsc_wall(&pbeg2,&pend2,&newdir,&interact->part1->dir, \
				&interact->param->half_len[0],orientin,positive,rcmz,&interact->param->rcut,&partbeg,&partend);
			if ( ((interact->param->geotype[0] == PSC)||(interact->param->geotype[0] == CHPSC)) ) {
				line2 = psc_wall(&pbeg2,&pend2,&newdir,&interact->part1->dir, \
					 positive,&interact->param->rcut,&partbeg,&partend);
			}
			//printf("line2: %d beg2 %f %f end2 %f %f \n",line2,pbeg2.x,pbeg2.y,pend2.x,pend2.y);
			
			/*calculate area*/
			if (extra == 0) {
				/*thish should only happen when there is PSC interacting only with end*/
				if (line1 == 0) {
					if (line2==0) {
						/*circle around middle-pbeg*/
						area = PI*( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y));
					}
					else{
						/* circle around middle-pbeg minus circle segment*/
						a = AVER(pbeg2.x,pend2.x);
						b = AVER(pbeg2.y,pend2.y);
						c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
						r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); /*radius squared*/
						area =  r2*(PI-acos(sqrt(c/r2))) + sqrt(r2*c-c*c);
					}	
				} else {
					if (line2==0) {
						/* circle around middle-pbeg minus circle segment*/
						a = AVER(pbeg1.x,pend1.x);
						b = AVER(pbeg1.y,pend1.y);
						c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
						r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); /*radius squared*/
						area =  r2*(PI-acos(sqrt(c/r2))) + sqrt(r2*c-c*c);
					} else {
						//area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pend,&pend1,&pbeg1,NULL,NULL); /* B B2 E2 E E1 B1 */
						/*circle minus two circle segments*/
						r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); /*radius squared*/
						area = r2*PI;
						a = AVER(pbeg1.x,pend1.x);
						b = AVER(pbeg1.y,pend1.y);
						c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
						area +=  -r2*acos(sqrt(c/r2)) + sqrt(r2*c-c*c);
						a = AVER(pbeg2.x,pend2.x);
						b = AVER(pbeg2.y,pend2.y);
						c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
						area +=  -r2*acos(sqrt(c/r2)) + sqrt(r2*c-c*c);
					}
				}
	
			} else {
				b = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend.y)- (pextr2.x-pend.x)*(pextr4.y-pextr2.y));/*pend on 42*/ 
				c = fabs((pextr1.x-pextr3.x)*(pextr3.y-pend.y)- (pextr3.x-pend.x)*(pextr1.y-pextr3.y));/*pend on 13*/
				if ( ( b< ZEROTOL) || ( c< ZEROTOL) ) 
					countend = FALSE;
				else 
					countend = TRUE;
				if (line1 == 0) {
					if (line2 == 0) {
						if ( countend ) {
							area = areaeightpoints(&pbeg,&pextr2,&pextr4,&pend,&pextr3,&pextr1,NULL,NULL);/* B 2 4 E 3 1 */
						} else 
							area = areaeightpoints(&pbeg,&pextr2,&pextr4,&pextr3,&pextr1,NULL,NULL,NULL);/* B 2 4 3 1 */
					} else {
						a = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend2.y)- (pextr2.x-pend2.x)*(pextr4.y-pextr2.y));
						if ( a< ZEROTOL) /*pend2 on 42*/ { 
							if ( countend ) {
								area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pextr1,NULL); /* B B2 E2 4 E 3 1 */
							} else {
								area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pextr1,NULL,NULL); /* B B2 E2 4 3 1 */
							}
						} else { 
							a = fabs((pextr1.x-pextr3.x)*(pextr3.y-pend2.y)- (pextr3.x-pend2.x)*(pextr1.y-pextr3.y));
							if ( a< ZEROTOL) /*pend2 on 13*/ {
								area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr1,NULL,NULL,NULL,NULL); /* B B2 E2 1 */
							} else { /*pend2 on 34 or on begining sphere  of psc*/
								if (line2 == 2) {
									if ( countend ) {
										area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pextr1,NULL); /* B B2 E2 4 E 3 1 */
									} else {
										area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pextr1,NULL,NULL); /* B B2 E2 4 3 1 */
									}
								} else {
									if ( countend ) {
										area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pend,&pextr3,&pextr1,NULL,NULL); /* B B2 E2 E 3 1 */
									} else {
										area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr3,&pextr1,NULL,NULL,NULL); /* B B2 E2 3 1 */
									}
								}
							}
						}
					}
				} else { 
					a = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend1.y)- (pextr2.x-pend1.x)*(pextr4.y-pextr2.y));
					if ( a< ZEROTOL) /*pend1 on 42*/ { 
						if (line2 == 0) {
							area = areaeightpoints(&pbeg,&pextr2,&pend1,&pbeg1,NULL,NULL,NULL,NULL); /* B 2 E1 B1 */
						} else {
							area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pend1,&pbeg1,NULL,NULL,NULL); /* B B2 E2 E1 B1 */
						}
					} else { 
						a = fabs((pextr1.x-pextr3.x)*(pextr3.y-pend1.y)- (pextr3.x-pend1.x)*(pextr1.y-pextr3.y));
						if ( a< ZEROTOL) /*pend1 on 13*/ {
							if (line2 == 0) {
								if (countend) {
									area =  areaeightpoints(&pbeg,&pextr2,&pextr4,&pend,&pextr3,&pend1,&pbeg1,NULL);  /* B 2 4 E 3 E1 B1 */
								} else {
									area =  areaeightpoints(&pbeg,&pextr2,&pextr4,&pextr3,&pend1,&pbeg1,NULL,NULL);  /* B 2 4 3 E1 B1 */
								}
							} else {
								a = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend2.y)- (pextr2.x-pend2.x)*(pextr4.y-pextr2.y));
								if ( a< ZEROTOL) /*pend2 on 42*/ { 
									if (countend)
										area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pend1,&pbeg1); /* B B2 E2 4 E 3 E1 B1 */
									else 
										area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 4 3 E1 B1 */
								} else {
									a = fabs((pextr3.x-pextr1.x)*(pextr1.y-pend2.y)- (pextr1.x-pend2.x)*(pextr3.y-pextr1.y));
									if ( a< ZEROTOL) /*pend2 on 31*/ {
										area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pend1,&pbeg1,NULL,NULL,NULL); /* B B2 E2 E1 B1 */
									} else { /*pend2 close to 34 or on begining sphere  of psc*/
										if (line2 == 2) {
											if (countend)
												area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pend1,&pbeg1); /* B B2 E2 4 E 3 E1 B1 */
											else
												area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 4 3 E1 B1 */
										} else {
											if (countend)
												area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pend,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 E 3 E1 B1 */
											else
												area = areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr3,&pend1,&pbeg1,NULL,NULL); /* B B2 E2 3 E1 B1 */
										}
									}
								}
							}
						} else {/*pend1 close to 34 or on beging sphere for psc*/
							if (line2 == 0) {
								if (line1 ==2) {
									if (countend)
										area = areaeightpoints(&pbeg,&pextr2,&pextr4,&pend,&pextr3,&pend1,&pbeg1,NULL); /* B 2 4 E 3 E1 B1*/
									else {
										area = areaeightpoints(&pbeg,&pextr2,&pextr4,&pextr3,&pend1,&pbeg1,NULL,NULL); /* B 2 4 3 E1 B1*/
									}
								} else {
									if (countend)
										area = areaeightpoints(&pbeg,&pextr2,&pextr4,&pend,&pend1,&pbeg1,NULL,NULL); /* B 2 4 E E1 B1*/
									else {
										area = areaeightpoints(&pbeg,&pextr2,&pextr4,&pend1,&pbeg1,NULL,NULL,NULL); /* B 2 4 E1 B1*/
									}
								}
		
							} else { 
								a = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend2.y)- (pextr2.x-pend2.x)*(pextr4.y-pextr2.y));
								if ( a< ZEROTOL) /* pend2 on 42 */ { 
									if (countend)
										area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pend1,&pbeg1,NULL);  /* B B2 E2 4 E E1 B1 */
									else
										area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend1,&pbeg1,NULL,NULL);  /* B B2 E2 4 E1 B1 */
								} else { /*pend1 and pend2 close to 34 or on beging sphere for psc*/
									if (line2 == 2) {
										if (line1 == 2) {
											if (countend)
												area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pend1,&pbeg1); /* B B2 E2 4 E 3 E1 B1 */
											else
												area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 4 3 E1 B1 */
										} else {
											if (countend)
												area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pend1,&pbeg1,NULL); /* B B2 E2 4 E E1 B1 */
											else
												area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend1,&pbeg1,NULL,NULL); /* B B2 E2 4 E1 B1 */
										}
										
										
									} else {
										if (line1 == 2) {
											if (countend)
												area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pend,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 E 3 E1 B1 */
											else
												area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pextr3,&pend1,&pbeg1,NULL,NULL); /* B B2 E2 3 E1 B1 */
										} else {
											if (countend)
												area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pend,&pend1,&pbeg1,NULL,NULL); /* B B2 E2 E E1 B1 */
											else
												area =  areaeightpoints(&pbeg,&pbeg2,&pend2,&pend1,&pbeg1,NULL,NULL,NULL); /* B B2 E2 E1 B1 */
										}
									}
								}
							}
						}
					}
				} /*extra != 0*/
				if ((interact->param->geotype[0] == PSC)||(interact->param->geotype[0] == CHPSC)) {
					if (line1==2) {
						/* add circle segment*/
						a = AVER(pextr1.x,pend1.x); /*end to cutoff - pextr1 ,pend1 */
						b = AVER(pextr1.y,pend1.y);
						c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
						r2 = ( (partbeg.x-pend1.x)*(partbeg.x-pend1.x) + (partbeg.y-pend1.y)*(partbeg.y-pend1.y)); /*radius squared*/
						area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
						a = AVER(pbeg.x,pbeg1.x); /* between beginings - pbeg ,pbeg1 */
						b = AVER(pbeg.y,pbeg1.y);
						c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
						area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
					} else {
						if (line1==0) {
							/* add circle segment*/
							a = AVER(pextr1.x,pbeg.x); /* begining to cutoff*/
							b = AVER(pextr1.y,pbeg.y);
							c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
							r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); /*radius squared*/
							area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
						}
					}
						if (line2==2) {
						/* add circle segment*/
						a = AVER(pextr3.x,pend2.x); /*end to cutoff - pextr3 ,pend2 */
						b = AVER(pextr3.y,pend2.y);
						c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
						r2 = ( (partbeg.x-pend2.x)*(partbeg.x-pend2.x) + (partbeg.y-pend2.y)*(partbeg.y-pend2.y)); /*radius squared*/
						area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
						a = AVER(pbeg.x,pbeg2.x); /* between beginings - pbeg ,pbeg2 */
						b = AVER(pbeg.y,pbeg2.y);
						c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
						area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
					} else {
						if (line2==0) {
							/* add circle segment*/
							a = AVER(pextr3.x,pbeg.x); /* begining to cutoff*/
							b = AVER(pextr3.y,pbeg.y);
							c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
							r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); /*radius squared*/
							area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
						}
					}
				}
			}
			/*area finished*/

			/*cm2 by average begining and end*/
			cm2.x = AVER(pbeg.x,pend.x);
			cm2.y = AVER(pbeg.y,pend.y);
			cm2.z = 0.0;
			/*length by size of end-benining*/
			//length2 = sqrt( (pend.x-pbeg.x)*(pend.x-pbeg.x)+(pend.y-pbeg.y)*(pend.y-pbeg.y) );
			inters.x = cm2.x - cm1.x;
			inters.y = cm2.y - cm1.y;
			inters.z = cm2.z - cm1.z;
			//printf("cm2 %f %f %f inters %f %f %f \n",cm2.x,cm2.y,cm2.z,inters.x,inters.y,inters.z);
			*ndist = sqrt(DOT(inters,inters));
			if (*ndist < interact->param->pdis) { 
				atrenergy = -interact->param->epsilon;
			}
			else {
				atrenergy= cos(PIH*(*ndist-interact->param->pdis)/interact->param->pswitch);
				atrenergy *= -atrenergy*interact->param->epsilon;
			}
			/*  scaling function1: dependence on the length of intersetions plus SCALING WITH AREA*/
			f0=(length1 + area / interact->param->sigma)*0.5;
			atrenergy *= f0;
			//printf(" %f   %f    %f %f %f %f %f %d %d %d \n",conf->particle[target].pos.z*conf->box.z,atrenergy, area, length1, length2,f0,ndist,extra,line1,line2);
			//printf("%f %f %f %f\n",pbeg.x,pbeg.y,pend.x,pend.y);
			//printf("%f %f %f %f %f %f\n",pbeg2.x,pend2.y,pextr2.x,pextr2.y,pextr1.x,pextr1.y);
		} else {
			atrenergy = 0.0;
		}

	} else {
		if (*ndist < interact->param->pdis) 
			atrenergy = -interact->param->epsilon;
		else  {
			atrenergy= cos(PIH*(*ndist-interact->param->pdis)/interact->param->pswitch);
			atrenergy *= -atrenergy*interact->param->epsilon;
		}
		/*add wall scaling wall area/ particle arear.. to reflect that we have a wall not sphere */
		atrenergy *= (interact->param->rcut*interact->param->rcut - (*ndist)*(*ndist))/(interact->param->sigma*interact->param->sigma) ;
	}

	
	//printf("%f   %f \n",conf->particle[target].pos.z*conf->box.z,atrenergy);
	return atrenergy;
}



/*..............................................................................*/

/* Initializes the array with the pointers to the energy function
 */
void init_intfce(double (* intfce[MAXT][MAXT])(struct interacts *), struct topo * topo){
	// NB
	// Fill in the names of the functions for calculating the 
	// interaction energy
	long geotype, other_geotype;
	int i, j;
	for(i = 0; i < MAXT; i++){
		for(j = 0; j < MAXT; j++){
			/* Initialize them as not existing */
			intfce[i][j] = &enoexist;
			geotype = topo->ia_params[i][j].geotype[0];
			other_geotype = topo->ia_params[i][j].geotype[1];
			
			if ( ( (geotype == CHCPSC || geotype == CPSC || geotype == TCHCPSC || geotype == TCPSC) &&
				    (other_geotype == CHPSC || other_geotype == PSC || other_geotype == TCHPSC || other_geotype == TPSC) ) ||
			      ( (geotype == CHPSC || geotype == PSC || geotype == TCHPSC || geotype == TPSC)  &&
				    (other_geotype == CHCPSC || other_geotype == CPSC || other_geotype == TCHCPSC || other_geotype == TCPSC) ) )  {
				intfce[i][j] = &e_psc_cpsc;
			}
			if ( (geotype == CHCPSC || geotype == CPSC || geotype == TCHCPSC || geotype == TCPSC) && 
				    (other_geotype == CHCPSC || other_geotype == CPSC || other_geotype == TCHCPSC || other_geotype == TCPSC) ){
				intfce[i][j] = &e_cpsc_cpsc;
			}
			if ( (geotype == CHPSC || geotype == PSC || geotype == TCHPSC || geotype == TPSC) && 
				    (other_geotype == CHPSC || other_geotype == PSC || other_geotype == TCHPSC || other_geotype == TPSC) ){
				intfce[i][j] = &e_psc_psc;
			}
						
			if(geotype == SCN || geotype == SPN
					|| other_geotype == SCN || other_geotype == SPN){
				intfce[i][j] = &e_spn_or_scn;
			}
			if((geotype == SCA && other_geotype == SCA) 
					|| (geotype == SPA && other_geotype == SPA)){
				intfce[i][j] = &e_2sca_or_2spa;
			}
			if((geotype == SCA && other_geotype == SPA) 
					|| (geotype == SPA && other_geotype == SCA)){
				intfce[i][j] = &e_spa_sca;
			}
			if(( (geotype == PSC || geotype == CHPSC || geotype == TCHPSC || geotype == TPSC) && other_geotype == SPA) 
					|| (geotype == SPA && (other_geotype == PSC||other_geotype == CHPSC || other_geotype == TCHPSC || other_geotype == TPSC) )){
				intfce[i][j] = &e_psc_spa;
			}
			if(( (geotype == CPSC ||geotype == CHCPSC || geotype == TCHCPSC || geotype == TCPSC) && other_geotype == SPA) 
					|| (geotype == SPA && (other_geotype == CPSC||other_geotype == CHCPSC || other_geotype == TCHCPSC || other_geotype == TCPSC)  )){
				intfce[i][j] = &e_cpsc_spa;
			}

		}
	}
}
/*..............................................................................*/

/*
   Compare energy change to temperature and based on Boltzmann probability
   return either 0 to accept or 1 to reject the move
 */

int movetry(double energyold, double energynew, double temperature)
{

	double ran2(long *);

	/*DEBUG   printf ("   Move trial:    %13.8lf %13.8lf %13.8lf %13.8lf\n",
	  energynew, energyold, temperature, ran2(&seed));*/
	if (energynew <= energyold ) 
		return 0;
	else if (exp(-1.0*(energynew-energyold)/temperature) > ran2(&seed)) 
		return 0;
	else 
		return 1;
}
/*..............................................................................*/


/*
 * Calculate the different energy contributions. This is a merge of the different
 * energy calculation functions (energyone, -chain, -all)
 * 0: all
 * 1: one
 * 2: chain
 */
double calc_energy(long target, double (* intfce[MAXT][MAXT])(struct interacts *), 
		int mode, struct topo * topo, struct conf * conf, struct sim * sim, int chainnum)
{
	long i=0,j=0;
	
	double paire(long, long, double (* intfce[MAXT][MAXT])(struct interacts *), 
			struct topo * topo, struct conf * conf);
//	double extere(long, struct topo * topo, struct conf * conf);
	double extere2(long, struct topo * topo, struct conf * conf);

	//DEBUG_SIM("Calculate the energy with mode %d", mode)
	double energy = 0;
	/* Calculates energy between particle "target" and the rest.  Returns energy */
	if(mode == 1){
		if (sim->pairlist_update) {
#ifdef OMP
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
			for (i = 0; i < sim->pairlist[target].num_pairs; i++){
				energy+= paire(target, sim->pairlist[target].pairs[i], intfce, topo, conf);
			}
		}

		else{
#ifdef OMP
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
			for (i = 0; i < target; i++) {
				energy+= paire(target, i, intfce, topo, conf);
			}
#ifdef OMP
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
			for (i = target + 1; i < topo->npart; i++) {
				energy+= paire(target, i, intfce, topo, conf);
			}
		}
		/*add interaction with external potential*/
		if (topo->exter.exist)
			energy+= extere2(target,topo,conf);
		
	}
	/*
	 * Calculates energy between particle "target" and the rest. skipping
	 * particles from the given chain -particles has to be sorted in chain!!
	 * so similar to energy one but with chain exception
	 */
	else if(mode == 2){
		//#ifdef OMP
		//#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
		//#endif
		for (i = 0; i < target; i++) {
			if (i != topo->chainlist[chainnum][j]) {
				energy+= paire(target, i, intfce, topo, conf);
			}
			else {
				j++;
			}
		}
		j++;
		//#ifdef OMP
		//#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
		//#endif
		for (i = target + 1; i < topo->npart; i++) {
			if (i != topo->chainlist[chainnum][j]) {
				energy+= paire(target, i, intfce, topo, conf);
			}
			else {
				j++;
			}
		}
		/*add interaction with external potential*/
		if (topo->exter.exist)
			energy+= extere2(target,topo,conf);
	}
	
	/* Calculates energy between all pairs. Returns energy */
	else if(mode == 0){
#ifdef OMP
#pragma omp parallel for private(i,j) reduction (+:energy) schedule (dynamic)
#endif
		for (i = 0; i < topo->npart - 1; i++) {
			for (j = i + 1; j < topo->npart; j++) {
				energy+= paire(i, j, intfce, topo, conf);
			}
			/*for every particle add interaction with external potential*/
			if (topo->exter.exist)
				energy+= extere2(i,topo,conf);
		}
		/*add interaction of last particle with external potential*/
		if (topo->exter.exist)
			energy+= extere2(topo->npart-1,topo,conf);
	}
	else {
		fprintf(stderr, "ERROR: Wrong mode (%d) was given to calc_energy!", mode);
		return 0.0;
	}
	//  DEBUG_SIM("Will return energy from calc_energy")
	//printf("energymove %f\n",energy);
	return energy;
}
/*..............................................................................*/

/*
   Checks for overlaps between particle "target" and the rest.  Returns 1 if overlap
   detected, 0 otherwise.
 */
int forbidden(long npart, struct particles *particle,
		long target, struct vector box, struct ia_param ia_params[MAXT][MAXT])
{
	long test;
	int overlap(struct particles, struct particles, struct vector,struct ia_param [MAXT][MAXT]);

	for (test=0; test<npart; test++) {
		if (test != target) {
			if ( overlap(particle[target], particle[test], box, ia_params) ) {
				return 1;
			}
		}
	}

	return 0;
}

/*..............................................................................*/

/*
   Checks for overlaps between all pairs of particles. Returns 1 if overlap
   detected, 0 otherwise.
 */

int checkall(long npart, struct particles *particle,
		struct vector box, struct ia_param ia_params[MAXT][MAXT])
{
	long i, j;
	int overlap(struct particles, struct particles, struct vector, 
			struct ia_param [MAXT][MAXT]);

	for (i=0; i<npart-1; i++) {
		for (j=i+1; j<npart; j++) {
			if ( overlap(particle[i], particle[j], box, ia_params) ) {
				return 1;
			}
		}
	}

	return 0;
}


/*..............................................................................*/

/*
   Optimize the maximum displacement within the specified limits and resets the
   acceptance counters to zero.
 */

void optimizestep(struct disp *x, double hi, double lo)
{
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

/*..............................................................................*/

/*
   Optimize the maximum rotation within the specified limits and resets the
   acceptance counters to zero. Rotation is given by cos of angle
   larger rotation = smaller cos
 */

void optimizerot(struct disp *x, double hi, double lo)
{
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


/*................................................................................*/

/*
   Accumulate a value into the statistics and update the mean and rms values.
 */
void accumulate(struct stat *q, double x)
{
	(*q).sum += x;
	(*q).sum2 += x*x;
	(*q).samples++;
	(*q).mean = (*q).sum / (*q).samples;
	(*q).rms = sqrt(fabs((*q).sum2 / (*q).samples -
				(*q).sum * (*q).sum / (*q).samples / (*q).samples));
}

void printeqstat(struct disp *dat, double scale, int length)
{
	int i;

	for (i=0;i<length;i++) {
		if (RATIO(dat[i]) > 0)
			printf ("   TYPE %d           %.6lf  /  %.6lf\n", i, dat[i].mx/scale,RATIO(dat[i]));
	}
}

int memoryalloc(struct conf * conf)
{

	printf ("Allocating memory...\n");
	conf->particle = malloc( sizeof(struct particles)*MAXN);
	if(conf->particle == NULL){
		return 1;
	}

	return 0;
}
int memorydealloc(struct conf * conf, struct topo * topo, struct sim * sim)
{
	int dealloc_pairlist(struct topo * topo, struct sim * sim);
	printf ("Deallocating memory...\n");
	if (conf->particle != NULL) 
	    free(conf->particle);
	conf->particle = NULL;
	if (sim->clusterlist != NULL) 
	    free(sim->clusterlist);
	if (sim->clustersenergy != NULL) 
	    free(sim->clustersenergy); 
	if(topo->switchlist){
	    free(topo->switchlist);
	}
	if (sim->pairlist_update) {
		if(dealloc_pairlist(topo, sim)){
			return 1;
		}
	}
	return 0;
}
/*............................................................................*/


/**
 * nice malloc, which does the error checking for us 
 */
void * xmalloc (size_t num){
	void *new = malloc (num);
	if (!new){
		fprintf(stderr, "Couldn't allocate any memory!\n");
		exit(1);
	}
	return new;
}
/*............................................................................*/


/* ***********************   GEOMETRICAL FUNCTIONS  ****************************  */
/*.........................PATCHY SPOHEROCYLINDERS INTERACTION....................*/
/*................................................................................*/

/*
   Calculate intersections of sc2 with a patch of sc1 and return them in
 */
int psc_intersect(struct particles * part1, struct particles * part2,
		double halfl1, double halfl2, struct vector r_cm, double intersections[5], double rcut,
		struct ia_param * param, int which, int patchnum)

{
	int intrs;
	double a, b, c, d, e, x1, x2, rcut2;
	struct vector cm21, vec1, vec2, vec3, vec4;

	struct vector vec_crossproduct(struct vector, struct vector);
	struct vector vec_sub(struct vector, struct vector);
	struct vector vec_create(double, double, double);
	struct vector vec_scale(struct vector, double);
	struct vector vec_perpproject(struct vector*, struct vector*);
	struct quat quat_create(struct vector, double, double);
	void vec_rotate(struct vector *, struct quat);
	void normalise(struct vector *);
	int find_intersect_plane(struct particles *, struct particles *, double,
			struct vector, struct vector, double, double, double *);
	int test_intrpatch(struct particles *, struct vector, double, double, double *,int);


	intrs=0;
	rcut2=rcut*rcut;
	/*1- do intersections of spherocylinder2 with patch of spherocylinder1 at 
	  cut distance C*/
	/*1a- test intersection with half planes of patch and look how far they are 
	  from spherocylinder. If closer then C  we got itersection*/

	/* plane1 */
	/* find intersections of part2 with plane by par1 and patchsides[0] */
	intrs+=find_intersect_plane(part1,part2,halfl2,r_cm,part1->patchsides[0+2*patchnum],rcut,param->pcanglsw[which+2*patchnum],intersections);
	//	    printf("plane1 %d\n", intrs);
	/* plane2 */
	/* find intersections of part2 with plane by par1 and patchsides[1] */
	intrs+=find_intersect_plane(part1,part2,halfl2,r_cm,part1->patchsides[1+2*patchnum],rcut,param->pcanglsw[which+2*patchnum],intersections);

	if ( (intrs == 2 ) && (param->pcanglsw[which+2*patchnum] <0) ) {
		fprintf (stderr, "ERROR: Patch is larger than 180 degrees and we are getting two segments - this hasnot been programed yet.\n\n");
		exit (1);
	}
	//	    printf("plane2 %d\n", intrs);

	/*1b- test intersection with cylinder - it is at distance C*/
	if (intrs < 2 )  {
		cm21=vec_scale(r_cm,-1.0);
		vec1=vec_crossproduct(cm21,part1->dir);
		vec2=vec_crossproduct(part2->dir,part1->dir);
		a = DOT(vec2,vec2);
		b = 2*DOT(vec1,vec2);
		c = -rcut*rcut + DOT(vec1,vec1);
		d = b*b - 4*a*c;
		if ( d >= 0) { /*there is intersection with infinite cylinder */
			x1 = (-b+sqrt(d))*0.5/a;/*parameter on line of SC2 determining intersection*/
			if ((x1 >=halfl2) || (x1 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
			else {
				/* vectors from center os sc1 to intersection with infinite cylinder*/
				vec1.x=part2->dir.x*x1-r_cm.x;
				vec1.y=part2->dir.y*x1-r_cm.y;
				vec1.z=part2->dir.z*x1-r_cm.z;
				e = DOT(part1->dir,vec1);
				if ((e >=halfl1) || (e <= -halfl1)) intrs+=0; /*intersection is outside sc1*/
				else {
					intrs+=test_intrpatch(part1,vec1,param->pcanglsw[which+2*patchnum],x1,intersections,patchnum);
				}
			}
			if ( d > 0 ){
				x2 = (-b-sqrt(d))*0.5/a;/*parameter on line of SC2 determining intersection*/
				if ((x2 >=halfl2) || (x2 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
				else {
					vec2.x = part2->dir.x*x2-r_cm.x;
					vec2.y = part2->dir.y*x2-r_cm.y;
					vec2.z = part2->dir.z*x2-r_cm.z;
					e = DOT(part1->dir,vec2);
					if ((e >=halfl1) || (e <= -halfl1)) intrs+=0; /*intersection is outside sc1*/
					else {
						intrs+=test_intrpatch(part1,vec2,param->pcanglsw[which+2*patchnum],x2,intersections,patchnum);
					}
				}
			}
		}
	}
	//	    printf ("cylinder %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
	/*1c- test intersection with spheres at the end - it is at distace C*/
	if (intrs < 2 )  {
		/*centers of spheres*/
		/*relative to the CM of sc2*/
		vec1.x =  part1->dir.x*halfl1 - r_cm.x;
		vec1.y =  part1->dir.y*halfl1 - r_cm.y;
		vec1.z =  part1->dir.z*halfl1 - r_cm.z;
		vec2.x = -part1->dir.x*halfl1 - r_cm.x;
		vec2.y = -part1->dir.y*halfl1 - r_cm.y;
		vec2.z = -part1->dir.z*halfl1 - r_cm.z;

		/*sphere1*/
		a = DOT(part2->dir,part2->dir);
		b = 2.0*DOT(vec1,part2->dir);
		c = DOT(vec1,vec1)-rcut*rcut;
		d = b*b-4*a*c;
		if (d >= 0) { /*if d<0 there are no intersections*/
			x1= (-b + sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
			if ((x1 >=halfl2) || (x1 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
			else {
				vec3.x = part2->dir.x*x1-r_cm.x;
				vec3.y = part2->dir.y*x1-r_cm.y;
				vec3.z = part2->dir.z*x1-r_cm.z;
				e = DOT(part1->dir,vec3);
				if ((e >= halfl1) || (e <= -halfl1)) { /*if not intersection is inside sc1*/
					intrs+=test_intrpatch(part1,vec3,param->pcanglsw[which+2*patchnum],x1,intersections,patchnum);
				}
			}
			if ( d > 0) {
				x2= (-b - sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
				if ((x2 >=halfl2) || (x2 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
				else {
					vec4.x = part2->dir.x*x2 - r_cm.x;
					vec4.y = part2->dir.y*x2 - r_cm.y;
					vec4.z = part2->dir.z*x2 - r_cm.z;
					e = DOT(part1->dir,vec4);
					if ((e >=halfl1) || (e <= -halfl1)) { /*if not intersection is inside sc1*/
						intrs+=test_intrpatch(part1,vec4,param->pcanglsw[which+2*patchnum],x2,intersections,patchnum);
					}
				}
			}
		}
		//		printf ("sphere1 %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
		/*sphere2*/
		a = DOT(part2->dir,part2->dir);
		b = 2.0*DOT(vec2,part2->dir);
		c = DOT(vec2,vec2)-rcut*rcut;
		d = b*b-4*a*c;
		if (d >= 0) { /*if d<0 there are no intersections*/
			x1= (-b + sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
			if ((x1 >=halfl2) || (x1 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
			else {
				vec3.x = part2->dir.x*x1 - r_cm.x;
				vec3.y = part2->dir.y*x1 - r_cm.y;
				vec3.z = part2->dir.z*x1 - r_cm.z;
				e = DOT(part1->dir,vec3);
				if ((e >=halfl1) || (e <= -halfl1)) { /*if not intersection is inside sc1*/
					intrs+=test_intrpatch(part1,vec3,param->pcanglsw[which+2*patchnum],x1,intersections,patchnum);
				}
			}
			if ( d > 0 ) {
				x2= (-b - sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
				if ((x2 >=halfl2) || (x2 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
				else {
					vec4.x = part2->dir.x*x2 - r_cm.x;
					vec4.y = part2->dir.y*x2 - r_cm.y;
					vec4.z = part2->dir.z*x2 - r_cm.z;
					e = DOT(part1->dir,vec4);
					if ((e >=halfl1) || (e <= -halfl1)) { /*if not intersection is inside sc1*/
						intrs+=test_intrpatch(part1,vec4,param->pcanglsw[which+2*patchnum],x2,intersections,patchnum);
					}
				}
			}
		}
		//		printf ("sphere2 %d\n", intrs);
	}

	/*1d- if there is only one itersection shperocylinder ends within patch wedge 
	  set as second intersection end inside patch*/
	if (intrs < 2 )  {
		/*whole spherocylinder is in or all out if intrs ==0*/
		vec1.x = part2->dir.x*halfl2 - r_cm.x;
		vec1.y = part2->dir.y*halfl2 - r_cm.y;
		vec1.z = part2->dir.z*halfl2 - r_cm.z;
		/*vector from CM of sc1 to end of sc2*/
		/*check is is inside sc1*/
		a=DOT(vec1,part1->dir);
		vec3.x = vec1.x - part1->dir.x*a;
		vec3.y = vec1.y - part1->dir.y*a;
		vec3.z = vec1.z - part1->dir.z*a;
		b=DOT(vec3,vec3);
		d = fabs(a)-halfl1;
		if ( d <= 0) 
			c = b; /*is inside cylindrical part*/
		else 
			c = d*d + b; /*is inside caps*/
		/*c is distance squared from line or end to test if is inside sc*/
		if (c < rcut2) 
			intrs+=test_intrpatch(part1,vec1,param->pcanglsw[which+2*patchnum],halfl2,intersections,patchnum);
		if (intrs < 2 ) {
			vec2.x = -part2->dir.x*halfl2 - r_cm.x;
			vec2.y = -part2->dir.y*halfl2 - r_cm.y;
			vec2.z = -part2->dir.z*halfl2 - r_cm.z;
			/*check is is inside sc1*/
			a=DOT(vec2,part1->dir);
			vec4.x = vec2.x - part1->dir.x*a;
			vec4.y = vec2.y - part1->dir.y*a;
			vec4.z = vec2.z - part1->dir.z*a;
			b=DOT(vec4,vec4);
			d = fabs(a) -halfl1;
			if (d <= 0) 
				c = b; /*is inside cylindrical part*/
			else 
				c = d*d + b; /*is inside caps*/
			/*c is distance squared from line or end to test if is inside sc*/
			if (c < rcut2) 
				intrs+=test_intrpatch(part1,vec2,param->pcanglsw[which+2*patchnum],-1.0*halfl2,intersections,patchnum);
		}
		//		    printf ("ends %d\n", intrs);
	}


	return intrs;
}

/*................................................................................*/
/*
   Find if vector vec has angular intersection with patch of sc1 
 */
int test_intrpatch(struct particles * part1, struct vector vec, double cospatch, 
	double ti, double intersections[5],int patchnum)
{
	double a;
	int i, intrs;

	struct vector vec_perpproject(struct vector*, struct vector*);
	void normalise(struct vector *);

	intrs=0;
	/*test if we have intersection*/
	/* do projection to patch plane*/
	vec=vec_perpproject(&vec,&part1->dir);
	normalise(&vec);
	/* test angle distance from patch*/
	a = DOT(part1->patchdir[patchnum],vec);
	if (a >= cospatch) {
		intrs=1;
		i=0;
		while (intersections[i] !=0) {
			if (ti == intersections[i]) 
				intrs=0; /* found intersection we already have -it is at boundary*/
			i++;
		}
		if (intrs > 0) 
			intersections[i]=ti;
	}

	return intrs;
}

/*................................................................................*/

/*
   Find intersections of SC and plane defined by vector w_vec.and returns number of them
 */
int find_intersect_plane(struct particles * part1, struct particles * part2, double halfl2, 
		struct vector r_cm, struct vector w_vec, double rcut, double cospatch, double intersections[5])
{
	int i, intrs;
	double a, c, d, ti, disti;
	struct vector nplane, d_vec;

	void normalise(struct vector *);
	struct vector vec_crossproduct(struct vector, struct vector);

	nplane=vec_crossproduct(part1->dir,w_vec);
	normalise(&nplane);
	normalise(&w_vec);
	a =  DOT(nplane, part2->dir);
	if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
	else {
		ti = DOT(nplane,r_cm)/a;
		if ((ti  > halfl2 ) || (ti < -halfl2)) intrs=0; /* there is no intersection plane sc is too short*/
		else {
			d_vec.x = ti * part2->dir.x - r_cm.x; /*vector from intersection point to CM*/
			d_vec.y = ti * part2->dir.y - r_cm.y;
			d_vec.z = ti * part2->dir.z - r_cm.z;
			c = DOT (d_vec, w_vec);
			if ( c * cospatch < 0) intrs=0; /* the intersection in plane is on other side of patch */
			else {
				d = fabs(DOT (d_vec, part1->dir)) - halfl2;
				if (d <= 0) disti = c*c; /*is inside cylinder*/
				else disti = d*d + c*c; /*is inside patch*/
				if (disti > rcut*rcut) intrs=0; /* the intersection is outside sc */
				else {
					intrs=1;
					i=0;
					while (intersections[i] !=0) {
						if (ti == intersections[i]) intrs=0; /* found intersection we already have -it is at boundary*/
						i++;
					}
					if (intrs > 0) {
					    intersections[i]=ti;
					}
				}
			}
		}

	}
	return intrs;
}


/*CPSC................................................................................*/

/*
   Calculate intersections of sc2 with a patch of sc1 and return them in this works for cylindrical psc -CPSC
 */
int cpsc_intersect(struct particles * part1, struct particles * part2,
		double halfl1, double halfl2, struct vector r_cm, double intersections[5], double rcut,
		struct ia_param * param, int which, int patchnum)

{
	int intrs;
	double a, b, c, d, e, x1, x2, rcut2;
	struct vector cm21, vec1, vec2, vec3, vec4;

	struct vector vec_crossproduct(struct vector, struct vector);
	struct vector vec_sub(struct vector, struct vector);
	struct vector vec_create(double, double, double);
	struct vector vec_scale(struct vector, double);
	struct vector vec_perpproject(struct vector*, struct vector*);
	struct quat quat_create(struct vector, double, double);
	void vec_rotate(struct vector *, struct quat);
	void normalise(struct vector *);
	int find_intersect_planec(struct particles *, struct particles *, double,
			struct vector, struct vector, double, double, double *);
	int test_intrpatch(struct particles *, struct vector, double, double, double *, int);


	intrs=0;
	rcut2=rcut*rcut;
	/*1- do intersections of spherocylinder2 with patch of spherocylinder1 at 
	  cut distance C*/
	/*1a- test intersection with half planes of patch and look how far they are 
	  from spherocylinder. If closer then C  we got itersection*/

	/* plane1 */
	/* find intersections of part2 with plane by par1 and part1->patchsides[0] */
	intrs+=find_intersect_planec(part1,part2,halfl2,r_cm,part1->patchsides[0+2*patchnum],rcut,param->pcanglsw[which+2*patchnum],intersections);
	//	    printf("plane1 %d\n", intrs);
	/* plane2 */
	/* find intersections of part2 with plane by par1 and part1->patchsides[1] */
	intrs+=find_intersect_planec(part1,part2,halfl2,r_cm,part1->patchsides[1+2*patchnum],rcut,param->pcanglsw[which+2*patchnum],intersections);

	if ( (intrs == 2 ) && (param->pcanglsw[which+2*patchnum] < 0) ) {
		fprintf (stderr, "ERROR: Patch is larger than 180 degrees and we are getting two segments - this hasnot been programed yet.\n\n");
		exit (1);
	}
	//	    printf("plane2 %d\n", intrs);

	/*1b- test intersection with cylinder - it is at distance C*/
	if (intrs < 2 )  {
		cm21=vec_scale(r_cm,-1.0);
		vec1=vec_crossproduct(cm21,part1->dir);
		vec2=vec_crossproduct(part2->dir,part1->dir);
		a = DOT(vec2,vec2);
		b = 2*DOT(vec1,vec2);
		c = -rcut*rcut + DOT(vec1,vec1);
		d = b*b - 4*a*c;
		if ( d >= 0) { /*there is intersection with infinite cylinder */
			x1 = (-b+sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
			if ((x1 >=halfl2) || (x1 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
			else {
				/* vectors from center os sc1 to intersection with infinite cylinder*/
				vec1.x=part2->dir.x*x1-r_cm.x;
				vec1.y=part2->dir.y*x1-r_cm.y;
				vec1.z=part2->dir.z*x1-r_cm.z;
				e = DOT(part1->dir,vec1);
				if ((e >=halfl1) || (e <= -halfl1)) intrs+=0; /*intersection is outside sc1*/
				else {
					intrs+=test_intrpatch(part1,vec1,param->pcanglsw[which+2*patchnum],x1,intersections,patchnum);
				}
			}
			if ( d > 0 ){
				x2 = (-b-sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
				if ((x2 >=halfl2) || (x2 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
				else {
					vec2.x = part2->dir.x*x2-r_cm.x;
					vec2.y = part2->dir.y*x2-r_cm.y;
					vec2.z = part2->dir.z*x2-r_cm.z;
					e = DOT(part1->dir,vec2);
					if ((e >=halfl1) || (e <= -halfl1)) intrs+=0; /*intersection is outside sc1*/
					else {
						intrs+=test_intrpatch(part1,vec2,param->pcanglsw[which+2*patchnum],x2,intersections,patchnum);
					}
				}
			}
		}
	}
	//	    printf ("cylinder %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
	/*1c- test intersection with plates at the end - it is at distace C and in wedge*/
	if (intrs < 2 )  {
		a =  DOT(part1->dir, part2->dir);
		if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
		else {
			/*plane cap1*/
			vec1.x= r_cm.x + halfl1*part1->dir.x;
			vec1.y= r_cm.y + halfl1*part1->dir.y;
			vec1.z= r_cm.z + halfl1*part1->dir.z;
			x1 = DOT(part1->dir,vec1)/a; /*parameter on line of SC2 determining intersection*/
			if ((x1 > halfl2 ) || (x1 < -halfl2)) intrs+=0; /* there is no intersection plane sc is too short*/
			else {
				vec2.x = x1*part2->dir.x - vec1.x; /*vector from ENDPOINT to intersection point */
				vec2.y = x1*part2->dir.y - vec1.y;
				vec2.z = x1*part2->dir.z - vec1.z;
				b = DOT (vec2, vec2);
				if (b > rcut*rcut) intrs+=0;   /* the intersection is outside sc */
				else {
					intrs+=test_intrpatch(part1,vec2,param->pcanglsw[which+2*patchnum],x1,intersections,patchnum);
				}
			}
			//		    printf ("plane cap1 %d %f\n", intrs, x1);
			/*plane cap2*/
			vec1.x= r_cm.x - halfl1*part1->dir.x;
			vec1.y= r_cm.y - halfl1*part1->dir.y;
			vec1.z= r_cm.z - halfl1*part1->dir.z;
			x2 = DOT(part1->dir,vec1)/a; /*parameter on line of SC2 determining intersection*/
			if ((x2  > halfl2 ) || (x2 < -halfl2)) intrs+=0; /* there is no intersection plane sc is too short*/
			else {
				vec2.x = x2*part2->dir.x - vec1.x; /*vector from ENDPOINT to intersection point */
				vec2.y = x2*part2->dir.y - vec1.y;
				vec2.z = x2*part2->dir.z - vec1.z;
				b = DOT (vec2, vec2);
				if (b > rcut*rcut) intrs+=0;   /* the intersection is outside sc */
				else {
					intrs+=test_intrpatch(part1,vec2,param->pcanglsw[which+2*patchnum],x2,intersections,patchnum);
				}
			}
			//		    printf ("plane cap2 %d %f\n", intrs,x2);

		}
	}

	/*1d- if there is only one itersection shperocylinder ends within patch wedge 
	  set as second intersection end inside patch*/
	if (intrs < 2 )  {
		/*whole spherocylinder is in or all out if intrs ==0*/
		vec1.x = part2->dir.x*halfl2 - r_cm.x;
		vec1.y = part2->dir.y*halfl2 - r_cm.y;
		vec1.z = part2->dir.z*halfl2 - r_cm.z;
		/*vector from CM of sc1 to end of sc2*/
		/*check is is inside sc1*/
		a=DOT(vec1,part1->dir);
		vec3.x = vec1.x - part1->dir.x*a;
		vec3.y = vec1.y - part1->dir.y*a;
		vec3.z = vec1.z - part1->dir.z*a;
		b=DOT(vec3,vec3);
		d = fabs(a)-halfl1;
		if ( d <= 0) { /*is in cylindrical part*/
			/*c is distance squared from line or end to test if is inside sc*/
			if (b < rcut2) intrs+=test_intrpatch(part1,vec1,param->pcanglsw[which+2*patchnum],halfl2,intersections,patchnum);
		}
		if (intrs < 2 ) {
			vec2.x = -part2->dir.x*halfl2 - r_cm.x;
			vec2.y = -part2->dir.y*halfl2 - r_cm.y;
			vec2.z = -part2->dir.z*halfl2 - r_cm.z;
			/*check is is inside sc1*/
			a=DOT(vec2,part1->dir);
			vec4.x = vec2.x - part1->dir.x*a;
			vec4.y = vec2.y - part1->dir.y*a;
			vec4.z = vec2.z - part1->dir.z*a;
			b=DOT(vec4,vec4);
			d = fabs(a) -halfl1;
			if (d <= 0) {
				/*c is distance squared from line or end to test if is inside sc*/
				if (b < rcut2) intrs+=test_intrpatch(part1,vec2,param->pcanglsw[which+2*patchnum],-1.0*halfl2,intersections,patchnum);
			}
		}
				//    printf ("ends %d\n", intrs);
	}


	return intrs;
}

/*CPSC................................................................................*/

/*
   Find intersections of plane defined by vector w_vec.and returns number of them - for cylindrical psc -CPSC
 */
int find_intersect_planec(struct particles * part1, struct particles * part2, double halfl, 
		struct vector r_cm, struct vector w_vec, double rcut, double cospatch, double intersections[5])
{
	int i, intrs=0;
	double a, c, d, ti, disti;
	struct vector nplane, d_vec;

	void normalise(struct vector *);
	struct vector vec_crossproduct(struct vector, struct vector);

	nplane=vec_crossproduct(part1->dir,w_vec);
	normalise(&nplane);
	normalise(&w_vec);
	a =  DOT(nplane, part2->dir);
	if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
	else {
		ti = DOT(nplane,r_cm)/a;
		if ((ti  > halfl ) || (ti < -halfl)) intrs=0; /* there is no intersection plane sc is too short*/
		else {
			d_vec.x = ti*part2->dir.x - r_cm.x; /*vector from intersection point to CM*/
			d_vec.y = ti*part2->dir.y - r_cm.y;
			d_vec.z = ti*part2->dir.z - r_cm.z;
			c = DOT (d_vec, w_vec);
			if ( c *cospatch < 0) intrs=0; /* the intersection in plane is on other side of patch */
			else {
				d = fabs(DOT (d_vec, part1->dir)) - halfl;
				if (d <= 0) {
					disti= c*c; /*is inside cylinder*/
					if (disti > rcut*rcut) intrs=0; /* the intersection is outside sc */
					else {
						intrs=1;
						i=0;
						while (intersections[i] !=0) {
							if (ti == intersections[i]) intrs=0; /* found intersection we already have -it is at boundary*/
							i++;
						}
						if (intrs > 0) intersections[i]=ti;
					}
				}
			}
		}
	}
	return intrs;
}


/*..................................................................................*/

/*
   Find projection of cpsc on plane (0,0,1) including cutoff and return 
   vector to its begining and end and cm
 */

int psc_wall(struct vector* pbeg, struct vector* pend, struct vector* projectdir, 
		struct vector* partdir, BOOL* positive, double * cutdist, 
		struct vector *partbeg, struct vector *partend)
{
	struct vector vec1;
	double k,x1,x2,y1,y2,a,b,c,e,d;

	void projectinz( struct vector* vec1, struct vector * projectdir, struct vector* result);
	void normalise(struct vector*);
	
	if (( (*positive)&& (projectdir->z > 0) ) || ( (!(*positive))&& (projectdir->z < 0) ))
		return 0;
	if ( fabs(partbeg->z) > (*cutdist) )
		return 0;
	/* we might have interacting segment*/
	x2 = 0.0;
	y2 = 0.0;
	/*begining point*/ 
	/*if begining projected along particle direction is within cutoff */
	if (fabs(partdir->z) > ZEROTOL2) {
		projectinz(partbeg,partdir,pbeg);
		a=0;
	}
	else {
		/*we need some starting point*/
		vec1.x = 2.0*partbeg->x - partend->x;
		vec1.y = 2.0*partbeg->y - partend->y;
		vec1.z = 2.0*partbeg->z - partend->z;
		projectinz(&vec1,projectdir,pbeg);
		a=1;
	}
	if (partdir->z != 0) {
		b = fabs(partbeg->z / partdir->z);
	} else {
		b = (*cutdist)+1.0;
	}


	if ( (b > (*cutdist)) || (a==1)) {
		/*else beginig is at sphere, find intersections with sphere of cutoff radius*/
		if ( fabs(projectdir->z) > ZEROTOL2) {
			projectinz(partbeg,projectdir,pend);
		} else {
			pend->x = pbeg->x + projectdir->x;
			pend->y = pbeg->y + projectdir->y;
		}
		if (pend->y == pbeg->y) {
			y1=pbeg->y;
			y2=pbeg->y;
			a=sqrt( (*cutdist)*(*cutdist) - partbeg->z*partbeg->z - (pbeg->y-partbeg->y)*(pbeg->y-partbeg->y) );
			x1 = partbeg->x + a;
			x2 = partbeg->x - a;
			if (pend->x > pbeg->x) {/*select the right intersection*/
				pbeg->x = x2;	
				x2 = x1;
			} else {
				pbeg->x = x1;
			}
			pbeg->y = y1;
		} else {
			k = (pend->x - pbeg->x)/ (pend->y - pbeg->y);
			a = k*k +1;
			b = partbeg->y + k*k*pbeg->y - k*pbeg->x + k*partbeg->x;
			c = partbeg->y*partbeg->y + partbeg->z*partbeg->z - (*cutdist)*(*cutdist) + (k*pbeg->y - pbeg->x + partbeg->x)*(k*pbeg->y - pbeg->x + partbeg->x);
			e = b*b-a*c;
			if (e < 0) {
				return 0; /*tehre might be no intersection with sphere*/
			}
			d = sqrt(e);
			if (pend->y > pbeg->y) {/*select the right intersection*/
				y1 = (b - d ) /a;
				y2 = (b + d ) /a;
			}
			else {
				y1 = (b + d ) /a;
				y2 = (b - d ) /a;
			}
			x1 = k * (y1 - pbeg->y) + pbeg->x;
			x2 = k * (y2 - pbeg->y) + pbeg->x;
			pbeg->x = x1;
			pbeg->y = y1;
			pbeg->z = 0.0;
		}
	}
	//printf("pscwall beg %f %f \n",pbeg->x,pbeg->y);
	
	/*end point*/
	a = -(*cutdist) * projectdir->z; /*z coordinate of point where projection is in cut distance*/
	//printf("sphere end %f %f   ",a,partend->z);
	if ( ((partend->z < a)&&(*positive)) || ((a < partend->z)&&(!(*positive))) ){ 
		/*end is within cut off - second sphere*/
		/*if this is the case vec1 is end of pherocylinder and pend is its projection*/
		if (projectdir->z != 0) {
			projectinz(partend,projectdir,pend);
		} else {
			pend->x = pbeg->x + projectdir->x;
			pend->y = pbeg->y + projectdir->y;
		}
		if (pend->y == pbeg->y) {
			y1=pend->y;
			y2=pend->y;
			a=sqrt( (*cutdist)*(*cutdist) - partend->z*partend->z - (pend->y-partend->y)*(pend->y-partend->y) );
			x1 = partend->x + a;
			x2 = partend->x - a;
			if (pbeg->x > pend->x) {/*select the right intersection*/
				pend->x = x2;	
			} else {
				pend->x = x1;
			}
			pend->y = y1;
		} else {
			k = (pbeg->x - pend->x)/ (pbeg->y - pend->y);
			a = k*k +1;
			b = partend->y + k*k*pend->y - k*pend->x + k*partend->x;
			c = partend->y*partend->y + partend->z*partend->z - (*cutdist)*(*cutdist) + (k*pend->y - pend->x + partend->x)*(k*pend->y - pend->x + partend->x);
			e = b*b-a*c;
			if (e < 0) {
				return 0; /*there might be no intersection with sphere*/
			}
			d = sqrt(e);
			if (pbeg->y > pend->y) {/*select the right intersection*/
				y1 = (b - d ) /a;
				y2 = (b + d ) /a;
			}
			else {
				y1 = (b + d ) /a;
				y2 = (b - d ) /a;
			}
			x1 = k * (y1 - pend->y) + pend->x;
			x2 = k * (y2 - pend->y) + pend->x;
			pend->x = x1;
			pend->y = y1;
			pend->z = 0.0;
		}
	} else {
		if ( ((partbeg->z < a)&&(*positive)) || ((a < partbeg->z)&&(!(*positive))) ) {
			/*end is at cutoff going through cylindrical part*/
			//printf("cylinder ");
			b = (a - partbeg->z)/ partdir->z;
			vec1.x = partbeg->x + b * partdir->x; 
			vec1.y = partbeg->y + b * partdir->y;
			vec1.z = a;
			projectinz(&vec1,projectdir,pend);
		} else {
			/* also projected end is within the same sphere as begining- no contribution from cylinder*/
			if (x2 == 0.0 ) {
				//printf("sphere beg ");
				if (projectdir->z != 0) {
					projectinz(partbeg,projectdir,pend);
				} else {
					pend->x = pbeg->x + projectdir->x;
					pend->y = pbeg->y + projectdir->y;
				}
				
				if (pend->y == pbeg->y) {
					y1=pbeg->y;
					y2=pbeg->y;
					a=sqrt( (*cutdist)*(*cutdist) - partbeg->z*partbeg->z - (pbeg->y-partbeg->y)*(pbeg->y-partbeg->y) );
					x1 = partbeg->x + a;
					x2 = partbeg->x - a;
					if (pend->x > pbeg->x) {/*select the right intersection*/
						pend->x = x1;	
					} else {
						pend->x = x2;
					}
					pend->y = y1;
				} else {
					k = (pend->x - pbeg->x)/ (pend->y - pbeg->y);
					a = k*k +1;
					b = partbeg->y + k*k*pbeg->y - k*pbeg->x + k*partbeg->x;
					c = partbeg->y*partbeg->y + partbeg->z*partbeg->z - (*cutdist)*(*cutdist) + (k*pbeg->y - pbeg->x + partbeg->x)*(k*pbeg->y - pbeg->x + partbeg->x);
					e = b*b-a*c;
					if (e < 0) {
						return 0; /*tehre might be no intersection with sphere*/
					}
					d = sqrt(e);
					if (pend->y > pbeg->y) {/*select the right intersection*/
						y1 = (b - d ) /a;
						y2 = (b + d ) /a;
					}
					else {
						y1 = (b + d ) /a;
						y2 = (b - d ) /a;
					}
					x1 = k * (y1 - pbeg->y) + pbeg->x;
					x2 = k * (y2 - pbeg->y) + pbeg->x;
					pend->x = x1;
					pend->y = y1;
					pend->z = 0.0;
				}
			} else {
				pend->x = x2;
				pend->y = y2;
				pend->z = 0.0;
			}
			return 2; /*line end is on sphere of particle begining = no cylindrical cutoff*/
		}
	}

		
	return 1;
}

int cpsc_wall(struct vector* pbeg, struct vector* pend, struct vector* projectdir, 
		struct vector* partdir, double* halfl, BOOL* orientin,
		BOOL* positive, double* rcmz, double * cutdist, 
		struct vector *partbeg, struct vector *partend)
{
	struct vector vec1;
	double a;

	void projectinz( struct vector* vec1, struct vector * projectdir, struct vector* result);

	if (( (*positive)&& (projectdir->z >= 0) ) || ( (!(*positive))&& (projectdir->z <= 0) ))
		return 0;
	/*if projected closer point beoynd cutoff no interaction*/
	/*project begining of spherocylinder*/
	vec1.x = partbeg->x;
	vec1.y = partbeg->y;		
	vec1.z = partbeg->z;
	
	if (-vec1.z/projectdir->z < (*cutdist) ) {
		projectinz(&vec1,projectdir,pbeg);	
	} else {
		return 0;
	}

	/* we have interacting segment*/
	if (-partend->z/projectdir->z < (*cutdist) ) {
		/*whole segment interacts*/
		vec1.z = partend->z;	
	} else {
		vec1.z = -(*cutdist)*projectdir->z;
	}
	if (partdir->z != 0.0)
	    a = (vec1.z - (*rcmz)) / partdir->z;
	else {
	    if (*orientin) 
		a = -(*halfl);
	    else
		a = (*halfl);
	}
	vec1.x = partdir->x * a;
	vec1.y = partdir->y * a;
	projectinz(&vec1,projectdir,pend);

	return 1;

}

int cutprojectatwall(struct vector* pextr1, struct vector* pextr2, struct vector* pextr3, struct vector* pextr4, 
					 struct vector* projectdir, struct vector* partdir, double * cutdist, 
					 struct vector *partbeg, struct vector *partend, struct vector *pend, 
					 double *cuttoproject, BOOL* orientin)
{
	double y1,y2,O2z,det,a,b,dirydirz,dir2x,dir2y,dir2z,dirzldiry;
	
	void projectinz( struct vector* vec1, struct vector * projectdir, struct vector* result);

	dirydirz = partdir->y * partdir->z;
	dir2x = partdir->x * partdir->x;
	dir2y = partdir->y * partdir->y;
	dir2z = partdir->z * partdir->z;
	a = 1/(dir2x+dir2y);

	
	if (partdir->x != 0) {
		O2z = partbeg->z * partbeg->z;
		b=dir2y*dir2z*O2z - (dir2x+dir2y) * (O2z*(dir2x+dir2z)- (*cutdist)*(*cutdist)*dir2x);
		if (b < 0 ) {
			/*no cutoff from cylindrical part*/
			return 0;
		}
		det = sqrt(b);
		y1 = partbeg->y + (dirydirz*partbeg->z + det )*a;
		y2 = partbeg->y + (dirydirz*partbeg->z - det )*a;
	
		if (( (partdir->x > 0)&&(!(*orientin)) ) || ( (partdir->x < 0)&&(*orientin) ))  {
			pextr1->y = y1;
			pextr2->y = y2;
		} else {
			pextr1->y = y2;
			pextr2->y = y1;
		}
		pextr1->x = partbeg->x  + (partbeg->z*partdir->z -  (pextr1->y - partbeg->y)*partdir->y) / partdir->x;
		pextr2->x = partbeg->x  + (partbeg->z*partdir->z -  (pextr2->y - partbeg->y)*partdir->y) / partdir->x;

		O2z = partend->z * partend->z;
		b= dir2y*dir2z*O2z - (dir2x+dir2y) * (O2z*(dir2x+dir2z)- (*cutdist)*(*cutdist)*dir2x);
		if (b >= 0) { /*we have intersections from end*/
			det = sqrt(b);
			y1 = partend->y + (dirydirz * partend->z + det )*a;
			y2 = partend->y + (dirydirz * partend->z - det )*a;
			//printf("det %f y1 %f y2 %f \n", det,y1,y2);	
			if (( (partdir->x > 0)&&(!(*orientin)) ) || ( (partdir->x < 0)&&(*orientin) ))  {
				pextr3->y = y1;
				pextr4->y = y2;
			} else {
				pextr3->y = y2;
				pextr4->y = y1;
			}
			pextr3->x = partend->x  + (partend->z*partdir->z -  (pextr3->y - partend->y)*partdir->y) / partdir->x;
			pextr4->x = partend->x  + (partend->z*partdir->z -  (pextr4->y - partend->y)*partdir->y) / partdir->x;
		} else { 
			/*no intersection at the end the cutoff intersects the plane 
			 in the perpendicular projection of line segemnt, so we have to use that point */
			if (partdir->z == 0) {
				fprintf (stderr, "\nERROR: Something went wrong in calculation of projection.\n\n");
				exit (1);
			} else {
				a = ((*cuttoproject) - partbeg->z)/ partdir->z;
				//if ( projectdir->y * partdir->x < 0  ) 
				pextr3->x = partbeg->x + a * partdir->x; 
				pextr3->y = partbeg->y + a * partdir->y;
				pextr3->z = (*cuttoproject);
				//printf("before proj %f %f dir %f %f %f ",pextr3->x,pextr3->y,projectdir->x,projectdir->y,projectdir->z);
				projectinz(pextr3,projectdir,pextr4);
				pextr3->x = pextr4->x;  
				pextr3->y = pextr4->y; 
				pextr3->z = 0.0;
				//printf("after proj %f %f \n",pextr3->x,pextr3->y);
				return 2;
			}
		}
	} else {
		if (partdir->y != 0) {
			dirzldiry = partdir->z/partdir->y;
			y1 = partbeg->y + partbeg->z * dirzldiry;
			det = sqrt( (*cutdist)*(*cutdist) - partbeg->z * partbeg->z * (1+dirzldiry*dirzldiry) );
			if (( (partdir->y > 0)&&(!(*orientin)) ) || ( (partdir->y < 0)&&(*orientin) ))  {
				pextr1->x = partbeg->x  + det;
				pextr2->x = partbeg->x  - det;
			} else {
				pextr1->x = partbeg->x  - det;
				pextr2->x = partbeg->x  + det;
			}
			pextr1->y = y1;
			pextr2->y = y1;
			
			y1 = partend->y + partend->z * dirzldiry;
			b = (*cutdist)*(*cutdist) - partend->z * partend->z * (1+dirzldiry*dirzldiry); 
			if (b >= 0) { /*we have intersections from end*/
				det = sqrt(b);
				if (( (partdir->y > 0)&&(!(*orientin)) ) || ( (partdir->y < 0)&&(*orientin) )) {
					pextr3->x = partend->x  + det;
					pextr4->x = partend->x  - det;
				} else {
					pextr3->x = partend->x  - det;
					pextr4->x = partend->x  + det;
				}
				pextr3->y = y1;
				pextr4->y = y1;
			} else {
				/*no intersection at the end the cutoff intersects the plane 
				 in the perpendicular projection of line segemnt, so we have to use that point */
				if (partdir->z == 0) {
					fprintf (stderr, "\nERROR: Something went wrong in calculation of projection.\n\n");
					exit (1);
				} else {
					a = ((*cutdist) - partbeg->z)/ partdir->z;
					y1 = a * partdir->y + partbeg->y;
					if ( projectdir->x * partdir->y > 0  ) {
						pextr3->x = a * partdir->x  + partbeg->x; 
						pextr3->y = y1;
						pextr4->x = pend->x;  
						pextr4->y = pend->y; 
					}else {
						pextr3->x = pend->x;  
						pextr3->y = pend->y; 
						pextr4->x = a * partdir->x  + partbeg->x; 
						pextr4->y = y1;
					}
				}
			}

			
		} else {
			return 0; /* if perpendicular to plane we don't have any intersections*/
		}

		
	}
	
	return 1;
}


/*project a point in project direction to z plane z=0*/
void projectinz(struct vector* vec1, struct vector* projectdir,struct vector * projection)	
{
	projection->x = vec1->x - vec1->z * projectdir->x/projectdir->z;
	projection->y = vec1->y - vec1->z * projectdir->y/projectdir->z;
	projection->z = 0;

}

/*calculates area defined by four points in z=0 plane */
double areafourpoints(struct vector * pbeg, struct vector * pend, struct vector * pbeg1, struct vector * pend1 )
{
	double area =0.0;
	struct vector vec1,vec2;

	/*area by four points... two half vector cross product  
	  |(pbegining1-pbegining)x(pend-pbegining)|/2  */
	vec1.x = pbeg1->x - pbeg->x; 
	vec1.y = pbeg1->y - pbeg->y;  
	vec2.x = pend->x - pbeg->x;
	vec2.y = pend->y - pbeg->y;
	//printf("a: %f %f %f %f \n",vec1.x,vec2.y,vec1.y,vec2.x);
	area += fabs(vec1.x*vec2.y - vec1.y*vec2.x)*0.5;
	/* + |(pend-pend1)x(pbegining1-pend1)|/2*/
	vec1.x = pend->x - pend1->x; 
	vec1.y = pend->y - pend1->y; 
	vec2.x = pbeg1->x - pend1->x;
	vec2.y = pbeg1->y - pend1->y;
	area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;

	return area;
}

/*calculates area defined by six points in z=0  plane */
double areaeightpoints(struct vector * p1, struct vector * p2, struct vector * p3, struct vector * p4, 
					 struct vector * p5, struct vector * p6,struct vector * p7, struct vector * p8)
{
	double area =0.0;
	struct vector vec1,vec2;
	
	/*area by half vector cross product  
	 |(pbegining-pbegining)x(pend-pbegining)|/2  */
	vec1.x = p2->x - p1->x; 
	vec1.y = p2->y - p1->y;  
	vec2.x = p3->x - p2->x;
	vec2.y = p3->y - p2->y;
	area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
	//printf("3");
	
	if (p4 != NULL) {
		vec1.x = p3->x - p1->x; 
		vec1.y = p3->y - p1->y; 
		vec2.x = p4->x - p3->x;
		vec2.y = p4->y - p3->y;
		area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
		//printf("4");
		
		if (p5 != NULL) {
			vec1.x = p4->x - p1->x; 
			vec1.y = p4->y - p1->y; 
			vec2.x = p5->x - p4->x;
			vec2.y = p5->y - p4->y;
			area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
			//printf("5");
			
			if (p6 != NULL) {
				vec1.x = p5->x - p1->x; 
				vec1.y = p5->y - p1->y; 
				vec2.x = p6->x - p5->x;
				vec2.y = p6->y - p5->y;
				area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
				//printf("6");
			
				if (p7 != NULL) {
					vec1.x = p6->x - p1->x; 
					vec1.y = p6->y - p1->y; 
					vec2.x = p7->x - p6->x;
					vec2.y = p7->y - p6->y;
					area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
					//printf("7");
										
					if (p8 != NULL) {
						vec1.x = p7->x - p1->x; 
						vec1.y = p7->y - p1->y; 
						vec2.x = p8->x - p7->x;
						vec2.y = p8->y - p7->y;
						area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
						//printf("8");
					}
				}
			}
		}
	}
		
	return area;
}

/*..............................................................................*/
/*........................INPUT OUTPUT..........................................*/
/*..............................................................................*/

/*..............................................................................*/

/**
 * convert string num into two integers
 */
void readii(char * num, int value[2]){
    char *end, *num2;
    void trim (char *);
    
    value[0] = strtol(num, &num2, 10);
    trim(num2);
    if ((int)strlen(num2) > 0) 
	value[1] = strtol(num2, &end, 10);
    else {
	value[1] =0;
	return;
	}
    if(*end){
        fprintf(stderr, "Could not convert %s into two integers\n", num);
        exit(1);
    }
 
    return;
}

/**
 * convert string num into integer
 */
int readi(char * num){
    char *end;
    int i = strtol(num, &end, 10);
    if(*end){
        fprintf(stderr, "Could not convert %s into integer\n", num);
        exit(1);
    }
    return (int) i;
}

/**
 * convert string num into long
 */
long readl(char * num){
    char *end;
    long i = strtol(num, &end, 10);
    if(*end){
        fprintf(stderr, "Could not convert %s into long\n", num);
        exit(1);
    }
    return i;
}

/**
 * convert string num into double
 */
double readd(char * num){
    char *end;
    double i = strtod(num, &end);
    if(*end){
        fprintf(stderr, "Could not convert %s into double\n", num);
        exit(1);
    }
    return i;
}
/*
   Reads the run parameters from the external file "options".  See the end of the
   code for a template. All comments starting with '#' are stripped out.  The 
   options are summarised on standard output and checked for validity of range.
 */
void read_options(struct sim* sim,char filename[30])
{
	int i;
	int num_options = -1;
	double transmx, rotmx, chainmmx, chainrmx;
	double angle, chain_angle;
	char *id, *value, *tokLine, *line;
	FILE *infile;

	void strip_comment (char *);
	void trim (char *);
	void readii(char * num, int value[2]);
	double readd(char *);
	long readl(char *);
	int readi(char *);

    /* for new options add before the last line */
    Option options[] = {
        {"write_cluster",       Long,   FALSE, &sim->write_cluster},
        {"adjust",              Long,   FALSE, &sim->adjust},
        {"movie",               Long,   FALSE, &sim->movie},
        {"nequil",              Long,   FALSE, &sim->nequil},
        {"nsweeps",             Long,   FALSE, &sim->nsweeps},
		{"nrepchange",          Long,   FALSE, &sim->nrepchange},
        {"paramfrq",            Long,   FALSE, &sim->paramfrq},
        {"report",              Long,   FALSE, &sim->report},
        {"seed",                Long,   FALSE, &seed},
        {"pairlist_update",     Int,    FALSE, &sim->pairlist_update},
        {"ptype",               Int,    FALSE, &sim->ptype},
        {"wlm",                 Int2,   FALSE, &sim->wlm},
        {"wlmtype",             Int,    FALSE, &sim->wl.wlmtype},
        {"press",               Double, FALSE, &sim->press},
		{"paralpress",          Double, FALSE, &sim->paralpress},
        {"edge_mx",             Double, FALSE, &sim->edge.mx},
        {"shave",               Double, FALSE, &sim->shave},
        {"chainprob",           Double, FALSE, &sim->chainprob},
        {"switchprob",          Double, FALSE, &sim->switchprob},
        {"temper",              Double, FALSE, &sim->temper},
		{"paraltemper",         Double, FALSE, &sim->paraltemper},
        {"transmx",             Double, FALSE, &transmx},
        {"rotmx",               Double, FALSE, &rotmx},
        {"chainmmx",            Double, FALSE, &chainmmx},
        {"chainrmx",            Double, FALSE, &chainrmx},
        {"last",                Int,    FALSE, NULL}
    };
    while(options[++num_options].var != NULL)
        ;

	/*--- 1. Read in values ---*/
	size_t line_size = (STRLEN + 1) * sizeof(char); 
	line = (char *) malloc(line_size);

	infile = fopen(filename, "r");
	if (infile == NULL) {
		fprintf (stderr, "\nERROR: Could not open options file.\n\n");
		exit (1);
	}

    while(getline(&line, &line_size, infile) != -1){
        // strip comments
        strip_comment(line);
        trim(line);
        if(strlen(line) == 0){
            continue;
        }
        // tokenize
        tokLine = line;
        id = strtok(tokLine, "=");
        if(id == NULL){
            fprintf(stderr, "error parsing Configuration line (%s)", line);
            free(line);
            exit(1);
        }
        trim(id);
        tokLine = NULL;
        value = strtok(tokLine, "=");
        trim(value);
        if(value == NULL){
            fprintf(stderr, "error parsing Configuration line (%s)", line);
            free(line);
            exit(1);
        }
        //printf("id: %s; value: %s\n", id, value);
        for(i = 0; i < num_options; i++){
            if(strcmp(id, options[i].id) == 0){
		if(options[i].type == Int2){
		    readii(value,*((int (*)[2]) options[i].var));
                    options[i].set = TRUE;
                    break;
                }
                if(options[i].type == Int){
                    *((int *) options[i].var) = readi(value);
                    options[i].set = TRUE;
                    break;
                }
                else if(options[i].type == Long){
                    *((long *) options[i].var) = readl(value);
                    options[i].set = TRUE;
                    break;
                }
                else if(options[i].type == Double){
                    *((double *) options[i].var) = readd(value);
                    options[i].set = TRUE;
                    break;
                }
                else {
                    fprintf(stderr, "Could not determine type of %s!\n", id);
                    free(line);
                    exit(1);
                }
            } 
        }
        if(i == num_options){
            fprintf(stderr, "Unknown identifier %s!\nWill procede.\n", id);
        }
    }
	fclose (infile);
    free(line);

    /* Check, wheter all options have been readin */
    for(i = 0; i < num_options; i++){
        if(!options[i].set){
            fprintf(stderr, "option '%s' is not set!\n", options[i].id);
            exit(1);
        }
    }

	/*--- 2. Summarize results on standard output ---*/
	/* Density of close-packed spherocylinders */
	//   rho_cp = 2.0/(sqrt(2.0) + *length * sqrt(3.0));

	printf (" Pressure coupling type:                             %d\n", sim->ptype);
	printf (" Pressure:                                           %.8lf\n", sim->press);
	printf (" Replica exchange pressure:                          %.8lf\n", sim->paralpress);
	printf (" Average volume change attempts per sweep:           %.8lf\n", sim->shave);
	printf (" Equilibration sweeps:                               %ld\n", sim->nequil);
	printf (" Sweeps between step size adjustments:               %ld\n", sim->adjust);
	printf (" Production sweeps:                                  %ld\n", sim->nsweeps);
	printf (" Sweeps between statistics samples:                  %ld\n", sim->paramfrq);
	printf (" Sweeps between statistics reports:                  %ld\n", sim->report);
	printf (" Average chain move attempts per sweep:              %.8lf\n", sim->chainprob);
	printf (" Initial maximum displacement:                       %.8lf\n", transmx);
	printf (" Inititial maximum angular change (degrees):         %.8lf\n", rotmx);
	printf (" Inititial maximum box edge change:                  %.8lf\n", sim->edge.mx);
	printf (" Initial maximum chain displacement:                 %.8lf\n", chainmmx);
	printf (" Inititial maximum chain angular change (degrees):   %.8lf\n", chainrmx);
	printf (" Temperature in kT/e:                                %.8lf\n", sim->temper);
	printf (" Parallel tempering temperature in kT/e:             %.8lf\n", sim->paraltemper);
	printf (" Sweeps between replica exchange:                    %ld\n", sim->nrepchange);
	printf (" Wang-Landau method:                                 %d %d\n", sim->wlm[0],sim->wlm[1]);
	printf (" Calculate the Wang-Landau method for atom type:     %d\n", sim->wl.wlmtype);
	printf (" Average type switch attempts per sweep:             %.8lf\n", sim->switchprob);
	printf (" Number of Sweeps per pairlist update:               %d\n", sim->pairlist_update);
	printf (" Random number seed:                                 %ld\n", seed);
	printf (" Number of sweeps per writing out cluster info:      %ld\n", sim->write_cluster);


	if (sim->movie > 0) {
		printf (" Sweeps between movie frames:                      %ld\n", sim->movie);
	} else {
		printf (" No movie\n");
	}
	printf ("\n");

	if(sim->pairlist_update){
		printf(" A pairlist will be generated every %d steps. This is a greedy"
				" algorithm; make sure you don't have big chains etc.!\n", 
				sim->pairlist_update);
	}

	/*--- 3. Validity checks ---*/
	if (rotmx < 0.0 || rotmx > 180) {
		fprintf (stderr, "ERROR: Maximum orientation change must be in range 0 to 180.\n\n");
		exit (1);
	}
	if (chainrmx < 0.0 || chainrmx > 180) {
		fprintf (stderr, "ERROR: Maximum orientation change for chains must be in range 0 to 180.\n\n");
		exit (1);
	}
	if ( (sim->ptype <0) || (sim->ptype>3) ) {
		fprintf (stderr, "ERROR: Unknown pressure coupling %d. Program only knows: 0 - anisotropic coupling, \
				1 - isotropic coupling, 2 - isotropic in xy z=const, 3 - isotropic xy V=const.\n\n",sim->ptype);
		exit (1);
	}
	if ( (sim->wlm[0] <0) || (sim->wlm[0] > 7) || (sim->wlm[1] <0) || (sim->wlm[1] > 7)  ) {
		fprintf (stderr, "ERROR: Unknown Wang-Landau method %d %d. Program only knows: 0 - none, \
				1 - z-direction od 1st particle, 2 - pore in membrane, 3 - zorientation of 0th particle,\
				4 - distance of fist two particles, 5 - pore around z-axis above CM,\
				6 - pore around z-axis above 0th particle, 7 - number of particles in contact \n\n",sim->wlm[0],sim->wlm[1]);
		exit (1);
	}
	if ( (sim->wlm[0] == 0) && (sim->wlm[1] > 0)  ) {
		fprintf (stderr, "ERROR: Wang-Landau method has to be set for first order parameter and then for second order parameter\n\n");
		exit (1);
	}
	if ( (sim->wlm[0] == 2) || (sim->wlm[0] == 5) || (sim->wlm[0] == 6)  ) { 
		if(sim->wl.wlmtype < 1){
			fprintf (stderr, "ERROR: Atom type for the Wang-Landau Method (%d) was false defined.\n\n",sim->wl.wlmtype);
			exit (1);
		}
		if ( (sim->wlm[1] == 2) || (sim->wlm[1] == 5) || (sim->wlm[1] == 6) ) { 
			fprintf (stderr, "ERROR: Simulaneous use of two pore order parameters has not been implemented yet.\n\n");
			exit (1);
		}
	}

	/* we store maximum rotation as half angle - useful for quaterions*/
	angle = rotmx / 180.0 * PIH *0.5;
	rotmx = cos((rotmx)/180.0*PIH);
	chain_angle = chainrmx / 180.0 * PIH;
	chainrmx = cos((chainrmx)/180.0*PIH);
	sim->edge.mx *= 2.0;   /* The full range is -maxl to +maxl, i.e. spanning 2*maxl */
	transmx *= 2.0;   /* The full range is -maxr to +maxr, i.e. spanning 2*maxr */
	chainmmx *= 2.0;   /* The full range is -maxr to +maxr, i.e. spanning 2*maxr */

	for (i=0;i<MAXT;i++) {
		sim->trans[i].mx = transmx;
		sim->rot[i].mx = rotmx;
		sim->rot[i].angle = angle;
	}
	for (i=0;i<MAXMT;i++) {
		sim->chainm[i].mx = chainmmx;
		sim->chainr[i].mx = chainrmx;
		sim->chainr[i].angle = chain_angle;
	}
	
	//parallel tempering
#ifdef MPI
	if ( (sim->temper != sim->paraltemper) && (sim->mpinprocs <2) ) {
		printf("ERROR: Paralllel tempering at single core does not work.\n\n");
		exit(1);
	}
	sim->dtemp = (sim->paraltemper - sim->temper )/(sim->mpinprocs-1);
	sim->temper += sim->dtemp * sim->mpirank;
	if ( (sim->press != sim->paralpress) && (sim->mpinprocs <2) ) {
		printf("ERROR: Pressure replica exchange at single core does not work.\n\n");
		exit(1);
	}
	sim->dpress = (sim->paralpress - sim->press )/(sim->mpinprocs-1);
	sim->press += sim->dpress * sim->mpirank;
	seed += sim->mpirank;
	sim->mpiexch.mx = sim->dtemp;
	sim->mpiexch.angle = sim->dpress;
#endif

}

/*..............................................................................*/

/*
   Used by read_options to read a long integer with error checking.
   NOT USED ANYMORE
 */

long read_long(FILE *infile, char *error) {

	char *gotline;
	char line[500];
	int fields;
	long value;

	gotline = fgets(line, sizeof(line), infile);
	fields = sscanf(line, "%ld", &value);
	if (gotline == NULL || fields != 1) {
		fprintf (stdout, "\nERROR reading %s from options file.\n\n", error);
		exit (1);
	}

	return value;
}

/*
   Used by read_options to read a long integer with error checking.
   NOT USED ANYMORE
 */

int read_int(FILE *infile, char *error) {

	char *gotline;
	char line[500];
	int fields;
	int value;

	gotline = fgets(line, sizeof(line), infile);
	fields = sscanf(line, "%d", &value);
	if (gotline == NULL || fields != 1) {
		fprintf (stdout, "\nERROR reading %s from options file.\n\n", error);
		exit (1);
	}

	return value;
}



/*..............................................................................*/

/*
   Used by read_options to read a double precision with error checking.
   NOT USED ANYMORE
 */

double read_double(FILE *infile, char *error) {

	char *gotline;
	char line[500];
	int fields;
	double value;

	gotline = fgets(line, sizeof(line), infile);
	fields = sscanf(line, "%le", &value);
	if (gotline == NULL || fields != 1) {
		fprintf (stdout, "\nERROR reading %s from options file.\n\n", error);
		exit (1);
	}

	return value;
}

/*..............................................................................*/

/****************************************************************************
 * CONFIG INITIALIZATION
 *****************************************************************************/
/*
   Reads in the initial configuration from the file "config.init".  Each line
   contains the three components of the position vector and three components of
   the direction vector and three components of patch direction for a spherocylinder.  
   The direction vector is normalised
   after being read in.  The configuration is checked for particle overlaps.
 */

void init_config(struct topo * topo, struct conf * conf, struct sim * sim, char filename[30])
{
	int err,fields,tmp_type;
	long i,j,current,first;
	FILE * infile;
	char * line, line2[STRLEN];
	size_t line_size = (STRLEN + 1) * sizeof(char); 
	line = (char *) malloc(line_size);
	struct particles chorig[MAXCHL];

	int overlap(struct particles, struct particles, struct vector, 
			struct ia_param [MAXT][MAXT]);
	void normalise(struct vector *);
	void ortogonalise(struct vector *, struct vector);
	void usepbc(struct vector *, struct vector);
	double anint(double);
	void strip_comment (char *);
	void trim (char *);
	void aftercommand(char *, char *, char);

	double maxlength = 0;
	for(i = 0; i < MAXT; i++){
		if(maxlength < topo->ia_params[i][i].len[0])
			maxlength = topo->ia_params[i][i].len[0];
	}


	infile = fopen(filename, "r");
	if (infile == NULL) {
		fprintf (stderr, "\nERROR: Could not open config.init file.\n\n");
		exit (1);
	}

	if(getline(&line, &line_size, infile) == -1){
		fprintf (stderr, "ERROR: Could not read box size.\n\n");
		exit (1);
	}
	strip_comment(line);
        trim(line);
	if (sscanf(line, "%le %le %le", &(conf->box.x), &(conf->box.y), &(conf->box.z)) != 3) {
		if(getline(&line, &line_size, infile) == -1){
		      fprintf (stderr, "ERROR: Could not read box size.\n\n");
		      exit (1);
		}
		aftercommand(line2,line,BOXSEP);
		strip_comment(line2);
		trim(line2);
		if (sscanf(line2, "%le %le %le", &(conf->box.x), &(conf->box.y), &(conf->box.z)) != 3) {
		      fprintf (stderr, "ERROR: Could not read box size.\n\n");
		      exit (1);
		}
	}
	if (conf->box.x < maxlength * 2.0 + 2.0) {
		printf ("WARNING: x box length is less than two spherocylinders long.\n\n");
	}
	if (conf->box.y < maxlength * 2.0 + 2.0) {
		printf ("WARNING: y box length is less than two spherocylinders long.\n\n");
	}
	if (conf->box.z < maxlength * 2.0 + 2.0) {
		printf ("WARNING: z box length is less than two spherocylinders long.\n\n");
	}

	DEBUG_INIT("Position of the particle");
	for (i=0; i < topo->npart; i++) {
		if(getline(&line, &line_size, infile) == -1){
			break;
		}
		strip_comment(line);
		trim(line);
		fields = sscanf(line, "%le %le %le %le %le %le %le %le %le %d",
				&conf->particle[i].pos.x, &conf->particle[i].pos.y, &conf->particle[i].pos.z,
				&conf->particle[i].dir.x, &conf->particle[i].dir.y, &conf->particle[i].dir.z,
				&conf->particle[i].patchdir[0].x, &conf->particle[i].patchdir[0].y, &conf->particle[i].patchdir[0].z,
				&conf->particle[i].switched);
		conf->particle[i].patchdir[1].x = conf->particle[i].patchdir[1].y = conf->particle[i].patchdir[1].z =0;
		conf->particle[i].chdir[0].x = conf->particle[i].chdir[0].y = conf->particle[i].chdir[0].z =0;
		conf->particle[i].chdir[1].x = conf->particle[i].chdir[1].y = conf->particle[i].chdir[1].z =0;
		DEBUG_INIT("Line: %s\nNumber of Fields: %d", line, fields);
		if (fields == 9){
			conf->particle[i].switched = 0;
			fprintf(stdout, "WARNING: Particle %ld is assumed to be not switched!\n", i+1);
			fields++;
		}
		if (fields != 10) {
			fprintf (stderr, "ERROR: Could not read coordinates for particle %ld.\n \
					Did you specify box size at the begining?\n\n", i+1);
            free(line);
			exit (1);
		}
		/* Scale position vector to the unit cube */
		usepbc(&conf->particle[i].pos, conf->box );

		conf->particle[i].pos.x /= conf->box.x;
		conf->particle[i].pos.y /= conf->box.y;
		conf->particle[i].pos.z /= conf->box.z;

		if ((topo->ia_params[conf->particle[i].type][conf->particle[i].type].geotype[0]<SP)&&( DOT(conf->particle[i].dir, conf->particle[i].dir) < ZEROTOL )) {
			//DEBUG_INIT("Geotype = %d < %d", conf->particle[i].geotype,SP);
			fprintf (stderr,
					"ERROR: Null direction vector supplied for particle %ld.\n\n", i+1);
            free(line);
			exit (1);
		} else {
			normalise(&conf->particle[i].dir);
		}

		if ((topo->ia_params[conf->particle[i].type][conf->particle[i].type].geotype[0]<SP)&&( DOT(conf->particle[i].patchdir[0], conf->particle[i].patchdir[0]) < ZEROTOL )) {
			fprintf (stderr,
					"ERROR: Null patch vector supplied for particle %ld.\n\n", i+1);
            free(line);
			exit (1);
		} else {
			ortogonalise(&conf->particle[i].patchdir[0],conf->particle[i].dir);
			normalise(&conf->particle[i].patchdir[0]);
		}
		// Switch the type
		if(conf->particle[i].switched){
			if(conf->particle[i].switchtype == 0){
				fprintf(stderr, "ERROR: Particle %ld switched even though it has no switchtype", i);
                free(line);
				exit(1);
			}
			tmp_type = conf->particle[i].type;
			conf->particle[i].type = conf->particle[i].switchtype;
			conf->particle[i].switchtype = tmp_type;
		}

		DEBUG_INIT("%ld:\t%lf\t%lf\t%lf", i, conf->particle[i].pos.x, conf->particle[i].pos.y, conf->particle[i].pos.z);

	}
	free(line);
	/*Make chains WHOLE*/
	for (i=0;i<topo->chainnum;i++){
	    j=0;
	    current = topo->chainlist[i][0];
	    first = current;
	    chorig[0].pos = conf->particle[first].pos;
	    while (current >=0 ) {
			/*shift the chain particle by first one*/
			conf->particle[current].pos.x -= chorig[0].pos.x;
			conf->particle[current].pos.y -= chorig[0].pos.y;
			conf->particle[current].pos.z -= chorig[0].pos.z;
			/*put it in orig box*/
			conf->particle[current].pos.x -=  anint(conf->particle[current].pos.x);
			conf->particle[current].pos.y -=  anint(conf->particle[current].pos.y);
			conf->particle[current].pos.z -=  anint(conf->particle[current].pos.z);
			//printf("ant: %f %f %f\n",conf->particle[current].pos.x,conf->particle[current].pos.y,conf->particle[current].pos.z);
			/*shot it back*/
			conf->particle[current].pos.x += chorig[0].pos.x;
			conf->particle[current].pos.y += chorig[0].pos.y;
			conf->particle[current].pos.z += chorig[0].pos.z;
			//printf("posstart: %f %f %f\n",conf->particle[current].pos.x,conf->particle[current].pos.y,conf->particle[current].pos.z);
			j++;
			current = topo->chainlist[i][j];
	    }
	}
	
	err = 0;
	//for (i=0; i < topo->npart-1; i++) {
	//    for (j=i+1; j < topo->npart; j++) {
	//        if ( overlap(conf->particle[i], conf->particle[j], conf->box, topo->ia_params) ) {
	//            fprintf (stderr,
	//                    "ERROR: Overlap in initial coniguration between particles %ld and %ld.\n",
	//                    i+1, j+1);
	//            err = 1;
	//        }
	//    }
	//}
	if (err) {
		printf ("\n");
		exit (1);
	}

	fclose (infile);

	fflush (stdout);
}
/*..............................................................................*/

/****************************************************************************
 * TOPOLOGY INITIALIZATION
 *****************************************************************************/
/* 
   Create lists for chain operations: Connectivity list where it is written for each sc
   with which sc it is connected. The order is important because spherocylinders have direction
   First is interacting tail then head. Chain list where particles are assigned to chains to 
   which they belong
 */

void init_top(struct topo * topo, struct conf * conf, struct sim * sim,char filename[30])
{
	long i,j,k,mol,maxch,maxpart;
	FILE *infile;
	char *pline=NULL, *dummy=NULL, *sysnames[MAXN];
	char line[STRLEN], keystr[STRLEN], molname[STRLEN];
	unsigned size;
	long  *sysmoln /*[MAXN]*/;
	BOOL exclusions[MAXT][MAXT];

	char *fgets2(char *, int , FILE *);
	void strip_comment (char *);
	void trim(char *);
	int continuing(char *);
	void upstring (char *);
	int filltypes(char **, struct topo * topo);
	int fillexter(char **, struct topo * topo);
	int fillexclusions(char **, BOOL (*exclusions)[MAXT][MAXT]);
	void beforecommand(char *, char *,  char);
	int fillmol(char *, char *, struct molecule * molecules, struct topo * topo);
	int fillsystem(char *, char *[MAXN], long **);
	void initparams(struct topo * topo);
	void genparampairs(struct topo * topo, BOOL (*exclusions)[MAXT][MAXT]);
	int topdealoc(char **, char *[MAXN], long **, struct molecule *);
	struct molecule molecules[MAXMT];

	if ((infile = fopen(filename, "r")) == NULL) {
		fprintf (stderr, "\nTOPOLOGY ERROR: Could not open top.init file.\n\n");
		exit (1);
	}
	fprintf (stdout, "Initialize chainlist...\n");
	fflush(stdout);
	sysmoln = malloc( sizeof(long)*MAXN);
	if(sysmoln == NULL){
		fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sysmoln");
		exit(1);
	}
	struct particles tmp_particles[MAXN];
	for (i=0;i<MAXN;i++) {
		if (i < MAXMT) {
			topo->chainparam[i].bond1eq = -1;
			topo->chainparam[i].bond1c = -1;
			topo->chainparam[i].bond2eq = -1;
			topo->chainparam[i].bond2c = -1;
			topo->chainparam[i].bonddc = -1;
            topo->chainparam[i].angle1eq = -1;
			topo->chainparam[i].angle1c = -1;
			topo->chainparam[i].angle2eq = -1;
			topo->chainparam[i].angle2c = -1;
			molecules[i].name = NULL;
			molecules[i].type = malloc(sizeof(long)*MAXN);
			molecules[i].switchtype = malloc(sizeof(long)*MAXN);
			molecules[i].delta_mu = malloc(sizeof(double)*MAXN);
			for (j=0;j<MAXN;j++) {
				molecules[i].type[j] = -1;
			}
		}
		for (j = 0; j < MAXCHL; j++){
			topo->chainlist[i][j] = -1;
		}
		sysnames[i]=NULL;

	}
	for (i=0;i<MAXT;i++) {
		for (j=0;j<MAXT;j++) {
		      exclusions[i][j]=FALSE;
		}
	}
	topo->exter.exist = FALSE;
	topo->exter.thickness = 0.0;
	topo->exter.epsilon = 0.0;
	topo->exter.attraction = 0.0;
	topo->exter.sqmaxcut = 0.0;
	for(i = 0; i < MAXT; i++){
		for(j = 0; j < MAXT; j++){
			for(k = 0; k < 2; k++){
				topo->ia_params[i][j].geotype[k] = 0;
			}
		}
	}
	fprintf (stdout, "Reading topology...\n");
	fflush(stdout);
	molname[0] = ' ';
	initparams(topo);
	pline=malloc((size_t)STRLEN);
	while (fgets2(line,STRLEN-2,infile) != NULL) {
		strcpy(pline,line);
		if (!pline) fprintf (stderr, "\nTOPOLOGY ERROR: Empty line in topology.\n\n");
		/* build one long line from several fragments */
		while (continuing(line) && (fgets2(line,STRLEN-1,infile) != NULL)) {
			size=strlen(pline)+strlen(line)+1;
			free(pline);
			pline=malloc((size_t)size);
			strcat(pline,line);
		}
		/* skip trailing and leading spaces and comment text */
		strip_comment (pline);
		trim (pline);
		/* if there is something left... */
		if ((int)strlen(pline) > 0) {
			// get the [COMMAND] key
			if (pline[0] == OPENKEY) {
				pline[0] = ' ';
				beforecommand(keystr,pline,CLOSEKEY);
				upstring (keystr);
			} else {
				//DEBUG		fprintf (stdout, "Topology read type:%s, %s \n",keystr,pline);
				if (!strcmp(keystr,"TYPES")) {
					fflush(stdout);
					if (!filltypes(&pline, topo)) {
						DEBUG_INIT("Something went wrong with filltypes");
						fprintf (stderr, "\nTOPOLOGY ERROR: in reading types\n\n");
						topdealoc(&pline,sysnames,&sysmoln, molecules);
						exit (1);
					}
					DEBUG_INIT("back in init_top");
				} else{
					if (!strcmp(keystr,"MOLECULES")){
						DEBUG_INIT("Let's go to the molecules");
						if (molname[0] == ' ') {
							beforecommand(molname,pline,SEPARATOR);
							i=0;
							while (molecules[i].name != NULL)
								i++;
							DEBUG_INIT("in the middle of getting to fillmol");
							molecules[i].name = malloc(strlen(molname)+1);
							strcpy(molecules[i].name, molname);
							fprintf (stdout, "Topology read for molecule: %s \n",molname);
						} 
						if (!fillmol(molname, pline, molecules, topo)) {
							fprintf (stderr, "\nTOPOLOGY ERROR: in reading molecules\n\n");
							topdealoc(&pline,sysnames,&sysmoln, molecules);
							exit (1);
						}
						if ((dummy = strchr (pline,CLOSEMOL)) != NULL) 
							molname[0] = ' ';
					} else {
						if (!strcmp(keystr,"SYSTEM")) {
							if (!fillsystem(pline,sysnames,&sysmoln)) {
								fprintf (stderr, "\nTOPOLOGY ERROR: in reading system\n\n");
								topdealoc(&pline,sysnames,&sysmoln, molecules);
								exit (1);
							}
						} else {
							if (!strcmp(keystr,"EXTER")) {
								fflush(stdout);
								if (!fillexter(&pline, topo)) {
									DEBUG_INIT("Something went wrong with external potential");
									fprintf (stderr, "\nTOPOLOGY ERROR: in reading external potential\n\n");
									topdealoc(&pline,sysnames,&sysmoln, molecules);
									exit (1);
								}
							} else {
								if (!strcmp(keystr,"EXCLUDE")) {
									fflush(stdout);
									if (!fillexclusions(&pline,&exclusions)) {
										DEBUG_INIT("Something went wrong with exclusions potential");
										fprintf (stderr, "\nTOPOLOGY ERROR: in reading exclusions\n\n");
										topdealoc(&pline,sysnames,&sysmoln, molecules);
										exit (1);
									}
								} else {
									fprintf (stderr, "\nTOPOLOGY ERROR: invalid keyword:%s.\n\n", keystr);
									topdealoc(&pline,sysnames,&sysmoln, molecules);
									exit (1);
								}
							}
						}
					}
				}
			}
		}
	}
	/*we have sucessfully read topology*/
	if (pline !=NULL) free(pline);
	pline=NULL;
	fclose (infile);
	fflush (stdout);

	/*fill ia_params combinations*/
	fprintf (stdout, "\nTopology succesfully read. Generating pair interactions...\n");
	genparampairs(topo,&exclusions);
	double maxlength = 0;
	for(i = 0; i < MAXT; i++){
		if(maxlength < topo->ia_params[i][i].len[0])
			maxlength = topo->ia_params[i][i].len[0];
	}
	topo->sqmaxcut += maxlength+2;
	topo->sqmaxcut *= 1.1;
	topo->maxcut = topo->sqmaxcut;
	topo->sqmaxcut = topo->sqmaxcut*topo->sqmaxcut;
	topo->exter.sqmaxcut +=  maxlength;
	topo->exter.sqmaxcut *= topo->exter.sqmaxcut*1.1;

	/*TODO fill chain list and maxch, park particle type*/
	fprintf (stdout, "Generating chainlist...\n");
	maxch=0;
	maxpart=0;
	i=0;
	while (sysnames[i]!=NULL) {
		mol=0;
		while (strcmp(molecules[mol].name,sysnames[i])) {
			mol++;
			if (molecules[mol].name == NULL) {
				fprintf (stderr, "TOPOLOGY ERROR: molecules %s is not defined.\n\n",sysnames[i]);
				topdealoc(&pline,sysnames,&sysmoln, molecules);
				exit(1);
			}
		}
		for (j=0;j<sysmoln[i];j++) {
			//DEBUG	    fprintf (stdout, "molnames %s sysname %s sysnum %ld \n",molnames[mol],sysnames[i],sysmoln[i]);
			k=0;
			while (molecules[mol].type[k] != -1) {
				tmp_particles[maxpart].type = molecules[mol].type[k];
				tmp_particles[maxpart].switchtype = molecules[mol].switchtype[k];
				tmp_particles[maxpart].delta_mu = molecules[mol].delta_mu[k];
				tmp_particles[maxpart].chaint = mol;
				tmp_particles[maxpart].chainn = maxch;
				if (k > MAXCHL) {
					fprintf (stderr, "TOPOLOGY ERROR: more particles in chan (%ld) than allowed(%d).\n",k,MAXCHL);
					fprintf (stderr, "Change MAXCHL in source and recompile the program. \n\n");
					topdealoc(&pline,sysnames,&sysmoln, molecules);
					exit(1);
				}
				if (molecules[mol].type[1] != -1) {
					topo->chainlist[maxch][k] = maxpart;
				}
				k++;
				maxpart++;
				if (maxpart > MAXN) {
					fprintf (stderr, "TOPOLOGY ERROR: more particles(%ld) than allowed(%d).\n",maxpart,MAXN);
					fprintf (stderr, "Change MAXN in source and recompile the program. \n\n");
					topdealoc(&pline,sysnames,&sysmoln, molecules);
					exit(1);
				}
			}
			if (molecules[mol].type[1] != -1) {
				maxch++;
			}
		}
		i++;
	}
	topo->npart = maxpart;
	
	/* write the particles from the temporary to the "permanent" conf */
	conf->particle = malloc(sizeof(struct particles) * topo->npart);
	if(conf->particle == NULL){
		fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for conf->particle");
		exit(1);
	}

	for(i = 0; i < topo->npart; i++){
		conf->particle[i].type       = tmp_particles[i].type;
		conf->particle[i].switchtype = tmp_particles[i].switchtype;
		conf->particle[i].delta_mu   = tmp_particles[i].delta_mu;
		conf->particle[i].chaint     = tmp_particles[i].chaint;
		conf->particle[i].chainn     = tmp_particles[i].chainn;
	}

	/* Initialize the clusterlist */
	sim->clusterlist = malloc(sizeof(long) * topo->npart);
	if(sim->clusterlist == NULL){
		fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sim->clusterlist!");
		exit(1);
	}
	sim->clustersenergy = malloc(sizeof(double) * topo->npart);
	if(sim->clustersenergy== NULL){
		fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sim->clustersenergy!");
		exit(1);
	}
	sim->clusters = NULL;

	/* get all the particles with switch type */
	long switchlist[topo->npart];
	long n_switch_part = 0;
	for(i = 0; i < topo->npart; i++){
		if(conf->particle[i].type != conf->particle[i].switchtype){
			switchlist[n_switch_part] = i;
			n_switch_part++;
		}
	}
	topo->n_switch_part = n_switch_part; 
	if (n_switch_part == 0 && sim->switchprob > 0){
		fprintf(stderr, "TOPOLOGY WARNING: No switchable particles found, but probability for a switch is not zero!\n");
		sim->switchprob = 0;
		fprintf(stderr, "TOPOLOGY WARNING: We changed Switch Probability to zero in this run!\n");
	}
	topo->switchlist=NULL;
	if (n_switch_part > 0){
		topo->switchlist = malloc(sizeof(long) * n_switch_part);
		for(i = 0; i < n_switch_part; i++){
			topo->switchlist[i] = switchlist[i];
			//DEBUG
			//printf("%ld is in switchlist\n", switchlist[i]);
		}
	}

	j = 0;
	while (topo->chainlist[j][0] >= 0) {
		j++;
	}
	topo->chainnum = j;
	if (topo->chainnum != maxch) {
		fprintf (stderr, "TOPOLOGY ERROR: Maximum number of chains(%ld) does not agree with number of chains (%ld)\n\n",maxch,topo->chainnum);
		topdealoc(&pline,sysnames,&sysmoln, molecules);
		exit (1);
	}
	k=0;
	/*clear connectivity and then fill it from chain list*/
	fprintf (stdout, "Generating connectivity...\n");
	for (i=0; i<MAXN; i++) {
		topo->conlist[i][0] = -1;
		topo->conlist[i][1] = -1;
		topo->conlist[i][2] = -1;
		topo->conlist[i][3] = -1;
	}
	conf->sysvolume = 0;
	for (i=0; i<maxpart; i++) {
		for (j=0; j<MAXCHL; j++) {
			if (topo->chainlist[i][j] >= 0) {
				k = topo->chainlist[i][j];
				if ((j+1 < MAXCHL)&&(topo->chainlist[i][j+1] >= 0)) 
					topo->conlist[k][1] = topo->chainlist[i][j+1]; /*if there is a next particle fill it to head bond*/
				if (j > 0) 
					topo->conlist[k][0] = topo->chainlist[i][j-1]; /*if this is not first particle fill tail bond*/
				if ((j+2 < MAXCHL)&& (topo->chainlist[i][j+2] >= 0)) 
					topo->conlist[k][3] = topo->chainlist[i][j+2]; /*if there is a second next particle fill it second neighbour*/
				if (j > 1) 
					topo->conlist[k][2] = topo->chainlist[i][j-2]; /*if this is not second or first particle fill second tail bond*/
			} 
		}
		conf->sysvolume += topo->ia_params[conf->particle[i].type][conf->particle[i].type].volume;
	}
	/*DEBUG
	  for (i=0; i<MAXN; i++) {
	  for (j=0; j<MAXCHL; j++) {
	  fprintf (stderr, " %d",chainlist[i][j]);
	  }
	  fprintf (stderr, " \n");
	  }
	  for (i=0; i<MAXN; i++) {
	  printf (" %ld %ld %ld %ld\n",conlist[i][0],conlist[i][1],conlist[i][2],conlist[i][3]);
	  }
	 */
	
	//  Mark particles as not switched
	for(i = 0; i < maxpart; i++){
		conf->particle[i].switched = 0;
	}
	topdealoc(&pline,sysnames,&sysmoln, molecules);
	DEBUG_INIT("Finished with reading the topology");

	/* Parallel tempering check */
	#ifdef MPI
		// probability to switch replicas = exp ( -0.5 * dT*dT * N / (1 + dT) )
		printf("Probability to switch replicas is roughly: %f\n",exp(-0.5 * maxpart * sim->dtemp * sim->dtemp / (1.0 + sim->dtemp)) );
	#endif
}


/*..........................................................................*/
/*dealocting memory for init_top*/
int topdealoc(char **pline,char *sysnames[MAXN], long **sysmoln, struct molecule * molecules)
{
	long i;

	if ((*pline) != NULL) free((*pline));
	(*pline)=NULL;
	if ((*sysmoln) != NULL) free((*sysmoln));
	(*sysmoln)=NULL;
	for (i=0;i<MAXN;i++) {
		if (i < MAXMT) {
			free(molecules[i].name);
			free(molecules[i].type);
			free(molecules[i].switchtype);
			free(molecules[i].delta_mu);
		}
		if ((sysnames[i]) != NULL) free(sysnames[i]);
		sysnames[i]=NULL;
	}

	return 0;
}

/* initiate vectors of a single particle*/
void int_partvec(long target, struct ia_param * ia_parami, struct conf * conf )
{
	struct quat quatrot;
	
	struct quat quat_create(struct vector, double, double);
	void vec_rotate(struct vector *, struct quat);
	void normalise(struct vector *);
	void ortogonalise(struct vector *,struct vector);
	
	if ( (ia_parami->geotype[0] == SCA) || (ia_parami->geotype[0] == SCN) ){
	    /*SCA and SCN are isotropic... nothing to initialize*/
	    return;
	}
	normalise (&conf->particle[target].dir);
	ortogonalise(&conf->particle[target].patchdir[0],conf->particle[target].dir);
	/*calculate patch sides*/
	if ( (ia_parami->geotype[0] == PSC) || (ia_parami->geotype[0] == CPSC) 
	  || (ia_parami->geotype[0] == TPSC) || (ia_parami->geotype[0] == TCPSC)  ){
		/* rotate patch vector by half size of patch*/
		conf->particle[target].patchsides[0] = conf->particle[target].patchdir[0];
		quatrot=quat_create(conf->particle[target].dir, ia_parami->pcoshalfi[0], ia_parami->psinhalfi[0]);
		vec_rotate(&(conf->particle[target].patchsides[0]),quatrot);
		/*second side*/
		conf->particle[target].patchsides[1] = conf->particle[target].patchdir[0];
		quatrot=quat_create(conf->particle[target].dir, ia_parami->pcoshalfi[0], -1.0*ia_parami->psinhalfi[0]);
		vec_rotate(&(conf->particle[target].patchsides[1]),quatrot);
	}
	/*calculate second patchdir*/
	if ( (ia_parami->geotype[0] == TPSC) || (ia_parami->geotype[0] == TCPSC) || 
	  (ia_parami->geotype[0] == TCHPSC) || (ia_parami->geotype[0] == TCHCPSC)){
		conf->particle[target].patchdir[1] = conf->particle[target].patchdir[0];
		quatrot=quat_create(conf->particle[target].dir, ia_parami->csecpatchrot[0], ia_parami->ssecpatchrot[0]);
		vec_rotate(&(conf->particle[target].patchdir[1]),quatrot);
		ortogonalise(&conf->particle[target].patchdir[1],conf->particle[target].dir);
	}
	/*calculate second patch sides*/
	if ( (ia_parami->geotype[0] == TPSC) || (ia_parami->geotype[0] == TCPSC) ){
		/* rotate patch vector by half size of patch*/
		conf->particle[target].patchsides[2] = conf->particle[target].patchdir[1];
		quatrot=quat_create(conf->particle[target].dir, ia_parami->pcoshalfi[2], ia_parami->psinhalfi[2]);
		vec_rotate(&(conf->particle[target].patchsides[2]),quatrot);
		/*second side*/
		conf->particle[target].patchsides[3] = conf->particle[target].patchdir[1];
		quatrot=quat_create(conf->particle[target].dir, ia_parami->pcoshalfi[2], -1.0*ia_parami->psinhalfi[2]);
		vec_rotate(&(conf->particle[target].patchsides[3]),quatrot);
	}
	
	/*calculate chdir vector*/
	if ( (ia_parami->geotype[0] == CHPSC) || (ia_parami->geotype[0] == CHCPSC) 
	  || (ia_parami->geotype[0] == TCHPSC) || (ia_parami->geotype[0] == TCHCPSC)){
		conf->particle[target].chdir[0] = conf->particle[target].dir;
		quatrot = quat_create(conf->particle[target].patchdir[0], ia_parami->chiral_cos[0], ia_parami->chiral_sin[0]);
		vec_rotate(&(conf->particle[target].chdir[0]), quatrot);	
		/* rotate patch vector by half size of patch*/
		conf->particle[target].patchsides[0] = conf->particle[target].patchdir[0];
		quatrot=quat_create(conf->particle[target].chdir[0], ia_parami->pcoshalfi[0], ia_parami->psinhalfi[0]);
		vec_rotate(&(conf->particle[target].patchsides[0]),quatrot);
		/*second side*/
		conf->particle[target].patchsides[1] = conf->particle[target].patchdir[0];
		quatrot=quat_create(conf->particle[target].chdir[0], ia_parami->pcoshalfi[0], -1.0*ia_parami->psinhalfi[0]);
		vec_rotate(&(conf->particle[target].patchsides[1]),quatrot);
	}
	/*calculate chdir vector for seond patch*/
	if ( (ia_parami->geotype[0] == TCHPSC) || (ia_parami->geotype[0] == TCHCPSC) ){
		conf->particle[target].chdir[1] = conf->particle[target].dir;
		quatrot = quat_create(conf->particle[target].patchdir[1], ia_parami->chiral_cos[0], ia_parami->chiral_sin[0]);
		vec_rotate(&(conf->particle[target].chdir[1]), quatrot);	
		/* rotate patch vector by half size of patch to get sides*/
		conf->particle[target].patchsides[2] = conf->particle[target].patchdir[1];
		quatrot=quat_create(conf->particle[target].chdir[1], ia_parami->pcoshalfi[2], ia_parami->psinhalfi[2]);
		vec_rotate(&(conf->particle[target].patchsides[2]),quatrot);
		/*second side*/
		conf->particle[target].patchsides[3] = conf->particle[target].patchdir[1];
		quatrot=quat_create(conf->particle[target].chdir[1], ia_parami->pcoshalfi[2], -1.0*ia_parami->psinhalfi[2]);
		vec_rotate(&(conf->particle[target].patchsides[3]),quatrot);
	}	
	
}

/* calculate vectors on particles for speedup*/
void partvecinit(struct topo * topo, struct sim * sim, struct conf * conf )
{
	long i;
	void int_partvec(long target, struct ia_param *, struct conf * conf );
	
	for(i = 0; i < topo->npart; i++){
	    if ( topo->ia_params[conf->particle[i].type][conf->particle[i].type].geotype[0]  < SP)
			int_partvec(i,&(topo->ia_params[conf->particle[i].type][conf->particle[i].type]),conf);
	}

}

/*generate interations pairs*/

void genparampairs(struct topo * topo, BOOL (*exclusions)[MAXT][MAXT])
{
	int i,j,k;
	int a[2];
	int len;
	double length = 0; // The length of a PSC, currently only one is allow, ie implemented

	for (i=0;i<MAXT;i++) {
		for (j=0;j<MAXT;j++) {
			if (i!=j) {
				if((topo->ia_params[j][j].geotype[0] != 0) && (topo->ia_params[i][i].geotype[0] != 0)){
					a[0] = i;
					a[1] = j;
					for(k = 0; k < 2; k++){
						topo->ia_params[i][j].geotype[k] = topo->ia_params[a[k]][a[k]].geotype[0];
						topo->ia_params[i][j].len[k] = topo->ia_params[a[k]][a[k]].len[0];
						if (topo->ia_params[a[k]][a[k]].len[0] > 0){
							if (length == 0){
								length = topo->ia_params[a[k]][a[k]].len[0];
							}
							else if (length > 0){
								if (length != topo->ia_params[a[k]][a[k]].len[0]){
									fprintf(stderr, "Error: ");
									fprintf(stderr, "Different lengths for spherocylinders have not been implemented yet!\n");
									fprintf(stderr, "\tCheck the length of type %d!\n", a[k]); 
									exit(1);
								}
							}
						}
						topo->ia_params[i][j].half_len[k] = topo->ia_params[a[k]][a[k]].half_len[0];
						/* Handle angles only, when geotype is a patchs sphero cylinder */
						if(topo->ia_params[i][j].geotype[k] >= PSC && topo->ia_params[i][j].geotype[k] < SP){
							topo->ia_params[i][j].pangl[k] = topo->ia_params[a[k]][a[k]].pangl[0];
							topo->ia_params[i][j].panglsw[k] = topo->ia_params[a[k]][a[k]].panglsw[0];
							topo->ia_params[i][j].pcangl[k] = cos(topo->ia_params[i][j].pangl[k]/2.0/180*PI);
							topo->ia_params[i][j].pcanglsw[k] = cos((topo->ia_params[i][j].pangl[k]/2.0+topo->ia_params[i][j].panglsw[k])/180*PI);
							topo->ia_params[i][j].pcoshalfi[k] = cos((topo->ia_params[i][j].pangl[k]/2.0+topo->ia_params[i][j].panglsw[k])/2.0/180*PI);
							topo->ia_params[i][j].psinhalfi[k] = sqrt(1.0 - topo->ia_params[i][j].pcoshalfi[k] * topo->ia_params[i][j].pcoshalfi[k]);
							
						}
						
						/* Only when the PSC is chiral */
						if( (topo->ia_params[i][j].geotype[k] == CHCPSC) || (topo->ia_params[i][j].geotype[k] == CHPSC) \
						   || (topo->ia_params[i][j].geotype[k] == TCHCPSC) || (topo->ia_params[i][j].geotype[k] == TCHPSC) ){
							topo->ia_params[i][j].chiral_cos[k] = topo->ia_params[a[k]][a[k]].chiral_cos[0];
							topo->ia_params[i][j].chiral_sin[k] = topo->ia_params[a[k]][a[k]].chiral_sin[0];
						}
						/* Information of two patches */
						if( (topo->ia_params[i][j].geotype[k] == TCPSC) || (topo->ia_params[i][j].geotype[k] == TPSC) \
						   || (topo->ia_params[i][j].geotype[k] == TCHCPSC) || (topo->ia_params[i][j].geotype[k] == TCHPSC) ){
							topo->ia_params[i][j].csecpatchrot[k] = topo->ia_params[a[k]][a[k]].csecpatchrot[0];
							topo->ia_params[i][j].ssecpatchrot[k] = topo->ia_params[a[k]][a[k]].ssecpatchrot[0];
							
							topo->ia_params[i][j].pangl[k+2] = topo->ia_params[a[k]][a[k]].pangl[2];
							topo->ia_params[i][j].panglsw[k+2] = topo->ia_params[a[k]][a[k]].panglsw[2];
							topo->ia_params[i][j].pcangl[k+2] = cos(topo->ia_params[i][j].pangl[k+2]/2.0/180*PI);
							topo->ia_params[i][j].pcanglsw[k+2] = cos((topo->ia_params[i][j].pangl[k+2]/2.0+topo->ia_params[i][j].panglsw[k+2])/180*PI);
							topo->ia_params[i][j].pcoshalfi[k+2] = cos((topo->ia_params[i][j].pangl[k+2]/2.0+topo->ia_params[i][j].panglsw[k+2])/2.0/180*PI);
							topo->ia_params[i][j].psinhalfi[k+2] = sqrt(1.0 - topo->ia_params[i][j].pcoshalfi[k+2] * topo->ia_params[i][j].pcoshalfi[k+2]);
							
						}
					}
					len = strlen(topo->ia_params[i][i].name);
					strncpy(topo->ia_params[i][j].name, topo->ia_params[i][i].name, len + 1);
					len = strlen(topo->ia_params[i][i].other_name);
					strncpy(topo->ia_params[i][j].other_name, topo->ia_params[i][i].other_name, len + 1);
					topo->ia_params[i][j].sigma = AVER(topo->ia_params[i][i].sigma,topo->ia_params[j][j].sigma);
					topo->ia_params[i][j].epsilon = sqrt(topo->ia_params[i][i].epsilon *  topo->ia_params[j][j].epsilon);
					topo->ia_params[i][j].pswitch = AVER(topo->ia_params[i][i].pswitch,topo->ia_params[j][j].pswitch);

					topo->ia_params[i][j].rcutwca = (topo->ia_params[i][j].sigma)*pow(2.0,1.0/6.0);
					// Averaging of the flat part of attraction
					topo->ia_params[i][j].pdis = AVER(topo->ia_params[i][i].pdis - topo->ia_params[i][i].rcutwca, \
							topo->ia_params[j][j].pdis - topo->ia_params[j][j].rcutwca) + topo->ia_params[i][j].rcutwca;
					topo->ia_params[i][j].rcut = topo->ia_params[i][j].pswitch+topo->ia_params[i][j].pdis;

					// if not non-attractive == if attractive
					if (!((topo->ia_params[i][j].geotype[0] % 10 == 0) || (topo->ia_params[i][j].geotype[1] % 10 == 0))){
						if (topo->ia_params[i][j].rcutwca >  topo->ia_params[i][j].rcut){
							fprintf(stderr, "Error: Repulsive cutoff is larger than the attractive cutoff!\n");
							fprintf(stderr, "       between %d and %d: %lf > %lf\n", i, j, topo->ia_params[i][j].rcutwca, topo->ia_params[i][j].rcut);
						}
					}

					if ( topo->ia_params[i][j].rcutwca > topo->sqmaxcut )
						topo->sqmaxcut = topo->ia_params[i][j].rcutwca;
					if ( topo->ia_params[i][j].rcut > topo->sqmaxcut )
						topo->sqmaxcut = topo->ia_params[i][j].rcut;
					
				}
			}
		}
		/*filling interaction with external potential*/
		if( (topo->exter.exist) && (topo->ia_params[i][i].geotype[0] != 0)){
			/*use everything like for given particles except distance and attraction, which is generated as for other interactions*/
			topo->exter.interactions[i] = topo->ia_params[i][i];
			topo->exter.interactions[i].sigma = AVER(topo->ia_params[i][i].sigma, topo->exter.thickness);
			topo->exter.interactions[i].rcutwca = (topo->exter.interactions[i].sigma)*pow(2.0,1.0/6.0);
			topo->exter.interactions[i].epsilon = sqrt(topo->ia_params[i][i].epsilon *  topo->exter.epsilon);
			topo->exter.interactions[i].pswitch = AVER(topo->ia_params[i][i].pswitch, topo->exter.attraction);
			topo->exter.interactions[i].pdis = AVER(topo->ia_params[i][i].pdis - topo->ia_params[i][i].rcutwca, 0.0) + topo->exter.interactions[i].rcutwca;
			topo->exter.interactions[i].rcut = topo->exter.interactions[i].pswitch + topo->exter.interactions[i].pdis;
			if (topo->exter.interactions[i].rcut > topo->exter.sqmaxcut ) topo->exter.sqmaxcut = topo->exter.interactions[i].rcut;
		}
	}
	for (i=0;i<MAXT;i++) {
		for (j=0;j<MAXT;j++) {
			if ( (*exclusions)[i][j] )
				topo->ia_params[i][j].epsilon = 0.0;
		}
	}
}

/*initialize parameters for interactions*/

void initparams(struct topo * topo)
{
	int i,j,k;

	for (i=0;i<MAXT;i++) {
		for (j=0;j<MAXT;j++) {
			for(k = 0; k < 2; k++){
				topo->ia_params[i][j].geotype[k] = 0;
				topo->ia_params[i][j].len[k] = 0.0;
				topo->ia_params[i][j].half_len[k] = 0.0;
				topo->ia_params[i][j].chiral_cos[k] = 0.0;
				topo->ia_params[i][j].chiral_sin[k] = 0.0;
				topo->ia_params[i][j].csecpatchrot[k] = 0.0;
				topo->ia_params[i][j].ssecpatchrot[k] = 0.0;
			}
			for(k = 2; k < 4; k++){
				topo->ia_params[i][j].pangl[k] = 0.0;
				topo->ia_params[i][j].panglsw[k] = 0.0;
				topo->ia_params[i][j].pcangl[k] = 0.0;
				topo->ia_params[i][j].pcanglsw[k] = 0.0;
				topo->ia_params[i][j].pcoshalfi[k] = 0.0;
				topo->ia_params[i][j].psinhalfi[k] = 0.0;
			}
			topo->ia_params[i][j].sigma = 0.0;
			topo->ia_params[i][j].epsilon = 0.0;
			topo->ia_params[i][j].rcutwca = 0.0;
			topo->ia_params[i][j].pdis = 0.0;
			topo->ia_params[i][j].pswitch = 0.0;
			topo->ia_params[i][j].rcut = 0.0;
			topo->ia_params[i][j].volume = 0.0;
			topo->ia_params[i][j].pvolscale = 0.0;
		}
	}
	topo->sqmaxcut = 0;
}
/*...........................................................................*/

/*filling the system parameters*/
int fillsystem(char *pline, char *sysnames[MAXN], long **sysmoln)
{
	int i,fields;
	char zz[STRLEN];

	void trim (char *);

	trim(pline);
	if (!pline) {
		fprintf (stderr, "TOPOLOGY ERROR: obtained empty line in fil system.\n\n");
		return 0;
	}
	i=0;
	while (sysnames[i]!=NULL) i++;
	fields = sscanf(pline, "%s %ld", zz, &(*sysmoln)[i]);
	sysnames[i]=malloc(strlen(zz)+1);
	strcpy(sysnames[i],zz);

	if (fields != 2) {
		fprintf (stderr, "TOPOLOGY ERROR: failed reading system from (%s).\n\n", pline);
		return 0;
	}
	if ((*sysmoln)[i] < 1) {
		fprintf (stderr, "TOPOLOGY ERROR: cannot have %ld number of molecules.\n\n", (*sysmoln)[i]);
		return 0;
	}
	fprintf (stdout, "system: %s %ld\n",sysnames[i],(*sysmoln)[i]);
	return 1;
}

/*filling the parameters for molecules*/
int fillmol(char *molname, char *pline, struct molecule * molecules, struct topo * topo)
{
	DEBUG_INIT("fillmol just has been called!");
	char str[STRLEN],str2[STRLEN],molcommand[STRLEN],molparams[STRLEN];
	int i,j,fields;
	double bondk,bonddist;

	void trim (char *);
	void upstring(char *);
	void beforecommand(char *, char *, char);
	void aftercommand(char *, char *, char);

	beforecommand(str2, pline, CLOSEMOL);
	aftercommand(str, str2, OPENMOL);
	trim(str);
	if (strlen(str) == 0) return 1;
	beforecommand(molcommand,str,SEPARATOR);
	aftercommand(molparams,str,SEPARATOR);
	trim(molcommand);
	trim(molparams);
	upstring (molcommand);
	DEBUG_INIT("molcommand: %s", molcommand);
	DEBUG_INIT("molparams: %s", molparams);
	i=0;
	while (strcmp(molecules[i].name, molname)) 
		i++;
	j=0;
	while (molecules[i].type[j] != -1)
		j++;
	if (!strcmp(molcommand,"PARTICLES")) {
		fprintf (stdout, "particle %d: \t", j + 1);
		fields =  sscanf(molparams,"%ld %ld %lf",molecules[i].type + j, molecules[i].switchtype + j, molecules[i].delta_mu + j);
		fprintf (stdout, "%ld ",molecules[i].type[j]);
		if (fields == 1){
			(molecules[i].switchtype[j]) = (molecules[i].type[j]);
			(molecules[i].delta_mu[j]) = 0;
			fields = 3;
		} else{
			fprintf(stdout, "(with switchtype: %ld and delta_mu: %lf)", molecules[i].switchtype[j], molecules[i].delta_mu[j]);
		}
		if (fields != 3) {
			fprintf (stderr, "TOPOLOGY ERROR: could not read a pacticle.\n\n");
			return 0;
		}
		fflush(stdout);
		if (molecules[i].type[j] < 0) {
			fprintf (stderr, "TOPOLOGY ERROR: pacticles include negative type.\n\n");
			return 0;
		}
		if (molecules[i].type[j] > MAXT) {
			fprintf (stderr, "TOPOLOGY ERROR: pacticles include type out of range 0-%ld.\n\n",(long)MAXT);
			return 0;
		}
		fprintf (stdout, "\n");
		return 1;
	}
	if (!strcmp(molcommand,"BOND1")) {
		fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
		if (fields < 2) {
			fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for bond1, should be 2.\n\n");
			return 0;
		}
		if (bonddist < 0) {
			fprintf (stderr, "TOPOLOGY ERROR: bonddist cannot be negative: %f \n\n",bonddist);
			return 0;
		}
		topo->chainparam[i].bond1c = bondk;
		topo->chainparam[i].bond1eq = bonddist;
		fprintf (stdout, "bond1: %f %f \n",topo->chainparam[i].bond1c,topo->chainparam[i].bond1eq);
		return 1;
	}
	if (!strcmp(molcommand,"BOND2")) {
		fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
		if (fields < 2) {
			fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for bond2, should be 2.\n\n");
			return 0;
		}
		if (bonddist < 0) {
			fprintf (stderr, "TOPOLOGY ERROR: bonddist cannot be negative: %f \n\n",bonddist);
			return 0;
		}
		topo->chainparam[i].bond2c = bondk;
		topo->chainparam[i].bond2eq = bonddist;
		fprintf (stdout, "bond2: %f %f \n",topo->chainparam[i].bond2c,topo->chainparam[i].bond2eq);
		return 1;
	}
    if (!strcmp(molcommand,"BONDD")) {
		fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
		if (fields < 2) {
			fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for bondd, should be 2.\n\n");
			return 0;
		}
		if (bonddist < 0) {
			fprintf (stderr, "TOPOLOGY ERROR: bonddist cannot be negative: %f \n\n",bonddist);
			return 0;
		}
		topo->chainparam[i].bonddc = bondk;
		topo->chainparam[i].bonddeq = bonddist;
		fprintf (stdout, "bondd: %f %f \n",topo->chainparam[i].bonddc,topo->chainparam[i].bonddeq);
		return 1;
	}
    
    if (!strcmp(molcommand,"ANGLE1")) {
		fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
		if (fields < 2) {
			fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for angle1, should be 2.\n\n");
			return 0;
		}
		if (bonddist < 0) {
			fprintf (stderr, "TOPOLOGY ERROR: equilibrium angle cannot be negative: %f \n\n",bonddist);
			return 0;
		}
		topo->chainparam[i].angle1c = bondk;
		topo->chainparam[i].angle1eq = bonddist/180.0*PI;
		fprintf (stdout, "angle1: %f %f \n",topo->chainparam[i].angle1c,topo->chainparam[i].angle1eq);
		return 1;
	}
    if (!strcmp(molcommand,"ANGLE2")) {
		fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
		if (fields < 2) {
			fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for angle2, should be 2.\n\n");
			return 0;
		}
		if (bonddist < 0) {
			fprintf (stderr, "TOPOLOGY ERROR: equilibrium angle cannot be negative: %f \n\n",bonddist);
			return 0;
		}
		topo->chainparam[i].angle2c = bondk;
		topo->chainparam[i].angle2eq = bonddist/180.0*PI;
		fprintf (stdout, "angle2: %f %f \n",topo->chainparam[i].angle2c,topo->chainparam[i].angle2eq);
		return 1;
	}
    
	fprintf (stderr, "TOPOLOGY ERROR: unknown parameter: %s.\n\n",molcommand);
	return 0;
}

/* Converts the geometrical type string into a number */
int convert_geotype(char * geotype){
	if (strcmp(geotype, "CPSC") == 0)
		return CPSC;
	if (strcmp(geotype, "CHCPSC") == 0)
		return CHCPSC;
	if (strcmp(geotype, "SCA") == 0)
		return SCA;
	if (strcmp(geotype, "PSC") == 0)
		return PSC;
	if (strcmp(geotype, "CHPSC") == 0)
		return CHPSC;
	if (strcmp(geotype, "TCPSC") == 0)
		return TCPSC;
	if (strcmp(geotype, "TCHCPSC") == 0)
		return TCHCPSC;
	if (strcmp(geotype, "TPSC") == 0)
		return TPSC;
	if (strcmp(geotype, "TCHPSC") == 0)
		return TCHPSC;
	if (strcmp(geotype, "SPN") == 0)
		return SPN;
	if (strcmp(geotype, "SPA") == 0)
		return SPA;

	return 0;
}

/*filling the parameters of external potentail - wall. Returns 1 on succes.*/
int fillexter(char **pline, struct topo * topo)
{
	int fields;
	
	double param[3];
	/* 0: thickness
	 * 1: epsilon
	 * 2: attraction 
	 */
	char typestr[STRLEN], paramstr[STRLEN];

	void trim (char *);
	void beforecommand(char *, char *, char);
	void aftercommand(char *, char *, char);

	beforecommand(typestr, *pline, SEPARATOR);
	aftercommand(paramstr, *pline, SEPARATOR);
	fields = sscanf(paramstr, "%le %le %le", &param[0], &param[1], &param[2]);
	if (fields >3) {
		fprintf (stderr, "TOPOLOGY ERROR: too many parameters for external potential. We have \
				thickness, epsilon, and attraction distance so far.\n\n");
		return 0;
	}
	if (fields >0) {
		topo->exter.exist = TRUE;
		topo->exter.thickness = param[0];
		fprintf(stdout, "External potential with thickness: %le ",topo->exter.thickness);
		if (fields >1) {
			topo->exter.epsilon = param[1];
			fprintf(stdout, "epsilon: %le ",topo->exter.epsilon);
			if (fields >2) {
				topo->exter.attraction = param[2];
				fprintf(stdout, "and range of attraction: %le ",topo->exter.attraction);
			}
		}
	} else{
		topo->exter.exist = FALSE;
		fprintf(stdout, "No external potential ");
	}

	fprintf(stdout, " \n");
	DEBUG_INIT("Finished filling external potential");
	return 1;
}

/*filling pair for which we exlude attraction interaction. Returns 1 on succes.*/
int fillexclusions(char **pline, BOOL (*exlusions)[MAXT][MAXT])
{
	long num1,num2;
	char *pline1, *pline2;

	void trim (char *);

	num1 = strtol(*pline, &pline2, 10);
	trim(pline2);
	if ((int)strlen(pline2) > 0) {
	    num2 = strtol(pline2, &pline1, 10);
	    trim(pline1);
	    (*exlusions)[num1][num2]=TRUE;
	    (*exlusions)[num2][num1]=TRUE;
	    fprintf(stderr, "Exclusions %ld %ld \n", num1, num2);
	} else {
	    fprintf(stderr, "Error in readin Topology exclusions, probably there is not even number of types \n");
	    return 0;
	}
	while ((int)strlen(pline1) > 0) {
	  num1 = strtol(pline1, &pline2, 10);
	  trim(pline2);
	  if ((int)strlen(pline2) > 0) {
	    num2 = strtol(pline2, &pline1, 10);
	    trim(pline1);
	    (*exlusions)[num1][num2]=TRUE;
	    (*exlusions)[num2][num1]=TRUE;
	    fprintf(stderr, "Exclusions %ld %ld \n", num1, num2);
	  } else {
	    fprintf(stderr, "Error in readin Topology exclusions, probably there is not even number of types \n");
	    return 0;
	  }
	}
	return 1;
}

/*filing the parameters for types from given strings. Returns 1 on succes.*/
int filltypes(char **pline, struct topo * topo)
{
	int type;
	int geotype_i;
	int fields;
	char name[SMSTR];
	char geotype[SMSTR];

	double param[11];
	/* 0: epsilon
	 * 1: sigma
	 * 2: attraction dist
	 * 3: sttraction switch
	 * 4: patch angle
	 * 5: patch switch
	 * 6: length
	 * 7(optional): second patche rotation
	 * 8(optional): second patch angle
	 * 9(optional): second patch angle switch
	 * +1: chirality
	 */
	char typestr[STRLEN], paramstr[STRLEN];

	void trim (char *);
	void beforecommand(char *, char *, char);
	void aftercommand(char *, char *, char);

	beforecommand(typestr, *pline, SEPARATOR);
	aftercommand(paramstr, *pline, SEPARATOR);
	fields = sscanf(paramstr, "%s %d %s %le %le %le %le %le %le %le %le %le %le %le", name, &type, geotype, &param[0], &param[1], &param[2], &param[3], &param[4], &param[5], &param[6], &param[7], &param[8], &param[9], &param[10]);
	fields -= 5; // number of parameter fields => I am too lazy to adjust everywhere below the numbers
	//DEBUG    fprintf (stdout, "Topology read geotype: %ld with parameters fields %d, str:%s and %s in pline %s\n",geotype,fields,geotypestr,paramstr,pline);


	geotype_i = convert_geotype(geotype);
	if(!geotype_i){
		fprintf(stderr, "TOPOLOGY ERROR: Unknown GEOTYPE: %s!", geotype);
		return 0;
	}
	DEBUG_INIT("geotype_i: %d; fields = %d", geotype_i, fields);
	if (( (geotype_i == SCN) || (geotype_i == SPN) ) && (fields != 0)) {
		fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 1.\n\n", geotype);
		return 0;
	}
	if (( (geotype_i == SCA) || (geotype_i == SPA)) && (fields != 2)) {
		fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 3.\n\n", geotype);
		return 0;
	}
	if (( (geotype_i == PSC) || (geotype_i == CPSC) ) && (fields != 5)) {
		fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 5.\n\n", geotype);
		return 0;
	}
	if (( (geotype_i == CHCPSC) || (geotype_i == CHCPSC) )&& ( fields != 6)) {
		fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 6.\n\n", geotype);
		return 0;
	}
	if (( (geotype_i == TPSC) || (geotype_i == TCPSC) ) && (fields != 8)) {
		fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 8.\n\n", geotype);
		return 0;
	}
	if (( (geotype_i == TCHCPSC) || (geotype_i == TCHCPSC) )&& ( fields != 9)) {
		fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 9.\n\n", geotype);
		return 0;
	}
	
	if ((geotype_i < 0) || (geotype_i > (MAXT + 10))) {
		fprintf (stderr, "TOPOLOGY ERROR: geotype (%s) is out of range: 0 - %d.\n\n", geotype, MAXT + 10);
		return 0;
	}

	strcpy(topo->ia_params[type][type].name, name);
	strcpy(topo->ia_params[type][type].other_name, name);
	topo->ia_params[type][type].geotype[0] = geotype_i;
	topo->ia_params[type][type].geotype[1] = geotype_i;
	topo->ia_params[type][type].epsilon = param[0];
	topo->ia_params[type][type].sigma = param[1];
	topo->ia_params[type][type].rcutwca = (topo->ia_params[type][type].sigma)*pow(2.0,1.0/6.0);
	fprintf(stdout, "Topology read of %d: %s (geotype: %s, %d) with parameters %lf %lf", type, name, geotype, geotype_i, topo->ia_params[type][type].epsilon, topo->ia_params[type][type].sigma);
	if (fields > 0) {
		topo->ia_params[type][type].pdis = param[2];
		topo->ia_params[type][type].pswitch = param[3];
		topo->ia_params[type][type].rcut = topo->ia_params[type][type].pswitch+topo->ia_params[type][type].pdis;
		fprintf(stdout, " %f %f",topo->ia_params[type][type].pdis,topo->ia_params[type][type].pswitch);
	}
	if (fields > 2) {
		int i;
		for(i = 0; i < 2; i++){
			topo->ia_params[type][type].len[i] = param[6];
			topo->ia_params[type][type].half_len[i] = param[6] / 2;
			topo->ia_params[type][type].pangl[i] = param[4];
			topo->ia_params[type][type].panglsw[i] = param[5];
			topo->ia_params[type][type].pcangl[i] = cos(param[4]/2.0/180*PI);                 // C1
			topo->ia_params[type][type].pcanglsw[i] = cos((param[4]/2.0+param[5])/180*PI);    // C2
			//topo->ia_params[type][type].pcangl[i] = topo->ia_params[type][type].pcangl[i];
			//topo->ia_params[type][type].pcanglsw[i] =	topo->ia_params[type][type].pcanglsw[i];
			topo->ia_params[type][type].pcoshalfi[i] = cos((param[4]/2.0+param[5])/2.0/180*PI);
			topo->ia_params[type][type].psinhalfi[i] = sqrt(1.0 - topo->ia_params[type][type].pcoshalfi[i] * topo->ia_params[type][type].pcoshalfi[i]);
		}
		fprintf(stdout, " %f %f", topo->ia_params[type][type].pangl[0], topo->ia_params[type][type].panglsw[0]);
	}
	if(fields == 6){
		int i;
		for(i = 0; i < 2; i++){
			topo->ia_params[type][type].chiral_cos[i] = cos(param[7] / 360 * PI);
			topo->ia_params[type][type].chiral_sin[i] = sqrt(1 - topo->ia_params[type][type].chiral_cos[i] * topo->ia_params[type][type].chiral_cos[i]);
			fprintf(stdout, " %f ", param[7]);
		}
	}
	if ((fields == 8)||(fields == 9)) {
		int i;
		for(i = 0; i < 2; i++){
			topo->ia_params[type][type].csecpatchrot[i] = cos(param[7] / 360 * PI);
			topo->ia_params[type][type].ssecpatchrot[i] = sqrt(1 - topo->ia_params[type][type].csecpatchrot[i] * topo->ia_params[type][type].csecpatchrot[i]);
			//fprintf(stdout, " %f %f", topo->ia_params[type][type].csecpatchrot[0], topo->ia_params[type][type].ssecpatchrot[0]);
			
			topo->ia_params[type][type].pangl[i+2] = param[8];
			topo->ia_params[type][type].panglsw[i+2] = param[9];
			topo->ia_params[type][type].pcangl[i+2] = cos(param[8]/2.0/180*PI);                 // C1
			topo->ia_params[type][type].pcanglsw[i+2] = cos((param[8]/2.0+param[9])/180*PI);    // C2
			//topo->ia_params[type][type].pcangl[i] = topo->ia_params[type][type].pcangl[i];
			//topo->ia_params[type][type].pcanglsw[i] = topo->ia_params[type][type].pcanglsw[i];
			topo->ia_params[type][type].pcoshalfi[i+2] = cos((param[8]/2.0+param[9])/2.0/180*PI);
			topo->ia_params[type][type].psinhalfi[i+2] = sqrt(1.0 - topo->ia_params[type][type].pcoshalfi[i+2] * topo->ia_params[type][type].pcoshalfi[i+2]);
		}
		fprintf(stdout, " %f  %f %f", param[7], topo->ia_params[type][type].pangl[2], topo->ia_params[type][type].panglsw[2]);
	}
	if(fields == 9){
		int i;
		for(i = 0; i < 2; i++){
			topo->ia_params[type][type].chiral_cos[i] = cos(param[10] / 360 * PI);
			topo->ia_params[type][type].chiral_sin[i] = sqrt(1 - topo->ia_params[type][type].chiral_cos[i] * topo->ia_params[type][type].chiral_cos[i]);
			fprintf(stdout, " %f ", param[9]);
		}
	}

	// Volume
	if (geotype_i < SP)
		topo->ia_params[type][type].volume = 4.0/3.0*PI*pow((topo->ia_params[type][type].sigma)/2.0,3.0) + PI/2.0*topo->ia_params[type][type].len[0]*pow((topo->ia_params[type][type].sigma)/2.0,2.0) ;
	else 
		topo->ia_params[type][type].volume = 4.0/3.0*PI*pow((topo->ia_params[type][type].sigma)/2.0,3.0);
	if ( topo->ia_params[type][type].rcutwca > topo->sqmaxcut )
		topo->sqmaxcut = topo->ia_params[type][type].rcutwca;
	if ( topo->ia_params[type][type].rcut > topo->sqmaxcut )
		topo->sqmaxcut = topo->ia_params[type][type].rcut;
	fprintf(stdout, " \n");
	DEBUG_INIT("Finished filltypes");
	return 1;
}

/************************************************
 * String Manipulation stuff for parsing files
 ************************************************/
/* return string that goes before comand character*/
void beforecommand(char *str,char *pline,char commandc)
{
	char *dummy;

	void trim(char *);

	strcpy(str,pline);
	if ((dummy = strchr (str,commandc)) != NULL) (*dummy) = 0;
	trim (str);
}

/* return string that goes after command character */
void aftercommand(char *str, char *pline,char commandc)
{
	char *dummy;
	int i;

	void trim(char *);

	strcpy(str,pline);
	if ((dummy = strchr (str,commandc)) != NULL) {
		i=0;
		while( (*dummy) != str[i]) {
			str[i] = ' ';
			i++;
		}
		str[i] = ' ';
	}
	trim (str);
}

/* reads a string from stream of max length n */
char *fgets2(char *line, int n, FILE *stream)
{
	char *c;

	if (fgets(line,n,stream)==NULL) {
		return NULL;
	}
	if ((c=strchr(line,'\n'))!=NULL)
		*c=0;

	return line;
}

/* remove comments */
void strip_comment (char *line)
{
	char *c;

	if (!line) return;
	/* search for a comment mark and replace it by a zero */
	if ((c = strchr(line,COMMENTSIGN)) != NULL) (*c) = 0;
}

/*test is there is still something left in string*/
int continuing(char *s)
{
	int sl;

	void rtrim (char *str);

	rtrim(s);
	sl = strlen(s);
	if ((sl > 0) && (s[sl-1] == CONTINUE)) {
		s[sl-1] = 0;
		return 1; /*true*/
	} else return 0; /*false*/
}

/*make strin uppercase*/
void upstring (char *str)
{
	int i;
	for (i=0; (i < (int)strlen(str)); i++) str[i] = toupper(str[i]);
}

/*trim string from left*/
void ltrim (char *str)
{
	char *tr;
	int c;

	if (!str) return;
	tr = strdup (str);
	c  = 0;
	while ((tr[c] == ' ') || (tr[c] == '\n') || (tr[c] == '\t')) c++;
	strcpy (str,tr+c);
	free (tr);
}

/*trim string from right*/
void rtrim (char *str)
{
	int nul;

	if (!str) return;
	nul = strlen(str)-1;
	while ((nul > 0) && ((str[nul] == ' ') || (str[nul] == '\t') || (str[nul] == '\n')) ) {
		str[nul] = '\0';
		nul--;
	}
}

/*trim strin from left and right*/
void trim (char *str)
{
	void ltrim (char *str);
	void rtrim (char *str);

	ltrim (str);
	rtrim (str);
}

/**
 *  Dumps a configuration to the supplied file handle.
 */
void draw(FILE *outfile, /*struct vector box, long npart,
			   struct particles *particle,*/ struct conf * conf, struct topo * topo)
{
	long i;

	double anint(double);

	//fprintf (outfile, "%15.8le %15.8le %15.8le\n", box.x, box.y, box.z);
	for (i = 0; i < topo->npart; i++) {
		fprintf (outfile, "%15.8le %15.8le %15.8le   %15.8le %15.8le %15.8le   %15.8le %15.8le %15.8le %d\n",
				conf->box.x * ((conf->particle[i].pos.x) - anint(conf->particle[i].pos.x)),
				conf->box.y * ((conf->particle[i].pos.y) - anint(conf->particle[i].pos.y)),
				conf->box.z * ((conf->particle[i].pos.z) - anint(conf->particle[i].pos.z)),
				conf->particle[i].dir.x, conf->particle[i].dir.y, conf->particle[i].dir.z,
				conf->particle[i].patchdir[0].x, conf->particle[i].patchdir[0].y, conf->particle[i].patchdir[0].z,
				conf->particle[i].switched);
	}
}
/*............................................................................*/


/****************************************************************************/
/* Pairlist stuf */
/****************************************************************************/
/** 
 * Initializes the pairlist and allocates memory
 */
void init_pairlist(struct topo * topo, struct sim * sim){
	printf("\nAllocating memory for pairlist...\n");
	sim->pairlist = xmalloc(sizeof(struct pairs) * topo->npart);
	// Highest guess: Every particle interacts with the others
	// TODO: Make it more sophisticated
	long i;
	for(i = 0; i < topo->npart; i++){
		sim->pairlist[i].pairs = malloc(sizeof(long) * topo->npart);
		sim->pairlist[i].num_pairs = 0;
	}
}
/*............................................................................*/

/** 
 * Cleans up: deallocates the memory for the pairlist 
 */
int dealloc_pairlist(struct topo * topo, struct sim * sim){
	long i;
	if(sim->pairlist != NULL){
		for(i = 0; i < topo->npart; i++){
			if(sim->pairlist[i].pairs != NULL){
				free(sim->pairlist[i].pairs);
			}
		}
		free(sim->pairlist);
	}
	return 0;
}
/*............................................................................*/

/**
 * Generates a pairlist with a very basic alogrithm 
 */
void gen_simple_pairlist(struct topo * topo, struct sim * sim, struct conf * conf){
	struct vector r_cm;
	double r_cm2;
	double max_dist;
	// Set the pairlist to zero
	//DEBUG_INIT("Gen Pairlist")
	long i, j;
	for(i = 0; i < topo->npart; i++){
		//DEBUG_INIT("%ld", i);
		sim->pairlist[i].num_pairs = 0;
	}
	long nj = topo->npart;
	long ni = nj - 1;
	for(i = 0; i < ni; i++){
		for(j = i + 1; j < nj; j++){
			r_cm.x = conf->particle[i].pos.x - conf->particle[j].pos.x;
			r_cm.y = conf->particle[i].pos.y - conf->particle[j].pos.y;
			r_cm.z = conf->particle[i].pos.z - conf->particle[j].pos.z;
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
			max_dist = AVER(sim->trans[conf->particle[i].type].mx, \
					sim->trans[conf->particle[j].type].mx);
			max_dist *= (1 + sim->pairlist_update) * 2;
			max_dist += topo->maxcut;
			max_dist *= max_dist; /* squared */

			if (r_cm2 <= max_dist){
				sim->pairlist[i].pairs[sim->pairlist[i].num_pairs++] = j;
				sim->pairlist[j].pairs[sim->pairlist[j].num_pairs++] = i;
			} 
		}
	}
	////Check for too many pairs
	//for(i = 0; i < topo->npart; i++){
	//    //if (sim->pairlist.list[i].num_pairs >= topo->npart)
	//    if (sim->pairlist[i].num_pairs >= topo->npart){
	//        fprintf(stderr, "ERROR: Too many pairs for particle %ld!!!\n", i);
	//        exit(1);
	//    }
	//}
}
/*.............................................................................*/

/**
 * Interface for the generation of the pairlist. Define other pairlist 
 * algorithms above.
 */
void gen_pairlist(struct topo * topo, struct sim * sim, struct conf * conf){
	gen_simple_pairlist(topo, sim, conf);
}
/*.............................................................................*/

/**
 * Print out the pairlist 
 */
void print_pairlist(FILE * stream, struct sim * sim, struct topo * topo){
	long i, j;
	for (i = 0; i < topo->npart; i++){
		fprintf(stream, "%ld (%ld):", i, sim->pairlist[i].num_pairs);
		for(j = 0; j < sim->pairlist[i].num_pairs; j++){
			fprintf(stream, " %ld", sim->pairlist[i].pairs[j]);
		}
		fprintf(stream, "\n");
	}
}
/*..........................................................................*/


/****************************************************************************/
/* Cluster statistics stuf */
/****************************************************************************/
/**
 * determines, wheter two particles are in the same cluster
 */
int same_cluster(struct topo * topo, struct conf * conf, 
		long fst, long snd, double (* intfce[MAXT][MAXT])(struct interacts *) ){

  
	/*if two particles are bonded they belong to the same cluster*/
	if ( ((topo->chainparam[conf->particle[fst].chaint]).bond1c >= 0) ||
	    ((topo->chainparam[conf->particle[fst].chaint]).bonddc >= 0) ){
		if ( (snd == topo->conlist[fst][1]) || (snd == topo->conlist[fst][0]) ) {
		  return TRUE;
		}
	}
	if ( ((topo->chainparam[conf->particle[snd].chaint]).bond1c >= 0) ||
	    ((topo->chainparam[conf->particle[snd].chaint]).bonddc >= 0) ){
		if ( (fst == topo->conlist[snd][1]) || (fst == topo->conlist[snd][0]) ) {
		  return TRUE;
		}
	}
	
	
	/*cluster is made of particles closer tna some distance*/
/*	struct vector image(struct vector r1, struct vector r2, struct vector box);
	struct vector r_cm = image(conf->particle[fst].pos,
			conf->particle[snd].pos,
			conf->box);
	double dist2 = DOT(r_cm, r_cm);
	* TODO: Make it much more efficient => define cluster_dist!!! *
	if(dist2 > topo->ia_params[conf->particle[fst].type][conf->particle[snd].type].sigma * topo->ia_params[conf->particle[fst].type][conf->particle[snd].type].sigma*4.0){
		return FALSE;
	}
	else {
		return TRUE;
	}*/

	/*cluster is made of attractively interacting particles*/
	double paire(long, long, double (* intfce[MAXT][MAXT])(struct interacts *), 
			struct topo * topo, struct conf * conf);

	
	if(paire(fst, snd, intfce, topo, conf) > -0.10 ){
		return FALSE;
	}
	else {
		return TRUE;
	}
}
/*............................................................................*/

/** 
 * generate the clusterlist 
 */
int gen_clusterlist(struct topo * topo, struct sim * sim, struct conf * conf, double (* intfce[MAXT][MAXT])(struct interacts *)){
	int change = TRUE; /* does it still change? */
	//long neighbour;
	long i, j, fst, snd, tmp, minnumber, maxnumber;
	
	int same_cluster(struct topo * topo, struct conf * conf, 
		long fst, long snd, double (* intfce[MAXT][MAXT])(struct interacts *));
	

	// Set clusterindex to the corresponding index
	for( i = 0; i < topo->npart; i++){
		sim->clusterlist[i] = i;
	}

	// Start determining the cluster
	while(change){
		change = FALSE;
		for(i = 0; i < topo->npart; i++){
			/*If nore pairlist go over all pairs*/
			maxnumber = topo->npart;
			minnumber = i ;
			if (sim->pairlist_update) {
			    maxnumber = sim->pairlist[i].num_pairs;
			    minnumber=0;
			}
			/* Go over pairs to see if they are in the cluster */
			for(j = minnumber; j < maxnumber; j++){
				fst = i;
				snd = j;
				if (sim->pairlist_update) {
				  snd = sim->pairlist[i].pairs[j];
				}
				/*do cluster analysis only for spherocylinders*/
				if ( (topo->ia_params[conf->particle[fst].type][conf->particle[snd].type].geotype[0] < SP) && \
					(topo->ia_params[conf->particle[fst].type][conf->particle[snd].type].geotype[1] < SP) ) {
				
					/* if they are close to each other */
					if(same_cluster(topo, conf, fst, snd, intfce)){
						if(fst > snd){
							tmp = snd;
							snd = fst;
							fst = tmp;
						}

						if(sim->clusterlist[fst] < sim->clusterlist[snd]){
							sim->clusterlist[snd] = sim->clusterlist[fst];
							change = TRUE;
							break;
							/* => will eventually start the i loop from new */
						}
						if(sim->clusterlist[snd] < sim->clusterlist[fst]){
							sim->clusterlist[fst] = sim->clusterlist[snd];
							change = TRUE;
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
/*............................................................................*/

/**
 * sort the clusterlist
 */
int sort_clusterlist(struct topo * topo, struct sim * sim){
	long cluster_indices[topo->npart];   /* holds the different cluster indices.
						(currently too much memory) */
	long num_cluster = 0;                /* number of clusters, temporary needed */
	long i, j;

	/* how many clusters are there? */
	long max_index = -1;
	for(i = 0; i < topo->npart; i++){
		if(max_index < sim->clusterlist[i]){
			max_index = sim->clusterlist[i];
			cluster_indices[num_cluster++] = max_index;
		}
	}

	/* free the memory from the old clusters */
	if(sim->clusters){
		for(i = 0; i < sim->num_cluster; i++){
			if(sim->clusters[i].particles){
				free(sim->clusters[i].particles);
			}
		}
		free(sim->clusters);
	}

	/* Allocate memory for the clusters */
	sim->clusters = xmalloc(sizeof(struct cluster) * num_cluster);
	for(i = 0; i < num_cluster; i++){
		/* allocate maximal space for all the clusters */
		sim->clusters[i].particles = xmalloc(sizeof(long) * topo->npart); 
		sim->clusters[i].npart = 0;
	}

	/* fill in the particles belonging to one cluster */
	for(i = 0; i < num_cluster; i++){
		for(j = 0; j < topo->npart; j++){
			if(sim->clusterlist[j] == cluster_indices[i]){
				sim->clusters[i].particles[sim->clusters[i].npart++] = j;
			}
		}
	}
	sim->num_cluster = num_cluster;

	/* Find the biggest size */
	sim->max_clust = 0;
	for(i = 0; i < num_cluster; i++){
		if(sim->clusters[i].npart > sim->max_clust){
			sim->max_clust = sim->clusters[i].npart;
		}
	}
	/* Set the statistics to zero */
	sim->clusterstat = xmalloc(sizeof(long) * sim->max_clust);
	for(i = 0; i < sim->max_clust; i++){
		sim->clusterstat[i] = 0;
	}
	/* Do the statistics */
	for(i = 0; i < num_cluster; i++){
		sim->clusterstat[sim->clusters[i].npart - 1]++;
	}

	return 0;
}
/*............................................................................*/


/**
 * calculate energies of clusters 
 * */
int calc_clusterenergies(struct topo * topo, struct sim * sim, struct conf * conf, 
			 double (* intfce[MAXT][MAXT])(struct interacts *)){
	long i,j,k;
	
	double paire(long, long, double (* intfce[MAXT][MAXT])(struct interacts *), 
			struct topo * topo, struct conf * conf);
	
	for(i = 0; i < sim->num_cluster; i++){
		sim->clustersenergy[i]=0.0;
		for(j = 0; j < sim->clusters[i].npart; j++){
			for(k = j+1; k < sim->clusters[i].npart; k++){
				sim->clustersenergy[i]+= paire(sim->clusters[i].particles[j], sim->clusters[i].particles[k], intfce, topo, conf);
			}
		}
	}
	return 0;
}
/*............................................................................*/

/**
 * print the clusterlist 
 * */
int print_clusterlist(FILE * stream, BOOL decor, 
		struct topo * topo, struct sim * sim, struct conf * conf){
	long i;
	if(decor){
		fprintf(stream, "\n"
				"-----------------------------------------------------\n"
				"  The Cluster List\n"
				"  (Index starts with 1)\n"
				"-----------------------------------------------------\n");
	}
	
	for(i = 0; i < topo->npart; i++){
		fprintf(stream,"%3ld %3ld %8.4lf %8.4lf %8.4lf", i + 1, 
				sim->clusterlist[i] + 1, 
				conf->particle[i].pos.x, 
				conf->particle[i].pos.y, 
				conf->particle[i].pos.z);
		fprintf(stream,"\n");
	}
	if(decor){
		fprintf(stream,"-----------------------------------------------------\n");
	}
	fflush(stream);
	return 0;
}
/*............................................................................*/

/**
 * print the clusters 
 * */
int print_clusters(FILE * stream, BOOL decor, struct sim * sim){
	long i, j;
	
			
	if(decor){
		fprintf(stream, "\n"
				"-----------------------------------------------------\n"
				"  The Clusters\n"
				"  (Index starts with 1)\n"
				"-----------------------------------------------------\n");
	}
	for(i = 0; i < sim->num_cluster; i++){
		fprintf(stream, "%3ld(%f):", i + 1,sim->clustersenergy[i]);
		for(j = 0; j < sim->clusters[i].npart; j++){
			fprintf(stream, "%5ld", sim->clusters[i].particles[j] + 1);
		}
		fprintf(stream, "\n");
	}
	if(decor){
		fprintf(stream,"---------------------------------------------------\n");
	}
	fflush(stream);
	return 0;
}
/*............................................................................*/

/**
 * print a statistics for the clusters
 */
int print_clusterstat(FILE * stream, BOOL decor, struct sim * sim){
	long i;
	if(decor){
		fprintf(stream, "\n"
				"-----------------------------------------------------\n"
				"   Cluster Distribution\n"
				"-----------------------------------------------------\n");
	}
	for(i = 0; i < sim->max_clust; i++){
		fprintf(stream, "%5ld\t%5ld\n", i + 1, sim->clusterstat[i]);
	}
	if(decor){
		fprintf(stream, "--------------------------------------------------\n");
	}
	fflush(stream);
	return 0;
}
/*............................................................................*/

/**
 * Alternative way of printing the cluster statistics: everything is on
 * one line. First monomers, then dimers etc.
 */
int print_clstat_oneline(FILE * stream, long sweep, struct sim * sim){
	long i;
	fprintf(stream, "%ld: ", sweep);
	for(i = 0; i < sim->max_clust; i++){
		fprintf(stream, "%5ld\t", sim->clusterstat[i]);
	}
	fprintf(stream, "\n");
	fflush(stream);
	return 0;
}

/**
 * write out all the cluster stat in files, if file name is given
 */
int write_cluster(FILE * cl_stat, FILE * cl, FILE * cl_list, BOOL decor, long sweep, 
	    struct sim * sim, struct topo * topo, struct conf * conf, double (* intfce[MAXT][MAXT])(struct interacts *)){
	
	int gen_clusterlist(struct topo * topo, struct sim * sim, struct conf * conf, double (* intfce[MAXT][MAXT])(struct interacts *));
	int sort_clusterlist(struct topo * topo, struct sim * sim);
	int print_clusters(FILE * stream, BOOL decor, struct sim * sim);
    int calc_clusterenergies(struct topo * topo, struct sim * sim, struct conf * conf,
                             double (* intfce[MAXT][MAXT])(struct interacts *));
	
	gen_clusterlist(topo, sim, conf, intfce);
	sort_clusterlist(topo, sim);
  calc_clusterenergies(topo, sim, conf, intfce);
	if(cl_stat){
		if(decor == FALSE){
			// if no decor, this means usually into a file. Hence print info
			// about number of line per frame
			fprintf(cl_stat, "Sweep: %ld | Maximal size: %ld\n", 
					sweep, sim->max_clust);
		}
		print_clusterstat(cl_stat, decor, sim);
		/*
		   print_clstat_oneline(cl_stat, sweep, sim);
		 */
	}
	if(cl){
		if(decor == FALSE){
			fprintf(cl, "Sweep: %ld | Number of clusters: %ld\n", 
					sweep, sim->num_cluster);
		}
		print_clusters(cl, decor, sim);
	}
	if(cl_list){
		if(decor == FALSE){
			fprintf(cl_list, "Sweep: %ld | Number of particles: %ld\n", 
					sweep, topo->npart);
		}
		print_clusterlist(cl, decor, topo, sim, conf);
	}
	return 0;
}
/*............................................................................*/


/****************************************************************************/
/* Wang-Landau stuf */
/****************************************************************************/
/*
   Initiate Wang-Landau calculation.
 */

int wlinit(struct wls *wl, char filename[30])
{
	long i,length,fields=0;
	double field[5];
	FILE *infile;
	char line[STRLEN];

	int wlend(struct wls *);
	void trim(char *);
	void strip_comment(char *);

	infile = fopen(filename, "r");
	if (infile == NULL) {
		fprintf (stderr, "\nERROR: Could not open %s file.\n\n",filename);
		return 1;
	}
	length=0;
	while (fgets2(line,STRLEN-2,infile) != NULL) {
		strip_comment (line);
		trim (line);
		/* if there is something left... */
		if ((int)strlen(line) > 0) {
			length++;
		}
	}
	length--; /*there is alpha at the first line*/
	(*wl).weights = malloc( sizeof(double) * length );
	(*wl).hist = malloc( sizeof(long) * length );
	(*wl).length[1] = 0;
	(*wl).dorder[1] = 0;

	fseek(infile,0,SEEK_SET);
	i=0;
	while (fgets2(line,STRLEN-2,infile) != NULL) {
		strip_comment (line);
		trim (line);
		/* if there is something left... */
		if ((int)strlen(line) > 0) {
			if (i == 0) {
				if (sscanf(line, "%le",&(*wl).alpha)!= 1) {
					fprintf (stderr, "ERROR: Could not read alpha at the begining.\n\n");
					wlend(wl);
					return 1;
				} else i++;
			} else {
				fields = sscanf(line, "%le %le %le %le",&field[0],&field[1],&field[2],&field[3]);
				if ( fields == 3 ) {
					if (i==1)
						(*wl).minorder[0] = field[0];
					(*wl).weights[i-1] = field[1];
					(*wl).hist[i-1] = field[2];
					(*wl).length[0]++;
					i++;
				} else if (fields == 4 ) {  
					if (i==1) {
						(*wl).minorder[0] = field[0];
						(*wl).minorder[1] = field[1];
					}
					if ( (*wl).minorder[1] == field[1] )
					      (*wl).length[0]++;
					(*wl).weights[i-1] = field[2];
					(*wl).hist[i-1] = field[3];
					i++;
				} else {
					fprintf (stderr, "ERROR: Could not read order parameter at line %ld.\n\n", i);
					wlend(wl);
					return 1;
				} 
			}
		}
	}
	
	if (fields == 4 ) { 
		(*wl).length[1] = length / (*wl).length[0];
		(*wl).dorder[1] = (field[1] - (*wl).minorder[1])/((*wl).length[1]-1);
	}
	(*wl).dorder[0] = (field[0] - (*wl).minorder[0])/((*wl).length[0]-1);
	if  ( ( (i-1) != (*wl).length[0] ) && (fields==3) )  {
		fprintf (stderr, "ERROR: In reading order parameters length %ld does not fit number of lines %ld.\n\n", (*wl).length[0],i-1);
		wlend(wl);
		return 1;
	}
	if  ( ( (i-1) != (*wl).length[0]*(*wl).length[1] ) && (fields==4) )  {
		fprintf (stderr, "ERROR: In reading order parameters lengths %ld %ld does not fit number of lines %ld.\n\n", (*wl).length[0],(*wl).length[1],i-1);
		wlend(wl);
		return 1;
	}
	/*DEBUG*/ 
	printf("Wang-Landau method init:\n");
	printf("alpha: %f\n",(*wl).alpha);
	/*int j=0;
	if ((*wl).length[1] == 0) {
		for (i=0; i<(*wl).length[0]; i++) {
			printf ("%15.8le %15.8le %ld \n",(*wl).minorder[0] + i * (*wl).dorder[0], (*wl).weights[i], (*wl).hist[i]);
		}
	} else {
		for (j=0; j<(*wl).length[1]; j++) {
			for (i=0; i<(*wl).length[0]; i++) {
				printf ("%15.8le %15.8le %15.8le %ld \n",(*wl).minorder[0] + i * (*wl).dorder[0], (*wl).minorder[1]+j*(*wl).dorder[1], (*wl).weights[i+(*wl).length[0]*j], (*wl).hist[i+(*wl).length[0]*j]);
			}
			printf (" \n");
		}
	}*/
	fclose(infile);
	fflush(stdout); 
	/**/
	return 0;
}
int wlwrite(struct wls *wl, char filename[30])
{
	long i,j;
	FILE *outfile;

	outfile = fopen(filename, "w");
	if (outfile == NULL) {
		fprintf (stderr, "\nERROR: Could not open %s file.\n\n",filename);
		return 1;
	}
	fprintf (outfile, "%15.8le \n",(*wl).alpha);
	if ((*wl).length[1] == 0) {
		for (i=0; i<(*wl).length[0]; i++) {
			fprintf (outfile, "%15.8le %15.8le %ld \n",(*wl).minorder[0] + i * (*wl).dorder[0], (*wl).weights[i], (*wl).hist[i]);
		}
	} else {
		for (j=0; j<(*wl).length[1]; j++) {
			for (i=0; i<(*wl).length[0]; i++) {
				fprintf (outfile, "%15.8le %15.8le %15.8le %ld \n",(*wl).minorder[0] + i * (*wl).dorder[0], (*wl).minorder[1]+j*(*wl).dorder[1], (*wl).weights[i+(*wl).length[0]*j], (*wl).hist[i+(*wl).length[0]*j]);
			}
			fprintf (outfile, " \n");
		}
	}
	
	fflush(outfile);
	fclose(outfile);

	return 0;
}
int wlend(struct wls *wl)
{
	free((*wl).weights);
	free((*wl).hist);

	return 0;
}

void wlreject(struct sim *sim, long oldlength)
{
	int mesh_cpy(struct meshs *, struct meshs *);
	int longarray_cpy (long **target, long **source,long,long);

	if ( sim->wlm[0] > 0 ) {
		sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]] -= sim->wl.alpha;
		sim->wl.hist[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]]++;
		if ( (sim->wlm[0] == 2) || (sim->wlm[1] == 2) )   
			mesh_cpy(&sim->wl.mesh,&sim->wl.origmesh);
		if ( (sim->wlm[0] == 5) || (sim->wlm[1] == 5)||(sim->wlm[0] == 6) || (sim->wlm[1] == 6) ) {
			longarray_cpy(&sim->wl.radiushole,&sim->wl.radiusholeold,sim->wl.radiusholemax,oldlength);
			sim->wl.radiusholemax = oldlength;
		}
		sim->wl.partincontact = sim->wl.partincontactold;
		
	}
}

void wlaccept(int wlm,struct wls *wl)
{
	int i;
	
	if ( wlm > 0 ) {
		for (i=0;i<2;i++)
			(*wl).currorder[i] = (*wl).neworder[i];
		(*wl).weights[ (*wl).currorder[0] + (*wl).currorder[1] * (*wl).length[0]] -= (*wl).alpha;
		(*wl).hist[ (*wl).currorder[0] + (*wl).currorder[1] * (*wl).length[0]]++;
	}
}

/*..............................................................................*/
/*........................NAMETIC ORDER.........................................*/
/*..............................................................................*/

/*
   Calculates the instantaneous value of the nematic order parameter for the
   specified configuration.  The nematic director is determined by diagonalisation
   of the tensor order parameter Q (see Allen & Tildesley p305).  The order
   parameter is the corresponding eigenvalue.  However, it is equivalent to take
   minus two times the middle eigenvalue (see Eppenga & Frenkel, Mol Phys vol.
   52, p.1303-1334 [1984]), and this is more reliable for comparing the isotropic
   phase.  This is the approach taken in this implementation.

   Routines from Numerical Recipes are used to perform the diagonalisation.  Note
   that these routines expect an n*n matrix to be stored in elements [1...n][1...n],
   rather than [0...n-1][0...n-1], so the arrays must be declared with one more
   element in each dimension.
 */

double nematic(long npart, struct particles *p)
{
	double q[4][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	double d[4], e[4];
	long i;

	void tred2(double [4][4], double [4], double [4]);
	void tqli(double [4], double [4]);

	for (i=0; i<npart; i++) {
		q[1][1] += p[i].dir.x * p[i].dir.x;
		q[1][2] += p[i].dir.x * p[i].dir.y;
		q[1][3] += p[i].dir.x * p[i].dir.z;
		q[2][1] += p[i].dir.y * p[i].dir.x;
		q[2][2] += p[i].dir.y * p[i].dir.y;
		q[2][3] += p[i].dir.y * p[i].dir.z;
		q[3][1] += p[i].dir.z * p[i].dir.x;
		q[3][2] += p[i].dir.z * p[i].dir.y;
		q[3][3] += p[i].dir.z * p[i].dir.z;
	}
	q[1][1] = (q[1][1] * 3.0 / npart - 1.0) / 2.0;
	q[1][2] = (q[1][2] * 3.0 / npart      ) / 2.0;
	q[1][3] = (q[1][3] * 3.0 / npart      ) / 2.0;
	q[2][1] = (q[2][1] * 3.0 / npart      ) / 2.0;
	q[2][2] = (q[2][2] * 3.0 / npart - 1.0) / 2.0;
	q[2][3] = (q[2][3] * 3.0 / npart      ) / 2.0;
	q[3][1] = (q[3][1] * 3.0 / npart      ) / 2.0;
	q[3][2] = (q[3][2] * 3.0 / npart      ) / 2.0;
	q[3][3] = (q[3][3] * 3.0 / npart - 1.0) / 2.0;

	tred2 (q, d, e);
	tqli (d, e);

	/* Sort eigenvalues */
	if (d[1] > d[2]) { d[0]=d[1]; d[1]=d[2]; d[2]=d[0]; }
	if (d[2] > d[3]) { d[0]=d[2]; d[2]=d[3]; d[3]=d[0]; }
	if (d[1] > d[2]) { d[0]=d[1]; d[1]=d[2]; d[2]=d[0]; }

	return -2.0*d[2];
}

/*..............................................................................*/

/*
   Returns the coefficient of the Fourier series term with period boxlength/n
   in the z direction.  The coefficients of the sine and cosine terms are added
   in quadrature and returned, making the result independent of phase shifts in
   the z direction.  A significantly non-zero value indicates layering of the
   particles in the z direction with periodicity boxlength/n.
 */

double smectic(long npart, struct particles *p, long n)
{
	double a, b;
	double omega = 8.0*n*atan(1.0);
	long i;

	a = b = 0.0;

	for (i=0; i<npart; i++) {
		a += cos(omega * p[i].pos.z);
		b += sin(omega * p[i].pos.z);
	}
	a /= (double)npart;
	b /= (double)npart;

	return sqrt(a*a + b*b);
}


/*..............................................................................*/
/*........................Z ORDER PARAMETER.....................................*/
long z_order(struct wls *wl, struct conf * conf,int wli)
{
	//    printf("%f %ld\n",particle[0].pos.z * box.z,lround(particle[0].pos.z * box.z / wl.dorder[wli] - wl.minorder[wli]));
	/* Because older C compilators do not know lround we can use ceil as well
	   return lround(particle[0].pos.z * box.z / wl.dorder[wli] - wl.minorder[wli]);*/
	/*
	printf("%f ",conf->particle[0].pos.z );
	printf("%f ",conf->syscm.z);
	printf("%f ",conf->box.z);
	printf("%f ", wl->minorder[wli]);
	printf("%f \n", wl->dorder[wli] );*/
	
	return (long) ceil( ((conf->particle[0].pos.z - conf->syscm.z) * conf->box.z- wl->minorder[wli]) / wl->dorder[wli]  );
}
/*..............................................................................*/
/*........................2 particles distance.....................................*/
long twopartdist(struct wls *wl, struct conf * conf, int wli)
{
	struct vector r_cm;
	
	r_cm.x = conf->particle[0].pos.x - conf->particle[1].pos.x;
	r_cm.y = conf->particle[0].pos.y - conf->particle[1].pos.y;
	r_cm.z = conf->particle[0].pos.z - conf->particle[1].pos.z;
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
	
	return (long) ceil( (( sqrt(r_cm.x*r_cm.x + r_cm.y*r_cm.y) ) - wl->minorder[wli]) / wl->dorder[wli]  );
}

/*..............................................................................*/
/*........................alignment ORDER PARAMETER.....................................*/
double alignment_order(struct conf * conf, struct topo * topo)
{
	double sumdot=0;
	long i,j;
	struct vector r_cm;
	struct vector image(struct vector, struct vector, struct vector);
	
	for (i = 0; i < topo->npart - 1; i++) {
		for (j = i + 1; j < topo->npart; j++) {
			r_cm = image(conf->particle[i].pos, conf->particle[j].pos, conf->box);
			if ( DOT(r_cm,r_cm) < 1.5*1.5 ) {
				sumdot+= DOT(conf->particle[i].dir,conf->particle[j].dir);
			}
		}
	}
	
	return sumdot;
}


/*..............................................................................*/
/*........................HOLE IN MESH-MEMBRANE ORDER PARAM.....................*/


/* return change in order parameter when one particle moves*/
long meshorder_moveone(struct vector oldpos, struct vector newpos, struct meshs *mesh,
		long npart, long target, struct conf * conf, struct sim * sim, int wli)
{
	int change;
	int nx,ny,ox,oy; /* position in mesh */
	double resid;

	void mesh_fill(struct meshs *, long , struct particles *, struct sim * sim);
	int mesh_findholes(struct meshs *);
	int mesh_addpart(double, double, int **, int [2]);
	int mesh_removepart(double, double, int **, int [2]);

	if ( conf->particle[target].type != sim->wl.wlmtype ) 
		return sim->wl.currorder[wli];

	nx = (int) (INBOX(newpos.x,resid) * (*mesh).dim[0]);
	ny = (int) (INBOX(newpos.y,resid) * (*mesh).dim[1]);
	ox = (int) (INBOX(oldpos.x,resid) * (*mesh).dim[0]);
	oy = (int) (INBOX(oldpos.y,resid) * (*mesh).dim[1]);
	if ( (nx == ox) && (ny == oy) ) return sim->wl.currorder[wli]; /* particle stayed in the same mesh bin*/

	change = mesh_addpart(newpos.x,newpos.y,&(*mesh).data,(*mesh).dim);
	if (change) {
		change = mesh_removepart(oldpos.x,oldpos.y,&(*mesh).data,(*mesh).dim);
	}
	if ( !change ) {
		/* fill the mesh with particles*/
		mesh_fill(mesh,npart,conf->particle, sim);
		return (long) (mesh_findholes(mesh) - sim->wl.minorder[wli]);
	}
	return sim->wl.currorder[wli];
}

/* return change in order parameter when chain moves*/
long meshorder_movechain(long chain[MAXN], struct meshs *mesh, 
		long npart, struct conf * conf, struct sim * sim,struct particles chorig[MAXCHL], int wli)
{
	long i,current;
	int change;

	void mesh_fill(struct meshs *, long , struct particles *, struct sim * sim);
	int mesh_findholes(struct meshs *);
	int mesh_addpart(double, double, int **, int [2]);
	int mesh_removepart(double, double, int **, int [2]);


	change= 1;
	i = 0;
	current = chain[0];
	while ( (current >=0 ) && (change) ) {
		if ( conf->particle[current].type == sim->wl.wlmtype )
			change = mesh_addpart(conf->particle[current].pos.x, conf->particle[current].pos.y, &(*mesh).data, (*mesh).dim);
		i++;
		current = chain[i];
	}
	i = 0;
	current = chain[0];
	while ( (current >=0 ) && (change) ) {
		if ( conf->particle[current].type == sim->wl.wlmtype )
			change = mesh_removepart(chorig[i].pos.x, chorig[i].pos.y, &(*mesh).data, (*mesh).dim);
		i++;
		current = chain[i];
	}

	if ( !change ) {
		/* fill the mesh with particles*/
		mesh_fill(mesh,npart,conf->particle, sim);
		return (long) (mesh_findholes(mesh) - sim->wl.minorder[wli]);
	}
	return sim->wl.currorder[wli];
}

/* filling the mesh */
void mesh_fill(struct meshs *mesh, long npart, struct particles *particle, struct sim * sim)
{
	long i;

	int mesh_addpart(double posx, double posy, int **mesh, int dim[2]);

	for ( i=0; i<((*mesh).dim[0] * (*mesh).dim[1]); i++) {
		(*mesh).data[i] = 0;
	}

	for (i=0; i<npart; i++) {
		/*calculate position of particle on mesh and add it to all where it belongs */
		if (particle[i].type == sim->wl.wlmtype)
			mesh_addpart(particle[i].pos.x,particle[i].pos.y, &(*mesh).data, (*mesh).dim);
	}

}

/* add particle on coordinates posx posy to mesh return 0 if it was placed on empty spot*/
int mesh_addpart(double posx, double posy, int **mesh, int dim[2])
{
	int i, square[9], onhole;
	double resid;

	void mesh_square(int , int , int [2], int (*)[9]);

	onhole = 1;
	mesh_square( (int) (INBOX(posx,resid) * dim[0]), (int) (INBOX(posy,resid) * dim[1]) , dim, &square);
	for(i=0;i<9;i++) {
		if ( (square[i] >= dim[0]*dim[1])||(square[i] <0) ) {
			printf ("Error: trying to write to %d\n",square[i]);
			printf ("%d %d and  %d\n", (int) (INBOX(posx,resid) * dim[0]), (int) (INBOX(posy,resid) * dim[1]),i );
			fflush(stdout);
		}
		if ( ((*mesh)[ square[i] ]) >= 0 ) onhole =  0;
		(*mesh)[ square[i] ]--;
	}
	return onhole;
}

/* remove particle on coordinates posx posy from mesh and return 0 if there is a empty spot now*/
int mesh_removepart(double posx, double posy, int **mesh, int dim[2])
{
	int i, square[9];
	double resid;

	void mesh_square(int , int , int [2], int (*)[9]);

	mesh_square((int) (INBOX(posx,resid) * dim[0]), (int) (INBOX(posy,resid) * dim[1]) , dim, &square);
	for(i=0;i<9;i++) {
		//DEBUG	    if (square[i] >= dim[0]*dim[1]) printf ("Error: trying to write to %d\n",square[i]);
		(*mesh)[ square[i] ]++;
		if ( ((*mesh)[ square[i] ]) == 0 ) return 0;
	}

	return 1;
}

void mesh_square(int x, int y, int dim[2], int (*square)[9])
{
	int a,b;

	b=y;
	(*square)[0] = x + dim[0]*b;
	a = x-1; 
	if ( a<0 ) a = dim[0]-1;
	(*square)[1] = a + dim[0]*b;
	a = x+1; 
	if ( a==dim[0] ) a = 0;
	(*square)[2] = a + dim[0]*b;

	b = y-1; 
	if ( b<0 ) b = dim[1]-1;
	(*square)[3] = x + dim[0]*b;
	a = x-1; 
	if ( a<0 ) a = dim[0]-1;
	(*square)[4] = a + dim[0]*b;
	a = x+1; 
	if ( a==dim[0] ) a = 0;
	(*square)[5] = a + dim[0]*b;

	b = y+1; 
	if ( b==dim[1] ) b = 0;
	(*square)[6] = x + dim[0]*b;
	a = x-1; 
	if ( a<0 ) a = dim[0]-1;
	(*square)[7] = a + dim[0]*b;
	a = x+1; 
	if ( a==dim[0] ) a = 0;
	(*square)[8] = a + dim[0]*b;

}

void mesh_neighbors(int pos, int dim[2], int neighbors[4])
{
	int x,y,a;

	x = pos % dim[0];
	y = pos / dim[0];

	a = x-1; 
	if ( a<0 ) a = dim[0]-1;
	neighbors[0] = a + dim[0]*y;
	a = x+1; 
	if ( a==dim[0] ) a = 0;
	neighbors[1] = a + dim[0]*y;

	a = y-1; 
	if ( a<0 ) a = dim[1]-1;
	neighbors[2] = x + dim[0]*a;
	a = y+1; 
	if ( a==dim[1] ) a = 0;
	neighbors[3] = x + dim[0]*a;

}

/* returns the number of holes and a list of mesh points belonging to each of them */
int mesh_findholes(struct meshs *mesh)
{
	int i,j, k, n, size, li, maxsize;
	int neighbors[4];

	void mesh_neighbors(int, int [2], int [4]);

	n=0;
	maxsize = 0;
	for (i=0;i<((*mesh).dim[0] * (*mesh).dim[1]);i++) {
		(*mesh).tmp[i] = 0;
		if ( (*mesh).data[i] > 0 ) (*mesh).data[i] = 0;
	}
	i=0;
	// go through all mesh points
	while ( i < ((*mesh).dim[0] * (*mesh).dim[1]) ) {
		// test if mesh point is occupied 
		if ( (*mesh).data[i] != 0 ) { i++; }
		else {
			// mesh point is free, create a new cluster
			n++;
			(*mesh).data[i] = n;
			// start new cluster, put mesh point as first element, and set list pointer on first element
			//DEBUG      if (n >= mesh.dim[0]*mesh.dim[1]) printf ("Error: trying to write to sizes position %d\n",n);
			size = 1;
			(*mesh).tmp[0] = i;
			li = 0;
			// go through all elements of the cluster
			while ( li < size ) {
				//go through all neighbors
				j =  (*mesh).tmp[li];
				mesh_neighbors(j, (*mesh).dim, neighbors);
				for ( k=0; k<4; k++ ) {
					// test if status is free and append it to the cluster
					if ( (*mesh).data[ neighbors[k] ] == 0 ) {
						(*mesh).data[ neighbors[k] ] = n;
						// append mesh point as element in the list
						(*mesh).tmp[size] = neighbors[k];
						size++;
					}
					if (  (*mesh).data[ neighbors[k] ] > 0 &&  (*mesh).data[ neighbors[k] ]<n ) {
						fprintf(stderr,"Error: Mesh cluster out of range, propably going infinite through pbc."); 
						fflush(stderr);
					}
				}
				li++;
			}
			if (size > maxsize) maxsize = size;
		}

	}

	return maxsize;
}

int mesh_init(struct meshs *mesh, double meshsize, long npart, struct conf * conf, struct sim * sim)
{
	//  int i;
	int maxsize,length;

	void mesh_fill(struct meshs *, long , struct particles *, struct sim * sim);
	int mesh_findholes(struct meshs *);

	(*mesh).dim[0] = (int)(conf->box.x/meshsize);
	(*mesh).dim[1] = (int)(conf->box.y/meshsize);
	if ( (*mesh).data != NULL ) free((*mesh).data);
	if ( (*mesh).tmp != NULL ) free((*mesh).tmp);
	length = (*mesh).dim[0] * (*mesh).dim[1];
	(*mesh).data = malloc( sizeof(int)* (length));
	(*mesh).tmp = malloc( sizeof(int)* (length+1));

	/* fill the mesh with particles*/
	mesh_fill(mesh, npart,conf->particle, sim);
	/* perfrom hole cluster algorithm */
	maxsize = mesh_findholes(mesh);

	/*DEBUG    printf("maxsize: %d\n",maxsize);
	  printf("mesh:\n");
	  for (i=0;i<mesh.dim[0]*mesh.dim[1];i++) {
	  printf("%d ",mesh.data[i]);
	  if ( ((i+1) % mesh.dim[0]) == 0) printf("\n");
	  }*/
	return maxsize;
}

void mesh_print (struct meshs *mesh)
{
	int i;

	int mesh_findholes(struct meshs *);

	printf("mesh:\n");
	for (i=0;i<(*mesh).dim[0] * (*mesh).dim[1];i++) {
		printf("%d ",(*mesh).data[i]);
		if ( ((i+1) % (*mesh).dim[0]) == 0) printf("\n");
	}
	printf("hole %d:\n", mesh_findholes(mesh) );
	printf("\n");
}

int mesh_cpy (struct meshs *target, struct meshs *source)
{
	if ( (*target).data != NULL) {
		if ( ((*target).dim[0] == (*source).dim[0]) && ((*target).dim[1] == (*source).dim[1]) ) {
			memcpy((*target).data,(*source).data, sizeof(int)* ((*target).dim[0] * (*target).dim[1])  );
			return 0;
		} else {
			free ((*target).data);
			if ( (*source).dim[0] * (*source).dim[1] > (*target).dim[0] * (*target).dim[1] ) {
				if ((*target).tmp != NULL ) free ((*target).tmp);
				(*target).tmp = malloc( sizeof(int)* ((*source).dim[0] * (*source).dim[1] + 1));
			}
		}
	}
	(*target).dim[0] = (*source).dim[0];
	(*target).dim[1] = (*source).dim[1];
	(*target).data = malloc( sizeof(int)* ((*target).dim[0] * (*target).dim[1]));
	if ((*target).tmp == NULL ) (*target).tmp = malloc( sizeof(int)* ((*source).dim[0] * (*source).dim[1] + 1));
	memcpy((*target).data,(*source).data, sizeof(int)* ((*target).dim[0] * (*target).dim[1])  );
	return 0;
}

int mesh_end(struct meshs *mesh)
{
	/* free allocated memory */
	if ( (*mesh).data!= NULL ) free((*mesh).data);
	if ( (*mesh).tmp!= NULL ) free((*mesh).tmp);
	return 0;
}



/*..............................................................................*/
/*........................RADIUS HOLE IN CENTER MEMBRANE ORDER PARAM............*/

/*return current bin of free radius*/
long radiushole_order(struct sim * sim)
{
	long i;
	
	for (i=0;i<sim->wl.radiusholemax-3;i++){
	    if ((sim->wl.radiushole[i] >0 ) && (sim->wl.radiushole[i+1] >0 ) && (sim->wl.radiushole[i+2] >0 ) && (sim->wl.radiushole[i+3] >0 ))
	      return i-1;
	}
	return -100;
}

/*return order of given radius */
long radiushole_position(double radius, struct sim * sim, int wli)
{
  return (long) ceil( ( radius - sim->wl.minorder[wli]) / sim->wl.dorder[wli]  );
}

/* return change in order parameter when one particle moves*/
long radiusholeorder_moveone(struct vector *oldpos, struct conf *conf, struct sim * sim, long target,int wli, struct vector *position)
{
	long nr,or; /* position in radiushole */
	double rx,ry,z;
	BOOL oz,nz;

	long radiushole_position(double radius, struct sim * sim,int);
	long radiushole_order(struct sim *sim);
	double anint(double);
	void radiushole_print (long *radiushole, long length);

	if ( conf->particle[target].type != sim->wl.wlmtype ) 
		return sim->wl.currorder[wli];
	
	z=conf->particle[target].pos.z - position->z; /*if above position*/
	if (z-anint(z) < 0) nz = FALSE;
	else nz=TRUE;
	z=oldpos->z - position->z;  /*if above position*/
	if (z-anint(z) < 0) oz = FALSE;
	else oz=TRUE;
	if ( !(nz) && !(oz) ) 
		return sim->wl.currorder[wli];
	
	rx = conf->box.x * (conf->particle[target].pos.x - anint(conf->particle[target].pos.x));
	ry = conf->box.y * (conf->particle[target].pos.y - anint(conf->particle[target].pos.y));
	nr = radiushole_position(sqrt(rx*rx+ry*ry),sim,wli);
	
	if (nr < 0)
	    return -100;
	/*particle move over radius bins*/
	if (nz) {
		sim->wl.radiushole[nr]++;
	}
	if (oz) {
		rx = conf->box.x * (oldpos->x - anint(oldpos->x));
		ry = conf->box.y * (oldpos->y - anint(oldpos->y));
		or = radiushole_position(sqrt(rx*rx+ry*ry),sim,wli);
		sim->wl.radiushole[or]--;  
		if ( sim->wl.radiushole[or] < 0 ) {
			printf ("Error(single particle move): trying to make number of beads in radiuspore smaller than 0 at position %ld\n",or);
			radiushole_print(sim->wl.radiushole,sim->wl.radiusholemax);
			fflush(stdout); 
		}
		if (sim->wl.radiushole[or] ==0)
			return radiushole_order(sim);
	}
	
	if ( (nz) && (sim->wl.radiushole[nr] ==1) )  {
		return radiushole_order(sim);
	} 
	
	return sim->wl.currorder[wli];
}

/* return change in order parameter when chain moves*/
long radiusholeorder_movechain(long chain[MAXN], struct conf * conf, struct sim * sim,struct particles chorig[MAXCHL],int wli, struct vector *position)
{
	long i,current,nr;
	double rx,ry,z;
	BOOL change=FALSE;

	long radiushole_position(double radius, struct sim * sim,int);
	long radiushole_order(struct sim *sim);
	double anint(double);
	void radiushole_print (long *radiushole, long length);

	i = 0;
	rx=0;
	current = chain[0];
	while (current >=0 ) {
		if ( conf->particle[current].type == sim->wl.wlmtype ) {
			z=conf->particle[current].pos.z - position->z; /*if above system CM*/
			if (z-anint(z) > 0) {
				rx = conf->box.x * (conf->particle[current].pos.x - anint(conf->particle[current].pos.x));
				ry = conf->box.y * (conf->particle[current].pos.y - anint(conf->particle[current].pos.y));
				nr = radiushole_position(sqrt(rx*rx+ry*ry),sim,wli);
				if (nr < 0)
				    return -100;
				sim->wl.radiushole[nr]++;
				if ( sim->wl.radiushole[nr] == 1 ) change = TRUE;
			}
		}
		i++;
		current = chain[i];
	}
	i = 0;
	current = chain[0];
	while (current >=0 )  {
		if ( conf->particle[current].type == sim->wl.wlmtype ) {
			z=chorig[i].pos.z - position->z; /*if above system CM*/
			if (z-anint(z) > 0) {
				rx = conf->box.x * (chorig[i].pos.x - anint(chorig[i].pos.x));
				ry = conf->box.y * (chorig[i].pos.y - anint(chorig[i].pos.y));
				nr = radiushole_position(sqrt(rx*rx+ry*ry),sim,wli);
				sim->wl.radiushole[nr]--;
				if ( sim->wl.radiushole[nr] < 0 ) {
					printf ("Error (chainmove): trying to make number of beads in radiuspore smaller than 0 at position %ld\n",nr);
					radiushole_print(sim->wl.radiushole,sim->wl.radiusholemax);
					fflush(stdout); 
				}
				if ( sim->wl.radiushole[nr] == 0 ) change = TRUE;
			}
		}
		i++;
		current = chain[i];
	}

	if ( change ) {
		return radiushole_order(sim);
	}
	return sim->wl.currorder[wli];
}

/* filling the radiushole above vec*/
long radiushole_all(struct topo *topo, struct conf *conf, struct sim * sim,int wli, struct vector *position)
{
	long i,nr,radiusholemax;
	double rx,ry,z;
	
	long radiushole_position(double radius, struct sim * sim,int);
	long radiushole_order(struct sim *sim);
	double anint(double);
	
	radiusholemax = radiushole_position(sqrt(conf->box.x*conf->box.x+conf->box.y*conf->box.y),sim,wli);
	if ( radiusholemax > sim->wl.radiusholemax ) {
		if (sim->wl.radiushole != NULL) 
			free(sim->wl.radiushole); 
		sim->wl.radiushole = malloc( sizeof(long)* (radiusholemax));
		sim->wl.radiusholemax = radiusholemax;
	}
		
	for (i=0;i<radiusholemax;i++) {
		sim->wl.radiushole[i] = 0;
	}

	for (i=0; i< topo->npart; i++) {
		/*calculate position of particle from z axis, and add it in array */
		if ( conf->particle[i].type == sim->wl.wlmtype ) {
			z=conf->particle[i].pos.z - (*position).z; /*if above position*/
			if (z-anint(z) > 0) {
				rx = conf->box.x * (conf->particle[i].pos.x - anint(conf->particle[i].pos.x));
				ry = conf->box.y * (conf->particle[i].pos.y - anint(conf->particle[i].pos.y));
				nr = radiushole_position(sqrt(rx*rx+ry*ry),sim,wli);
				if (nr < 0)
				    return -100;
				sim->wl.radiushole[nr]++;
			}
		}
	}

	return radiushole_order(sim);
}


void radiushole_print (long *radiushole, long length)
{
	long i;

	printf("radiushole:\n");
	for (i=0;i<length;i++) {
		printf("%ld ",radiushole[i]);
	}
	printf("\n");
}

int longarray_cpy (long **target, long **source, long targetlength, long sourcelength)
{
		
	/*if ( (*target) != NULL) {
		if ( targetlength == sourcelength ) {
			memcpy((*target),(*source), sizeof(long)*(sourcelength));
			return 0;
		} else {
			free(*target);
		}
	}*/
	if ( (*target) != NULL) 
		(*target) = (long*) realloc((*target), sizeof(long)*(sourcelength));
	else
		(*target) = malloc( sizeof(long)*(sourcelength));
	memcpy((*target),(*source), sizeof(long)*(sourcelength));
	
	return 0;
}

/*..............................................................................*/
/*  ............................... particles in contact .....................   */


/*return order for particles in contact */
long contparticles_order(struct sim * sim, int wli)
{
	return (long) ceil( ( sim->wl.partincontact - sim->wl.minorder[wli]) / sim->wl.dorder[wli]  );
}

/*returns if particle is in contact*/
BOOL particleinncontact (struct vector *vec, struct conf *conf)
{
	double x,y,z;
	
	double anint(double);

	x = vec->x - conf->particle[0].pos.x;
	y = vec->y - conf->particle[0].pos.y;
	z = vec->z - conf->particle[0].pos.z;
	x = conf->box.x * (x - anint(x));
	y = conf->box.y * (y - anint(y));
	z = conf->box.z * (z - anint(z));
	if ( x*x + y*y + z*z < WL_CONTACTS) {
		return TRUE;
	}
	else {
		return FALSE;
	}
	
}

/* return change in number of particles in contact when one particle moves*/
long contparticles_moveone(struct vector *oldpos, struct conf *conf, struct sim * sim, long target,int wli)
{
	long contparticles_order(struct sim * sim, int wli);
	BOOL particleinncontact (struct vector *vec, struct conf *conf);
	
	if ( conf->particle[target].type != sim->wl.wlmtype ) 
		return sim->wl.currorder[wli];
	
	if ( particleinncontact (&(conf->particle[target].pos),conf) )
		sim->wl.partincontact++;
	if ( particleinncontact (oldpos,conf) )
		sim->wl.partincontact--;
 
	return contparticles_order(sim,wli);
}

/* return change in order parameter when chain moves*/
long contparticles_movechain(long chain[MAXN], struct conf * conf, struct sim * sim,struct particles chorig[MAXCHL],int wli)
{
	long i,current;
	
	long contparticles_order(struct sim * sim, int wli);
	BOOL particleinncontact (struct vector *vec, struct conf *conf);
	
	i = 0;
	current = chain[0];
	while (current >=0 ) {
		if ( conf->particle[current].type == sim->wl.wlmtype ) {
			if ( particleinncontact (&(conf->particle[current].pos),conf) )
				sim->wl.partincontact++;
		}
		i++;
		current = chain[i];
	}
	i = 0;
	current = chain[0];
	while (current >=0 )  {
		if ( conf->particle[current].type == sim->wl.wlmtype ) {
			if ( particleinncontact (&(chorig[i].pos),conf) )
				sim->wl.partincontact--;
		}
		i++;
		current = chain[i];
	}
	
	return contparticles_order(sim,wli);
}

/* filling all particles in the contact */
long contparticles_all(struct topo *topo, struct conf *conf, struct sim * sim,int wli)
{
	long i;
	
	long contparticles_order(struct sim * sim, int wli);
	BOOL particleinncontact (struct vector *vec, struct conf *conf);
	
	
	sim->wl.partincontact = 0;
	
	for (i=1; i< topo->npart; i++) {
		/*calculate position of particle and add it if in contact */
		if ( conf->particle[i].type == sim->wl.wlmtype ) {
			if ( particleinncontact (&(conf->particle[i].pos),conf) )
				sim->wl.partincontact++;
		}
	}
	
	return contparticles_order(sim,wli);
}


/*..............................................................................*/
/*........................GEOMETRIC STUFF.......................................*/
/*..............................................................................*/


/*..............................................................................*/

/*
   Find closest distance between line segments and return its vector
   gets orientations and lengths of line segments and the vector connecting
   their center os masses (from vec1 to vec2)
 */
// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

struct vector mindist_segments(struct vector dir1, double halfl1,
		struct vector dir2, double halfl2, struct vector r_cm)
{
	struct vector u,v,w,vec;
	double a,b,c,d,e,D,sc,sN,sD,tc,tN,tD;

	struct vector vec_scale(struct vector, double);

	u = vec_scale(dir1,2.0*halfl1); //S1.P1 - S1.P0;
	v = vec_scale(dir2,2.0*halfl2); //S2.P1 - S2.P0;
	w.x = dir2.x*halfl2 - dir1.x*halfl1 - r_cm.x;
	w.y = dir2.y*halfl2 - dir1.y*halfl1 - r_cm.y;
	w.z = dir2.z*halfl2 - dir1.z*halfl1 - r_cm.z; //S1.P0 - S2.P0;
	a = DOT(u,u);        // always >= 0
	b = DOT(u,v);
	c = DOT(v,v);        // always >= 0
	d = DOT(u,w);
	e = DOT(v,w);
	D = a*c - b*b;       // always >= 0
	sc = D;
	sN = D;
	sD = D;      // sc = sN / sD, default sD = D >= 0
	tc = D;
	tN = D;
	tD = D;      // tc = tN / tD, default tD = D >= 0

	// compute the line parameters of the two closest points
	if (D < 0.00000001) { // the lines are almost parallel
		sN = 0.0;        // force using point P0 on segment S1
		sD = 1.0;        // to prevent possible division by 0.0 later
		tN = e;
		tD = c;
	}
	else {                // get the closest points on the infinite lines
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
			sN = 0.0;
			tN = e;
			tD = c;
		}
		else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}

	if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
		tN = 0.0;
		// recompute sc for this edge
		if (-d < 0.0)
			sN = 0.0;
		else if (-d > a)
			sN = sD;
		else {
			sN = -d;
			sD = a;
		}
	}
	else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
		tN = tD;
		// recompute sc for this edge
		if ((-d + b) < 0.0)
			sN = 0;
		else if ((-d + b) > a)
			sN = sD;
		else {
			sN = (-d + b);
			sD = a;
		}
	}
	// finally do the division to get sc and tc
	if (fabs(sN) < 0.00000001) sc = 0.0 ;
	else sc = sN / sD;
	if (fabs(tN) < 0.00000001) tc = 0.0 ;
	else tc = tN / tD;

	// get the difference of the two closest points
	//Vector = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
	vec.x = u.x*sc + w.x - v.x*tc;
	vec.y = u.y*sc + w.y - v.y*tc;
	vec.z = u.z*sc + w.z - v.z*tc;

	return vec;
}


/*..............................................................................*/

/*
   Find closest distance between line segment and point and return it as vector 
   (from point to closest segment point)
   Function gets orientation and length of line segments and the vector connecting
   their center os masses (from segment to point)
 */

struct vector mindist_segmentpoint(struct vector dir1, double length, struct vector r_cm)
{
	struct vector vec;
	double c,d,halfl;

	halfl=length*0.5;
	c = DOT(dir1,r_cm); 
	if (c >= halfl) d = halfl;
	else {
		if (c > -halfl) d = c;
		else d = -halfl;
	}

	vec.x = - r_cm.x + dir1.x * d;
	vec.y = - r_cm.y + dir1.y * d;
	vec.z = - r_cm.z + dir1.z * d;

	return vec;
}


/*..............................................................................*/

/*
   Determines whether two particles overlap.
   Returns 1 if there is an overlap, 0 if not.
 */

int overlap(struct particles part1, struct particles part2,
		struct vector box, struct ia_param ia_params[MAXT][MAXT])
{
	double b, c, d, e, f;   /* Coefficients in distance quadratic */
	double boundary;     /* Half length of central boundary zone of quadratic */
	double det;
	double halfl;        /* Half length of cylinder */
	double s0, t0;       /* det times location of min separation of infinite lines */
	double ss, tt;       /* Location of min separation of line segments */
	struct vector r_cm;  /* Vector between centres of mass */
	double dist;         /* Distance between particles*/
	struct vector distvec; /* Distance vector between particles*/

	double linemin(double, double);
	struct vector image(struct vector, struct vector, struct vector);

	r_cm = image(part1.pos, part2.pos, box);
	if ((part1.type >= SP) && (part2.type >= SP)) { /*we have two spheres - most common, do nothing*/
		dist=sqrt(DOT(r_cm,r_cm));
	} else {
		if ((ia_params[part1.type][part2.type].geotype[0] < SP) && (ia_params[part1.type][part2.type].geotype[1] < SP)) { /*we have two spherocylinders*/
			/*finding closes contact between them*/
			b = -DOT(part1.dir, part2.dir);
			d =  DOT(part1.dir, r_cm);
			e = -DOT(part2.dir, r_cm);
			f =  DOT(r_cm, r_cm);
			det = 1.0 - b*b;
			//halfl = length / 2.0;
			// Just take the mean
			halfl = ia_params[part1.type][part2.type].half_len[0] = ia_params[part1.type][part2.type].half_len[1];
			halfl /= 2;
			boundary = det * halfl;
			/* Location of smallest separation of the infinite lines */
			s0 = b*e - d;
			t0 = b*d - e;
			/* Location of smallest separation of line segments */
			if (s0 >= boundary) {
				if (t0 >= boundary) {
					/* Region 2 */
					if ( d + halfl + halfl*b < 0.0 ) {
						ss = halfl;
						tt = linemin( -ss*b - e, halfl );
					} else {
						tt = halfl;
						ss = linemin( -tt*b - d, halfl );
					}
				} else if (t0 >= -boundary) {
					/* Region 1 */
					ss = halfl;
					tt = linemin( -ss*b - e, halfl );
				} else {
					/* Region 8 */
					if ( d + halfl - halfl*b < 0.0 ) {
						ss = halfl;
						tt = linemin( -ss*b - e, halfl );
					} else {
						tt = -halfl;
						ss = linemin( -tt*b - d, halfl );
					}
				}
			} else if (s0 >= -boundary) {
				if (t0 >= boundary) {
					/* Region 3 */
					tt = halfl;
					ss = linemin( -tt*b - d, halfl );
				} else if (t0 >= -boundary) {
					/* Region 0 */
					ss = s0/det;
					tt = t0/det;
				} else {
					/* Region 7 */
					tt = -halfl;
					ss = linemin( -tt*b - d, halfl );
				}
			} else {
				if (t0 >= boundary) {
					/* Region 4 */
					if ( d - halfl + halfl*b > 0.0 ) {
						ss = -halfl;
						tt = linemin( -ss*b - e, halfl );
					} else {
						tt = halfl;
						ss = linemin( -tt*b - d, halfl );
					}
				} else if (t0 >= -boundary) {
					/* Region 5 */
					ss = -halfl;
					tt = linemin( -ss*b - e, halfl );
				} else {
					/* Region 6 */
					if ( d - halfl - halfl*b > 0.0 ) {
						ss = -halfl;
						tt = linemin( -ss*b - e, halfl );
					} else {
						tt = -halfl;
						ss = linemin( -tt*b - d, halfl );
					}
				}
			}
			/*ss snd tt are Location of min separation of line segments */
			dist=sqrt(f + ss*ss + tt*tt + 2.0*(ss*d + tt*e + ss*tt*b));
		} else {
			if (ia_params[part1.type][part2.type].geotype[0] < SP) { /*We have one spherocylinder -it is first one*/
				//halfl=length/2;/*finding closest vector from sphyrocylinder to sphere*/
				halfl=ia_params[part1.type][part2.type].half_len[0];/*finding closest vector from sphyrocylinder to sphere*/
				c = DOT(part1.dir,r_cm);
				if (c >= halfl) d = halfl;
				else {
					if (c > -halfl) d = c;
					else d = -halfl;
				}
				distvec.x = - r_cm.x + part1.dir.x * d;
				distvec.y = - r_cm.y + part1.dir.y * d;
				distvec.z = - r_cm.z + part1.dir.z * d;
				dist=sqrt(DOT(distvec,distvec));
			} else { /*lst option first one is sphere second one spherocylinder*/
				//halfl=length/2; /*finding closest vector from sphyrocylinder to sphere*/
				halfl=ia_params[part1.type][part2.type].half_len[1];/*finding closest vector from sphyrocylinder to sphere*/
				c = DOT(part2.dir,r_cm); 
				if (c >= halfl) d = halfl;
				else {
					if (c > -halfl) d = c;
					else d = -halfl;
				}
				distvec.x = r_cm.x - part2.dir.x * d;
				distvec.y = r_cm.y - part2.dir.y * d;
				distvec.z = r_cm.z - part2.dir.z * d;
				dist=sqrt(DOT(distvec,distvec));
			}
		}
	}

	/* Overlap exists if smallest separation is less than diameter of cylinder */
	if (dist < ia_params[part1.type][part2.type].sigma*0.5 ) {
		return 1;
	} else {
		return 0;
	}
}


/*..............................................................................*/

double linemin(double criterion, double halfl)
{
	if      (criterion >=  halfl) { return  halfl; }
	else if (criterion >= -halfl) { return  criterion; }
	else                          { return -halfl; }
}



/*..............................................................................*/
/*........................SOME USEFUL MATH......................................*/
/*..............................................................................*/

/*
   ran2 from Numerical Recipes.
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX 

/*..............................................................................*/

/*
   From Numerical Recipes.  Simplified to deal specifically with 3*3 matrices
   (stored as elements [1...3][1...3] or a 4*4 array).
 */

void tred2(double a[4][4], double d[4], double e[4])
{
	int l, k, j, i;
	double scale, hh, h, g, f;

	for (i=3; i>=2; i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++) scale += fabs(a[i][k]);
			if (scale == 0.0) e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
					/* a[j][i]=a[i][j]/h; */
					g=0.0;
					for (k=1;k<=j;k++) g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++) g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++) a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else e[i]=a[i][l];
		d[i]=h;
	}
	/* d[1]=0.0; */
	e[1]=0.0;
	for (i=1; i<=3; i++) {
		/* l=i-1;
		   if (d[i]) {
		   for (j=1;j<=l;j++) {
		   g=0.0;
		   for (k=1;k<=l;k++) g += a[i][k]*a[k][j];
		   for (k=1;k<=l;k++) a[k][j] -= g*a[k][i];
		   }
		   } */
		d[i]=a[i][i];
		/* a[i][i]=1.0;
		   for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0; */
	}
}  

/*..............................................................................*/

/*
   From Numerical Recipes.  Simplified to deal specifically with 3*3 matrices
   (stored as elements [1...3][1...3] or a 4*4 array).
 */

#define NRANSI
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void tqli(double d[4], double e[4])
{
	double pythag(double a, double b);
	int m, l, iter, i;
	/* int k; */
	double s, r, p, g, f, dd, c, b;

	for (i=2; i<=3; i++) e[i-1] = e[i];
	e[3] = 0.0;
	for (l=1; l<=3; l++) {
		iter = 0;
		do {
			for (m=l; m<=3-1; m++) {
				dd = fabs(d[m]) + fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) {
					fprintf(stderr, "Too many iterations in tqli\n");
					exit (2);
				}
				g = (d[l+1] - d[l]) / (2.0*e[l]);
				r = pythag(g, 1.0);
				g = d[m] - d[l] + e[l] / (g + SIGN(r,g));
				s = c = 1.0;
				p = 0.0;
				for (i=m-1; i>=l; i--) {
					f = s * e[i];
					b = c * e[i];
					e[i+1] = (r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f/r;
					c = g/r;
					g = d[i+1] - p;
					r = (d[i] - g)*s + 2.0*c*b;
					d[i+1] = g+(p=s*r);
					g = c*r - b;
					/* for (k=1; k<=3; k++) {
					   f = z[k][i+1];
					   z[k][i+1] = s*z[k][i]+c*f;
					   z[k][i] = c*z[k][i]i - s*f;
					   } */
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}
}

#undef NRANSI

/*..............................................................................*/

/*
   From Numerical Recipes.  Used by tqli.
 */

#define NRANSI
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

double pythag(double a, double b)
{
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

#undef NRANSI

/*..............................................................................*/

/*
   Normalise a vector to have unit length.  For speed during heavy use, it is
   not checked that the supplied vector has non-zero length.
 */

void normalise(struct vector *u)
{
	double tot;

	tot = sqrt( DOT(*u,*u) );
	if (tot !=0.0) {
		tot=1/tot;
		(*u).x *= tot;
		(*u).y *= tot;
		(*u).z *= tot;
	}
}


/*
   Returns the vector pointing from the centre of mass of particle 2 to the
   centre of mass of the closest image of particle 1.
 */

struct vector image(struct vector r1, struct vector r2, struct vector box)
{
	struct vector r12;

	double anint(double);

	r12.x = r1.x - r2.x;
	r12.y = r1.y - r2.y;
	r12.z = r1.z - r2.z;

	r12.x = box.x * (r12.x - anint(r12.x));
	r12.y = box.y * (r12.y - anint(r12.y));
	r12.z = box.z * (r12.z - anint(r12.z));

	return r12;
}


/*
   Returns the nearest integer to its argument as a double precision number. e.g.
   anint(-0.49) = 0.0 and anint(-0.51) = -1.0. Equivalent to the Fortran intrinsic
   ANINT.
 */

double anint(double arg)
{
	if (arg < 0) {
		return (double)( (long)(arg-0.5) );
	} else {
		return (double)( (long)(arg+0.5) );
	}
}

/*..............................................................................*/

/*
   Returns an evenly distributed random unit vector of unit length.  See Allen &
   Tildesley p349 or Frenkel & Smit p410.
   RANDOM VECTOR ON UNIT SPHERE
 */

struct vector ranvec(void)
{
	double a, b, xi1, xi2;
	struct vector unit;

	double ran2(long *);

	do {
		xi1 = 1.0 - 2.0*ran2(&seed);
		xi2 = 1.0 - 2.0*ran2(&seed);
		a = xi1*xi1 + xi2*xi2;
	} while (a > 1.0);

	b = 2.0 * sqrt(1.0 - a);
	unit.x = xi1 * b;
	unit.y = xi2 * b;
	unit.z = 1.0 - 2.0*a;

	return unit;
}

/**
 * returns a point randomly and evenly distributed inside of a unit sphere
 */
struct vector ranvecsph(void)
{
	struct vector ranvec;

	double ran2(long *);
	do{
		ranvec.x = 2 * ran2(&seed) - 1.0;
		ranvec.y = 2 * ran2(&seed) - 1.0;
		ranvec.z = 2 * ran2(&seed) - 1.0;
	} while(ranvec.x*ranvec.x +
			ranvec.y*ranvec.y +
			ranvec.z*ranvec.z >= 1);

	//printf("%lf\t%lf\t%lf\n", ranvec.x,ranvec.y,ranvec.z);
	return ranvec;
}


/****  some useful math *******/

struct vector vec_create(double x, double y, double z)
{
	struct vector newvec;

	newvec.x=x;
	newvec.y=y;
	newvec.z=z;

	return newvec;
}

struct vector vec_createarr(double a[3])
{
	struct vector newvec;

	newvec.x=a[0];
	newvec.y=a[1];
	newvec.z=a[2];

	return newvec;
}


double vec_dotproduct(struct vector A,struct vector B)
{
	double dp;

	dp = A.x*B.x + A.y*B.y + A.z*B.z;
	return dp;
}

/* vector projection of vector A to direction of B*/

struct vector vec_project(struct vector* A,struct vector* B)
{
	double dp;
	struct vector pr;

	dp = A->x*B->x + A->y*B->y + A->z*B->z;

	pr.x=B->x*dp;
	pr.y=B->y*dp;
	pr.z=B->z*dp;
	return pr;
}

void ortogonalise(struct vector *A, struct vector B)
{
	double dp;
	double vec_dotproduct(struct vector A,struct vector B);

	dp=vec_dotproduct(*A,B);

	(*A).x -= B.x * dp;
	(*A).y -= B.y * dp;
	(*A).z -= B.z * dp;
}


/* vector projection of vector A perpendicular to direction of B*/

struct vector vec_perpproject(struct vector *A,struct vector *B)
{
	struct vector pp;
	double dp;

	struct vector vec_project(struct vector *, struct vector*);

	dp=DOT((*A),(*B));

	pp.x = A->x - B->x*dp;
	pp.y = A->y - B->y*dp;
	pp.z = A->z - B->z*dp;
	//    fprintf (stderr, "pp x: %.8f y: %.8f z: %.8f \n",pp.x,pp.y,pp.z);
	return pp;
}


/*  returns a vector perpendicular to A
    nothing special about the vector except that it's one of the perpendicular options and is normalized
 */

struct vector vec_perp(struct vector A)
{
	double ratio,x,y;
	struct vector somevector;

	struct vector vec_create(double, double, double);
	struct vector vec_normalize(struct vector);
	void normalise(struct vector *);
	struct vector vec_crossproduct(struct vector, struct vector);

	x=A.x;
	y=A.y;
	if (x == 0) x=1;
	else {
		if (y == 0) y=1;
		else {
			ratio=y/x;
			y=x*ratio*2;
		}
	}
	somevector= vec_create(x, y, A.z);
	normalise(&somevector);
	return vec_crossproduct(A,somevector);
}

/* Perform the multiplication of a matrix A and a vector B where A is the
   first argument and B is the second argument.  The routine will
   return AxB*/

struct vector matrix_vec_multiply(double A[3][3],struct vector B)
{
	int i;
	double vecarr[3];
	struct vector AB,RA;

	struct vector vec_createarr(double[3]);
	double vec_dotproduct(struct vector,struct vector);

	for (i=0;i<3;i++) {
		/* index the row vector from A*/
		RA=vec_createarr(A[i]);
		/* Now find the dot product of this row with B*/
		vecarr[i]=vec_dotproduct(RA,B);
	}
	AB=vec_createarr(vecarr);
	return AB;
}


/* Distance between two vectors*/
double vec_distance(struct vector vec1,struct vector vec2)
{
	double sum;

	sum= (vec1.x-vec2.x)*(vec1.x-vec2.x)+(vec1.y-vec2.y)*(vec1.y-vec2.y)+(vec1.z-vec2.z)*(vec1.z-vec2.z);

	return pow(sum,0.5);
}

/* Vector size */
double vec_size(struct vector vec)
{
	double size;

	size=sqrt(vec.x*vec.x+ vec.y*vec.y+ vec.z*vec.z);

	return size;
}


/* Normalize a vector*/

struct vector vec_normalize(struct vector vec)
{
	double mag;
	struct vector newvec;
	double vec_size(struct vector);

	mag= vec_size (vec);
	mag=1/mag;
	newvec.x=vec.x*mag;
	newvec.y=vec.y*mag;
	newvec.z=vec.z*mag;

	return newvec;
}

/* Scale a vector */

struct vector vec_scale(struct vector vec, double scale)
{
	vec.x=vec.x*scale;
	vec.y=vec.y*scale;
	vec.z=vec.z*scale;

	return vec;
}

/* cross_product*/

struct vector vec_crossproduct(struct vector A,struct vector B)
{
	struct vector cp;

	cp.x=( A.y*B.z - A.z*B.y);
	cp.y=( -A.x*B.z + A.z*B.x);
	cp.z=( A.x*B.y - A.y*B.x);
	return cp;
}

/* addition of vectors*/
	inline
struct vector vec_sum(struct vector A,struct vector B)
{
	struct vector C;

	C.x=(A.x + B.x);
	C.y=(A.y + B.y);
	C.z=(A.z + B.z);
	return C;
}

/* subtraction of vectors*/
	inline
struct vector vec_sub(struct vector A,struct vector B)
{
	struct vector C;

	C.x=(A.x - B.x);
	C.y=(A.y - B.y);
	C.z=(A.z - B.z);
	return C;
}
/* asign vlues of vector A by values in vector B*/
	inline
void vec_asign(struct vector *A, struct vector B)
{
	(*A).x=B.x;
	(*A).y=B.y;
	(*A).z=B.z;

}

/* generate random unit vector*/
struct vector vec_random(void)
{
	struct vector newvec;
	struct vector ranvec(void);

	newvec=ranvec();

	return newvec;
}

/*generate random unit quaternion*/
struct quat quat_random(void)
{
	double cosv, sinv;
	struct quat newquat;
	struct vector newaxis;
	struct vector ranvec(void);

	/* generate quaternion for rotation*/
	newaxis = ranvec(); /*random axes for rotation*/
	cosv = cos(PIH * ran2(&seed) );
	if (ran2(&seed) <0.5) sinv = sqrt(1.0 - cosv*cosv);
	else sinv = -sqrt(1.0 - cosv*cosv);
	newquat.w=cosv;
	newquat.x=newaxis.x*sinv;
	newquat.y=newaxis.y*sinv;
	newquat.z=newaxis.z*sinv;

	return newquat;

}
/* Create quaternion for rotation around vector "vec" of angle in degrees "angle" 
   function need cos of half angle and its sin*/
struct quat quat_create(struct vector vec, double vc, double vs)
{
	struct quat newquat;

	newquat.w=vc;
	newquat.x=vec.x*vs;
	newquat.y=vec.y*vs;
	newquat.z=vec.z*vs;

	return newquat;
}

/*rotate vector with quaternion*/
void vec_rotate(struct vector *vec, struct quat quat)
{
	double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz;

	/*    t1 = quat.w * quat.w; */
	t2 =  quat.w * quat.x;
	t3 =  quat.w * quat.y;
	t4 =  quat.w * quat.z;
	t5 = -quat.x * quat.x;
	t6 =  quat.x * quat.y;
	t7 =  quat.x * quat.z;
	t8 = -quat.y * quat.y;
	t9 =  quat.y * quat.z;
	t10 = -quat.z * quat.z;
	newx = 2.0 * ( (t8+t10)*(*vec).x + (t6-t4)*(*vec).y + (t3+t7)*(*vec).z ) + (*vec).x;
	newy = 2.0 * ( (t4+t6)*(*vec).x + (t5+t10)*(*vec).y + (t9-t2)*(*vec).z ) + (*vec).y;
	newz = 2.0 * ( (t7-t3)*(*vec).x + (t2+t9)*(*vec).y + (t5+t8)*(*vec).z ) + (*vec).z;

	(*vec).x = newx;
	(*vec).y = newy;
	(*vec).z = newz;

}

/* rotate spherocylinder by quaternion of random axis and angle smaller than
   maxcos(cosine of angle half), we do everything on site for speed */
void psc_rotate(struct particles *psc, double max_angle,int geotype)
{
	double vc, vs, t2, t3, t4, t5, t6, t7, t8, t9, t10;
	double d1, d2, d3, d4, d5, d6, d7, d8, d9 , newx, newy, newz;
	int k,m;
	struct quat newquat;
	struct vector newaxis;
	struct vector ranvec(void);

	/* generate quaternion for rotation*/
	newaxis = ranvec(); /*random axes for rotation*/
	//    maxcos = cos(maxorient/2/180*PI);
	// vc = maxcos + ran2(&seed)*(1-maxcos); /*cos of angle must be bigger than maxcos and smaller than one*/
	vc = cos(max_angle * ran2(&seed) );
	if (ran2(&seed) <0.5) vs = sqrt(1.0 - vc*vc);
	else vs = -sqrt(1.0 - vc*vc); /*randomly choose orientation of direction of rotation clockwise or counterclockwise*/
	newquat.w=vc;
	newquat.x=newaxis.x*vs;
	newquat.y=newaxis.y*vs;
	newquat.z=newaxis.z*vs;

	/* do quaternion rotation*/
	t2 =  newquat.w * newquat.x;
	t3 =  newquat.w * newquat.y;
	t4 =  newquat.w * newquat.z;
	t5 = -newquat.x * newquat.x;
	t6 =  newquat.x * newquat.y;
	t7 =  newquat.x * newquat.z;
	t8 = -newquat.y * newquat.y;
	t9 =  newquat.y * newquat.z;
	t10 = -newquat.z * newquat.z;

	d1 = t8 + t10;
	d2 = t6 - t4;
	d3 = t3 + t7;
	d4 = t4 + t6;
	d5 = t5 + t10;
	d6 = t9 - t2;
	d7 = t7 - t3;
	d8 = t2 + t9;
	d9 = t5 + t8;

	/*rotate spherocylinder direction vector*/
	newx = 2.0 * ( d1*psc->dir.x + d2*psc->dir.y + d3*psc->dir.z ) + psc->dir.x;
	newy = 2.0 * ( d4*psc->dir.x + d5*psc->dir.y + d6*psc->dir.z ) + psc->dir.y;
	newz = 2.0 * ( d7*psc->dir.x + d8*psc->dir.y + d9*psc->dir.z ) + psc->dir.z;
	psc->dir.x = newx;
	psc->dir.y = newy;
	psc->dir.z = newz;
	
	m=1;
	if ( (geotype != SCN) && (geotype != SCA) ) {
	    if ( (geotype == TPSC) || (geotype == TCPSC) || (geotype == TCHPSC) || (geotype == TCHCPSC) )
		m=2;
	    for (k=0;k<m;k++) {
		/*rotate patch direction vector*/
		newx = 2.0 * ( d1*psc->patchdir[k].x + d2*psc->patchdir[k].y + d3*psc->patchdir[k].z ) + psc->patchdir[k].x;
		newy = 2.0 * ( d4*psc->patchdir[k].x + d5*psc->patchdir[k].y + d6*psc->patchdir[k].z ) + psc->patchdir[k].y;
		newz = 2.0 * ( d7*psc->patchdir[k].x + d8*psc->patchdir[k].y + d9*psc->patchdir[k].z ) + psc->patchdir[k].z;
		psc->patchdir[k].x = newx;
		psc->patchdir[k].y = newy;
		psc->patchdir[k].z = newz;
	
		/*rotate patch sides vectors*/
		newx = 2.0 * ( d1*psc->patchsides[0+2*k].x + d2*psc->patchsides[0+2*k].y + d3*psc->patchsides[0+2*k].z ) + psc->patchsides[0+2*k].x;
		newy = 2.0 * ( d4*psc->patchsides[0+2*k].x + d5*psc->patchsides[0+2*k].y + d6*psc->patchsides[0+2*k].z ) + psc->patchsides[0+2*k].y;
		newz = 2.0 * ( d7*psc->patchsides[0+2*k].x + d8*psc->patchsides[0+2*k].y + d9*psc->patchsides[0+2*k].z ) + psc->patchsides[0+2*k].z;
		psc->patchsides[0+2*k].x = newx;
		psc->patchsides[0+2*k].y = newy;
		psc->patchsides[0+2*k].z = newz;
		newx = 2.0 * ( d1*psc->patchsides[1+2*k].x + d2*psc->patchsides[1+2*k].y + d3*psc->patchsides[1+2*k].z ) + psc->patchsides[1+2*k].x;
		newy = 2.0 * ( d4*psc->patchsides[1+2*k].x + d5*psc->patchsides[1+2*k].y + d6*psc->patchsides[1+2*k].z ) + psc->patchsides[1+2*k].y;
		newz = 2.0 * ( d7*psc->patchsides[1+2*k].x + d8*psc->patchsides[1+2*k].y + d9*psc->patchsides[1+2*k].z ) + psc->patchsides[1+2*k].z;
		psc->patchsides[1+2*k].x = newx;
		psc->patchsides[1+2*k].y = newy;
		psc->patchsides[1+2*k].z = newz;
	    }	
	}
	
	m=1;
	if ( (geotype == CHPSC) || (geotype == CHCPSC) || (geotype == TCHPSC) || (geotype == TCHCPSC) ) {
	    if ( (geotype == TCHPSC) || (geotype == TCHCPSC) ) 
		m=2;
	    for (k=0;k<m;k++) {
		/*rotate chiral direction vector*/
		newx = 2.0 * ( d1*psc->chdir[k].x + d2*psc->chdir[k].y + d3*psc->chdir[k].z ) + psc->chdir[k].x;
		newy = 2.0 * ( d4*psc->chdir[k].x + d5*psc->chdir[k].y + d6*psc->chdir[k].z ) + psc->chdir[k].y;
		newz = 2.0 * ( d7*psc->chdir[k].x + d8*psc->chdir[k].y + d9*psc->chdir[k].z ) + psc->chdir[k].z;
		psc->chdir[k].x = newx;
		psc->chdir[k].y = newy;
		psc->chdir[k].z = newz;
	    }
	}

}

/*returns a position of center of mass of system*/
void masscenter(long npart, struct ia_param ia_params[MAXT][MAXT], struct conf * conf)
{
	long i;

	double anint(double);

	conf->syscm.x = 0;
	conf->syscm.y = 0;
	conf->syscm.z = 0;
	for (i=0; i<npart; i++) {
		/*using periodic boundary conditions*/
		conf->syscm.x += (conf->particle[i].pos.x - anint(conf->particle[i].pos.x) ) * 
			ia_params[conf->particle[i].type][conf->particle[i].type].volume;
		conf->syscm.y += (conf->particle[i].pos.y - anint(conf->particle[i].pos.y) ) * 
			ia_params[conf->particle[i].type][conf->particle[i].type].volume;
		conf->syscm.z += (conf->particle[i].pos.z - anint(conf->particle[i].pos.z) ) * 
			ia_params[conf->particle[i].type][conf->particle[i].type].volume;
	}
	
	conf->syscm.x /= conf->sysvolume;
	conf->syscm.y /= conf->sysvolume;
	conf->syscm.z /= conf->sysvolume;
	return;
}

/* rotate cluster of particles by quaternion of random axis and angle smaller than 
   maxcos(cosine of angle half), we do everything on site for speed */

void cluster_rotate(long target, struct vector gc, double max_angle, struct topo * topo, struct conf * conf)
{ 
	long current,i;
	double vc,vs;
	//double quatsize;
	
	struct quat newquat;
	struct vector newaxis;

	struct vector ranvec(void);
	void vec_rotate(struct vector *, struct quat);

	// create rotation quaternion
	newaxis = ranvec(); /*random axes for rotation*/
	//    maxcos = cos(maxorient/2/180*PI);
	//vc = maxcos + ran2(&seed)*(1-maxcos); /*cos of angle must be bigger than maxcos and smaller than one*/
	vc = cos(max_angle * ran2(&seed) );
	if (ran2(&seed) <0.5) vs = sqrt(1.0 - vc*vc);
	else vs = -sqrt(1.0 - vc*vc); /*randomly choose orientation of direction of rotation clockwise or counterclockwise*/
	newquat.w=vc;
	newquat.x=newaxis.x*vs;
	newquat.y=newaxis.y*vs;
	newquat.z=newaxis.z*vs;
	//quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);
	
	//shift position to geometrical center
	i=0;
	current = topo->chainlist[target][0];
	while (current >=0 ) {
		//shift position to geometrical center
		conf->particle[current].pos.x -= gc.x;
		conf->particle[current].pos.y -= gc.y;
		conf->particle[current].pos.z -= gc.z;
		//scale things by box not to have them distorted
		conf->particle[current].pos.x *= conf->box.x;
		conf->particle[current].pos.y *= conf->box.y;
		conf->particle[current].pos.z *= conf->box.z;
		//do rotation
		vec_rotate(&conf->particle[current].pos, newquat);
		vec_rotate(&conf->particle[current].dir, newquat);
		vec_rotate(&conf->particle[current].patchdir[0], newquat);
		vec_rotate(&conf->particle[current].patchdir[1], newquat);
		vec_rotate(&conf->particle[current].chdir[0], newquat);
		vec_rotate(&conf->particle[current].chdir[1], newquat);
		vec_rotate(&conf->particle[current].patchsides[0], newquat);
		vec_rotate(&conf->particle[current].patchsides[1], newquat);
		vec_rotate(&conf->particle[current].patchsides[2], newquat);
		vec_rotate(&conf->particle[current].patchsides[3], newquat);
		//sclae back
		conf->particle[current].pos.x /= conf->box.x;
		conf->particle[current].pos.y /= conf->box.y;
		conf->particle[current].pos.z /= conf->box.z;		
		//shift positions back
		conf->particle[current].pos.x += gc.x;
		conf->particle[current].pos.y += gc.y;
		conf->particle[current].pos.z += gc.z;
		i++;
		current = topo->chainlist[target][i];
	}

}

/* put the particle in the original box using periodic boundary conditions
   in our system the particle positions are scaled by box size so to get them
   into original obx is to get htem between 0 and 1 and then scale this back 
   by size of box*/
void origbox(struct vector *pos,struct vector box)
{
	double anint(double);

	(*pos).x = box.x * ((*pos).x - anint((*pos).x));
	(*pos).y = box.y * ((*pos).y - anint((*pos).y));
	(*pos).z = box.z * ((*pos).z - anint((*pos).z));
}

/* use of periodic boundary conditions*/
void usepbc(struct vector *pos,struct vector pbc)
{

	do {
		(*pos).x += pbc.x;
	} while ((*pos).x < 0.0);
	do {
		(*pos).x -= pbc.x;
	} while ((*pos).x > pbc.x);

	do {
		(*pos).y += pbc.y;
	} while ((*pos).y < 0.0);
	do {
		(*pos).y -= pbc.y;
	} while ((*pos).y > pbc.y);

	do {
		(*pos).z += pbc.z;
	} while ((*pos).z < 0.0);
	do {
		(*pos).z -= pbc.z;
	} while ((*pos).z > pbc.z);

}

/*..............................................................................*/
/*.......................TEMPLATE FILES.........................................*/
/*..............................................................................*/

/*
# Template for the "options" file. Options start with an '#'.
# Pressure couplings:
#   0 = anisotropic coupling, 1 = isotropic coupling, 2 = isotropic in xy z=const, 3 = isotropic 
#   xy and keep Volume constant
# Wang-Landau method: (with constant decrease of bias addition by factor of 2, until less than WL_ALPHATOL)
#   O = none, 1 = z-direction of 1st paticle, 2 = hole in xyplane, 3 = z-orientation of 0th particle 
#   4 = distance of first two particles, 5 = pore around z axis and above CM, 6 = pore around z axis and above 0th particle
#   7 = number of particles in contact (within distance sqrt(WL_CONTACTS))
ptype = 1               #  Pressure coupling type (0-anisotropic xyz, 1-isotropic xyz, 2 - isotropic in xy z=const, 3 - isotropic in xy and V=const)
press = 1               #  Pressure
paralpress = 1          #  Parallel pressure for replica exchange
shave = 0               #  Average number of volume change attempts per sweep (usually 1)
nequil = 0              #  Number of equilibration sweeps
adjust = 0              #  Number of equilibration sweeps between step size adjustments
nsweeps  = 1000000      #  Number of production sweeps
paramfrq = 1000000      #  Number of sweeps between order parameter samples
report   = 1000000      #  Number of sweeps between statistics reports
nrepchange = 1000       #  Number of sweeps between replica exchanges
movie    = 100000       #  Number of sweeps between movie frames (0 = no movie)
chainprob = 0.0         #  Probability of chain move attempts per sweep ( 0.25/number of particles in chain)
transmx = 0.212         #  Initial maximum displacement
rotmx = 7.5             #  Initial maximum orientation change (degrees)
edge_mx = 0.0           #  Initial maximum box length change
chainmmx = 0.0          #  Initial maximum chain displacement
chainrmx = 0.0          #  Initial maximum chain rotation change (degrees)
temper = 1.0            #  Temperature in units kT/e
paraltemper = 1.5       #  Temperature for parallel tempering in kT/e 
wlm = 0                 #  Wang-Landau method
wlmtype = 0             #  For which atomic type (from top.init) should the Wang-Landau method be calculated?
switchprob = 0.0016     #  Probability of type switch attempts per sweep
pairlist_update = 8     #  Number of sweeps after which the pairlist should be updated
seed = 1                #  Random number seed
write_cluster = 10000   #  Number of sweeps per writing out cluster info
# End of the file

 */


/*
   Example of 'Config.init' file, but you must delete comments... there are only number in configuration file
#box
10.0 10.0 10.0
#particles (x,y,z) (direction_x,direction_y, direction_z) (patchdirection_x,patchdirection_y,patchdirection_z) (switched)

 */


/*
   Template for the topology file 'top.init'. ( "\\" is symbol for line continue,
   "#" is symbol for comment, "[" is starting sign for keyword, "]" is ending sign 
   for  kyeword ) There are three keywords, types, molecules, and system. They
   should be given in this order.
TYPES:
spherocylinders
SC - purely repulsive spherocylinder with WCA potential on closest distance
SCA - isotropic cos^2 potential is acting isotropicaly dependent only on
closest distance between spherocylinders.. 
PSC - Attractive potential in limited to an angular wedge on spherocylinder. Patch
goes all the way through, making also hemispherical caps on end attractive
CPSC - Attractive potential in limited to an angular wedge on cylindrical part
of spherocylinders. The hemispherical caps on ends are repulsive
spheres
(T)(CH)PSC - T adds second patch, CH - adds chirality
SP - purely repulsive shpere with WCA potential on closest distance
SPA - isotropic cos^2 potential is acting isotropicaly dependent only on
closest distance between obejcts


[Types]
# NAME  NUMBER    GEOTYPE  EPSILON  SIGMA    ATTRACTION_DIST   ATTRACTION_SWITCH PATCH_ANGLE PATCH_SWITCH   SC_LENGTH  (Optional second patch: PATCH_ROTATION PATCH_ANGLE PATCH_SWITCH )CHIRAL_ANGLE
Prot1   1         PSC      1        1.2      1.346954458       1.0               80.0        5.0            3
Prot2   2         PSC      1        1.2      1.346954458       1.0               170.0       5.0            3
Prot3   3         CHCPSC   1        1.2      1.346954458       1.0               170.0       5.0            3             10
Prot4   4         TCHCPSC  1        1.2      1.346954458       1.0               170.0       5.0            3             90.0    90.0     5.0    10
[Molecules]
# Molecules letter
#  bond1 - harmonic bond between nearest neighbours (end points for spherocylinders) (first constant then eq distance)
#  bond2 - harmonic bond between second nearest neighbours (their center of mass) (first constant then eq distance)
#  bondd - directional harmonic bond between nearest neighbours (end point of the second spherocylinder is attached to the point of bondlength extension of the first spherocylinder) (first constant then eq distance) 
#  angle1 - angle between two spherocylinders -nearest neighbours (first constant then eq degrees 0-180.0)
#  angle2 - angle between two spherocylinder patches -nearest neighbours (first constant then eq degrees 0-180.0)
#  particles - types as they go in chain in molecule
A: {
#what:       TYPE     SWITCHTYPE  DELTA_MU
particles:   1        2           0.5
particles:   2
}
B: {
particles:  1
particles:  2         1           0.3
}
[System]
A 2
B 2
[EXTER]
# wall interaction
# THICKNESS  EPSILON   ATTRACTION_SWITCH 
5.0        1.0           1.0
[EXCLUDE]
#set pair types for which attraction will be excluded  (reversepair is automaticaly added)
 1  2
 1  3
 */

