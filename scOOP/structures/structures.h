/** @file structures.h*/

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "macros.h"
#include "../mc/math_calc.h"
#include <vector>


/**
 * @brief contains all the particles of one cluster
 */
typedef struct{
    long npart;
    long * particles;
} Cluster;

/**
 * @brief This structure is for io only
 */
typedef struct{
    char * name;                  // The name of the molecule
    int * type;                   // The type of the particle
    long * switchtype;            // The switchtype of the particle
    double * delta_mu;            // The chemical potential for the switch
} Molecule;

/**
 * @brief Holds the type of a variable in struct option
 */
typedef enum {
    Int,
    Int2,
    Long,
    Double
} Type;

typedef struct {        // for reading in the options
    const char *id;     // The name of the value in the option file
    Type type;          // The type (int, double or long)
    bool set;           // Wheter the variable has been set
    void *var;          // The variable
} Option;

typedef struct{
    // input files
    char configurationInFile[30];
    char topologyInFile[30];
    char optionsfile[30];
    char wlinfile[30];
    // output files
    char configurationoutfile[30];
    char topologyOutFile[30];
    char moviefile[30];
    char wloutfile[30];
    char statfile[30];
    char clusterfile[30];
    char clusterstatfile[30];
    char energyfile[30];
} FileNames;

typedef struct {               /* Define step size and acceptance ratio statistics */
    double mx;               /* Maximum value displacement, cos(angle), etc.  */
    double angle;            /* Maximum angle, since in .mx cos(angle) is saved */
    long acc;                /* Number of accepted steps */
    long rej;                /* Number of rejected steps */
    double oldrmsd;          /* Averaged mx value in previous equilibration round */
    double oldmx;            /* Change in mx in last equlibrium step */
} Disp;

typedef struct{      /*Parameters for inner interaction in chains*/
    char * name;

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

    double mu;                      // chemical potential, specific to each molecule type
    double lnThermalWavelengh;      // ln (de Broglie thermal wavelenght), specific to each type
    bool muVTmove;                  // true - attempt muVT move for that type
    int particleTypes[MAXCHL];      // -1 - not a single particle, 0..40 single particle type

} Chainparams;



typedef struct{      /* Contatins properties and parameters of particle types */
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
    double csecpatchrot[2];     /* Cosine of Rotation of second patches in 2psc models*/
    double ssecpatchrot[2];     /* Sine of Rotation of second patches in 2psc models*/
    double volume;              /* Volume of particle for geometrical center calculations*/
    double pvolscale;           /* Scale of patch volume size*/
    double len[2];              /* Length of the PSC */
    double half_len[2];         /* Half length of the PSC */
    double chiral_cos[2];       /* Coctains the cosinus for the chiral rotation of the patch */
    double chiral_sin[2];       /* Contains the sinus for the chiral rotation of the patch */
} Ia_param;

typedef struct{
    bool exist;                           /* existence of external potential*/
    double thickness;                     /* external wall thicnkess*/
    double epsilon;                       /* depth of attraction*/
    double attraction;                    /* distance of attraction*/
    double sqmaxcut;                      /* distance when nothing can interact*/
    Ia_param interactions[MAXT];   /* Interaction parameters with particle types generated from above params*/
} Exters;

typedef struct{                /* It would be nice, if this struct would contain all the topo stuff in the end*/
    long * switchlist;                   /* List containing the number of all the particles with switchtypes */
    long n_switch_part;                  /* number of particles with switchtype */
    double sqmaxcut;                     /* square of distance over which even spherocylinders cannot interact (distance between CM) */
    double maxcut;                       /* distance over which even spherocylinders cannot interact (distance between CM) */

    /** @brief parameters for chains */
    Chainparams chainparam[MAXMT];
    int molTypeCount;

    Ia_param ia_params[MAXT][MAXT];      /* parametrization of particles for all interations*/
    Exters exter;                        /* external potential - wall */

    /** @brief List of particles in chain */
    long chainlist[MAXN][MAXCHL];

    /** @brief Number of chains */
    long chainCount;

    /** @brief list of particles with same type, size [MAXMT][MAXN], stack size insufficient */
    long** typeList;

    /** @brief number of particles with same type */
    long typeCount[MAXMT];

} Topo;







#endif // STRUCTURES_H
