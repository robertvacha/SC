/** @file structures.h*/

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "macros.h"
#include "../mc/math_calc.h"
#include <vector>
#include <cstdlib>
#include "moleculeparams.h"

#include <cstdio>

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
class Molecule {
public:
    char * name;        ///< \brief The name of the molecule
    int * type;         ///< \brief The type of the particle
    long * switchtype;  ///< \brief The switchtype of the particle
    double * delta_mu;  ///< \brief The chemical potential for the switch

    Molecule() {
        name = NULL;
        type = (int*) malloc(sizeof(long)*MAXN);
        switchtype = (long int*) malloc(sizeof(long)*MAXN);
        delta_mu = (double*) malloc(sizeof(double)*MAXN);
        for (int j=0;j<MAXN;j++) {
            type[j] = -1;
        }
    }

    ~Molecule() {
        free(name);
        free(type);
        free(switchtype);
        free(delta_mu);
    }
};




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
    bool set;           // Whether the variable has been set
    void *var;          // The variable
} Option;




typedef struct{
    // input files
    char configurationInFileMuVTChains[30];
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




/**
  * @brief Define step size and acceptance ratio statistics
  **/
class Disp {
public:
    double mx;          ///< \brief Maximum value displacement, cos(angle), etc.
    double angle;       ///< \brief Maximum angle, since in .mx cos(angle) is saved
    long acc;           ///< \brief Number of accepted steps
    long rej;           ///< \brief Number of rejected steps
    double oldrmsd;     ///< \brief Averaged mx value in previous equilibration round
    double oldmx;       ///< \brief Change in mx in last equlibrium step

    Disp() : mx(0.0), angle(0.0), acc(0), rej(0), oldrmsd(0.0), oldmx(0.0) {}
};




/**
  * @brief Contatins properties and parameters of particle types
  **/
class Ia_param{
public:
    char name[SMSTR];           ///< \brief The name of the particle type
    char other_name[SMSTR];     ///< \brief The name of the particle type
    int geotype[2];             ///< \brief The geometrical type:
                                /// spherocylinder (0-repulsive, 1-isotropic, 2-patchy, 3-cylindrical)
                                /// or sphere (0-repulsive, 1-isotropic)
    double sigma;               ///< \brief Repulsion wca
    double epsilon;             ///< \brief Repulsion strength
    double pdis;                ///< \brief Interaction distance of patch
    double pswitch;             ///< \brief Switch of distance of patch
    double pangl[4];            ///< \brief angular size of patch as was specifid in input
    double panglsw[4];          ///< \brief angular size of patchswitch as was specifid in input
    double pcangl[4];           ///< \brief cosine of half size angle - rotation from patch direction to side
    double pcanglsw[4];         ///< \brief cosine of half size angle plus switch - rotation from patch direction to side
    double rcut;                ///< \brief Cutoff for attraction
    double rcutwca;             ///< \brief Cutoff for repulsion
    double pcoshalfi[4];        ///< \brief Cosine of half angle going to side of interaction
    double psinhalfi[4];        ///< \brief Sine of half angle going to side of interaction -useful for quaterion rotation
    double csecpatchrot[2];     ///< \brief Cosine of Rotation of second patches in 2psc models
    double ssecpatchrot[2];     ///< \brief Sine of Rotation of second patches in 2psc models
    double volume;              ///< \brief Volume of particle for geometrical center calculations
    double pvolscale;           ///< \brief Scale of patch volume size
    double len[2];              ///< \brief Length of the PSC
    double half_len[2];         ///< \brief Half length of the PSC
    double chiral_cos[2];       ///< \brief Coctains the cosinus for the chiral rotation of the patch
    double chiral_sin[2];       ///< \brief Contains the sinus for the chiral rotation of the patch

    Ia_param() {
        for(int k = 0; k < 2; k++){
            geotype[k] = 0;
            len[k] = 0.0;
            half_len[k] = 0.0;
            chiral_cos[k] = 0.0;
            chiral_sin[k] = 0.0;
            csecpatchrot[k] = 0.0;
            ssecpatchrot[k] = 0.0;
        }
        for(int k = 0; k < 4; k++){
            pangl[k] = 0.0;
            panglsw[k] = 0.0;
            pcangl[k] = 0.0;
            pcanglsw[k] = 0.0;
            pcoshalfi[k] = 0.0;
            psinhalfi[k] = 0.0;
        }
        sigma = 0.0;
        epsilon = 0.0;
        rcutwca = 0.0;
        pdis = 0.0;
        pswitch = 0.0;
        rcut = 0.0;
        volume = 0.0;
        pvolscale = 0.0;
    }
};




class Exters{
public:
    bool exist;                    ///< \brief existence of external potential
    double thickness;              ///< \brief external wall thicnkess
    double epsilon;                ///< \brief depth of attraction
    double attraction;             ///< \brief distance of attraction
    double sqmaxcut;               ///< \brief distance when nothing can interact
    Ia_param interactions[MAXT];   ///< \brief Interaction parameters with particle types generated from above params

    Exters() {
        exist = false;
        thickness = 0.0;
        epsilon = 0.0;
        attraction = 0.0;
        sqmaxcut = 0.0;
    }
};



#endif // STRUCTURES_H
