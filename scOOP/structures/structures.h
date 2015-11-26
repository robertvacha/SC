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

class Molecule : public vector<int > {};

/**
 * @brief This structure is for io only
 */
class MolIO {
public:
    char * name;        ///< \brief The name of the molecule
    int * type;         ///< \brief The type of the particle
    long * switchtype;  ///< \brief The switchtype of the particle
    double * delta_mu;  ///< \brief The chemical potential for the switch
    long * count;       ///< \brief count of molecules

    MolIO() {
        name = NULL;
        type = (int*) malloc(sizeof(long)*MAXMT);
        switchtype = (long int*) malloc(sizeof(long)*MAXMT);
        delta_mu = (double*) malloc(sizeof(double)*MAXMT);
        for (int j=0;j<MAXMT;j++) {
            type[j] = -1;
        }
    }

    ~MolIO() {
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




class FileNames{
public:
    FileNames(){}
    // for MPI
    // input files
    char configurationPool[30];
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

    void initMPIRank(int rank) {
        // in files
        sprintf(configurationPool, "%dpool", rank);
        sprintf(configurationInFile, "%dconfig.init", rank);
        // topology in file -> change only for grand  canonical paralel tempering
        // options file is the same
        sprintf(wlinfile, "%dwl.dat", rank);

        // out files
        sprintf(configurationoutfile, "%dconfig.last", rank);
        // topology out file -> change only for grand  canonical paralel tempering
        sprintf(moviefile, "%dmovie", rank);
        sprintf(wloutfile, "%dwl-new.dat", rank);
        sprintf(statfile, "%dstat.dat", rank);
        sprintf(clusterfile, "%dcluster.dat", rank);
        sprintf(clusterstatfile, "%dcluster_stat.dat", rank);
        sprintf(energyfile, "%denergy.dat", rank);
    }
};




class StatsGrand {
public:
    StatsGrand() : delAcc(0), delRej(0), insAcc(0), insRej(0),
    muVtAverageParticles(0), muVtSteps(0) {}

    int delAcc;
    int delRej;
    int insAcc;
    int insRej;

    unsigned long long int muVtAverageParticles;
    unsigned int muVtSteps;

#ifdef ENABLE_MPI
    MPI_Datatype* defDataType(MPI_Datatype* MPI_stat) {

        MPI_Aint     dispstart;

        MPI_Datatype mpiexdataType[6] = {MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_LONG, MPI_INT};
        int          mpiexdataLen[6]  = {1,          1,          1,          1,          1,        1};

        MPI_Aint     mpiexdata[6];
        MPI_Address( this, &dispstart);
        MPI_Address( &(this->delAcc), &mpiexdata[0]);
        MPI_Address( &(this->delRej), &mpiexdata[1]);
        MPI_Address( &(this->insAcc), &mpiexdata[2]);

        MPI_Address( &(this->insRej), &mpiexdata[3]);
        MPI_Address( &(this->muVtAverageParticles), &mpiexdata[4]);
        MPI_Address( &(this->muVtSteps), &mpiexdata[5]);

        for (int i=0; i <6; i++)
            mpiexdata[i] -= dispstart;

        MPI_Type_struct(6, mpiexdataLen, mpiexdata, mpiexdataType, MPI_stat);
        MPI_Type_commit( MPI_stat);

        return MPI_stat;
    }
#endif
};

/**
  * @brief Define step size and acceptance ratio statistics
  **/
class Disp {
public:
    double mx;          ///< \brief Maximum value displacement, cos(angle), etc.
    double angle;       ///< \brief Maximum angle, since in .mx cos(angle) is saved
    double oldrmsd;     ///< \brief Averaged mx value in previous equilibration round
    double oldmx;       ///< \brief Change in mx in last equlibrium step
    long acc;           ///< \brief Number of accepted steps
    long rej;           ///< \brief Number of rejected steps

    Disp() : mx(0.0), angle(0.0), oldrmsd(0.0), oldmx(0.0), acc(0), rej(0) {}

double ratio() {
    if(acc + rej > 0) {
        return 1.0*acc/(acc+rej);
    }
    else return 0.0;
}

#ifdef ENABLE_MPI
    MPI_Datatype* defDataType(MPI_Datatype* MPI_stat) {

        MPI_Aint     dispstart;

        MPI_Datatype mpiexdataType[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_LONG};
        int          mpiexdataLen[6]  = {1,          1,          1,          1,          1,        1};

        MPI_Aint     mpiexdata[6];
        MPI_Address( this, &dispstart);
        MPI_Address( &(this->mx), &mpiexdata[0]);
        MPI_Address( &(this->angle), &mpiexdata[1]);
        MPI_Address( &(this->oldrmsd), &mpiexdata[2]);

        MPI_Address( &(this->oldmx), &mpiexdata[3]);
        MPI_Address( &(this->acc), &mpiexdata[4]);
        MPI_Address( &(this->rej), &mpiexdata[5]);

        for (int i=0; i <6; i++)
            mpiexdata[i] -= dispstart;

        MPI_Type_struct(6, mpiexdataLen, mpiexdata, mpiexdataType, MPI_stat);
        MPI_Type_commit( MPI_stat);

        return MPI_stat;
    }
#endif
};


class Statistics {
public:
    Statistics (){}

    StatsGrand grand[MAXMT];

    Disp edge;                  ///< \brief Maximum box length change and statistics
    Disp rot[MAXT];             ///< \brief Maximum rotation and statistics
    Disp trans[MAXT];           ///< \brief Maximum translation  and statistics
    Disp chainm[MAXMT];         ///< \brief Maximum translation for chain  and statistics
    Disp chainr[MAXMT];         ///< \brief Maximum rotation for chain and statistics
    Disp mpiexch;               ///< \brief MPI statistics

#ifdef ENABLE_MPI
    MPI_Datatype* defDataType(MPI_Datatype* MPI_stat) {

        StatsGrand g;
        Disp d;

        MPI_Datatype MPI_DISP, MPI_GRAND;

        d.defDataType(&MPI_DISP);
        g.defDataType(&MPI_GRAND);

        MPI_Aint     dispstart;

        MPI_Datatype mpiexdataType[7] = {MPI_GRAND, MPI_DISP, MPI_DISP, MPI_DISP, MPI_DISP, MPI_DISP, MPI_DISP};
        int          mpiexdataLen[7]  = {MAXMT,     1,        MAXT,     MAXT,     MAXMT,    MAXMT,    1};

        MPI_Aint     mpiexdata[7];
        MPI_Address( this, &dispstart);
        MPI_Address( &(this->grand), &mpiexdata[0]);

        MPI_Address( &(this->edge), &mpiexdata[1]);
        MPI_Address( &(this->rot), &mpiexdata[2]);
        MPI_Address( &(this->trans), &mpiexdata[3]);
        MPI_Address( &(this->chainm), &mpiexdata[4]);
        MPI_Address( &(this->chainr), &mpiexdata[5]);
        MPI_Address( &(this->mpiexch), &mpiexdata[6]);

        for (int i=0; i <7; i++)
            mpiexdata[i] -= dispstart;

        MPI_Type_struct(7, mpiexdataLen, mpiexdata, mpiexdataType, MPI_stat);
        MPI_Type_commit( MPI_stat);

        MPI_Type_free(&MPI_DISP);
        MPI_Type_free(&MPI_GRAND);

        return MPI_stat;
    }
#endif
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
    double parallel;            ///< \brief additional epsilon directional interactions parallel(>0), isotpropic(0), or antiparallel(<0)

    bool exclude;
    
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
        parallel = 0.0;
        exclude = false;
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




class MpiExchangeData {        // extra type for mpi communication
public:
    MpiExchangeData() {
        wantedTemp = -1;
        mpiRank = -1;
        pseudoMpiRank = -1;
        for (int wli=0;wli<2;wli++) {
            wl_order[wli] = 0;
        }
    }

    Vector box;         // box of configuration
    Vector syscm;       // system CM of configuration
    double energy;	    // energy of configuration
    double volume;      // volume of configuration
    double temperature;
    double wantedTemp;
    long radiusholemax; // size of array for WL
    long wl_order[2];   // wang-landau order parameter
    int accepted;       // BOOL if accepted
    int partNum[MAXMT];        // Number of particles
    int pseudoMpiRank;  // Pseudo-rank of the process - basically same as temperature
    int mpiRank;        // mpiRank of this process


#ifdef ENABLE_MPI
    MPI_Datatype* defDataType(MPI_Datatype* MPI_exchange, MPI_Datatype* MPI_vector2) {
        Vector vec;

        MPI_Aint     dispstart;

        MPI_Datatype vecType[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
        int           vecLen[3] = {1,          1,          1};

        MPI_Aint     vecData[3];
        MPI_Address( &vec, &dispstart);
        MPI_Address( &(vec.x), &vecData[0]);
        MPI_Address( &(vec.y), &vecData[1]);
        MPI_Address( &(vec.z), &vecData[2]);

        for (int i=0; i <3; i++)
            vecData[i] -= dispstart;

        MPI_Type_struct( 3, vecLen, vecData, vecType, MPI_vector2);
        MPI_Type_commit( MPI_vector2);

        //
        // Prepare MpiExchangeData structure in MPI
        //
        //                                box           syscm         energy      volume      temperature wantedtemp  radius..  wl_order, acc,     partNum, pseu,    mpiRank
        MPI_Datatype mpiexdataType[12] = {*MPI_vector2, *MPI_vector2, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_LONG, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
        int          mpiexdataLen[12]  = {1,            1,            1,          1,          1,          1,          1,        2,        1,       MAXMT,       1,       1};

        MPI_Aint     mpiexdata[12];
        MPI_Address( this, &dispstart);
        MPI_Address( &(this->box), &mpiexdata[0]);
        MPI_Address( &(this->syscm), &mpiexdata[1]);

        MPI_Address( &(this->energy), &mpiexdata[2]);
        MPI_Address( &(this->volume), &mpiexdata[3]);
        MPI_Address( &(this->temperature), &mpiexdata[4]);
        MPI_Address( &(this->wantedTemp), &mpiexdata[5]);

        MPI_Address( &(this->radiusholemax), &mpiexdata[6]);
        MPI_Address( &(this->wl_order), &mpiexdata[7]);

        MPI_Address( &(this->accepted), &mpiexdata[8]);
        MPI_Address( &(this->partNum), &mpiexdata[9]);
        MPI_Address( &(this->pseudoMpiRank), &mpiexdata[10]);
        MPI_Address( &(this->mpiRank), &mpiexdata[11]);

        for (int i=0; i <12; i++)
            mpiexdata[i] -= dispstart;

        MPI_Type_struct(12, mpiexdataLen, mpiexdata, mpiexdataType, MPI_exchange);
        MPI_Type_commit( MPI_exchange);

        return MPI_exchange;
    }
#endif
};



class Neighbors {
public:
    Neighbors() : neighborCount(0), neighborID(NULL) {}

    long neighborCount; ///< \brief The number of neighbors (pairs eq. num_pairs)
    long * neighborID;  ///< \brief The particle indexes of the neighbors
};

#ifdef ENABLE_MPI
extern MPI_Datatype MPI_vector, MPI_Particle, MPI_exchange;
#endif

#endif // STRUCTURES_H
