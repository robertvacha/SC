/** @file structures.h*/

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "macros.h"
#include "../mc/math_calc.h"
#include <vector>
#include <cstdlib>
#include "moleculeparams.h"

#include <cstdio>
#include <array>


class Molecule : public vector<int > {
public:
    std::string info(){
        std::ostringstream o;
        o << "[";
        for ( unsigned int i = 0; i < this->size()-1; i++ ){
            o << this->operator [](i) << ", ";
        }
        o << this->operator [](this->size()-1) << "]";
        return o.str();
    }
};

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
    Double,
    Tuple
} Type;




typedef struct {        // for reading in the options
    const char *id;     // The name of the value in the option file
    Type type;          // The type (int, double or long)
    bool set;           // Whether the variable has been set
    void *var;          // The variable
} Option;




class FileNames{
public:
    FileNames(int rank=0){
        sprintf(configurationPool, "pool");
        sprintf(configurationInFile, "config.init");
        sprintf(configurationoutfile, "config.last");
        sprintf(optionsfile, "options");
        sprintf(topologyInFile, "top.init");
        sprintf(topologyOutFile, "top.last");
        sprintf(moviefile, "movie");
        sprintf(wlinfile, "wl.dat");
        sprintf(wloutfile, "wl-new.dat");
        sprintf(statfile, "stat.dat");
        sprintf(clusterfile, "cluster.dat");
        sprintf(clusterstatfile, "cluster_stat.dat");
        sprintf(energyfile, "energy.dat");

#ifdef ENABLE_MPI
        initMPI(rank);
#endif
    }

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

    void initMPI(int rank) {
        FILE *infile;
        initMPIRank(rank);

        //test if there is a specific input configuration for mpi run

        infile = fopen(configurationInFile, "r");
        if (infile != NULL)
            fclose (infile);
        else  sprintf(configurationInFile, "config.init");

        infile = fopen(topologyInFile, "r");
        if (infile != NULL)
            fclose (infile);
        else  sprintf(topologyInFile, "top.init");

        //test if there is a specific input wang-landau for mpi run

        infile = fopen(wlinfile, "r");
        if (infile != NULL)
            fclose (infile);
        else  sprintf(wlinfile, "wl.dat");
    }

    void initMPIRank(int rank) {
        // in files
        sprintf(configurationPool, "%dpool", rank);
        sprintf(configurationInFile, "%dconfig.init", rank);
        sprintf(topologyInFile, "%dtop.init", rank);

        // topology in file -> change only for grand  canonical paralel tempering
        // options file is the same
        sprintf(wlinfile, "%dwl.dat", rank);

        // out files
        sprintf(configurationoutfile, "%dconfig.last", rank);
        sprintf(topologyOutFile, "%dtop.last", rank);

        // topology out file -> change only for grand  canonical paralel tempering
        sprintf(moviefile, "%dmovie", rank);

        if(false /* Multiple walkers Wang-Landau */) {
            sprintf(wloutfile, "%dwl-new.dat", rank);
        }

        sprintf(statfile, "%dstat.dat", rank);
        sprintf(clusterfile, "%dcluster.dat", rank);
        sprintf(clusterstatfile, "%dcluster_stat.dat", rank);
        sprintf(energyfile, "%denergy.dat", rank);
    }
};




/**
  * @brief Contatins properties and parameters of particle types
  **/
class Ia_param{
public:
    int type=-1;
    string name;           ///< \brief The name of the particle type
    string other_name;     ///< \brief The name of the particle type
     std::array<int, 2> geotype = {};   ///< \brief The geometrical type:
                                        /// spherocylinder (0-repulsive, 1-isotropic, 2-patchy, 3-cylindrical)
                                        /// or sphere (0-repulsive, 1-isotropic)
    double sigma        = 0.0;  ///< \brief Repulsion wca
    double sigmaSq      = 0.0;  ///< \brief Repulsion wca
    double epsilon      = 0.0;  ///< \brief Repulsion strength
    double A            = 0.0;  ///< \brief A = 4 * epsilon * sigma^12
    double B            = 0.0;  ///< \brief A = 4 * epsilon * sigma^6

    double pdis         = 0.0;  ///< \brief Interaction distance of patch
    double pdisSq       = 0.0;  ///< \brief Interaction distance of patch squared
    double pswitch      = 0.0;  ///< \brief Switch of distance of patch
    double pswitchINV   = 0.0;  ///< \brief Inverted Switch of distance of patch

    double rcut         = 0.0;  ///< \brief Cutoff for attraction
    double rcutSq       = 0.0;  ///< \brief Cutoff for attraction squared
    double rcutwca      = 0.0;  ///< \brief Cutoff for repulsion
    double rcutwcaSq    = 0.0;  ///< \brief Cutoff for repulsion squared

    std::array<double, 4> pcoshalfi     = {};   ///< \brief Cosine of half angle going to side of interaction
    std::array<double, 4> psinhalfi     = {};   ///< \brief Sine of half angle going to side of interaction -useful for quaterion rotation
    std::array<double, 2> csecpatchrot  = {};   ///< \brief Cosine of Rotation of second patches in 2psc models
    std::array<double, 2> ssecpatchrot  = {};   ///< \brief Sine of Rotation of second patches in 2psc models

    std::array<double, 4> pangl         = {};   ///< \brief angular size of patch as was specifid in input
    std::array<double, 4> panglsw       = {};   ///< \brief angular size of patchswitch as was specifid in input
    std::array<double, 4> pcangl        = {};   ///< \brief cosine of half size angle - rotation from patch direction to side
    std::array<double, 4> pcanglsw      = {};   ///< \brief cosine of half size angle plus switch - rotation from patch direction to side

    std::array<double, 2> len           = {};   ///< \brief Length of the PSC
    std::array<double, 2> half_len      = {};   ///< \brief Half length of the PSC
    std::array<double, 2> chiral_cos    = {};   ///< \brief Coctains the cosinus for the chiral rotation of the patch
    std::array<double, 2> chiral_sin    = {};   ///< \brief Contains the sinus for the chiral rotation of the patch

    double patchRot = 0.0;
    double chiral = 0.0;

    double volume = 0.0;              ///< \brief Volume of particle for geometrical center calculations
    double pvolscale = 0.0;           ///< \brief Scale of patch volume size
    double parallel = 0.0;            ///< \brief additional epsilon directional interactions parallel(>0), isotpropic(0), or antiparallel(<0)
    bool exclude = false;

    bool operator==(const Ia_param& o) const {
        return (this->geotype == o.geotype) && (this->sigma == o.sigma)
                && (this->sigmaSq == o.sigmaSq) && (this->sigmaSq == o.sigmaSq) && (this->epsilon == o.epsilon)
                && (this->A == o.A) && (this->B == o.B) && (this->pdis == o.pdis) && (this->pdisSq == o.pdisSq)
                && (this->pswitch == o.pswitch) && (this->pswitchINV == o.pswitchINV) && (this->pangl == o.pangl)
                && (this->panglsw == o.panglsw) && (this->pcangl == o.pcangl) && (this->pcanglsw == o.pcanglsw)
                && (this->rcut == o.rcut) && (this->rcutSq == o.rcutSq) && (this->rcutwca == o.rcutwca) && (this->rcutwcaSq == o.rcutwcaSq)
                && (this->pcoshalfi == o.pcoshalfi) && (this->psinhalfi == o.psinhalfi) && (this->csecpatchrot == o.csecpatchrot)
                && (this->ssecpatchrot == o.ssecpatchrot) && (this->volume == o.volume) && (this->pvolscale == o.pvolscale)
                && (this->len == o.len) && (this->half_len == o.half_len) && (this->chiral_cos == o.chiral_cos) && (this->chiral_sin == o.chiral_sin)
                && (this->parallel == o.parallel) && (this->exclude == o.exclude);
    }

    string convertGeotype(int type) const {
        if(type == SCN)     return "SCN";
        if(type == SCA)     return "SCA";
        if(type == PSC)     return "PSC";
        if(type == CPSC)    return "CPSC";
        if(type == CHPSC)   return "CHPSC";
        if(type == CHCPSC)  return "CHCPSC";
        if(type == TPSC)    return "TPSC";
        if(type == TCPSC)   return "TCPSC";
        if(type == TCHPSC)  return "TCHPSC";
        if(type == TCHCPSC) return "TCHCPSC";
        if(type == SPN)     return "SPN";
        if(type == SPA)     return "SPA";
        return string();
    }

    string toString() {
        stringstream ss;

        ss << name << " " << type << " " << convertGeotype(geotype[0]) << " " << epsilon << " " << sigma;
        if(geotype[0] == SPN) return ss.str();

        ss << " " << pdis << " " << pswitch;
        if(geotype[0] == SPA) return ss.str();

        ss << " " << pangl[0] << " " << panglsw[0] << " " << len[0] << " " << parallel;
        if(geotype[0] == PSC || geotype[0] == CPSC) return ss.str();
        if(geotype[0] == CHPSC || geotype[0] == CHCPSC) {
            ss << " " << chiral;
            return ss.str();
        }

        ss << " " << patchRot << " " << pangl[2] << " " << panglsw[2];
        if(geotype[0] == TPSC || geotype[0] == TCPSC) return ss.str();
        if(geotype[0] == TCHPSC || geotype[0] == TCHCPSC) {
            ss << " " << chiral;
            return ss.str();
        }
        return ss.str();
    }
};

class Exters{
public:
    bool exist          = false;    ///< \brief existence of external potential
    double thickness    = 0.0;      ///< \brief external wall thicnkess
    double epsilon      = 0.0;      ///< \brief depth of attraction
    double attraction   = 0.0;      ///< \brief distance of attraction
    double sqmaxcut     = 0.0;      ///< \brief distance when nothing can interact
    std::array<Ia_param, MAXT>  interactions;   ///< \brief Interaction parameters with particle types generated from above params

    bool operator==(Exters& o) const {
        return (this->exist == o.exist) && (this->thickness == o.thickness) && (this->epsilon == o.epsilon)
                && (this->attraction == o.attraction) && (this->sqmaxcut == o.sqmaxcut) && (this->interactions == o.interactions);
    }

    string toString() {
        stringstream ss;

        ss << "[EXTER]" << endl;
        ss << thickness << " " << epsilon << " " << attraction << endl;

        return ss.str();
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

    double press;       // pressure
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
        //                                press       energy      volume      temperature wantedtemp  radius..  wl_order, acc,     partNum, pseu,    mpiRank
        MPI_Datatype mpiexdataType[11] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_LONG, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
        int          mpiexdataLen[11]  = {1,          1,          1,          1,          1,          1,        2,        1,       MAXMT,       1,       1};

        MPI_Aint     mpiexdata[11];
        MPI_Address( this, &dispstart);
        MPI_Address( &(this->press), &mpiexdata[0]);
        MPI_Address( &(this->energy), &mpiexdata[1]);
        MPI_Address( &(this->volume), &mpiexdata[2]);
        MPI_Address( &(this->temperature), &mpiexdata[3]);
        MPI_Address( &(this->wantedTemp), &mpiexdata[4]);

        MPI_Address( &(this->radiusholemax), &mpiexdata[5]);
        MPI_Address( &(this->wl_order), &mpiexdata[6]);

        MPI_Address( &(this->accepted), &mpiexdata[7]);
        MPI_Address( &(this->partNum), &mpiexdata[8]);
        MPI_Address( &(this->pseudoMpiRank), &mpiexdata[9]);
        MPI_Address( &(this->mpiRank), &mpiexdata[10]);

        for (int i=0; i <11; i++)
            mpiexdata[i] -= dispstart;

        MPI_Type_struct(11, mpiexdataLen, mpiexdata, mpiexdataType, MPI_exchange);
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
