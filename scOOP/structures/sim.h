/** @file sim.h*/

#ifndef SIM_H
#define SIM_H

#include "structures.h"
#include "../mc/wanglandau.h"

/**
 * @brief Should contain mostly all the simulation options and variables
 */
class Sim {
public:
    double press;               ///< \brief Pressure
    double paralpress;          ///< \brief Parallel pressure for replica exachnge
    double dpress;		        ///< \brief Pressure change for replica exchange
    double shave;               ///< \brief Average number of volume changes to attempt per sweep
    double shprob;              ///< \brief Probability of attempting a volume change
    double chainprob;           ///< \brief Average number of chain move attempt per sweep
    double switchprob;          ///< \brief Average number of type switch attempt per sweep
    int pairlist_update;        ///< \brief Number of sweep per updating the pairlist
    double temper;              ///< \brief Temperature
    double paraltemper;         ///< \brief Temperature for parallel tempering
    double dtemp;               ///< \brief Temprature step
    int ptype;                  ///< \brief Type of pressure coupling
    long adjust;                ///< \brief Number of sweeps between step size adjustments
    long movie;                 ///< \brief Number of sweeps between movie frames
    long nequil;                ///< \brief Number of equilibration sweeps
    long nsweeps;               ///< \brief Number of production sweeps
    long paramfrq;              ///< \brief Number of sweeps between order parameter samples
    long report;                ///< \brief Number of sweeps between statistics reports
    //long terms;                 ///< \brief Number of Fourier terms as smectic order parameters
    long nrepchange;            ///< \brief Number of sweeps between replica exchanges
    long nGrandCanon;           ///< \brief Number of sweeps between particle insert/delete
    int wlm[2];                 ///< \brief Wang landau method (wl)
    Disp edge;                  ///< \brief Maximum box length change and statistics
    Disp rot[MAXT];             ///< \brief Maximum rotation and statistics
    Disp trans[MAXT];           ///< \brief Maximum translation  and statistics
    Disp chainm[MAXMT];         ///< \brief Maximum translation for chain  and statistics
    Disp chainr[MAXMT];         ///< \brief Maximum rotation for chain and statistics
    Disp mpiexch;               ///< \brief MPI statistics
    long write_cluster;         ///< \brief Number of sweeps per writing out cluster info
    long * clusterlist;         ///< \brief clusterlist[i] = cluster index of particle i
    Cluster * clusters;         ///< \brief informations about the single clusters
    double *clustersenergy;     ///< \brief list of energies of clusters
    long num_cluster;           ///< \brief number of single clusters
    long * clusterstat;         ///< \brief Statistics about the seize of cluster
    long max_clust;             ///< \brief maximal clustersize
    WangLandau wl;              ///< \brief Wang landau data
    int mpirank;                ///< \brief MPI number for given process
    int mpinprocs;              ///< \brief MPI number of processes

    Sim() {}
};

#endif // SIM_H
