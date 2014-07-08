/** @file sim.h*/

#ifndef SIM_H
#define SIM_H

#include "structures.h"
#include "../mc/wanglandau.h"

typedef struct{                 /* Should contain mostly all the simulation options and variables. */
    double press;               /* Pressure */
    double paralpress;          /* Parallel pressure for replica exachnge*/
    double dpress;		        /* Pressure change for replica exchange*/
    double shave;               /* Average number of volume changes to attempt per sweep */
    double shprob;              /* Probability of attempting a volume change */
    double chainprob;           /* Average number of chain move attempt per sweep */
    double switchprob;          /* Average number of type switch attempt per sweep */
    int pairlist_update;        /* Number of sweep per updating the pairlist */
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
    long nGrandCanon;           // Number of sweeps between particle insert/delete
    int wlm[2];                 /* Wang landau method (wl) */
    Disp edge;           /* Maximum box length change and statistics */
    Disp rot[MAXT];      /* Maximum rotation and statistics */
    Disp trans[MAXT];    /* Maximum translation  and statistics*/
    Disp chainm[MAXMT];  /* Maximum translation for chain  and statistics*/
    Disp chainr[MAXMT];  /* Maximum rotation for chain and statistics */
    Disp mpiexch;        /* MPI statistics*/
    long write_cluster;         /* Number of sweeps per writing out cluster info */
    long * clusterlist;         /* clusterlist[i] = cluster index of particle i */
    Cluster * clusters;  /* informations about the single clusters */
    double *clustersenergy;     /* list of energies of clusters*/
    long num_cluster;           /* number of single clusters */
    long * clusterstat;         /* Statistics about the seize of cluster */
    long max_clust;             /* maximal clustersize */
    WangLandau wl;              /* Wang landau data */
    int mpirank;                /* MPI number for given process*/
    int mpinprocs;              /* MPI number of processes */

} Sim;

#endif // SIM_H
