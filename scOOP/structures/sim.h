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
    long nClustMove;            ///< \brief Number of sweeps between cluster moves

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
    ~Sim() {
        printf ("Deallocating Sim...\n");

        if (clusterlist != NULL)
            free(clusterlist);

        if (clustersenergy != NULL)
            free(clustersenergy);

        /*if (pairlist_update) {
            if(deallocPairlist()) {
                return 1;
            }
        }*/
    }

    void printEqStat() {

        printf ("   Equilibrated maximum displacement / acceptance ratio:            \n");
        printEqStat(trans,2.0,MAXT);

        printf ("   Equilibrated maximum rotation / acceptance ratio:                       \n");
        printEqStat(rot,1.0,MAXT);

        printf ("   Equilibrated maximum box length change / acceptance ratio:              \n");
        printf ("                     %.6e  /  %.6e\n", edge.mx/2.0,RATIO(edge));

        printf ("   Equilibrated maximum displacement of chain / acceptance ratio:   \n");
        printEqStat(chainm,2.0,MAXMT);

        printf ("   Equilibrated maximum rotation of chain / acceptance ratio:              \n");
        printEqStat(chainr,1.0,MAXMT);
        printf ("\n");
    }

    void printEqStat(Disp *dat, double scale, int length) {
        for(int i=0; i<length; i++) {
            if (RATIO(dat[i]) > 0)
                printf ("   TYPE %d           %.6f  /  %.6f\n", i, dat[i].mx/scale,RATIO(dat[i]));
        }
    }

    void info() {
        printf (" Pressure coupling type:                             %d\n", ptype);
        printf (" Pressure:                                           %.8f\n", press);
        printf (" Replica exchange pressure:                          %.8f\n", paralpress);
        printf (" Average volume change attempts per sweep:           %.8f\n", shave);
        printf (" Equilibration sweeps:                               %ld\n", nequil);
        printf (" Sweeps between step size adjustments:               %ld\n", adjust);
        printf (" Production sweeps:                                  %ld\n", nsweeps);
        printf (" Sweeps between statistics samples:                  %ld\n", paramfrq);
        printf (" Sweeps between statistics reports:                  %ld\n", report);
        printf (" Average chain move attempts per sweep:              %.8f\n", chainprob);
        printf (" Inititial maximum box edge change:                  %.8f\n", edge.mx);
        printf (" Temperature in kT/e:                                %.8f\n", temper);
        printf (" Parallel tempering temperature in kT/e:             %.8f\n", paraltemper);
        printf (" Sweeps between replica exchange:                    %ld\n", nrepchange);
        printf (" Sweeps between Grand-Canonical move:                %ld\n", nGrandCanon);
        printf (" Sweeps between Cluster moves:                       %ld\n", nClustMove);
        printf (" Wang-Landau method:                                 %d %d\n", wl.wlm[0],wl.wlm[1]);
        printf (" Calculate the Wang-Landau method for atom type:     %d\n", wl.wlmtype);
        printf (" Average type switch attempts per sweep:             %.8f\n", switchprob);
        printf (" Number of Sweeps per pairlist update:               %d\n", pairlist_update);
        printf (" Number of sweeps per writing out cluster info:      %ld\n", write_cluster);
    }

private:
    int deallocPairlist() {  // deprecated, done in memoryDealloc()
        /*long i;
        if(sim.pairlist != NULL){
            for(i = 0; i < topo.npart; i++){
                if(sim.pairlist[i].pairs != NULL){
                    free(sim.pairlist[i].pairs);
                }
            }
            free(sim.pairlist);
        }*/
        return 0;
    }
};

#endif // SIM_H
