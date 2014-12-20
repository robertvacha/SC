/** @file updater.h*/

#ifndef UPDATER_H
#define UPDATER_H

#include "printStat.h"
#include "movecreator.h"
#include "../structures/Conf.h"

class Updater
{
public:
    Updater(Sim* sim, Conf* conf, FileNames* files) :
        sim(sim), conf(conf), files(files),
        calcEnergy(sim, conf), move(sim, conf, &calcEnergy) {}

private:
    Sim* sim;                  // Should contain the simulation options.
    Conf* conf;                // Should contain fast changing particle and box(?) information
    FileNames* files;

    TotalEnergyCalculator calcEnergy;
    MoveCreator move;

    long nsweeps;
    long adjust;
    long paramfrq;
    long report;

public:
    /**
     * @brief simulate
     * @param nsweeps
     * @param adjust
     * @param paramfrq
     * @param report
     */
    void simulate(long nsweeps, long adjust, long paramfrq, long report);

private:
    void openFilesClusterStatistics(FILE **cl_stat, FILE **cl, FILE **cl_list, FILE **ef, FILE **statf);
    void initValues(long &next_adjust, long &next_calc, long &next_dump, long &next_frame);

    /**
     * @brief optimizestep Optimize the maximum displacement within the specified limits and resets the
       acceptance counters to zero.
     * @param x
     * @param hi
     * @param lo
     */
    void optimizeStep(Disp *x, double hi, double lo);

    /**
     * @brief optimizerot Optimize the maximum rotation within the specified limits and resets the
       acceptance counters to zero. Rotation is given by cos of angle
       larger rotation = smaller cos
     * @param x
     * @param hi
     * @param lo
     */
    void optimizeRot(Disp *x, double hi, double lo);

    /**
     * @brief alignment_order alignment ORDER PARAMETER
     * @return
     */
    double alignmentOrder();

    /**
     * @brief gen_pairlist Interface for the generation of the pairlist. Define other pairlist algorithms above.
     */
    void genPairList();

    /**
     * @brief gen_simple_pairlist Generates a pairlist with a very basic alogrithm
     */
    void genSimplePairList();


    /****************************************************************************/
    /* Cluster statistics stuf                                                  */
    /****************************************************************************/


    /**
     * @brief same_cluster determines, wheter two particles are in the same cluster
     * @param fst
     * @param snd
     * @return
     */
    int sameCluster(long fst, long snd);

    /**
     * @brief gen_clusterlist generate the clusterlist
     * @return
     */
    int genClusterList();

    /**
     * @brief sort_clusterlist sort the clusterlist
     * @return
     */
    int sortClusterList();

    /**
     * @brief calc_clusterenergies calculate energies of clusters
     * @return
     */
    int calcClusterEnergies();

    /**
     * @brief write_cluster write out all the cluster stat in files, if file name is given
     * @param cl_stat
     * @param cl
     * @param cl_list
     * @param decor
     * @param sweep
     * @return
     */
    int writeCluster(FILE * cl_stat, FILE * cl, FILE * cl_list, bool decor, long sweep);

};

#endif // UPDATER_H
