/** @file sampler.h*/

#ifndef SAMPLER_H
#define SAMPLER_H

#include "../structures/sim.h"

//
namespace printStat {

    /**
     * @brief print_clusterlist print the clusterlist
     * @param stream
     * @param decor
     * @param topo
     * @param sim
     * @param conf
     * @return
     */
    int printClusterList(FILE * stream, bool decor, Sim * sim, Conf * conf);

    /**
     * @brief print_clusters print the clusters
     * @param stream
     * @param decor
     * @param sim
     * @return
     */
    int printClusters(FILE * stream, bool decor, Sim * sim);

    /**
     * @brief print_clusterstat print a statistics for the clusters
     * @param stream
     * @param decor
     * @param sim
     * @return
     */
    int printClusterStat(FILE * stream, bool decor, Sim * sim);

    /**
     * @brief print_pairlist Print out the pairlist
     * @param stream
     * @param sim
     * @param topo
     */
    void printPairList(FILE * stream, Conf *conf);

    /**
     * @brief printeqstat
     * @param dat
     * @param scale
     * @param length
     */
    void printEqStat(Disp *dat, double scale, int length);

    /**
     * @brief draw Dumps a configuration to the supplied file handle.
     * @param outfile
     * @param conf
     * @param topo
     */
    void draw(FILE *outfile, Conf * conf);

}


#endif // SAMPLER_H
