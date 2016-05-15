#ifndef CLUST_H
#define CLUST_H

#include <algorithm>

#include "../structures/Conf.h"
#include "../structures/sim.h"
#include "totalenergycalculator.h"
#include "../structures/structures.h"

/**
 * @brief contains all the particles of one cluster
 */
typedef struct{
    long npart;
    long * particles;
} Cluster;

class Clusters {
public:
    Conf* conf;
    Sim* sim;
    TotalEnergyCalculator* calcEnergy;
    FileNames* files;

    long * clusterstat;         ///< \brief Statistics about the size of cluster

    Cluster * clusters;         ///< \brief informations about the single clusters

    long * clusterlist;         ///< \brief clusterlist[i] = cluster index of particle i
    double *clustersenergy;     ///< \brief list of energies of clusters
    long num_cluster;           ///< \brief number of single clusters
    long max_clust;             ///< \brief maximal clustersize

    Clusters(Conf* conf, Sim* sim, TotalEnergyCalculator* calcEnergy, FileNames* files) : conf(conf),
    sim(sim), calcEnergy(calcEnergy), files(files) {

        clusterstat = (long int*) malloc(sizeof(long) * max_clust);
        clusterlist = (long int*) malloc(sizeof(long) * MAXN);
        if(clusterlist == NULL){
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for clusterlist!");
            exit(1);
        }
        clustersenergy = (double*) malloc(sizeof(double) * MAXN);
        if(clustersenergy== NULL){
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sim->clustersenergy!");
            exit(1);
        }
        clusters = NULL;
    }

    ~Clusters() {
        if (clusterstat != NULL)
            free(clusterstat);

        if (clusterlist != NULL)
            free(clusterlist);

        if (clustersenergy != NULL)
            free(clustersenergy);
    }

    /**
     * @brief write_cluster write out all the cluster stat in files, if file name is given
     * @param cl_stat
     * @param cl
     * @param cl_list
     * @param decor
     * @param sweep
     * @return
     */
    int writeCluster(bool decor, long sweep);

private:

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
    int calcClusterEnergies() {
        for(int i = 0; i < num_cluster; i++) {
            clustersenergy[i]=0.0;
            for(int j = 0; j < clusters[i].npart; j++) {
                for(int k = j+1; k < clusters[i].npart; k++) {
                    clustersenergy[i]+= calcEnergy->p2p(clusters[i].particles[j], // particle 1
                                                       clusters[i].particles[k]);        // particle 2
                }
            }
        }
        return 0;
    }

    int printClusterList(FILE *stream, bool decor);


    int printClusters(FILE *stream, bool decor);


    int printClusterStat(FILE *stream, bool decor);
};


#endif // CLUST_H
