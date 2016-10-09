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

    int printClusterList(FILE *stream, bool decor) {
        if(decor){
            fprintf(stream, "\n"
                    "-----------------------------------------------------\n"
                    "  The Cluster List\n"
                    "  (Index starts with 1)\n"
                    "-----------------------------------------------------\n");
        }

        for(int i=0; i < (long)conf->pvec.size(); i++){
            fprintf(stream,"%3d %3ld %8.4lf %8.4f %8.4f", i + 1,
                    clusterlist[i] + 1,
                    conf->pvec[i].pos.x,
                    conf->pvec[i].pos.y,
                    conf->pvec[i].pos.z);
            fprintf(stream,"\n");
        }
        if(decor){
            fprintf(stream,"-----------------------------------------------------\n");
        }
        fflush(stream);
        return 0;
    }


    int printClusters(FILE *stream, bool decor) {
        if(decor){
            fprintf(stream, "\n"
                    "-----------------------------------------------------\n"
                    "  The Clusters\n"
                    "  (Index starts with 1)\n"
                    "-----------------------------------------------------\n");
        }
        for(int i = 0; i < num_cluster; i++){
            fprintf(stream, "%3d(%f):", i + 1,clustersenergy[i]);
            for(int j = 0; j < clusters[i].npart; j++){
                fprintf(stream, "%5ld", clusters[i].particles[j] + 1);
            }
            fprintf(stream, "\n");
        }
        if(decor){
            fprintf(stream,"---------------------------------------------------\n");
        }
        fflush(stream);
        return 0;
    }


    int printClusterStat(FILE *stream, bool decor) {
        if(decor){
            fprintf(stream, "\n"
                    "-----------------------------------------------------\n"
                    "   Cluster Distribution\n"
                    "-----------------------------------------------------\n");
        }
        for(int i=0; i < max_clust; i++){
            fprintf(stream, "%5d\t%5ld\n", i + 1, clusterstat[i]);
        }
        if(decor){
            fprintf(stream, "--------------------------------------------------\n");
        }
        fflush(stream);
        return 0;
    }
};


#endif // CLUST_H
