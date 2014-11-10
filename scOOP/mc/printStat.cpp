#include "printStat.h"


int printStat::printClusterList(FILE *stream, bool decor, Sim *sim, Conf *conf) {
    if(decor){
        fprintf(stream, "\n"
                "-----------------------------------------------------\n"
                "  The Cluster List\n"
                "  (Index starts with 1)\n"
                "-----------------------------------------------------\n");
    }

    for(int i=0; i < (long)conf->pvec.size(); i++){
        fprintf(stream,"%3d %3ld %8.4lf %8.4f %8.4f", i + 1,
                sim->clusterlist[i] + 1,
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


int printStat::printClusters(FILE *stream, bool decor, Sim *sim) {
    if(decor){
        fprintf(stream, "\n"
                "-----------------------------------------------------\n"
                "  The Clusters\n"
                "  (Index starts with 1)\n"
                "-----------------------------------------------------\n");
    }
    for(int i = 0; i < sim->num_cluster; i++){
        fprintf(stream, "%3d(%f):", i + 1,sim->clustersenergy[i]);
        for(int j = 0; j < sim->clusters[i].npart; j++){
            fprintf(stream, "%5ld", sim->clusters[i].particles[j] + 1);
        }
        fprintf(stream, "\n");
    }
    if(decor){
        fprintf(stream,"---------------------------------------------------\n");
    }
    fflush(stream);
    return 0;
}


int printStat::printClusterStat(FILE *stream, bool decor, Sim *sim) {
    if(decor){
        fprintf(stream, "\n"
                "-----------------------------------------------------\n"
                "   Cluster Distribution\n"
                "-----------------------------------------------------\n");
    }
    for(int i=0; i < sim->max_clust; i++){
        fprintf(stream, "%5d\t%5ld\n", i + 1, sim->clusterstat[i]);
    }
    if(decor){
        fprintf(stream, "--------------------------------------------------\n");
    }
    fflush(stream);
    return 0;
}


void printStat::printPairList(FILE *stream, Conf *conf) {
    for (unsigned long i = 0; i < conf->neighborList.size(); i++){
        fprintf(stream, "%ld (%ld):", i, conf->neighborList[i].neighborCount);
        for(long j = 0; j < conf->neighborList[i].neighborCount; j++){
            fprintf(stream, " %ld", conf->neighborList[i].neighborID[j]);
        }
        fprintf(stream, "\n");
    }
}
