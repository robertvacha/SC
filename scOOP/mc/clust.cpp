#include "clust.h"





int Clusters::writeCluster(bool decor, long sweep) {

    FILE *cl_stat = NULL;
    FILE *cl = NULL;

    cl_stat = fopen(files->clusterstatfile, "a");
    cl = fopen(files->clusterfile, "a");

    genClusterList();
    sortClusterList();
    calcClusterEnergies();
    if(cl_stat){
        if(decor == false){
            // if no decor, this means usually into a file. Hence print info
            // about number of line per frame
            fprintf(cl_stat, "Sweep: %ld | Maximal size: %ld\n", sweep, max_clust);
        }
        printClusterStat(cl_stat, decor);
        /*
           print_clstat_oneline(cl_stat, sweep, sim);
         */
    }
    if(cl){
        if(decor == false){
            fprintf(cl, "Sweep: %ld | Number of clusters: %ld\n",
                    sweep, num_cluster);
        }
        printClusters(cl, decor);
    }
    FILE *cl_list = NULL; // CL_LIST NEVER USED
    if(cl_list){
        if(decor == false){
            fprintf(cl_list, "Sweep: %ld | Number of particles: %ld\n", sweep, (long)conf->pvec.size());
        }
        printClusterList(cl, decor);
    }

    fclose(cl_stat);
    fclose(cl);
    //fclose(cl_list);

    return 0;
}

int Clusters::sameCluster(long fst, long snd) {

    ConList conFst = conf->pvec.getConlist(fst);
    ConList conSnd = conf->pvec.getConlist(snd);

    /*if two particles are bonded they belong to the same cluster*/
    if ( ((topo.moleculeParam[conf->pvec[fst].molType]).bond1c >= 0) ||
         ((topo.moleculeParam[conf->pvec[fst].molType]).bonddc >= 0) ){
        if ( (&conf->pvec[snd] == conFst.conlist[1]) || (&conf->pvec[snd] == conFst.conlist[0]) ) {
            return true;
        }
    }
    if ( ((topo.moleculeParam[conf->pvec[snd].molType]).bond1c >= 0) ||
         ((topo.moleculeParam[conf->pvec[snd].molType]).bonddc >= 0) ){
        if ( (&conf->pvec[fst] == conSnd.conlist[1]) || (&conf->pvec[fst] == conSnd.conlist[0]) ) {
            return false;
        }
    }

    /*cluster is made of particles closer tna some distance*/
    /*	struct vector2 image(struct vector2 r1, struct vector2 r2, struct vector2 geo.box);
    struct vector2 r_cm = image(conf->pvec[fst].pos,
            conf->pvec[snd].pos,
            conf->geo.box);
    double dist2 = DOT(r_cm, r_cm);
    * TODO: Make it much more efficient => define cluster_dist!!! *
    if(dist2 > topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].sigma * topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].sigma*4.0){
        return false;
    }
    else {
        return true;
    }*/

    /*cluster is made of attractively interacting particles*/
    /*Here we add particle into cluster if particles fst and snd are from same molecule*/
    Molecule fstMol = conf->pvec.getMolOfPart(fst);
    if( std::find(fstMol.begin(), fstMol.end(), snd) != fstMol.end() ){
        return true;
    }

    if(calcEnergy->p2p(fst, snd) > -0.10 ) {
        return false;
    }
    else {
        return true;
    }
}

int Clusters::genClusterList() {
    bool change = true; /* does it still change? */
    //long neighbour;
    long i, j, fst, snd, tmp, minnumber, maxnumber;

    // Set clusterindex to the corresponding index
    for( i = 0; i < (long)conf->pvec.size(); i++){
        clusterlist[i] = i;
    }

    // Start determining the cluster
    while(change){
        change = false;
        for(i = 0; i < (long)conf->pvec.size(); i++){
            /*If nore pairlist go over all pairs*/
            maxnumber = (long)conf->pvec.size();
            minnumber = i ;
            if (sim->pairlist_update) {
                maxnumber = conf->neighborList[i].neighborCount;
                //maxnumber = sim->pairlist[i].num_pairs; // del after
                minnumber=0;
            }
            /* Go over pairs to see if they are in the cluster */
            for(j = minnumber; j < maxnumber; j++){
                fst = i;
                snd = j;
                if (sim->pairlist_update) {
                    snd = conf->neighborList[i].neighborID[j];
                    //snd = sim->pairlist[i].pairs[j];
                }
                /*do cluster analysis only for spherocylinders*/
                //if ( (topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].geotype[0] < SP) && \
                //   (topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].geotype[1] < SP) ) {

                /* if they are close to each other */
                if(sameCluster(fst, snd)){
                    if(fst > snd){
                        tmp = snd;
                        snd = fst;
                        fst = tmp;
                    }

                    if(clusterlist[fst] < clusterlist[snd]){
                        clusterlist[snd] = clusterlist[fst];
                        change = true;
                        break;
                        /* => will eventually start the i loop from new */
                    }
                    if(clusterlist[snd] < clusterlist[fst]){
                        clusterlist[fst] = clusterlist[snd];
                        change = true;
                        break;
                        /* => will eventually start the i loop from new */
                    }
                }
                //}
            }
            if(change){
                break;
            }
        }
    }

    return 0;
}

int Clusters::sortClusterList() {
    long cluster_indices[(long)conf->pvec.size()];   /* holds the different cluster indices.
                        (currently too much memory) */
    long num_cluster = 0;                /* number of clusters, temporary needed */

    /* how many clusters are there? */
    long max_index = -1;
    for(int i = 0; i < (long)conf->pvec.size(); i++){
        if(max_index < clusterlist[i]){
            max_index = clusterlist[i];
            cluster_indices[num_cluster++] = max_index;
        }
    }

    // free the memory from the old clusters */
    if(clusters){
        for(int i = 0; i < this->num_cluster; i++){
            if(clusters[i].particles){
                free(clusters[i].particles);
            }
        }
        free(clusters);
        clusters = NULL;
    }

    /* Allocate memory for the clusters */
    clusters = (Cluster*) malloc(sizeof(Cluster) * num_cluster);

    if (!clusters){
        fprintf(stderr, "Couldn't allocate any memory!\n");
        exit(1);
    }

    for(int i = 0; i < num_cluster; i++){
        /* allocate maximal space for all the clusters */
        clusters[i].particles = (long int*) malloc(sizeof(long) * (long)conf->pvec.size());

        if (!clusters[i].particles){
            fprintf(stderr, "Couldn't allocate any memory!\n");
            exit(1);
        }

        clusters[i].npart = 0;
    }

    /* fill in the particles belonging to one cluster */
    for(int i = 0; i < num_cluster; i++){
        for(int j = 0; j < (long)conf->pvec.size(); j++){
            if(clusterlist[j] == cluster_indices[i]){
                clusters[i].particles[clusters[i].npart++] = j;
            }
        }
    }
    this->num_cluster = num_cluster;

    /* Find the biggest size */
    max_clust = 0;
    for(int i = 0; i < num_cluster; i++){
        if(clusters[i].npart > max_clust){
            max_clust = clusters[i].npart;
        }
    }
    /* Set the statistics to zero */
    clusterstat = (long int*) realloc( clusterstat, sizeof(long) * max_clust); // OLD, MISTAKE? memmory dont have to be 0, no free
    memset(clusterstat, 0, sizeof(long) * max_clust);

    if (!clusterstat){
        fprintf(stderr, "Couldn't allocate any memory!\n");
        exit(1);
    }

    for(int i = 0; i < max_clust; i++) {
        clusterstat[i] = 0;
    }
    /* Do the statistics */
    for(int i = 0; i < num_cluster; i++){
        clusterstat[clusters[i].npart - 1]++;
    }

    return 0;
}
