#include "totalenergycalculator.h"

#include <stdlib.h>
#include <iostream>

extern Topo topo;


double TotalEnergyCalculator::operator ()(int target, int mode, int chainnum) {

    //DEBUG_SIM("Calculate the energy with mode %d", mode)

    if(mode == 1)
        return oneToAll(target);

    if(mode == 2)
        return chainToAll(target, chainnum);

    if(mode == 0)
        return allToAll();

    fprintf(stderr, "ERROR: Wrong mode (%d) was given to calc_energy!", mode);
    return 0.0;
}


double TotalEnergyCalculator::chainToAll(int target, int chainnum) {

    double energy=0.0;
    long i,j=0;

#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
    for (i = 0; i < target; i++) {
        if (i != conf->pvecGroupList.getChain(chainnum, j)) {
            energy+= pairE[getThreadNum()](&conf->pvec[target], &conf->conlist[target], &conf->pvec[i], &conf->conlist[i]);
        } else {
            j++;
        }
    }
    j++;

#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
    for (i = target + 1; i < (long)conf->pvec.size(); i++) {
        if (i != conf->pvecGroupList.getChain(chainnum, j)) {
            energy+= pairE[getThreadNum()](&conf->pvec[target], &conf->conlist[target], &conf->pvec[i], &conf->conlist[i]);
        } else {
            j++;
        }
    }

    //add interaction with external potential
    if (topo.exter.exist)
        energy+= extere2(target);

    return energy;
}


double TotalEnergyCalculator::chainToAll(vector<Particle>::iterator chain,
                                         vector<ConList>::iterator con, int size) {
    double energy=0.0;
    bool notChain = true;
    int i;

#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
    for(i=0; i<size; i++) {
        if (topo.exter.exist)
            energy += exterE[getThreadNum()].extere2(&*(chain+i));

        for(unsigned int j=0; j < conf->pvec.size(); j++) {
            notChain=true;
            for(int c=0; c<size; c++)
                if(&conf->pvec[j] == &*(chain+c) )
                    notChain=false;

            if(notChain)
                energy += pairE[getThreadNum()](&*(chain+i), &*(con+i), &conf->pvec[j], &conf->conlist[j]);
        }
    }

    return energy;
}


double TotalEnergyCalculator::oneToAll(int target) {
    double energy=0.0;
    long i;

    if (sim->pairlist_update) {
#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
        for (i = 0; i < conf->neighborList[target].neighborCount; i++){
           energy += (pairE[getThreadNum()])(&conf->pvec[target],
                             &conf->conlist[target],
                             &conf->pvec[ conf->neighborList[target].neighborID[i] ],
                             &conf->conlist[ conf->neighborList[target].neighborID[i] ]);
        }
    } else { // no neighborList
#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
        for (i = 0; i < (long)conf->pvec.size(); i++) {
            if(target != i) {
                energy += (pairE[getThreadNum()])(&conf->pvec[target], &conf->conlist[target], &conf->pvec[i], &conf->conlist[i]);
            }
        }
    }
    //add interaction with external potential
    if (topo.exter.exist)
        energy += extere2(target);

    return energy;
}


double TotalEnergyCalculator::oneToAll(Particle *target, ConList* conlist, Neighbors* neighborList) {
    double energy=0.0;
    unsigned long i;

    if (sim->pairlist_update) {
#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
        for (i = 0; i < (unsigned long)neighborList->neighborCount; i++){
           energy += (pairE[getThreadNum()])(target, conlist,
                                             &conf->pvec[ neighborList->neighborID[i] ],
                                             &conf->conlist[ neighborList->neighborID[i] ]);
        }
    } else {
#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
        for (i = 0; i < conf->pvec.size(); i++) {
            if(*target != conf->pvec[i]) {
                energy += (pairE[getThreadNum()])(target, conlist, &conf->pvec[i], &conf->conlist[i]);
            }
        }
    }

    //add interaction with external potential
    if (topo.exter.exist)
        energy += exterE[getThreadNum()].extere2(target);

    return energy;
}

double TotalEnergyCalculator::chainInner(vector<Particle >::iterator chain,
                                         vector<ConList>::iterator con, int size) {
    double energy = 0.0;

    for (int i=0; i<size-1; i++) {
        for (int j=i+1; j<size; j++)
            energy += (pairE[getThreadNum()])(&*(chain+i), &*(con+i), &*(chain+j), &*(con+j));

        //for every particle add interaction with external potential
        if (topo.exter.exist)
            energy += exterE[getThreadNum()].extere2(&*(chain+i));
    }

    //add interaction of last particle with external potential
    if (topo.exter.exist)
        energy+= exterE[getThreadNum()].extere2(&*(chain+size-1));

    return energy;
}

double TotalEnergyCalculator::allToAll() {
    double energy=0.0;
    unsigned long i;

#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
    for (i = 0; i < conf->pvec.size() - 1; i++) {

        for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
            energy += (pairE[getThreadNum()])(&conf->pvec[i], &conf->conlist[i], &conf->pvec[j], &conf->conlist[j]);
        }

        //for every particle add interaction with external potential
        if (topo.exter.exist)
            energy += extere2(i);
    }

    //add interaction of last particle with external potential
    if (topo.exter.exist)
        energy+= extere2(conf->pvec.size() - 1);

    return energy;
}
