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
    ConList conlist = conf->pvec.getConlist(target);

#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
    for (i = 0; i < target; i++) {
        if (i != conf->pvec.getChainPart(chainnum, j)) {
            energy+= pairE[getThreadNum()](&conf->pvec[target], &conf->pvec[i], &conlist);
        } else {
            j++;
        }
    }
    j++;

#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
    for (i = target + 1; i < (long)conf->pvec.size(); i++) {
        if (i != conf->pvec.getChainPart(chainnum, j)) {
            energy+= pairE[getThreadNum()](&conf->pvec[target], &conf->pvec[i],  &conlist);
        } else {
            j++;
        }
    }

    //add interaction with external potential
    if (topo.exter.exist)
        energy+= extere2(target);

    return energy;
}


double TotalEnergyCalculator::mol2others(vector<Particle>& mol) {
    double energy=0.0;

#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
    for(unsigned int i=0; i<mol.size(); i++) {
        if (topo.exter.exist)
            energy += exterE[getThreadNum()].extere2(&mol[i]);

        for(unsigned int j=0; j < conf->pvec.size(); j++) {
            energy += pairE[getThreadNum()](&mol[i], &conf->pvec[j]);
        }
    }

    return energy;
}

double TotalEnergyCalculator::mol2others(Molecule& mol) {
    double energy=0.0;

#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
    for(unsigned int j=0; j<mol.size(); j++) {
        for (long i = 0; i < mol[0]; i++)
            energy+= pairE[getThreadNum()](&conf->pvec[mol[j]], &conf->pvec[i]);

        for (long i = mol[mol.size()-1] + 1; i < (long)conf->pvec.size(); i++)
            energy+= pairE[getThreadNum()](&conf->pvec[mol[j]], &conf->pvec[i]);

        //add interaction with external potential
        if (topo.exter.exist)
            energy+= extere2(mol[j]);
    }

    return energy;
}


double TotalEnergyCalculator::oneToAll(int target) {
    double energy=0.0;
    long i;
    ConList conlist = conf->pvec.getConlist(target);
    if (sim->pairlist_update) {
#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
        for (i = 0; i < conf->neighborList[target].neighborCount; i++){
           energy += (pairE[getThreadNum()])(&conf->pvec[target],
                             &conf->pvec[ conf->neighborList[target].neighborID[i] ],
                             &conlist);
        }
    } else { // no neighborList
#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
        for (i = 0; i < (long)conf->pvec.size(); i++) {
            if(target != i) {
                energy += (pairE[getThreadNum()])(&conf->pvec[target], &conf->pvec[i], &conlist);
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
           energy += (pairE[getThreadNum()])(target, &conf->pvec[ neighborList->neighborID[i] ], conlist );
        }
    } else {
#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
        for (i = 0; i < conf->pvec.size(); i++) {
            if(*target != conf->pvec[i]) {
                energy += pairE[getThreadNum()](target, &conf->pvec[i], conlist);
                /*std::cout.precision(15);
                cout << pairE[getThreadNum()](target, conlist, &conf->pvec[i], &conf->conlist[i])
                        <<" == "<< pairE[getThreadNum()](&conf->pvec[i], &conf->conlist[i], target, conlist) << endl;
                assert((float)pairE[getThreadNum()](target, conlist, &conf->pvec[i], &conf->conlist[i])
                        == (float)pairE[getThreadNum()](&conf->pvec[i], &conf->conlist[i], target, conlist) );*/
            }
        }
    }

    //add interaction with external potential
    if (topo.exter.exist)
        energy += exterE[getThreadNum()].extere2(target);

    return energy;
}

double TotalEnergyCalculator::chainInner(vector<Particle >::iterator chain,
                                         int size, vector<ConList>::iterator con) {
    double energy = 0.0;

    for (int i=0; i<size-1; i++) {
        for (int j=i+1; j<size; j++)
            energy += (pairE[getThreadNum()])(&*(chain+i), &*(chain+j),  &*(con+i));

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
    if(conf->pvec.empty()) return 0.0;
    double energy=0.0;
    unsigned long i;
    ConList conlist;

#ifdef OMP1
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
    for (i = 0; i < conf->pvec.size() - 1; i++) {
        conlist = conf->pvec.getConlist(i);
        for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
            energy += (pairE[getThreadNum()])(&conf->pvec[i], &conf->pvec[j], &conlist);
        }

        //for every particle add interaction with external potential
        if (topo.exter.exist)
            energy += extere2(i);
    }
    //add interaction of last particle with external potential
    if (topo.exter.exist && !conf->pvec.empty())
        exterE[getThreadNum()].extere2(&conf->pvec.back());
    return energy;
}
