/** @file totalenergycalculator.h*/

#ifndef TOTALENERGYCALCULATOR_H
#define TOTALENERGYCALCULATOR_H

#include "pairenergycalculator.h"
#include "externalenergycalculator.h"
#include "../structures/sim.h"

using namespace std;

class TotalEnergyCalculator
{
private:
    vector<PairEnergyCalculator> pairE;
    vector<ExternalEnergyCalculator> exterE;

#ifdef OMP
    inline int getThreadNum() {
        return omp_get_thread_num();
    }
#else
    inline int getThreadNum() {return 0;}
#endif

    Sim* sim;
    Conf* conf;

public:
    TotalEnergyCalculator(Sim * sim, Conf * conf): sim(sim), conf(conf) {
        int threadCount = 1;
#ifdef OMP
        threadCount = 8;
        omp_set_num_threads(threadCount);
        cout << "Number of threads: " << threadCount << endl;
#endif

        for(int i=0; i < threadCount; i++) {
            pairE.push_back(PairEnergyCalculator(&conf->geo));
            exterE.push_back(ExternalEnergyCalculator(&conf->geo.box));
            pairE[i].initIntFCE();
        }
    }

    /**
     * @brief calc_energy Calculate the different energy contributions.
     *
     * This is a merge of the different energy calculation functions (energyone, -chain, -all)
     * 0: all
     * 1: one
     * 2: chain
     *
     * @param target
     * @param mode
     * @param chainnum
     * @return
     */
    double operator() (int target, int mode, int chainnum);

    /**
     * @brief Calculates energy between all pairs. Returns energy
     * @return
     */
    double allToAll();

    /**
     * @brief p2p
     * @param part1
     * @param part2
     * @param conlist of first particle
     * @return
     */
    double p2p(int part1, int part2) {
        ConList conlist = conf->pvec.getConlist(part1);
        return pairE[getThreadNum()](&conf->pvec[part1], &conf->pvec[part2], &conlist);
    }

    /**
     * @brief Calculates energy between particle "target" and the rest skipping particles from the given chain
              -particles has to be sorted in chain!!
              similar to oneToAll, but with chain exception
     */
    double chainToAll(int target, int chainnum);

    double mol2others(vector<Particle> &mol);

    double mol2others(Molecule &mol);

    /**
     * @brief Calculates energy between particle "target" and the rest
     */
    double oneToAll(int target);

    /**
     * @brief Calculates energy between particle "target" and the rest
     */
    double oneToAll(Particle* target, ConList* conlist, Neighbors* neighborList);

    /**
     * @brief Calculates inner chain energy
     */
    double chainInner(vector<Particle >& chain, vector<ConList>& con);

    /**
     * @brief extere2    Calculates interaction of target particle and external field version 2
                         calculate projection of spherocylinder in direction of patch and calculate
                         interacting line segment within cutoff
     * @param target
     * @return
     */
    double extere2(int target) {
        return exterE[getThreadNum()].extere2(&conf->pvec[target]);
    }
};

#endif // TOTALENERGYCALCULATOR_H
