/** @file totalenergycalculator.h*/

#ifndef TOTALENERGYCALCULATOR_H
#define TOTALENERGYCALCULATOR_H

#include <stdlib.h> // exit()

#include <iostream>

#include "pairenergycalculator.h"
#include "externalenergycalculator.h"
#include "../structures/sim.h"

using namespace std;

class TotalEnergyCalculator
{
public:
    PairEnergyCalculator pairE;
    ExternalEnergyCalculator exterE;

private:
    Topo* topo;
    Sim* sim;
    Conf* conf;

public:
    TotalEnergyCalculator(Topo * topo, Sim * sim, Conf * conf):
        pairE(topo, &conf->box), exterE(&topo->exter, &conf->box),
        topo(topo), sim(sim), conf(conf) {pairE.initIntFCE();}

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
     * @brief Calculates energy between particle "target" and the rest skipping particles from the given chain
              -particles has to be sorted in chain!!
              similar to oneToAll, but with chain exception
     */
    double chainToAll(int target, int chainnum);

    double chainToAll(vector<Particle>::iterator chain, vector<ConList>::iterator con, int size);

    /**
     * @brief Calculates energy between particle "target" and the rest
     */
    double oneToAll(int target);

    /**
     * @brief Calculates energy between particle "target" and the rest
     */
    double oneToAll(Particle* target, ConList* conlist);

    /**
     * @brief Calculates energy between particle "target" and the rest
     */
    double chainInner(vector<Particle >::iterator chain, vector<ConList>::iterator con, int size);

    /**
     * @brief Calculates energy between all pairs. Returns energy
     * @return
     */
    double allToAll();

    /**
     * @brief extere2    Calculates interaction of target particle and external field version 2
                         calculate projection of spherocylinder in direction of patch and calculate
                         interacting line segment within cutoff
     * @param target
     * @return
     */
    double extere2(int target) {
        return exterE.extere2(&conf->pvec[target]);
    }
};

#endif // TOTALENERGYCALCULATOR_H
