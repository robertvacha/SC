/** @file totalenergycalculator.h*/

#ifndef TOTALENERGYCALCULATOR_H
#define TOTALENERGYCALCULATOR_H

#include "pairenergycalculator.h"
#include "paire.h"
#include "externalenergycalculator.h"
#include "../structures/sim.h"

using namespace std;

class TotalEnergyCalculator
{
private:
    PairEnergyCalculator pairEControl;
    PairE pairE;

    //PairE pairEControl;
    //PairEnergyCalculator pairE;

    ExternalEnergyCalculator exterE;

    bool pairListUpdate;
    Conf* conf;

public:
    TotalEnergyCalculator(Sim * sim, Conf * conf): pairEControl(PairEnergyCalculator(&conf->geo)), pairE(PairE(&conf->geo)),
                                                   //pairEControl(PairE(&conf->geo)), pairE(PairEnergyCalculator(&conf->geo)),
                                                   exterE(ExternalEnergyCalculator(&conf->geo.box)),
                                                   pairListUpdate(sim->pairlist_update), conf(conf) {}

    /************************************************************************************************/
    /*                                                                                              */
    /*  ENERGY FUNCTIONS:                                                                           */
    /*  always provide 3 functions, one utilizing energy matrix, one which doesnt for testing       */
    /*  and third for update of part of energy matrix                                               */
    /*                                                                                              */
    /*  Further diverge all function but the ones for testing on the usage of pairlist              */
    /*                                                                                              */
    /************************************************************************************************/

    /**
     * @brief Calculates energy between all pairs, generates energy matrix. Returns energy
     * @return
     */
    double allToAll(std::vector< std::vector<double> >* energyMatrix);

    /**
     * @brief allToAll - Calculation with the use of energy matrix
     * @return
     */
    double allToAll();

    /**
     * @brief allToAllBasic - Base calculation with no optimization
     * @return
     */
    double allToAllBasic();

    ////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Calculates energy between particle target and rest, generates energy matrix. Returns energy
     * @return
     */
    double oneToAll(int target, vector<double> *changes);

    /**
     * @brief Calculates energy between particle target and rest using energy matrix Returns energy
     * @return
     */
    double oneToAll(int target);

    /**
     * @brief oneToAllBasic - Base calculation with no optimization
     * @param target
     * @return
     */
    double oneToAllBasic(int target);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief mol2others - Calculates energy between Molecule mol and rest using energy matrix Returns energy
     * @param mol
     * @return
     */
    double mol2others(Molecule &mol);

    /**
     * @brief mol2others - Calculates energy between Molecule mol and rest, generates energy matrix change, Returns energy
     * @param mol
     * @param changes
     * @return
     */
    double mol2others(Molecule &mol, vector<double> *changes);

    /**
     * @brief mol2othersBasic - Base calculation with no optimization
     * @param mol
     * @return
     */
    double mol2othersBasic(Molecule &mol);

    ///////////////////////////////////////////////////////////////////////////////////////////////



    /********************************************************************************/
    /*                            LIMITED USE FUNCTIONS                             */
    /********************************************************************************/

    /**
     * @brief p2p
     * @param part1
     * @param part2
     * @return
     */
    double p2p(int part1, int part2) {
        ConList conlist = conf->pvec.getConlist(part1);
        return pairE(&conf->pvec[part1], &conf->pvec[part2], &conlist);
    }

    double p2p(Particle* part1, int part2) {
        ConList conlist = conf->pvec.getConlist(part2);
        return pairE(part1, &conf->pvec[part2], &conlist);
    }

    /**
     * @brief mol2others - used on vector of particles which arent inserted into main pvec
     * @param mol
     * @return
     */
    double mol2others(vector<Particle> &mol);

    /**
     * @brief Calculates inner chain energy
     */
    double chainInner(vector<Particle >& chain);

    double chainInner(Molecule& chain);

    /**
     * @brief Calculates energy between particle "target" and the rest, conlist and neighbor list can be null, for GrandCanonical mainly
     */
    double oneToAll(Particle* target, ConList* conlist=NULL, Neighbors* neighborList=NULL);


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




    /********************************************************************************/
    /*                              LEGACY FUNCTIONS                                */
    /********************************************************************************/

    /**
     * @brief calc_energy Calculate the different energy contributions. Legacy function
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
              Legacy Function
     */
    double chainToAll(int target, int chainnum);


};

#endif // TOTALENERGYCALCULATOR_H
