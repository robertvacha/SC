#ifndef CONF_H
#define CONF_H

#include "structures.h"
#include "particle.h"

/**
 * @brief Configuration of the system
 */
class Conf {
public:
    /** @brief Main store of all particles, grouped by Molecular types */
    std::vector<Particle* > particleStore;

    /** @brief Box size */
    Vector box;

    /** @brief Something like total mass -> for each += massOfParticle */
    double sysvolume;

    /** @brief System center of mass */
    Vector syscm;

    long muVTpartCount;

    Conf(){muVTpartCount =0;}

public:
    /**
     * @brief massCenter
     * @param topo
     */
    void massCenter(Topo* topo);

    /**
     * @brief partvecinit calculate vectors on particles for speedup
     */
    void partVecInit(Topo *topo);

    /**
     * @brief Determines whether two particles overlap.
     * @param part1
     * @param part2
     * @param ia_params 1 if there is an overlap, 0 if not.
     * @return
     */
    int overlap(Particle* part1, Particle* part2, Ia_param ia_params[MAXT][MAXT]);

    double linemin(double criterion, double halfl)  {
        if      (criterion >=  halfl) { return  halfl; }
        else if (criterion >= -halfl) { return  criterion; }
        else                          { return -halfl; }
    }

    /**
     * @brief forbidden Checks for overlaps between particle "target" and the rest.
     * @param target
     * @param ia_params
     * @return Returns true if overlap detected, false otherwise.
     */
    bool overlapAll(Particle* target, Ia_param ia_params[MAXT][MAXT]);

    /**
     * @brief checkall Checks for overlaps between all pairs of particles.
     * @param ia_params
     * @return  Returns 1 if overlap detected, 0 otherwise.
     */
    int checkall(Ia_param ia_params[MAXT][MAXT]);


};


#endif // CONF_H
