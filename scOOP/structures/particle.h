#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vector.h"
#include "structures.h"

/**
 * @brief Define a particle
 */
class Particle {
public:  
    /** @brief Position vector */
    Vector pos;

    /** @brief Unit direction vector of axis */
    Vector dir;

    /** @brief Vector defining orientation of patch */
    Vector patchdir[2];

    /** @brief Vector defining sides of patch */
    Vector patchsides[4];

    /** @brief Direction for chirality - keep in memory to increase speed */
    Vector chdir[2];



    /** @brief Chain type 0-100, given sequentialy from 0, effectively means molecule type */
    long chaint;

    /** @brief Chain number */
    long chainIndex;

    /** @brief Type of the particle 0 - 40 */
    int type;

    /** @brief 0: in initial stat; 1: in the switched stat */
    int switched;


    /** @brief With which kind of particle do you want to switch? */
    int switchtype;

    /** @brief Chemical potential for the switch */
    double delta_mu;


    /** @brief The number of neighbors (pairs eq. num_pairs) */
    long neighborCount;

    /** @brief The particle numbers of the neighbors */
    long * neighborID;

    /** @brief Connectivity list, we have connection to tail and head and secon neighbours so far */
    long conlist[4];


public:
    /**
     * @brief int_partvec initiate vectors of a single particle
     * @param target
     * @param ia_parami
     */
    void init(Ia_param * ia_parami);
};

#endif // PARTICLE_H
