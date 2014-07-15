#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vector.h"
#include "structures.h"

/**
 * @brief Define a particle
 */
class Particle {
public:    
    Vector pos;             ///< \brief Position vector
    Vector dir;             ///< \brief Unit direction vector of axis
    Vector patchdir[2];     ///< \brief Vector defining orientation of patch
    Vector patchsides[4];   ///< \brief Vector defining sides of patch
    Vector chdir[2];        ///< \brief Direction for chirality - keep in memory to increase speed

    long molType;           ///< \brief Molecule type 0-100, given sequentialy from 0
    long chainIndex;        ///< \brief Chain number, only for Molecules of two or more particles
    int type;               ///< \brief Type of the particle 0 - 40
    int switched;           ///< \brief 0: in initial stat; 1: in the switched stat
    int switchtype;         ///< \brief With which kind of particle do you want to switch?
    double delta_mu;        ///< \brief Chemical potential for the switch

    long neighborCount; ///< \brief The number of neighbors (pairs eq. num_pairs)
    long * neighborID;  ///< \brief The particle numbers of the neighbors
    long conlist[4];    ///< \brief Connectivity list, we have connection to tail and head and secon neighbours so far

public:
    /**
     * @brief int_partvec initiate vectors of a single particle
     * @param target
     * @param ia_parami
     */
    void init(Ia_param * ia_parami);
};

#endif // PARTICLE_H
