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
    //long chainIndex;        ///< \brief Chain number, only for Molecules of two or more particles
    int type;               ///< \brief Type of the particle 0 - 40
    int switchtype;         ///< \brief With which kind of particle do you want to switch?
    double delta_mu;        ///< \brief Chemical potential for the switch
    int switched;           ///< \brief 0: in initial stat; 1: in the switched stat

public:
    /**
     * @brief int_partvec initiate vectors of a single particle
     * @param target
     * @param ia_parami
     */
    void init(Ia_param * ia_parami);

    void random(int molType,int type, Vector& box) {
        pos.random();
        pos.x *= box.x;
        pos.y *= box.y;
        pos.z *= box.z;

        dir.random();

        // init new Particle
        this->type = type;
        patchdir[0].random();

        //chainIndex = -1;
        this->molType = molType;
        delta_mu = 0;
        switchtype = 0;
        switched = 0;
    }
};

#endif // PARTICLE_H
