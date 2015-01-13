#ifndef PARTICLE_H
#define PARTICLE_H

#include "structures.h"
#include "topo.h"

extern Topo topo;

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

    void init(Vector pos, Vector dir, Vector patchDir, int molType, int type) {
        this->pos = pos;
        this->dir = dir;
        this->patchdir[0] = patchDir;

        // init new Particle
        this->type = type;

        //chainIndex = -1;
        this->molType = molType;
        delta_mu = 0;
        switchtype = 0;
        switched = 0;
    }

    inline bool operator== (Particle& other) {
        if(pos == other.pos) return true;
        else return false;
    }

    inline bool operator!= (Particle& other) {
        if(pos != other.pos) return true;
        else return false;
    }

    string info() {
        std::ostringstream o;
        o << "pos:" << pos.info() << "\ndir:" << dir.info() << "\npatchdir:" << patchdir[0].info();
        return o.str();

    }

    /**
     * @brief psc_rotate rotate spherocylinder by quaternion of random axis and angle smaller than
       maxcos(cosine of angle half), we do everything on site for speed
     * @param part
     * @param angle (radians)
     * @param clockWise
     */
    void pscRotate(Vector axis, double angle, bool clockWise) {
        int geotype = topo.ia_params[type][type].geotype[0]; // ASK ROBERT -> difference [0] [1]
        double vc, vs, t2, t3, t4, t5, t6, t7, t8, t9, t10;
        double d1, d2, d3, d4, d5, d6, d7, d8, d9;
        double newx,newy,newz;
        int k,m;

        vc = cos(angle);

        if (clockWise) vs = sqrt(1.0 - vc*vc);
        else vs = -sqrt(1.0 - vc*vc); //randomly choose orientation of direction of rotation clockwise or counterclockwise

        Quat newquat(vc, axis.x*vs, axis.y*vs, axis.z*vs);

        // do quaternion rotation
        t2 =  newquat.w * newquat.x;
        t3 =  newquat.w * newquat.y;
        t4 =  newquat.w * newquat.z;
        t5 = -newquat.x * newquat.x;
        t6 =  newquat.x * newquat.y;
        t7 =  newquat.x * newquat.z;
        t8 = -newquat.y * newquat.y;
        t9 =  newquat.y * newquat.z;
        t10 = -newquat.z * newquat.z;

        d1 = t8 + t10;
        d2 = t6 - t4;
        d3 = t3 + t7;
        d4 = t4 + t6;
        d5 = t5 + t10;
        d6 = t9 - t2;
        d7 = t7 - t3;
        d8 = t2 + t9;
        d9 = t5 + t8;

        //rotate spherocylinder direction vector2
        newx = 2.0 * ( d1*dir.x + d2*dir.y + d3*dir.z ) + dir.x;
        newy = 2.0 * ( d4*dir.x + d5*dir.y + d6*dir.z ) + dir.y;
        newz = 2.0 * ( d7*dir.x + d8*dir.y + d9*dir.z ) + dir.z;
        this->dir.x = newx;
        this->dir.y = newy;
        this->dir.z = newz;

        m=1;
        if ( (geotype != SCN) && (geotype != SCA) ) {
            if ( (geotype == TPSC) || (geotype == TCPSC) || (geotype == TCHPSC) || (geotype == TCHCPSC) )
            m=2;
            for (k=0;k<m;k++) {
            //rotate patch direction vector2
            newx = 2.0 * ( d1*this->patchdir[k].x + d2*this->patchdir[k].y + d3*this->patchdir[k].z ) + this->patchdir[k].x;
            newy = 2.0 * ( d4*this->patchdir[k].x + d5*this->patchdir[k].y + d6*this->patchdir[k].z ) + this->patchdir[k].y;
            newz = 2.0 * ( d7*this->patchdir[k].x + d8*this->patchdir[k].y + d9*this->patchdir[k].z ) + this->patchdir[k].z;
            this->patchdir[k].x = newx;
            this->patchdir[k].y = newy;
            this->patchdir[k].z = newz;

            //rotate patch sides vector2s
            newx = 2.0 * ( d1*patchsides[0+2*k].x + d2*patchsides[0+2*k].y + d3*patchsides[0+2*k].z ) + patchsides[0+2*k].x;
            newy = 2.0 * ( d4*patchsides[0+2*k].x + d5*patchsides[0+2*k].y + d6*patchsides[0+2*k].z ) + patchsides[0+2*k].y;
            newz = 2.0 * ( d7*patchsides[0+2*k].x + d8*patchsides[0+2*k].y + d9*patchsides[0+2*k].z ) + patchsides[0+2*k].z;
            this->patchsides[0+2*k].x = newx;
            this->patchsides[0+2*k].y = newy;
            this->patchsides[0+2*k].z = newz;
            newx = 2.0 * ( d1 * patchsides[1+2*k].x + d2 * patchsides[1+2*k].y + d3 *patchsides[1+2*k].z ) + patchsides[1+2*k].x;
            newy = 2.0 * ( d4 * patchsides[1+2*k].x + d5 * patchsides[1+2*k].y + d6 *patchsides[1+2*k].z ) + patchsides[1+2*k].y;
            newz = 2.0 * ( d7 * patchsides[1+2*k].x + d8 * patchsides[1+2*k].y + d9 * patchsides[1+2*k].z ) + patchsides[1+2*k].z;
            this->patchsides[1+2*k].x = newx;
            this->patchsides[1+2*k].y = newy;
            this->patchsides[1+2*k].z = newz;
            }
        }

        m=1;
        if ( (geotype == CHPSC) || (geotype == CHCPSC) || (geotype == TCHPSC) || (geotype == TCHCPSC) ) {
            if ( (geotype == TCHPSC) || (geotype == TCHCPSC) )
            m=2;
            for (k=0;k<m;k++) {
            //rotate chiral direction vector2
            newx = 2.0 * ( d1*this->chdir[k].x + d2*this->chdir[k].y + d3*this->chdir[k].z ) + this->chdir[k].x;
            newy = 2.0 * ( d4*this->chdir[k].x + d5*this->chdir[k].y + d6*this->chdir[k].z ) + this->chdir[k].y;
            newz = 2.0 * ( d7*this->chdir[k].x + d8*this->chdir[k].y + d9*this->chdir[k].z ) + this->chdir[k].z;
            this->chdir[k].x = newx;
            this->chdir[k].y = newy;
            this->chdir[k].z = newz;
            }
        }
    }
};

#endif // PARTICLE_H
