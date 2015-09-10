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
    Vector patchdir[2];     ///< \brief Vector defining orientation of patches
    Vector patchsides[4];   ///< \brief Vector defining sides of patch
    Vector chdir[2];        ///< \brief Direction for chirality - keep in memory to increase speed

    int molType;           ///< \brief Molecule type 0-100, given sequentialy from 0
    //long chainIndex;        ///< \brief Chain number, only for Molecules of two or more particles
    int type;               ///< \brief Type of the particle 0 - 40
    int switchtype;         ///< \brief With which kind of particle do you want to switch?
    double delta_mu;        ///< \brief Chemical potential for the switch
    int switched;           ///< \brief 0: in initial stat; 1: in the switched stat

    Particle() {}
    Particle(Vector pos, Vector dir, Vector patchDir, int molType, int type) {
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

    /**
     * @brief int_partvec initiate vectors of a single particle
     * @param target
     * @param ia_parami
     */
    void init(Ia_param * ia_parami);

    inline bool operator== (Particle* other) {
        if(this == other) return true;
        else return false;
    }

    inline bool operator!= (Particle* other) {
        if(this != other) return true;
        else return false;
    }

    string info() {
        std::ostringstream o;
        o << "pos:" << pos.info() << "\ndir:" << dir.info() << "\npatchdir1:" << patchdir[0].info() << "\npatchdir2:" << patchdir[1].info()<<endl;
        return o.str();

    }

    /**
     * @brief testInit - dir, patchdir and chdir must be unit vectors
     * @return
     */
    bool testInit() {
        if(dir.dot(dir) > 1.000001 || dir.dot(dir) < 0.999999)
            return false;
        if(pos.x < 0.0 || pos.x > 1.0)
            return false;
        if(pos.y < 0.0 || pos.y > 1.0)
            return false;
        if(pos.z < 0.0 || pos.z > 1.0)
            return false;
        return true;
    }

    inline void rotateRandom(double max_angle, int geotype) {
        pscRotate(max_angle * ran2() , geotype, Vector::getRandomUnitSphere());
    }

    /**
     * @brief psc_rotate rotate spherocylinder by quaternion of random axis and angle smaller than
       maxcos(cosine of angle half), we do everything on site for speed
     * @param part
     * @param angle (radians)
     * @param clockWise
     */
    void pscRotate(double angle, int geotype, Vector newaxis, int clockwise=2) {
        double vc, vs, t2, t3, t4, t5, t6, t7, t8, t9, t10;
        double d1, d2, d3, d4, d5, d6, d7, d8, d9 , newx, newy, newz;
        int k,m;

        /* generate quaternion for rotation*/
        //    maxcos = cos(maxorient/2/180*PI);
        // vc = maxcos + ran2(&seed)*(1-maxcos); /*cos of angle must be bigger than maxcos and smaller than one*/

        vc = cos(angle);

        // GNU compiler wont keep the order of calls for ran2() when taken as an argument for clockwise/counterclockwise
        if ( (ran2()<0.5 && clockwise==2) || clockwise==1 )
            vs = sqrt(1.0 - vc*vc);
        else vs = -sqrt(1.0 - vc*vc); /*randomly choose orientation of direction of rotation clockwise or counterclockwise*/

        Quat newquat(vc, newaxis.x*vs, newaxis.y*vs, newaxis.z*vs);

        /* do quaternion rotation*/
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

        /*rotate spherocylinder direction vector2*/
        newx = 2.0 * ( d1*this->dir.x + d2*this->dir.y + d3*this->dir.z ) + this->dir.x;
        newy = 2.0 * ( d4*this->dir.x + d5*this->dir.y + d6*this->dir.z ) + this->dir.y;
        newz = 2.0 * ( d7*this->dir.x + d8*this->dir.y + d9*this->dir.z ) + this->dir.z;
        this->dir.x = newx;
        this->dir.y = newy;
        this->dir.z = newz;

        m=1;
        if ( (geotype != SCN) && (geotype != SCA) ) {
            if ( (geotype == TPSC) || (geotype == TCPSC) || (geotype == TCHPSC) || (geotype == TCHCPSC) )
                m=2;
            for (k=0;k<m;k++) {
                /*rotate patch direction vector2*/
                newx = 2.0 * ( d1*this->patchdir[k].x + d2*this->patchdir[k].y + d3*this->patchdir[k].z ) + this->patchdir[k].x;
                newy = 2.0 * ( d4*this->patchdir[k].x + d5*this->patchdir[k].y + d6*this->patchdir[k].z ) + this->patchdir[k].y;
                newz = 2.0 * ( d7*this->patchdir[k].x + d8*this->patchdir[k].y + d9*this->patchdir[k].z ) + this->patchdir[k].z;
                this->patchdir[k].x = newx;
                this->patchdir[k].y = newy;
                this->patchdir[k].z = newz;

                /*rotate patch sides vector2s*/
                newx = 2.0 * ( d1*this->patchsides[0+2*k].x + d2*this->patchsides[0+2*k].y + d3*this->patchsides[0+2*k].z ) + this->patchsides[0+2*k].x;
                newy = 2.0 * ( d4*this->patchsides[0+2*k].x + d5*this->patchsides[0+2*k].y + d6*this->patchsides[0+2*k].z ) + this->patchsides[0+2*k].y;
                newz = 2.0 * ( d7*this->patchsides[0+2*k].x + d8*this->patchsides[0+2*k].y + d9*this->patchsides[0+2*k].z ) + this->patchsides[0+2*k].z;
                this->patchsides[0+2*k].x = newx;
                this->patchsides[0+2*k].y = newy;
                this->patchsides[0+2*k].z = newz;
                newx = 2.0 * ( d1*this->patchsides[1+2*k].x + d2*this->patchsides[1+2*k].y + d3*this->patchsides[1+2*k].z ) + this->patchsides[1+2*k].x;
                newy = 2.0 * ( d4*this->patchsides[1+2*k].x + d5*this->patchsides[1+2*k].y + d6*this->patchsides[1+2*k].z ) + this->patchsides[1+2*k].y;
                newz = 2.0 * ( d7*this->patchsides[1+2*k].x + d8*this->patchsides[1+2*k].y + d9*this->patchsides[1+2*k].z ) + this->patchsides[1+2*k].z;
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
                /*rotate chiral direction vector2*/
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
