/** @file movecreator.h*/

#ifndef MOVECREATOR_H
#define MOVECREATOR_H

#include "totalenergycalculator.h"
#include "wanglandau.h"

extern Topo topo;


class MoveCreator
{
public:
    MoveCreator(Sim* sim, Conf* conf, TotalEnergyCalculator* calcEnergy)
        :  sim(sim), conf(conf), wl(conf, sim), calcEnergy(calcEnergy) {

        try{
            insert.reserve(MAXCHL);
        } catch(std::bad_alloc& bad) {
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for insert, see MoveCreator::constructor\n");
            exit(1);
        }
    }

private:
    Sim* sim;                  // Should contain the simulation options.
    Conf* conf;                // Should contain fast changing particle and box(?) information  

    std::vector<Particle > insert;  // insert vector for muVt -> malloc() error when in muVTmove()
    Particle chorig[MAXCHL];

public:
    WangLandau wl;
    TotalEnergyCalculator* calcEnergy;

    /**
     * @brief particlemove
     * @return
     */
    double particleMove();

    /**
     * @brief Function to select clusters and move/rotate whole cluster
     * @param
     * @return
     */
    double clusterMove();

    int isInCluster(double *list, int size, double value);

    /**
     * @brief Function to produce one geometrical cluster move ()
     * @param
     * @return
     */
    double clusterMoveGeom(long target);

    double printClustersConf();

    /**
     * @brief switchtypemove This is an attempt to switch a type
     * @return
     */
    double switchTypeMove();

    /**
     * @brief chainmove
     * @return
     */
    double chainMove();

    /**
     * @brief pressuremove
     * @return
     */
    double pressureMove();

    /**
     * @brief replicaexchangemove for canocnical, grandcanonical and isobaric-isotermal ensembles
     *
     * We are switching statistics and Hamiltonians. Configurations remains the same
     *
     * Data to switch:
     *
     * Statistics, pressure, temperature, Wang-Landau data
     *
     * Data to stay:
     *
     * Conf, Geometry, MoleculeParams, Ia_param
     *
     *
     *
     * @param sweep
     * @return
     */
    double replicaExchangeMove(long sweep );

    /**
     * @brief muVTMove
     * @return
     */
    double muVTMove();

private:
    int getRandomMuVTType();

    /**
     * @brief partdisplace
     * @param target
     * @return
     */
    double partDisplace(long target);

    /**
     * @brief partrotate
     * @param target
     * @return
     */
    double partRotate(long target);

    /**
     * @brief partAxialRotate(), take particle and rotate it (anti or clockwise) by random angle around particle axis
     * @param target
     * @return edriftchange
     */
    double partAxialRotate(long target);

    /**
     * @brief chaindisplace
     * @param target
     * @return
     */
    double chainDisplace(long target);

    /**
     * @brief chainrotate
     * @param target
     * @return
     */
    double chainRotate(long target);


    /**
     * @brief clusterRotate rotate cluster of Particle by quaternion of random axis and angle smaller than
       maxcos(cosine of angle half), we do everything on site for speed
     * @param begin
     * @param size
     * @param max_angle
     */
    void clusterRotate(vector<int >& cluster, double max_angle);

    void clusterRotate(vector<Particle >& cluster, double max_angle);

    inline Vector clusterCM(vector<int >& cluster);

    inline Vector clusterCM(vector<Particle >& cluster);




    /**
     * @brief movetry Compare energy change to temperature and based on Boltzmann probability
                        return either 0 to accept or 1 to reject the move
     * @param energyold
     * @param energynew
     * @param temperature
     * @return false - move accepted, True - move rejected
     */
    inline bool moveTry(double energyold, double energynew, double temperature) {
        /*DEBUG   printf ("   Move trial:    %13.8f %13.8f %13.8f %13.8f\n",
      energynew, energyold, temperature, ran2(&seed));*/
        if (energynew <= energyold ) {
            return false;
        } else {
            if (exp(-1.0*(energynew-energyold)/temperature) > ran2()) {
                return false;
            } else {
                return true;
            }
        }
    }

    /*
     *  MESH METHODS
     */

    /**
     * @brief meshorder_moveone return change in order parameter when one particle moves
     * @param oldpos
     * @param newpos
     * @param mesh
     * @param npart
     * @param target
     * @param wli
     * @return
     */
    long meshOrderMoveOne(Vector oldpos, Vector newpos, Mesh* mesh, long npart, long target, int wli);

    /**
     * @brief meshorder_movechain return change in order parameter when chain moves
     * @param chain
     * @param mesh
     * @param npart
     * @param chorig
     * @param wli
     * @return
     */
    long meshOrderMoveChain(Molecule chain, Mesh* mesh, long npart, Particle chorig[MAXCHL], int wli);




    /*
     *  BACKWARD COMPATIBILITY FUNCTIONS
     */

    /*void pscRotate(Particle *psc, double max_angle, int geotype) {
        double vc, vs, t2, t3, t4, t5, t6, t7, t8, t9, t10;
        double d1, d2, d3, d4, d5, d6, d7, d8, d9 , newx, newy, newz;
        int k,m;
        Vector newaxis;

        // generate quaternion for rotation
        newaxis.randomUnitSphere(); // random axes for rotation
        //    maxcos = cos(maxorient/2/180*PI);
        // vc = maxcos + ran2(&seed)*(1-maxcos); //cos of angle must be bigger than maxcos and smaller than one

        vc = cos(max_angle * ran2() );

        if (ran2() <0.5) vs = sqrt(1.0 - vc*vc);
        else vs = -sqrt(1.0 - vc*vc); //randomly choose orientation of direction of rotation clockwise or counterclockwise

        Quat newquat(vc, newaxis.x*vs, newaxis.y*vs, newaxis.z*vs);

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

        /rotate spherocylinder direction vector2
        newx = 2.0 * ( d1*psc->dir.x + d2*psc->dir.y + d3*psc->dir.z ) + psc->dir.x;
        newy = 2.0 * ( d4*psc->dir.x + d5*psc->dir.y + d6*psc->dir.z ) + psc->dir.y;
        newz = 2.0 * ( d7*psc->dir.x + d8*psc->dir.y + d9*psc->dir.z ) + psc->dir.z;
        psc->dir.x = newx;
        psc->dir.y = newy;
        psc->dir.z = newz;

        m=1;
        if ( (geotype != SCN) && (geotype != SCA) ) {
            if ( (geotype == TPSC) || (geotype == TCPSC) || (geotype == TCHPSC) || (geotype == TCHCPSC) )
                m=2;
            for (k=0;k<m;k++) {
                //rotate patch direction vector2
                newx = 2.0 * ( d1*psc->patchdir[k].x + d2*psc->patchdir[k].y + d3*psc->patchdir[k].z ) + psc->patchdir[k].x;
                newy = 2.0 * ( d4*psc->patchdir[k].x + d5*psc->patchdir[k].y + d6*psc->patchdir[k].z ) + psc->patchdir[k].y;
                newz = 2.0 * ( d7*psc->patchdir[k].x + d8*psc->patchdir[k].y + d9*psc->patchdir[k].z ) + psc->patchdir[k].z;
                psc->patchdir[k].x = newx;
                psc->patchdir[k].y = newy;
                psc->patchdir[k].z = newz;

                //rotate patch sides vector2s
                newx = 2.0 * ( d1*psc->patchsides[0+2*k].x + d2*psc->patchsides[0+2*k].y + d3*psc->patchsides[0+2*k].z ) + psc->patchsides[0+2*k].x;
                newy = 2.0 * ( d4*psc->patchsides[0+2*k].x + d5*psc->patchsides[0+2*k].y + d6*psc->patchsides[0+2*k].z ) + psc->patchsides[0+2*k].y;
                newz = 2.0 * ( d7*psc->patchsides[0+2*k].x + d8*psc->patchsides[0+2*k].y + d9*psc->patchsides[0+2*k].z ) + psc->patchsides[0+2*k].z;
                psc->patchsides[0+2*k].x = newx;
                psc->patchsides[0+2*k].y = newy;
                psc->patchsides[0+2*k].z = newz;
                newx = 2.0 * ( d1*psc->patchsides[1+2*k].x + d2*psc->patchsides[1+2*k].y + d3*psc->patchsides[1+2*k].z ) + psc->patchsides[1+2*k].x;
                newy = 2.0 * ( d4*psc->patchsides[1+2*k].x + d5*psc->patchsides[1+2*k].y + d6*psc->patchsides[1+2*k].z ) + psc->patchsides[1+2*k].y;
                newz = 2.0 * ( d7*psc->patchsides[1+2*k].x + d8*psc->patchsides[1+2*k].y + d9*psc->patchsides[1+2*k].z ) + psc->patchsides[1+2*k].z;
                psc->patchsides[1+2*k].x = newx;
                psc->patchsides[1+2*k].y = newy;
                psc->patchsides[1+2*k].z = newz;
            }
        }

        m=1;
        if ( (geotype == CHPSC) || (geotype == CHCPSC) || (geotype == TCHPSC) || (geotype == TCHCPSC) ) {
            if ( (geotype == TCHPSC) || (geotype == TCHCPSC) )
                m=2;
            for (k=0;k<m;k++) {
                //rotate chiral direction vector2
                newx = 2.0 * ( d1*psc->chdir[k].x + d2*psc->chdir[k].y + d3*psc->chdir[k].z ) + psc->chdir[k].x;
                newy = 2.0 * ( d4*psc->chdir[k].x + d5*psc->chdir[k].y + d6*psc->chdir[k].z ) + psc->chdir[k].y;
                newz = 2.0 * ( d7*psc->chdir[k].x + d8*psc->chdir[k].y + d9*psc->chdir[k].z ) + psc->chdir[k].z;
                psc->chdir[k].x = newx;
                psc->chdir[k].y = newy;
                psc->chdir[k].z = newz;
            }
        }
    }*/

    /*void clusterRotate(long target, Vector gc, double max_angle) {
        long current,i;
        double vc,vs;
        //double quatsize;
        Vector newaxis;

        // create rotation quaternion
        newaxis.randomUnitSphere(); //random axes for rotation
        //    maxcos = cos(maxorient/2/180*PI);
        //vc = maxcos + ran2(&seed)*(1-maxcos); //cos of angle must be bigger than maxcos and smaller than one
        vc = cos(max_angle * ran2() );
        if (ran2() <0.5) vs = sqrt(1.0 - vc*vc);
        else vs = -sqrt(1.0 - vc*vc); //randomly choose orientation of direction of rotation clockwise or counterclockwise

        Quat newquat(vc, newaxis.x*vs, newaxis.y*vs, newaxis.z*vs);

        //quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);

        //shift position to geometrical center
        i=0;
        current = conf->pvec.getChainPart(target,0);
        while (current >=0 ) {
            //shift position to geometrical center
            conf->pvec[current].pos.x -= gc.x;
            conf->pvec[current].pos.y -= gc.y;
            conf->pvec[current].pos.z -= gc.z;
            //scale things by geo.box not to have them distorted
            conf->pvec[current].pos.x *= conf->geo.box.x;
            conf->pvec[current].pos.y *= conf->geo.box.y;
            conf->pvec[current].pos.z *= conf->geo.box.z;
            //do rotation
            conf->pvec[current].pos.rotate(newquat);
            conf->pvec[current].dir.rotate(newquat);
            conf->pvec[current].patchdir[0].rotate(newquat);
            conf->pvec[current].patchdir[1].rotate(newquat);
            conf->pvec[current].chdir[0].rotate(newquat);
            conf->pvec[current].chdir[1].rotate(newquat);
            conf->pvec[current].patchsides[0].rotate(newquat);
            conf->pvec[current].patchsides[1].rotate(newquat);
            conf->pvec[current].patchsides[2].rotate(newquat);
            conf->pvec[current].patchsides[3].rotate(newquat);
            //sclae back
            conf->pvec[current].pos.x /= conf->geo.box.x;
            conf->pvec[current].pos.y /= conf->geo.box.y;
            conf->pvec[current].pos.z /= conf->geo.box.z;
            //shift positions back
            conf->pvec[current].pos.x += gc.x;
            conf->pvec[current].pos.y += gc.y;
            conf->pvec[current].pos.z += gc.z;
            i++;
            current = conf->pvec.getChainPart(target,i);
        }
    }*/
};

#endif // MOVECREATOR_H
