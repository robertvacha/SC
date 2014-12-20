/** @file movecreator.h*/

#ifndef MOVECREATOR_H
#define MOVECREATOR_H

#include "totalenergycalculator.h"

extern Topo topo;

class MoveCreator
{
public:
    MoveCreator(Sim* sim, Conf* conf, TotalEnergyCalculator* calcEnergy)
        :  sim(sim), conf(conf), calcEnergy(calcEnergy) {

        try{
            insert.reserve(MAXCHL);
            conlist.reserve(MAXCHL);
        } catch(std::bad_alloc& bad) {
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for insert or conlist, see MoveCreator::constructor\n");
            exit(1);
        }
    }

private:
    Sim* sim;                  // Should contain the simulation options.
    Conf* conf;                // Should contain fast changing particle and box(?) information

    std::vector<Particle > insert;  // insert vector for muVt -> malloc() error when in muVTmove()
    std::vector<ConList > conlist;  // insert vector for muVt -> malloc() error when in muVTmove()

public:
    TotalEnergyCalculator* calcEnergy;

    /**
     * @brief particlemove
     * @return
     */
    double particleMove();

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
     * @brief replicaexchangemove
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
     * @brief psc_rotate rotate spherocylinder by quaternion of random axis and angle smaller than
       maxcos(cosine of angle half), we do everything on site for speed
     * @param psc
     * @param max_angle
     * @param geotype
     */
    void pscRotate(Particle *psc, double max_angle,int geotype);




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
     * @brief cluster_rotate rotate cluster of Particle by quaternion of random axis and angle smaller than
       maxcos(cosine of angle half), we do everything on site for speed
     * @param target
     * @param gc
     * @param max_angle
     */
    void clusterRotate(long target, Vector gc, double max_angle);

    void clusterRotate(vector<Particle>::iterator begin, unsigned int size, double max_angle);

    inline Vector clusterCM(vector<Particle>::iterator begin, unsigned int size) {
        double chainVolume=0.0;
        Vector cluscm(0.0, 0.0, 0.0);

        for(vector<Particle >::iterator it=begin; it!=begin+size; ++it) {
            cluscm.x += it->pos.x * topo.ia_params[it->type][it->type].volume;
            cluscm.y += it->pos.y * topo.ia_params[it->type][it->type].volume;
            cluscm.z += it->pos.z * topo.ia_params[it->type][it->type].volume;

            chainVolume += topo.ia_params[it->type][it->type].volume;
        }

        cluscm.x /= chainVolume;
        cluscm.y /= chainVolume;
        cluscm.z /= chainVolume;

        return cluscm;
    }




    /**
     * @brief movetry Compare energy change to temperature and based on Boltzmann probability
                        return either 0 to accept or 1 to reject the move
     * @param energyold
     * @param energynew
     * @param temperature
     * @return
     */
    int moveTry(double energyold, double energynew, double temperature);




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
    long meshOrderMoveChain(vector<int> chain, Mesh* mesh, long npart, Particle chorig[MAXCHL], int wli);

};

#endif // MOVECREATOR_H
