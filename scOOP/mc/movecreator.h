/** @file movecreator.h*/

#ifndef MOVECREATOR_H
#define MOVECREATOR_H


#include "totalenergycalculator.h"
#include "../structures/Conf.h"

#include "randomGenerator.h"

class MoveCreator
{
public:
    MoveCreator(Topo* topo, Sim* sim, Conf* conf, TotalEnergyCalculator* calcEnergy)
        : topo(topo), sim(sim), conf(conf), calcEnergy(calcEnergy) {

        try{
            insert.reserve(MAXCHL);
            conlist.reserve(MAXCHL);
        } catch(std::bad_alloc& bad) {
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for insert or conlist, see MoveCreator::constructor\n");
            exit(1);
        }
    }

private:
    Topo* topo;                // will maybe contain all the topo stuff in future
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

    /**
     * @brief radiushole_all filling the radiushole above vec
     * @param wli
     * @param position
     * @return
     */
    long radiusholeAll(int wli, Vector *position);

    /**
     * @brief contParticle_all filling all Particle in the contact
     * @param wli
     * @return
     */
    long contParticlesAll(int wli);



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

    void clusterRotate(vector<Particle>::iterator begin, unsigned int size);


    /**
     * @brief radiushole_position return order of given radius
     * @param radius
     * @param wli
     * @return
     */
    long radiusholePosition(double radius, int wli) {
        return (long) ceil( ( radius - sim->wl.minorder[wli]) / sim->wl.dorder[wli]  );
    }

    /**
     * @brief radiushole_order return current bin of free radius
     * @return
     */
    long radiusholeOrder();

    /**
     * @brief radiusholeorder_moveone return change in order parameter when one particle moves
     * @param oldpos
     * @param target
     * @param wli
     * @param position
     * @return
     */
    long radiusholeOrderMoveOne(Vector *oldpos, long target, int wli, Vector *position);

    /**
     * @brief radiusholeorder_movechain return change in order parameter when chain moves
     * @param chain
     * @param chorig
     * @param wli
     * @param position
     * @return
     */
    long radiusholeOrderMoveChain(vector<int> chain, Particle chorig[MAXCHL], int wli, Vector *position);

    /**
     * @brief radiushole_print
     * @param radiushole
     * @param length
     */
    void radiusholePrint(long *radiushole, long length);




    /**
     * @brief particleinncontact returns if particle is in contact
     * @param vec
     * @return
     */
    bool particlesInContact(Vector *vec);

    /**
     * @brief contParticle_order return order for Particle in contact
     * @param wli
     * @return
     */
    inline long contParticlesOrder(int wli) {
        return (long) ceil( ( sim->wl.partincontact - sim->wl.minorder[wli]) / sim->wl.dorder[wli]  );
    }

    /**
     * @brief contParticle_moveone return change in number of Particle in contact when one particle moves
     * @param oldpos
     * @param target
     * @param wli
     * @return
     */
    long contParticlesMoveOne(Vector *oldpos, long target,int wli);

    /**
     * @brief contParticle_movechain return change in order parameter when chain moves
     * @param chain
     * @param conf
     * @param sim
     * @param chorig
     * @param wli
     * @return
     */
    long contParticlesMoveChain(vector<int> chain, Particle chorig[MAXCHL], int wli);

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
