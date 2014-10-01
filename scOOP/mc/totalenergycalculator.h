/** @file totalenergycalculator.h*/

#ifndef TOTALENERGYCALCULATOR_H
#define TOTALENERGYCALCULATOR_H

#include <stdlib.h> // exit()

#include <iostream>

#include "pairenergycalculator.h"
#include "../structures/sim.h"

using namespace std;

class TotalEnergyCalculator
{
public:
    PairEnergyCalculator pairE;

private:
    Topo* topo;
    Sim* sim;
    Conf* conf;

    double dist;            /* closest distance */
    Vector distvec;         /* vector2 of closes distance */
    Particle* part1;       /* particle 1 */
    Particle* part2;       /* particle 2 */
    Vector box;             /* box size */
    Ia_param * param;       /* interaction parameters */
    Vector r_cm;            /* vector2 connecting center of masses */
    double distcm;          /* distance between center of masses */
    double dotrcm;          /* square size of r_cm*/
    double contt;           /* closest point on spherocylinder to sphere */
    bool positive;
    bool orientin;
    double orient;
    double rcmz;            /* z distance between*/
    double interendz;       /* z coordinate of interaction end*/
    Vector project;         /*vector2 for projection down to plane */

public:
    TotalEnergyCalculator(Topo * topo, Sim * sim, Conf * conf): pairE(topo, conf->box, conf->conlist),
        topo(topo), sim(sim), conf(conf) {pairE.initIntFCE();}

    /**
     * @brief calc_energy Calculate the different energy contributions.
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
     */
    double chainToAll(int target, int chainnum);

    /**
     * @brief Calculates energy between particle "target" and the rest
     */
    double oneToAll(int target);

    /**
     * @brief Calculates energy between particle "target" and the rest
     */
    double oneToAll(Particle* target, int conlistTarget);

    /**
     * @brief Calculates energy between all pairs. Returns energy
     * @return
     */
    double allToAll();


private:
    /**
     * @brief extere2    Calculates interaction of target particle and external field version 2
                         calculate projection of spherocylinder in direction of patch and calculate
                         interacting line segment within cutoff
     * @param target
     * @return
     */
    double extere2(int target);
    double extere2(Particle* target, int overload);

    /**
     * @brief exter2ClosestDist
     * @param interact
     */
    void exter2ClosestDist();

    /**
     * @brief exter2_atre
     * @param interact
     * @param ndist
     * @param numpatch
     * @param halfl
     * @return
     */
    double exter2Atre(double* ndist, int numpatch, double halfl);

    /**
     * @brief exter_atre    Calculates interaction of target particle and external field
                            calculate projection of patch of spherocylinder on wall
                            evaluate intersection area and calculate interaction from that
     * @param ndist
     * @param numpatch
     * @param halfl
     * @return
     *
     * @deprecated
     */
    double exterAtre(double *ndist,int numpatch,double halfl);

    /**
     * @brief areaeightpoints calculates area defined by six points in z=0  plane
     * @param p1
     * @param p2
     * @param p3
     * @param p4
     * @param p5
     * @param p6
     * @param p7
     * @param p8
     * @return
     */
    double areaEightPoints(Vector * p1, Vector * p2, Vector * p3, Vector * p4,
                         Vector * p5, Vector * p6,Vector * p7, Vector * p8);

    /**
     * @brief psc_wall Find projection of cpsc on plane (0,0,1) including cutoff and return
                        vector2 to its begining and end and cm
     * @param pbeg
     * @param pend
     * @param projectdir
     * @param partdir
     * @param cutdist
     * @param partbeg
     * @param partend
     * @return
     */
    int pscWall(Vector* pbeg, Vector* pend, Vector* projectdir, Vector* partdir,
                 double * cutdist, Vector *partbeg, Vector *partend);

    /**
     * @brief cpsc_wall
     * @param pbeg
     * @param pend
     * @param projectdir
     * @param partdir
     * @param halfl
     * @param cutdist
     * @param partbeg
     * @param partend
     * @return
     */
    int cpscWall(Vector* pbeg, Vector* pend, Vector* projectdir, Vector* partdir,
                  double* halfl, double * cutdist, Vector *partbeg, Vector *partend);

    /**
     * @brief cutprojectatwall
     * @param pextr1
     * @param pextr2
     * @param pextr3
     * @param pextr4
     * @param projectdir
     * @param partdir
     * @param cutdist
     * @param partbeg
     * @param partend
     * @param pend
     * @param cuttoproject
     * @return
     */
    int cutProjectAtWall(Vector* pextr1, Vector* pextr2, Vector* pextr3, Vector* pextr4,
                         Vector* projectdir, Vector* partdir, double * cutdist,
                         Vector *partbeg, Vector *partend, Vector *pend, double *cuttoproject);

    /**
     * @brief projectinz project a point in project direction to z plane z=0
     * @param vec1
     * @param projectdir
     * @param projection
     */
    inline void projectinZ(Vector* vec1, Vector* projectdir,Vector * projection) {
        projection->x = vec1->x - vec1->z * projectdir->x/projectdir->z;
        projection->y = vec1->y - vec1->z * projectdir->y/projectdir->z;
        projection->z = 0;
    }

};

#endif // TOTALENERGYCALCULATOR_H
