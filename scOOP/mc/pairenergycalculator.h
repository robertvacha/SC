/** @file pairenergycalculator.h*/

#ifndef PAIRENERGYCALCULATOR_H
#define PAIRENERGYCALCULATOR_H

/*
 *  Calculates energy of a pair of Particle
 *  One particle stored as a internal state
 */

#include <stdlib.h>
#include "../structures/topo.h"
#include <stdio.h>
#include "math_calc.h"
#include "../structures/particle.h"

class PairEnergyCalculator
{
private:
    Particle* part1;
    Particle* part2;
    int num2;

    /*
     *  Parameters pased to functions of interaction calculation
     */
    double dist;                    // closest distance
    Vector distvec;                 // vector of closes distance
    Vector r_cm;                    // Vector between centres of mass from part2 to part1
    double distcm;                  // distance between center of masses
    double dotrcm;                  // square size of r_cm
    double contt;                   // closest point on spherocylinder to sphere

    Topo* topo;
    Vector box; // box size

    double (PairEnergyCalculator::*intFCE[MAXT][MAXT])();


public:
    PairEnergyCalculator(Topo* topo, Vector box) : topo(topo), box(box) {}

    double operator() (Particle* part1, Particle* part2, int num2);

    /**
     * @brief init_intfce Initializes the array with the pointers to the energy function
     * @param topo
     */
    void initIntFCE();

private:

    /**
     * @brief Determines angle energy between spherocylinders
     */
    double angleEnergy();

    /**
     * @brief Determines bond energy
     */
    double bondEnergy();

    /**
     * @brief harmonic    function for calculation of harmonic potential
     * @param aktualvalue
     * @param eqvalue
     * @param springconst
     * @return
     */
    double harmonicPotential(double aktualvalue, double eqvalue, double springconst) {
        return springconst*(aktualvalue-eqvalue)*(aktualvalue-eqvalue)*0.5;
    }

    /**
     * @brief enoexist    Indicates not yet programmed interaction
     * @return
     */
    double eNoExist();

    /**
     * @brief e_psc_psc    Determines total energy of two spherocylinders type PSC PSC
     * @return energy
     */
    double ePscPsc();

    /**
     * @brief eattractive_psc_psc    Determines attractive energy of two spherocylinders type PSC PSC
     * @param patchnum1
     * @param patchnum2
     * @return
     */
    double eattractivePscPsc(int patchnum1,int patchnum2);

    /**
     * @brief e_cpsc_cpsc    Determines total energy of two spherocylinders of type 3 -cylindrical psc -CPSC
     * @return
     */
    double eCpscCpsc();

    /**
     * @brief eattractive_cpsc_cpsc    Determines attractive energy of two spherocylinders of type 3 -cylindrical psc -CPSC
     * @param patchnum1
     * @param patchnum2
     * @return
     */
    double eattractiveCpscCpsc(int patchnum1, int patchnum2);

    /**
     * @brief e_psc_cpsc    Determines total energy of spherocylinders type PSC and CPSC
     * @return
     */
    double ePscCpsc();

    /**
     * @brief eattractive_psc_cpsc    Determines attractive energy of spherocylinders type PSC and CPSC
     * @param patchnum1
     * @param patchnum2
     * @return
     */
    double eattractivePscCpsc(int patchnum1,int patchnum2);

    /**
     * @brief e_spa_sca    Determines total energy of spherocylinder type 1 and sphere type 11
     * @return
     */
    double eSpaSca();

    /**
     * @brief e_psc_spa    Determines total energy of spherocylinder type 2 and sphere type 11
     * @return
     */
    double ePscSpa();

    /**
     * @brief eattractive_psc_spa    Determines attractive energy of spherocylinder type 2 and sphere type 11
     * @param patchnum1
     * @return
     */
    double eattractivePscSpa(int patchnum1);

    /**
     * @brief e_cpsc_spa    Determines total energy of spherocylinder type 3 and sphere type 11
     * @return
     */
    double eCpscSpa();

    /**
     * @brief eattractive_cpsc_spa    Determines attractive energy of spherocylinder type 3 and sphere type 11
     * @param patchnum1
     * @return
     */
    double eattractiveCpscSpa(int patchnum1);

    /**
     * @brief e_2sca_or_2spa    Determines total energy of two spherocylinders type 11
     * @return
     */
    double e2ScaOr2Spa();

    /**
     * @brief e_spn_or_scn    Determines total energy with purely repulsive types
     * @return
     */
    double eSpnOrScn();

    /**
     * @brief closestdist closest distance calculation
     */
    void closestDist();

    /**
     * @brief fanglscale  a = r_ij * n_i
     * @param a
     * @param which
     * @return
     */
    double fanglScale(double a, int which);

    /**
     * @brief erepulsive    Determines repulsive energy of two spherocylinders
     * @return
     */
    double eRepulsive();

    /**
     * @brief mindist_segments Find closest distance between line segments and return its vector
       gets orientations and lengths of line segments and the vector connecting
       their center os masses (from vec1 to vec2)
     * @param dir1
     * @param halfl1
     * @param dir2
     * @param halfl2
     * @param r_cm
     * @return
     */
    Vector minDistSegments(double halfl1, double halfl2, Vector r_cm);

    /**
     * @brief psc_intersect Calculate intersections of sc2 with a patch of sc1 and return them in
     * @param halfl1
     * @param halfl2
     * @param r_cm
     * @param intersections
     * @param which
     * @param patchnum
     * @return
     */
    int pscIntersect(Particle* part1, Particle* part2, double halfl1, double halfl2, Vector r_cm,
                      double intersections[5], int which, int patchnum);

    /**
     * @brief test_intrpatch Find if vector vec has angular intersection with patch of sc1
     * @param part1
     * @param vec
     * @param cospatch
     * @param ti
     * @param intersections
     * @param patchnum
     * @return
     */
    int testIntrPatch(Particle * part1, Vector vec, double cospatch,
        double ti, double intersections[5],int patchnum);

    /**
     * @brief find_intersect_plane Find intersections of SC and plane defined by vector w_vec.and returns number of them
     * @param part1
     * @param part2
     * @param halfl2
     * @param r_cm
     * @param w_vec
     * @param rcut
     * @param cospatch
     * @param intersections
     * @return
     */
    int findIntersectPlane(Particle * part1, Particle * part2, double halfl2,
            Vector r_cm, Vector w_vec, double cospatch, double intersections[5]);

    /**
     * @brief cpsc_intersect Calculate intersections of sc2 with a patch of sc1
                             and return them in this works for cylindrical psc -CPSC
     * @param part1
     * @param part2
     * @param halfl1
     * @param halfl2
     * @param r_cm
     * @param intersections
     * @param which
     * @param patchnum
     * @return
     */
    int cpscIntersect(Particle * part1, Particle * part2, double halfl1, double halfl2, Vector r_cm,
                       double intersections[5], int which, int patchnum);

    /**
     * @brief find_intersect_planec Find intersections of plane defined by vector w_vec.
                                    and returns number of them - for cylindrical psc -CPSC
     * @param part1
     * @param part2
     * @param halfl
     * @param r_cm
     * @param w_vec
     * @param cospatch
     * @param intersections
     * @return
     */
    int findIntersectPlanec(Particle * part1, Particle * part2, double halfl,
            Vector r_cm, Vector w_vec, double cospatch, double intersections[5]);

    /**
     * @brief vec_perpproject vector projection of vector A perpendicular to direction of B
     * @param A
     * @param B
     * @return
     */
    Vector vecPerpProject(Vector *A,Vector *B);
};



#endif // PAIRENERGYCALCULATOR_H
