/** @file externalenergycalculator.h*/

#ifndef EXTERNALENERGYCALCULATOR_H
#define EXTERNALENERGYCALCULATOR_H

#include "../structures/sim.h"

class ExternalEnergyCalculator
{
private:
    Vector* box;

    double dist;            /* closest distance */
    Particle* part1;       /* particle 1 */
    Ia_param * param;       /* interaction parameters */
    bool positive;
    bool orientin;
    double orient;
    double rcmz;            /* z distance between*/
    double interendz;       /* z coordinate of interaction end*/
    Vector project;         /*vector2 for projection down to plane */

public:
    ExternalEnergyCalculator(Vector* box) : box(box) {}

    double extere2(Particle* target);

private:
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
};

#endif // EXTERNALENERGYCALCULATOR_H
