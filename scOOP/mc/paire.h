/** @file paire.h*/

#ifndef PAIRE_H
#define PAIRE_H

#include "../structures/topo.h"
#include "../structures/Conf.h"

inline double fanglScale(double a, double pcangl, double pcanglsw) { // TODO for different types
    if (a <= pcanglsw)
        return 0.0;
    else {
        if (a >= pcangl)
            return 1.0;
        else {
            return 0.5 - ((pcanglsw + pcangl)*0.5 - a ) / (pcangl - pcanglsw);
        }
    }
    assert(false);
    return 0.0;
}

/**
 * @brief closestDist - between Spherocylinder and sphere
 * @param r_cm - distance of center masses
 * @param dir1 - direction vector of spherocylinder
 * @param halfl - half lenght of the spherocylinder
 * @return distance squared
 */
inline double closestDist(const Vector&& r_cm, const Vector& dir1, double halfl, double& contt, Vector& distvec) {
    double d;

    contt = DOT(dir1,r_cm);
    d = -halfl;
    if (contt >= halfl)
        d = halfl;
    else {
        if (contt > -halfl)
            d = contt;
    }
    distvec.x = - r_cm.x + dir1.x * d;
    distvec.y = - r_cm.y + dir1.y * d;
    distvec.z = - r_cm.z + dir1.z * d;

    return DOT(distvec, distvec);
}



class EPatchToSphere {
public:
    virtual double operator() (double dist, double contt, Vector& distvec, const Ia_param& iaParam, const Vector& p1Dir, const Patch p1P)  {
        return 0.0;
    }

    virtual void getGeoTypes(const Ia_param& iaParam, bool& chiral, bool& sec)  {
        chiral=false; sec=false;
    }

    virtual bool isFar(double contt, double halfl) {
        return false;
    }
};

class EPatch {
public:
    virtual double operator() (const Ia_param& iaParam, const Vector& p1Dir, const Vector& p2Dir, const Patch p1P, const Patch p2P, const Vector& r_cm,
                               int patchnum1, int patchnum2) {
        return 0.0;
    }

    virtual void getGeoTypes(const Ia_param& iaParam, bool& firstCH, bool& secondCH, bool& firstT, bool& secondT) {
        firstCH = false; secondCH = false; firstT = false; secondT = false;
    }
};

class EPotential {
public:
    virtual double operator() (double dist, const Particle* part1, const Particle* part2) {
        return 0.0;
    }
};

class EAngle {
public:
    virtual double operator() (double dist, const Particle* part1, const Particle* part2, const ConList* conlist) {
        return 0.0;
    }
};

class EBond {
public:
    virtual double operator() (double dist, const Particle* part1, const Particle* part2, const ConList* conlist) {
        return 0.0;
    }
};

class EBasic {
public:
    virtual double operator() (double dist, const Vector& r_cm, const Particle* part1, const Particle* part2, const ConList* conlist) {
        fprintf (stderr, "ERROR: We have not programed interaction of types %d and %d\n", part1->type,part2->type);
        return 0.0;
    }
};




class HarmonicSp : public EBond {
public:
    double operator() (double dist, const Particle* part1, const Particle* part2, const ConList* conlist) override {
        assert(conlist != nullptr);

        // interaction with nearest neighbours
        if ((topo.moleculeParam[part1->molType]).bond1c >= 0)
            if (part2 == conlist->conlist[1] || part2 == conlist->conlist[0])
                return harmonicPotential(dist,topo.moleculeParam[part1->molType].bond1eq,topo.moleculeParam[part1->molType].bond1c);

        // interaction with second nearest neighbours
        if (topo.moleculeParam[part1->molType].bond2c >= 0)
            if (part2 == conlist->conlist[2] || part2 == conlist->conlist[3])
                return harmonicPotential(dist,topo.moleculeParam[part1->molType].bond2eq,topo.moleculeParam[part1->molType].bond2c);

        // interaction with nearest neighbours - direct harmonic bond
        if ((topo.moleculeParam[part1->molType]).bonddc > 0)
            if (part2 == conlist->conlist[1] || part2 == conlist->conlist[0])
                return harmonicPotential(dist,topo.moleculeParam[part1->molType].bonddeq,topo.moleculeParam[part1->molType].bonddc);

        // interaction with nearest neighbours - direct harmonic bond hinged
        if ((topo.moleculeParam[part1->molType]).bondhc > 0)
            if (part2 == conlist->conlist[1] || part2 == conlist->conlist[0])
                return harmonicPotential(dist,topo.moleculeParam[part1->molType].bondheq,topo.moleculeParam[part1->molType].bondhc);

        return 0.0;
    }

private:
    /**
     * @brief harmonic    function for calculation of harmonic potential
     * @param aktualvalue
     * @param eqvalue
     * @param springconst
     * @return
     */
    inline double harmonicPotential(double aktualvalue, double eqvalue, double springconst) {
        return springconst*(aktualvalue-eqvalue)*(aktualvalue-eqvalue)*0.5;
    }
};




class Lj : public EPotential {
public:
    double operator() (double dist, const Particle* part1, const Particle* part2) override {
        return topo.ia_params[part1->type][part2->type].A * pow(dist, -12) - topo.ia_params[part1->type][part2->type].B * pow(dist, -6);
    }
};

/**
 * @brief The WcaTrunc class - shifted truncated WCA potential (repulsive only)
 */
class WcaTrunc : public EPotential {
public:
    double operator() (double dist, const Particle* part1, const Particle* part2) override {
        if (dist > topo.ia_params[part1->type][part2->type].rcutwca)
            return 0.0;
        return topo.ia_params[part1->type][part2->type].A * pow(dist, -12) - topo.ia_params[part1->type][part2->type].B * pow(dist, -6) + topo.ia_params[part1->type][part2->type].epsilon;
    }
};

class WcaTruncSq : public EPotential {
public:
    double operator() (double distSq, const Particle* part1, const Particle* part2) override {
        if (distSq > topo.ia_params[part1->type][part2->type].rcutwcaSq)
            return 0.0;
        else
            return topo.ia_params[part1->type][part2->type].epsilon + topo.ia_params[part1->type][part2->type].A * pow(distSq, -6)- topo.ia_params[part1->type][part2->type].B * pow(distSq, -3);
    }
};

/**
 * @brief The WcaCos2 class - shifted truncated WCA potential + cos^2 potential
 *
 * Decreasing distance, energy:
 * too large : 0
 * below rcut -> attraction fce
 * below pdis -> -epsilon
 * below rcutWCA -> -epsilon + repulsive fce
 *
 * 0        rcutWCA         pdis        rcut        infinity
 */
class WcaCos2 : public EPotential {
public:
    double operator() (double dist, const Particle* part1, const Particle* part2) override {
        double e = 0.0;
        if(dist > topo.ia_params[part1->type][part2->type].rcut ||
                (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
                topo.ia_params[part1->type][part2->type].exclude)
            return 0.0;

        if(dist > topo.ia_params[part1->type][part2->type].pdis) {
            // taylor cos, precise to ~10^(-7)
            e = PIH * (dist - topo.ia_params[part1->type][part2->type].pdis) * topo.ia_params[part1->type][part2->type].pswitchINV;
            e *= e;
            e = ( 1 - e + e*e*(1.0/3.0) - e*e*e*(2.0/45.0) + e*e*e*e*(1.0/315.0) - e*e*e*e*e*(2.0/14175.0)
                    + e*e*e*e*e*e*(2.0/467775.0) - e*e*e*e*e*e*e*(4.0/42567525) + e*e*e*e*e*e*e*e*(1.0/638512875) ) * -topo.ia_params[part1->type][part2->type].epsilon;
        } else {
            e = -topo.ia_params[part1->type][part2->type].epsilon;
        }

        // WCA repulsion
        if (dist > topo.ia_params[part1->type][part2->type].rcutwca) {
            return e;
        }
        return topo.ia_params[part1->type][part2->type].A * pow(dist, -12) - topo.ia_params[part1->type][part2->type].B * pow(dist, -6);
    }
};

class WcaCos2Full : public EPotential {
public:
    double operator() (double dist, const Particle* part1, const Particle* part2) override {
        double e = 0.0;
        if(dist > topo.ia_params[part1->type][part2->type].rcut ||
                (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
                topo.ia_params[part1->type][part2->type].exclude)
            return 0.0;

        if(dist > topo.ia_params[part1->type][part2->type].pdis) {
            // full precision cos
            e = cos(PIH*(dist - topo.ia_params[part1->type][part2->type].pdis) * topo.ia_params[part1->type][part2->type].pswitchINV);
            e *= -e * topo.ia_params[part1->type][part2->type].epsilon;
        } else {
            e = -topo.ia_params[part1->type][part2->type].epsilon;
        }

        // WCA repulsion
        if (dist > topo.ia_params[part1->type][part2->type].rcutwca) {
            return e;
        }
        return topo.ia_params[part1->type][part2->type].A * pow(dist, -12) - topo.ia_params[part1->type][part2->type].B * pow(dist, -6);
    }
};




class Common {
public:
    inline double closestDist(const Vector& r_cm, const Vector& dir1, const Vector& dir2, const Ia_param& iaParam) {
        Vector distvec = minDistSegments(dir1, dir2, iaParam.half_len[0], iaParam.half_len[1], r_cm);
        return DOT(distvec,distvec);
    }

    Vector minDistSegments(const Vector &segA, const Vector &segB, double halfl1, double halfl2, const Vector &r_cm);

    inline double atrE(const Ia_param& iaParam, const Vector& p1Dir, const Vector& p2Dir, const Patch& p1P, const Patch& p2P, const Vector& r_cm,
                int patchnum1, int patchnum2, double& S1, double& S2, double& T1, double& T2) {
        Vector vec1, vec2, vec_intrs, vec_mindist;
        double atrenergy = 0, ndist, a, paral, f1, f2;

        // 3 - scaling function1: dependence on the length of intersetions
        double v1 = fabs(S1-S2);
        double v2 = fabs(T1-T2);
        double f0 = 0.5*(v1+v2);

        // 4a - with two intersection pieces calculate vector between their CM -this is for angular orientation
        vec1 = ((S1+S2)*0.5) * p1Dir;
        vec2 = ((T1+T2)*0.5) * p2Dir;
        vec_intrs.x = vec2.x - vec1.x - r_cm.x; //vec_intrs should be from sc1 to sc2
        vec_intrs.y = vec2.y - vec1.y - r_cm.y;
        vec_intrs.z = vec2.z - vec1.z - r_cm.z;

        // 4b - calculate closest distance attractive energy from it
        vec_mindist = minDistSegments(p1Dir, p2Dir, v1, v2, vec_intrs);
        ndist = sqrt(DOT(vec_mindist, vec_mindist));

        //
        // POTENTIAL FUNCTION
        //
        if (ndist < iaParam.pdis)
            atrenergy = -iaParam.epsilon;
        else {
            atrenergy = cos(PIH*(ndist - iaParam.pdis) / iaParam.pswitch);
            atrenergy *= -atrenergy * iaParam.epsilon;
        }

        //
        // 5 - scaling function2: angular dependence of patch1
        //
        vec1 = vec_intrs;
        vec1 = vec1.perpProject(p1Dir);
        a = DOT(vec1, p1P.dir) / vec1.size();
        f1 = fanglScale(a, iaParam.pcangl[0+2*patchnum1], iaParam.pcanglsw[0+2*patchnum1]);

        //
        // 6 - scaling function3: angular dependence of patch2
        //
        vec1 = -1.0 * vec_intrs;
        vec1 = vec1.perpProject(p2Dir);
        a = DOT(vec1, p2P.dir) / vec1.size();
        f2 = fanglScale(a, iaParam.pcangl[1+2*patchnum2], iaParam.pcanglsw[1+2*patchnum2]);

        //
        // 7 - add scaling increased if particles are parallel or antiparallel
        //
        paral = 1.0;
        if( iaParam.parallel != 0.0)
            paral = scparallel(iaParam.parallel, p1Dir, p2Dir);

        //7- put it all together, CPSC + PSC -> dont have PARAL E
        atrenergy *= f0 * f1 * f2 * paral;

        return atrenergy;
    }




    int cpscIntersect(const Vector& p1Dir, const Vector& p2Dir, const Patch& p1P, const Vector& r_cm,
                      double& in1, double& in2, const double pcanglsw, double rcutSq, double halfl1, double halfl2) {
        double a, b, c, d, x1, x2, disti = 0.0, ti;
        Vector vec1, vec2;
        double intersections[2] = {0};

        // 1- do intersections of spherocylinder2 with patch of spherocylinder1 at cut distance C
        // 1a- test intersection with half planes of patch and look how far they are from spherocylinder. If closer then C  we got itersection

        // plane1, find intersections of part2 with plane by par1 and part1->patchsides[0]
        if(findIntersectPlaneUni(p1Dir, p2Dir, halfl2, r_cm, p1P.sides[0], pcanglsw, ti, c, d)) {
            if (d <= 0) {
                disti = c*c; // is inside cylinder
                if (disti <= rcutSq) {// the intersection is outside sc
                    intersections[0] = ti;
                }
            }
        }

        // plane2, find intersections of part2 with plane by par1 and part1->patchsides[1]
        if( findIntersectPlaneUni(p1Dir, p2Dir, halfl2, r_cm, p1P.sides[1], pcanglsw, ti, c, d) ) {
            if (d <= 0) {
                disti = c*c; // is inside cylinder
                if (disti <= rcutSq) { // the intersection is outside sc
                    if(intersections[0] == 0.0) {
                        intersections[0] = ti;
                    } else {
                        if(ti != intersections[0]) {
                            intersections[1] = ti;
                        }
                    }
                }
            }
        }

        if ( intersections[1] != 0.0 ) {
            if(pcanglsw < 0)
                fprintf (stderr, "ERROR: Patch is larger than 180 degrees and we are getting two segments - this hasnot been programed yet.\n\n");

            in1 = intersections[0];
            in2 = intersections[1];

            return 2;
        } else {
            // 1b- test intersection with cylinder - it is at distance C
            testIntrC(p1Dir, p2Dir, p1P, r_cm, pcanglsw, rcutSq, halfl1, halfl2, intersections);
        }

        // 1c- test intersection with plates at the end - it is at distace C and in wedge
        if ( intersections[1] == 0.0 )  {
            a =  DOT(p1Dir, p2Dir);
            if (a == 0.0)
                ; // there is no intersection plane and sc are paralel
            else {
                // plane cap1
                vec1.x= r_cm.x + halfl1*p1Dir.x;
                vec1.y= r_cm.y + halfl1*p1Dir.y;
                vec1.z= r_cm.z + halfl1*p1Dir.z;
                x1 = DOT(p1Dir,vec1)/a; // parameter on line of SC2 determining intersection
                if ((x1 > halfl2 ) || (x1 < -halfl2))
                    ; // there is no intersection plane sc is too short
                else {
                    vec2.x = x1*p2Dir.x - vec1.x; // vector from ENDPOINT to intersection point
                    vec2.y = x1*p2Dir.y - vec1.y;
                    vec2.z = x1*p2Dir.z - vec1.z;
                    b = DOT (vec2, vec2);
                    if (b > rcutSq)
                        ;   // the intersection is outside sc
                    else {
                        testIntrPatch(p1Dir,p1P.dir,vec2, pcanglsw,x1,intersections);
                    }
                }
                //		    printf ("plane cap1 %d %f\n", intrs, x1);
                //plane cap2
                vec1.x= r_cm.x - halfl1*p1Dir.x;
                vec1.y= r_cm.y - halfl1*p1Dir.y;
                vec1.z= r_cm.z - halfl1*p1Dir.z;
                x2 = DOT(p1Dir,vec1)/a; // parameter on line of SC2 determining intersection
                if ((x2  > halfl2 ) || (x2 < -halfl2))
                    ; // there is no intersection plane sc is too short
                else {
                    vec2.x = x2*p2Dir.x - vec1.x; // vector from ENDPOINT to intersection point
                    vec2.y = x2*p2Dir.y - vec1.y;
                    vec2.z = x2*p2Dir.z - vec1.z;
                    b = DOT (vec2, vec2);
                    if (b > rcutSq)
                        ;   // the intersection is outside sc
                    else {
                        testIntrPatch(p1Dir,p1P.dir,vec2, pcanglsw,x2,intersections);
                    }
                }
            }
        } else {
            in1 = intersections[0];
            in2 = intersections[1];

            return 2;
        }

        // 1d- if there is only one itersection shperocylinder ends within patch wedge set as second intersection end inside patch
        if ( intersections[1] == 0.0 )  {
            // whole spherocylinder is in or all out if intrs ==0
            vec1.x = p2Dir.x*halfl2 - r_cm.x;
            vec1.y = p2Dir.y*halfl2 - r_cm.y;
            vec1.z = p2Dir.z*halfl2 - r_cm.z;
            // vector from CM of sc1 to end of sc2
            // check is is inside sc1
            a=DOT(vec1,p1Dir);
            vec2.x = vec1.x - p1Dir.x*a;
            vec2.y = vec1.y - p1Dir.y*a;
            vec2.z = vec1.z - p1Dir.z*a;
            b=DOT(vec2,vec2);
            d = fabs(a)-halfl1;
            if ( d <= 0) { // is in cylindrical part
                // c is distance squared from line or end to test if is inside sc
                if (b < rcutSq)
                    testIntrPatch(p1Dir,p1P.dir,vec1, pcanglsw,halfl2,intersections);
            }
            if ( intersections[1] == 0.0 ) {
                vec1.x = -p2Dir.x*halfl2 - r_cm.x;
                vec1.y = -p2Dir.y*halfl2 - r_cm.y;
                vec1.z = -p2Dir.z*halfl2 - r_cm.z;
                // check is is inside sc1
                a=DOT(vec1,p1Dir);
                vec2.x = vec1.x - p1Dir.x*a;
                vec2.y = vec1.y - p1Dir.y*a;
                vec2.z = vec1.z - p1Dir.z*a;
                b=DOT(vec2,vec2);
                d = fabs(a) -halfl1;
                if (d <= 0) {
                    // c is distance squared from line or end to test if is inside sc
                    if (b < rcutSq)
                        testIntrPatch(p1Dir,p1P.dir,vec1, pcanglsw, -1.0*halfl2,intersections);
                }
            }
        }
        in1 = intersections[0];
        in2 = intersections[1];

        return (intersections[1] == 0.0) ? 0 : 2;
    }

    int pscIntersect(const Vector& p1Dir, const Vector& p2Dir, const Patch& p1P, const Vector& r_cm,
                     double& in1, double& in2, const double pcanglsw, double rcutSq, double halfl1, double halfl2) {
        double c, d;
        double ti, disti;
        Vector vec1, vec2;
        double intersections[2] = {0};

        //1- do intersections of spherocylinder2 with patch of spherocylinder1 at cut distance C
        //1a- test intersection with half planes of patch and look how far they are from spherocylinder. If closer then C  we got itersection
        // plane1 // find intersections of part2 with plane by par1 and patchsides[0]
        if(findIntersectPlaneUni(p1Dir, p2Dir, halfl2, r_cm, p1P.sides[0], pcanglsw, ti, c, d)) {
            if (d <= 0)
                disti = c*c; // is inside cylinder
            else
                disti = d*d + c*c; // is inside patch

            if (disti <= rcutSq) {// the intersection is outside sc
                intersections[0] = ti;
            }
        }

        // plane2 // find intersections of part2 with plane by par1 and patchsides[1]
        if( findIntersectPlaneUni(p1Dir, p2Dir, halfl2, r_cm, p1P.sides[1], pcanglsw, ti, c, d) ) {
            if (d <= 0)
                disti = c*c; // is inside cylinder
            else
                disti = d*d + c*c; // is inside patch

            if (disti <= rcutSq) { // the intersection is outside sc
                if(intersections[0] == 0.0) {
                    intersections[0] = ti;
                } else {
                    if(ti != intersections[0]) {
                        intersections[1] = ti;
                    }
                }
            }
        }

        if ( intersections[1] != 0.0 ) {
            if(pcanglsw < 0)
                fprintf (stderr, "ERROR: Patch is larger than 180 degrees and we are getting two segments - this hasnot been programed yet.\n\n");

            in1 = intersections[0];
            in2 = intersections[1];

            return 2;
        } else {
            // 1b - test intersection with cylinder - it is at distance C
            testIntrC(p1Dir, p2Dir, p1P, r_cm, pcanglsw, rcutSq, halfl1, halfl2, intersections);
        }

        // 1c - test intersection with spheres at the end - it is at distace C
        if ( intersections[1] == 0.0 )  {
            //centers of spheres, relative to the CM of sc2
            vec1.x =  p1Dir.x * halfl1 - r_cm.x;
            vec1.y =  p1Dir.y * halfl1 - r_cm.y;
            vec1.z =  p1Dir.z * halfl1 - r_cm.z;
            vec2.x = -p1Dir.x * halfl1 - r_cm.x;
            vec2.y = -p1Dir.y * halfl1 - r_cm.y;
            vec2.z = -p1Dir.z * halfl1 - r_cm.z;

            calcIntersections(p1Dir, p2Dir, p1P, r_cm, intersections, pcanglsw, halfl1, halfl2, 2.0*DOT(vec1, p2Dir), DOT(vec1,vec1) - rcutSq);
            calcIntersections(p1Dir, p2Dir, p1P, r_cm, intersections, pcanglsw, halfl1, halfl2, 2.0*DOT(vec2, p2Dir), DOT(vec2,vec2) - rcutSq);
        } else {
            in1 = intersections[0];
            in2 = intersections[1];

            return 2;
        }

        // 1d - if there is only one itersection shperocylinder ends within patch wedge set as second intersection end inside patch
        if ( intersections[1] == 0.0 )  {

            testIntrA(p1Dir, p2Dir, p1P, r_cm, pcanglsw, rcutSq, halfl1, halfl2, intersections);

            if ( intersections[1] == 0.0 ) {
                testIntrA(p1Dir, p2Dir, p1P, r_cm, pcanglsw, rcutSq, halfl1, -halfl2, intersections);
            }
        }
        in1 = intersections[0];
        in2 = intersections[1];

        return (intersections[1] == 0.0) ? 0 : 2;
    }

private:


    void calcIntersections(const Vector& p1Dir, const Vector& p2Dir, const Patch& p1P, const Vector& r_cm,
             double intersections[], double pcanglsw, double halfl1, double halfl2, double b, double c) {
        int intrs=0;
        Vector vec1;
        double x1,e;

        double d = b*b-4*c;

        if (d >= 0) { //if d<0 there are no intersections

            d = sqrt(d);
            x1 = (-b + d)*0.5; //parameter on line of SC2 determining intersection

            if ((x1 >=halfl2) || (x1 <= -halfl2)) {
                intrs+=0; //intersection is outside sc2
            } else {
                vec1.x = p2Dir.x*x1 - r_cm.x;
                vec1.y = p2Dir.y*x1 - r_cm.y;
                vec1.z = p2Dir.z*x1 - r_cm.z;
                e = DOT(p1Dir, vec1);
                if ((e >= halfl1) || (e <= -halfl1)) { //if not intersection is inside sc1
                    testIntrPatch(p1Dir, p1P.dir, vec1, pcanglsw, x1, intersections);
                }
            }
            if ( d > 0) {

                x1= (-b - d)*0.5; //parameter on line of SC2 determining intersection

                if ((x1 >=halfl2) || (x1 <= -halfl2))
                    intrs+=0; //intersection is outside sc2
                else {
                    vec1.x = p2Dir.x*x1 - r_cm.x;
                    vec1.y = p2Dir.y*x1 - r_cm.y;
                    vec1.z = p2Dir.z*x1 - r_cm.z;
                    e = DOT(p1Dir, vec1);
                    if ((e >=halfl1) || (e <= -halfl1)) { //if not intersection is inside sc1
                        testIntrPatch(p1Dir, p1P.dir, vec1, pcanglsw, x1, intersections);
                    }
                }
            }
        }
        return;
    }

    inline void testIntrA(const Vector& p1Dir, const Vector& p2Dir, const Patch& p1P, const Vector& r_cm,
                          const double pcanglsw, double rcutSq, double halfl1, double halfl2, double intersections[]) {
        double a, b, c, d;
        Vector vec1, vec2;

        //whole spherocylinder is in or all out if intrs ==0
        vec1.x = p2Dir.x*halfl2 - r_cm.x;
        vec1.y = p2Dir.y*halfl2 - r_cm.y;
        vec1.z = p2Dir.z*halfl2 - r_cm.z;

        //vector from CM of sc1 to end of sc2, check if is inside sc1
        a=DOT(vec1,p1Dir);
        vec2.x = vec1.x - p1Dir.x*a;
        vec2.y = vec1.y - p1Dir.y*a;
        vec2.z = vec1.z - p1Dir.z*a;
        b=DOT(vec2,vec2);
        d = fabs(a)-halfl1;

        if ( d <= 0)
            c = b; //is inside cylindrical part
        else
            c = d*d + b; //is inside caps

        //c is distance squared from line or end to test if is inside sc
        if (c < rcutSq)
            testIntrPatch(p1Dir, p1P.dir, vec1, pcanglsw, halfl2, intersections);
    }

    inline void testIntrC(const Vector& p1Dir, const Vector& p2Dir, const Patch& p1P, const Vector& r_cm,
                  const double pcanglsw, double rcutSq, double halfl1, double halfl2, double intersections[]) {
        Vector vec1,vec2;
        double a,b,c,d,e,x1;

        vec1 = vecCrossProduct(-1.0 * r_cm,p1Dir);
        vec2 = vecCrossProduct(p2Dir,p1Dir);

        a = DOT(vec2,vec2);
        b = 2*DOT(vec1,vec2);
        c = -rcutSq + DOT(vec1,vec1);

        d = b*b - 4*a*c;

        if ( d >= 0) { //there is intersection with infinite cylinder

            d = sqrt(d);
            a = 0.5/a;
            x1 = (-b + d) * a;//parameter on line of SC2 determining intersection

            if ( (x1 >= halfl2) || (x1 <= -halfl2) ) {
                ; //intersection is outside sc2
            } else {
                // vectors from center os sc1 to intersection with infinite cylinder
                vec1.x = p2Dir.x*x1-r_cm.x;
                vec1.y = p2Dir.y*x1-r_cm.y;
                vec1.z = p2Dir.z*x1-r_cm.z;
                e = DOT(p1Dir, vec1);
                if ((e >=halfl1) || (e <= -halfl1))
                    ; //intersection is outside sc1
                else {
                    testIntrPatch(p1Dir, p1P.dir, vec1, pcanglsw, x1, intersections);
                }
            }

            if ( d > 0 ){
                x1 = (-b - d) * a;//parameter on line of SC2 determining intersection

                if ( (x1 >= halfl2) || (x1 <= -halfl2) ) {
                    ; //intersection is outside sc2
                } else {
                    // vectors from center os sc1 to intersection with infinite cylinder
                    vec1.x = p2Dir.x*x1-r_cm.x;
                    vec1.y = p2Dir.y*x1-r_cm.y;
                    vec1.z = p2Dir.z*x1-r_cm.z;
                    e = DOT(p1Dir, vec1);
                    if ((e >=halfl1) || (e <= -halfl1))
                        ; //intersection is outside sc1
                    else {
                        testIntrPatch(p1Dir, p1P.dir, vec1, pcanglsw, x1, intersections);
                    }
                }
            }
        }
        return;
    }



    inline void testIntrPatch(const Vector& dir, const Vector& patchdir, Vector vec, double cospatch,
                                            double ti, double intersections[]) { // test if we have intersection

        vec = vec.perpProject(dir); // do projection to patch plane
        if (DOT(patchdir,vec) / vec.size() >= cospatch) { // test angle distance from patch
            if(intersections[0] == 0) {
                intersections[0] = ti;
                return;
            }
            if(intersections[1] == 0 && intersections[0] != ti) {
                intersections[1] = ti;
                return;
            }
        }
        return;
    }


    inline double scparallel(double epsilonparallel, const Vector& dir1, const Vector& dir2){
        double cosa=DOT(dir1,dir2);

        if ((epsilonparallel>0 && cosa>0) || (epsilonparallel<0 && cosa<0))
            return 1.0 + epsilonparallel*cosa;
        else
            return 1.0;
    }


    bool findIntersectPlaneUni(const Vector& dirA, const Vector& dirB, double halfl,
                                                 const Vector& r_cm, Vector w_vec, double cospatch, double& ti, double& c, double& d) {
        assert(fabs(w_vec.size() - 1.0) < 1e-13 || !(cout << "w_vec isnt normalised" << endl));

        double a;
        Vector d_vec, nplane = vecCrossProduct(dirA, w_vec);
        a =  DOT(nplane, dirB);
        c = 1.0, d = 1.0;

        if (a == 0.0) {
            return false; // there is no intersection plane and sc are paralel
        } else {
            ti = DOT(nplane,r_cm)/a;
            if ((ti  > halfl ) || (ti < -halfl)) {
                return false; // there is no intersection plane sc is too short
            } else {
                d_vec.x = ti * dirB.x - r_cm.x; // vector from intersection point to CM
                d_vec.y = ti * dirB.y - r_cm.y;
                d_vec.z = ti * dirB.z - r_cm.z;

                c = DOT (d_vec, w_vec);

                if ( c * cospatch < 0)  {
                    return false; // the intersection in plane is on other side of patch
                } else {
                    d = fabs(DOT (d_vec, dirA)) - halfl;
                    return true;
                }
            }
        }
        return false;
    }
};






class Psc : public EPatch, public Common {
public:
    double operator() (const Ia_param& iaParam, const Vector& p1Dir, const Vector& p2Dir, const Patch p1P, const Patch p2P, const Vector& r_cm,
                       int patchnum1, int patchnum2) override {
        double T1, T2, S1, S2;
        Vector vec1;

        // 1 - do intersections of spherocylinder2 with patch of spherocylinder1 at cut distance C
        if( 2 > pscIntersect( p1Dir, p2Dir, p1P, r_cm, T1, T2, iaParam.pcanglsw[2*patchnum1], iaParam.rcutSq, iaParam.half_len[0], iaParam.half_len[1]) ) // No intersection
            return 0.0; // sc is all outside patch, attractive energy is 0

        // 2 - now do the same oposite way psc1 in patch of psc2
        vec1 = -1.0*r_cm;
        if( 2 > pscIntersect( p2Dir, p1Dir, p2P, vec1, S1, S2, iaParam.pcanglsw[2*patchnum2+1], iaParam.rcutSq, iaParam.half_len[1], iaParam.half_len[0]))
            return 0.0; //sc is all outside patch, attractive energy is 0

        return atrE(iaParam, p1Dir, p2Dir, p1P, p2P, r_cm, patchnum1, patchnum2, S1, S2, T1, T2);
    }

    void getGeoTypes(const Ia_param& iaParam, bool& firstCH, bool& secondCH, bool& firstT, bool& secondT) override {
        firstCH = false; secondCH = false; firstT = false; secondT = false;

        if ( (iaParam.geotype[0] == CHPSC) || (iaParam.geotype[0] == TCHPSC) )
            firstCH = true;
        if ( (iaParam.geotype[1] == CHPSC) || (iaParam.geotype[1] == TCHPSC) )
            secondCH = true;
        if ( (iaParam.geotype[0] == TPSC)  || (iaParam.geotype[0] == TCHPSC) )
            firstT = true;
        if ( (iaParam.geotype[1] == TPSC)  || (iaParam.geotype[1] == TCHPSC) )
            secondT = true;
    }
};

class CPsc : public EPatch, public Common {
public:
    double operator() (const Ia_param& iaParam, const Vector& p1Dir, const Vector& p2Dir, const Patch p1P, const Patch p2P, const Vector& r_cm,
                       int patchnum1, int patchnum2) override {
        double T1, T2, S1, S2;
        Vector vec1;

        // 1 - do intersections of spherocylinder2 with patch of spherocylinder1 at cut distance C
        if( 2 > cpscIntersect( p1Dir, p2Dir, p1P, r_cm, T1, T2, iaParam.pcanglsw[2*patchnum1], iaParam.rcutSq, iaParam.half_len[0], iaParam.half_len[1])) // No intersection
            return 0.0; // sc is all outside patch, attractive energy is 0

        // 2 - now do the same oposite way psc1 in patch of psc2
        vec1 = -1.0*r_cm;
        if( 2 > cpscIntersect(   p2Dir, p1Dir, p2P, vec1, S1, S2, iaParam.pcanglsw[2*patchnum2+1], iaParam.rcutSq, iaParam.half_len[1], iaParam.half_len[0]) )
            return 0.0; //sc is all outside patch, attractive energy is 0

        return atrE(iaParam, p1Dir, p2Dir, p1P, p2P, r_cm, patchnum1, patchnum2, S1, S2, T1, T2);
    }

    void getGeoTypes(const Ia_param& iaParam, bool& firstCH, bool& secondCH, bool& firstT, bool& secondT) override {
        firstCH = false; secondCH = false; firstT = false; secondT = false;

        if ( (iaParam.geotype[0] == CHCPSC) || (iaParam.geotype[0] == TCHCPSC) )
            firstCH = true;
        if ( (iaParam.geotype[1] == CHCPSC) || (iaParam.geotype[1] == TCHCPSC) )
            secondCH = true;
        if ( (iaParam.geotype[0] == TCPSC)  || (iaParam.geotype[0] == TCHCPSC) )
            firstT = true;
        if ( (iaParam.geotype[1] == TCPSC)  || (iaParam.geotype[1] == TCHCPSC) )
            secondT = true;
    }
};

class Scn : public EPatch, public Common {
public:
    double operator() (const Ia_param& iaParam, const Vector& p1Dir, const Vector& p2Dir, const Patch p1P, const Patch p2P, const Vector& r_cm,
                       int patchnum1, int patchnum2) override {
        return 0.0;
    }

    void getGeoTypes(const Ia_param& iaParam, bool& firstCH, bool& secondCH, bool& firstT, bool& secondT) override {
        firstCH = false; secondCH = false; firstT = false; secondT = false;
    }
};

class Sca : public EPatch, public Common {
public:
    double operator() (const Ia_param& iaParam, const Vector& p1Dir, const Vector& p2Dir, const Patch p1P, const Patch p2P, const Vector& r_cm,
                       int patchnum1, int patchnum2) override {
        double dist = closestDist(r_cm, p1Dir, p2Dir, iaParam); // redundant, but in practise SCA never used
        dist = sqrt(dist);
        if (dist > iaParam.rcutwca)
            return 0.0;
        return iaParam.A * pow(dist, -12) - iaParam.B * pow(dist, -6) + iaParam.epsilon;
    }

    void getGeoTypes(const Ia_param& iaParam, bool& firstCH, bool& secondCH, bool& firstT, bool& secondT) override {
        firstCH = false; secondCH = false; firstT = false; secondT = false;
    }
};

class PscCPsc : public EPatch, public Common {
public:
    double operator() (const Ia_param& iaParam, const Vector& p1Dir, const Vector& p2Dir, const Patch p1P, const Patch p2P, const Vector& r_cm,
                       int patchnum1, int patchnum2) override {
        double T1, T2, S1, S2;
        Vector vec1 = -1.0 * r_cm;;

        bool first =false;
        if ( (iaParam.geotype[0] == PSC)||(iaParam.geotype[0] == CHPSC)||(iaParam.geotype[0] == TPSC)||(iaParam.geotype[0] == TCHPSC) ){
            first = true;
        }

        // 1- do intersections of spherocylinder2 with patch of spherocylinder1 at. cut distance C
        // 2- now do the same oposite way psc1 in patch of psc2
        if (first) {
            if( 2 > pscIntersect( p1Dir, p2Dir, p1P, r_cm, T1, T2, iaParam.pcanglsw[2*patchnum1], iaParam.rcutSq, iaParam.half_len[0], iaParam.half_len[1]) )
                return 0.0;
            if( 2 > cpscIntersect(p2Dir, p1Dir, p2P, vec1, S1, S2, iaParam.pcanglsw[2*patchnum2+1], iaParam.rcutSq, iaParam.half_len[1], iaParam.half_len[0]))
                return 0.0;
       } else {
            if( 2 > cpscIntersect(   p1Dir, p2Dir, p1P, r_cm, T1, T2, iaParam.pcanglsw[2*patchnum1], iaParam.rcutSq, iaParam.half_len[0], iaParam.half_len[1]))
                return 0.0;
            if( 2 > pscIntersect(p2Dir, p1Dir, p2P, vec1, S1, S2, iaParam.pcanglsw[2*patchnum2+1], iaParam.rcutSq, iaParam.half_len[1], iaParam.half_len[0]))
                return 0.0;
        }

        return atrE(iaParam, p1Dir, p2Dir, p1P, p2P, r_cm, patchnum1, patchnum2, S1, S2, T1, T2);
    }

    void getGeoTypes(const Ia_param& iaParam, bool& firstCH, bool& secondCH, bool& firstT, bool& secondT) override {
        firstCH = false; secondCH = false; firstT = false; secondT = false;

        if ((iaParam.geotype[0] == CHPSC) || (iaParam.geotype[0] == CHCPSC)||
          (iaParam.geotype[0] == TCHPSC) || (iaParam.geotype[0] == TCHCPSC) )
            firstCH = true;
        if ((iaParam.geotype[1] == CHPSC) || (iaParam.geotype[1] == CHCPSC)||
          (iaParam.geotype[1] == TCHPSC) || (iaParam.geotype[1] == TCHCPSC) )
            secondCH = true;
        if ( (iaParam.geotype[0] == TCPSC) || (iaParam.geotype[0] == TCHCPSC) ||
             (iaParam.geotype[0] == TPSC) || (iaParam.geotype[0] == TCHPSC) )
            firstT = true;
        if ( (iaParam.geotype[1] == TCPSC) || (iaParam.geotype[1] == TCHCPSC) ||
             (iaParam.geotype[1] == TPSC) || (iaParam.geotype[1] == TCHPSC) )
            secondT = true;
    }
};




class ScaSpa : public EPatchToSphere  {
public:
    double operator() (double dist, double contt, Vector& distvec, const Ia_param& iaParam, const Vector& p1Dir, const Patch p1P) override {

        double atrenergy, b, f0, halfl;

        if (dist < iaParam.pdis)
            atrenergy = -iaParam.epsilon;
        else {
            atrenergy = cos(PIH*(dist - iaParam.pdis) / iaParam.pswitch);
            atrenergy *= -atrenergy * iaParam.epsilon ;
        }

        //scaling function for the length of spherocylinder within cutoff
        halfl = iaParam.half_len[0];
        b = sqrt(iaParam.rcutSq - dist*dist);
        if ( contt + b > halfl )
            f0 = halfl;
        else
            f0 = contt + b;
        if ( contt - b < -halfl )
            f0 -= -halfl;
        else
            f0 -= contt - b;

        atrenergy *= f0;

        return atrenergy;
    }
};

class PscSpa : public EPatchToSphere  {
public:
    double operator() (double dist, double contt, Vector& distvec, const Ia_param& iaParam, const Vector& p1Dir, const Patch p1P) override {

        double atrenergy, a, b, f0, halfl;
        Vector vec1;
        int which;

        if (dist < iaParam.pdis)
            atrenergy = -iaParam.epsilon;
        else {
            atrenergy = cos(PIH*(dist - iaParam.pdis) / iaParam.pswitch);
            atrenergy *= -atrenergy * iaParam.epsilon ;
        }

        // scaling function: angular dependence of patch1
        which = 0;
        vec1 = distvec.perpProject(p1Dir);
        a = DOT(vec1, p1P.dir) / vec1.size();
        halfl = iaParam.half_len[0];

        //scaling function for the length of spherocylinder within cutoff
        b = sqrt(iaParam.rcutSq - dist*dist);
        if ( contt + b > halfl )
            f0 = halfl;
        else
            f0 = contt + b;
        if ( contt - b < -halfl )
            f0 -= -halfl;
        else
            f0 -= contt - b;

        atrenergy *= fanglScale(a, iaParam.pcangl[which], iaParam.pcanglsw[which])*(f0);

        return atrenergy;
    }

    void getGeoTypes(const Ia_param& iaParam, bool& chiral, bool& sec) override {
        chiral=false; sec=false;

        if ( (iaParam.geotype[0] == CHPSC) || (iaParam.geotype[0] == TCHPSC) )
            chiral = true;

        if ( (iaParam.geotype[0] == TPSC) || (iaParam.geotype[0] == TCHPSC) )
            sec = true;
    }
};

class CPscSpa : public PscSpa  {
public:
    void getGeoTypes(const Ia_param& iaParam, bool& chiral, bool& sec) override {
        chiral=false; sec=false;

        if ( (iaParam.geotype[0] == CHCPSC) || (iaParam.geotype[0] == TCHCPSC) )
            chiral = true;

        if ( (iaParam.geotype[0] == TCPSC) || (iaParam.geotype[0] == TCHCPSC) )
            sec = true;
    }

    virtual bool isFar(double contt, double halfl) override {
        return (contt > halfl) || (contt < -halfl);
    }
};




template <typename PatchE, typename RepE, typename BondE, typename AngleE>
class MixSpSc : public EBasic {
public:
    PatchE patchE;
    RepE repE;
    BondE bondE;
    AngleE angleE;

    double operator() (double dist, const Vector& r_cm, const Particle* part1, const Particle* part2, const ConList* conlist) {

        double atrenergy, repenergy = 0.0, contt=0.0;
        Vector distvec;
        bool isp1Spc = (topo.ia_params[part1->type][part2->type].geotype[0] < SP);
        const Particle* spc = (isp1Spc) ? part1 : part2;
        const Ia_param& iaParam = (isp1Spc) ? topo.ia_params[part1->type][part2->type] : topo.ia_params[part2->type][part1->type];

        double distSq = closestDist((isp1Spc) ? std::move(r_cm) : -1.0*r_cm, spc->dir, iaParam.half_len[0], contt, distvec);

        // WCA
        if (distSq < iaParam.rcutwcaSq)
            repenergy = iaParam.epsilon + iaParam.A * pow(distSq, -6)- iaParam.B * pow(distSq, -3);

        if ( ( distSq > iaParam.rcutSq ) || (iaParam.epsilon == 0.0 ) || iaParam.exclude || patchE.isFar(contt, iaParam.half_len[0]) ) {
            atrenergy = 0.0;
        } else {
            bool chiral=false, sec=false;

            patchE.getGeoTypes(iaParam, chiral, sec);

            if (chiral) {
                distSq = closestDist((isp1Spc) ? std::move(r_cm) : -1.0*r_cm, (chiral) ? spc->chdir[0] : spc->dir, iaParam.half_len[0], contt, distvec);
            }

            dist = sqrt(distSq);
            atrenergy = patchE(dist, contt, distvec, iaParam, (chiral) ? spc->chdir[0] : spc->dir,
                                           (Patch(spc->patchdir[0], spc->patchsides[0], spc->patchsides[1])) );

            //addition of interaction of second patches
            if(sec) {
                distSq = closestDist((isp1Spc) ? std::move(r_cm) : -1.0*r_cm, (chiral) ? spc->chdir[1] : spc->dir, iaParam.half_len[0], contt, distvec);
                dist = sqrt(distSq);
                atrenergy += patchE(dist, contt, distvec, iaParam, (chiral) ? spc->chdir[1] : spc->dir,
                                               (Patch(spc->patchdir[1], spc->patchsides[2], spc->patchsides[3])) );
            }
        }
        return repenergy+atrenergy;
    }
};

/**
 *  @brief Inteactions of spheres
 */
template <typename PotentialE, typename BondE, typename AngleE>
class Sphere : public EBasic {
    PotentialE potE;
    BondE bondE;
    AngleE angleE;
public:
    double operator() (double dist, const Vector& r_cm, const Particle* part1, const Particle* part2, const ConList* conlist) {
        double e = 0.0;
        if(conlist != nullptr) {
            e = bondE(dist, part1, part2, conlist);
            e += angleE(dist, part1, part2, conlist);
        }
        return e + potE(dist, part1, part2);
    }
};

/**
 *  @brief Inteactions of spherocylinders
 */
template <typename PatchE, typename RepE, typename BondE, typename AngleE>
class SpheroCylinder : public EBasic {
public:
    PatchE patchE;
    RepE repE;
    BondE bondE;
    AngleE angleE;

    double operator() (double dist, const Vector& r_cm, const Particle* part1, const Particle* part2, const ConList* conlist) {
        double abE = 0.0, distSq = 0.0, atrenergy = 0.0, repenergy = 0.0;

        if(conlist != nullptr) {
            abE = bondE(dist, part1, part2, conlist);
            abE += angleE(dist, part1, part2, conlist);
        }

        distSq = patchE.closestDist(r_cm, part1->dir, part2->dir, topo.ia_params[part1->type][part2->type]);
        repenergy = repE(distSq, part1, part2); // WCA repulsion from SC shape

        if ( ( distSq > topo.ia_params[part1->type][part2->type].rcutSq ) ||
             (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
             topo.ia_params[part1->type][part2->type].exclude ) { // cutoff or not interacting
            atrenergy = 0.0;
        } else {
            bool firstCH, secondCH, firstT, secondT;
            patchE.getGeoTypes(topo.ia_params[part1->type][part2->type], firstCH, secondCH, firstT, secondT);

            atrenergy = patchE(topo.ia_params[part1->type][part2->type],
                    (firstCH) ? part1->chdir[0] : part1->dir,
                    (secondCH) ? part2->chdir[0] : part2->dir,
                    Patch(part1->patchdir[0], part1->patchsides[0], part1->patchsides[1]),
                    Patch(part2->patchdir[0], part2->patchsides[0], part2->patchsides[1]), r_cm, 0,0);

            // addition of interaction of second patches
            if(firstT) { // part1 has second patch
                atrenergy += patchE(topo.ia_params[part1->type][part2->type],
                        (firstCH) ? part1->chdir[1] : part1->dir,
                        (secondCH) ? part2->chdir[0] : part2->dir,
                    Patch(part1->patchdir[1], part1->patchsides[2], part1->patchsides[3]),
                    Patch(part2->patchdir[0], part2->patchsides[0], part2->patchsides[1]), r_cm, 1,0);
            }

            if(secondT) { // part2 has second patch
                atrenergy += patchE(topo.ia_params[part1->type][part2->type],
                        (firstCH) ? part1->chdir[0] : part1->dir,
                        (secondCH) ? part2->chdir[1] : part2->dir,
                    Patch(part1->patchdir[0], part1->patchsides[0], part1->patchsides[1]),
                    Patch(part2->patchdir[1], part2->patchsides[2], part2->patchsides[3]), r_cm, 0,1);
            }

            if(firstT && secondT) { // part1 and part2 has second patch
                atrenergy += patchE(topo.ia_params[part1->type][part2->type],
                        (firstCH) ? part1->chdir[1] : part1->dir,
                        (secondCH) ? part2->chdir[1] : part2->dir,
                    Patch(part1->patchdir[1], part1->patchsides[2], part1->patchsides[3]),
                    Patch(part2->patchdir[1], part2->patchsides[2], part2->patchsides[3]), r_cm, 1,1);
            }
        }

#ifdef EXTRA_HYDROPHOBIC_ALL_BODY_ATTRACTION
        double extraAttr = 0;
        double extraEpsilon = -1.0 * E_ISO; //  -> isotropic hydrophobic interaction Epsilon

        double extraInteractionSwitch = 0.2; //  -> isotropic hydrophobic interaction Switch

        double distRcm = sqrt(dotrcm);

        if (distRcm > topo.ia_params[part1->type][part2->type].pdis+extraInteractionSwitch) {
            extraAttr = 0.0;
        } else {
            if (distRcm < topo.ia_params[part1->type][part2->type].pdis)
                extraAttr = extraEpsilon;
            else {
                extraAttr = cos(PIH*(distRcm - topo.ia_params[part1->type][part2->type].pdis)/extraInteractionSwitch);
                extraAttr *= extraAttr * extraEpsilon ;
            }
            // cos between two SC
            extraAttr *= (part1->dir.dot(part2->dir));
        }
        atrenergy += extraAttr;
#endif

        return abE + repenergy + atrenergy;
    }
};




class PairE {
    GeoBase* pbc;                   // box size
    EBasic* eFce[MAXT][MAXT];
public:
    PairE(GeoBase* pbc) : pbc(pbc) { initIntFCE(); }

    double operator() (Particle* part1, Particle* part2, ConList* conlist=NULL) {
        Vector r_cm = pbc->image(&part1->pos, &part2->pos);
        double dotrcm = r_cm.dot(r_cm);

        if (dotrcm > topo.sqmaxcut)
            return 0.0;  // distance so far that even spherocylinders cannot be within cutoff

        double dist = sqrt(dotrcm);

        return (*eFce[part1->type][part2->type])(dist, r_cm, part1, part2, conlist); // on fast fce (SPA, SPN) function call slows the sim ~4%
    }

private:

    /**
     * @brief init_intfce Initializes the array with energy functors
     */
    void initIntFCE();
};

#endif // PAIRE_H
