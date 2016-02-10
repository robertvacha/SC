#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "particle.h"

class GeoBase {
public:
    Vector box;                             ///< \brief Box size */

    GeoBase(){}

    virtual void usePBC(Particle *pos) = 0;
    virtual Vector image(Vector* r1, Vector* r2) = 0;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
    virtual bool boundaryOverlap(Vector *pos) {return false;}       // all - for insertion because of uniformity
#pragma GCC diagnostic pop
    virtual Vector randomPos() =0;
    virtual double volume() = 0;
};

class Cuboid: public GeoBase {
public:

    Cuboid(){}
    Cuboid(Vector box) {
        this->box = box;       
    }

    void info() {
        cout << "Box: " << this->box.info() << endl;
    }

    inline double volume() {
        return box.x*box.y*box.z;
    }

    /**
     * @brief Use of periodic boundary conditions on Particle
     * @param particle
     * Function return particle back in central box
     */ 
    void usePBC(Particle *part) {
        usePBC(part->pos);
    }

    /**
     * @brief use of periodic boundary conditions range 0 - box.particular
     * @param part
     */
    void usePBC2(Particle *part) {
        while (part->pos.x < 0.0) {
            part->pos.x += box.x;
        }
        while (part->pos.x > box.x) {
            part->pos.x -= box.x;
        }
        while (part->pos.y < 0.0) {
            part->pos.y += box.y;
        }
        while (part->pos.y > box.y) {
            part->pos.y -= box.y;
        }
        while (part->pos.z < 0.0) {
            part->pos.z += box.z;
        }
        while (part->pos.z > box.z) {
            part->pos.z -= box.z;
        }
    }


     /**
      * @brief usePBC
     * Function take point in reduce coordinate space and then put point (vector) in central box
     * @param vec - reference to vector in real coordinates
     */
    void usePBC( Vector &vec ) {
        while ( vec.x < 0.0 ) {
            vec.x += 1.0;
        }
        while ( vec.x > 1.0 ) {
            vec.x -= 1.0;
        }
        while ( vec.y < 0.0 ) {
            vec.y += 1.0;
        }
        while ( vec.y > 1.0 ) {
            vec.y -= 1.0;
        }
        while ( vec.z < 0.0 ) {
            vec.z += 1.0;
        }
        while ( vec.z > 1.0 ) {
            vec.z -= 1.0;
        }
    }


    /**
     * @brief Returns the vector pointing from the centre of mass of particle 2 to the
       centre of mass of the closest image of particle 1.
       BECAREFULL r_cm vector is in real values ... in other words when you use image finction
       in code u have to normalize it by box size!!
     * @param r1
     * @param r2
     * @param box
     * @return
     */
    inline Vector image(Vector* r1, Vector* r2) {
        Vector r_cm( r1->x - r2->x, r1->y - r2->y, r1->z - r2->z );

        Vector temp(r_cm.x + 6755399441055744.0, r_cm.y + 6755399441055744.0, r_cm.z + 6755399441055744.0);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
        r_cm.x = box.x * (r_cm.x - static_cast<double>(reinterpret_cast<int&>(temp.x) ) );
        r_cm.y = box.y * (r_cm.y - static_cast<double>(reinterpret_cast<int&>(temp.y) ) );
        r_cm.z = box.z * (r_cm.z - static_cast<double>(reinterpret_cast<int&>(temp.z) ) );
#pragma GCC diagnostic pop
        return r_cm;
    }

    virtual Vector randomPos() {
        return Vector(ran2(), ran2(), ran2());
    }

};

/**
 * @brief The Wedge class
 * Solid rotational boundaries - inner and outer radius - around Z axis
 * Unit Z pbc
 * rotational XY pbc
 *
 * CHAINS NOT SUPPORTED
 *
 */
class Wedge : public GeoBase {
private:
    double angleSin;
    double angleCos;
    double angleRad;

    inline void rotateClockWise(Vector* pos) { // OK
        double x = pos->x;
        pos->x = x * angleCos + pos->y * angleSin;
        pos->y = pos->y * angleCos - x * angleSin;
    }

    inline void rotateCounterClock(Vector* pos) { // OK
        double x = pos->x;
        pos->x = x * angleCos - pos->y * angleSin;
        pos->y = pos->y * angleCos + x * angleSin;
    }

public:
    double innerR; // diameter = 2 * radius
    double outerR;
    double angleDeg;

    Wedge(){cout << "NO PARAMETERS GIVEN !!!" << endl;}
    Wedge(double z, double angle, double outerR, double innerR) : angleDeg(angle) {
        this->box = Vector(outerR, outerR, z);
        cout << "Wedge: box: " << this->box.info() << " Angle:" << angleDeg
             << " Inner radius: " << innerR << " Outer radius: " << outerR << endl;
        if(innerR > outerR) {
            cout << "Mixed Inner and outer radius!!!" << endl;
            exit(1);
        }
        if(0.0 <= angle && angle > 90) {
            cout << "Angle: " << angle << " must be between 0 and 90 degrees!!!" << endl;
            exit(1);
        }
        if(2*innerR*innerR < topo.sqmaxcut) {
            cout << "inner radius too small for interactions of species, see constructor Wedge" << endl;
            exit(1);
        }
        this->innerR = innerR/box.x;
        this->outerR = outerR/box.x;
        angleRad = angle*PI / 180.0;
        angleSin = sin(angleRad);
        angleCos = cos(angleRad);
    }

    inline double volume() {
        return box.x*box.x*box.z*angleRad*0.5*(outerR*outerR - innerR*innerR);
    }

    inline bool boundaryOverlap(Vector *pos) { // OK
        double radius = pos->x*pos->x + pos->y*pos->y;
        if(radius > outerR*outerR || radius < innerR*innerR)
            return true;
        // check plane YZ (angle 0) -> x go negative -> rotate clockwise
        if(pos->x < 0.0)
            return true;
        // check angle
        if(atan2(pos->x, pos->y) > angleRad)
            return true;
        return false;
    }

    virtual void usePBC(Particle *part) {
        // check plane YZ (angle 0) -> x go negative -> rotate clockwise
        if(part->pos.x < 0.0) {
            rotateClockWise(&part->pos);
            //part->pscRotate(angleRad*0.5, topo.ia_params[part->type][part->type].geotype[0], Vector(0,0,1), false);
        }
        // check angle
        if(atan2(part->pos.x, part->pos.y) > angleRad) {
            rotateCounterClock(&part->pos);
            //part->pscRotate(angleRad*0.5, topo.ia_params[part->type][part->type].geotype[0], Vector(0,0,1), true);
        }

        /*do { // unit lenght scaling for z axis
            (part->pos).z += box.z;
        } while ((part->pos).z < 0.0);
        do {
            (part->pos).z -= box.z;
        } while ((part->pos).z > box.z);*/
    }

    virtual Vector image(Vector* r1, Vector* r2) {
        Vector r_cm, r_cm2 = *r1;
        r_cm.x = r1->x - r2->x;
        r_cm.y = r1->y - r2->y;
        r_cm.z = r1->z - r2->z;

        double temp =  r_cm.z + 6755399441055744.0;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
        r_cm.z = box.z * (r_cm.z - static_cast<double>(reinterpret_cast<int&>(temp) ) );
#pragma GCC diagnostic pop
        if(atan2(r_cm2.x, r_cm2.y) < angleRad*0.5)
            rotateClockWise(&r_cm2);
        else rotateCounterClock(&r_cm2);

        r_cm2.x = r_cm2.x - r2->x;
        r_cm2.y = r_cm2.y - r2->y;
        r_cm2.z = r_cm2.z - r2->z;

        temp =  r_cm2.z + 6755399441055744.0;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
        r_cm2.z = box.z * (r_cm2.z - static_cast<double>(reinterpret_cast<int&>(temp) ) );
#pragma GCC diagnostic pop
        if(r_cm.dot(r_cm) < r_cm2.dot(r_cm2) )
            return r_cm;
        else return r_cm2;
    }

    virtual Vector randomPos() {
        Vector pos;
        pos.randomUnitCube();
        while(boundaryOverlap(&pos)) {
            pos.randomUnitCube();
        }
        return pos;
    }

};

#endif // GEOMETRY_H
