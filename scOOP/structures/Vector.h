/** @file Vector.h*/

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include "quaternion.h"
#include <sstream>
#include "../mc/randomGenerator.h"

class Quat;

class Vector {
public:
    double x, y, z;
    Vector(){}
    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    std::string info() {
        std::ostringstream o;
        o << "(" <<x << ", " << y << ", " << z <<")";
        return o.str();
    }

    /**
     * @brief Return a norm of Vector
     * @return
     */
    inline double size() {
        return sqrt( pow(this->x,2) + pow(this->y,2) + pow(this->z,2));
    }

    /**
     * @brief normalise Normalise a vector to have unit length.  For speed during heavy use, it is
       not checked that the supplied vector has non-zero length.
     */
    inline void normalise() {
        double tot = size();
        if (tot !=0.0) {
            tot = 1.0 / tot;
            x *= tot;
            y *= tot;
            z *= tot;
        }
    }

    /**
     * @brief dotProduct == Scallar product
     * @param other
     * @return
     */
    inline double dot(Vector& other)  {
        return x*other.x + y*other.y + z*other.z;
    }

    inline void scale(double scale) {
        x=x*scale; y=y*scale, z=z*scale;
    }

    inline Vector operator- (const Vector& o) const {
        return Vector(x-o.x, y-o.y,z-o.z);
    }

    inline void operator-= (const Vector& o) {
        x-=o.x;
        y-=o.y;
        z-=o.z;
    }

    inline void operator+= (const Vector& o) {
        x+=o.x;
        y+=o.y;
        z+=o.z;
    }

    inline bool operator==(Vector& o) {
        if(x == o.x && y == o.y && z == o.z) return true;
        else return false;
    }

    inline bool operator!= (Vector& o) {
        if(o.x != x || o.y != y || o.z != z ) return true;
        return false;
    }

    inline Vector operator* (double scale) {
        return Vector(this->x*scale, this->y*scale, this->z*scale);
    }

    inline void operator*= (double scale) {
        this->x*=scale, this->y*=scale, this->z*=scale;
    }

    friend inline Vector operator* (double,Vector&);

    inline Vector cross(Vector& B) {
        return Vector(this->y*B.z - this->z*B.y, -this->x*B.z + this->z*B.x, this->x*B.y - this->y*B.x);
    }

    inline void ortogonalise(Vector& B) {
        double dp(dot(B));    this->x -= dp*B.x;    this->y -= dp*B.y;    this->z -= dp*B.z;
    }

    // using Vector&, double, double instead of Quat -> circular include
    inline void rotate(Vector& vec, double vc, double vs) {
        double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz, qw(vc), qx(vec.x * vs), qy(vec.y * vs), qz(vec.z * vs);

        /*    t1 = quat.w * quat.w; */
        t2 =  qw * qx;
        t3 =  qw * qy;
        t4 =  qw * qz;
        t5 = -qx * qx;
        t6 =  qx * qy;
        t7 =  qx * qz;
        t8 = -qy * qy;
        t9 =  qy * qz;
        t10 = -qz * qz;

        newx = 2.0 * ( (t8+t10) * x + (t6-t4)  * y + (t3+t7) * z ) + x;
        newy = 2.0 * ( (t4+t6)  * x + (t5+t10) * y + (t9-t2) * z ) + y;
        newz = 2.0 * ( (t7-t3)  * x + (t2+t9)  * y + (t5+t8) * z ) + z;

        x = newx;
        y = newy;
        z = newz;
    }

    inline void rotate(Quat& q) {
        double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz;

        /*    t1 = quat.w * quat.w; */
        t2 =  q.w * q.x;
        t3 =  q.w * q.y;
        t4 =  q.w * q.z;
        t5 = -q.x * q.x;
        t6 =  q.x * q.y;
        t7 =  q.x * q.z;
        t8 = -q.y * q.y;
        t9 =  q.y * q.z;
        t10 = -q.z * q.z;

        newx = 2.0 * ( (t8+t10) * x + (t6-t4)  * y + (t3+t7) * z ) + x;
        newy = 2.0 * ( (t4+t6)  * x + (t5+t10) * y + (t9-t2) * z ) + y;
        newz = 2.0 * ( (t7-t3)  * x + (t2+t9)  * y + (t5+t8) * z ) + z;

        x = newx;
        y = newy;
        z = newz;
    }

    /**
     * @brief ranvec    Returns an evenly distributed random unit vector2 of unit length.
                        See Allen & Tildesley p349 or Frenkel & Smit p410.
     * @return    RANDOM vector2 ON UNIT SPHERE
     */
    inline void randomUnitSphere() {
        double a, b, xi1, xi2;

        do {
            xi1 = 1.0 - 2.0*ran2();
            xi2 = 1.0 - 2.0*ran2();

            a = xi1*xi1 + xi2*xi2;
        } while (a > 1.0);

        b = 2.0 * sqrt(1.0 - a);

        x = xi1 * b;
        y = xi2 * b;
        z = 1.0 - 2.0*a;
    }

    inline void randomUnitCube() {
        x = ran2();
        y = ran2();
        z = ran2();
    }

    static inline Vector getRandomUnitCube() {
        return Vector(ran2(), ran2(), ran2());
    }

    static inline Vector getRandomUnitSphere() {
        Vector vec;
        vec.randomUnitSphere();
        return vec;
    }
};

inline Vector operator* (double scale, Vector& vec) {
    return Vector(vec.x*scale, vec.y*scale, vec.z*scale);
}




#endif // VECTOR_H
