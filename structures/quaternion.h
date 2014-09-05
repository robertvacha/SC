/** @file quaternion.h*/

#ifndef QUATERNION_H
#define QUATERNION_H

class Quat{
public:
    double w,x,y,z;

    Quat(double w,double x, double y, double z): w(w), x(x), y(y), z(z) {}
};

#endif // QUATERNION_H
