/** @file quaternion.h*/

#ifndef QUATERNION_H
#define QUATERNION_H

#include <iostream>

using namespace std;

class Quat{
public:
    double w,x,y,z;

    Quat(double w,double x, double y, double z): w(w), x(x), y(y), z(z) {}

    void print() {
        cout <<"( "<<w<<", "<<x<<", "<<y<<", "<<z<<" )"<< endl;
    }
};

#endif // QUATERNION_H
