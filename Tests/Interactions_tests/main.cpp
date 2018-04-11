#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <cstdlib>
#include <algorithm>
#include <random>

using namespace std;

class Vector{
public:
    double x,y,z;
    int type;
    Vector(){}
    Vector(double x, double y, double z, int type=0): x(x), y(y), z(z), type(type) {}

    bool operator==(const Vector& o) {
        if(x != o.x) return false;
        if(y != o.y) return false;
        if(z != o.z) return false;
        /*const double aprox = 0.0000000000001;
        if(x < o.x+aprox && x > o.x - aprox)
            return false;
        if(y < o.y+aprox && y > o.y - aprox)
            return false;
        if(z < o.z+aprox && z > o.z - aprox)
            return false;*/

        return true;
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

    bool isAproxSame(const Vector& o) {
        const double aprox = 0.000001;
        if(!(x < o.x+aprox && x > o.x - aprox))
            return false;
        if(!(y < o.y+aprox && y > o.y - aprox))
            return false;
        if(!(z < o.z+aprox && z > o.z - aprox))
            return false;

        return true;
    }

    double size() {
        return sqrt(x*x + y*y + z*z);
    }

    Vector operator*(double a) {
        return Vector(x*a,y*a,z*a, type);
    }

    void operator*=(double a) {
        x*=a;
        y*=a;
        z*=a;
    }

    double dot(Vector& o) {
        return this->x*o.x + this->y*o.y + this->z*o.z;
    }

    Vector operator+(Vector& o) {
        return Vector(x+o.x, y+o.y, z+o.z);
    }

    Vector operator-(Vector& o) {
        return Vector(x-o.x, y-o.y, z-o.z);
    }

    double dist(Vector& o) {
        return sqrt( (this->x-o.x)*(this->x-o.x) + (this->y-o.y)*(this->y-o.y) + (this->z-o.z)*(this->z-o.z) );
    }

    inline Vector cross(const Vector& B) const {
        return Vector(this->y*B.z - this->z*B.y, -this->x*B.z + this->z*B.x, this->x*B.y - this->y*B.x);
    }

    bool isNeighbor(const Vector& o) {
        Vector vec(o.x - x, o.y - y, o.z - z);
        double dist = vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
        return (dist < 1.00000000000001 && dist >  0.999999999999999 );
    }

    inline void rotate(Vector& axis, double angle) {
        double cosAngle = cos(angle);
        double sinAngle = sin(angle);
        double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz;
        double qw = cosAngle, qx = (axis.x * sinAngle), qy = (axis.y * sinAngle), qz = (axis.z * sinAngle);

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
};

std::ostream& operator<<(std::ostream& os, const Vector& vec) {
  os << vec.x << " " << vec.y << " " << vec.z;
  return os;
}

vector< Vector> beads;



void fibonacci_sphere(int samples, int type) {
    const double PI = 3.141592653589793;
    double offset = 2.0/samples;
    double increment = PI * (3.0 - sqrt(5.0));
    double x,y,z,r,phi;

    for(int i=0; i<samples; ++i) {
        z = ((i * offset) - 1) + (offset / 2);
        r = sqrt(1 - z*z);

        phi = i * increment;

        x = cos(phi) * r;
        y = sin(phi) * r;

        beads.push_back(Vector(x,y,z,type));
    }
}



int main(int argc, char* argv[]) {

    if(argc != 6){
        cout << "Parameters: x, y, z, samples, angle_increment" << endl;
        cout << "x, y, z - coordinates" << endl;
        cout << "samples - number of points to generate on sphere (Fibonacci algorithm)" << endl;
        cout << "angle_increment - in radians, sampling full 2pi angle" << endl;
        cout << "\nOutput - SC program coorninates for spherocylinders" << endl;
        cout << "Last line - number of particles written" << endl;
        exit(1);
    }

    //
    // params: x, y ,z ,samples , angle_increment
    //
    int num = 0;
    const int x = atoi(argv[1]);
    const int y = atoi(argv[2]);
    const int z = atoi(argv[3]);
    const int samples = atoi(argv[4]);
    string str2(argv[5]);

    double angle = std::stod(str2);

    fibonacci_sphere(samples, 0);

    for(int i=0; i<samples; ++i) {
        Vector dir = beads[i];

        for(double j=0; j<6.28318530718; j+=angle) {
            Vector par = beads[i];

            par.z = - (dir.x+dir.y) / dir.z;
            par.rotate(dir, j);

            cout << x << " " << y << " " << z << " " << dir.x << " " << dir.y << " " << dir.z << " " << par.x << " " << par.y << " " << par.z << " 0 0" << endl;
            ++num;
        }
    }
    cout << num << endl;

    return 0.0;
}
