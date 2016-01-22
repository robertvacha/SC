#include "analysis.h"

void analyzeCur(double &r1, double &r2, double &fi, int& R1i, int& R1j, int& R2i, int& R2j, int &cen, Conf* conf, Sim* sim) {
    const double deg = 180.0/PI;

    int result = 0;
    int size = 19; // number of spc in sides (hexagon, each line has 19 spc)

    // calculate center of mass
    Vector center(0.0, 0.0, 0.0);
    for(unsigned int i=0; i<conf->pvec.size(); i++) {
        center += conf->pvec[i].pos;
    }
    center.scale(1.0/conf->pvec.size());

    // find particle closest to center of mass
    int partID = 0;
    for(unsigned int i=0; i<conf->pvec.size(); i++) {
        if((center - conf->pvec[i].pos).size() < (center - conf->pvec[partID].pos).size() )
            partID = i;
    }

    cen = partID;
    center = conf->pvec[partID].pos;
    // partID - center of plane -> for curvature calculation

    // find particle positions at edges -> doesnt work, both as distance or energy calc -> calculate between all
    //
    // CONDITION -> 90deg between (center of structure to midpoint of line between two points Vector) AND line between two points Vector

    //
    // Find principle curvatures
    //
    double dist;
    double phi;
    double alfa;
    double R, R_min=999999.9, R_min2 = 999999.9, angle;
    int index;
    Vector temp, c2c;
    Vector ref;
    Vector plane1, plane2, ab, cd;
    Vector midPoint;
    bool c2cOK = false;
    for(unsigned int i = 0; i<conf->pvec.size(); i++) { // search edge particles
        dist = 999.0;
        index = -1;

        //
        // Smallest distance from center of two points to center of mass
        //    AND 90 deg
        //    AND minimal distance between points!!!!
        //    AND center of structure closest point to middle of line between points
        //
        //  Smallest distance so that we are describing the curvature,
        // (imagine a plane, arbitralily choose 3 points. Here the smallest distance will be when the points are in line,
        //  however if there is a curvature, then the 3 points with minimal distance from center of two points to center of mass of structure should descibe curvature
        //
        for(unsigned int j = 0; j<conf->pvec.size(); j++) {
            c2cOK = true;
            temp = (conf->pvec[j].pos - conf->pvec[i].pos); // vector i to j
            temp.scale(0.5);          // half of it
            c2c = ((conf->pvec[i].pos + temp) - center);
            midPoint = conf->pvec[i].pos + temp;


            for(unsigned int k = 0; k<conf->pvec.size(); k++) {
                if(c2c.size() > (midPoint-conf->pvec[k].pos).size()) {
                    c2cOK = false;
                    break;
                }
            }

            if( deg*asin(temp.cross(c2c).size() / (temp.size() * c2c.size())) > 88
                    && temp.size() > 5 // temp == half the distance between points
                    && (c2c.size() < dist)
                    && c2cOK) {
                dist = c2c.size(); // center of two points to mass center distance
                index = j;
            }
        }
        // we have distance and two ID: index and i



        // calculate phi: form dist(d) and temp (x)
        temp = (conf->pvec[index].pos - conf->pvec[i].pos);
        temp.scale(0.5);
        phi = atan((temp).size() / dist);
        alfa = phi*-2 + PI;

        // calculate principal curvature, dist == d (from proof)
        R = temp.size() / sin(alfa);

        // calculate second curvature plain normal vector
        if(R < R_min) { // get R_max
            R_min = R;
            ab = conf->pvec[i].pos - conf->pvec[index].pos;
            R1i = i;
            R1j = index;
            ref = (conf->pvec[i].pos + temp) - center; // center of line to center of structure vector
            plane1 = temp.cross( ref );
        }
    }

    //
    // Repeat process and find second max shifted at least by 45 deg (arbitralily chosen)
    //
    for(unsigned int i = 0; i<conf->pvec.size(); i++) {
        dist = 999.0;
        index = -1;

        //
        // Smallest distance from center of two points on edges to center of mass
        //
        for(unsigned int j = 0; j<conf->pvec.size(); j++) {
            c2cOK = true;
            temp = (conf->pvec[j].pos - conf->pvec[i].pos); // vector i to j
            temp.scale(0.5);          // half of it
            c2c = ((conf->pvec[i].pos + temp) - center);
            midPoint = conf->pvec[i].pos + temp;
            cd = conf->pvec[i].pos - conf->pvec[j].pos;


            for(unsigned int k = 0; k<conf->pvec.size(); k++) {
                if(c2c.size() > (midPoint-conf->pvec[k].pos).size()) {
                    c2cOK = false;
                    break;
                }
            }

            if( deg*asin(temp.cross(c2c).size() / (temp.size() * c2c.size())) > 88
                    && temp.size() > 5 // temp == half the distance between points
                    && (c2c.size() < dist)
                    && c2cOK
                    && deg*asin(ab.cross(cd).size() / (ab.size() * cd.size())) > 80) {
                dist = c2c.size(); // center of two points to mass center distance
                index = j;
            }
        }

        // we have distance and two ID: index and i

        // calculate phi: form dist(d) and temp (x)
        temp = (conf->pvec[index].pos - conf->pvec[i].pos);
        temp.scale(0.5);
        phi = atan((temp).size() / dist);
        alfa = phi*-2 + PI;

        // calculate principal curvature, dist == d (from proof)
        R = temp.size() / sin(alfa);

        // calculate second curvature plain normal vector
        ref = (conf->pvec[i].pos + temp) - center;
        plane2 = temp.cross( ref );
        if(R < R_min2 && (conf->pvec[i].pos - conf->pvec[index].pos).size() > 14.0 && deg*asin(plane1.cross(plane2).size() / (plane1.size() * plane2.size())) > 45 ) { // get R_max
            R_min2 = R;
            R2i = i;
            R2j = index;
            angle = deg*asin(plane1.cross(plane2).size() / (plane1.size() * plane2.size()));
        }
    }


    r1 = R_min;
    r2 = R_min2;
    fi = angle;
    //cout << R_max << " " << R_max2 << " " << angle << endl;

    //
    //  Proof for transformation
    //
    /*const double R=1.0;
    double x;
    double y;
    double d;
    double alfa;
    double phi;
    for(int i=1; i<100; i++) {
        x = (double)i*R/100.0;
        y = sqrt(R*R-x*x);
        d = R-y;
        phi = atan(x/(R-y));
        alfa = asin(x/R);

        cout << x << " " << y << " " << R-y << " " << deg*alfa << " " << deg*phi << " " << deg*2*atan(1/tan(phi)) << " " << deg*-2*phi +180  << endl;
    }*/
    //
    // alfa = 2*atan(1/tan(phi)) = -2*phi + 180
    //

}
