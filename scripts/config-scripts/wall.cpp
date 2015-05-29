#include <iostream>
#include <cstdlib>
#include <fstream>
#include <math.h>

using namespace std;

int main() {
    ofstream out("config.init");


    const int edge = 14;
    double equiDist=0.5;
    const double boxX = 20.0;
    const double boxY = 20.0;
    const double boxZ = 20.0;

    double height = sqrt( (equiDist*equiDist)-(equiDist*0.25*equiDist));
    const double x = boxX/2;
    double y = ( boxY-(equiDist*(edge+edge-1)) ) *0.75;
    double z = boxZ*0.5 + (edge-1)*height;

    int rowNum = edge;
    int count =0;

    if(out.is_open() ) {

        out <<  " " << boxX << "  " << boxY << "  " << boxZ << "\n";

        for(int i=0; i<(edge+edge-1); i++) { // vyska
            for(int j=0; j< rowNum; j++) { // sirka
                if(x>=boxX-1.0 || y>=boxY-1.0 ||z>=boxZ-1.0 || x<= 1.0 || y<= 1.0 || z<= 1.0) {
                    cout << "Too big\n";
                    exit(1);
                }
                out << " " << x << " " << y << " " << z ;
                out << "   1.0  0.0  0.0   0.0 1.0 0.0 0\n";

                y+=0.5;
                count++;

                cout << " .";
            }
            cout << endl;
            y-= 0.5*rowNum;

            if(i<(edge-1)) {
                    rowNum++;
                    y -= 0.25;
            } else {
                rowNum--;
                y += 0.25;
            }

            z -= height;

        }

        out.close();
    }

    cout << "count: " << count << endl;

    return 0;
}

