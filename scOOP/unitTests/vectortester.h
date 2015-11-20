#ifndef VECTORTESTER_H
#define VECTORTESTER_H

class VectorTester
{
public:
    VectorTester();

    void test() {
        double R=1.0, h=0.5, r=0.5;
        double S = PI*R*(2*h+r);
        double S_rem = 3*PI - S;
        int in=0;

        ran2.setSeed(5000);
        Vector test[(int)1e6];
        for(int i=0; i<1e6; i++) {
            test[i] = Vector::getRandomUnitConeUniform(Vector(1,0,0), 90);
            if(test[i].x >= 0.5 || test[i].x <= -0.5 )
                in++;
        }
        cout << "in cone: " << in << ", outside: " << 1e6-in << endl;
        cout << "surface cone: " << S << ", rem: " << S_rem << endl;
        exit(0);
    }
};

#endif // VECTORTESTER_H
