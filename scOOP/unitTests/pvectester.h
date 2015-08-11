#ifndef PVECTESTER_H
#define PVECTESTER_H

#include "../structures/Conf.h"
#include <vector>

using namespace std;

class pVecTester
{
public:
    pVecTester() {}

    bool test() {
        ParticleVector pvec;
        vector<Particle> part;

        part.push_back(Particle(Vector(1.0, 1.0, 1.0), Vector::getRandomUnitSphere(), Vector::getRandomUnitSphere(), 1, 10));

        pvec.molTypeCount = 1;
        pvec.first[0] = 0;
        pvec.first[1] = 0;
        pvec.first[2] = 0;
        pvec.calcChainCount();

        for(int i=0; i<10; i++) {
            pvec.insertMolecule(part);
            part[0].pos.x += 1.0;
        }

        // pvec of 10 particles initialized
        for(int i=0; i<10; i++) {
            if(pvec.getMolecule(i,1,1)[0] != i) {
                cout << "getMolecule FAILED" << endl;
                return false;
            }
        }

        return true;

    }
};

#endif // PVECTESTER_H
