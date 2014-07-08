/** @file particlestore.h*/

#ifndef PARTICLESTORE_H
#define PARTICLESTORE_H

#include <vector>
#include <cstdio>

#include "particle.h"
#include "macros.h"

class ParticleStore // tested ok
{
public:
    std::vector<Particle* > store[MAXMT];

    /** @brief How many particles per molecule for each type */
    int partInMol[MAXMT];

    ParticleStore() {
        for(int i=0; i<MAXMT; i++) partInMol[i] = 0;
    }

    ~ParticleStore() {
        std::vector<Particle*>::iterator iter;
        int i=0;
        while(partInMol[i] != 0) {
            for(iter = store[i].begin(); iter != store[i].end(); iter++){
                delete (*iter)->neighborID;
                delete *iter;
            }
            i++;
        }
    }

    /*int particleCount() {
        int count=0;
        for(unsigned int i=0; i < store.size(); i++) {
            count += store[i].size();
        }
        return count;
    }

    int moleculeCount() {
        int count=0;
        for(unsigned int i=0; i < store.size(); i++) {
            count += store[i].size()/(i+1);
        }
        return count;
    }

    Particle* getParticleAt(unsigned int num) {
        int i=0;
        while(num >= store[i].size()) {
            num -= store[i].size();
            i++;
        }
        return store[i][num];
    }

    Particle* getMoleculeAt(unsigned int num) {
        int i=0;
        while(num >= store[i].size()/(i+1)) {
            num -= store[i].size()/(i+1);
            i++;
        }
        return store[i][(num)*(i+1)];
    }

    void pushParticle(Particle* part,int chainsize=1) {
        store[chainsize-1].push_back(part);
    }

    void deleteParticleAt(unsigned int num);

    void deleteMoleculeAt(unsigned int num);*/



};

#endif // PARTICLESTORE_H
