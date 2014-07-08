#include "particlestore.h"

/*void ParticleStore::deleteParticleAt(unsigned int num) {
    unsigned int i=0;
    while(num >= store[i].size()) {
        num -= store[i].size();
        i++;
    }

    unsigned int molBegin = num/(i+1);
    molBegin *= (i+1);

    for(int j=i; j >= 0; j--) {
        if((molBegin+j) != num) {
            store[i-1].push_back(store[i][molBegin+j]);
            store[i][molBegin+j] = store[i][store[i].size()-1];
            store[i][store[i].size()-1] = NULL;
            store[i].pop_back();
        } else {
            delete store[i][num];
            store[i][num] = store[i][store[i].size()-1];
            store[i][store[i].size()-1] = NULL;
            store[i].pop_back();
        }
    }

}

void ParticleStore::deleteMoleculeAt(unsigned int num) {
    int i=0;
    while(num >= store[i].size()/(i+1)) {
        num -= store[i].size()/(i+1);
        i++;
    }
    for(int j=0; j<=i; j++) {
        delete store[i][(num)*(i+1) + j];
        store[i][num*(i+1) + j] = store[i][store[i].size()-1];
        store[i][store[i].size()-1] = NULL;
        store[i].pop_back();
    }
}*/
