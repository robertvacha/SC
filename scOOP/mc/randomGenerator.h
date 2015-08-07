/** @file randomGenerator.h*/

#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include "dSFMT-src-2.2.3/dSFMT.h"
#include <iostream>

using namespace std;

class randomBase {
public:
    randomBase() {}

    virtual void setSeed(long int seed) =0;

    /**
     * @return random number in close-open range [0,1)
     */
    virtual double operator()() =0;
};

class Ran2 : public randomBase
{
public:
    Ran2() : seed(13) {
        cout << "Using RAN2 random generator" << endl;
    }

    void setSeed(long int seed) {this->seed=seed;}
    int getSeed() {return seed;}

    /**
     * @brief ran2 from Numerical Recipes.
     * @param idum
     * @return
     */
    double operator()();

private:
    long int seed;
};

class Dsfmt : public randomBase
{
public:
    Dsfmt() {
        cout << "Using double precision Mersenne twister" << endl;
        setSeed(13);
    }

    void setSeed(long int num) {
        dsfmt_init_gen_rand(&seed, num);
    }

    double operator()() {
        return dsfmt_genrand_close1_open2(&seed) - 1.0;
    }

private:
    dsfmt_t seed;
};


#endif // RANDOMGENERATOR_H
