/** @file ran2.h*/

#ifndef RAN2_H
#define RAN2_H

class Ran2
{
public:
    Ran2(long int seed): seed(seed){}

    long int seed;

    /**
     * @brief ran2 ran2 from Numerical Recipes.
     * @param idum
     * @return
     */
    double operator()();
};

#endif // RAN2_H
