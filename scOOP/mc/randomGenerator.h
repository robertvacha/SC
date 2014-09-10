/** @file randomGenerator.h*/

#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

class randomBase {
public:
    randomBase() {}

    virtual void setSeed(long int seed) =0;
    virtual int getSeed() =0;

    /**
     * @return random number in close-open range [0,1)
     */
    virtual double operator()() =0;
};

class Ran2 : public randomBase
{
public:
    Ran2() : seed(13) {}

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

/*#if defined(MERSENNETWISTER)
  typedef Faunus::RandomTwister<double,std::mt19937> slump;
#else*/
  typedef Ran2 ranGen;
//#endif
  extern ranGen ran2;


#endif // RANDOMGENERATOR_H
