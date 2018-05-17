/** @file inicializer.h*/

#ifndef INICIALIZER_H
#define INICIALIZER_H

/*
 *  Initialize Topology, Inicial configuration,
 */

#include "../structures/sim.h"

using namespace std;

class Inicializer
{  
public:
    bool poolConfig;
private:
    Sim* sim;                  // Should contain the simulation options.
    Conf* conf;                // Should contain fast changing particle and box(?) information
    FileNames* files;

public:
    Inicializer(Sim* sim, Conf *conf, FileNames* files):
        poolConfig(false), sim(sim), conf(conf), files(files) {}

    /*
     *  INPUT
     */

    /**
     * @brief Config initialization

       Reads in the initial configuration from the file "config.init".  Each line
       contains the three components of the position vector and three components of
       the direction vector and three components of patch direction for a spherocylinder.
       The direction vector is normalised
       after being read in.  The configuration is checked for particle overlaps.

       @return True - successful, False - unsuccessful
     */
    bool initConfig(FILE** infile, std::vector<Particle > &pvec, bool scale_to_box = true);

    /// @brief test if simulation contains Chains, sets probability of chain move to 0 if no chains
    void testChains();

    /// @brief Initializes the pairlist and allocates memory
    void initNeighborList();

public:

    void initGroupLists();

    void setParticlesParams() {
        setParticlesParamss(&conf->pvec, topo.system);
        setParticlesParamss(&conf->pool, topo.pool);
    }

    void readTopoFile();

    /**
     * @brief setParticlesParamss - Generate System and Pool according to what is stored in System class
     * @param molecules
     * @param pvec
     * @param system
     */
    void setParticlesParamss(std::vector<Particle >* pvec, System& system);

    /**
     * @brief xmalloc nice malloc, which does the error checking for us
     * @param num
     * @return
     */
    void* xMalloc (size_t num);

    /**
     * @brief Converts the geometrical type string into a number
     * @param geotype
     * @return
     */
    int convertGeotype(char * geotype);
};

#endif // INICIALIZER_H
