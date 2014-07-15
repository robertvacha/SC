/** @file inicializer.h*/

#ifndef INICIALIZER_H
#define INICIALIZER_H

/*
 *  Initialize Topology, Inicial configuration,
 */

#include "mygetline.h"
#include <iostream>

#include "simlib.h"
#include "../structures/sim.h"

using namespace std;

class Inicializer
{
public:
    Inicializer(Topo *topo, Sim* sim, Conf *conf, FileNames* files):
        topo(topo), sim(sim), conf(conf), files(files) {}

private:
    long seed;
    Topo* topo;                // will maybe contain all the topo stuff in future
    Sim* sim;                  // Should contain the simulation options.
    Conf* conf;                // Should contain fast changing particle and box(?) information
    FileNames* files;

public:

    long int getSeed() {return seed;}

    /*
     *  INPUT
     */

    /**
     * @brief Reads the run parameters from the external file "options".  See the end of the
       code for a template. All comments starting with '#' are stripped out.  The
       options are summarised on standard output and checked for validity of range.
     */
    void readOptions();




    /**
     * @brief Inicialization of topology

       Create lists for chain operations: Connectivity list where it is written for each sc
       with which sc it is connected. The order is important because spherocylinders have direction
       First is interacting tail then head. Chain list where particles are assigned to chains to
       which they belong

     */
    void initTop();

    /**
     * @brief Config initialization

       Reads in the initial configuration from the file "config.init".  Each line
       contains the three components of the position vector and three components of
       the direction vector and three components of patch direction for a spherocylinder.
       The direction vector is normalised
       after being read in.  The configuration is checked for particle overlaps.

     */
    void initConfig();

    /**
     * @brief test if simulation contains Chains, sets probability of chain move to 0 if no chains
     */
    void testChains();

    /**
     * @brief Sets names of "write files"
     */
    void initWriteFiles();

    /**
     * @brief Initializes the pairlist and allocates memory
     */
    void initPairlist();

    /**
     * @brief Paralel tempering(Replica exchange move) initialization
     */
    void initMPI();

private:

    void readTopoFile(Molecule *molecules, long *sysmoln, char *sysnames[], bool exclusions[][MAXT]);

    void openTopoFile(FILE* infile);

    void allocSysmoln(long* sysmoln);

    /**
     * @brief xmalloc nice malloc, which does the error checking for us
     * @param num
     * @return
     */
    void* xMalloc (size_t num);

    /**
     * @brief dealocating memory for initTop
     * @param pline
     * @param sysnames
     * @param sysmoln
     * @param molecules
     * @return
     */
    int topDealoc(char *sysnames[MAXN], long **sysmoln);

    /**
     * @brief filling pair for which we exlude attraction interaction. Returns 1 on succes.
     * @param pline
     * @param (*exlusions)[][]
     * @return
     */
    int fillExclusions(char **pline, bool exlusions[][MAXT]);

    /**
     * @brief filling the system parameters
     * @param pline
     * @param sysnames
     * @param sysmoln
     * @return
     */
    int fillSystem(char *pline, char *sysnames[MAXN], long **sysmoln);

    /**
     * @brief filing the parameters for types from given strings. Returns 1 on succes.
     * @param pline
     * @return
     */
    int fillTypes(char **pline);

    /**
     * @brief Converts the geometrical type string into a number
     * @param geotype
     * @return
     */
    int convertGeotype(char * geotype);

    /**
     * @brief filling the parameters of external potentail - wall. Returns 1 on succes.
     * @param pline
     * @return
     */
    int fillExter(char **pline);

    /**
     * @brief filling the parameters for molecules
     * @param molname
     * @param pline
     * @param molecules
     * @return
     */
    int fillMol(char *molname, char *pline, Molecule * molecules);

    /**
     * @brief generate interations pairs
     * @param (*exlusions)[][]
     */
    void genParamPairs(bool (*exclusions)[MAXT][MAXT]);

    /**
     * @brief use of periodic boundary conditions
     * @param pos
     * @param pbc
     */
    void usePBC(Vector *pos,Vector pbc);

    /**
     * @brief convert string num into two integers
     * @param num
     * @param value
     */
    void readii2(char * num, int value[2]);

    /**
     * @brief convert string num into double
     * @param num
     * @return
     */
    double readd2(char * num);

    /**
     * @brief convert string num into long
     * @param num
     * @return
     */
    long readl2(char * num);

    /**
     * @brief convert string num into integer
     * @param num
     * @return
     */
    int readi2(char * num);
};

#endif // INICIALIZER_H
