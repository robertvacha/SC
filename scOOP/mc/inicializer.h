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

    MolIO* molecules;    ///< @brief List of AtomType parameters read from init.top

    char *sysnames[MAXN];   ///< @brief List of MoleculeType names of system
    char *poolNames[MAXN];  ///< @brief List of MoleculeType names of pool

    long  *sysmoln /*[MAXN]*/;
    long  *poolMolNum /*[MAXN]*/;


public:
    Inicializer(Sim* sim, Conf *conf, FileNames* files):
        poolConfig(false), sim(sim), conf(conf), files(files) {

        for(int i=0; i<MAXN; i++) {
            sysnames[i] = NULL;
            poolNames[i] = NULL;
        }

        molecules = new MolIO[MAXMT];

        sysmoln = (long int*) malloc( sizeof(long)*MAXN);
        if(sysmoln == NULL){
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sysmoln");
            exit(1);
        }

        poolMolNum = (long int*) malloc( sizeof(long)*MAXN);
        if(poolMolNum == NULL){
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for poolMolNum");
            exit(1);
        }
    }

    /*
     *  INPUT
     */

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

       @return True - successful, False - unsuccessful
     */
    bool initConfig(FILE** infile, std::vector<Particle > &pvec);

    /// @brief test if simulation contains Chains, sets probability of chain move to 0 if no chains
    void testChains();

    /// @brief Initializes the pairlist and allocates memory
    void initNeighborList();

private:

    void topDealoc() {
        delete[] molecules;

        if (sysmoln != NULL) free(sysmoln);
            sysmoln=NULL;

        if (poolMolNum != NULL) free(poolMolNum);
            poolMolNum=NULL;

        for (int i=0;i<MAXN;i++) {
            if ((sysnames[i]) != NULL) free(sysnames[i]);
                sysnames[i]=NULL;
        }

        for (int i=0;i<MAXN;i++) {
            if ((poolNames[i]) != NULL) free(poolNames[i]);
                poolNames[i]=NULL;
        }
    }

    void initGroupLists();

    void initSwitchList();

    void setParticlesParams() {
        setParticlesParamss(molecules, sysmoln, sysnames, &conf->pvec);
        setParticlesParamss(molecules, poolMolNum, poolNames, &conf->pool);
    }

    void readTopoFile(bool exclusions[][MAXT]);

    void setParticlesParamss(MolIO *molecules, long  *sysmoln, char **sysnames, std::vector<Particle >* pvec);

    /**
     * @brief xmalloc nice malloc, which does the error checking for us
     * @param num
     * @return
     */
    void* xMalloc (size_t num);

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
    int fillSystem(char *pline, char *sysnames[MAXN], long **sysmoln, char* name);

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
    int fillMol(char *molname, char *pline, MolIO *molecules);
};

#endif // INICIALIZER_H
