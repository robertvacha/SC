#ifndef TOPO_H
#define TOPO_H

#include "macros.h"
#include "structures.h"


/* It would be nice, if this struct would contain all the topo stuff in the end*/
class Topo
{
public:
    double sqmaxcut;    ///< \brief square of distance over which even spherocylinders cannot interact (distance between CM)
    double maxcut;      ///< \brief distance over which even spherocylinders cannot interact (distance between CM)

    MoleculeParams moleculeParam[MAXMT];   ///< \brief parameters for Molecules

    Ia_param ia_params[MAXT][MAXT];     ///< \brief parametrization of particles for all interations
    Exters exter;                       ///< \brief external potential - wall

    //
    //  METHODS
    //

    Topo():sqmaxcut(0) {}
    ~Topo() {
        printf ("Deallocating Topo...\n");

        for(int i=0; i<MAXMT; i++) {
            free(moleculeParam[i].name);
        }
    }

    /**
     * @brief generate interations pairs
     * @param (*exlusions)[][]
     */
    void genParamPairs(bool (*exclusions)[MAXT][MAXT]);

    void genTopoParams();
};

#endif // TOPO_H
