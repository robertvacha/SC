/** @file wanglandau.h*/

#ifndef WANGLANDAU_H
#define WANGLANDAU_H

#include "stdio.h"
#include "mesh.h"
#include "simlib.h"

/**
 * @brief The WangLandau class - Wang landau method (wl)
 */
class WangLandau
{
public:
    WangLandau() {}

    void setTopoConf(Topo* topo, Conf* conf) {
        this->topo = topo;
        this->conf = conf;
    }

    Topo* topo;
    Conf* conf;

    double *weights;           ///< \brief Array of weights for wl method
    long   *hist;              ///< \brief Array of histogram for wl method
    long   length[2];          ///< \brief Length of above arrays
    double dorder[2];          ///< \brief Increments of order parameter
    double minorder[2];        ///< \brief Minimum order parameter
    double alpha;              ///< \brief Current modifier of weights
    long   currorder[2];       ///< \brief Walue of current order parameter
    long   neworder[2];        ///< \brief wl order parameter in new step
    long   max;                ///< \brief wl maximum of histogram
    long   min;                ///< \brief wl minimum of histogram
    double wmin;               ///< \brief weights minimum
    int    wlmdim;             ///< \brief Dimwnsionality of wang landau
    int    wlmtype;            ///< \brief Atom type for the Wang landau method (wl)
    double wl_meshsize;        ///< \brief Size of mesh bin for hole order paremeter
    Mesh mesh;                 ///< \brief Mesh for hole order
    Mesh origmesh;             ///< \brief Mesh store for rejected moves
    long* radiushole;          ///< \brief Array for hole radius around origin
    long* radiusholeold;       ///< \brief Array for hole radius around origin-bigmove
    long  radiusholemax;       ///< \brief Size of array for hole radius
    long  partincontact;       ///< \brief Number of particles in contact
    long  partincontactold;    ///< \brief Number of particles in contact - old for move

public:

    /**
     * @brief z_order
     * @param wli
     * @return
     */
    long zOrder(int wli);

    /**
     * @brief twopartdist
     * @param wli
     * @return 2 particles distance
     */
    long twoPartDist(int wli);

    /**
     * @brief wlwrite
     * @param filename
     * @return
     */
    int write(char filename[30]);


    /**
     * @brief wlreject
     * @param oldlength
     */
    void reject(long oldlength, int wlm[2]);

    /**
     * @brief wlaccept
     * @param wlm
     */
    void accept(int wlm);

    /**
     * @brief wlend
     * @return
     */
    int end();

    /**
     * @brief wlinit Initiate Wang-Landau calculation.
     * @param filename
     * @return
     */
    int initCalc(char filename[30]);

    inline char* wlfgets(char *line, int n, FILE *stream) {
        char *c;
        if (fgets(line,n,stream)==NULL)     return NULL;
        if ((c=strchr(line,'\n'))!=NULL)    *c=0;
        return line;
    }

};

#endif // WANGLANDAU_H
