/** @file wanglandau.h*/

#ifndef WANGLANDAU_H
#define WANGLANDAU_H


#include "stdio.h"
#include "mesh.h"
#include "simlib.h"

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

    /* Wang landau method (wl) */
    double *weights;        /* Array of weights for wl method */
    long   *hist;           /* Array of histogram for wl method */
    long   length[2];       /* Length of above arrays */
    double dorder[2];       /* Increments of order parameter */
    double minorder[2];     /* Minimum order parameter */
    double alpha;           /* Current modifier of weights */
    long   currorder[2];    /* Walue of current order parameter*/
    long   neworder[2];     /* wl order parameter in new step */
    long   max;             /* wl maximum of histogram */
    long   min;             /* wl minimum of histogram */
    double wmin;            /* weights minimum */
    int    wlmdim;          /* Dimwnsionality of wang landau */
    int    wlmtype;         /* Atom type for the Wang landau method (wl) */
    double wl_meshsize;         /* Size of mesh bin for hole order paremeter*/
    Mesh mesh;          /* Mesh for hole order */
    Mesh origmesh;      /* Mesh store for rejected moves */
    long * radiushole;          /* Array for hole radius around origin */
    long * radiusholeold;       /* Array for hole radius around origin-bigmove */
    long   radiusholemax;         /* Size of array for hole radius*/
    long   partincontact;         /* Number of particles in contact */
    long   partincontactold;      /* Number of particles in contact - old for move*/

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
