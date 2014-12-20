/** @file wanglandau.h*/

#ifndef WANGLANDAU_H
#define WANGLANDAU_H

#include "mesh.h"
#include "simlib.h"
#include <cstring>

extern long int test[30];
extern bool cond;

/**
 * @brief The WangLandau class - Wang landau method (wl)
 */
class WangLandau
{
public:
    WangLandau() {}

    void setConf(Conf* conf) {
        this->conf = conf;
    }

    Conf* conf;

    int wlm[2];                ///< \brief Wang landau method (wl)

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
     * @brief init Initialization of wang-landaou method
     * @param wlinfile
     */
    void init(char wlinfile[30]);

    /**
     * @brief endWangLandau Finalizing of wang-landaou method
     * @param wloutfile
     */
    void endWangLandau(char wloutfile[30]);

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

    /**
     * @brief radiushole_all filling the radiushole above vec
     * @param wli
     * @param position
     * @return
     */
    long radiusholeAll(int wli, Vector *position);

    /**
     * @brief contParticle_all filling all Particle in the contact
     * @param wli
     * @return
     */
    long contParticlesAll(int wli);

    /**
     * @brief radiusholeorder_moveone return change in order parameter when one particle moves
     * @param oldpos
     * @param target
     * @param wli
     * @param position
     * @return
     */
    long radiusholeOrderMoveOne(Vector *oldpos, long target, int wli, Vector *position);

    /**
     * @brief contParticle_moveone return change in number of Particle in contact when one particle moves
     * @param oldpos
     * @param target
     * @param wli
     * @return
     */
    long contParticlesMoveOne(Vector *oldpos, long target,int wli);

    /**
     * @brief radiusholeorder_movechain return change in order parameter when chain moves
     * @param chain
     * @param chorig
     * @param wli
     * @param position
     * @return
     */
    long radiusholeOrderMoveChain(vector<int> chain, Particle chorig[MAXCHL], int wli, Vector *position);

    /**
     * @brief contParticle_movechain return change in order parameter when chain moves
     * @param chain
     * @param conf
     * @param sim
     * @param chorig
     * @param wli
     * @return
     */
    long contParticlesMoveChain(vector<int> chain, Particle chorig[MAXCHL], int wli);

private:
    /**
     * @brief radiushole_position return order of given radius
     * @param radius
     * @param wli
     * @return
     */
    long radiusholePosition(double radius, int wli) {
        return (long) ceil( ( radius - minorder[wli]) / dorder[wli]  );
    }

    /**
     * @brief radiushole_order return current bin of free radius
     * @return
     */
    long radiusholeOrder();

    /**
     * @brief radiushole_print
     * @param radiushole
     * @param length
     */
    void radiusholePrint(long *radiushole, long length);

    /**
     * @brief contParticle_order return order for Particle in contact
     * @param wli
     * @return
     */
    inline long contParticlesOrder(int wli) {
        return (long) ceil( ( partincontact - minorder[wli]) / dorder[wli]  );
    }

    /**
     * @brief particleinncontact returns if particle is in contact
     * @param vec
     * @return
     */
    bool particlesInContact(Vector *vec);

};

#endif // WANGLANDAU_H
