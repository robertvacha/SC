/** @file wanglandau.h*/

#ifndef WANGLANDAU_H
#define WANGLANDAU_H

#include "../structures/sim.h"
#include "mesh.h"
#include "simlib.h"
#include <cstring>

/**
 * @brief The WangLandau class - Wang landau method (wl)
 * Option 1: position in box respective to center of mass
 * Option 3: cos(angle of orientation to z axis)
 */
class WangLandau
{
#ifdef ENABLE_MPI
    MPI_Win _winInt;
    MPI_Win _winDouble;
    MPI_Win _winShared;
    MPI_Aint _sizeInt;
    MPI_Aint _sizeDouble;
    MPI_Aint _sizeShared;
#endif
public:
    const int mpirank;

    WangLandau(Conf* conf, Sim* sim) : conf(conf), mpirank(sim->mpirank), mesh(true), origmesh(false) {
        wlmtype = sim->wlmtype;
        wlm[0] = sim->wlm[0];
        wlm[1] = sim->wlm[1];
    }

    Conf* conf;

    int wlm[2];                ///< \brief Wang landau method (wl)

    double *shared_weights;           ///< \brief Array of weights for wl method
    long   *shared_hist;              ///< \brief Array of histogram for wl method
    double *shared_A_min_wmin;             ///< \brief Alpha, min, wmin

    long   length[2];          ///< \brief Length of above arrays
    double dorder[2];          ///< \brief Increments of order parameter
    double minorder[2];        ///< \brief Minimum order parameter
    double alpha;              ///< \brief Current modifier of weights

    long   currorder[2];       ///< \brief Walue of current order parameter
    long   neworder[2];        ///< \brief wl order parameter in new step
    long   max;                ///< \brief wl maximum of histogram
    long   min;                ///< \brief wl minimum of histogram
    double wmin;               ///< \brief weights minimum
    int    wlmdim;             ///< \brief Dimensionality of wang landau
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

    bool update(double temper) {
        if(mpirank == 0) {
            min = shared_hist[0];
            max = shared_hist[0]; // only local
            int j;
            for (int i=0;i < length[0];i++) {
                j=0;
                if ( shared_hist[i+j*length[0]] > max ) max = shared_hist[i+j*length[0]];
                if ( shared_hist[i+j*length[0]] < min ) min = shared_hist[i+j*length[0]];
                for (j=1;j < length[1];j++) {
                    if ( shared_hist[i+j*length[0]] > max ) max = shared_hist[i+j*length[0]];
                    if ( shared_hist[i+j*length[0]] < min ) min = shared_hist[i+j*length[0]];
                }
            }
            if ( min > WL_MINHIST ) {
                if ( temper * log(max/min) < WL_GERR ) {
                    /*DEBUG
                      for (i=1;i<wl.length;i++) {
                      printf (" %15.8e %15ld %15.8f\n",wl.weights[i],wl.hist[i],pvec[0].pos.z);
                      fflush(stdout);
                      }
                     */
                    if ( alpha < WL_ALPHATOL)
                        return true;
                    alpha/=2;
                    printf("%f \n", alpha);
                    fflush (stdout);
                    wmin = shared_weights[0];

                    for (int i=0;i < length[0];i++) {
                        j=0;
                        shared_hist[i+j*length[0]] = 0;
                        shared_weights[i+j*length[0]] -= wmin;
                        for (j=1;j < length[1];j++) {
                            shared_hist[i+j*length[0]] = 0;
                            shared_weights[i+j*length[0]] -= wmin;
                        }
                    }

                }
            }
            shared_A_min_wmin[0] = alpha;
            shared_A_min_wmin[1] = min;
            shared_A_min_wmin[2] = wmin;
        }
#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
        if(mpirank != 0){
            alpha = shared_A_min_wmin[0];
            min = shared_A_min_wmin[1];
            wmin = shared_A_min_wmin[2];
        }
#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
        return false;
    }

    inline double runRot(int& reject, int target) {
        double wlener = 0.0;

        for(int wli=0;wli<wlmdim;wli++) {
            switch(wlm[wli]) {
            case 3:     zOrient(wli, target);   break;
            default:    def(wli);               break;
            }
            if ( (neworder[wli] < 0) || (neworder[wli] >= length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += shared_weights[neworder[0]+neworder[1]*length[0]] - shared_weights[currorder[0]+currorder[1]*length[0]];
        }

        return wlener;
    }

    double runChainRot(int& reject, Particle orig[MAXCHL], double& radiusholemax_orig, Molecule& target) {
        double wlener = 0.0;

        for (int wli=0;wli<wlmdim;wli++) {
            switch (wlm[wli]) {
            case 1:
                if (target[0] == 0) zDir(wli);
                else def(wli);
                break;
            case 2:     holeXYPlane(wli, orig, target);                 break;
            case 3:     zOrient(wli, target[0]);                        break;
            case 4:     neworder[wli] = twoPartDist(wli);               break;
            case 5:     poreZCM(wli, radiusholemax_orig);               break;
            case 6:     poreZ0(wli, radiusholemax_orig,target, orig);   break;
            case 7:     partInContact(wli, target, orig);               break;
            default:    def(wli);                                       break;
            }
            if ( (neworder[wli] < 0) || (neworder[wli] >= length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += shared_weights[neworder[0]+neworder[1]*length[0]] - shared_weights[currorder[0]+currorder[1]*length[0]];
        }

        return wlener;
    }

    double runPress(int& reject, double& radiusholemax_orig, bool case2=false) {
        double wlener = 0.0;

        for (int wli=0;wli<wlmdim;wli++) {
            switch (wlm[wli]) {
                case 1: // press case 2 z constant, goto default
                    if(case2)
                        def(wli);
                    else
                        zDir(wli);
                    break;
                case 2:     holeXYPlane(wli);                   break;
                case 4:     neworder[wli] = twoPartDist(wli);   break;
                case 5:     poreZCM(wli, radiusholemax_orig);   break;
                case 6:     poreZ0(wli, radiusholemax_orig);    break;
                case 7:     partInContact(wli);                 break;
                default:    def(wli);                           break;
            }
            if ( (neworder[wli] < 0) || (neworder[wli] >= length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener = shared_weights[neworder[0]+neworder[1]*length[0]] - shared_weights[currorder[0]+currorder[1]*length[0]];
        }

        return wlener;
    }

    double run(int& reject, Particle orig[MAXCHL], Vector& cm, double& radiusholemax_orig, Molecule& target) { // for part and chain displace
        double wlener = 0.0;
        for(int wli=0;wli<wlmdim;wli++) {
            switch (wlm[wli]) {
            case 1:     zDir(wli, cm);                                  break;
            case 2:     holeXYPlane(wli, orig, target);                 break;
            case 4:     neworder[wli] = twoPartDist(wli);               break;
            case 5:     poreZCM(wli, cm, radiusholemax_orig);           break;
            case 6:     poreZ0(wli, radiusholemax_orig, target, orig);  break;
            case 7:     partInContact(wli, target, orig);               break;
            default:    def(wli);                                       break;
            }
            if ( (neworder[wli] < 0) || (neworder[wli] >= length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += shared_weights[neworder[0]+neworder[1]*length[0]] - shared_weights[currorder[0]+currorder[1]*length[0]];
        }

        return wlener;
    }


    double runSwitch(int& reject, double& radiusholemax_orig) {
        double wlener = 0.0;
        for (int wli=0;wli < wlmdim; wli++) {
            switch (wlm[wli]) {
            //case 1: sim->wl.neworder = z_order(&sim->wl, conf,wli); break;
            case 2:     holeXYPlane(wli);                   break;
            //case 4: sim->wl.neworder = twopartdist(&sim->wl,conf,wli); break;
            case 5:     poreZCM(wli, radiusholemax_orig);   break;
            case 6:     poreZ0(wli, radiusholemax_orig);    break;
            case 7:     partInContact(wli);                 break;
            default:    def(wli);                           break;
            }
            if ( ( neworder[wli] < 0) || (neworder[wli] >= length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += shared_weights[neworder[0]+neworder[1]*length[0]] - shared_weights[currorder[0]+currorder[1]*length[0]];
        }

        return wlener;
    }

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
     * @brief wlwrite
     * @param filename
     * @return
     */
    int write(char filename[30]);

    /**
     * @brief wlreject
     * @param oldlength
     */
    void reject(long oldlength, int wlm[2]) {
        if ( wlm[0] > 0 ) {
            shared_weights[currorder[0]+currorder[1]*length[0]] -= alpha;
            shared_hist[currorder[0]+currorder[1]*length[0]]++;
            if ( (wlm[0] == 2) || (wlm[1] == 2) )
                this->mesh = this->origmesh;
            if ( (wlm[0] == 5) || (wlm[1] == 5)||(wlm[0] == 6) || (wlm[1] == 6) ) {
                longarrayCpy(&radiushole,&radiusholeold,radiusholemax,oldlength);
                radiusholemax = oldlength;
            }
            partincontact = partincontactold;
        }
    }

    /**
     * @brief wlaccept
     * @param wlm
     */
    inline void accept(int wlm) {
        if ( wlm > 0 ) {
            for (int i=0;i<2;i++)
                currorder[i] = neworder[i];
            shared_weights[ currorder[0] + currorder[1] * length[0]] -= alpha;
            shared_hist[ currorder[0] + currorder[1] * length[0]]++;
        }
    }

private:

    // case:1
    inline void zDir(int wli, Vector& cm) {
        conf->syscm.x += cm.x / conf->sysvolume;
        conf->syscm.y += cm.y / conf->sysvolume;
        conf->syscm.z += cm.z / conf->sysvolume;
        neworder[wli] = zOrder(wli);
    }

    inline void zDir(int wli) {
        neworder[wli] = zOrder(wli);
    }

    // case:2
    inline void holeXYPlane(int wli, Particle orig[MAXCHL], Molecule& target) {
        origmesh = mesh;
        neworder[wli] = meshOrderMoveMolecule(target, orig, &mesh, wli);
        assert(target.size() == 1 || meshOrderMoveChain(target, &mesh, conf->pvec.size(), orig, wli) == neworder[wli]); // test chain
        //assert(target.size() > 1 || ); // test monomers
    }

    inline void holeXYPlane(int wli) {
        origmesh = mesh;
        neworder[wli] = (long) (mesh.meshInit(wl_meshsize, conf->pvec.size(), wlmtype, conf->geo.box, &conf->pvec) - minorder[wli]);
    }


    // case: 3
    inline void zOrient(int wli, int target) {
        if (target == 0) neworder[wli] = (long) floor( (conf->pvec[0].dir.z - minorder[wli])/ dorder[wli] );
        else neworder[wli] = currorder[wli];
    }

    //case:5
    inline void poreZCM(int wli, Vector& cm, double& radiusholemax_orig) {
        radiusholemax_orig = radiusholemax;
        conf->syscm.x += cm.x / conf->sysvolume;
        conf->syscm.y += cm.y / conf->sysvolume;
        conf->syscm.z += cm.z / conf->sysvolume;
        longarrayCpy(&radiusholeold,&radiushole,radiusholemax,radiusholemax);
        neworder[wli] = radiusholeAll(wli,&(conf->syscm));
    }

    inline void poreZCM(int wli, double& radiusholemax_orig) {
        radiusholemax_orig = radiusholemax;
        longarrayCpy(&radiusholeold,&radiushole,radiusholemax,radiusholemax);
        neworder[wli] = radiusholeAll(wli,&(conf->syscm));
    }

    // case:6
    inline void poreZ0(int wli, double& radiusholemax_orig, Molecule& target, Particle orig[MAXCHL]) {
        radiusholemax_orig = radiusholemax;
        longarrayCpy(&radiusholeold, &radiushole, radiusholemax, radiusholemax);
        if ( target[0] == 0 )
            neworder[wli] = radiusholeAll(wli,&(conf->pvec[0].pos));
        else {
            neworder[wli] = radiusholeOrderMoveMolecule(target, orig,wli,&(conf->pvec[0].pos));
        }
    }

    inline void poreZ0(int wli, double& radiusholemax_orig) {
        radiusholemax_orig = radiusholemax;
        longarrayCpy(&radiusholeold, &radiushole, radiusholemax, radiusholemax);
        neworder[wli] = radiusholeAll(wli,&(conf->pvec[0].pos));
    }

    // case:7
    inline void partInContact(int wli, Molecule& target, Particle orig[MAXCHL]) {
        partincontactold = partincontact;
        if ( target[0] == 0 )
            neworder[wli] = contParticlesAll(wli);
        else {
            neworder[wli] = contParticlesMoveMolecule(target,orig,wli);
        }
    }

    inline void partInContact(int wli) {
        partincontactold = partincontact;
        neworder[wli] = contParticlesAll(wli);
    }


    inline void def(int wli) {
        neworder[wli] = currorder[wli];
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief twopartdist
     * @param wli
     * @return 2 particles distance
     */
    inline long twoPartDist(int wli) {
        Vector r_cm = conf->geo.image(&conf->pvec[0].pos, &conf->pvec[1].pos);
        return (long) ceil( (( sqrt(r_cm.x*r_cm.x + r_cm.y*r_cm.y) ) - minorder[wli]) / dorder[wli]  );
    }


    /**
     * @brief radiusholeorder_movechain return change in order parameter when chain moves
     * @param chain
     * @param chorig
     * @param wli
     * @param position
     * @return
     */
    long radiusholeOrderMoveMolecule(Molecule molecule, Particle chorig[MAXCHL], int wli, Vector *position) {
        long nr,oor; // position in radiushole
        double rx,ry,z;
        bool oz,nz, change=false;;

        if(molecule.size() == 1) {

            if ( conf->pvec[molecule[0]].type != wlmtype )
                return currorder[wli];

            z=conf->pvec[molecule[0]].pos.z - position->z; //if above position

            if (z-anInt(z) < 0)
                nz = false;
            else
                nz=true;

            z=chorig[0].pos.z - position->z;  /*if above position*/

            if (z-anInt(z) < 0)
                oz = false;
            else
                oz=true;

            if ( !(nz) && !(oz) )
                return currorder[wli];

            rx = conf->geo.box.x * (conf->pvec[molecule[0]].pos.x - anInt(conf->pvec[molecule[0]].pos.x));
            ry = conf->geo.box.y * (conf->pvec[molecule[0]].pos.y - anInt(conf->pvec[molecule[0]].pos.y));
            nr = radiusholePosition(sqrt(rx*rx+ry*ry),wli);

            if (nr < 0)
                return -100;
            //particle move over radius bins*/
            if (nz) {
                radiushole[nr]++;
            }
            if (oz) {
                rx = conf->geo.box.x * (chorig[0].pos.x - anInt(chorig[0].pos.x));
                ry = conf->geo.box.y * (chorig[0].pos.y - anInt(chorig[0].pos.y));
                oor = radiusholePosition(sqrt(rx*rx+ry*ry),wli);
                radiushole[oor]--;
                if ( radiushole[oor] < 0 ) {
                    printf ("Error(single particle move): trying to make number of beads in radiuspore smaller than 0 at position %ld\n",oor);
                    radiusholePrint(radiushole,radiusholemax);
                    fflush(stdout);
                }
                if (radiushole[oor] ==0)
                    return radiusholeOrder();
            }

            if ( (nz) && (radiushole[nr] ==1) )  {
                return radiusholeOrder();
            }

            return currorder[wli];
        }

        rx=0;
        for(unsigned int i=0; i<molecule.size(); i++ ) {
            if ( conf->pvec[molecule[i]].type == wlmtype ) {
                z=conf->pvec[molecule[i]].pos.z - position->z; // if above system CM
                if (z-anInt(z) > 0) {
                    rx = conf->geo.box.x * (conf->pvec[molecule[i]].pos.x - anInt(conf->pvec[molecule[i]].pos.x));
                    ry = conf->geo.box.y * (conf->pvec[molecule[i]].pos.y - anInt(conf->pvec[molecule[i]].pos.y));
                    nr = radiusholePosition(sqrt(rx*rx+ry*ry),wli);
                    if (nr < 0)
                        return -100;
                    radiushole[nr]++;
                    if ( radiushole[nr] == 1 ) change = true;
                }
            }
        }
        for(unsigned int i=0; i<molecule.size(); i++ ) {
            if ( conf->pvec[molecule[i]].type == wlmtype ) {
                z=chorig[i].pos.z - position->z; // if above system CM
                if (z-anInt(z) > 0) {
                    rx = conf->geo.box.x * (chorig[i].pos.x - anInt(chorig[i].pos.x));
                    ry = conf->geo.box.y * (chorig[i].pos.y - anInt(chorig[i].pos.y));
                    nr = radiusholePosition(sqrt(rx*rx+ry*ry),wli);
                    radiushole[nr]--;
                    if ( radiushole[nr] < 0 ) {
                        printf ("Error (chainmove): trying to make number of beads in radiuspore smaller than 0 at position %ld\n",nr);
                        radiusholePrint(radiushole,radiusholemax);
                        fflush(stdout);
                    }
                    if ( radiushole[nr] == 0 ) change = true;
                }
            }
        }

        if ( change ) {
            return radiusholeOrder();
        }
        return currorder[wli];
    }


    /**
     * @brief contParticle_movechain return change in order parameter when chain moves
     * @param chain
     * @param conf
     * @param sim
     * @param chorig
     * @param wli
     * @return
     */
    long contParticlesMoveMolecule(Molecule molecule, Particle chorig[MAXCHL], int wli) {

        for(unsigned int i=0; i<molecule.size(); i++ ) {
            if ( conf->pvec[molecule[i]].type == wlmtype ) {
                if ( particlesInContact (&(conf->pvec[molecule[i]].pos)) )
                    partincontact++;
            }
        }
        for(unsigned int i=0; i<molecule.size(); i++ ) {
            if ( conf->pvec[molecule[i]].type == wlmtype ) {
                if ( particlesInContact (&(chorig[i].pos)) )
                    partincontact--;
            }
        }

        return contParticlesOrder(wli);
    }

    /**
     * @brief meshOrderMoveMolecule
     * @param target
     * @param chorig
     * @param mesh
     * @param wli
     * @return
     */
    long meshOrderMoveMolecule(Molecule& target, Particle chorig[], Mesh *mesh, int wli);



    /**
     * @brief z_order
     * @param wli
     * @return
     */
    inline long zOrder(int wli) {
        return (long) ceil( ((conf->pvec[0].pos.z - conf->syscm.z) * conf->geo.box.z- minorder[wli]) / dorder[wli]  );
        //return (long) ceil( ((conf->pvec[0].pos.z ) * conf->geo.box.z- minorder[wli]) / dorder[wli]  ); // for testing
    }

    /**
     * @brief wlend
     * @return
     */
    int end() {
#ifdef ENABLE_MPI
        MPI_Win_free(&_winDouble);
        MPI_Win_free(&_winInt);
        MPI_Win_free(&_winShared);
#else
        free(shared_weights);
        shared_weights = NULL;
        free(shared_hist);
        delete shared_A_min_wmin;
        shared_hist = NULL;
        shared_A_min_wmin = NULL;
#endif
        return 0;
    }

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
    long contParticlesAll(int wli) {
        partincontact = 0;

        for (int i=1; i< (long)conf->pvec.size(); i++) {
            /*calculate position of particle and add it if in contact */
            if ( conf->pvec[i].type == wlmtype ) {
                if ( particlesInContact (&(conf->pvec[i].pos)) )
                    partincontact++;
            }
        }

        return contParticlesOrder(wli);
    }


    /**
     * @brief radiushole_position return order of given radius
     * @param radius
     * @param wli
     * @return
     */
    inline long radiusholePosition(double radius, int wli) {
        return (long) ceil( ( radius - minorder[wli]) / dorder[wli]  );
    }

    /**
     * @brief radiushole_order return current bin of free radius
     * @return
     */
    inline long radiusholeOrder() {
        for (int i=0;i<radiusholemax-3;i++){
            if ((radiushole[i] >0 ) && (radiushole[i+1] >0 ) && (radiushole[i+2] >0 ) && (radiushole[i+3] >0 ))
                return i-1;
        }
        return -100;
    }

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
    inline bool particlesInContact(Vector *vec) {
        double x,y,z;

        x = vec->x - conf->pvec[0].pos.x;
        y = vec->y - conf->pvec[0].pos.y;
        z = vec->z - conf->pvec[0].pos.z;

        x = conf->geo.box.x * (x - anInt(x));
        y = conf->geo.box.y * (y - anInt(y));
        z = conf->geo.box.z * (z - anInt(z));

        if ( x*x + y*y + z*z < WL_CONTACTS) {
            return true;
        }
        else {
            return false;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                          Deprecated, for debug only                                 ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    long meshOrderMoveOne(Vector oldpos, Vector newpos, Mesh *mesh, long npart, long target, int wli);

    long meshOrderMoveChain(Molecule chain, Mesh *mesh, long npart, Particle chorig[], int wli);

};

#endif // WANGLANDAU_H
