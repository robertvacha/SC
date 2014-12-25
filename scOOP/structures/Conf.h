#ifndef CONF_H
#define CONF_H

#include "topo.h"
#include "structures.h"
#include "moleculeparams.h"
#include "particle.h"
#include <assert.h>
#include <map>

using namespace std;

class ConList {
public:
    ConList() {
        conlist[0]=NULL;
        conlist[1]=NULL;
        conlist[2]=NULL;
        conlist[3]=NULL;
    }
    Particle* conlist[4];    ///< \brief Connectivity list, we have connection to tail and head and secon neighbours so far
};


class GroupList {
public:
    /// @brief first Index(of pvec) of first particle of molecule type, array over molecular types
    /// on [molTypeCount] == size of vector
    int first[MAXMT+1]; ///<
    int molSize[MAXMT]; ///< \brief Number of particles per molecule of moleculeType, array over molecular types

    /// \brief chainIndex of first chain of molType
    /// on [molTypeCount] == number of all chains
    int firstChain[MAXMT+1];

    int molTypeCount;   ///< \brief Count of molecular types in use

    GroupList() {
        for(int i=0; i<MAXMT; i++) {
            first[i] = -1; molSize[i] = -1; firstChain[i]=0;
        }
    }

    inline int getChainCount() {return firstChain[molTypeCount];}

    int vecSize() {return first[molTypeCount];}

    void calcChainCount() {
        firstChain[molTypeCount]=0;
        for(int molType=0; molType<=molTypeCount; molType++) {
            if(molSize[molType] > 1) {
                firstChain[molType] = firstChain[molTypeCount];
                firstChain[molTypeCount] += molCountOfType(molType);
            } else {
                firstChain[molType] = firstChain[molTypeCount];
            }
        }
    }

    inline int getChain(int chainN, int pos) {
        for(int type=0; type<=molTypeCount; type++) {
            if(firstChain[type] > chainN) {
                type--;
                if(molSize[type] <= pos) return -1;
                chainN -= firstChain[type];
                return first[type]+chainN*molSize[type] + pos;
            }
        }
        assert(false && "ChainN not found");
        return -1;
    }

    vector<int> getChain(int chainN) {
        vector<int> ret;
        for(int type=0; type<=molTypeCount; type++) {
            if(firstChain[type] > chainN) {
                type--;
                chainN -= firstChain[type];
                chainN *= molSize[type];
                for(int i=0; i<molSize[type]; i++)
                    ret.push_back(first[type]+chainN + i);
                return ret;
            }
        }
        assert(false && "ChainN not found");
        return ret;
    }

    void info() {
        cout << "number of types: " << molTypeCount << endl;
        cout << "first:";
        for(int i=0; i<=molTypeCount; i++)
             cout << " " << first[i];
        cout << endl;
    }

    /**
     * @brief Converts molID of molType to pvec Index
     * molID -> starts from 0 for each molType
     * @return Index of first particle of molecule
     */
    inline int getStoreIndex(int molType, int molID) const {
        return first[molType] + molSize[molType]*molID;
    }

    inline int getInsertIndex(int molType) const {
        return first[molType+1];
    }

    /**
     * @param molType Type of molecule
     * @return Number of molecules a given type
     */
    int molCountOfType(int molType) {
        assert( ( first[molType+1] -  first[molType]) /  molSize[molType] >= 0);
        return ( first[molType+1] -  first[molType]) /  molSize[molType];
    }

    void insertMolecule(int molType) {
        for(int i= molType+1; i<=molTypeCount; i++)
            first[i] += molSize[molType];

        if(molSize[molType] > 1)
            calcChainCount();

        assert(checkConsistency());
    }

    void deleteMolecule(int molType) {
        for(int i=molType+1; i<=molTypeCount; i++)
            first[i] -= molSize[molType];

        if(molSize[molType] > 1)
            calcChainCount();

        assert(checkConsistency());
    }

#ifndef NDEBUG
    int checkConsistency() {
        for(int i=0; i<molTypeCount; i++) {
            if(first[i+1] >= first[i])return 1;
        }
        return 0;
    }
#endif

};

class PbcBase {
public:
    Vector* box;

    PbcBase(){}
    PbcBase(Vector* box) : box(box) {}

    virtual void usePBC(Vector *pos) = 0;
    virtual Vector image(Vector* r1, Vector* r2) = 0;
    virtual bool boundaryOverlap(Vector *pos) {return false;} // default - no solid boundaries
};

class Cuboid: public PbcBase {
public:

    Cuboid(){cout << "NO PARAMETERS GIVEN !!!" << endl;}
    Cuboid(Vector* box) : PbcBase(box) {}

    /**
     * @brief use of periodic boundary conditions
     * @param pos
     * @param pbc
     */
    void usePBC(Vector *pos) {
        do {
            (*pos).x += box->x;
        } while ((*pos).x < 0.0);
        do {
            (*pos).x -= box->x;
        } while ((*pos).x > box->x);

        do {
            (*pos).y += box->y;
        } while ((*pos).y < 0.0);
        do {
            (*pos).y -= box->y;
        } while ((*pos).y > box->y);

        do {
            (*pos).z += box->z;
        } while ((*pos).z < 0.0);
        do {
            (*pos).z -= box->z;
        } while ((*pos).z > box->z);
    }

    /**
     * @brief Returns the vector pointing from the centre of mass of particle 2 to the
       centre of mass of the closest image of particle 1.
     * @param r1
     * @param r2
     * @param box
     * @return
     */
    inline Vector image(Vector* r1, Vector* r2) {
        /*double x = r1->x - r2->x, y = r1->y - r2->y, z = r1->z - r2->z;
        return Vector( box->x * (x - anInt(x)),
                      box->y * (y - anInt(y)),
                      box->z * (z - anInt(z)) );*/

        Vector r_cm;
        r_cm.x = r1->x - r2->x;
        r_cm.y = r1->y - r2->y;
        r_cm.z = r1->z - r2->z;

        Vector temp(r_cm.x + 6755399441055744.0, r_cm.y + 6755399441055744.0, r_cm.z + 6755399441055744.0);

        r_cm.x = box->x * (r_cm.x - static_cast<double>(reinterpret_cast<int&>(temp.x) ) );
        r_cm.y = box->y * (r_cm.y - static_cast<double>(reinterpret_cast<int&>(temp.y) ) );
        r_cm.z = box->z * (r_cm.z - static_cast<double>(reinterpret_cast<int&>(temp.z) ) );

        return r_cm;
    }

};

/**
 * @brief The Wedge class
 * Solid rotational boundaries - inner and outer radius - around Z axis
 * Unit Z pbc
 * rotational XY pbc
 *
 * CHAINS NOT SUPPORTED
 *
 */
class Wedge : public PbcBase {
private:
    double angleSin;
    double angleCos;
    double angleRad;

    inline void rotateClockWise(Vector* pos) {
        pos->x = pos->x * angleCos - pos->y * (angleSin);
        pos->y = pos->x * angleSin + pos->y * angleCos;
    }

    inline void rotateCounterClock(Vector* pos) {
        pos->x = pos->x * angleCos + pos->y * (angleSin);
        pos->y = pos->y * angleCos - pos->x * angleSin;
    }

public:
    double innerR; // diameter = 2 * radius
    double outerR;
    double angleDeg;

    Wedge(){cout << "NO PARAMETERS GIVEN !!!" << endl;}
    Wedge(Vector* box, double angle, double outerR, double innerR) : PbcBase(box), angleDeg(angle) {
        if(box->x != box->y) {
            cout << "Box size X and Y must be the same for WEDGE BOX!!!" << endl;
            exit(1);
        }
        if(innerR > outerR) {
            cout << "Mixed Inner and outer radius!!!" << endl;
            exit(1);
        }
        this->innerR = innerR/box->x;
        this->outerR = outerR/box->x;
        angleRad = angle*PI / 180.0;
        angleSin = sin(angleRad);
        angleCos = cos(angleRad);
    }

    inline bool boundaryOverlap(Vector *pos) {
        // check if between inner and outer radius
        double radius = pos->x*pos->x + pos->y*pos->y;
        if(radius < outerR && radius > innerR)
            return false;
        // check plane YZ (angle 0) -> x go negative -> rotate clockwise
        if(pos->x < 0.0)
            return true;
        // check angle
        if(PI - atan2(pos->x, pos->y) > angleRad)
            return true;
        return false;
    }

    virtual void usePBC(Vector *pos) {
        // check plane YZ (angle 0) -> x go negative -> rotate clockwise
        if(pos->x < 0.0)
            rotateClockWise(pos);
        // check angle
        if(PI - atan2(pos->x, pos->y) > angleRad)
            rotateCounterClock(pos);

        do {
            (*pos).z += box->z;
        } while ((*pos).z < 0.0);
        do {
            (*pos).z -= box->z;
        } while ((*pos).z > box->z);
    }

    virtual Vector image(Vector* r1, Vector* r2) {
        Vector r_cm, r_cm2 = *r1;
        r_cm.x = r1->x - r2->x;
        r_cm.y = r1->y - r2->y;
        r_cm.z = r1->z - r2->z;

        if ( r_cm.z < 0  )
            r_cm.z = box->z * (r_cm.z - (double)( (long)(r_cm.z-0.5) ) );
        else
            r_cm.z = box->z * (r_cm.z - (double)( (long)(r_cm.z+0.5) ) );

        if(PI - atan2(r_cm2.x, r_cm2.y) < angleRad*0.5)
            rotateClockWise(&r_cm2);
        else rotateCounterClock(&r_cm2);

        r_cm2.x = r_cm2.x - r2->x;
        r_cm2.y = r_cm2.y - r2->y;
        r_cm2.z = r_cm2.z - r2->z;

        if ( r_cm2.z < 0  )
            r_cm2.z = box->z * (r_cm2.z - (double)( (long)(r_cm2.z-0.5) ) );
        else
            r_cm2.z = box->z * (r_cm2.z - (double)( (long)(r_cm2.z+0.5) ) );

        if(r_cm.dot(r_cm) < r_cm2.dot(r_cm2) )
            return r_cm;
        else return r_cm2;
    }

};

/**
 * @brief Configuration of the system
 */
class Conf {
public:  
    std::vector<Particle > pvec;  ///< \brief Main store of all particles, grouped by Molecular types
    std::vector<Neighbors > neighborList;
    std::vector<ConList > conlist;
    std::vector<Particle > pool; ///< \brief Store for chains for muVT insert of chain

    Vector box;                             ///< \brief Box size */
    double sysvolume;                       ///< \brief Something like total mass -> for each += massOfParticle
    Vector syscm;                           ///< \brief System center of mass

    GroupList pvecGroupList;
    GroupList poolGroupList;

    bool pairlist_update;

#ifdef WEDGE
    Wedge pbc;
#else
    Cuboid pbc;
#endif


public:
    /**
     * @brief Conf Constructor, initializing variables
     */
    Conf() : pairlist_update(false) {
        try{
            pvec.reserve(MAXN);
            conlist.reserve(MAXN);
            pool.reserve(MAXN);
            neighborList.reserve(MAXN);
        } catch(std::bad_alloc& bad) {
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for pvec, pool, conlist, neighborlist, conf inicializer");
            exit(1);
        }
    }

    ~Conf() {
        printf ("Deallocating Conf...\n");

        for(unsigned int i=0; i < neighborList.size(); i++){
            delete neighborList[i].neighborID;
            neighborList[i].neighborID = NULL;
        }

        /*if (conf.particle != NULL) // left just for peace of mind, before Particles* --> std::vector<Particles* >
            free(conf.particle);
        conf.particle = NULL;*/
    }

    double distSq(Particle* part1, Particle* part2) {
        Vector r_cm = pbc.image(&part1->pos, &part2->pos);
        return DOT(r_cm,r_cm);
    }

    double dist(Particle* part1, Particle* part2) {
        return sqrt(distSq(part1,part2));
    }

    void recalcConList();

    /**
     * @brief draw Dumps a configuration to the supplied file handle.
     * @param outfile
     */
    void draw(FILE *outfile) {
        //fprintf (outfile, "%15.8e %15.8e %15.8e\n", box.x, box.y, box.z);
#ifdef TESTING
        for (unsigned int i=0; i < pvec.size(); i++) {
            fprintf (outfile, "%15.7e %15.7e %15.7e   %15.7e %15.7e %15.7e   %15.7e %15.7e %15.7e %d\n",
                    box.x * ((pvec[i].pos.x) - anInt(pvec[i].pos.x)),
                    box.y * ((pvec[i].pos.y) - anInt(pvec[i].pos.y)),
                    box.z * ((pvec[i].pos.z) - anInt(pvec[i].pos.z)),
                    pvec[i].dir.x, pvec[i].dir.y, pvec[i].dir.z,
                    pvec[i].patchdir[0].x, pvec[i].patchdir[0].y, pvec[i].patchdir[0].z,
                    pvec[i].switched);
        }
#else
        for (unsigned int i=0; i < pvec.size(); i++) {
            fprintf (outfile, "%15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e %d\n",
                    box.x * ((pvec[i].pos.x) - anInt(pvec[i].pos.x)),
                    box.y * ((pvec[i].pos.y) - anInt(pvec[i].pos.y)),
                    box.z * ((pvec[i].pos.z) - anInt(pvec[i].pos.z)),
                    pvec[i].dir.x, pvec[i].dir.y, pvec[i].dir.z,
                    pvec[i].patchdir[0].x, pvec[i].patchdir[0].y, pvec[i].patchdir[0].z,
                    pvec[i].switched);
        }
#endif
    }

    /**
     * @brief Adds molecule to pvec, ensures Particle order, grouplist, conlist
     * @param molecule
     */
    void addMolecule(std::vector<Particle>* molecule);

    /**
     * @brief removeMolecule
     * @param molType
     * @param molID
     */
    void removeMolecule(int target, int size);

    std::vector<Particle> getRandomPoolConf(int molType);

    /**
     * @brief massCenter
     */
    void massCenter();

    /**
     * @brief partvecinit calculate vectors on particles for speedup
     */
    void partVecInit();

    /**
     * @brief Determines whether two particles overlap.
     * @param part1
     * @param part2
     * @param ia_params 1 if there is an overlap, 0 if not.
     * @return
     */
    int overlap(Particle* part1, Particle* part2, Ia_param ia_params[MAXT][MAXT]);

    double linemin(double criterion, double halfl)  {
        if      (criterion >=  halfl) { return  halfl; }
        else if (criterion >= -halfl) { return  criterion; }
        else                          { return -halfl; }
    }

    /**
     * @brief forbidden Checks for overlaps between particle "target" and the rest.
     * @param target
     * @param ia_params
     * @return Returns true if overlap detected, false otherwise.
     */
    bool overlapAll(Particle* target, Ia_param ia_params[MAXT][MAXT]);

    /**
     * @brief checkall Checks for overlaps between all pairs of particles.
     * @param ia_params
     * @return  Returns 1 if overlap detected, 0 otherwise.
     */
    int checkall(Ia_param ia_params[MAXT][MAXT]);
};


#endif // CONF_H
