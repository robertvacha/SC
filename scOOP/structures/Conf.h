#ifndef CONF_H
#define CONF_H

#include "topo.h"
#include "structures.h"
#include "moleculeparams.h"
#include "particle.h"
#include <assert.h>
#include <map>

extern Topo topo;

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
        //assert(false && "ChainN not found");
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

class GeoBase {
public:
    Vector box;                             ///< \brief Box size */

    GeoBase(){}

    virtual void usePBC(Particle *pos) = 0;
    virtual Vector image(Vector* r1, Vector* r2) = 0;
    virtual bool boundaryOverlap(Vector *pos) {return false;}       // all - for insertion because of uniformity

    virtual Vector randomPos() =0;
    virtual double volume() = 0;
};

class Cuboid: public GeoBase {
public:

    Cuboid(){}
    Cuboid(Vector box) {
        this->box = box;
        cout << "Box: " << this->box.info() << endl;
    }

    inline double volume() {
        return box.x*box.y*box.z;
    }

    /**
     * @brief use of periodic boundary conditions
     * @param pos
     * @param pbc
     */
    void usePBC(Particle *part) {
        do {
            part->pos.x += box.x;
        } while (part->pos.x < 0.0);
        do {
            part->pos.x -= box.x;
        } while (part->pos.x > box.x);

        do {
            part->pos.y += box.y;
        } while (part->pos.y < 0.0);
        do {
            part->pos.y -= box.y;
        } while (part->pos.y > box.y);

        do {
            part->pos.z += box.z;
        } while (part->pos.z < 0.0);
        do {
            part->pos.z -= box.z;
        } while (part->pos.z > box.z);
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
        Vector r_cm( r1->x - r2->x, r1->y - r2->y, r1->z - r2->z );

        Vector temp(r_cm.x + 6755399441055744.0, r_cm.y + 6755399441055744.0, r_cm.z + 6755399441055744.0);

        r_cm.x = box.x * (r_cm.x - static_cast<double>(reinterpret_cast<int&>(temp.x) ) );
        r_cm.y = box.y * (r_cm.y - static_cast<double>(reinterpret_cast<int&>(temp.y) ) );
        r_cm.z = box.z * (r_cm.z - static_cast<double>(reinterpret_cast<int&>(temp.z) ) );

        return r_cm;
    }

    virtual Vector randomPos() {
        return Vector(ran2(), ran2(), ran2());
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
class Wedge : public GeoBase {
private:
    double angleSin;
    double angleCos;
    double angleRad;
    Vector scale;       // scale unitCube -> higher prob of insert into wedge
    Vector offset;      // offset scaled unitCube to wedge

    inline void rotateClockWise(Vector* pos) { // OK
        double x = pos->x;
        pos->x = x * angleCos + pos->y * angleSin;
        pos->y = pos->y * angleCos - x * angleSin;
    }

    inline void rotateCounterClock(Vector* pos) { // OK
        double x = pos->x;
        pos->x = x * angleCos - pos->y * angleSin;
        pos->y = pos->y * angleCos + x * angleSin;
    }

public:
    double innerR; // diameter = 2 * radius
    double outerR;
    double angleDeg;

    Wedge(){cout << "NO PARAMETERS GIVEN !!!" << endl;}
    Wedge(double z, double angle, double outerR, double innerR) : angleDeg(angle) {
        this->box = Vector(outerR, outerR, z);
        cout << "Wedge: box: " << this->box.info() << " Angle:" << angleDeg
             << " Inner radius: " << innerR << " Outer radius: " << outerR << endl;
        if(innerR > outerR) {
            cout << "Mixed Inner and outer radius!!!" << endl;
            exit(1);
        }
        if(0.0 <= angle && angle > 90) {
            cout << "Angle: " << angle << " must be between 0 and 90 degrees!!!" << endl;
            exit(1);
        }
        if(2*innerR*innerR < topo.sqmaxcut) {
            cout << "inner radius too small for interactions of species, see constructor Wedge" << endl;
            exit(1);
        }
        this->innerR = innerR/box.x;
        this->outerR = outerR/box.x;
        angleRad = angle*PI / 180.0;
        angleSin = sin(angleRad);
        angleCos = cos(angleRad);

        //cout << "Sin: " << angleSin << " Cos: " << angleCos << endl;

        scale = Vector(angleSin, angleCos, 0.0);
        offset = Vector(0.0, 1.0 - angleCos, 0.0);
    }

    inline double volume() {
        return box.x*box.x*box.z*angleRad*0.5*(outerR*outerR - innerR*innerR);
    }

    inline bool boundaryOverlap(Vector *pos) { // OK
        double radius = pos->x*pos->x + pos->y*pos->y;
        if(radius > outerR*outerR || radius < innerR*innerR)
            return true;
        // check plane YZ (angle 0) -> x go negative -> rotate clockwise
        if(pos->x < 0.0)
            return true;
        // check angle
        if(atan2(pos->x, pos->y) > angleRad)
            return true;
        return false;
    }

    virtual void usePBC(Particle *part) {
        // check plane YZ (angle 0) -> x go negative -> rotate clockwise
        if(part->pos.x < 0.0) {
            rotateClockWise(&part->pos);
            part->pscRotate(Vector(0,0,1), angleRad, true);
        }
        // check angle
        if(atan2(part->pos.x, part->pos.y) > angleRad) {
            rotateCounterClock(&part->pos);
            part->pscRotate(Vector(0,0,1), angleRad, false);
        }

        /*do { // unit lenght scaling for z axis
            (*pos).z += box->z;
        } while ((*pos).z < 0.0);
        do {
            (*pos).z -= box->z;
        } while ((*pos).z > box->z);*/
    }

    virtual Vector image(Vector* r1, Vector* r2) {
        Vector r_cm, r_cm2 = *r1;
        r_cm.x = r1->x - r2->x;
        r_cm.y = r1->y - r2->y;
        r_cm.z = r1->z - r2->z;

        double temp =  r_cm.z + 6755399441055744.0;
        r_cm.z = box.z * (r_cm.z - static_cast<double>(reinterpret_cast<int&>(temp) ) );

        if(atan2(r_cm2.x, r_cm2.y) < angleRad*0.5)
            rotateClockWise(&r_cm2);
        else rotateCounterClock(&r_cm2);

        r_cm2.x = r_cm2.x - r2->x;
        r_cm2.y = r_cm2.y - r2->y;
        r_cm2.z = r_cm2.z - r2->z;

        temp =  r_cm2.z + 6755399441055744.0;
        r_cm2.z = box.z * (r_cm2.z - static_cast<double>(reinterpret_cast<int&>(temp) ) );

        if(r_cm.dot(r_cm) < r_cm2.dot(r_cm2) )
            return r_cm;
        else return r_cm2;
    }

    virtual Vector randomPos() {
        Vector pos;
        pos.randomUnitCube();
        pos.x = pos.x * scale.x;
        pos.y = pos.y * scale.y + offset.y;
        while(boundaryOverlap(&pos)) {
            pos.randomUnitCube();
            pos.x = pos.x * scale.x;
            pos.y = pos.y * scale.y + offset.y;
        }
        return pos;
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

    double sysvolume;                       ///< \brief Something like total mass -> for each += massOfParticle
    Vector syscm;                           ///< \brief System center of mass

    GroupList pvecGroupList;
    GroupList poolGroupList;

    bool pairlist_update;

#ifdef WEDGE
    Wedge geo;
#else
    Cuboid geo;
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

    inline double distSq(Particle* part1, Particle* part2) {
        Vector r_cm = geo.image(&part1->pos, &part2->pos);
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
            fprintf (outfile, "%15.6e %15.6e %15.6e   %15.6e %15.6e %15.6e   %15.6e %15.6e %15.6e %d\n",
                    geo.box.x * ((pvec[i].pos.x) - anInt(pvec[i].pos.x)),
                    geo.box.y * ((pvec[i].pos.y) - anInt(pvec[i].pos.y)),
                    geo.box.z * ((pvec[i].pos.z) - anInt(pvec[i].pos.z)),
                    pvec[i].dir.x, pvec[i].dir.y, pvec[i].dir.z,
                    pvec[i].patchdir[0].x, pvec[i].patchdir[0].y, pvec[i].patchdir[0].z,
                    pvec[i].switched);
        }
#else
#ifdef WEDGE
        for (unsigned int i=0; i < pvec.size(); i++) {
            fprintf (outfile, "%15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e %d\n",
                    geo.box.x * ((pvec[i].pos.x)),
                    geo.box.y * ((pvec[i].pos.y)),
                    geo.box.z * ((pvec[i].pos.z) - anInt(pvec[i].pos.z)),
                    pvec[i].dir.x, pvec[i].dir.y, pvec[i].dir.z,
                    pvec[i].patchdir[0].x, pvec[i].patchdir[0].y, pvec[i].patchdir[0].z,
                    pvec[i].switched);
        }
#else
        for (unsigned int i=0; i < pvec.size(); i++) {
            fprintf (outfile, "%15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e %d\n",
                    geo.box.x * ((pvec[i].pos.x) - anInt(pvec[i].pos.x)),
                    geo.box.y * ((pvec[i].pos.y) - anInt(pvec[i].pos.y)),
                    geo.box.z * ((pvec[i].pos.z) - anInt(pvec[i].pos.z)),
                    pvec[i].dir.x, pvec[i].dir.y, pvec[i].dir.z,
                    pvec[i].patchdir[0].x, pvec[i].patchdir[0].y, pvec[i].patchdir[0].z,
                    pvec[i].switched);
        }
#endif
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
