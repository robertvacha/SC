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
        Vector r_cm;
        r_cm.x = part1->pos.x - part2->pos.x;
        r_cm.y = part1->pos.y - part2->pos.y;
        r_cm.z = part1->pos.z - part2->pos.z;

        if ( r_cm.x < 0  )
            r_cm.x = box.x * (r_cm.x - (double)( (long)(r_cm.x-0.5) ) );
        else
            r_cm.x = box.x * (r_cm.x - (double)( (long)(r_cm.x+0.5) ) );
        if ( r_cm.y < 0  )
            r_cm.y = box.y * (r_cm.y - (double)( (long)(r_cm.y-0.5) ) );
        else
            r_cm.y = box.y * (r_cm.y - (double)( (long)(r_cm.y+0.5) ) );
        if ( r_cm.z < 0  )
            r_cm.z = box.z * (r_cm.z - (double)( (long)(r_cm.z-0.5) ) );
        else
            r_cm.z = box.z * (r_cm.z - (double)( (long)(r_cm.z+0.5) ) );

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
        for (unsigned int i=0; i < pvec.size(); i++) {
            fprintf (outfile, "%15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e %d\n",
                    box.x * ((pvec[i].pos.x) - anInt(pvec[i].pos.x)),
                    box.y * ((pvec[i].pos.y) - anInt(pvec[i].pos.y)),
                    box.z * ((pvec[i].pos.z) - anInt(pvec[i].pos.z)),
                    pvec[i].dir.x, pvec[i].dir.y, pvec[i].dir.z,
                    pvec[i].patchdir[0].x, pvec[i].patchdir[0].y, pvec[i].patchdir[0].z,
                    pvec[i].switched);
        }
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
