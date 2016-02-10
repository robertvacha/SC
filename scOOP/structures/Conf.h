#ifndef CONF_H
#define CONF_H

#include "topo.h"
#include "structures.h"
#include "moleculeparams.h"
#include "particle.h"
#include <assert.h>
#include <map>
#include "geometry.h"

extern Topo topo;

using namespace std;

class ConList {
public:
    ConList() : conlist() {}
    Particle* conlist[4];    ///< \brief Connectivity list, we have connection to tail and head and secon neighbours so far

    inline Particle* operator[] (int j) {
        return conlist[j];
    }
};


/**
 * @brief The ParticleVector class
 */
class ParticleVector : public std::vector<Particle > {
public:
    /// @brief first Index(of pvec) of first particle of molecule type, array over molecular types
    /// on [molTypeCount] == size of vector
    int first[MAXMT+1];

    /// \brief chainIndex of first chain of molType
    /// on [molTypeCount] == number of all chains
    int chainCount[MAXMT+1];

    int molTypeCount;   ///< \brief Count of molecular types in use

    /***********************************************************************************************************/

    ParticleVector() {
        for(int i=0; i<MAXMT; i++) {
            first[i] = 0; chainCount[i]=0;
        }
    }

    inline int getChainCount() {return chainCount[molTypeCount];}

    void calcChainCount() {
        chainCount[molTypeCount]=0;
        for(int molType=0; molType<=molTypeCount; molType++) {
            if(topo.moleculeParam[molType].molSize() > 1) {
                chainCount[molType] = chainCount[molTypeCount];
                chainCount[molTypeCount] += molCountOfType(molType);
            } else {
                chainCount[molType] = chainCount[molTypeCount];
            }
        }
    }

    /**
     * @brief Get the global vector index of a particle based on chainN of a molecule and the pos position in this molecule
     * @param chainN - sequential number of all molecules bigger than 1 particle
     * @param pos - the sequential number of a particle within this molecule
     * @return index of this particle in the global vector
     */
    inline int getChainPart(int chainN, int pos) {
        for(int type=0; type<=molTypeCount; type++) {
            if(chainCount[type] > chainN) {
                type--;
                if(topo.moleculeParam[type].molSize() <= pos) return -1;
                chainN -= chainCount[type];
                return first[type]+chainN*topo.moleculeParam[type].molSize() + pos;
            }
        }
        assert(false && "ChainN not found");
        return -1;
    }

    inline ConList getConlist(int part1) {
        ConList conlist; // all NULL already

        if(topo.moleculeParam[(*this)[part1].molType].isAtomic())  // particle not in any chain
            return conlist;

        int pos = (part1 - first[(*this)[part1].molType]) % topo.moleculeParam[ (*this)[part1].molType ].molSize(); // position within chain

        if(pos > 0)
            conlist.conlist[0] = &this->operator [](part1-1);
        if(pos+1 < topo.moleculeParam[ (*this)[part1].molType ].molSize())
        conlist.conlist[1] = &this->operator [](part1+1);
        if(pos > 1)
            conlist.conlist[2] = &this->operator [](part1-2);
        if(pos+2 < topo.moleculeParam[ (*this)[part1].molType ].molSize())
            conlist.conlist[3] = &this->operator [](part1+2);

        return conlist;
    }

    inline ConList getConlist(int part1, int pos) {
        ConList conlist; // all NULL already

        if(pos > 0)
            conlist.conlist[0] = &this->operator [](part1-1);
        if(pos+1 < topo.moleculeParam[(*this)[part1].molType].molSize())
        conlist.conlist[1] = &this->operator [](part1+1);
        if(pos > 1)
            conlist.conlist[2] = &this->operator [](part1-2);
        if(pos+2 < topo.moleculeParam[(*this)[part1].molType].molSize())
            conlist.conlist[3] = &this->operator [](part1+2);

        return conlist;
    }


    Molecule getChain(int chainN) {
        Molecule ret;
        for(int type=0; type<=molTypeCount; type++) {
            if(chainCount[type] > chainN) {
                type--;
                chainN -= chainCount[type];
                chainN *= topo.moleculeParam[type].molSize();
                for(int i=0; i < topo.moleculeParam[type].molSize(); i++)
                    ret.push_back(first[type]+chainN + i);
                return ret;
            }
        }
        assert(false && "ChainN not found");
        return ret;
    }

    /**
     * @brief getMolOfPart - Get Molecule of Particle
     * @param index
     * @return Molecule
     */
    Molecule getMolOfPart(int index) {
        int i=0;
        while(index >= first[i]) {
            i++;
        }
        i--; // i is now the molecule type
        index -= first[i];//index sequnce number of particle in moltype
        index /= topo.moleculeParam[i].molSize(); // index sequnce number of molecule in type
        return getMolecule(index, i, topo.moleculeParam[i].molSize());
    }

    /**
     * @brief getMolecule
     * @param molN - sequence number of molecule in given moltype
     * @param molType - molecular type
     * @param molSize - number of particles in 1 molecule
     * @return
     */
    Molecule getMolecule(int molN, int molType, int molSize) {
        Molecule ret;
        molN *= molSize;
        molN += first[molType];

        for(int i=0; i < molSize; i++)
            ret.push_back(molN + i);

        assert(this->operator [](ret[0]).molType == molType); // did we pick correct molType
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
     * @param molType
     * @param molN - sequence number of molecule in given moltype
     * @return Index of first particle of molecule
     */
    inline int getStoreIndex(int molType, int molN) const {
        return first[molType] + topo.moleculeParam[molType].molSize()*molN;
    }

    inline int getInsertIndex(int molType) const {
        return first[molType+1];
    }

    /**
     * @param molType Type of molecule
     * @return Number of molecules a given type
     */
    int molCountOfType(int molType) {
        assert( (first[molType+1] - first[molType]) / topo.moleculeParam[molType].molSize() >= 0);
        return ( first[molType+1] - first[molType]) / topo.moleculeParam[molType].molSize();
    }



    long switchPartCount() {
        int size=0;
        for(int i=0; i<molTypeCount; i++) {
            // size += number of swichtypes per molType * number of molecules of moltype
            size += topo.moleculeParam[i].switchCount*( first[i+1] - first[i]) / topo.moleculeParam[i].molSize();
        }
        return size;
    }

    /**
     * @brief getSwitchPart
     * @param target - sequnce number of chosen particle with switchtype from all particles with switchtypes
     * @param seq - sequence number of type in chain
     * @return
     */
    int getSwitchPart(int target, int& seq) {
        int count;
        for(int i=0; i<molTypeCount; i++) {
            count = topo.moleculeParam[i].switchCount*( first[i+1] - first[i]) / topo.moleculeParam[i].molSize(); // number of particles with switchtype per molType
            if(target > count)
                target -= count;
            else {
                seq = target % topo.moleculeParam[i].switchCount; // sequence number of type in chain
                target /= topo.moleculeParam[i].switchCount;          // sequence number of molecule
                // target = first index of type (==offset) + target molecule index + particle sequence
                target = first[i] + target*topo.moleculeParam[i].molSize() + topo.moleculeParam[i].switchTypeSeq(seq);
                break;
            }
        }
        return target;
    }

    void insertMolecule(std::vector<Particle>& molecule) {
        insert(begin()+first[molecule[0].molType+1], molecule.begin(), molecule.end()); // insert actual data at the end of moltype (== begining of next type)

        for(int i= molecule[0].molType+1; i<=molTypeCount; i++)     // move next molTypes, molSize() == number of particles in molecule
            first[i] += topo.moleculeParam[molecule[0].molType].molSize();

        if(topo.moleculeParam[molecule[0].molType].molSize() > 1) // recalculate Particle vector, TODO do this nicer
            calcChainCount();
    }

    void removeMolecule(Molecule& mol) {
        for(int i=this->operator [](mol[0]).molType+1; i<=molTypeCount; i++) // move next molTypes, molSize() == number of particles in molecule
            first[i] -= topo.moleculeParam[ (*this)[mol[0]].molType ].molSize();

        // mol[0] - index of first particle of molecule
        // (*this)[ mol[0] ] - access Particle
        // (*this)[mol[0]].molType - moltype
        // topo.moleculeParam[molType].molSize() ...
        if(topo.moleculeParam[ (*this)[mol[0]].molType ].molSize() > 1) // recalculate Particle vector, TODO do this nicer
            calcChainCount();

        erase(begin()+mol[0], begin()+mol[0]+mol.size()); // delete actual data
    }


#ifndef NDEBUG
    bool checkConsistency() {
        for(int i=0; i<molTypeCount; i++) {

            if(first[i] > first[i+1]) // current entry of first cant be larger than the next one
                return false;
        }
        return true;
    }
#endif
};


/**
 * @brief Configuration of the system
 */
class Conf {
public:  
    ParticleVector pvec;  ///< \brief Main store of all particles, grouped by Molecular types
    std::vector<Neighbors > neighborList;
    //std::vector<ConList > conlist;
    ParticleVector pool; ///< \brief Store for chains for muVT insert of chain

    double sysvolume;                       ///< \brief Something like total mass -> for each += massOfParticle
    Vector syscm;                           ///< \brief System center of mass

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

    /**
     * @brief makeMoleculeWhole
     * Function move all particles in Molecule (*mol) to be directly next to first particle in molecule.
     * In other words function restore Molecule to be in one pice if molecule were broken by using PBC.
     * @param *mol
     */
    void makeMoleculeWhole( Molecule *mol ){ // BACHA POKUD JE mol->size() == 0 TAK SE SAHA NA NEPLATNOU PAMET !!!
        Vector r_cm;
        for ( std::vector<int>::iterator it = mol->begin()+1 ; it != mol->end() ; ++it ){
            r_cm = geo.image( &pvec[(*it)].pos, &pvec[(*(it-1))].pos );
            r_cm.x/=geo.box.x;//Image function return distance in real space so we need to get back in internal coordinates [0:1]
            r_cm.y/=geo.box.y;
            r_cm.z/=geo.box.z;
            pvec[(*it)].pos = pvec[(*(it-1))].pos + r_cm;
        }
    }

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
            fprintf (outfile, "%15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e %d %d\n",
                     geo.box.x * pvec[i].pos.x,
                     geo.box.y * pvec[i].pos.y,
                     geo.box.z * ((pvec[i].pos.z) - anInt(pvec[i].pos.z)),
                     pvec[i].dir.x, pvec[i].dir.y, pvec[i].dir.z,
                     pvec[i].patchdir[0].x, pvec[i].patchdir[0].y, pvec[i].patchdir[0].z,
                    pvec[i].switched,
                    pvec[i].molType);
        }
    #else
        for (unsigned int i=0; i < pvec.size(); i++) {
            fprintf (outfile, "%15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e %d %d\n",
                     geo.box.x * ((pvec[i].pos.x) - anInt(pvec[i].pos.x)),
                     geo.box.y * ((pvec[i].pos.y) - anInt(pvec[i].pos.y)),
                     geo.box.z * ((pvec[i].pos.z) - anInt(pvec[i].pos.z)),
                     pvec[i].dir.x, pvec[i].dir.y, pvec[i].dir.z,
                     pvec[i].patchdir[0].x, pvec[i].patchdir[0].y, pvec[i].patchdir[0].z,
                    pvec[i].switched,
                    pvec[i].molType);
        }
    #endif
#endif
    }

    /**
     * @brief Adds molecule to pvec, ensures Particle order, grouplist, conlist
     * @param molecule
     */
    void insertMolecule(std::vector<Particle> &molecule);

    /**
     * @brief removeMolecule
     * @param molType
     * @param molID
     */
    void removeMolecule(Molecule &target);

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

    void info() {
        cout << "conf:\n";
        for(unsigned int i=0; i< pvec.size(); i++)
            cout << pvec[i].molType <<"=" << pvec[i].type << ",";
        cout << endl;
    }
};


#endif // CONF_H
