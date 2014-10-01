#ifndef CONF_H
#define CONF_H

#include "topo.h"
#include "structures.h"
#include "moleculeparams.h"
#include "particle.h"
#include <assert.h>
#include <map>

#include <iostream>

/**
 * @brief The GroupList class - simple grouplist
 */
class GroupList {
public:
    int first[MAXMT];   ///< \brief Index(of pvec) of first particle of molecule type, array over molecular types
    int molSize[MAXMT]; ///< \brief Number of particles per molecule of moleculeType, array over molecular types
    int molTypeCount;   ///< \brief Count of molecular types in use
    int vecSize;        ///< \brief Size of vector

    GroupList() {
        for(int i=0; i<MAXMT; i++) {
            first[i] = -1; molSize[i] = -1;
        }
    }

    /**
     * @brief Converts molID of molType to pvec Index
     * molID -> starts from 0 for each molType
     * @return Index of first particle of molecule
     */
    int getStoreIndex(int molType, int molID) {
        return first[molType] + molSize[molType]*molID;
    }

    int getInsertIndex(int molType) {
        if(first[molType+1] != -1)
            return first[molType+1];
        return vecSize;
    }

    /**
     * @param molType Type of molecule
     * @return Number of molecules a given type
     */
    int molCountOfType(int molType) {
        if(first[molType+1] != -1)
            return ( first[molType+1] -  first[molType]) /  molSize[molType];
        return (vecSize -  first[molType]) /  molSize[molType];
    }

    void insertMolecule(int molType) {
        for(int i= molType+1; i<molTypeCount; i++)
            first[i] += molSize[molType];
        vecSize += molSize[molType];
    }

    void deleteMolecule(int molType) {
        for(int i=molType+1; i<molTypeCount; i++)
            first[i] -= molSize[molType];
        vecSize -= molSize[molType];
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

    // chainlist - molecules of 1 particle on included
    long chainlist[MAXN][MAXCHL];       ///< \brief List of particles in chain
    long chainCount;                    ///< \brief Number of chains

    bool pairlist_update;

public:
    /**
     * @brief Conf Constructor, initializing variables
     */
    Conf() {
        try{
            pvec.reserve(MAXN);
            conlist.reserve(MAXN);
            pool.reserve(MAXN);
            neighborList.reserve(MAXN);
        } catch(std::bad_alloc& bad) {
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for pvec, pool, conlist, neighborlist, conf inicializer");
            exit(1);
        }

        for (int i=0;i<MAXN;i++) {
            for (int j = 0; j < MAXCHL; j++){
                chainlist[i][j] = -1;
            }
        }
        chainCount=0;
    }

    /**
     * @brief Adds molecule to pvec, ensures Particle order, changes chainlist, grouplist
     * @param molecule
     * @param topo
     */
    void addMolecule(std::vector<Particle>* molecule);

    /**
     * @brief removeMolecule
     * @param molType
     * @param molID
     */
    void removeMolecule(int target, int size);

    /**
     * @brief massCenter
     * @param topo
     */
    void massCenter(Topo* topo);

    /**
     * @brief partvecinit calculate vectors on particles for speedup
     */
    void partVecInit(Topo *topo);

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
