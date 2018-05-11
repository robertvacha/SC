#ifndef MOLECULEPARAMS_H
#define MOLECULEPARAMS_H

#include "macros.h"
#include <cassert>

/**
 * @brief Parameters for interactions of molecules(mainly inner interaction in chains), each molType has a MoleculeParams instance
 */
class MoleculeParams
{
public:
    char * name;        ///< \brief name of the molecule
    int molType;

    double bond1eq;     ///< \brief Equilibrium distance of harmonic bond between nearest neighbours
    double bond1c;      ///< \brief Spring constant for harmonic bond between nearest neighbours

    double bond2eq;     ///< \brief Equilibrium distance of harmonic bond between second nearest neighbours
    double bond2c;      ///< \brief Spring constant for harmonic bond between second nearest neighbours

    double bonddeq;     ///< \brief Equilibrium distance of directional harmonic bond between the nearest neighbours
    double bonddc;      ///< \brief Spring constant for directional harmonic bond between the nearest neighbours

    double bondheq;     ///< \brief Equilibrium distance of hinge like harmonic bond between the nearest neighbours
    double bondhc;      ///< \brief Spring constant for hinge like harmonic bond between the nearest neighbours

    double angle1eq;    ///< \brief Equilibrium angle between two spherocylinders -neerest neighbours
    double angle1c;     ///< \brief Spring constant angle between two spherocylinders -nearest neighbours

    double angle2eq;    ///< \brief Equilibrium angle between two spherocylinder patches -nearest neighbours
    double angle2c;     ///< \brief Spring constant for angle between two spherocylinder patches -nearest neighbours

    double activity;                ///< \brief activity, specific to each molecule type
    double chemPot;                 ///< \brief mu + 3*ln(A), specific to each type
    std::vector<int> particleTypes;    ///< \brief 0..40 particle type
    std::vector<int> switchTypes;
    std::vector<double> deltaMu;
    int switchCount;        ///< \brief count of particles with defined switchtype

public:
    MoleculeParams() : name(NULL), bond1eq(-1.0), bond1c(-1.0), bond2eq(-1.0), bond2c(-1.0), bonddeq(-1), bonddc(-1.0), bondheq(-1.0), bondhc(-1.0),
                       angle1eq(-1.0), angle1c(-1.0), angle2eq(-1.0), angle2c(-1.0), activity(-1.0), chemPot(-1.0) {
        particleTypes.reserve(MAXCHL);
    }

    /**
     * @brief operator ==
     * @param o
     * @return
     */
    bool operator==(const MoleculeParams& o) const {
        //
        // operator== on std::vector: Compares the contents of two containers.
        //   Checks if the contents of lhs and rhs are equal, that is, whether lhs.size() == rhs.size() and each element in lhs has equivalent element in rhs at the same position.
        //
        return (this->particleTypes == o.particleTypes) && (this->switchTypes == o.switchTypes)
                && (this->deltaMu == o.deltaMu) && (this->molType == o.molType)
                && (this->bond1eq == o.bond1eq) && (this->bond1c == o.bond1c)
                && (this->bond2eq == o.bond2eq) && (this->bond2c == o.bond2c)
                && (this->bonddeq == o.bonddeq) && (this->bonddc == o.bonddc)
                && (this->bondheq == o.bondheq) && (this->bondhc == o.bondhc)
                && (this->angle1eq == o.angle1eq) && (this->angle1c == o.angle1c)
                && (this->angle2eq == o.angle2eq) && (this->angle2c == o.angle2c)
                && (this->activity == o.activity) && (this->chemPot == o.chemPot)
                && (this->switchCount == o.switchCount);
    }

    /**
     * @return True if Molecule of only one particle
     */
    bool isAtomic() {return (particleTypes.size() == 1);}

    bool isGrandCanonical() {return chemPot!=-1.0;}

    /**
     * @brief molSize
     * @return Number of particles in molecule
     */
    inline int molSize() {
        return particleTypes.size();
    }

    inline int switchTypeSeq(int seq) {
        for(int i=0; i<(int)switchTypes.size(); i++) {
            if(switchTypes[i] != -1)
                seq--;
            if(seq == -1)
                return i;
        }
        assert(false && "This should never happen!!!");
        return -1;
    }

    void info() {
        cout << name << ", activity: " << activity << ", atoms:\n";
        for(unsigned int i=0; i<particleTypes.size(); i++)
            cout << particleTypes[i] << ", ";
        cout << endl;
    }
};

#endif // MOLECULEPARAMS_H
