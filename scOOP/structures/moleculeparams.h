#ifndef MOLECULEPARAMS_H
#define MOLECULEPARAMS_H

#include "macros.h"
#include <cassert>
#include <array>

/**
 * @brief Parameters for interactions of molecules(mainly inner interaction in chains), each molType has a MoleculeParams instance
 */
class MoleculeParams
{
public:
    static const int maxChainParticles = 100;

    string name;        ///< \brief name of the molecule
    int molType = -1;

    double bond1eq  = -1.0;     ///< \brief Equilibrium distance of harmonic bond between nearest neighbours
    double bond1c   = -1.0;      ///< \brief Spring constant for harmonic bond between nearest neighbours
    double bond2eq  = -1.0;     ///< \brief Equilibrium distance of harmonic bond between second nearest neighbours
    double bond2c   = -1.0;      ///< \brief Spring constant for harmonic bond between second nearest neighbours
    double bonddeq  = -1.0;     ///< \brief Equilibrium distance of directional harmonic bond between the nearest neighbours
    double bonddc   = -1.0;      ///< \brief Spring constant for directional harmonic bond between the nearest neighbours
    double bondheq  = -1.0;     ///< \brief Equilibrium distance of hinge like harmonic bond between the nearest neighbours
    double bondhc   = -1.0;      ///< \brief Spring constant for hinge like harmonic bond between the nearest neighbours

    double angle1eq = -1.0;    ///< \brief Equilibrium angle between two spherocylinders -neerest neighbours
    double angle1c  = -1.0;     ///< \brief Spring constant angle between two spherocylinders -nearest neighbours
    double angle2eq = -1.0;    ///< \brief Equilibrium angle between two spherocylinder patches -nearest neighbours
    double angle2c  = -1.0;     ///< \brief Spring constant for angle between two spherocylinder patches -nearest neighbours

    double activity = -1.0;                ///< \brief activity, specific to each molecule type
    double chemPot  = -1.0;                 ///< \brief mu + 3*ln(A), specific to each type

    std::vector<int> particleTypes = {};    ///< \brief 0..40 particle type
    std::array<int, maxChainParticles> switchTypes = {};
    std::array<double, maxChainParticles> deltaMu = {};
    int switchCount = 0;        ///< \brief count of particles with defined switchtype

public:
    MoleculeParams() { particleTypes.reserve(MAXCHL); switchTypes.fill(-1); }

    string toString() {
        stringstream ss;

        ss << name << ":{" << endl;
        ss << ((bond1c == -1.0) ? string() : ( string("bond1: ")+to_string(bond1c)+string(" ")+to_string(bond1eq)+string("\n")) );
        ss << ((bond2c == -1.0) ? string() : ( string("bond2: ")+to_string(bond2c)+string(" ")+to_string(bond2eq)+string("\n")) );
        ss << ((bonddc == -1.0) ? string() : ( string("bondd: ")+to_string(bonddc)+string(" ")+to_string(bonddeq)+string("\n")) );
        ss << ((bondhc == -1.0) ? string() : ( string("bondh: ")+to_string(bondhc)+string(" ")+to_string(bondheq)+string("\n")) );

        ss << ((angle1c == -1.0) ? string() : ( string("bond1: ")+to_string(angle1c)+string(" ")+to_string(angle1eq))+string("\n") );
        ss << ((angle2c == -1.0) ? string() : ( string("bond1: ")+to_string(angle2c)+string(" ")+to_string(angle2eq))+string("\n") );

        ss << ((activity == -1.0) ? string() : ( string("activity: ")+to_string(activity))+string("\n"));
        for(int j=0; j < particleTypes.size(); ++j) {
            ss << "particles: " << particleTypes[j]
                 << " " << ( ( switchTypes[j] == -1 ) ? string() : ( std::to_string(switchTypes[j]) ) )
                 << " " << ( ( deltaMu[j] == 0.0 ) ? string() : (std::to_string(deltaMu[j]) ) )
                 << endl;
        }
        ss << "}";

        return ss.str();
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
                && (this->switchCount == o.switchCount) && (this->name == o.name);
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
