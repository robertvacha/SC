#ifndef MOLECULEPARAMS_H
#define MOLECULEPARAMS_H

#include "macros.h"

/**
 * @brief Parameters for interactions of molecules(mainly inner interaction in chains), each molType has a MoleculeParams instance
 */
class MoleculeParams
{
public:
    char * name;        ///< \brief name of the molecule
    double bond1eq;     ///< \brief Equilibrium distance of harmonic bond between nearest neighbours
    double bond1c;      ///< \brief Spring constant for harmonic bond between nearest neighbours
    double bond2eq;     ///< \brief Equilibrium distance of harmonic bond between second nearest neighbours
    double bond2c;      ///< \brief Spring constant for harmonic bond between second nearest neighbours
    double bonddeq;     ///< \brief Equilibrium distance of directional harmonic bond between the nearest neighbours
    double bonddc;      ///< \brief Spring constant for directional harmonic bond between the nearest neighbours
    double angle1eq;    ///< \brief Equilibrium angle between two spherocylinders -neerest neighbours
    double angle1c;     ///< \brief Spring constant angle between two spherocylinders -nearest neighbours
    double angle2eq;    ///< \brief Equilibrium angle between two spherocylinder patches -nearest neighbours
    double angle2c;     ///< \brief Spring constant for angle between two spherocylinder patches -nearest neighbours

    // For muVT enseble
    double activity;                ///< \brief activity, specific to each molecule type
    double chemPot;                 ///< \brief mu + 3*ln(A), specific to each type
    std::vector<int> particleTypes; ///< \brief 0..40 particle type, -1: no particle, What particle types make a molecule

    int delAcc;
    int delRej;
    int insAcc;
    int insRej;

    unsigned long long int muVtAverageParticles;
    unsigned int muVtSteps;

public:
    MoleculeParams() : bond1eq(-1.0), bond1c(-1.0), bond2eq(-1.0), bond2c(-1.0), /*bonddeq(-1),*/ bonddc(-1.0),
                       angle1eq(-1.0), angle1c(-1.0), angle2eq(-1.0), angle2c(-1.0),
                       activity(-1.0), chemPot(-1.0),
                       delAcc(0), delRej(0), insAcc(0), insRej(0),
                       muVtAverageParticles(0), muVtSteps(0) {
        particleTypes.reserve(MAXCHL);
    }

    /**
     * @return True if Molecule of only one particle
     */
    bool isAtomic() {return (particleTypes.size() == 1);}

    bool isGrandCanonical() {return chemPot!=-1.0;}

    int molSize() {
        return particleTypes.size();
    }

};

#endif // MOLECULEPARAMS_H
