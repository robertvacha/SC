#ifndef MOLECULEPARAMS_H
#define MOLECULEPARAMS_H

#include "macros.h"


/**
 * @brief Parameters for inner interaction in chains
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
    double mu;                  ///< \brief chemical potential, specific to each molecule type
    double lnThermalWavelengh;  ///< \brief ln(de Broglie thermal wavelenght), specific to each type
    bool muVTmove;              ///< \brief true - attempt muVT move for that type
    int particleTypes[MAXCHL];  ///< \brief 0..40 single particle type of each particle, -1: no particle

public:
    MoleculeParams() {
        bond1eq = -1;
        bond1c = -1;
        bond2eq = -1;
        bond2c = -1;
        bonddc = -1;
        angle1eq = -1;
        angle1c = -1;
        angle2eq = -1;
        angle2c = -1;
        mu = -1;
        muVTmove = false;
        lnThermalWavelengh = -1;
        for(int j=0; j<MAXCHL; j++)
            particleTypes[j] = -1;
    }

    /**
     * @return True if Molecule of only one particle
     */
    bool isAtomic() {return (particleTypes[1] == -1);}

    int molSize() {
        int size = 0;
        for(int i=0; i<MAXCHL; i++) if(particleTypes[i] != -1) size++;
        return size;
    }

};

#endif // MOLECULEPARAMS_H
