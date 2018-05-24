/** @file totalenergycalculator.h*/

#ifndef TOTALENERGYCALCULATOR_H
#define TOTALENERGYCALCULATOR_H

#include <iomanip>
#include <limits>

#include "pairenergycalculator.h"
#include "paire.h"
#include "externalenergycalculator.h"
#include "../structures/sim.h"

using namespace std;

template<typename TFloat>
class EnergyMatrixF {
private:
    std::vector< std::vector<TFloat> > eMatrix; /// \brief Diagonal matrix[i][j], i > j
    std::vector< std::vector<TFloat> > eMatrix2;
    Conf* conf;

public:
    typedef TFloat TEMFloat;
    typedef std::vector< std::vector<TFloat> > TEmVector;

    EnergyMatrixF(Conf* conf) : conf(conf) {
        energyMatrix = &eMatrix;
        energyMatrixTrial = &eMatrix2;

        try{
            changes.reserve(conf->pvec.size());

            energyMatrix->reserve(conf->pvec.size()+1024);
            energyMatrix->resize(conf->pvec.size());
            for(unsigned int i=0; i < conf->pvec.size(); i++) {
                energyMatrix->operator [](i).resize(i);
            }
            mcout.get() << "Alloc of EM " << conf->pvec.size() * sizeof(TFloat) * conf->pvec.size() / (2 * 1024 * 1024) << " MB" << endl;

            energyMatrixTrial->reserve(conf->pvec.size()+1024);
            energyMatrixTrial->resize(conf->pvec.size());
            for(unsigned int i=0; i < conf->pvec.size(); i++) {
                energyMatrixTrial->operator [](i).resize(i);
            }
            mcout.get() << "Alloc of trial EM " << conf->pvec.size() * sizeof(TFloat) * conf->pvec.size() / (2 * 1024 * 1024) << " MB" << endl;
        } catch(std::bad_alloc& bad) {
            cerr << "\nTOPOLOGY ERROR: Could not allocate memory for Energy matrix" << endl;
        }
    }

    vector<double> changes;

    std::vector< std::vector<TFloat> >* energyMatrix;
    std::vector< std::vector<TFloat> >* energyMatrixTrial;

    void swapEMatrices() {
        std::swap(energyMatrix, energyMatrixTrial);
    }

    void fixEMatrixSingle(bool pairlist_update, int target) {
        // FIX NEIGHBORLIST
        if(pairlist_update) {
            for(unsigned int i=0; i < changes.size(); i++) {
                if(target > conf->neighborList[target].neighborID[i])
                    (*energyMatrix)[target][conf->neighborList[target].neighborID[i]] = changes[i];
                else
                    energyMatrix->operator [](conf->neighborList[target].neighborID[i])[target] = changes[i];
            }
        } else {
            assert(changes.size() == conf->pvec.size());
            for(unsigned int i=0; i < conf->pvec.size(); i++) {
                if(target > i)
                    energyMatrix->operator [](target)[i] = changes[i];
                if(target < i)
                    energyMatrix->operator [](i)[target] = changes[i];
            }
        }
    }

    void fixEMatrixChain(bool pairlist_update, Molecule chain) {
        // FIX ENERGY MATRIX
        if(chain.empty())
            return;
        if(pairlist_update) {
            vector<double>::iterator it = changes.begin();
            for(unsigned int i=0; i < chain.size(); i++) {
                for(int j=0; j < conf->neighborList[chain[i]].neighborCount; j++) {

                    if ( chain.back() < conf->neighborList[chain[i]].neighborID[j]) {
                        (*energyMatrix)[  conf->neighborList[chain[i]].neighborID[j]  ] [  chain[i]  ] = *it;
                        it++;
                    }
                    if ( chain[0] > conf->neighborList[chain[i]].neighborID[j] ) {
                        (*energyMatrix)[  chain[i]  ][  conf->neighborList[chain[i]].neighborID[j]  ] = *it;
                        it++;
                    }
                }
            }
        } else {
            for(unsigned int i=0; i<chain.size(); i++) { // for all particles in molecule
                // for cycles => pair potential with all particles except those in molecule
                for (int j = 0; j < chain[0]; j++) { // pair potential with all particles from 0 to the first one in molecule
                    energyMatrix->operator [](chain[i])[j] = changes[i * conf->pvec.size() + j];
                }

                for (int j = chain[chain.size()-1] + 1; j < (long)conf->pvec.size(); j++) {
                    energyMatrix->operator [](j)[chain[i]] = changes[i * conf->pvec.size() + j];
                }

            }
        }

    }

    void resizeEMatrix() {
        unsigned int oldSize = energyMatrix->size();
        energyMatrix->resize(conf->pvec.size());
        for(unsigned int i=oldSize; i<conf->pvec.size(); ++i) {
            energyMatrix->operator [](i).resize(i);
        }

        energyMatrixTrial->resize(conf->pvec.size());
        for(unsigned int i=oldSize; i<conf->pvec.size(); ++i) {
            energyMatrixTrial->operator [](i).resize(i);
        }
    }
};


typedef EnergyMatrixF<double> EnergyMatrix;

class EMResize {};

template<typename pairEFce>
class TotalE {
public:
    pairEFce pairE;

    ExternalEnergyCalculator exterE;

    bool pairListUpdate;
    Conf* conf;

    TotalE(Sim * sim, Conf * conf) : pairE(pairEFce(&conf->geo)), exterE(ExternalEnergyCalculator(&conf->geo.box)),
                                     pairListUpdate(sim->pairlist_update), conf(conf) {}

    virtual void update() {}
    virtual void update(int target) {}
    virtual void update(Molecule target) {}
    virtual void update(EMResize a) {}

    virtual void initEM() {}

    /**
     * @brief Calculates energy between all pairs, generates energy matrix. Returns energy
     * @return
     */
    virtual double allToAllTrial() {return 0.0;}

    /**
     * @brief allToAll - Calculation with the use of energy matrix
     * @return
     */
    virtual double allToAll() {return 0.0;}

    /**
     * @brief Calculates energy between particle target and rest, generates energy matrix. Returns energy
     * @return
     */
    virtual double oneToAllTrial(int target) {return 0.0;}

    /**
     * @brief Calculates energy between particle target and rest using energy matrix Returns energy
     * @return
     */
    virtual double oneToAll(int target) {return 0.0;}

    /**
     * @brief mol2others - Calculates energy between Molecule mol and rest using energy matrix Returns energy
     * @param mol
     * @return
     */
    virtual double mol2others(Molecule &mol) {return 0.0;}

    /**
     * @brief mol2others - Calculates energy between Molecule mol and rest, generates energy matrix change, Returns energy
     * @param mol
     * @param changes
     * @return
     */
    virtual double mol2othersTrial(Molecule &mol) {return 0.0;}

    /**
     * @brief p2p
     * @param part1
     * @param part2
     * @return
     */
    double p2p(int part1, int part2) {
        ConList conlist = conf->pvec.getConlist(part1);
        return pairE(&conf->pvec[part1], &conf->pvec[part2], &conlist);
    }

    double p2p(Particle* part1, int part2) {
        ConList conlist = conf->pvec.getConlist(part2);
        return pairE(part1, &conf->pvec[part2], &conlist);
    }

    /**
     * @brief mol2others - used on vector of particles which arent inserted into main pvec
     * @param mol
     * @return
     */
    double mol2others(vector<Particle> &mol) {
        ConList conlist; // empty conlist
        double energy=0.0;
        unsigned int i=0;

        for(i=0; i<mol.size(); i++) {
            if (topo.exter.exist)
                energy += this->exterE.extere2(&mol[i]);

            for(unsigned int j=0; j < conf->pvec.size(); j++)
                energy += pairE(&mol[i], &conf->pvec[j], &conlist);
        }
        return energy;
    }

    /**
     * @brief Calculates inner chain energy
     */
    double chainInner(ParticleVector& chain) {
        double energy = 0.0;
        ConList conlist;

        for (unsigned int i=0; i<chain.size(); i++) {
            for (unsigned int j=i+1; j<chain.size(); j++) {
                conlist = chain.getConlist( i );
                energy += (pairE)(&chain[i], &chain[j],  &conlist);
            }
        }

        return energy;
    }

    double chainInner(Molecule& chain) {
        double energy = 0.0;
        ConList conlist;

        for (unsigned int i=0; i<chain.size(); i++) {
            for (unsigned int j=i+1; j<chain.size(); j++) {
                conlist = conf->pvec.getConlist(chain[i]);
                energy += (pairE)(&conf->pvec[chain[i]], &conf->pvec[chain[j]],  &conlist);
            }
        }
        return energy;
    }

    /**
     * @brief Calculates energy between particle "target" and the rest, conlist and neighbor list can be null, for GrandCanonical mainly
     */
    double oneToAll(Particle* target, ConList* conlist=ConList(), Neighbors* neighborList=NULL) {
        double energy=0.0;
        unsigned long i;

        if (pairListUpdate && neighborList!=NULL) {
            for (i = 0; i < (unsigned long)neighborList->neighborCount; i++){
                energy += (pairE)(target, &conf->pvec[ neighborList->neighborID[i] ], conlist );
            }
        } else {
            for (i = 0; i < conf->pvec.size(); i++) {
                if(target != &conf->pvec[i]) {
                    energy += pairE(target, &conf->pvec[i], conlist);
                }
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(target);

        return energy;
    }
};

template<typename pairEFce>
class TotalEMatrix : public TotalE<pairEFce>
{
public:
    EnergyMatrix eMat;

    using TotalE<pairEFce>::mol2others;
    using TotalE<pairEFce>::conf;
    using TotalE<pairEFce>::pairE;
    using TotalE<pairEFce>::pairListUpdate;

    TotalEMatrix(Sim * sim, Conf * conf) : TotalE<pairEFce>(sim, conf), eMat(EnergyMatrix(conf)) {
        mcout.get() << "Simulation uses Energy Matrix acceleration" << endl;
    }

    void initEM() override {
        allToAll(eMat.energyMatrix);
    }

    void update() override {
        eMat.swapEMatrices();
    }

    void update(EMResize a) override {
        eMat.resizeEMatrix();
    }

    void update(int target) override {
        eMat.fixEMatrixSingle(pairListUpdate, target);
    }

    void update(Molecule target) {
        eMat.fixEMatrixChain(pairListUpdate, target);
    }

    double allToAllTrial() override {
        return allToAll(eMat.energyMatrixTrial);
    }

    double allToAll() override {
        if(conf->pvec.empty()) return 0.0;

        double energy = 0.0;
        for (unsigned int i = 1; i < conf->pvec.size(); ++i) {
            for (unsigned long j = 0; j < i; ++j) {
                energy += eMat.energyMatrix->operator [](i)[j];
            }
        }

        if (topo.exter.exist) //for every particle add interaction with external potential
            for (unsigned int i = 0; i < conf->pvec.size(); ++i)
                energy += this->exterE.extere2(&conf->pvec[i]);

        return energy;
    }

    double oneToAll(int target) override {
        assert(target >=0 && target < conf->pvec.size());
        double energy=0.0;

        if (pairListUpdate) {
            for (long i = 0; i < conf->neighborList[target].neighborCount; i++) {
                if(target > conf->neighborList[target].neighborID[i])
                    energy += eMat.energyMatrix->operator [](target)[conf->neighborList[target].neighborID[i]];
                else
                    energy += eMat.energyMatrix->operator [](conf->neighborList[target].neighborID[i])[target];
            }
        } else {
            for (long i = 0; i < (long)conf->pvec.size(); i++) {
                if(target > i) {
                    energy += eMat.energyMatrix->operator [](target)[i];
                }
                if(target < i) {
                    energy += eMat.energyMatrix->operator [](i)[target];
                }
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(&conf->pvec[target]);

        return energy;
    }

    double oneToAllTrial(int target) override {
        assert(target >= 0 && target < conf->pvec.size());

        double energy=0.0;
        long i;
        ConList conlist = conf->pvec.getConlist(target);

        if (pairListUpdate) {
            eMat.changes.resize( conf->neighborList[target].neighborCount );
            for (i = 0; i < conf->neighborList[target].neighborCount; i++) {
                eMat.changes[i] = pairE(&conf->pvec[target], &conf->pvec[ conf->neighborList[target].neighborID[i] ], &conlist);
                energy += eMat.changes[i];
            }
        } else {
            eMat.changes.resize(conf->pvec.size());
            eMat.changes[target] = 0.0;
            for (i = 0; i < target; i++) {
                eMat.changes[i] = pairE(&conf->pvec[target], &conf->pvec[i], &conlist);
                energy += eMat.changes[i];
            }
            for (i = target+1; i < (long)conf->pvec.size(); i++) {
                eMat.changes[i] = pairE(&conf->pvec[target], &conf->pvec[i], &conlist);
                energy += eMat.changes[i];
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(&conf->pvec[target]);

        return energy;

    }

    double mol2others(Molecule &mol) override {
        double energy=0.0;
        long i = 0;

        if (pairListUpdate) {

            for(unsigned int j=0; j<mol.size(); ++j) { // for all particles in molecule
                for (i = 0; i < conf->neighborList[mol[j]].neighborCount; ++i) {
                    // exclude interactions within molecule atoms
                    if( mol[0] > conf->neighborList[mol[j]].neighborID[i] || mol.back() < conf->neighborList[mol[j]].neighborID[i] ) {

                        if(mol[j] > conf->neighborList[mol[j]].neighborID[i])
                            energy += (*eMat.energyMatrix)[ mol[j] ][ conf->neighborList[mol[j]].neighborID[i] ];
                        else
                            energy += (*eMat.energyMatrix)[ conf->neighborList[mol[j]].neighborID[i] ][ mol[j] ];
                    }
                }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }
        } else {

            for(unsigned int j=0; j<mol.size(); ++j) { // for all particles in molecule

                for (i = 0; i < mol[0]; i++)  // pair potential with all particles from 0 to the first one in molecule
                    energy += (*eMat.energyMatrix)[ mol[j] ][ i ];

                for (long i = mol.back() + 1; i < (long)conf->pvec.size(); i++)
                    energy += (*eMat.energyMatrix)[ i ][ mol[j] ];

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }
        }
        return energy;
    }

    double mol2othersTrial(Molecule &mol) {
        ConList conlist;
        eMat.changes.clear();
        assert(eMat.changes.empty());

        double energy=0.0;
        long i = 0;

        if (pairListUpdate) {

            //
            // SAVE the changes to energy to easily swap them later -> instance when oneToAll(target, false, some vector) used
            //
            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < conf->neighborList[mol[j]].neighborCount; i++) {
                    if(mol[0] > conf->neighborList[mol[j]].neighborID[i] || mol.back() < conf->neighborList[mol[j]].neighborID[i]) {
                        eMat.changes.push_back( pairE(&conf->pvec[mol[j]], &conf->pvec[conf->neighborList[mol[j]].neighborID[i]], &conlist) );
                        energy += eMat.changes.back();
                    }
                }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }

        } else {

            eMat.changes.resize( mol.size() * conf->pvec.size() );
            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < mol[0]; i++) { // pair potential with all particles from 0 to the first one in molecule
                    eMat.changes[j*conf->pvec.size() + i] = pairE(&conf->pvec[mol[j]], &conf->pvec[i], &conlist);
                    energy += eMat.changes[j*conf->pvec.size() + i];
                }

                for (long i = mol[mol.size()-1] + 1; i < (long)conf->pvec.size(); i++) {
                    eMat.changes[j*conf->pvec.size() + i] = pairE(&conf->pvec[mol[j]], &conf->pvec[i], &conlist);
                    energy += eMat.changes[j*conf->pvec.size() + i];
                }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }
        }
        return energy;
    }

private:
    double allToAll(EnergyMatrix::TEmVector* energyMatrix) {
        if(conf->pvec.empty()) return 0.0;
        double energy=0.0;
        ConList conlist;
        assert(energyMatrix != NULL);

        for (unsigned int i = 1; i < conf->pvec.size(); ++i) {
            conlist = conf->pvec.getConlist(i);
            for (unsigned long j = 0; j < i; ++j) {
                (*energyMatrix)[i][j] = (pairE)(&conf->pvec[i], &conf->pvec[j], &conlist);
                energy += (*energyMatrix)[i][j];
            }
        }

        if (topo.exter.exist) //for every particle add interaction with external potential
            for (unsigned int i = 0; i < conf->pvec.size(); ++i)
                energy += this->exterE.extere2(&conf->pvec[i]);

        return energy;
    }
};

// PairList not included - need to be done for benchmarks
template<typename pairEFce>
class TotalEFull : public TotalE<pairEFce>
{
public:
    using TotalE<pairEFce>::mol2others;
    using TotalE<pairEFce>::conf;
    using TotalE<pairEFce>::pairE;
    using TotalE<pairEFce>::pairListUpdate;

    TotalEFull(Sim * sim, Conf * conf) : TotalE<pairEFce>(sim, conf) {}

    double allToAllTrial() override {
        return allToAll();
    }

    double allToAll() override {
        if(conf->pvec.empty()) return 0.0;

        double energy=0.0;
        ConList conlist;
        for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
            conlist = conf->pvec.getConlist(i);
            for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
                energy += (pairE)(&conf->pvec[i], &conf->pvec[j], &conlist);
            }
        }

        if (topo.exter.exist && !conf->pvec.empty()) //add interaction with external potential
            for (unsigned int i = 0; i < conf->pvec.size(); i++)
                energy += this->exterE.extere2(&conf->pvec[i]);

        return energy;
    }

    double oneToAllTrial(int target) override {
        return oneToAll(target);
    }

    double oneToAll(int target) override {
        ConList conlist = conf->pvec.getConlist(target);

        double energy=0.0;

        if (pairListUpdate) {
            for (int i = 0; i < conf->neighborList[target].neighborCount; i++) {
                energy += pairE(&conf->pvec[target], &conf->pvec[ conf->neighborList[target].neighborID[i] ], &conlist);
            }
        } else {
            for (int i = 0; i < (int)conf->pvec.size(); i++) {
                if(target != i)
                    energy += (pairE)(&conf->pvec[target], &conf->pvec[i], &conlist);
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(&conf->pvec[target]);

        return energy;
    }

    double mol2othersTrial(Molecule &mol) override {
        return mol2others(mol);
    }

    double mol2others(Molecule &mol) override {
        ConList conlist; // Empty conlist, looping over particles that arent in the molecule
        double energy=0.0;
        long i = 0;

        if (pairListUpdate) {
            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < conf->neighborList[mol[j]].neighborCount; i++)
                    if(mol[0] > conf->neighborList[mol[j]].neighborID[i] || mol.back() < conf->neighborList[mol[j]].neighborID[i]) {
                        energy += pairE(&conf->pvec[mol[j]], &conf->pvec[conf->neighborList[mol[j]].neighborID[i]], &conlist);
                    }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }
        } else {
            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < mol[0]; i++) // pair potential with all particles from 0 to the first one in molecule
                    energy += pairE(&conf->pvec[mol[j]], &conf->pvec[i], &conlist);

                for (i = mol[mol.size()-1] + 1; i < (long)conf->pvec.size(); i++)
                    energy += pairE(&conf->pvec[mol[j]], &conf->pvec[i], &conlist);

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }
        }
        return energy;
    }
};

template<typename pairEFce>
class TotalEFullSymetry : public TotalE<pairEFce>
{
public:
    using TotalE<pairEFce>::mol2others;
    using TotalE<pairEFce>::conf;
    using TotalE<pairEFce>::pairE;
    using TotalE<pairEFce>::pairListUpdate;
    using TotalE<pairEFce>::eMat;

    TotalEFullSymetry(Sim * sim, Conf * conf) : TotalE<pairEFce>(sim, conf) {}

    bool test(double a, double b) {
        return fabs(a - b) > 1e-10;
    }

    double allToAll(EnergyMatrix::TEmVector* energyMatrix) override {
        return allToAll();
    }

    double allToAll() override {
        if(conf->pvec.empty()) return 0.0;

        double energy=0.0;
        ConList conlist;
        for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
            conlist = conf->pvec.getConlist(i);
            for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
                energy += (pairE)(&conf->pvec[i], &conf->pvec[j], &conlist);
            }

            if (topo.exter.exist) //for every particle add interaction with external potential
                energy += this->exterE.extere2(&conf->pvec[i]);
        }

        if (topo.exter.exist && !conf->pvec.empty()) //add interaction of last particle with external potential
            energy += this->exterE.extere2(&conf->pvec.back());

        return energy;
    }

    double oneToAll(int target, vector<double> *changes) override {
        return oneToAll(target);
    }

    double oneToAll(int target) override {
        ConList conlist = conf->pvec.getConlist(target);
        ConList conlist2;

        double energy=0.0, temp, tempSym;

        for (int i = 0; i < (int)conf->pvec.size(); i++) {
            if(target != i) {
                temp = (pairE)(&conf->pvec[target], &conf->pvec[i], &conlist);
                conlist2 = conf->pvec.getConlist(i);
                tempSym = (pairE)(&conf->pvec[i], &conf->pvec[target], &conlist2);
                if( test(temp, tempSym) )
                    cout << std::setprecision(10) << "Error, energy: " << temp << " symetrical: " << tempSym << endl;
                energy += temp;
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(&conf->pvec[target]);

        return energy;
    }

    double mol2others(Molecule &mol, vector<double> *changes) override {
        return mol2others(mol);
    }

    double mol2others(Molecule &mol) override {
        double energy=0.0;
        long i = 0;

        for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
            for (i = 0; i < mol[0]; i++) // pair potential with all particles from 0 to the first one in molecule
                energy += pairE(&conf->pvec[mol[j]], &conf->pvec[i]);

            for (i = mol[mol.size()-1] + 1; i < (long)conf->pvec.size(); i++)
                energy += pairE(&conf->pvec[mol[j]], &conf->pvec[i]);

            if (topo.exter.exist) //add interaction with external potential
                energy += this->exterE.extere2(&conf->pvec[mol[j]]);
        }
        return energy;
    }
};

template<typename TotE, typename TotEControl>
class TestE : public TotE {
public:
    using TotE::mol2others;

    TotEControl control;

    TestE(Sim * sim, Conf * conf) : TotE(sim, conf), control(sim, conf)  {}

    bool test(double a, double b) {
        return fabs(a - b) > 1e-7;
    }

    double allToAllTrial() override {
        double e = 0.0, eControl = 0.0;

        e = TotE::allToAllTrial();
        eControl = control.allToAllTrial();

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error(e), allToAll energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double allToAll() override {
        double e = 0.0, eControl = 0.0;

        e = TotE::allToAll();
        eControl = control.allToAll();

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error(), allToAll energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double oneToAllTrial(int target) override {
        double e = 0.0, eControl = 0.0;

        e = TotE::oneToAllTrial(target);
        eControl = control.oneToAllTrial(target);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error(c), oneToAll(change) energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double oneToAll(int target) override {
        double e = 0.0, eControl = 0.0;

        e = TotE::oneToAll(target);
        eControl = control.oneToAll(target);

        if(test(e,eControl)) {
            cout << std::setprecision(10) << "Error, oneToAll() energy: " << e << " control: " << eControl << ", itself calculated: " << TotE::oneToAllTrial(target) << endl;
        }

        return e;
    }

    double mol2others(Molecule &mol, vector<double> *changes) {
        double e = 0.0, eControl = 0.0;

        e = TotE::mol2others(mol, changes);
        eControl = control.mol2others(mol, changes);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, mol2others energy: " << e << " control: " << eControl << endl;

        return e;
    }


    double mol2othersTrial(Molecule &mol) override {
        double e = 0.0, eControl = 0.0;

        e = TotE::mol2othersTrial(mol);
        eControl = control.mol2othersTrial(mol);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, mol2others energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double mol2others(Molecule &mol) override {
        double e = 0.0, eControl = 0.0;

        e = TotE::mol2others(mol);
        eControl = control.mol2others(mol);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, mol2others energy: " << e << " control: " << eControl << endl;

        return e;
    }
};

template<typename pairEFce>
class TotalEHardSphereFull : public TotalE<pairEFce>
{
public:
    using TotalE<pairEFce>::mol2others;
    using TotalE<pairEFce>::conf;
    using TotalE<pairEFce>::pairE;
    using TotalE<pairEFce>::pairListUpdate;
    using TotalE<pairEFce>::eMat;

    TotalEHardSphereFull(Sim * sim, Conf * conf) : TotalE<pairEFce>(sim, conf) {}

    inline double allToAll(EnergyMatrix::TEmVector* energyMatrix) override {
        return allToAll();
    }

    inline double allToAll() override { // just to check a valid config
        if(conf->pvec.empty()) return 0.0;

        Vector r_cm;
        double dotrcm;
        for(unsigned int i=0; i < conf->pvec.size()-1; ++i) {
            for(unsigned int j= i+1; j < conf->pvec.size(); ++j) {
                r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[j].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].sigmaSq)
                    return std::numeric_limits<double>::infinity();
            }
        }

        return 0.0;
    }

    inline double oneToAll(int target, vector<double> *changes) override {
        Vector r_cm;
        double dotrcm;

        if (pairListUpdate) {
            for (int i = 0; i < conf->neighborList[target].neighborCount; i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[ conf->neighborList[target].neighborID[i] ].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].sigmaSq)
                    return std::numeric_limits<double>::infinity();
            }
        } else {
            for(unsigned int i=0; i < conf->pvec.size(); ++i) {
                if(i != target) {
                    r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[target].pos);
                    dotrcm = r_cm.dot(r_cm);
                    if(dotrcm < topo.ia_params[conf->pvec[i].type][conf->pvec[target].type].sigmaSq)
                        return std::numeric_limits<double>::infinity();
                }
            }
        }

        return 0.0;
    }

    inline double oneToAll(int target) override { // before move, always 0.0 in a valid config
        return 0.0;
    }

    inline double mol2others(Molecule &mol, vector<double> *changes) override {
        return 0.0;
    }

    inline double mol2others(Molecule &mol) override {
        return 0.0;
    }
};

template<typename pairEFce>
class TotalEHardSphereEMatrix : public TotalE<pairEFce>
{
public:
    using TotalE<pairEFce>::mol2others;
    using TotalE<pairEFce>::conf;
    using TotalE<pairEFce>::pairE;
    using TotalE<pairEFce>::pairListUpdate;
    using TotalE<pairEFce>::eMat;

    TotalEHardSphereEMatrix(Sim * sim, Conf * conf) : TotalE<pairEFce>(sim, conf) {}

    double allToAll(EnergyMatrix::TEmVector* energyMatrix) override {
        if(conf->pvec.empty()) return 0.0;
        Vector r_cm;
        double dotrcm;
        double energy=0.0;
        assert(energyMatrix != NULL);

        for (unsigned int i = 1; i < conf->pvec.size(); ++i) {
            for (unsigned long j = 0; j < i; ++j) {
                r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[j].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].sigmaSq) {
                    (*energyMatrix)[i][j] = std::numeric_limits<double>::infinity();
                } else {
                    (*energyMatrix)[i][j] = 0.0;
                }
                energy += (*energyMatrix)[i][j];
            }
        }
        return energy;
    }

    double allToAll() override {
        if(conf->pvec.empty()) return 0.0;
        double energy = 0.0;

        for (unsigned int i = 1; i < conf->pvec.size(); ++i) {
            for (unsigned long j = 0; j < i; ++j) {
                energy += eMat.energyMatrix->operator [](i)[j];
            }
        }
        return energy;
    }

    double oneToAll(int target, vector<double> *changes) override {
        assert(target >= 0 && target < conf->pvec.size());
        assert(changes != NULL);

        Vector r_cm;
        double dotrcm;
        double energy=0.0;
        long i;

        if (pairListUpdate) {
            changes->resize( conf->neighborList[target].neighborCount );
            for (i = 0; i < conf->neighborList[target].neighborCount; i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[conf->neighborList[target].neighborID[i]].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[conf->neighborList[target].neighborID[i]].type].sigmaSq) {
                    (*changes)[i] = std::numeric_limits<double>::infinity();
                    energy += std::numeric_limits<double>::infinity();
                } else
                    (*changes)[i] = 0.0;
            }
        } else {
            changes->resize(conf->pvec.size());
            (*changes)[target] = 0.0;
            for (i = 0; i < target; i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[i].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].sigmaSq) {
                    (*changes)[i] = std::numeric_limits<double>::infinity();
                    energy = std::numeric_limits<double>::infinity();
                }
                else
                    (*changes)[i] = 0.0;
            }
            for (i = target+1; i < (long)conf->pvec.size(); i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[i].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].sigmaSq) {
                    (*changes)[i] = std::numeric_limits<double>::infinity();
                    energy = std::numeric_limits<double>::infinity();
                }
                else
                    (*changes)[i] = 0.0;
            }
        }

        return energy;
    }

    double oneToAll(int target) override {
        assert(target >=0 && target < conf->pvec.size());

        double energy = 0.0;
        if (pairListUpdate) {
            for (long i = 0; i < conf->neighborList[target].neighborCount; i++) {
                if(target > conf->neighborList[target].neighborID[i])
                    energy += eMat.energyMatrix->operator [](target)[conf->neighborList[target].neighborID[i]];
                else
                    energy += eMat.energyMatrix->operator [](conf->neighborList[target].neighborID[i])[target];
            }
        } else {
            for(unsigned int i=0; i < conf->pvec.size(); ++i) {
                if(i != target) {
                    if( eMat.energyMatrix->operator [](target)[i] == std::numeric_limits<double>::infinity() )
                        return std::numeric_limits<double>::infinity();
                }
            }
        }

        return 0.0;
    }

    double mol2others(Molecule &mol) override {
        return 0.0;
    }

    double mol2others(Molecule &mol, vector<double> *changes) override {
        return 0.0;
    }
};

template<typename pairEFce>
class TotalELJFull : public TotalE<pairEFce>
{
public:
    using TotalE<pairEFce>::mol2others;
    using TotalE<pairEFce>::conf;
    using TotalE<pairEFce>::pairE;
    using TotalE<pairEFce>::pairListUpdate;
    using TotalE<pairEFce>::eMat;

    TotalELJFull(Sim * sim, Conf * conf) : TotalE<pairEFce>(sim, conf) {}

    inline double allToAll(EnergyMatrix::TEmVector* energyMatrix) override {
        return allToAll();
    }

    inline double allToAll() override { // just to check a valid config
        if(conf->pvec.empty()) return 0.0;

        Vector r_cm;
        double dotrcm, energy = 0.0;
        for(unsigned int i=0; i < conf->pvec.size()-1; ++i) {
            for(unsigned int j= i+1; j < conf->pvec.size(); ++j) {
                r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[j].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].sigmaSq*(6.25))
                    energy += topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].A * pow(dotrcm, -6)
                            - topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].B * pow(dotrcm, -3);
            }
        }

        return energy;
    }

    inline double oneToAll(int target, vector<double> *changes) override {
        Vector r_cm;
        double dotrcm, energy = 0.0;

        if (pairListUpdate) {
            for (int i = 0; i < conf->neighborList[target].neighborCount; i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[ conf->neighborList[target].neighborID[i] ].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].sigmaSq*(6.25))
                    energy += topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].A * pow(dotrcm, -6)
                            - topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].B * pow(dotrcm, -3);
            }
        } else {
            for(unsigned int i=0; i < conf->pvec.size(); ++i) {
                if(i != target) {
                    r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[target].pos);
                    dotrcm = r_cm.dot(r_cm);
                    if(dotrcm < topo.ia_params[conf->pvec[i].type][conf->pvec[target].type].sigmaSq*(6.25))
                        energy += topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].A * pow(dotrcm, -6)
                                - topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].B * pow(dotrcm, -3);
                }
            }
        }

        return energy;
    }

    inline double oneToAll(int target) override {
        Vector r_cm;
        double dotrcm, energy = 0.0;

        if (pairListUpdate) {
            for (int i = 0; i < conf->neighborList[target].neighborCount; i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[ conf->neighborList[target].neighborID[i] ].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].sigmaSq*(6.25))
                    energy += topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].A * pow(dotrcm, -6)
                            - topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].B * pow(dotrcm, -3);
            }
        } else {
            for(unsigned int i=0; i < conf->pvec.size(); ++i) {
                if(i != target) {
                    r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[target].pos);
                    dotrcm = r_cm.dot(r_cm);
                    if(dotrcm < topo.ia_params[conf->pvec[i].type][conf->pvec[target].type].sigmaSq*(6.25))
                        energy += topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].A * pow(dotrcm, -6)
                                - topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].B * pow(dotrcm, -3);
                }
            }
        }

        return energy;
    }

    inline double mol2others(Molecule &mol, vector<double> *changes) override {
        return 0.0;
    }

    inline double mol2others(Molecule &mol) override {
        return 0.0;
    }
};

template<typename pairEFce>
class TotalELJEMatrix : public TotalE<pairEFce>
{
public:
    using TotalE<pairEFce>::mol2others;
    using TotalE<pairEFce>::conf;
    using TotalE<pairEFce>::pairE;
    using TotalE<pairEFce>::pairListUpdate;
    using TotalE<pairEFce>::eMat;

    TotalELJEMatrix(Sim * sim, Conf * conf) : TotalE<pairEFce>(sim, conf) {}

    double allToAll(EnergyMatrix::TEmVector* energyMatrix) override {
        if(conf->pvec.empty()) return 0.0;
        Vector r_cm;
        double dotrcm;
        double energy=0.0;
        assert(energyMatrix != NULL);

        for (unsigned int i = 1; i < conf->pvec.size(); ++i) {
            for (unsigned long j = 0; j < i; ++j) {
                r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[j].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].sigmaSq*(6.25)) // 2.5 sigma cutoff
                    (*energyMatrix)[i][j] = topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].A * pow(dotrcm, -6)
                        - topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].B * pow(dotrcm, -3);
                else
                    (*energyMatrix)[i][j] = 0.0;
                energy += (*energyMatrix)[i][j];
            }
        }
        return energy;
    }

    double allToAll() override {
        if(conf->pvec.empty()) return 0.0;
        double energy = 0.0;

        for (unsigned int i = 1; i < conf->pvec.size(); ++i) {
            for (unsigned long j = 0; j < i; ++j) {
                energy += this->eMat.energyMatrix->operator [](i)[j];
            }
        }
        return energy;
    }

    double oneToAll(int target, vector<double> *changes) override {
        assert(target >= 0 && target < conf->pvec.size());
        assert(changes != NULL);

        Vector r_cm;
        double dotrcm;
        double energy=0.0;

        if (pairListUpdate) {
            changes->resize( conf->neighborList[target].neighborCount );
            for (int i = 0; i < conf->neighborList[target].neighborCount; i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[conf->neighborList[target].neighborID[i]].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[conf->neighborList[target].neighborID[i]].type].sigmaSq*(6.25)) {
                    (*changes)[i] = topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].A * pow(dotrcm, -6)
                                      - topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].B * pow(dotrcm, -3);
                    energy += (*changes)[i];
                } else
                    (*changes)[i] = 0.0;
            }
        } else {
            changes->resize(conf->pvec.size());
            (*changes)[target] = 0.0;
            for (int i = 0; i < target; i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[i].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].sigmaSq*(6.25)) {
                    (*changes)[i] = topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].A * pow(dotrcm, -6)
                                  - topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].B * pow(dotrcm, -3);
                    energy += (*changes)[i];
                }
                else
                    (*changes)[i] = 0.0;
            }
            for (int i = target+1; i < (long)conf->pvec.size(); i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[i].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].sigmaSq*(6.25)) {
                    (*changes)[i] = topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].A * pow(dotrcm, -6)
                                  - topo.ia_params[conf->pvec[target].type][conf->pvec[i].type].B * pow(dotrcm, -3);
                    energy += (*changes)[i];
                }
                else
                    (*changes)[i] = 0.0;
            }
        }

        return energy;
    }

    double oneToAll(int target) override {
        assert(target >=0 && target < conf->pvec.size());
        double energy = 0.0;

        if (pairListUpdate) {
            for (long i = 0; i < conf->neighborList[target].neighborCount; i++) {
                if(target > conf->neighborList[target].neighborID[i])
                    energy += (*eMat.energyMatrix)[ target ][ conf->neighborList[target].neighborID[i] ];
                else
                    energy += (*eMat.energyMatrix)[ conf->neighborList[target].neighborID[i] ][ target ];
            }
        } else {
            for(unsigned int i=0; i < conf->pvec.size(); ++i) {
                if(i != target) {
                    energy += (*eMat.energyMatrix)[target][i];
                }
            }
        }

        return energy;
    }

    double mol2others(Molecule &mol) override {
        return 0.0;
    }

    double mol2others(Molecule &mol, vector<double> *changes) override {
        return 0.0;
    }
};

//
//  Full-scope Total energy calculators
//
//typedef TotalE<PairE> TotalEnergyCalculator;                        // Empty Energy
//typedef TotalEFull<PairEnergyCalculator> TotalEnergyCalculator;         // Energy matrix optimization
typedef TotalEMatrix<PairE> TotalEnergyCalculator;                        // Energy matrix optimization
//typedef TotalEFull<PairE> TotalEnergyCalculator;                          // Full calculation
//typedef TestE< TotalEMatrix<PairE>, TotalEFull<PairE> > TotalEnergyCalculator;         // Test of Pair energy
//typedef TotalEFullSymetry<PairE> TotalEnergyCalculator;         // Test of Pair energy symetry

//
//  Specialized Total energy calculators
//
//typedef TotalEHardSphereFull<PairE> TotalEnergyCalculator;         // Hardspheres, full calculation
//typedef TotalEHardSphereEMatrix<PairE> TotalEnergyCalculator;         // Hardspheres, energy matrix calculation
//typedef TotalELJFull<PairE> TotalEnergyCalculator;         // LJ, full calculation //cut-off distance of rc = 2.5σ,
//typedef TotalELJEMatrix<PairE> TotalEnergyCalculator;         // LJ, energy matrix calculation //cut-off distance of rc = 2.5σ,



#endif // TOTALENERGYCALCULATOR_H
