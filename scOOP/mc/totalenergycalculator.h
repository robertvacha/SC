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

class EnergyMatrix {
private:
    std::vector< std::vector<double> > eMatrix;
    std::vector< std::vector<double> > eMatrix2;
    Conf* conf;

public:
    EnergyMatrix(Conf* conf) : conf(conf) {
        energyMatrix = &eMatrix;
        energyMatrixTrial = &eMatrix2;

        try{
            changes.reserve(1024 + conf->pvec.size());
            energyMatrix->reserve(1024 + conf->pvec.size());
            energyMatrix->resize(conf->pvec.size());
            for(unsigned int i=0; i < conf->pvec.size(); i++) {
                energyMatrix->operator [](i).reserve(1024 + conf->pvec.size());
                energyMatrix->operator [](i).resize(conf->pvec.size());
            }

            energyMatrixTrial->reserve(1024 + conf->pvec.size());
            energyMatrixTrial->resize(conf->pvec.size());
            for(unsigned int i=0; i < conf->pvec.size(); i++) {
                energyMatrixTrial->operator [](i).reserve(1024 + conf->pvec.size());
                energyMatrixTrial->operator [](i).resize(conf->pvec.size());
            }
        } catch(std::bad_alloc& bad) {
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for Energy matrix");
            exit(1);
        }
    }

    vector<double> changes;

    std::vector< std::vector<double> >* energyMatrix;
    std::vector< std::vector<double> >* energyMatrixTrial;

    void swapEMatrices() {
        std::vector< std::vector<double> >*  temp;
        temp = energyMatrix;
        energyMatrix = energyMatrixTrial;
        energyMatrixTrial = temp;
    }

    void fixEMatrixSingle(bool pairlist_update, int target) {

        // FIX NEIGHBORLIST
        if(pairlist_update) {
            for(unsigned int i=0; i < changes.size(); i++) {
                (*energyMatrix)[target][conf->neighborList[target].neighborID[i]] = changes[i];
                energyMatrix->operator [](conf->neighborList[target].neighborID[i])[target] = changes[i];
            }
        } else {
            assert(changes.size() == pvec.size());
            for(unsigned int i=0; i < conf->pvec.size(); i++) {
                energyMatrix->operator [](target)[i] = changes[i];
                energyMatrix->operator [](i)[target] = changes[i];
            }
        }
    }

    void fixEMatrixChain(bool pairlist_update, Molecule chain) {
        // FIX ENERGY MATRIX
        bool partOfChain;
        if(pairlist_update) {
            vector<double>::iterator it = changes.begin();
            for(unsigned int i=0; i < chain.size(); i++) {
                for(int j=0; j < conf->neighborList[chain[i]].neighborCount; j++) {
                    partOfChain = false;
                    for(unsigned int k=0; k<chain.size(); k++)
                        if(chain[k] == conf->neighborList[chain[i]].neighborID[j])
                            partOfChain = true;

                    if(!partOfChain) {
                        energyMatrix->operator [](chain[i])[conf->neighborList[chain[i]].neighborID[j]] = *it;
                        energyMatrix->operator [](conf->neighborList[chain[i]].neighborID[j])[chain[i]] = *it;
                        it++;
                    }
                }
            }
        } else {
            for(unsigned int i=0; i<chain.size(); i++) { // for all particles in molecule
                // for cycles => pair potential with all particles except those in molecule
                for (int j = 0; j < chain[0]; j++) { // pair potential with all particles from 0 to the first one in molecule
                    energyMatrix->operator [](chain[i])[j] = changes[i * conf->pvec.size() + j];
                    energyMatrix->operator [](j)[chain[i]] = changes[i * conf->pvec.size() + j];
                }

                for (int j = chain[chain.size()-1] + 1; j < (long)conf->pvec.size(); j++) {
                    energyMatrix->operator [](chain[i])[j] = changes[i * conf->pvec.size() + j];
                    energyMatrix->operator [](j)[chain[i]] = changes[i * conf->pvec.size() + j];
                }

            }
        }

    }

    void resizeEMatrix() {
        energyMatrix->resize(conf->pvec.size());
        for(unsigned int i=0; i<conf->pvec.size(); i++) {
            energyMatrix->operator [](i).resize(conf->pvec.size());
        }

        energyMatrixTrial->resize(conf->pvec.size());
        for(unsigned int i=0; i<conf->pvec.size(); i++) {
            energyMatrixTrial->operator [](i).resize(conf->pvec.size());
        }
    }
};

template<typename pairEFce>
class TotalE {
public:
    pairEFce pairE;

    ExternalEnergyCalculator exterE;

    EnergyMatrix eMat;

    bool pairListUpdate;
    Conf* conf;

    TotalE(Sim * sim, Conf * conf) : pairE(pairEFce(&conf->geo)),
                                                   exterE(ExternalEnergyCalculator(&conf->geo.box)), eMat(EnergyMatrix(conf)),
                                                   pairListUpdate(sim->pairlist_update), conf(conf) {}

    /**
     * @brief Calculates energy between all pairs, generates energy matrix. Returns energy
     * @param energyMatrix - We want to set either trial(for move) or actual matrix(re-initialization)
     * @return
     */
    virtual double allToAll(std::vector< std::vector<double> >* energyMatrix) {return 0.0;}

    /**
     * @brief allToAll - Calculation with the use of energy matrix
     * @return
     */
    virtual double allToAll() {return 0.0;}

    /**
     * @brief Calculates energy between particle target and rest, generates energy matrix. Returns energy
     * @return
     */
    virtual double oneToAll(int target, vector<double> *changes) {return 0.0;}

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
    virtual double mol2others(Molecule &mol, vector<double> *changes) {return 0.0;}

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
        double energy=0.0;
        unsigned int i=0;

        for(i=0; i<mol.size(); i++) {
            if (topo.exter.exist)
                energy += this->exterE.extere2(&mol[i]);

            for(unsigned int j=0; j < conf->pvec.size(); j++)
                energy += pairE(&mol[i], &conf->pvec[j]);
        }
        return energy;
    }

    /**
     * @brief Calculates inner chain energy
     */
    double chainInner(vector<Particle >& chain) {
        double energy = 0.0;

        vector<ConList> con;
        for(unsigned int j=0; j<chain.size(); j++) {
            con.push_back(ConList() );

            if (j > 0) //if this is not first particle fill tail bond
                con[j].conlist[0] = &chain[j-1];

            if ( j+1 < chain.size() ) //if there is a next particle fill it to head bond
                con[j].conlist[1] = &chain[j+1];

            if (j > 1) //if this is not second or first particle fill second tail bond
                con[j].conlist[2] = &chain[j-2];

            if ( j+2 < chain.size() ) //if there is a second next particle fill it second neighbour
                con[j].conlist[3] = &chain[j+2];
        }

        for (unsigned int i=0; i<chain.size(); i++)
            for (unsigned int j=i+1; j<chain.size(); j++)
                energy += (pairE)(&chain[i], &chain[j],  &con[i]);

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
    double oneToAll(Particle* target, ConList* conlist=NULL, Neighbors* neighborList=NULL) {
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
    using TotalE<pairEFce>::mol2others;
    using TotalE<pairEFce>::conf;
    using TotalE<pairEFce>::pairE;
    using TotalE<pairEFce>::pairListUpdate;
    using TotalE<pairEFce>::eMat;

    TotalEMatrix(Sim * sim, Conf * conf) : TotalE<PairE>(sim, conf) {}

    double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
        if(conf->pvec.empty()) return 0.0;
        double energy=0.0;
        ConList conlist;
        assert(energyMatrix != NULL);

        for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
            conlist = conf->pvec.getConlist(i);
            for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
                (*energyMatrix)[i][j] = (pairE)(&conf->pvec[i], &conf->pvec[j], &conlist);
                (*energyMatrix)[j][i] = (*energyMatrix)[i][j];
                energy += (*energyMatrix)[i][j];
            }

            if (topo.exter.exist) //for every particle add interaction with external potential
                energy += this->exterE.extere2(&conf->pvec[i]);
        }

        if (topo.exter.exist && !conf->pvec.empty()) //add interaction of last particle with external potential
            energy += this->exterE.extere2(&conf->pvec.back());

        return energy;
    }

    double allToAll() override {
        if(conf->pvec.empty()) return 0.0;

        double energy = 0.0;
        for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
            for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
                energy += eMat.energyMatrix->operator [](i)[j];
            }

            if (topo.exter.exist) //for every particle add interaction with external potential
                energy += this->exterE.extere2(&conf->pvec[i]);
        }

        if (topo.exter.exist && !conf->pvec.empty()) //add interaction of last particle with external potential
            energy += this->exterE.extere2(&conf->pvec.back());

        return energy;
    }

    double oneToAll(int target, vector<double> *changes) override {
        assert(target >= 0 && target < conf->pvec.size());
        assert(changes != NULL);
        eMat.changes.clear();
        assert(eMat.changes.empty());

        double energy=0.0;
        long i;
        ConList conlist = conf->pvec.getConlist(target);

        if (pairListUpdate) {
            for (i = 0; i < conf->neighborList[target].neighborCount; i++) {
                changes->push_back(pairE(&conf->pvec[target], &conf->pvec[ conf->neighborList[target].neighborID[i] ], (conlist.isEmpty) ? nullptr : &conlist));
                energy += changes->back();
            }
        } else {
            changes->resize(conf->pvec.size());
            (*changes)[target] = 0.0;
            for (i = 0; i < target; i++) {
                (*changes)[i] = pairE(&conf->pvec[target], &conf->pvec[i], &conlist);
                energy += (*changes)[i];
            }
            for (i = target+1; i < (long)conf->pvec.size(); i++) {
                (*changes)[i] = pairE(&conf->pvec[target], &conf->pvec[i], &conlist);
                energy += (*changes)[i];
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(&conf->pvec[target]);

        return energy;
    }

    double oneToAll(int target) override {
        assert(target >=0 && target < conf->pvec.size());
        double energy=0.0;

        if (pairListUpdate) {
            for (long i = 0; i < conf->neighborList[target].neighborCount; i++) {
                energy += eMat.energyMatrix->operator [](target)[conf->neighborList[target].neighborID[i]];
            }
        } else {
            for (long i = 0; i < (long)conf->pvec.size(); i++) {
                if(target != i) {
                    energy += eMat.energyMatrix->operator [](target)[i];
                }
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(&conf->pvec[target]);

        return energy;
    }

    double mol2others(Molecule &mol) override {
        double energy=0.0;
        long i = 0;
        bool partOfChain;

        if (pairListUpdate) {

            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < conf->neighborList[mol[j]].neighborCount; i++) {
                    partOfChain = false;
                    for(unsigned int k=0; k<mol.size(); k++)
                        if(mol[k] == conf->neighborList[mol[j]].neighborID[i])
                            partOfChain = true;

                    if(!partOfChain)
                        energy += eMat.energyMatrix->operator [](mol[j])[conf->neighborList[mol[j]].neighborID[i]];
                }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }
        } else {

            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < mol[0]; i++)  // pair potential with all particles from 0 to the first one in molecule
                    energy += eMat.energyMatrix->operator [](mol[j])[i];

                for (long i = mol[mol.size()-1] + 1; i < (long)conf->pvec.size(); i++)
                    energy += eMat.energyMatrix->operator [](mol[j])[i];

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }
        }
        return energy;
    }

    double mol2others(Molecule &mol, vector<double> *changes) override {
        eMat.changes.clear();
        assert(eMat.changes.empty());

        double energy=0.0;
        long i = 0;
        bool partOfChain;

        if (pairListUpdate) {

            //
            // SAVE the changes to energy to easily swap them later -> instance when oneToAll(target, false, some vector) used
            //
            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < conf->neighborList[mol[j]].neighborCount; i++) {
                    partOfChain = false;
                    for(unsigned int k=0; k<mol.size(); k++)
                        if(mol[k] == conf->neighborList[mol[j]].neighborID[i])
                            partOfChain = true;

                    if(!partOfChain) {
                        changes->push_back( pairE(&conf->pvec[mol[j]], &conf->pvec[conf->neighborList[mol[j]].neighborID[i]]) );
                        energy += changes->back();
                    }
                }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }

        } else {

            changes->resize( mol.size() * conf->pvec.size() );
            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < mol[0]; i++) { // pair potential with all particles from 0 to the first one in molecule
                    (*changes)[j*conf->pvec.size() + i] = pairE(&conf->pvec[mol[j]], &conf->pvec[i]);
                    energy += (*changes)[j*conf->pvec.size() + i];
                }

                for (long i = mol[mol.size()-1] + 1; i < (long)conf->pvec.size(); i++) {
                    (*changes)[j*conf->pvec.size() + i] = pairE(&conf->pvec[mol[j]], &conf->pvec[i]);
                    energy += (*changes)[j*conf->pvec.size() + i];
                }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&conf->pvec[mol[j]]);
            }
        }
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
    using TotalE<pairEFce>::eMat;

    TotalEFull(Sim * sim, Conf * conf) : TotalE<pairEFce>(sim, conf) {}

    double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
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

        double energy=0.0;

        for (int i = 0; i < (int)conf->pvec.size(); i++) {
            if(target != i)
                energy += (pairE)(&conf->pvec[target], &conf->pvec[i], &conlist);
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

    double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
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
class TestE : public TotalE<PairE> {
public:
    using TotalE<PairE>::mol2others;
    using TotalE<PairE>::eMat;

    TotE tested;
    TotEControl control;

    TestE(Sim * sim, Conf * conf) : TotalE<PairE>(sim, conf), tested(sim, conf), control(sim, conf)  {}

    bool test(double a, double b) {
        return fabs(a - b) > 1e-4;
    }

    double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
        double e = 0.0, eControl = 0.0;

        e = tested.allToAll(energyMatrix);
        eControl = control.allToAll(energyMatrix);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, allToAll energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double allToAll() override {
        double e = 0.0, eControl = 0.0;

        e = tested.allToAll();
        eControl = control.allToAll();

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, allToAll energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double oneToAll(int target, vector<double> *changes) override {
        double e = 0.0, eControl = 0.0;

        e = tested.oneToAll(target, changes);
        eControl = control.oneToAll(target, changes);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, oneToAll energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double oneToAll(int target) override {
        double e = 0.0, eControl = 0.0;

        e = tested.oneToAll(target);
        eControl = control.oneToAll(target);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, oneToAll energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double mol2others(Molecule &mol, vector<double> *changes) override {
        double e = 0.0, eControl = 0.0;

        e = tested.mol2others(mol, changes);
        eControl = control.mol2others(mol, changes);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, mol2others energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double mol2others(Molecule &mol) override {
        double e = 0.0, eControl = 0.0;

        e = tested.mol2others(mol);
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

    inline double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
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

    double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
        if(conf->pvec.empty()) return 0.0;
        Vector r_cm;
        double dotrcm;
        double energy=0.0;
        assert(energyMatrix != NULL);

        for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
            for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
                r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[j].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].sigmaSq) {
                    (*energyMatrix)[i][j] = std::numeric_limits<double>::infinity();
                    energy = std::numeric_limits<double>::infinity();
                } else
                    (*energyMatrix)[i][j] = 0.0;
                (*energyMatrix)[j][i] = (*energyMatrix)[i][j];
            }
        }
        return energy;
    }

    double allToAll() override {
        if(conf->pvec.empty()) return 0.0;

        for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
            for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
                if(conf->energyMatrix->operator [](i)[j] == std::numeric_limits<double>::infinity())
                    return std::numeric_limits<double>::infinity();
            }
        }
        return 0.0;
    }

    double oneToAll(int target, vector<double> *changes) override {
        assert(target >= 0 && target < conf->pvec.size());
        assert(changes != NULL);

        Vector r_cm;
        double dotrcm;
        double energy=0.0;
        long i;

        if (pairListUpdate) {
            for (i = 0; i < conf->neighborList[target].neighborCount; i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[conf->neighborList[target].neighborID[i]].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[conf->neighborList[target].neighborID[i]].type].sigmaSq) {
                    changes->push_back( std::numeric_limits<double>::infinity() );
                    energy = std::numeric_limits<double>::infinity();
                } else
                    changes->push_back( 0.0 );
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

        if (pairListUpdate) {
            for (long i = 0; i < conf->neighborList[target].neighborCount; i++) {
                if(conf->energyMatrix->operator [](target)[conf->neighborList[target].neighborID[i]] == std::numeric_limits<double>::infinity())
                    return std::numeric_limits<double>::infinity();
            }
        } else {
            for(unsigned int i=0; i < conf->pvec.size(); ++i) {
                if(i != target) {
                    if( conf->energyMatrix->operator [](target)[i] == std::numeric_limits<double>::infinity() )
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

    inline double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
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

    double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
        if(conf->pvec.empty()) return 0.0;
        Vector r_cm;
        double dotrcm;
        double energy=0.0;
        assert(energyMatrix != NULL);

        for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
            for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
                r_cm = conf->geo.image(&conf->pvec[i].pos, &conf->pvec[j].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].sigmaSq*(6.25))
                    (*energyMatrix)[i][j] = topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].A * pow(dotrcm, -6)
                        - topo.ia_params[conf->pvec[i].type][conf->pvec[j].type].B * pow(dotrcm, -3);
                else
                    (*energyMatrix)[i][j] = 0.0;
                (*energyMatrix)[j][i] = (*energyMatrix)[i][j];
                energy += (*energyMatrix)[i][j];
            }
        }
        return energy;
    }

    double allToAll() override {
        if(conf->pvec.empty()) return 0.0;
        double energy = 0.0;

        for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
            for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
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
            for (int i = 0; i < conf->neighborList[target].neighborCount; i++) {
                r_cm = conf->geo.image(&conf->pvec[target].pos, &conf->pvec[conf->neighborList[target].neighborID[i]].pos);
                dotrcm = r_cm.dot(r_cm);
                if(dotrcm < topo.ia_params[conf->pvec[target].type][conf->pvec[conf->neighborList[target].neighborID[i]].type].sigmaSq*(6.25)) {
                    changes->push_back( topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].A * pow(dotrcm, -6)
                                      - topo.ia_params[conf->pvec[target].type][conf->pvec[ conf->neighborList[target].neighborID[i] ].type].B * pow(dotrcm, -3) );
                    energy += changes->back();
                } else
                    changes->push_back( 0.0 );
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
                energy += this->eMat.energyMatrix->operator [](target)[conf->neighborList[target].neighborID[i]];
            }
        } else {
            for(unsigned int i=0; i < conf->pvec.size(); ++i) {
                if(i != target) {
                    energy += this->eMat.energyMatrix->operator [](target)[i];
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
//typedef TotalEMatrix<PairEnergyCalculator> TotalEnergyCalculator;         // Energy matrix optimization
typedef TotalEMatrix<PairE> TotalEnergyCalculator;                        // Energy matrix optimization
//typedef TotalEFull<PairE> TotalEnergyCalculator;                          // Full calculation
//typedef TestE<TotalEFull<PairE>, TotalEHardSphereEMatrix<PairE> > TotalEnergyCalculator;         // Test of Pair energy
//typedef TotalEFullSymetry<PairEnergyCalculator> TotalEnergyCalculator;         // Test of Pair energy symetry

//
//  Specialized Total energy calculators
//
//typedef TotalEHardSphereFull<PairE> TotalEnergyCalculator;         // Hardspheres, full calculation
//typedef TotalEHardSphereEMatrix<PairE> TotalEnergyCalculator;         // Hardspheres, energy matrix calculation
//typedef TotalELJFull<PairE> TotalEnergyCalculator;         // LJ, full calculation //cut-off distance of rc = 2.5σ,
//typedef TotalELJEMatrix<PairE> TotalEnergyCalculator;         // LJ, energy matrix calculation //cut-off distance of rc = 2.5σ,



#endif // TOTALENERGYCALCULATOR_H
