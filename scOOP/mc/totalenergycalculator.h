/** @file totalenergycalculator.h*/

#ifndef TOTALENERGYCALCULATOR_H
#define TOTALENERGYCALCULATOR_H

#include <iomanip>

#include "pairenergycalculator.h"
#include "paire.h"
#include "externalenergycalculator.h"
#include "../structures/sim.h"

using namespace std;

template<typename PairE>
class TotalE {
public:
    PairE pairE;

    ExternalEnergyCalculator exterE;

    bool pairListUpdate;
    Conf* conf;

    TotalE(Sim * sim, Conf * conf) : pairE(PairE(&conf->geo)),
                                                   exterE(ExternalEnergyCalculator(&conf->geo.box)),
                                                   pairListUpdate(sim->pairlist_update), conf(conf) {}

    /**
     * @brief Calculates energy between all pairs, generates energy matrix. Returns energy
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
        ConList conlist = this->conf->pvec.getConlist(part1);
        return this->pairE(&this->conf->pvec[part1], &this->conf->pvec[part2], &conlist);
    }

    double p2p(Particle* part1, int part2) {
        ConList conlist = this->conf->pvec.getConlist(part2);
        return this->pairE(part1, &this->conf->pvec[part2], &conlist);
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

            for(unsigned int j=0; j < this->conf->pvec.size(); j++)
                energy += this->pairE(&mol[i], &this->conf->pvec[j]);
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
                energy += (this->pairE)(&chain[i], &chain[j],  &con[i]);

        return energy;
    }

    double chainInner(Molecule& chain) {
        double energy = 0.0;
        ConList conlist;

        for (unsigned int i=0; i<chain.size(); i++) {
            for (unsigned int j=i+1; j<chain.size(); j++) {
                conlist = this->conf->pvec.getConlist(chain[i]);
                energy += (this->pairE)(&this->conf->pvec[chain[i]], &this->conf->pvec[chain[j]],  &conlist);
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

        if (this->pairListUpdate && neighborList!=NULL) {
            for (i = 0; i < (unsigned long)neighborList->neighborCount; i++){
                energy += (this->pairE)(target, &this->conf->pvec[ neighborList->neighborID[i] ], conlist );
            }
        } else {
            for (i = 0; i < this->conf->pvec.size(); i++) {
                if(target != &this->conf->pvec[i]) {
                    energy += this->pairE(target, &this->conf->pvec[i], conlist);
                }
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(target);

        return energy;
    }
};

template<typename PairE>
class TotalEMatrix : public TotalE<PairE>
{
public:
    using TotalE<PairE>::mol2others;

    TotalEMatrix(Sim * sim, Conf * conf) : TotalE<PairE>(sim, conf) {}

    double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
        if(this->conf->pvec.empty()) return 0.0;
        double energy=0.0;
        ConList conlist;
        assert(energyMatrix != NULL);

        for (unsigned int i = 0; i < this->conf->pvec.size() - 1; i++) {
            conlist = this->conf->pvec.getConlist(i);
            for (unsigned long j = i + 1; j < this->conf->pvec.size(); j++) {
                (*energyMatrix)[i][j] = (this->pairE)(&this->conf->pvec[i], &this->conf->pvec[j], &conlist);
                (*energyMatrix)[j][i] = (*energyMatrix)[i][j];
                energy += (*energyMatrix)[i][j];
            }

            if (topo.exter.exist) //for every particle add interaction with external potential
                energy += this->exterE.extere2(&this->conf->pvec[i]);
        }

        if (topo.exter.exist && !this->conf->pvec.empty()) //add interaction of last particle with external potential
            energy += this->exterE.extere2(&this->conf->pvec.back());

        return energy;
    }

    double allToAll() override {
        if(this->conf->pvec.empty()) return 0.0;

        double energy = 0.0;
        for (unsigned int i = 0; i < this->conf->pvec.size() - 1; i++) {
            for (unsigned long j = i + 1; j < this->conf->pvec.size(); j++) {
                energy += this->conf->energyMatrix->operator [](i)[j];
            }

            if (topo.exter.exist) //for every particle add interaction with external potential
                energy += this->exterE.extere2(&this->conf->pvec[i]);
        }

        if (topo.exter.exist && !this->conf->pvec.empty()) //add interaction of last particle with external potential
            energy += this->exterE.extere2(&this->conf->pvec.back());

        return energy;
    }

    double oneToAll(int target, vector<double> *changes) override {
        assert(target >= 0 && target < this->conf->pvec.size());
        assert(changes != NULL);

        double energy=0.0;
        long i;
        ConList conlist = this->conf->pvec.getConlist(target);

        if (this->pairListUpdate) {
            for (i = 0; i < this->conf->neighborList[target].neighborCount; i++){
                changes->push_back(this->pairE(&this->conf->pvec[target], &this->conf->pvec[ this->conf->neighborList[target].neighborID[i] ], (conlist.isEmpty) ? nullptr : &conlist));
                energy += changes->back();
            }
        } else {
            changes->resize(this->conf->pvec.size());
            (*changes)[target] = 0.0;
            for (i = 0; i < target; i++) {
                (*changes)[i] = this->pairE(&this->conf->pvec[target], &this->conf->pvec[i], &conlist);
                energy += (*changes)[i];
            }
            for (i = target+1; i < (long)this->conf->pvec.size(); i++) {
                (*changes)[i] = this->pairE(&this->conf->pvec[target], &this->conf->pvec[i], &conlist);
                energy += (*changes)[i];
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(&this->conf->pvec[target]);

        return energy;
    }

    double oneToAll(int target) override {
        assert(target >=0 && target < this->conf->pvec.size());
        double energy=0.0;

        if (this->pairListUpdate) {
            for (long i = 0; i < this->conf->neighborList[target].neighborCount; i++) {
                energy += this->conf->energyMatrix->operator [](target)[this->conf->neighborList[target].neighborID[i]];
            }
        } else {
            for (long i = 0; i < (long)this->conf->pvec.size(); i++) {
                if(target != i) {
                    energy += this->conf->energyMatrix->operator [](target)[i];
                }
            }
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(&this->conf->pvec[target]);

        return energy;
    }

    double mol2others(Molecule &mol) override {
        double energy=0.0;
        long i = 0;
        bool partOfChain;

        if (this->pairListUpdate) {

            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < this->conf->neighborList[mol[j]].neighborCount; i++) {
                    partOfChain = false;
                    for(unsigned int k=0; k<mol.size(); k++)
                        if(mol[k] == this->conf->neighborList[mol[j]].neighborID[i])
                            partOfChain = true;

                    if(!partOfChain)
                        energy += this->conf->energyMatrix->operator [](mol[j])[this->conf->neighborList[mol[j]].neighborID[i]];
                }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&this->conf->pvec[mol[j]]);
            }
        } else {

            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < mol[0]; i++)  // pair potential with all particles from 0 to the first one in molecule
                    energy += this->conf->energyMatrix->operator [](mol[j])[i];

                for (long i = mol[mol.size()-1] + 1; i < (long)this->conf->pvec.size(); i++)
                    energy += this->conf->energyMatrix->operator [](mol[j])[i];

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&this->conf->pvec[mol[j]]);
            }
        }
        return energy;
    }

    double mol2others(Molecule &mol, vector<double> *changes) override {
        double energy=0.0;
        long i = 0;
        bool partOfChain;

        if (this->pairListUpdate) {

            //
            // SAVE the changes to energy to easily swap them later -> instance when oneToAll(target, false, some vector) used
            //
            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < this->conf->neighborList[mol[j]].neighborCount; i++) {
                    partOfChain = false;
                    for(unsigned int k=0; k<mol.size(); k++)
                        if(mol[k] == this->conf->neighborList[mol[j]].neighborID[i])
                            partOfChain = true;

                    if(!partOfChain) {
                        changes->push_back( this->pairE(&this->conf->pvec[mol[j]], &this->conf->pvec[this->conf->neighborList[mol[j]].neighborID[i]]) );
                        energy += changes->back();
                    }
                }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&this->conf->pvec[mol[j]]);
            }

        } else {

            changes->resize( mol.size() * this->conf->pvec.size() );
            for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
                for (i = 0; i < mol[0]; i++) { // pair potential with all particles from 0 to the first one in molecule
                    (*changes)[j*this->conf->pvec.size() + i] = this->pairE(&this->conf->pvec[mol[j]], &this->conf->pvec[i]);
                    energy += (*changes)[j*this->conf->pvec.size() + i];
                }

                for (long i = mol[mol.size()-1] + 1; i < (long)this->conf->pvec.size(); i++) {
                    (*changes)[j*this->conf->pvec.size() + i] = this->pairE(&this->conf->pvec[mol[j]], &this->conf->pvec[i]);
                    energy += (*changes)[j*this->conf->pvec.size() + i];
                }

                if (topo.exter.exist) //add interaction with external potential
                    energy += this->exterE.extere2(&this->conf->pvec[mol[j]]);
            }
        }
        return energy;
    }
};


template<typename pairE>
class TotalEFull : public TotalE<pairE>
{
public:
    using TotalE<PairE>::mol2others;

    TotalEFull(Sim * sim, Conf * conf) : TotalE<pairE>(sim, conf) {}

    double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
        return allToAll();
    }

    double allToAll() override {
        if(this->conf->pvec.empty()) return 0.0;

        double energy=0.0;
        ConList conlist;
        for (unsigned int i = 0; i < this->conf->pvec.size() - 1; i++) {
            conlist = this->conf->pvec.getConlist(i);
            for (unsigned long j = i + 1; j < this->conf->pvec.size(); j++) {
                energy += (this->pairE)(&this->conf->pvec[i], &this->conf->pvec[j], &conlist);
            }

            if (topo.exter.exist) //for every particle add interaction with external potential
                energy += this->exterE.extere2(&this->conf->pvec[i]);
        }

        if (topo.exter.exist && !this->conf->pvec.empty()) //add interaction of last particle with external potential
            energy += this->exterE.extere2(&this->conf->pvec.back());

        return energy;
    }

    double oneToAll(int target, vector<double> *changes) override {
        return oneToAll(target);
    }

    double oneToAll(int target) override {
        ConList conlist = this->conf->pvec.getConlist(target);

        double energy=0.0;

        for (int i = 0; i < (int)this->conf->pvec.size(); i++) {
            if(target != i)
                energy += (this->pairE)(&this->conf->pvec[target], &this->conf->pvec[i], &conlist);
        }

        if (topo.exter.exist) //add interaction with external potential
            energy += this->exterE.extere2(&this->conf->pvec[target]);

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
                energy += this->pairE(&this->conf->pvec[mol[j]], &this->conf->pvec[i]);

            for (i = mol[mol.size()-1] + 1; i < (long)this->conf->pvec.size(); i++)
                energy += this->pairE(&this->conf->pvec[mol[j]], &this->conf->pvec[i]);

            if (topo.exter.exist) //add interaction with external potential
                energy += this->exterE.extere2(&this->conf->pvec[mol[j]]);
        }
        return energy;
    }
};

template<typename PairE, typename PairEControl>
class TestE : public TotalE<PairE> {
public:
    using TotalE<PairE>::mol2others;

    TotalEFull<PairE> tested;
    TotalEMatrix<PairEControl> control;

    TestE(Sim * sim, Conf * conf) : TotalE<PairE>(sim, conf), tested(sim, conf), control(sim, conf)  {}

    bool test(double a, double b) {
        return fabs(a - b) > 1e-7;
    }

    double allToAll(std::vector< std::vector<double> >* energyMatrix) override {
        double e = 0.0, eControl = 0.0;

        e = tested.allToAll(energyMatrix);
        eControl = control.allToAll(energyMatrix);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double allToAll() override {
        double e = 0.0, eControl = 0.0;

        e = tested.allToAll();
        eControl = control.allToAll();

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double oneToAll(int target, vector<double> *changes) override {
        double e = 0.0, eControl = 0.0;

        e = tested.oneToAll(target, changes);
        eControl = control.oneToAll(target, changes);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double oneToAll(int target) override {
        double e = 0.0, eControl = 0.0;

        e = tested.oneToAll(target);
        eControl = control.oneToAll(target);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double mol2others(Molecule &mol, vector<double> *changes) override {
        double e = 0.0, eControl = 0.0;

        e = tested.mol2others(mol, changes);
        eControl = control.mol2others(mol, changes);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, energy: " << e << " control: " << eControl << endl;

        return e;
    }

    double mol2others(Molecule &mol) override {
        double e = 0.0, eControl = 0.0;

        e = tested.mol2others(mol);
        eControl = control.mol2others(mol);

        if(test(e,eControl))
            cout << std::setprecision(10) << "Error, energy: " << e << " control: " << eControl << endl;

        return e;
    }
};

//typedef TotalEMatrix<PairEnergyCalculator> TotalEnergyCalculator; // Energy matrix optimization
typedef TotalEMatrix<PairE> TotalEnergyCalculator; // Energy matrix optimization
//typedef TotalEFull<PairE> TotalEnergyCalculator;        // Full calculation
//typedef TestE<PairE, PairEnergyCalculator> TotalEnergyCalculator;        // Test of Pair energy

#endif // TOTALENERGYCALCULATOR_H
