#include "totalenergycalculator.h"

#include <stdlib.h>
#include <iostream>
#include <iomanip>

extern Topo topo;


/************************************************************************************************/
/*                                                                                              */
/*  ENERGY FUNCTIONS:                                                                           */
/*  always provide 3 functions, one utilizing energy matrix, second which doesnt for testing    */
/*  and third for update of part of energy matrix                                               */
/*                                                                                              */
/*  Further diverge all function but the ones for testing on the usage of pairlist              */
/*                                                                                              */
/************************************************************************************************/


double TotalEnergyCalculator::allToAll(std::vector< std::vector<double> >* energyMatrix) {
    if(conf->pvec.empty()) return 0.0;
    double energy=0.0;
    ConList conlist;

    assert(energyMatrix != NULL);

    //
    //  DONT USE ENERGY MATRIX, CALCULATION FROM CONFIGURATION AND AT THE SAME TIME GENERATION OF ENERGY MATRIX
    //
    for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
        conlist = conf->pvec.getConlist(i);
        for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
            (*energyMatrix)[i][j] = (pairE)(&conf->pvec[i], &conf->pvec[j], &conlist);
            (*energyMatrix)[j][i] = (*energyMatrix)[i][j];
            energy += (*energyMatrix)[i][j];
        }

        //for every particle add interaction with external potential
        if (topo.exter.exist)
            energy += extere2(i);
    }


    //add interaction of last particle with external potential
    if (topo.exter.exist && !conf->pvec.empty())
        energy += exterE.extere2(&conf->pvec.back());

    return energy;
}

double TotalEnergyCalculator::allToAll() {
    if(conf->pvec.empty()) return 0.0;
    //
    //  CALCULATE ENERGIES FROM MATRIX
    //
    double energy = 0.0;
    for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
        for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
            energy += conf->energyMatrix->operator [](i)[j];
        }

        //for every particle add interaction with external potential
        if (topo.exter.exist)
            energy += extere2(i);
    }

    //add interaction of last particle with external potential
    if (topo.exter.exist && !conf->pvec.empty())
        energy += exterE.extere2(&conf->pvec.back());

#ifndef NDEBUG
    double energy2 = allToAllBasic();
    assert(energy+0.0000001 >= energy2 && energy-0.0000001 <= energy2);
#endif
    return energy;
}

double TotalEnergyCalculator::allToAllBasic() {
    if(conf->pvec.empty()) return 0.0;
    //
    //  STANDART CALCULATION FROM CONFIGURATION, TEST for energy Matrix
    //
    double energy=0.0;
    ConList conlist;
    for (unsigned int i = 0; i < conf->pvec.size() - 1; i++) {
        conlist = conf->pvec.getConlist(i);
        for (unsigned long j = i + 1; j < conf->pvec.size(); j++) {
            energy += (pairE)(&conf->pvec[i], &conf->pvec[j], &conlist);
        }

        //for every particle add interaction with external potential
        if (topo.exter.exist)
            energy += extere2(i);
    }

    //add interaction of last particle with external potential
    if (topo.exter.exist && !conf->pvec.empty())
        energy += exterE.extere2(&conf->pvec.back());

    return energy;
}


////////////////////////////////////////////////////////////////////////////////////////////////


double TotalEnergyCalculator::oneToAll(int target, vector<double>* changes) {
    assert(target >= 0 && target < conf->pvec.size());
    assert(changes != NULL);

    double energy=0.0;
    long i;
    ConList conlist = conf->pvec.getConlist(target);

    if (pairListUpdate) {
        for (i = 0; i < conf->neighborList[target].neighborCount; i++){
            changes->push_back(pairE(&conf->pvec[target], &conf->pvec[ conf->neighborList[target].neighborID[i] ], (conlist.isEmpty) ? nullptr : &conlist));
            //changes->push_back(pairEControl(&conf->pvec[target], &conf->pvec[ conf->neighborList[target].neighborID[i] ], &conlist));

            energy += changes->back();

#ifndef NDEBUG // check calculation
            double ener = pairEControl(&conf->pvec[target], &conf->pvec[ conf->neighborList[target].neighborID[i] ], &conlist);
            double ener2 = pairE(&conf->pvec[target], &conf->pvec[ conf->neighborList[target].neighborID[i] ], &conlist);
            if(fabs(ener - ener2) > 1e-6) {
                cout << "New vs Control: " << std::setprecision(20) << ener << " != " << ener2 << endl;
                exit(0);
            }
            assert(fabs(ener - ener2) < 1e-6 || !(cout << "New vs Control: " << ener << " != " << ener2 << endl) );
            assert(fabs(ener - changes->back()) < 1e-6 || !(cout << "New vs Stored: " << ener << " != " << changes->back() << endl) );
            assert(fabs(ener2 - changes->back()) < 1e-6 || !(cout << "Control vs Stored: " << ener2 << " != " << changes->back() << endl) );
#endif
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

    //add interaction with external potential
    if (topo.exter.exist)
        energy += extere2(target);

#ifndef NDEBUG
    double energy2 = oneToAllBasic(target);
    assert(energy+0.0000001 >= energy2 && energy-0.0000001 <= energy2);
#endif

    return energy;
}


double TotalEnergyCalculator::oneToAll(int target) {
    //
    // USE PRECOMPUTED ENERGY -> instance when oneToAll(target, true) used
    //
    assert(target >=0 && target < conf->pvec.size());
    double energy=0.0;
    long i;

    if (pairListUpdate) {
        for (i = 0; i < conf->neighborList[target].neighborCount; i++) {
            energy += conf->energyMatrix->operator [](target)[conf->neighborList[target].neighborID[i]];
#ifndef NDEBUG
            ConList conlist = conf->pvec.getConlist(target);
            double ener = pairE(&conf->pvec[target], &conf->pvec[ conf->neighborList[target].neighborID[i] ], &conlist);
            double ener2 = pairEControl(&conf->pvec[target], &conf->pvec[ conf->neighborList[target].neighborID[i] ], &conlist);
            double matrE = conf->energyMatrix->operator [](target)[conf->neighborList[target].neighborID[i]];
            assert(fabs(ener - matrE) < 1e-10 || !(cout << ener << " " << ener2 << " != " << matrE << endl) );
#endif
        }
    } else {
        for (i = 0; i < (long)conf->pvec.size(); i++) {
            if(target != i) {
                energy += conf->energyMatrix->operator [](target)[i];
#ifndef NDEBUG
                ConList conlist = conf->pvec.getConlist(target);
                double ener = pairE(&conf->pvec[target], &conf->pvec[i], &conlist);
                double matrE = conf->energyMatrix->operator [](target)[i];
                assert(ener + 0.0000001 >= matrE && ener - 0.0000001 <= matrE);
#endif
            }
        }
    }

    //add interaction with external potential
    if (topo.exter.exist)
        energy += extere2(target);

#ifndef NDEBUG
    double energy2 = oneToAllBasic(target);
    assert(energy+0.0000001 >= energy2 && energy-0.0000001 <= energy2);
#endif

    return energy;
}


double TotalEnergyCalculator::oneToAllBasic(int target) {
    ConList conlist = conf->pvec.getConlist(target);

    double energy=0.0;

    for (int i = 0; i < (int)conf->pvec.size(); i++) {
        if(target != i) {
            energy += (pairE)(&conf->pvec[target], &conf->pvec[i], &conlist);
        }
    }

    //add interaction with external potential
    if (topo.exter.exist)
        energy += extere2(target);

    return energy;
}


////////////////////////////////////////////////////////////////////////////////////////////////


double TotalEnergyCalculator::mol2others(Molecule &mol) {
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

                if(!partOfChain) {
                    energy += conf->energyMatrix->operator [](mol[j])[conf->neighborList[mol[j]].neighborID[i]];
#ifndef NDEBUG
                    double ener = pairE(&conf->pvec[mol[j]], &conf->pvec[conf->neighborList[mol[j]].neighborID[i]]);
                    double matrE = conf->energyMatrix->operator [](mol[j])[conf->neighborList[mol[j]].neighborID[i]];
                    double matrE2 = conf->energyMatrix->operator [](conf->neighborList[mol[j]].neighborID[i])[mol[j]];
                    assert(matrE == matrE2);
                    assert(ener + 0.0000001 >= matrE && ener - 0.0000001 <= matrE);
#endif
                }
            }

            //add interaction with external potential
            if (topo.exter.exist)
                energy += extere2(mol[j]);
        }
    } else {

        for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
            // for cycles => pair potential with all particles except those in molecule
            for (i = 0; i < mol[0]; i++) { // pair potential with all particles from 0 to the first one in molecule
                energy += conf->energyMatrix->operator [](mol[j])[i];
#ifndef NDEBUG
                ConList conlist = conf->pvec.getConlist(i);
                double ener = pairE(&conf->pvec[mol[j]], &conf->pvec[i], &conlist);
                double matrE = conf->energyMatrix->operator [](mol[j])[i];
                assert(ener + 0.0000001 >= matrE && ener - 0.0000001 <= matrE);
#endif
            }

            for (long i = mol[mol.size()-1] + 1; i < (long)conf->pvec.size(); i++) {
                energy += conf->energyMatrix->operator [](mol[j])[i];
#ifndef NDEBUG
                ConList conlist = conf->pvec.getConlist(i);
                double ener = pairE(&conf->pvec[mol[j]], &conf->pvec[i], &conlist);
                double matrE = conf->energyMatrix->operator [](mol[j])[i];
                assert(ener + 0.0000001 >= matrE && ener - 0.0000001 <= matrE);
#endif
            }

            //add interaction with external potential
            if (topo.exter.exist)
                energy += extere2(mol[j]);
        }
    }

#ifndef NDEBUG
    double energy2 = mol2othersBasic(mol);
    assert(energy+0.0000001 >= energy2 && energy-0.0000001 <= energy2);
#endif

    return energy;
}



double TotalEnergyCalculator::mol2others(Molecule &mol, vector<double> *changes) {
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

            //add interaction with external potential
            if (topo.exter.exist)
                energy += extere2(mol[j]);
        }

    } else {

        changes->resize( mol.size() * conf->pvec.size() );
        for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
            // for cycles => pair potential with all particles except those in molecule
            for (i = 0; i < mol[0]; i++) { // pair potential with all particles from 0 to the first one in molecule
                (*changes)[j*conf->pvec.size() + i] = pairE(&conf->pvec[mol[j]], &conf->pvec[i]);
                energy += (*changes)[j*conf->pvec.size() + i];
            }

            for (long i = mol[mol.size()-1] + 1; i < (long)conf->pvec.size(); i++) {
                (*changes)[j*conf->pvec.size() + i] = pairE(&conf->pvec[mol[j]], &conf->pvec[i]);
                energy += (*changes)[j*conf->pvec.size() + i];
            }

            //add interaction with external potential
            if (topo.exter.exist)
                energy += extere2(mol[j]);
        }
    }

#ifndef NDEBUG
    double energy2 = mol2othersBasic(mol);
    assert(energy+0.0000001 >= energy2 && energy-0.0000001 <= energy2);
#endif

    return energy;
}

double TotalEnergyCalculator::mol2othersBasic(Molecule &mol) {
    double energy=0.0;
    long i = 0;

    for(unsigned int j=0; j<mol.size(); j++) { // for all particles in molecule
        // for cycles => pair potential with all particles except those in molecule
        for (i = 0; i < mol[0]; i++) // pair potential with all particles from 0 to the first one in molecule
            energy += pairE(&conf->pvec[mol[j]], &conf->pvec[i]);

        for (i = mol[mol.size()-1] + 1; i < (long)conf->pvec.size(); i++)
            energy += pairE(&conf->pvec[mol[j]], &conf->pvec[i]);

        //add interaction with external potential
        if (topo.exter.exist)
            energy += extere2(mol[j]);
    }

    return energy;
}


////////////////////////////////////////////////////////////////////////////////////////////////


double TotalEnergyCalculator::mol2others(vector<Particle>& mol) {
    double energy=0.0;
    unsigned int i=0;

    //double test, test2;

    for(i=0; i<mol.size(); i++) {
        if (topo.exter.exist)
            energy += exterE.extere2(&mol[i]);

        for(unsigned int j=0; j < conf->pvec.size(); j++) {
            /*test=pairE(&mol[i], &conf->pvec[j]);
            test2 = pairE2(&mol[i], &conf->pvec[j]);
            if( test < 10000 && fabs( test - test2 ) > 0.00000001 ) {
                cout << std::setprecision(20) << fabs( test - test2 ) << endl;
                cout << std::setprecision(20) << test << " != " << test2 << endl;
            }*/
            energy += pairE(&mol[i], &conf->pvec[j]);
        }
    }

    return energy;
}


double TotalEnergyCalculator::chainInner(vector<Particle >& chain) {
    double energy = 0.0;

    //generate conlists
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


double TotalEnergyCalculator::chainInner(Molecule &chain) {
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


double TotalEnergyCalculator::oneToAll(Particle *target, ConList* conlist, Neighbors* neighborList) {
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
                /*std::cout.precision(15);
                cout << pairE[getThreadNum()](target, &conf->pvec[i], conlist) << endl;*/
            }
        }
    }

    //add interaction with external potential
    if (topo.exter.exist)
        energy += exterE.extere2(target);

    return energy;
}



/********************************************************************************/
/*                              LEGACY FUNCTIONS                                */
/********************************************************************************/


double TotalEnergyCalculator::operator ()(int target, int mode, int chainnum) {

    //DEBUG_SIM("Calculate the energy with mode %d", mode)

    if(mode == 1)
        return oneToAll(target);

    if(mode == 2)
        return chainToAll(target, chainnum);

    if(mode == 0)
        return allToAll();

    fprintf(stderr, "ERROR: Wrong mode (%d) was given to calc_energy!", mode);
    return 0.0;
}


double TotalEnergyCalculator::chainToAll(int target, int chainnum) {

    double energy=0.0;
    long i,j=0;
    ConList conlist = conf->pvec.getConlist(target);

    for (i = 0; i < target; i++) {
        if (i != conf->pvec.getChainPart(chainnum, j)) {
            energy+= pairE(&conf->pvec[target], &conf->pvec[i], &conlist);
        } else {
            j++;
        }
    }
    j++;

    for (i = target + 1; i < (long)conf->pvec.size(); i++) {
        if (i != conf->pvec.getChainPart(chainnum, j)) {
            energy+= pairE(&conf->pvec[target], &conf->pvec[i],  &conlist);
        } else {
            j++;
        }
    }

    //add interaction with external potential
    if (topo.exter.exist)
        energy+= extere2(target);

    return energy;
}
