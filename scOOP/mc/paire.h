/** @file paire.h*/

#ifndef PAIRE_H
#define PAIRE_H

#include "../structures/topo.h"
#include "../structures/Conf.h"

class PairE {
private:
    Particle* part1;
    Particle* part2;

    ConList* conlist;

    /*
     *  Parameters of interaction calculation
     */
    double dist;                    // closest distance
    Vector distvec;                 // vector of closes distance
    Vector r_cm;                    // Vector between centres of mass from part2 to part1
    double distcm;                  // distance between center of masses
    double dotrcm;                  // square size of r_cm
    double contt;                   // closest point on spherocylinder to sphere

    GeoBase* pbc;                   // box size

    double (PairE::*intFCE[MAXT][MAXT])();


public:
    PairE(GeoBase* pbc) : pbc(pbc) {}

    double operator() (Particle* part1, Particle* part2, ConList* conlist=NULL);

    /**
     * @brief init_intfce Initializes the array with the pointers to the energy function
     */
    void initIntFCE();

private:

    /**
     * @brief enoexist    Indicates not yet programmed interaction
     * @return
     */
    double eNoExist() {
        fprintf (stderr, "ERROR: We have not programed interaction of types %d and %d\n", part1->type,part2->type);
        return 0.0;
    }
};

#endif // PAIRE_H
