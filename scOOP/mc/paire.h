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
    PairE(GeoBase* pbc) : pbc(pbc) {
        initIntFCE();
    }

    double operator() (Particle* part1, Particle* part2, ConList* conlist=NULL) {
        double energy=0.0;

        this->part1 = part1;
        this->part2 = part2;
        this->conlist = conlist;

        // Placing interactin particle in unit pbc->box and finding vector connecting CM
        r_cm = pbc->image(&part1->pos, &part2->pos);

        dotrcm = DOT(r_cm,r_cm);

        if (dotrcm > topo.sqmaxcut) return 0.0;  // distance so far that even spherocylinders cannot be within cutoff

        /*contt = 0; // not used in SPA-SPA interaction
        distvec.x = 0;
        distvec.y = 0;
        distvec.z = 0;*/

        if(intFCE[part1->type][part2->type] == NULL)
            fprintf(stderr, "interaction function for type %d and %d not defined!\n", part1->type, part2->type);

        energy = (this->*intFCE[part1->type][part2->type])();

        //printf("num: %ld  %ld   e: %f dist: %f",num1,num2,energy,interact.dist);

        /*if(conlist != NULL) {
            energy += bondEnergy ();
            energy += angleEnergy ();
        }*/

        return energy;
    }

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

    inline double e2ScaOr2Spa() {
        double repenergy=0.0, atrenergy=0.0;

        dist = sqrt(dotrcm);
        distcm = dist;

        repenergy = eRepulsive();

        if ( ( dotrcm > topo.ia_params[part1->type][part2->type].rcutSq ) ||
             (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
             topo.ia_params[part1->type][part2->type].exclude ) {
            atrenergy = 0.0;
        } else {
            if (dotrcm < topo.ia_params[part1->type][part2->type].pdisSq)
                atrenergy = -topo.ia_params[part1->type][part2->type].epsilon;
            else  {
                atrenergy = cos(PIH*(dist - topo.ia_params[part1->type][part2->type].pdis) * topo.ia_params[part1->type][part2->type].pswitchINV);
                atrenergy *= -atrenergy*topo.ia_params[part1->type][part2->type].epsilon ;
            }
        }

    #ifdef LJ    // LJ
        double en6 = pow((topo.ia_params[part1->type][part2->type].sigma / dist),6);
        repenergy = 4*en6*(en6-1);
        atrenergy = 0.0;
    #endif

        return repenergy+atrenergy;
    }

    inline double eRepulsive() {
        // WCA repulsion
        if (dotrcm > topo.ia_params[part1->type][part2->type].rcutwcaSq) {
            return 0.0;
        } else {
            return topo.ia_params[part1->type][part2->type].A * pow(dotrcm, -6) - topo.ia_params[part1->type][part2->type].B * pow(dotrcm, -3) + topo.ia_params[part1->type][part2->type].epsilon;
        }
    }
};

#endif // PAIRE_H
