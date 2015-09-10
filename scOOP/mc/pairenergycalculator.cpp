#include "pairenergycalculator.h"

#include "math_calc.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

extern Topo topo;

double PairEnergyCalculator::operator ()(Particle *part1, Particle *part2, ConList* conlist) {

    double energy=0.0;

#ifdef OMP1 // switching of particles data occurs durin energy calc -> cant use with openMP
    Particle par1 = *part1;
    Particle par2 = *part2;
    this->part1 = &par1;
    this->part2 = &par2;
#else // advantageous performace wise for single core simulation
    this->part1 = part1;
    this->part2 = part2;
#endif

    this->conlist = conlist;

    /*Placing interactin particle in unit pbc->box and finding vector connecting CM*/
    r_cm = pbc->image(&part1->pos, &part2->pos); // explicit statement below for performance optimization*/
/*
    r_cm.x = part1->pos.x - part2->pos.x;
    r_cm.y = part1->pos.y - part2->pos.y;
    r_cm.z = part1->pos.z - part2->pos.z;

    if ( r_cm.x < 0  )
        r_cm.x = pbc->box.x * (r_cm.x - (double)( (long)(r_cm.x-0.5) ) );
    else
        r_cm.x = pbc->box.x * (r_cm.x - (double)( (long)(r_cm.x+0.5) ) );
    if ( r_cm.y < 0  )
        r_cm.y = pbc->box.y * (r_cm.y - (double)( (long)(r_cm.y-0.5) ) );
    else
        r_cm.y = pbc->box.y * (r_cm.y - (double)( (long)(r_cm.y+0.5) ) );
    if ( r_cm.z < 0  )
        r_cm.z = pbc->box.z * (r_cm.z - (double)( (long)(r_cm.z-0.5) ) );
    else
        r_cm.z = pbc->box.z * (r_cm.z - (double)( (long)(r_cm.z+0.5) ) );
*/

    dotrcm = DOT(r_cm,r_cm);

    // HARDSPHERE
    //if(dotrcm >= 1.0) return 0.0;
    //if(dotrcm < 1.0) return INFINITY;

    if (dotrcm > topo.sqmaxcut) return 0.0;  /* distance so far that even spherocylinders cannot be within cutoff  */

    contt = 0;
    distvec.x = 0;
    distvec.y = 0;
    distvec.z = 0;

    if(intFCE[part1->type][part2->type] == NULL)
        fprintf(stderr, "interaction function for type %d and %d not defined!\n", part1->type, part2->type);

    energy = (this->*intFCE[part1->type][part2->type])();

    //printf("num: %ld  %ld   e: %f dist: %f",num1,num2,energy,interact.dist);

    if(conlist != NULL) {
        energy += bondEnergy ();
        energy += angleEnergy ();
    }

    //printf("  e: %f\n",energy);
    return energy;
}

void PairEnergyCalculator::initIntFCE() {
    // NB
    // Fill in the names of the functions for calculating the
    // interaction energy
    long geotype, other_geotype;
    int i, j;
    for(i = 0; i < MAXT; i++){
        for(j = 0; j < MAXT; j++){
            /* Initialize them as not existing */
            intFCE[i][j] = &PairEnergyCalculator::eNoExist;
            geotype = topo.ia_params[i][j].geotype[0];
            other_geotype = topo.ia_params[i][j].geotype[1];

            if ( ( (geotype == CHCPSC || geotype == CPSC || geotype == TCHCPSC || geotype == TCPSC) &&
                    (other_geotype == CHPSC || other_geotype == PSC || other_geotype == TCHPSC || other_geotype == TPSC) ) ||
                  ( (geotype == CHPSC || geotype == PSC || geotype == TCHPSC || geotype == TPSC)  &&
                    (other_geotype == CHCPSC || other_geotype == CPSC || other_geotype == TCHCPSC || other_geotype == TCPSC) ) )  {
                intFCE[i][j] = &PairEnergyCalculator::ePscPsc;
            }
            if ( (geotype == CHCPSC || geotype == CPSC || geotype == TCHCPSC || geotype == TCPSC) &&
                    (other_geotype == CHCPSC || other_geotype == CPSC || other_geotype == TCHCPSC || other_geotype == TCPSC) ){
                intFCE[i][j] = &PairEnergyCalculator::eCpscCpsc;
            }
            if ( (geotype == CHPSC || geotype == PSC || geotype == TCHPSC || geotype == TPSC) &&
                    (other_geotype == CHPSC || other_geotype == PSC || other_geotype == TCHPSC || other_geotype == TPSC) ){
                intFCE[i][j] = &PairEnergyCalculator::ePscPsc;
            }

            if(geotype == SCN || geotype == SPN
                    || other_geotype == SCN || other_geotype == SPN){
                intFCE[i][j] = &PairEnergyCalculator::eSpnOrScn;
            }
            if((geotype == SCA && other_geotype == SCA)
                    || (geotype == SPA && other_geotype == SPA)){
                intFCE[i][j] = &PairEnergyCalculator::e2ScaOr2Spa;
            }
            if((geotype == SCA && other_geotype == SPA)
                    || (geotype == SPA && other_geotype == SCA)){
                intFCE[i][j] = &PairEnergyCalculator::eSpaSca;
            }
            if(( (geotype == PSC || geotype == CHPSC || geotype == TCHPSC || geotype == TPSC) && other_geotype == SPA)
                    || (geotype == SPA && (other_geotype == PSC||other_geotype == CHPSC || other_geotype == TCHPSC || other_geotype == TPSC) )){
                intFCE[i][j] = &PairEnergyCalculator::ePscSpa;
            }
            if(( (geotype == CPSC ||geotype == CHCPSC || geotype == TCHCPSC || geotype == TCPSC) && other_geotype == SPA)
                    || (geotype == SPA && (other_geotype == CPSC||other_geotype == CHCPSC || other_geotype == TCHCPSC || other_geotype == TCPSC)  )){
                intFCE[i][j] = &PairEnergyCalculator::eCpscSpa;
            }

        }
    }
}


double PairEnergyCalculator::angleEnergy() {
    double energy=0.0, currangle, halfl;
    Vector vec1, vec2;
    int * geotype = topo.ia_params[part1->type][part2->type].geotype;

    /*angle interaction with nearest neighbours -harmonic*/
    if ((topo.moleculeParam[part1->molType]).angle1c >= 0) {
        if (part2 == conlist->conlist[0]) {
            /*num1 is connected to num2 by tail*/
            if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                /*spheres do not have this interaction*/
                energy += 0.0;
            else {
                if (geotype[0] < SP)
                    vec1 = part1->dir;
                else {
                    halfl = topo.ia_params[part1->type][part2->type].half_len[1];
                    //sphere angle is defined versus the end of spherocylind.er
                    vec1.x = part2->pos.x - part2->dir.x * halfl / pbc->box.x;
                    vec1.y = part2->pos.y - part2->dir.y * halfl / pbc->box.y;
                    vec1.z = part2->pos.z - part2->dir.z * halfl / pbc->box.z;
                    vec1 = pbc->image(&vec1, &part1->pos);
                }
                if (geotype[1] < SP)
                    vec2 = part2->dir;
                else {
                    halfl = topo.ia_params[part1->type][part2->type].half_len[0];
                    vec2.x = part1->pos.x + part1->dir.x * halfl / pbc->box.x;
                    vec2.y = part1->pos.y + part1->dir.y * halfl / pbc->box.y;
                    vec2.z = part1->pos.z + part1->dir.z * halfl / pbc->box.z;
                    vec2 = pbc->image(&vec2, &part2->pos);
                }
                vec1.normalise();
                vec2.normalise();
                currangle = acos(DOT(vec1,vec2));
                energy += harmonicPotential(currangle,topo.moleculeParam[part1->molType].angle1eq,topo.moleculeParam[part1->molType].angle1c);
            }
        } else {
            if (part2 == conlist->conlist[1]) {
                /*num1 is connected to num2 by head*/
                if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                /*spheres do not have this interaction*/
                    energy += 0.0;
                else {
                    if (geotype[0] < SP)
                        vec1 = part1->dir;
                    else {
                        halfl = topo.ia_params[part1->type][part2->type].half_len[1];
                        //sphere angle is defined versus the end of spherocylinder
                        vec1.x = part2->pos.x + part2->dir.x * halfl / pbc->box.x;
                        vec1.y = part2->pos.y + part2->dir.y * halfl / pbc->box.y;
                        vec1.z = part2->pos.z + part2->dir.z * halfl / pbc->box.z;
                        vec1 = pbc->image(&vec1, &part1->pos);
                    }
                    if (geotype[1] < SP)
                        vec2 = part2->dir;
                    else {
                        halfl = topo.ia_params[part1->type][part2->type].half_len[0];
                        vec2.x = part1->pos.x - part1->dir.x * halfl / pbc->box.x;
                        vec2.y = part1->pos.y - part1->dir.y * halfl / pbc->box.y;
                        vec2.z = part1->pos.z - part1->dir.z * halfl / pbc->box.z;
                        vec2 = pbc->image(&vec2, &part2->pos);
                    }
                    vec1.normalise();
                    vec2.normalise();
                    currangle = acos(DOT(vec1,vec2));
                    energy += harmonicPotential(currangle,topo.moleculeParam[part2->molType].angle1eq,topo.moleculeParam[part2->molType].angle1c);
                }
            }
        }
    }

    /*interaction between the orientation of  spherocylinders patches -harmonic*/
    if (topo.moleculeParam[part1->molType].angle2c >= 0) {
        if (part2 == conlist->conlist[0]) {
            /*num1 is connected to num2 by tail*/
            if ( (geotype[0] < SP) && (geotype[1] < SP) ) {
                currangle = acos(DOT(part1->patchdir[0],part2->patchdir[0]) - DOT(part1->dir,part2->patchdir[0])  );
                energy += harmonicPotential(currangle,topo.moleculeParam[part1->molType].angle2eq,topo.moleculeParam[part1->molType].angle2c);
            } else {
                energy += 0.0;
            }
        } else {
            if (part2 == conlist->conlist[1]) {
                /*num1 is connected to num2 by head*/
                if ( (geotype[0] < SP) && (geotype[1] < SP) ) {
                    currangle = acos(DOT(part2->patchdir[0],part1->patchdir[0]) - DOT(part2->dir,part1->patchdir[0])  );
                    energy += harmonicPotential(currangle,topo.moleculeParam[part2->molType].angle2eq,topo.moleculeParam[part2->molType].angle2c);
                } else {
                    energy += 0.0;
                }
            }
        }
    }
    //    printf("angleener: %f\n",energy);
    return energy;
}

double PairEnergyCalculator::bondEnergy() {
    double energy=0.0, bondlength, halfl;
    Vector vec1, vec2, vecbond;
    int * geotype = topo.ia_params[part1->type][part2->type].geotype;

    /*interaction with nearest neighbours -harmonic*/
    if ((topo.moleculeParam[part1->molType]).bond1c >= 0) {

        if (part2 == conlist->conlist[1]) {
            /*num1 is connected to num2 by tail*/
            if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                energy = harmonicPotential(distcm,topo.moleculeParam[part1->molType].bond1eq,topo.moleculeParam[part1->molType].bond1c);
            else {
                if (geotype[0] < SP)
                    halfl =topo.ia_params[part1->type][part2->type].half_len[0];
                else
                    halfl = 0.0;

                vec1.x = part1->pos.x - part1->dir.x * halfl / pbc->box.x;
                vec1.y = part1->pos.y - part1->dir.y * halfl / pbc->box.y;
                vec1.z = part1->pos.z - part1->dir.z * halfl / pbc->box.z;

                if (geotype[1] < SP)
                    halfl =topo.ia_params[part1->type][part2->type].half_len[1];
                else
                    halfl = 0.0;

                vec2.x = part2->pos.x + part2->dir.x * halfl / pbc->box.x;
                vec2.y = part2->pos.y + part2->dir.y * halfl / pbc->box.y;
                vec2.z = part2->pos.z + part2->dir.z * halfl / pbc->box.z;
                vecbond = pbc->image(&vec1, &vec2);
                bondlength = sqrt(DOT(vecbond,vecbond));
                energy = harmonicPotential(bondlength,topo.moleculeParam[part1->molType].bond1eq,topo.moleculeParam[part1->molType].bond1c);
            }
        } else {
            if (part2 == conlist->conlist[0]) {
                /*num1 is connected to num2 by head*/
                if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                    energy = harmonicPotential(distcm,topo.moleculeParam[part1->molType].bond1eq,topo.moleculeParam[part1->molType].bond1c);
                else {
                    if (geotype[0] < SP)
                        halfl =topo.ia_params[part1->type][part2->type].half_len[0];
                    else
                        halfl = 0.0;
                    vec1.x = part1->pos.x + part1->dir.x * halfl / pbc->box.x;
                    vec1.y = part1->pos.y + part1->dir.y * halfl / pbc->box.y;
                    vec1.z = part1->pos.z + part1->dir.z * halfl / pbc->box.z;

                    if (geotype[0] < SP)
                        halfl =topo.ia_params[part1->type][part2->type].half_len[0];
                    else
                        halfl = 0.0;

                    vec2.x = part2->pos.x - part2->dir.x * halfl / pbc->box.x;
                    vec2.y = part2->pos.y - part2->dir.y * halfl / pbc->box.y;
                    vec2.z = part2->pos.z - part2->dir.z * halfl / pbc->box.z;
                    vecbond = pbc->image(&vec1, &vec2);
                    bondlength = sqrt(DOT(vecbond,vecbond));
                    energy = harmonicPotential(bondlength,topo.moleculeParam[part1->molType].bond1eq,topo.moleculeParam[part1->molType].bond1c);
                }
            }
        }
    }
    /*interaction with second nearest neighbours -harmonic*/
    if (topo.moleculeParam[part1->molType].bond2c >= 0) {
        if (part2 == conlist->conlist[2]) {
            /*num1 is connected to num2 by tail*/
            if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                energy = harmonicPotential(distcm,topo.moleculeParam[part1->molType].bond2eq,topo.moleculeParam[part1->molType].bond2c);
            else {
                vecbond = pbc->image(&part1->pos, &part2->pos);
                bondlength = sqrt(DOT(vecbond,vecbond));
                energy = harmonicPotential(bondlength,topo.moleculeParam[part1->molType].bond2eq,topo.moleculeParam[part1->molType].bond2c);
            }
        } else {
            if (part2 == conlist->conlist[3]) {
                /*num1 is connected to num2 by head*/
                if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                    energy = harmonicPotential(distcm,topo.moleculeParam[part1->molType].bond2eq,topo.moleculeParam[part1->molType].bond2c);
                else {
                    vecbond = pbc->image(&part1->pos, &part2->pos);
                    bondlength = sqrt(DOT(vecbond,vecbond));
                    energy = harmonicPotential(bondlength,topo.moleculeParam[part1->molType].bond2eq,topo.moleculeParam[part1->molType].bond2c);
                }
            }
        }
    }
    /*interaction with nearest neighbours - direct harmonic bond*/
    if ((topo.moleculeParam[part1->molType]).bonddc > 0) {

        if (part2 == conlist->conlist[1]) {
            /*num1 is connected to num2 by tail*/
            if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                energy = harmonicPotential(distcm,topo.moleculeParam[part1->molType].bonddeq,topo.moleculeParam[part1->molType].bonddc);
            else {
                if (geotype[0] < SP)
                    halfl =topo.ia_params[part1->type][part2->type].half_len[0];
                else
                    halfl = 0.0;
                vec1.x = part1->pos.x - part1->dir.x * halfl / pbc->box.x;
                vec1.y = part1->pos.y - part1->dir.y * halfl / pbc->box.y;
                vec1.z = part1->pos.z - part1->dir.z * halfl / pbc->box.z;
                if (geotype[1] < SP)
                    halfl =topo.ia_params[part1->type][part2->type].half_len[1];
                else
                    halfl = 0.0;
                vec2.x = part2->pos.x + part2->dir.x * (halfl + topo.moleculeParam[part1->molType].bonddeq) / pbc->box.x ;
                vec2.y = part2->pos.y + part2->dir.y * (halfl + topo.moleculeParam[part1->molType].bonddeq) / pbc->box.y ;
                vec2.z = part2->pos.z + part2->dir.z * (halfl + topo.moleculeParam[part1->molType].bonddeq) / pbc->box.z ;
                vecbond = pbc->image(&vec1, &vec2);
                bondlength = sqrt(DOT(vecbond,vecbond));
                energy = harmonicPotential(bondlength,0.0,topo.moleculeParam[part1->molType].bonddc);
            }
        } else {
            if (part2 == conlist->conlist[0]) {
                /*num1 is connected to num2 by head*/
                if ( (geotype[0] >= SP) && (geotype[1] >= SP) )
                    energy = harmonicPotential(distcm,topo.moleculeParam[part1->molType].bond1eq,topo.moleculeParam[part1->molType].bond1c);
                else {
                    if (geotype[0] < SP)
                        halfl =topo.ia_params[part1->type][part2->type].half_len[0];
                    else
                        halfl = 0.0;
                    vec1.x = part1->pos.x + part1->dir.x * (halfl + topo.moleculeParam[part1->molType].bonddeq) / pbc->box.x ;
                    vec1.y = part1->pos.y + part1->dir.y * (halfl + topo.moleculeParam[part1->molType].bonddeq) / pbc->box.y ;
                    vec1.z = part1->pos.z + part1->dir.z * (halfl + topo.moleculeParam[part1->molType].bonddeq) / pbc->box.z ;
                    if (geotype[0] < SP)
                        halfl =topo.ia_params[part1->type][part2->type].half_len[0];
                    else
                        halfl = 0.0;
                    vec2.x = part2->pos.x - part2->dir.x * halfl / pbc->box.x;
                    vec2.y = part2->pos.y - part2->dir.y * halfl / pbc->box.y;
                    vec2.z = part2->pos.z - part2->dir.z * halfl / pbc->box.z;
                    vecbond = pbc->image(&vec1, &vec2);
                    bondlength = sqrt(DOT(vecbond,vecbond));
                    energy = harmonicPotential(bondlength,0.0,topo.moleculeParam[part1->molType].bonddc);
                }
            }
        }
    }
    //printf("bondlength: %f\n",bondlength);
    //    printf("bondener: %f\n",energy);
    return energy;
}

double PairEnergyCalculator::eNoExist() {
    double energy=0.0;

    fprintf (stderr, "ERROR: We have not programed interaction of types %d and %d\n", part1->type,part2->type);
    //exit (1);

    return energy;
}

double PairEnergyCalculator::ePscPsc() {
    double atrenergy, repenergy;

    closestDist();
    repenergy = eRepulsive();
    if ( ( dist >topo.ia_params[part1->type][part2->type].rcut ) ||
         (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
         topo.ia_params[part1->type][part2->type].exclude ) {
        atrenergy = 0.0;
    } else {
        bool firstCH=false, secondCH=false;
        Vector olddir1 = part1->dir;
        Vector olddir2 = part2->dir;
        if ( (topo.ia_params[part1->type][part2->type].geotype[0] == CHPSC)||(topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) )
            firstCH = true;
        if ( (topo.ia_params[part1->type][part2->type].geotype[1] == CHPSC)||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHPSC) )
            secondCH = true;
        if (firstCH)
            part1->dir = part1->chdir[0];
        if (secondCH)
            part2->dir = part2->chdir[0];

        if ((firstCH) || (secondCH) ) {
            closestDist();
        }
        atrenergy = eattractivePscPsc(0,0);

        /*addition of interaction of second patches*/
        if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) ||
          (topo.ia_params[part1->type][part2->type].geotype[1] == TPSC) ||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHPSC) ) {
            bool firstT=false, secondT=false;
            if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) )
            firstT = true;
            if ( (topo.ia_params[part1->type][part2->type].geotype[1] == TPSC) ||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHPSC)  )
            secondT = true;

            if (firstT) {
            if (firstCH) {
                part1->dir = part1->chdir[1];
                closestDist();
            }
            atrenergy += eattractivePscPsc(1,0);
            }
            if ( (firstT) && (secondT) ) {
            if (secondCH) {
                part2->dir = part2->chdir[1];
                closestDist();
            }
            atrenergy += eattractivePscPsc(1,1);
            }
            if (secondT) {
            if (firstT && firstCH ) {
                part1->dir = part1->chdir[0];
                closestDist();
            }
            atrenergy += eattractivePscPsc(0,1);
            }
        }

        if (firstCH)
            part1->dir = olddir1;
        if (secondCH)
            part2->dir = olddir2;

    }
    return repenergy+atrenergy;
}

double PairEnergyCalculator::eattractivePscPsc(int patchnum1, int patchnum2) {
    int i, intrs;
    double atrenergy, ndist;
    double v1, v2, f0, f1, f2, T1, T2, S1, S2, a, paral;
    double intersections[5];
    Vector vec1, vec2, vec_intrs, vec_mindist;

   
    //interact->halfl =topo.ia_params[part1->type][part2->type].half_len[0];
    //DEBUG_SIM("halfl = %lf", interact->halfl);
    for(i=0 ;i<5; i++)
        intersections[i]=0;
    //cospatch =topo.ia_params.pcanglsw;
    //cospatchinr =topo.ia_params.pcangl;
    /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at.
      cut distance C*/
    //DEBUG_SIM("first intersection");
    intrs = pscIntersect(part1,part2,topo.ia_params[part1->type][part2->type].half_len[0],topo.ia_params[part1->type][part2->type].half_len[1], r_cm, intersections, 0, patchnum1);
    if (intrs <2){
        //DEBUG_SIM("No intersection :(");
        return 0.0; /*sc is all outside patch, attractive energy is 0*/
    }

    T1 = intersections[0]; /*points on sc2*/
    T2 = intersections[1];
    /*2- now do the same oposite way psc1 in patch of psc2*/
    for(i=0;i<5;i++)
        intersections[i]=0;
    //DEBUG_SIM("get vector");
    vec1 = -1.0*r_cm;
    //DEBUG_SIM("second intersection");
    intrs = pscIntersect(part2, part1,topo.ia_params[part1->type][part2->type].half_len[1],topo.ia_params[part1->type][part2->type].half_len[0], vec1, intersections, 1, patchnum2);
    if (intrs <2)
        return 0.0; /*sc is all outside patch, attractive energy is 0*/

    S1 = intersections[0]; /*points on sc1*/
    S2 = intersections[1];

    /*3- scaling function1: dependence on the length of intersetions*/
    // BUG FIX 23.9.2014
    v1=fabs(S1-S2); // v1=0.5*fabs(S1-S2);
    v2=fabs(T1-T2); // v2=0.5*fabs(T1-T2);
    f0=0.5*(v1+v2);

    /*4a- with two intersection pices calculate vector between their CM
      -this is for angular orientation*/
    vec1 = ((S1+S2)*0.5) * part1->dir;
    vec2 = ((T1+T2)*0.5) * part2->dir;
    vec_intrs.x = vec2.x - vec1.x - r_cm.x;
    vec_intrs.y = vec2.y - vec1.y - r_cm.y;
    vec_intrs.z = vec2.z - vec1.z - r_cm.z;
    /*vec_intrs should be from sc1 to sc2*/
    //fprintf (stderr, "segments_CM: %.8f %.8f %.8f \n",vec_intrs.x,vec_intrs.y,vec_intrs.z);

    /*4b - calculate closest distance attractive energy from it*/
    vec_mindist = minDistSegments(v1,v2,vec_intrs);
    //fprintf (stderr, "segments closest dist: %.8f %.8f %.8f \n",vec_mindist.x,vec_mindist.y,vec_mindist.z);
    ndist=sqrt(DOT(vec_mindist,vec_mindist));
    //dist=DOT(vec_intrs,vec_intrs);
    if (ndist <topo.ia_params[part1->type][part2->type].pdis)
        atrenergy = -topo.ia_params[part1->type][part2->type].epsilon;
    //atrenergy = -1.0;
    else {
        atrenergy = cos(PIH*(ndist-topo.ia_params[part1->type][part2->type].pdis)/topo.ia_params[part1->type][part2->type].pswitch);
        atrenergy *= -atrenergy*topo.ia_params[part1->type][part2->type].epsilon;
    }
    /*5- scaling function2: angular dependence of patch1*/
    vec1 = 1.0 * vec_intrs;
    //vec1=vecScale(vec_mindist,-1.0);
    vec1=vecPerpProject(&vec1, &part1->dir);
    vec1.normalise();
    a = DOT(vec1,part1->patchdir[patchnum1]);
    f1 = fanglScale(a, 0+2*patchnum1);

    /*6- scaling function3: angular dependence of patch2*/
    vec1 = -1.0 * vec_intrs;
    //vec1=vecScale(vec_mindist,1.0);
    vec1=vecPerpProject(&vec1, &part2->dir);
    vec1.normalise();
    a = DOT(vec1,part2->patchdir[patchnum2]);
    f2 = fanglScale(a, 1+2*patchnum2);
    //printf("v1: %f v2: %f f0: %f f1: %f f2: %f ener: %f\n",v1,v2,f0,f1,f2,atrenergy);

    //add scaling increased if particles are parallel or antiparallel
    paral=1.0;
    paral = scparallel(topo.ia_params[part1->type][part2->type].parallel,part1->dir,part2->dir);
    
    /*7- put it all together*/
    atrenergy *=f0*f1*f2*paral;
    //if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
    //	    fprintf (stderr, "attraction  %.8f \n",atrenergy);
    //	    exit(1);

    return atrenergy;
}

double PairEnergyCalculator::eCpscCpsc() {
    double atrenergy, repenergy;

    //DEBUG_SIM("do energy 33") ;
    closestDist();
    repenergy = eRepulsive();
    //DEBUG_SIM("got the rep. energy");
    if ( ( dist > topo.ia_params[part1->type][part2->type].rcut ) ||
         ( topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
         topo.ia_params[part1->type][part2->type].exclude ) {
        atrenergy = 0.0;
    } else {
        bool firstCH=false, secondCH=false;
        Vector olddir1 = part1->dir;
        Vector olddir2 = part2->dir;
        if ( (topo.ia_params[part1->type][part2->type].geotype[0] == CHCPSC)||(topo.ia_params[part1->type][part2->type].geotype[0] == TCHCPSC) )
            firstCH = true;
        if ( (topo.ia_params[part1->type][part2->type].geotype[1] == CHCPSC)||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHCPSC) )
            secondCH = true;
        if(firstCH)
            part1->dir = part1->chdir[0];
        if(secondCH)
            part2->dir = part2->chdir[0];

        if ((firstCH) || (secondCH) ) {
            closestDist();
        }
        atrenergy = eattractiveCpscCpsc(0,0);

        /*addition of interaction of second patches*/
        if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TCPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHCPSC) ||
          (topo.ia_params[part1->type][part2->type].geotype[1] == TCPSC) ||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHCPSC) ) {
            bool firstT=false, secondT=false;
            if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TCPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHCPSC) )
            firstT = true;
            if ( (topo.ia_params[part1->type][part2->type].geotype[1] == TCPSC) ||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHCPSC)  )
            secondT = true;

            if (firstT) {
            if (firstCH) {
                part1->dir = part1->chdir[1];
                closestDist();
            }
            atrenergy += eattractiveCpscCpsc(1,0);
            }
            if ( (firstT) && (secondT) ) {
            if (secondCH) {
                part2->dir = part2->chdir[1];
                closestDist();
            }
            atrenergy += eattractiveCpscCpsc(1,1);
            }
            if (secondT) {
            if (firstT && firstCH ) {
                part1->dir = part1->chdir[0];
                closestDist();
            }
            atrenergy += eattractiveCpscCpsc(0,1);

            }
        }

        if (firstCH)
            part1->dir = olddir1;
        if (secondCH)
            part2->dir = olddir2;

#ifdef EXTRA_HYDROPHOBIC_ALL_BODY_ATTRACTION
        double extraAttr = 0;
        double extraEpsilon = -1.0 * E_ISO; //  -> isotropic hydrophobic interaction Epsilon

        double extraInteractionSwitch = 0.2; //  -> isotropic hydrophobic interaction Switch

        double distRcm = sqrt(dotrcm);

        if (distRcm > topo.ia_params[part1->type][part2->type].pdis+extraInteractionSwitch) {
            extraAttr = 0.0;
        } else {
            if (distRcm < topo.ia_params[part1->type][part2->type].pdis)
                extraAttr = extraEpsilon;
            else {
                extraAttr = cos(PIH*(distRcm - topo.ia_params[part1->type][part2->type].pdis)/extraInteractionSwitch);
                extraAttr *= extraAttr * extraEpsilon ;
            }
            // cos between two SC
            extraAttr *= (part1->dir.dot(part2->dir));
        }
        atrenergy += extraAttr;
#endif

    }
    return repenergy+atrenergy;
}

double PairEnergyCalculator::eattractiveCpscCpsc(int patchnum1, int patchnum2) {
    int i, intrs;
    double atrenergy, v1, v2, f0, f1, f2, T1, T2, S1, S2, a, paral, ndist;
    double intersections[5];
    Vector vec1, vec2, vec_intrs, vec_mindist;

    //rcut =topo.ia_params[part1->type][part2->type].rcut; // dont use, just a reminder

//	interact->halfl =topo.ia_params[part1->type][part2->type].half_len[0];
    for(i=0;i<5;i++)
        intersections[i]=0;
    /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at.
      cut distance C*/
    intrs = cpscIntersect(part1,part2,topo.ia_params[part1->type][part2->type].half_len[0],
            topo.ia_params[part1->type][part2->type].half_len[1], r_cm, intersections, 0, patchnum1);
    if (intrs <2)
        return 0.0; /*sc is all outside patch, attractive energy is 0*/
    T1 = intersections[0]; /*points on sc2*/
    T2 = intersections[1];
    /*2- now do the same oposite way psc1 in patch of psc2*/
    for(i=0;i<5;i++)
        intersections[i]=0;
    vec1 = -1.0 * r_cm;
    intrs = cpscIntersect(part2,part1,topo.ia_params[part1->type][part2->type].half_len[1],
            topo.ia_params[part1->type][part2->type].half_len[0],vec1, intersections, 1, patchnum2);
    if (intrs <2)
        return 0.0; /*sc is all outside patch, attractive energy is 0*/
    S1=intersections[0]; /*points on sc1*/
    S2=intersections[1];

    /*3- scaling function1: dependence on the length of intersetions*/
    // BUG fix 23.9.2014
    v1=fabs(S1-S2); // v1=0.5*fabs(S1-S2);
    v2=fabs(T1-T2); // v2=0.5*fabs(T1-T2);
    f0=0.5*(v1+v2);

    /*4a- with two intersection pices calculate vector between their CM
      -this is for angular orientation*/
    vec1 = ((S1+S2)*0.5) * part1->dir;
    vec2 = ((T1+T2)*0.5) * part2->dir;
    vec_intrs.x =vec2.x - vec1.x - r_cm.x;
    vec_intrs.y =vec2.y - vec1.y - r_cm.y;
    vec_intrs.z =vec2.z - vec1.z - r_cm.z;
    /*vec_intrs should be from sc1 to sc2*/
    //	    fprintf (stderr, "segments_CM: %.8f %.8f %.8f \n",vec_intrs.x,vec_intrs.y,vec_intrs.z);

    /*4b - calculate closest distance attractive energy from it*/

    vec_mindist = minDistSegments(v1,v2,vec_intrs);
    ndist=sqrt(DOT(vec_mindist,vec_mindist));

    //dist=DOT(vec_intrs,vec_intrs);

    if (ndist <topo.ia_params[part1->type][part2->type].pdis)
        atrenergy = -topo.ia_params[part1->type][part2->type].epsilon;
    else {
        atrenergy = cos(PIH*(ndist-topo.ia_params[part1->type][part2->type].pdis)/topo.ia_params[part1->type][part2->type].pswitch);
        atrenergy *= -atrenergy*topo.ia_params[part1->type][part2->type].epsilon ;
    }

    /*5- scaling function2: angular dependence of patch1*/
    vec1 = 1.0 * vec_intrs;
    //vec1=vecScale(vec_mindist,-1.0);
    vec1=vecPerpProject(&vec1, &part1->dir);
    vec1.normalise();
    a = DOT(vec1,part1->patchdir[patchnum1]);
    f1 = fanglScale(a, 0+2*patchnum1);

    /*6- scaling function3: angular dependence of patch2*/
    vec1 = -1.0 * vec_intrs;
    //vec1=vecScale(vec_mindist,1.0);
    vec1=vecPerpProject(&vec1, &part2->dir);
    vec1.normalise();
    a = DOT(vec1, part2->patchdir[patchnum2]);
    f2 = fanglScale(a, 1+2*patchnum2);

    //add scaling increased if particles are parallel or antiparallel
    paral=1.0;
    paral = scparallel(topo.ia_params[part1->type][part2->type].parallel,part1->dir,part2->dir);
    
    /*7- put it all together*/
    atrenergy *=f0*f1*f2*paral;

    //if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
    //	    fprintf (stderr, "attraction  %.8f \n",atrenergy);
    //	    exit(1);

    return atrenergy;
}

double PairEnergyCalculator::ePscCpsc() {
    double atrenergy, repenergy;

    //DEBUG_SIM("do energy 23") ;
    closestDist();
    repenergy = eRepulsive();
    //DEBUG_SIM("got the rep. energy");
    if ( ( dist >topo.ia_params[part1->type][part2->type].rcut ) ||
         (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
         topo.ia_params[part1->type][part2->type].exclude ) {
        atrenergy = 0.0;
    } else {
        bool firstCH=false, secondCH=false;
        Vector olddir1 = part1->dir;
        Vector olddir2 = part2->dir;
        if ((topo.ia_params[part1->type][part2->type].geotype[0] == CHPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == CHCPSC)||
          (topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHCPSC) )
            firstCH = true;
        if ((topo.ia_params[part1->type][part2->type].geotype[1] == CHPSC) || (topo.ia_params[part1->type][part2->type].geotype[1] == CHCPSC)||
          (topo.ia_params[part1->type][part2->type].geotype[1] == TCHPSC) || (topo.ia_params[part1->type][part2->type].geotype[1] == TCHCPSC) )
            secondCH = true;
        if(firstCH)
            part1->dir = part1->chdir[0];
        if(secondCH)
            part2->dir = part2->chdir[0];

        if ((firstCH) || (secondCH) ) {
            closestDist();
        }
        atrenergy = eattractivePscCpsc(0,0);


        /*addition of interaction of second patches*/
        if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TCPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHCPSC) ||
              (topo.ia_params[part1->type][part2->type].geotype[0] == TPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) ||
              (topo.ia_params[part1->type][part2->type].geotype[1] == TCPSC) ||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHCPSC) ||
              (topo.ia_params[part1->type][part2->type].geotype[1] == TPSC) || (topo.ia_params[part1->type][part2->type].geotype[1] == TCHPSC) ) {
            bool firstT=false, secondT=false;
            if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TCPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHCPSC) ||
                (topo.ia_params[part1->type][part2->type].geotype[0] == TPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) )
            firstT = true;
            if ( (topo.ia_params[part1->type][part2->type].geotype[1] == TCPSC) || (topo.ia_params[part1->type][part2->type].geotype[1] == TCHCPSC) ||
                (topo.ia_params[part1->type][part2->type].geotype[1] == TPSC) || (topo.ia_params[part1->type][part2->type].geotype[1] == TCHPSC) )
            secondT = true;

            if (firstT) {
            if (firstCH) {
                part1->dir = part1->chdir[1];
                closestDist();
            }
            atrenergy += eattractivePscCpsc(1,0);
            }
            if ( (firstT) && (secondT) ) {
            if (secondCH) {
                part2->dir = part2->chdir[1];
                closestDist();
            }
            atrenergy += eattractivePscCpsc(1,1);
            }
            if (secondT) {
            if (firstT && firstCH ) {
                part1->dir = part1->chdir[0];
                closestDist();
            }
            atrenergy += eattractivePscCpsc(0,1);
            }
        }

        if (firstCH)
            part1->dir = olddir1;
        if (secondCH)
            part2->dir = olddir2;

    }
    return repenergy+atrenergy;
}

double PairEnergyCalculator::eattractivePscCpsc(int patchnum1, int patchnum2) {
    int i, intrs;
    double atrenergy, ndist;
    double v1, v2, f0, f1, f2, T1, T2, S1, S2, a;
    double intersections[5];
    Vector vec1, vec2, vec_intrs, vec_mindist;

    //rcut =topo.ia_params[part1->type][part2->type].rcut;

    //interact->halfl =topo.ia_params[part1->type][part2->type].half_len[0];
    //DEBUG_SIM("halfl = %lf", interact->halfl);
    for(i=0;i<5;i++)
        intersections[i]=0;
    bool first;
    if ( (topo.ia_params[part1->type][part2->type].geotype[0] == PSC)||(topo.ia_params[part1->type][part2->type].geotype[0] == CHPSC)||(topo.ia_params[part1->type][part2->type].geotype[0] == TPSC)||(topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) ){
        first = true;
    } else {
        first = false;
    }
    //cospatch =topo.ia_params.pcanglsw;
    //cospatchinr =topo.ia_params.pcangl;
    /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at.
      cut distance C*/
    //DEBUG_SIM("first intersection");
    if (first) {
        intrs = pscIntersect(part1,part2,topo.ia_params[part1->type][part2->type].half_len[0],topo.ia_params[part1->type][part2->type].half_len[1],r_cm, intersections, 0, patchnum1);
    } else {
        intrs = cpscIntersect(part1,part2,topo.ia_params[part1->type][part2->type].half_len[0],topo.ia_params[part1->type][part2->type].half_len[1],r_cm, intersections, 0, patchnum1);
    }
    //DEBUG_SIM("first intersection: done");
    if (intrs <2){
        //DEBUG_SIM("No intersection :(");
        return 0.0; /*sc is all outside patch, attractive energy is 0*/
    }

    T1 = intersections[0]; /*points on sc2*/
    T2 = intersections[1];
    /*2- now do the same oposite way psc1 in patch of psc2*/
    for(i=0;i<5;i++)
        intersections[i]=0;
    //DEBUG_SIM("get vector");
    vec1 = -1.0 * r_cm;
    //DEBUG_SIM("second intersection");
    if (first) {
        intrs = cpscIntersect(part2,part1,topo.ia_params[part1->type][part2->type].half_len[1],topo.ia_params[part1->type][part2->type].half_len[0],vec1, intersections, 1, patchnum2);
    } else {
        intrs = pscIntersect(part2,part1,topo.ia_params[part1->type][part2->type].half_len[1],topo.ia_params[part1->type][part2->type].half_len[0],vec1, intersections, 1, patchnum2);
    }
    if (intrs <2)
        return 0.0; /*sc is all outside patch, attractive energy is 0*/

    S1=intersections[0]; /*points on sc1*/
    S2=intersections[1];

    /*3- scaling function1: dependence on the length of intersetions*/
    // BUG 23.9.2014
    v1=fabs(S1-S2); // v1=0.5*fabs(S1-S2);
    v2=fabs(T1-T2); // v2=0.5*fabs(T1-T2);
    f0=0.5*(v1+v2);

    /*4a- with two intersection pices calculate vector between their CM
      -this is for angular orientation*/
    vec1 = ((S1+S2)*0.5) * part1->dir;
    vec2 = ((T1+T2)*0.5) * part2->dir;
    vec_intrs.x = vec2.x - vec1.x - r_cm.x;
    vec_intrs.y = vec2.y - vec1.y - r_cm.y;
    vec_intrs.z = vec2.z - vec1.z - r_cm.z;
    /*vec_intrs should be from sc1 to sc2*/
    //	    fprintf (stderr, "segments_CM: %.8f %.8f %.8f \n",vec_intrs.x,vec_intrs.y,vec_intrs.z);

    /*4b - calculate closest distance attractive energy from it*/
    vec_mindist = minDistSegments(v1,v2,vec_intrs);
    //	    fprintf (stderr, "segments closest dist: %.8f %.8f %.8f \n",vec_mindist.x,vec_mindist.y,vec_mindist.z);
    ndist=sqrt(DOT(vec_mindist,vec_mindist));
    //dist=DOT(vec_intrs,vec_intrs);
    if (ndist <topo.ia_params[part1->type][part2->type].pdis)
        atrenergy = -topo.ia_params[part1->type][part2->type].epsilon;
    //atrenergy = -1.0;
    else {
        atrenergy = cos(PIH*(ndist-topo.ia_params[part1->type][part2->type].pdis)/topo.ia_params[part1->type][part2->type].pswitch);
        atrenergy *= -atrenergy*topo.ia_params[part1->type][part2->type].epsilon ;
    }

    /*5- scaling function2: angular dependence of patch1*/
    vec1 = 1.0 * vec_intrs;
    //vec1=vecScale(vec_mindist,-1.0);
    vec1=vecPerpProject(&vec1, &part1->dir);
    vec1.normalise();
    a = DOT(vec1, part1->patchdir[patchnum1]);
    f1 = fanglScale(a, 0+2*patchnum1);

    /*6- scaling function3: angular dependence of patch2*/
    vec1 = -1.0 * vec_intrs;
    //vec1=vecScale(vec_mindist,1.0);
    vec1 = vecPerpProject(&vec1, &part2->dir);
    vec1.normalise();
    a = DOT(vec1,part2->patchdir[patchnum2]);
    f2 = fanglScale(a, 1+2*patchnum2);

    /*7- put it all together*/
    atrenergy *=f0*f1*f2;
    //if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
    //	    fprintf (stderr, "attraction  %.8f \n",atrenergy);
    //	    exit(1);

    return atrenergy;
}

double PairEnergyCalculator::eSpaSca() {
    double atrenergy, repenergy, b, f0, halfl;

    //DEBUG    printf ("do energy 111 \n\n");
    closestDist();
    repenergy = eRepulsive();

    if ( ( dist >topo.ia_params[part1->type][part2->type].rcut ) ||
         (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
         topo.ia_params[part1->type][part2->type].exclude ) {
        atrenergy = 0.0;
    } else {
        /*calculate closest distance attractive energy*/
        if (dist <topo.ia_params[part1->type][part2->type].pdis)
            atrenergy = -topo.ia_params[part1->type][part2->type].epsilon;
        else {
            atrenergy = cos(PIH*(dist -topo.ia_params[part1->type][part2->type].pdis)/topo.ia_params[part1->type][part2->type].pswitch);
            atrenergy *= -atrenergy*topo.ia_params[part1->type][part2->type].epsilon ;
        }

        /*scaling function for the length of spherocylinder within cutoff*/

        if (topo.ia_params[part1->type][part2->type].geotype [0] < SP)
            halfl =topo.ia_params[part1->type][part2->type].half_len[0];
        else
        halfl =topo.ia_params[part1->type][part2->type].half_len[1];
        b = sqrt(topo.ia_params[part1->type][part2->type].rcut*topo.ia_params[part1->type][part2->type].rcut - dist*dist);
        if ( contt + b > halfl )
            f0 = halfl;
        else
            f0 = contt + b;
        if ( contt - b < -halfl )
            f0 -= -halfl;
        else
            f0 -= contt - b;
        // BUG FIX 23.9.2014
        atrenergy *= (f0);
        //if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
        //fprintf (stderr, "attraction211  %.8f x: %.8f y: %.8f z: %.8f \n",atrenergy,vec1.x,vec1.y,vec1.z);
        //exit(1);
    }

    return repenergy+atrenergy;
}

double PairEnergyCalculator::ePscSpa() {
    double atrenergy, repenergy;

    //DEBUG_SIM("do energy 211") ;
    closestDist();
    repenergy = eRepulsive();
    //DEBUG_SIM("got the rep. energy");
    if ( ( dist >topo.ia_params[part1->type][part2->type].rcut ) ||
         (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
         topo.ia_params[part1->type][part2->type].exclude ) {
        atrenergy = 0.0;
    } else {
        bool firstCH=false, secondCH=false;
        Vector olddir1 = part1->dir;
        Vector olddir2 = part2->dir;
        if ( (topo.ia_params[part1->type][part2->type].geotype[0] == CHPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) )
            firstCH = true;
        if ( (topo.ia_params[part1->type][part2->type].geotype[1] == CHPSC) || (topo.ia_params[part1->type][part2->type].geotype[1] == TCHPSC) )
            secondCH = true;
        if(firstCH)
            part1->dir = part1->chdir[0];
        if(secondCH)
            part2->dir = part2->chdir[0];

        if ((firstCH) || (secondCH) ) {
            closestDist();
        }
        atrenergy = eattractivePscSpa(0);

        /*addition of interaction of second patches*/
        if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) ||
          (topo.ia_params[part1->type][part2->type].geotype[1] == TPSC) ||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHPSC) ) {
            bool firstT=false, secondT=false;
            if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHPSC) )
            firstT = true;
            if ( (topo.ia_params[part1->type][part2->type].geotype[1] == TPSC) ||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHPSC)  )
            secondT = true;

            if (firstT) {
            if (firstCH) {
                part1->dir = part1->chdir[1];
                closestDist();
            }
            atrenergy += eattractivePscSpa(1);
            }
            if (secondT) {
            if(secondCH) {
                part2->dir = part2->chdir[1];
                closestDist();
            }
            atrenergy += eattractivePscSpa(1);
            }
            if ( (firstT) && (secondT) ) {
              fprintf (stderr, "ERROR PSC should interact s SPA but got two PSC \n");
              //exit(1);
            }
        }

        if (firstCH)
            part1->dir = olddir1;
        if (secondCH)
            part2->dir = olddir2;
    }
    return repenergy+atrenergy;
}

double PairEnergyCalculator::eattractivePscSpa(int patchnum1) {
    double atrenergy, a, b, f0, halfl;
    Vector vec1;
    int which;

    /*calculate closest distance attractive energy*/
    if (dist <topo.ia_params[part1->type][part2->type].pdis)
        atrenergy = -topo.ia_params[part1->type][part2->type].epsilon;
    else {
        atrenergy = cos(PIH*(dist -topo.ia_params[part1->type][part2->type].pdis) /topo.ia_params[part1->type][part2->type].pswitch);
        atrenergy *= -atrenergy*topo.ia_params[part1->type][part2->type].epsilon ;
    }
    /*scaling function: angular dependence of patch1*/
    if (topo.ia_params[part1->type][part2->type].geotype[0] < SP) {
        which = 0;
        vec1=vecPerpProject(&distvec, &part1->dir);
        vec1.normalise();
        a = DOT(vec1, part1->patchdir[patchnum1]);
        halfl=topo.ia_params[part1->type][part2->type].half_len[0];
    } else {
        which = 1;
        vec1 = vecPerpProject(&distvec, &part2->dir);
        vec1.normalise();
        a = DOT(vec1, part2->patchdir[patchnum1]);
        halfl=topo.ia_params[part1->type][part2->type].half_len[1];
    }
    /*scaling function for the length of spherocylinder within cutoff*/

    b = sqrt(topo.ia_params[part1->type][part2->type].rcut*topo.ia_params[part1->type][part2->type].rcut - dist*dist);
    if ( contt + b > halfl )
        f0 = halfl;
    else
        f0 = contt + b;
    if ( contt - b < -halfl )
        f0 -= -halfl;
    else
        f0 -= contt - b;
    // Bug fix 23.9.2014, was f0+1.0 to f0
    atrenergy *= fanglScale(a, which)*(f0);
    //if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
    //fprintf (stderr, "attraction211  %.8f x: %.8f y: %.8f z: %.8f \n",atrenergy,vec1.x,vec1.y,vec1.z);
    //exit(1);

    return atrenergy;
}

double PairEnergyCalculator::eCpscSpa() {
    double atrenergy, repenergy, halfl;

    //DEBUG_SIM("do energy 311") ;
    closestDist();
    repenergy = eRepulsive();
    //DEBUG_SIM("got the rep. energy");

    if (topo.ia_params[part1->type][part2->type].geotype[0] < SP) {
        halfl=topo.ia_params[part1->type][part2->type].half_len[0];
    } else {
        halfl=topo.ia_params[part1->type][part2->type].half_len[1];
    }
    if ( ( dist >topo.ia_params[part1->type][part2->type].rcut ) ||
         (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
         topo.ia_params[part1->type][part2->type].exclude ||
         ( dist >topo.ia_params[part1->type][part2->type].rcut ) ||
         (contt > halfl) ||
         (contt < -halfl) )
        atrenergy = 0.0;
    else {
        bool firstCH=false, secondCH=false;
        Vector olddir1 = part1->dir;
        Vector olddir2 = part2->dir;
        if ( (topo.ia_params[part1->type][part2->type].geotype[0] == CHCPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHCPSC) )
            firstCH = true;
        if ( (topo.ia_params[part1->type][part2->type].geotype[1] == CHCPSC) || (topo.ia_params[part1->type][part2->type].geotype[1] == TCHCPSC) )
            secondCH = true;
        if(firstCH)
            part1->dir = part1->chdir[0];
        if(secondCH)
            part2->dir = part2->chdir[0];

        if ((firstCH) || (secondCH) ) {
            closestDist();
        }
        atrenergy = eattractiveCpscSpa(0);

        /*addition of interaction of second patches*/
        if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TCPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHCPSC) ||
          (topo.ia_params[part1->type][part2->type].geotype[1] == TCPSC) ||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHCPSC) ) {
            bool firstT=false, secondT=false;
            if ( (topo.ia_params[part1->type][part2->type].geotype[0] == TCPSC) || (topo.ia_params[part1->type][part2->type].geotype[0] == TCHCPSC) )
            firstT = true;
            if ( (topo.ia_params[part1->type][part2->type].geotype[1] == TCPSC) ||(topo.ia_params[part1->type][part2->type].geotype[1] == TCHCPSC)  )
            secondT = true;

            if (firstT) {
            if (firstCH) {
                part1->dir = part1->chdir[1];
                closestDist();
            }
            atrenergy += eattractiveCpscCpsc(1,0);
            }
            if (secondT) {
            if(secondCH) {
                part2->dir = part2->chdir[1];
                closestDist();
            }
            atrenergy += eattractiveCpscCpsc(0,1);
            }
            if ( (firstT) && (secondT) ) {
            fprintf (stderr, "ERROR PSC should interact s SPA but got two PSC \n");
            //exit(1);
            }
        }

        if (firstCH)
            part1->dir = olddir1;
        if (secondCH)
            part2->dir = olddir2;

    }
    return repenergy+atrenergy;
}

double PairEnergyCalculator::eattractiveCpscSpa(int patchnum1) {
    double atrenergy, a, b, f0, halfl;
    Vector vec1;
    int which;

    /*if it is in cylindrical part c>-halfl and c<halfl*/
    /*calculate closest distance attractive energy*/
    if (dist <topo.ia_params[part1->type][part2->type].pdis)
        atrenergy = -topo.ia_params[part1->type][part2->type].epsilon;
    else {
        atrenergy = cos(PIH*(dist -topo.ia_params[part1->type][part2->type].pdis)/topo.ia_params[part1->type][part2->type].pswitch);
        atrenergy *= -atrenergy*topo.ia_params[part1->type][part2->type].epsilon ;
    }
    /*scaling function: angular dependence of patch1*/
    if (topo.ia_params[part1->type][part2->type].geotype[0] < SP) {
        which = 0;
        vec1 = vecPerpProject(&distvec, &part1->dir);
        vec1.normalise();
        a = DOT(vec1, part1->patchdir[patchnum1]);
        halfl =topo.ia_params[part1->type][part2->type].half_len[0];
    } else {
        which = 1;
        vec1 = vecPerpProject(&distvec, &part2->dir);
        vec1.normalise();
        a = DOT(vec1,part2->patchdir[patchnum1]);
        halfl =topo.ia_params[part1->type][part2->type].half_len[1];
    }
    /*scaling function for the length of spherocylinder within cutoff*/
    b = sqrt(topo.ia_params[part1->type][part2->type].rcut*topo.ia_params[part1->type][part2->type].rcut - dist*dist);
    if ( contt + b > halfl )
        f0 = halfl;
    else
        f0 = contt + b;
    if ( contt - b < -halfl )
        f0 -= -halfl;
    else
        f0 -= contt - b;
    // BUG 23.9.2014 fix f0+1.0 to f0
    atrenergy *= fanglScale(a, which)*(f0);

    //if (atrenergy < 0) printf ("atraction %f\n",atrenergy);
    //fprintf (stderr, "attraction311  %.8f a: %.8f\n",atrenergy,a);
    //exit(1);

    return atrenergy;
}

double PairEnergyCalculator::e2ScaOr2Spa() {
    double repenergy=0.0, atrenergy=0.0;

    closestDist();

    repenergy = eRepulsive();
    if ( ( dist >topo.ia_params[part1->type][part2->type].rcut ) ||
         (topo.ia_params[part1->type][part2->type].epsilon == 0.0 ) ||
         topo.ia_params[part1->type][part2->type].exclude ) {
        atrenergy = 0.0;
    } else {
        if (dist <topo.ia_params[part1->type][part2->type].pdis)
            atrenergy = -topo.ia_params[part1->type][part2->type].epsilon;
        else  {
            atrenergy = cos(PIH*(dist - topo.ia_params[part1->type][part2->type].pdis)/topo.ia_params[part1->type][part2->type].pswitch);
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

double PairEnergyCalculator::eSpnOrScn() {
    closestDist();
    return eRepulsive();
}

    //add scaling increased if particles are parallel or antiparallel
double PairEnergyCalculator::scparallel(double epsilonparallel,Vector dir1,Vector dir2){
    double cosa;  
  
    cosa=DOT(dir1,dir2);
    if ((epsilonparallel>0 && cosa>0) || (epsilonparallel<0 && cosa<0)) return 1.0 + epsilonparallel*cosa;
    else return 1.0;
}

void PairEnergyCalculator::closestDist() {
    double c, d, halfl;

    //printf("we have %d %d  ",topo.ia_params[part1->type][part2->type].geotype[0],topo.ia_params[part1->type][part2->type].geotype[1] );
    if ((topo.ia_params[part1->type][part2->type].geotype[0] >= SP) && (topo.ia_params[part1->type][part2->type].geotype[1] >= SP)) { /*we have two spheres - most common, do nothing*/
        //printf("we have two spheres ");
        distvec = r_cm;
        dist = sqrt(dotrcm);
        distcm = dist;
    } else {
        if ((topo.ia_params[part1->type][part2->type].geotype[0] < SP) && (topo.ia_params[part1->type][part2->type].geotype[1] < SP)) { /*we have two spherocylinders*/
            distvec = minDistSegments(topo.ia_params[part1->type][part2->type].half_len[0],topo.ia_params[part1->type][part2->type].half_len[1], r_cm);
            dist = sqrt(DOT(distvec,distvec));
        } else {
            if (topo.ia_params[part1->type][part2->type].geotype[0] < SP) { /*We have one spherocylinder -it is first one*/
                halfl=topo.ia_params[part1->type][part2->type].half_len[0];/*finding closest vector from sphyrocylinder to sphere*/
                c = DOT(part1->dir,r_cm);
                if (c >= halfl) d = halfl;
                else {
                    if (c > -halfl) d = c;
                    else d = -halfl;
                }
                contt = c;
                distvec.x = - r_cm.x + part1->dir.x * d;
                distvec.y = - r_cm.y + part1->dir.y * d;
                distvec.z = - r_cm.z + part1->dir.z * d;
                dist=sqrt(DOT(distvec, distvec));
            } else { /*lst option first one is sphere second one spherocylinder*/
                halfl =topo.ia_params[part1->type][part2->type].half_len[1]; /*finding closest vector from sphyrocylinder to sphere*/
                c = DOT(part2->dir,r_cm);
                if (c >= halfl) d = halfl;
                else {
                    if (c > -halfl) d = c;
                    else d = -halfl;
                }
                contt = -c;
                distvec.x = r_cm.x - part2->dir.x * d;
                distvec.y = r_cm.y - part2->dir.y * d;
                distvec.z = r_cm.z - part2->dir.z * d;
                dist=sqrt(DOT(distvec,distvec));
            }
        }
    }
}

double PairEnergyCalculator::fanglScale(double a, int which) {
    double f;

    // TODO for different types
    if (a <=topo.ia_params[part1->type][part2->type].pcanglsw[which])
        f=0.0;
    else {
        if (a >=topo.ia_params[part1->type][part2->type].pcangl[which])
            f=1.0;
        else {
            f = 0.5 - ((topo.ia_params[part1->type][part2->type].pcanglsw[which]
                    +topo.ia_params[part1->type][part2->type].pcangl[which])*0.5 - a )
                    /(topo.ia_params[part1->type][part2->type].pcangl[which] -topo.ia_params[part1->type][part2->type].pcanglsw[which]);
        }
    }

    return f;
}

double PairEnergyCalculator::eRepulsive() {
    double repenergy, en6;

    /* WCA repulsion */
    if (dist >topo.ia_params[part1->type][part2->type].rcutwca) repenergy = 0.0;
    else {
            en6 = pow((topo.ia_params[part1->type][part2->type].sigma / dist),6);
            //repenergy = topo.ia_params[part1->type][part2->type].epsilon*(4*en6*(en6-1) + 1.0);
            repenergy = (4*en6*(en6-1) + 1.0);
    }
    //printf("repenergy: %f dist: %f\n",repenergy, dist);

    return repenergy;
}

// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
//
// FIXED parallel case by calculation both P1 on S0 and P0 on S1 and comparing
//
Vector PairEnergyCalculator::minDistSegments(double halfl1, double halfl2, Vector r_cm) {
    Vector u,v,w,vec;
    double a,b,c,d,e,D,sc,sN,sD,tc,tN,tD;
    bool paralel = false;

    // direction of lines
    u = (2.0*halfl1) * part1->dir; //S1.P1 - S1.P0;
    v = (2.0*halfl2) * part2->dir; //S2.P1 - S2.P0;

    w.x = part2->dir.x*halfl2 - part1->dir.x*halfl1 - r_cm.x;
    w.y = part2->dir.y*halfl2 - part1->dir.y*halfl1 - r_cm.y;
    w.z = part2->dir.z*halfl2 - part1->dir.z*halfl1 - r_cm.z; //S1.P0 - S2.P0;

    a = DOT(u,u);        // always >= 0
    b = DOT(u,v);
    c = DOT(v,v);        // always >= 0
    d = DOT(u,w);
    e = DOT(v,w);
    D = a*c - b*b;       // always >= 0
    sc = D;
    sN = D;
    sD = D;      // sc = sN / sD, default sD = D >= 0
    tc = D;
    tN = D;
    tD = D;      // tc = tN / tD, default tD = D >= 0

    // compute the linetopo.ia_paramseters of the two closest points
    if (D < 0.00000001) { // the lines are almost parallel
        paralel = true;
        sN = 0.0;        // force using point P0 on segment S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    } else {                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    if (fabs(sN) < 0.00000001) sc = 0.0 ;
    else sc = sN / sD;
    if (fabs(tN) < 0.00000001) tc = 0.0 ;
    else tc = tN / tD;

    // get the difference of the two closest points
    //Vector = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
    vec.x = u.x*sc + w.x - v.x*tc;
    vec.y = u.y*sc + w.y - v.y*tc;
    vec.z = u.z*sc + w.z - v.z*tc;

    if(paralel) { // force using point P1 on line S0 and compare the distance to point P0 on segment S1
        Vector vec2;
        // direction of lines
        // note u is v and v is u

        w.x = part1->dir.x*halfl1 - part2->dir.x*halfl2 + r_cm.x;
        w.y = part1->dir.y*halfl1 - part2->dir.y*halfl2 + r_cm.y;
        w.z = part1->dir.z*halfl1 - part2->dir.z*halfl2 + r_cm.z; //S1.P0 - S2.P0;

        // a is c
        // b is same
        // c is a
        d = DOT(v,w); // recalc
        e = DOT(u,w); // recalc
        D = a*c - b*b;
        sc = D;
        sN = D;
        sD = D;      // sc = sN / sD, default sD = D >= 0
        tc = D;
        tN = D;
        tD = D;      // tc = tN / tD, default tD = D >= 0

        // compute the linetopo.ia_paramseters of the two closest points
        if (D < 0.00000001) { // the lines are almost parallel
            sN = 0.0;        // force using point P0 on segment S1
            sD = 1.0;        // to prevent possible division by 0.0 later
            tN = e;
            tD = a;
        }

        if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
            tN = 0.0;
            // recompute sc for this edge
            if (-d < 0.0)
                sN = 0.0;
            else if (-d > c)
                sN = sD;
            else {
                sN = -d;
                sD = c;
            }
        } else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
            tN = tD;
            // recompute sc for this edge
            if ((-d + b) < 0.0)
                sN = 0;
            else if ((-d + b) > c)
                sN = sD;
            else {
                sN = (-d + b);
                sD = c;
            }
        }

        // finally do the division to get sc and tc
        if (fabs(sN) < 0.00000001) sc = 0.0 ;
        else sc = sN / sD;
        if (fabs(tN) < 0.00000001) tc = 0.0 ;
        else tc = tN / tD;

        // get the difference of the two closest points
        //Vector = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
        vec2.x = v.x*sc + w.x - u.x*tc;
        vec2.y = v.y*sc + w.y - u.y*tc;
        vec2.z = v.z*sc + w.z - u.z*tc;

        if(vec2.dot(vec2) < vec.dot(vec))
            return vec2;
    }

    return vec;
}


int PairEnergyCalculator::pscIntersect(Particle *part1, Particle *part2, double halfl1, double halfl2,
                                       Vector r_cm, double intersections[], int which, int patchnum) {
    int intrs;
    double a, b, c, d, e, x1, x2, rcut2;
    Vector cm21, vec1, vec2, vec3, vec4;

    intrs = 0;
    rcut2 =topo.ia_params[part1->type][part2->type].rcut *topo.ia_params[part1->type][part2->type].rcut;
    /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at
      cut distance C*/
    /*1a- test intersection with half planes of patch and look how far they are
      from spherocylinder. If closer then C  we got itersection*/

    /* plane1 */
    /* find intersections of part2 with plane by par1 and patchsides[0] */
    intrs += findIntersectPlane(part1,part2,halfl2,r_cm,part1->patchsides[0+2*patchnum],
           topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],intersections);
    //	    printf("plane1 %d\n", intrs);
    /* plane2 */
    /* find intersections of part2 with plane by par1 and patchsides[1] */
    intrs += findIntersectPlane(part1,part2,halfl2,r_cm,part1->patchsides[1+2*patchnum],
           topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],intersections);

    if ( (intrs == 2 ) && (topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum] <0) ) {
        fprintf (stderr, "ERROR: Patch is larger than 180 degrees and we are getting two segments - this hasnot been programed yet.\n\n");
        //exit (1);
    }
    //	    printf("plane2 %d\n", intrs);

    /*1b- test intersection with cylinder - it is at distance C*/
    if (intrs < 2 )  {
        cm21 = -1.0 * r_cm;
        vec1 = vecCrossProduct(&cm21,&part1->dir);
        vec2 = vecCrossProduct(&part2->dir,&part1->dir);
        a = DOT(vec2,vec2);
        b = 2*DOT(vec1,vec2);
        c = -rcut2 + DOT(vec1,vec1);
        d = b*b - 4*a*c;
        if ( d >= 0) { /*there is intersection with infinite cylinder */
            x1 = (-b+sqrt(d))*0.5/a;/*parameter on line of SC2 determining intersection*/
            if ((x1 >=halfl2) || (x1 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
            else {
                /* vectors from center os sc1 to intersection with infinite cylinder*/
                vec1.x = part2->dir.x*x1-r_cm.x;
                vec1.y = part2->dir.y*x1-r_cm.y;
                vec1.z = part2->dir.z*x1-r_cm.z;
                e = DOT(part1->dir,vec1);
                if ((e >=halfl1) || (e <= -halfl1)) intrs+=0; /*intersection is outside sc1*/
                else {
                    intrs+=testIntrPatch(part1,vec1,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x1,intersections,patchnum);
                }
            }
            if ( d > 0 ){
                x2 = (-b-sqrt(d))*0.5/a;/*parameter on line of SC2 determining intersection*/
                if ((x2 >=halfl2) || (x2 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
                else {
                    vec2.x = part2->dir.x*x2 - r_cm.x;
                    vec2.y = part2->dir.y*x2 - r_cm.y;
                    vec2.z = part2->dir.z*x2 - r_cm.z;
                    e = DOT(part1->dir,vec2);
                    if ((e >=halfl1) || (e <= -halfl1)) intrs+=0; /*intersection is outside sc1*/
                    else {
                        intrs+=testIntrPatch(part1,vec2,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x2,intersections,patchnum);
                    }
                }
            }
        }
    }
    //	    printf ("cylinder %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
    /*1c- test intersection with spheres at the end - it is at distace C*/
    if (intrs < 2 )  {
        /*centers of spheres*/
        /*relative to the CM of sc2*/
        vec1.x =  part1->dir.x * halfl1 - r_cm.x;
        vec1.y =  part1->dir.y * halfl1 - r_cm.y;
        vec1.z =  part1->dir.z * halfl1 - r_cm.z;
        vec2.x = -part1->dir.x * halfl1 - r_cm.x;
        vec2.y = -part1->dir.y * halfl1 - r_cm.y;
        vec2.z = -part1->dir.z * halfl1 - r_cm.z;

        /*sphere1*/
        a = DOT(part2->dir,part2->dir);
        b = 2.0*DOT(vec1,part2->dir);
        c = DOT(vec1,vec1)- rcut2;
        d = b*b-4*a*c;
        if (d >= 0) { /*if d<0 there are no intersections*/
            x1= (-b + sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >=halfl2) || (x1 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
            else {
                vec3.x = part2->dir.x*x1 - r_cm.x;
                vec3.y = part2->dir.y*x1 - r_cm.y;
                vec3.z = part2->dir.z*x1 - r_cm.z;
                e = DOT(part1->dir,vec3);
                if ((e >= halfl1) || (e <= -halfl1)) { /*if not intersection is inside sc1*/
                    intrs += testIntrPatch(part1,vec3,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x1,intersections,patchnum);
                }
            }
            if ( d > 0) {
                x2= (-b - sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >=halfl2) || (x2 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
                else {
                    vec4.x = part2->dir.x*x2 - r_cm.x;
                    vec4.y = part2->dir.y*x2 - r_cm.y;
                    vec4.z = part2->dir.z*x2 - r_cm.z;
                    e = DOT(part1->dir,vec4);
                    if ((e >=halfl1) || (e <= -halfl1)) { /*if not intersection is inside sc1*/
                        intrs += testIntrPatch(part1,vec4,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x2,intersections,patchnum);
                    }
                }
            }
        }
        //		printf ("sphere1 %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
        /*sphere2*/
        a = DOT(part2->dir,part2->dir);
        b = 2.0*DOT(vec2,part2->dir);
        c = DOT(vec2,vec2)- rcut2;
        d = b*b-4*a*c;
        if (d >= 0) { /*if d<0 there are no intersections*/
            x1= (-b + sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >=halfl2) || (x1 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
            else {
                vec3.x = part2->dir.x*x1 - r_cm.x;
                vec3.y = part2->dir.y*x1 - r_cm.y;
                vec3.z = part2->dir.z*x1 - r_cm.z;
                e = DOT(part1->dir,vec3);
                if ((e >=halfl1) || (e <= -halfl1)) { /*if not intersection is inside sc1*/
                    intrs+=testIntrPatch(part1,vec3,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x1,intersections,patchnum);
                }
            }
            if ( d > 0 ) {
                x2= (-b - sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >=halfl2) || (x2 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
                else {
                    vec4.x = part2->dir.x*x2 - r_cm.x;
                    vec4.y = part2->dir.y*x2 - r_cm.y;
                    vec4.z = part2->dir.z*x2 - r_cm.z;
                    e = DOT(part1->dir,vec4);
                    if ((e >=halfl1) || (e <= -halfl1)) { /*if not intersection is inside sc1*/
                        intrs+=testIntrPatch(part1,vec4,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x2,intersections,patchnum);
                    }
                }
            }
        }
        //		printf ("sphere2 %d\n", intrs);
    }

    /*1d- if there is only one itersection shperocylinder ends within patch wedge
      set as second intersection end inside patch*/
    if (intrs < 2 )  {
        /*whole spherocylinder is in or all out if intrs ==0*/
        vec1.x = part2->dir.x*halfl2 - r_cm.x;
        vec1.y = part2->dir.y*halfl2 - r_cm.y;
        vec1.z = part2->dir.z*halfl2 - r_cm.z;
        /*vector from CM of sc1 to end of sc2*/
        /*check is is inside sc1*/
        a=DOT(vec1,part1->dir);
        vec3.x = vec1.x - part1->dir.x*a;
        vec3.y = vec1.y - part1->dir.y*a;
        vec3.z = vec1.z - part1->dir.z*a;
        b=DOT(vec3,vec3);
        d = fabs(a)-halfl1;
        if ( d <= 0)
            c = b; /*is inside cylindrical part*/
        else
            c = d*d + b; /*is inside caps*/
        /*c is distance squared from line or end to test if is inside sc*/
        if (c < rcut2)
            intrs += testIntrPatch(part1,vec1,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],halfl2,intersections,patchnum);
        if (intrs < 2 ) {
            vec2.x = -part2->dir.x*halfl2 - r_cm.x;
            vec2.y = -part2->dir.y*halfl2 - r_cm.y;
            vec2.z = -part2->dir.z*halfl2 - r_cm.z;
            /*check is is inside sc1*/
            a = DOT(vec2,part1->dir);
            vec4.x = vec2.x - part1->dir.x*a;
            vec4.y = vec2.y - part1->dir.y*a;
            vec4.z = vec2.z - part1->dir.z*a;
            b = DOT(vec4,vec4);
            d = fabs(a) -halfl1;
            if (d <= 0)
                c = b; /*is inside cylindrical part*/
            else
                c = d*d + b; /*is inside caps*/
            /*c is distance squared from line or end to test if is inside sc*/
            if (c < rcut2)
                intrs += testIntrPatch(part1,vec2,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],-1.0*halfl2,intersections,patchnum);
        }
        //		    printf ("ends %d\n", intrs);
    }
    return intrs;
}

int PairEnergyCalculator::testIntrPatch(Particle *part1, Vector vec, double cospatch,
                                        double ti, double intersections[], int patchnum) {
    double a;
    int i, intrs;

    intrs=0;
    /*test if we have intersection*/
    /* do projection to patch plane*/
    vec = vecPerpProject(&vec,&part1->dir);
    vec.normalise();
    /* test angle distance from patch*/
    a = DOT(part1->patchdir[patchnum],vec);
    if (a >= cospatch) {
        intrs=1;
        i=0;
        while (intersections[i] !=0) {
            if (ti == intersections[i])
                intrs=0; /* found intersection we already have -it is at boundary*/
            i++;
        }
        if (intrs > 0)
            intersections[i]=ti;
    }

    return intrs;
}

int PairEnergyCalculator::findIntersectPlane(Particle *part1, Particle *part2, double halfl2,
                                             Vector r_cm, Vector w_vec, double cospatch, double intersections[]) {
    int i, intrs;
    double a, c, d, ti, disti;
    Vector d_vec, nplane = vecCrossProduct(&part1->dir,&w_vec);

    nplane.normalise();
    w_vec.normalise();
    a =  DOT(nplane, part2->dir);
    if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
    else {
        ti = DOT(nplane,r_cm)/a;
        if ((ti  > halfl2 ) || (ti < -halfl2)) intrs=0; /* there is no intersection plane sc is too short*/
        else {
            d_vec.x = ti * part2->dir.x - r_cm.x; /*vector from intersection point to CM*/
            d_vec.y = ti * part2->dir.y - r_cm.y;
            d_vec.z = ti * part2->dir.z - r_cm.z;
            c = DOT (d_vec, w_vec);
            if ( c * cospatch < 0) intrs=0; /* the intersection in plane is on other side of patch */
            else {
                d = fabs(DOT (d_vec, part1->dir)) - halfl2;
                if (d <= 0) disti = c*c; /*is inside cylinder*/
                else disti = d*d + c*c; /*is inside patch*/
                if (disti >topo.ia_params[part1->type][part2->type].rcut *
                       topo.ia_params[part1->type][part2->type].rcut) intrs=0; /* the intersection is outside sc */
                else {
                    intrs=1;
                    i=0;
                    while (intersections[i] !=0) {
                        if (ti == intersections[i]) intrs=0; /* found intersection we already have -it is at boundary*/
                        i++;
                    }
                    if (intrs > 0) {
                        intersections[i]=ti;
                    }
                }
            }
        }

    }
    return intrs;
}



Vector PairEnergyCalculator::vecPerpProject(Vector *A, Vector *B) {
    Vector pp;
    double dp;

    dp=DOT((*A),(*B));

    pp.x = A->x - B->x*dp;
    pp.y = A->y - B->y*dp;
    pp.z = A->z - B->z*dp;
    //    fprintf (stderr, "pp x: %.8f y: %.8f z: %.8f \n",pp.x,pp.y,pp.z);
    return pp;
}


int PairEnergyCalculator::cpscIntersect(Particle *part1, Particle *part2, double halfl1, double halfl2,
                                        Vector r_cm, double intersections[], int which, int patchnum) {
    int intrs;
    double a, b, c, d, e, x1, x2, rcut2;
    Vector cm21, vec1, vec2, vec3, vec4;

    intrs = 0;
    rcut2 =topo.ia_params[part1->type][part2->type].rcut *topo.ia_params[part1->type][part2->type].rcut;
    /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at
      cut distance C*/
    /*1a- test intersection with half planes of patch and look how far they are
      from spherocylinder. If closer then C  we got itersection*/

    /* plane1 */
    /* find intersections of part2 with plane by par1 and part1->patchsides[0] */
    intrs += findIntersectPlanec(part1,part2,halfl2,r_cm,part1->patchsides[0+2*patchnum],
           topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],intersections);
    //	    printf("plane1 %d\n", intrs);
    /* plane2 */
    /* find intersections of part2 with plane by par1 and part1->patchsides[1] */
    intrs += findIntersectPlanec(part1,part2,halfl2,r_cm,part1->patchsides[1+2*patchnum],
           topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],intersections);

    if ( (intrs == 2 ) && (topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum] < 0) ) {
        fprintf (stderr, "ERROR: Patch is larger than 180 degrees and we are getting two segments - this hasnot been programed yet.\n\n");
        //exit (1);
    }
    //	    printf("plane2 %d\n", intrs);

    /*1b- test intersection with cylinder - it is at distance C*/
    if (intrs < 2 )  {
        cm21 = -1.0 * r_cm;
        vec1 = vecCrossProduct(&cm21,&part1->dir);
        vec2 = vecCrossProduct(&part2->dir,&part1->dir);
        a = DOT(vec2,vec2);
        b = 2*DOT(vec1,vec2);
        c = -rcut2 + DOT(vec1,vec1);
        d = b*b - 4*a*c;
        if ( d >= 0) { /*there is intersection with infinite cylinder */
            x1 = (-b+sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >=halfl2) || (x1 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
            else {
                /* vectors from center os sc1 to intersection with infinite cylinder*/
                vec1.x=part2->dir.x*x1-r_cm.x;
                vec1.y=part2->dir.y*x1-r_cm.y;
                vec1.z=part2->dir.z*x1-r_cm.z;
                e = DOT(part1->dir,vec1);
                if ((e >=halfl1) || (e <= -halfl1)) intrs+=0; /*intersection is outside sc1*/
                else {
                    intrs += testIntrPatch(part1,vec1,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x1,intersections,patchnum);
                }
            }
            if ( d > 0 ){
                x2 = (-b-sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >=halfl2) || (x2 <= -halfl2)) intrs+=0; /*intersection is outside sc2*/
                else {
                    vec2.x = part2->dir.x*x2-r_cm.x;
                    vec2.y = part2->dir.y*x2-r_cm.y;
                    vec2.z = part2->dir.z*x2-r_cm.z;
                    e = DOT(part1->dir,vec2);
                    if ((e >=halfl1) || (e <= -halfl1)) intrs+=0; /*intersection is outside sc1*/
                    else {
                        intrs += testIntrPatch(part1,vec2,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x2,intersections,patchnum);
                    }
                }
            }
        }
    }
    //	    printf ("cylinder %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
    /*1c- test intersection with plates at the end - it is at distace C and in wedge*/
    if (intrs < 2 )  {
        a =  DOT(part1->dir, part2->dir);
        if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
        else {
            /*plane cap1*/
            vec1.x= r_cm.x + halfl1*part1->dir.x;
            vec1.y= r_cm.y + halfl1*part1->dir.y;
            vec1.z= r_cm.z + halfl1*part1->dir.z;
            x1 = DOT(part1->dir,vec1)/a; /*parameter on line of SC2 determining intersection*/
            if ((x1 > halfl2 ) || (x1 < -halfl2)) intrs+=0; /* there is no intersection plane sc is too short*/
            else {
                vec2.x = x1*part2->dir.x - vec1.x; /*vector from ENDPOINT to intersection point */
                vec2.y = x1*part2->dir.y - vec1.y;
                vec2.z = x1*part2->dir.z - vec1.z;
                b = DOT (vec2, vec2);
                if (b > rcut2) intrs+=0;   /* the intersection is outside sc */
                else {
                    intrs += testIntrPatch(part1,vec2,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x1,intersections,patchnum);
                }
            }
            //		    printf ("plane cap1 %d %f\n", intrs, x1);
            /*plane cap2*/
            vec1.x= r_cm.x - halfl1*part1->dir.x;
            vec1.y= r_cm.y - halfl1*part1->dir.y;
            vec1.z= r_cm.z - halfl1*part1->dir.z;
            x2 = DOT(part1->dir,vec1)/a; /*parameter on line of SC2 determining intersection*/
            if ((x2  > halfl2 ) || (x2 < -halfl2)) intrs+=0; /* there is no intersection plane sc is too short*/
            else {
                vec2.x = x2*part2->dir.x - vec1.x; /*vector from ENDPOINT to intersection point */
                vec2.y = x2*part2->dir.y - vec1.y;
                vec2.z = x2*part2->dir.z - vec1.z;
                b = DOT (vec2, vec2);
                if (b > rcut2) intrs+=0;   /* the intersection is outside sc */
                else {
                    intrs += testIntrPatch(part1,vec2,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],x2,intersections,patchnum);
                }
            }
            //		    printf ("plane cap2 %d %f\n", intrs,x2);

        }
    }

    /*1d- if there is only one itersection shperocylinder ends within patch wedge
      set as second intersection end inside patch*/
    if (intrs < 2 )  {
        /*whole spherocylinder is in or all out if intrs ==0*/
        vec1.x = part2->dir.x*halfl2 - r_cm.x;
        vec1.y = part2->dir.y*halfl2 - r_cm.y;
        vec1.z = part2->dir.z*halfl2 - r_cm.z;
        /*vector from CM of sc1 to end of sc2*/
        /*check is is inside sc1*/
        a=DOT(vec1,part1->dir);
        vec3.x = vec1.x - part1->dir.x*a;
        vec3.y = vec1.y - part1->dir.y*a;
        vec3.z = vec1.z - part1->dir.z*a;
        b=DOT(vec3,vec3);
        d = fabs(a)-halfl1;
        if ( d <= 0) { /*is in cylindrical part*/
            /*c is distance squared from line or end to test if is inside sc*/
            if (b < rcut2) intrs += testIntrPatch(part1,vec1,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],halfl2,intersections,patchnum);
        }
        if (intrs < 2 ) {
            vec2.x = -part2->dir.x*halfl2 - r_cm.x;
            vec2.y = -part2->dir.y*halfl2 - r_cm.y;
            vec2.z = -part2->dir.z*halfl2 - r_cm.z;
            /*check is is inside sc1*/
            a=DOT(vec2,part1->dir);
            vec4.x = vec2.x - part1->dir.x*a;
            vec4.y = vec2.y - part1->dir.y*a;
            vec4.z = vec2.z - part1->dir.z*a;
            b=DOT(vec4,vec4);
            d = fabs(a) -halfl1;
            if (d <= 0) {
                /*c is distance squared from line or end to test if is inside sc*/
                if (b < rcut2) intrs += testIntrPatch(part1,vec2,topo.ia_params[part1->type][part2->type].pcanglsw[which+2*patchnum],-1.0*halfl2,intersections,patchnum);
            }
        }
                //    printf ("ends %d\n", intrs);
    }

    return intrs;
}

int PairEnergyCalculator::findIntersectPlanec(Particle *part1, Particle *part2, double halfl,
                                              Vector r_cm, Vector w_vec, double cospatch, double intersections[]) {
    int i, intrs=0;
    double a, c, d, ti, disti;
    Vector nplane, d_vec;

    nplane=vecCrossProduct(&part1->dir,&w_vec);
    nplane.normalise();
    w_vec.normalise();
    a =  DOT(nplane, part2->dir);
    if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
    else {
        ti = DOT(nplane,r_cm)/a;
        if ((ti  > halfl ) || (ti < -halfl)) intrs=0; /* there is no intersection plane sc is too short*/
        else {
            d_vec.x = ti*part2->dir.x - r_cm.x; /*vector from intersection point to CM*/
            d_vec.y = ti*part2->dir.y - r_cm.y;
            d_vec.z = ti*part2->dir.z - r_cm.z;
            c = DOT (d_vec, w_vec);
            if ( c *cospatch < 0) intrs=0; /* the intersection in plane is on other side of patch */
            else {
                d = fabs(DOT (d_vec, part1->dir)) - halfl;
                if (d <= 0) {
                    disti= c*c; /*is inside cylinder*/
                    if (disti >topo.ia_params[part1->type][part2->type].rcut*topo.ia_params[part1->type][part2->type].rcut)
                        intrs=0; /* the intersection is outside sc */
                    else {
                        intrs=1;
                        i=0;
                        while (intersections[i] !=0) {
                            if (ti == intersections[i]) intrs=0; /* found intersection we already have -it is at boundary*/
                            i++;
                        }
                        if (intrs > 0) intersections[i]=ti;
                    }
                }
            }
        }
    }
    return intrs;
}

