#include "totalenergycalculator.h"


double TotalEnergyCalculator::operator ()(Particle* target, int mode, int chainnum) {
    long i=0,j=0;

    //DEBUG_SIM("Calculate the energy with mode %d", mode)
    double energy = 0;

    //
    // Calculates energy between particle "target" and the rest.  Returns energy
    //
    if(mode == 1) {
        pairE.setPrimaryParticle(target);
        if (sim->pairlist_update) {
#ifdef OMP
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
            for (i = 0; i < target->neighborCount; i++){
               energy += (pairE)(target->neighborID[i], &conf->particleStore[target->neighborID[i]]);
            }
        } else {

#ifdef OMP
#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
#endif
            for (i = 0; i < (long)conf->particleStore.size(); i++) {
                if(target != &conf->particleStore[i]) {
                    energy += (pairE)(i, &conf->particleStore[i]);
                }
            }
        }
        //add interaction with external potential
        if (topo->exter.exist)
            energy += extere2(target);
    } // END of MODE=1



    //
    // Calculates energy between particle "target" and the rest. skipping
    // particles from the given chain -particles has to be sorted in chain!!
    // so similar to energy one but with chain exception
    //
    else if(mode == 2){

        pairE.setPrimaryParticle(target);

        //#ifdef OMP
        //#pragma omp parallel for private(i) reduction (+:energy) schedule (dynamic)
        //#endif

        for (i = 0; i < (long)conf->particleStore.size(); i++) {
            if (i != conf->chainlist[chainnum][j]) {
                if(target != &conf->particleStore[i])
                energy += (pairE)(i, &conf->particleStore[i]);
            }
            else {
                j++;
            }
        }
        j++;

        //add interaction with external potential
        if (topo->exter.exist)
            energy+= extere2(target);

    // END of MODE=2
    } else if(mode == 0) { // Calculates energy between all pairs. Returns energy

#ifdef OMP
#pragma omp parallel for private(i,j) reduction (+:energy) schedule (dynamic)
#endif

        for (i = 0; i < (long)conf->particleStore.size() - 1; i++) {

            pairE.setPrimaryParticle(&conf->particleStore[i]);
            for (j = i + 1; j < (long)conf->particleStore.size(); j++) {
                energy += (pairE)(j, &conf->particleStore[j]);
            }

            //for every particle add interaction with external potential
            if (topo->exter.exist)
                energy += extere2(&conf->particleStore[i]);
        }

        //add interaction of last particle with external potential
        if (topo->exter.exist)
            energy+= extere2(&conf->particleStore.back());

    } else {
        fprintf(stderr, "ERROR: Wrong mode (%d) was given to calc_energy!", mode);
        return 0.0;
    }
    //  DEBUG_SIM("Will return energy from calc_energy")
    //printf("energymove %f\n",energy);
    return energy;
}

double TotalEnergyCalculator::extere2(Particle* target) {

    double repenergy=0.0,atrenergy=0.0;  /* energy*/

    double ndist;                        /* distance for CM of interacting line segment*/
    double halfl;

    Vector olddir;

    /* calcualte distance to center of mass*/
    if (  target->pos.z < 0  ) {
        rcmz = conf->box.z * (target->pos.z - (double)( (long)(target->pos.z - 0.5) ) );
    } else {
        rcmz = conf->box.z * (target->pos.z - (double)( (long)(target->pos.z + 0.5) ) );
    }

    project.x=0;
    project.y=0;

    if (rcmz < 0) {
        dist = - rcmz;
        positive = false;
        interendz = -1.0;
        project.z = 1.0;
    } else {
        dist = rcmz;
        positive = true;
        interendz = 1.0;
        project.z = -1.0;
    }

    dotrcm = rcmz * rcmz;
    if ( dotrcm > topo->exter.sqmaxcut) return 0.0;  /* distance so far that even spherocylinders cannot be within cutoff  */
    distvec.z = r_cm.z;
    distcm = dist;
    box = conf->box;
    part1 = target;
    param = &topo->exter.interactions[target->type];
    halfl = 0.5* topo->exter.interactions[target->type].len[0];
    ndist = dist;
    orientin = true;
    orient = 0.0;

    exter2ClosestDist();

    /* now we have closest distance so we can calculate repulsion*/

    /*
     *  Direct code insert to avoid struct interact
     */
    //repenergy = erepulsive(&interact);
    /*******************************************************************************/
    double en6;
    /* WCA repulsion */
    if (dist > topo->ia_params[part1->type][part2->type].rcutwca) {
        repenergy = 0.0;
    } else {
            en6 = pow((topo->ia_params[part1->type][part2->type].sigma / dist),6);
            repenergy = 4*en6*(en6-1) + 1.0;
    }
    //printf("repenergy: %f dist: %f\n", repenergy, dist);
    /***********************************************************************************/

    //printf("dist: %f", dist);

    /*save chiral stuff*/
    olddir = part1->dir;

    if ((param->geotype[0] == CHCPSC) || (param->geotype[0] == CHPSC)) {
        part1->dir = part1->chdir[0];
        exter2ClosestDist();
    }
    if (( dist > param->rcut ) || (param->epsilon == 0.0 ) ||
        ((part1->patchdir[0].z >0) && (positive))
            || ((part1->patchdir[0].z <0)&&(!positive)) ) {
        atrenergy = 0.0;
    } else {
        atrenergy = exter2Atre(&ndist,0,halfl);
    }
    if ((param->geotype[0] == TCPSC)||(param->geotype[0] == TPSC)||
      (param->geotype[0] == TCHCPSC)||(param->geotype[0] == TCHPSC)) {
        if ((param->geotype[0] == TCHCPSC)||(param->geotype[0] == TCHPSC)) {
            part1->dir = part1->chdir[1];
            exter2ClosestDist();
        }
        exter2ClosestDist();
        if (( dist > param->rcut ) || (param->epsilon == 0.0 ) ||
              ( (part1->patchdir[1].z >0) && (positive) )
                || ( (part1->patchdir[1].z <0) && (!(positive)) )  )
            atrenergy += 0.0;
        else {
            atrenergy += exter2Atre(&ndist,1,halfl);
        }

    }

    if ((param->geotype[0] == CHCPSC)||(param->geotype[0] == CHPSC)||
      (param->geotype[0] == TCHCPSC)||(param->geotype[0] == TCHPSC) ) {
        part1->dir = olddir;
    }

    //printf("%f   %f \n",conf->particleStore[target].pos.z*conf->box.z,repenergy+atrenergy);
    return repenergy+atrenergy;
}

void TotalEnergyCalculator::exter2ClosestDist() {

    if (rcmz < 0) {
        dist = -(rcmz);
        positive = false;
        interendz = -1.0;
        project.z = 1.0;
    } else {
        dist = (rcmz);
        positive = true;
        interendz = 1.0;
        project.z = -1.0;
    }
    /*psc closest is allways end closer to wall*/
    if (param->geotype[0] < SP ){
        /*calculate closest point distance*/
        if (part1->dir.z > 0) {
            if (positive) {
                orientin = false;
                orient = -1.0;
                dist = rcmz - part1->dir.z * param->half_len[0];
            } else {
                orientin = true;
                orient = 1.0;
                dist = -( rcmz + part1->dir.z * param->half_len[0]);
            }
        } else {
            if (positive) {
                orientin = true;
                orient = 1.0;
                dist = rcmz + part1->dir.z * param->half_len[0];
            } else {
                orientin = false;
                orient = -1.0;
                dist = -( rcmz -part1->dir.z * param->half_len[0]);
            }
        }
    }
}

double TotalEnergyCalculator::exter2Atre(double *ndist, int numpatch, double halfl) {
    Vector pbeg,pend;             /* projected spherocylinder begining and end*/
    double a,length1,length2, f0,f1;
    Vector cm1,cm2;               /* centrar of interacting segments */
    int line;
    Vector partbeg,partend;       /*closest and furthest point of particle*/
    Vector inters;
    double atrenergy=0.0;

    /*interaction with PATCHY SPHEROCYLINDERS*/
    if ((param->geotype[0] < SP)&&(param->geotype[0] > SCA)) {
        //printf("partdir: %f %f %f \n ",part1->dir.x,part1->dir.y,part1->dir.z);
        //printf("patchdir: %f %f %f \n ",part1->patchdir[0].x,part1->patchdir[0].y,part1->patchdir[0].z);
        /* calculate position of closest and furthest point (begining and end of spherocylinder)*/
        a = (orientin-0.5)*2;  /*if orientin a =1 else a=-1 */
        partbeg.x =  a * part1->dir.x * halfl;
        partbeg.y =  a * part1->dir.y * halfl;
        partbeg.z = rcmz + a * part1->dir.z *halfl;
        partend.x =  - a * part1->dir.x * halfl;
        partend.y =  - a * part1->dir.y * halfl;
        partend.z = rcmz - a * part1->dir.z * halfl;


        //printf("partbeg %f %f %f partend %f %f %f  \n",partbeg.x,partbeg.y,partbeg.z,partend.x,partend.y,partend.z);
        /*calculate interacting line segment and its cm of spherocylinder*/
        /*calculate end point z*/
        if ( (param->rcut - dist)/fabs(part1->dir.z) < 2.0*halfl ){
            /*if cutoff goes through spherocylinder the end point is at cutoff*/
            interendz *= param->rcut;
        } else {
            /*endpoint is at the end of spherocylinders*/
            interendz = partend.z;
        }
        /*calculate CM of interacting line segment of spherocylinder*/
        if (positive) {
            cm1.z = AVER(interendz,dist);
        } else {
            cm1.z = AVER(interendz,-dist);
        }
        if (part1->dir.z != 0.0 ) {
            a = (interendz - cm1.z ) / part1->dir.z;
            length1= -orient*2.0*a;
            a = a + orient*halfl;
        } else {
            a = 0.0;
            length1 = 2.0*halfl;
        }
        //printf("len1: %f rcm %f interz %f cutoff %f  \n",length1,rcmz, interendz,dist);
        cm1.x = part1->dir.x * a;
        cm1.y = part1->dir.y * a;

        /* we have interacting segment*/
        if ((param->geotype[0] == CPSC)||(param->geotype[0] == CHCPSC)) {
            /*CPSC type*/
            if ( ((interendz >= dist)&&(positive))
                 || ((interendz <= -dist)&&(!(positive))) ) {
                /*test if projection is not all out of interaction*/
                line = cpscWall(&pbeg,&pend,&project,&part1->dir, &param->half_len[0],&param->rcut,&partbeg,&partend);
                //printf("line: %d beg %f %f end %f %f  \n",line,pbeg.x,pbeg.y,pend.x,pend.y);
            } else {
                line = 0;
            }
        } else {
            /*PSC and CHPSC interaction with wall */
            line = pscWall(&pbeg,&pend,&project,&part1->dir,&param->rcut,&partbeg,&partend);
            //printf("line: %d beg %f %f end %f %f  \n",line,pbeg.x,pbeg.y,pend.x,pend.y);
        }

        if (line > 0) {
            /*cm2 by average begining and end*/
            cm2.x = AVER(pbeg.x,pend.x);
            cm2.y = AVER(pbeg.y,pend.y);
            cm2.z = 0.0;

            /*length by size of end-benining*/
            length2 = sqrt( (pend.x-pbeg.x)*(pend.x-pbeg.x)+(pend.y-pbeg.y)*(pend.y-pbeg.y) );

            inters.x = cm2.x - cm1.x;
            inters.y = cm2.y - cm1.y;
            inters.z = cm2.z - cm1.z;
            //printf("cm2 %f %f %f inters %f %f %f \n",cm2.x,cm2.y,cm2.z,inters.x,inters.y,inters.z);
            *ndist = sqrt(DOT(inters,inters));
            if (*ndist < param->pdis) {
                atrenergy = -param->epsilon;
            }
            else {
                atrenergy= cos(PIH*(*ndist-param->pdis)/param->pswitch);
                atrenergy *= -atrenergy*param->epsilon;
            }
            /*  scaling function1: dependence on the length of intersetions plus*/
            f0=(length1 + length2)*0.5;
            /*scaling with angle*/
            f1 = fabs(part1->patchdir[numpatch].z);
            atrenergy *= f0*f1;

            //printf(" %f   %f    %f %f %f %f %f \n",conf->particleStore[target].pos.z*conf->box.z,atrenergy, area, length1, length2,f0,ndist);
            //printf("%f %f %f %f\n",pbeg.x,pbeg.y,pend.x,pend.y);
        } else {
            atrenergy = 0.0;
        }

    } else {
        if (*ndist < param->pdis)
            atrenergy = -param->epsilon;
        else  {
            atrenergy= cos(PIH*(*ndist-param->pdis)/param->pswitch);
            atrenergy *= -atrenergy*param->epsilon;
        }
        /*add wall scaling wall area/ particle arear.. to reflect that we have a wall not sphere */
        atrenergy *= (param->rcut*param->rcut - (*ndist)*(*ndist))/(param->sigma*param->sigma) ;
    }

    return atrenergy;
}

double TotalEnergyCalculator::exterAtre(double *ndist, int numpatch, double halfl) {
    double area,a,b,c,r2;
    double atrenergy=0.0;  // energy
    bool countend;
    Vector cm1,cm2;               // centrar of interacting segments
    Vector pbeg,pend;             // projected spherocylinder begining and end
    Vector inters,newdir;
    Vector pbeg1,pend1,pbeg2,pend2,pextr1,pextr2,pextr3,pextr4; //additinal point of projected patch for calculation of area
    double length1, cuttoproject, f0;
    int line, line1, line2,extra;
    Vector partbeg,partend;        //closest and furthest point of particle

    //interaction with PATCHY SPHEROCYLINDER
    if ((param->geotype[0] < SP)&&(param->geotype[0] > SCA)) {
        //printf("partdir: %f %f %f \n ",part1->dir.x,part1->dir.y,part1->dir.z);
        //printf("patchdir: %f %f %f \n ",part1->patchdir[numpatch].x,part1->patchdir[numpatch].y,part1->patchdir[numpatch].z);
        // calculate position of closest and furthest point (begining and end of spherocylinder)
        a = (orientin-0.5)*2;  //if orientin a =1 else a=-1
        partbeg.x =  a * part1->dir.x * halfl;
        partbeg.y =  a * part1->dir.y * halfl;
        partbeg.z = rcmz + a * part1->dir.z * halfl;
        partend.x =  - a * part1->dir.x * halfl;
        partend.y =  - a * part1->dir.y * halfl;
        partend.z = rcmz - a * part1->dir.z * halfl;

        //printf("partbeg %f %f %f partend %f %f %f  \n",partbeg.x,partbeg.y,partbeg.z,partend.x,partend.y,partend.z);
        //calculate interacting line segment and its cm of spherocylinder
        //calculate end point z
        if ( (param->rcut - dist)/fabs(part1->dir.z) < halfl*2.0 ){
            //if cutoff goes through spherocylinder the end point is at cutoff
            interendz *= param->rcut;
        } else {
            //endpoint is at the end of spherocylinders
            interendz = partend.z;
        }
        //calculate CM of interacting line segment of spherocylinder
        if (positive) {
            cm1.z = AVER(interendz,dist);
        } else {
            cm1.z = AVER(interendz,-dist);
        }
        if (part1->dir.z != 0.0 ) {
            a = (interendz - cm1.z ) / part1->dir.z;
            length1= -orient*2.0*a;
            a = a + orient*halfl;
        } else {
            a = 0.0;
            length1 = 2.0*halfl;
        }
        //printf("len1: %f rcm %f interz %f cutoff %f  \n",length1,rcmz, interendz,dist);
        cm1.x = part1->dir.x * a;
        cm1.y = part1->dir.y * a;
        //calculate projection on wall as infinite line and make it interacting segment
        if (part1->patchdir[numpatch].z != 0) {
            cuttoproject = -param->rcut*part1->patchdir[numpatch].z; //z coordinate of point where projection is in cut distance
            if ( ((partend.z < cuttoproject)&&(positive)) || ((cuttoproject < partend.z)&&(!positive)) ){
                cuttoproject = partend.z;
            }
        } else {
            cuttoproject = partbeg.z;
        }
        //printf("cutproject %f \n",cuttoproject);
        //printf("cm1 %f %f %f \n",cm1.x, cm1.y,cm1.z );

        // we have interacting segment
        if ((param->geotype[0] == CPSC)||(param->geotype[0] == CHCPSC)) {
            //CPSC type
            if ( ((cuttoproject >= dist)&&(positive)) || ((cuttoproject <= -dist)&&(!(positive))) ){
                //test if projection is not all out of interaction
                line = cpscWall(&pbeg,&pend,&part1->patchdir[numpatch],&part1->dir, \
                        &param->half_len[0],&param->rcut,&partbeg,&partend);
                //printf("line: %d beg %f %f end %f %f  \n",line,pbeg.x,pbeg.y,pend.x,pend.y);
            } else {
                line = 0;
            }
        } else {
            //PSC and CHPSC interaction with wall
            line = pscWall(&pbeg,&pend,&part1->patchdir[numpatch],&part1->dir,&param->rcut,&partbeg,&partend);
            //printf("line: %d beg %f %f end %f %f  \n",line,pbeg.x,pbeg.y,pend.x,pend.y);
        }

        if (line > 0) {
            area = 0.0;
            //project cutoff boudaries
            if (line == 2 ) {
                //if projection end is on sphere of begining don't care about cylinder cutoff
                extra = 0;
            } else {
                extra = cutProjectAtWall(&pextr1, &pextr2, &pextr3, &pextr4, &part1->patchdir[numpatch],  \
                     &part1->dir, &param->rcut, &partbeg, &partend,&pend,&cuttoproject);
            }
            //printf("extr1: %d %f %f extr2 %f %f extr3 %f %f extr4 %f %f \n",extra,pextr1.x,pextr1.y,pextr2.x,pextr2.y,pextr3.x,pextr3.y,pextr4.x,pextr4.y);

            //project patch boundaries on the first side
            newdir=part1->patchsides[0+2*numpatch];
            line1 = cpscWall(&pbeg1,&pend1,&newdir,&part1->dir, \
                &param->half_len[0],&param->rcut,&partbeg,&partend);
            if ( ((param->geotype[0] == PSC)||(param->geotype[0] == CHPSC)) ) {
                line1 = pscWall(&pbeg1,&pend1,&newdir,&part1->dir,&param->rcut,&partbeg,&partend);
            }
            //printf("line1: %d beg1 %f %f end1 %f %f  \n",line1,pbeg1.x,pbeg1.y,pend1.x,pend1.y);

            //project patch boundaries on the second side
            newdir=part1->patchsides[1+2*numpatch];
            line2 = cpscWall(&pbeg2,&pend2,&newdir,&part1->dir,&param->half_len[0],&param->rcut,&partbeg,&partend);
            if ( ((param->geotype[0] == PSC)||(param->geotype[0] == CHPSC)) ) {
                line2 = pscWall(&pbeg2,&pend2,&newdir,&part1->dir,&param->rcut,&partbeg,&partend);
            }
            //printf("line2: %d beg2 %f %f end2 %f %f \n",line2,pbeg2.x,pbeg2.y,pend2.x,pend2.y);

            //calculate area
            if (extra == 0) {
                //thish should only happen when there is PSC interacting only with end
                if (line1 == 0) {
                    if (line2==0) {
                        //circle around middle-pbeg
                        area = PI*( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y));
                    }
                    else{
                        // circle around middle-pbeg minus circle segment
                        a = AVER(pbeg2.x,pend2.x);
                        b = AVER(pbeg2.y,pend2.y);
                        c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
                        r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); /*radius squared*/
                        area =  r2*(PI-acos(sqrt(c/r2))) + sqrt(r2*c-c*c);
                    }
                } else {
                    if (line2==0) {
                        // circle around middle-pbeg minus circle segment
                        a = AVER(pbeg1.x,pend1.x);
                        b = AVER(pbeg1.y,pend1.y);
                        c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); //height of triangle to segment
                        r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); //radius squared
                        area =  r2*(PI-acos(sqrt(c/r2))) + sqrt(r2*c-c*c);
                    } else {
                        //area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pend,&pend1,&pbeg1,NULL,NULL); // B B2 E2 E E1 B1
                        //circle minus two circle segments
                        r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); //radius squared
                        area = r2*PI;
                        a = AVER(pbeg1.x,pend1.x);
                        b = AVER(pbeg1.y,pend1.y);
                        c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); //height of triangle to segment
                        area +=  -r2*acos(sqrt(c/r2)) + sqrt(r2*c-c*c);
                        a = AVER(pbeg2.x,pend2.x);
                        b = AVER(pbeg2.y,pend2.y);
                        c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); //height of triangle to segment
                        area +=  -r2*acos(sqrt(c/r2)) + sqrt(r2*c-c*c);
                    }
                }

            } else {
                b = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend.y)- (pextr2.x-pend.x)*(pextr4.y-pextr2.y));//pend on 42
                c = fabs((pextr1.x-pextr3.x)*(pextr3.y-pend.y)- (pextr3.x-pend.x)*(pextr1.y-pextr3.y));//pend on 13
                if ( ( b< ZEROTOL) || ( c< ZEROTOL) )
                    countend = false;
                else
                    countend = true;
                if (line1 == 0) {
                    if (line2 == 0) {
                        if ( countend ) {
                            area = areaEightPoints(&pbeg,&pextr2,&pextr4,&pend,&pextr3,&pextr1,NULL,NULL);// B 2 4 E 3 1
                        } else
                            area = areaEightPoints(&pbeg,&pextr2,&pextr4,&pextr3,&pextr1,NULL,NULL,NULL);// B 2 4 3 1
                    } else {
                        a = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend2.y)- (pextr2.x-pend2.x)*(pextr4.y-pextr2.y));
                        if ( a< ZEROTOL) /*pend2 on 42*/ {
                            if ( countend ) {
                                area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pextr1,NULL); // B B2 E2 4 E 3 1
                            } else {
                                area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pextr1,NULL,NULL); // B B2 E2 4 3 1
                            }
                        } else {
                            a = fabs((pextr1.x-pextr3.x)*(pextr3.y-pend2.y)- (pextr3.x-pend2.x)*(pextr1.y-pextr3.y));
                            if ( a< ZEROTOL) /*pend2 on 13*/ {
                                area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr1,NULL,NULL,NULL,NULL); // B B2 E2 1
                            } else { /*pend2 on 34 or on begining sphere  of psc*/
                                if (line2 == 2) {
                                    if ( countend ) {
                                        area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pextr1,NULL); // B B2 E2 4 E 3 1
                                    } else {
                                        area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pextr1,NULL,NULL); // B B2 E2 4 3 1
                                    }
                                } else {
                                    if ( countend ) {
                                        area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pend,&pextr3,&pextr1,NULL,NULL); // B B2 E2 E 3 1
                                    } else {
                                        area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr3,&pextr1,NULL,NULL,NULL); // B B2 E2 3 1
                                    }
                                }
                            }
                        }
                    }
                } else {
                    a = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend1.y)- (pextr2.x-pend1.x)*(pextr4.y-pextr2.y));
                    if ( a< ZEROTOL) /*pend1 on 42*/ {
                        if (line2 == 0) {
                            area = areaEightPoints(&pbeg,&pextr2,&pend1,&pbeg1,NULL,NULL,NULL,NULL); /* B 2 E1 B1 */
                        } else {
                            area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pend1,&pbeg1,NULL,NULL,NULL); /* B B2 E2 E1 B1 */
                        }
                    } else {
                        a = fabs((pextr1.x-pextr3.x)*(pextr3.y-pend1.y)- (pextr3.x-pend1.x)*(pextr1.y-pextr3.y));
                        if ( a< ZEROTOL) /*pend1 on 13*/ {
                            if (line2 == 0) {
                                if (countend) {
                                    area =  areaEightPoints(&pbeg,&pextr2,&pextr4,&pend,&pextr3,&pend1,&pbeg1,NULL);  /* B 2 4 E 3 E1 B1 */
                                } else {
                                    area =  areaEightPoints(&pbeg,&pextr2,&pextr4,&pextr3,&pend1,&pbeg1,NULL,NULL);  /* B 2 4 3 E1 B1 */
                                }
                            } else {
                                a = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend2.y)- (pextr2.x-pend2.x)*(pextr4.y-pextr2.y));
                                if ( a< ZEROTOL) /*pend2 on 42*/ {
                                    if (countend)
                                        area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pend1,&pbeg1); /* B B2 E2 4 E 3 E1 B1 */
                                    else
                                        area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 4 3 E1 B1 */
                                } else {
                                    a = fabs((pextr3.x-pextr1.x)*(pextr1.y-pend2.y)- (pextr1.x-pend2.x)*(pextr3.y-pextr1.y));
                                    if ( a< ZEROTOL) /*pend2 on 31*/ {
                                        area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pend1,&pbeg1,NULL,NULL,NULL); /* B B2 E2 E1 B1 */
                                    } else { /*pend2 close to 34 or on begining sphere  of psc*/
                                        if (line2 == 2) {
                                            if (countend)
                                                area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pend1,&pbeg1); /* B B2 E2 4 E 3 E1 B1 */
                                            else
                                                area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 4 3 E1 B1 */
                                        } else {
                                            if (countend)
                                                area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pend,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 E 3 E1 B1 */
                                            else
                                                area = areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr3,&pend1,&pbeg1,NULL,NULL); /* B B2 E2 3 E1 B1 */
                                        }
                                    }
                                }
                            }
                        } else {/*pend1 close to 34 or on beging sphere for psc*/
                            if (line2 == 0) {
                                if (line1 ==2) {
                                    if (countend)
                                        area = areaEightPoints(&pbeg,&pextr2,&pextr4,&pend,&pextr3,&pend1,&pbeg1,NULL); /* B 2 4 E 3 E1 B1*/
                                    else {
                                        area = areaEightPoints(&pbeg,&pextr2,&pextr4,&pextr3,&pend1,&pbeg1,NULL,NULL); /* B 2 4 3 E1 B1*/
                                    }
                                } else {
                                    if (countend)
                                        area = areaEightPoints(&pbeg,&pextr2,&pextr4,&pend,&pend1,&pbeg1,NULL,NULL); /* B 2 4 E E1 B1*/
                                    else {
                                        area = areaEightPoints(&pbeg,&pextr2,&pextr4,&pend1,&pbeg1,NULL,NULL,NULL); /* B 2 4 E1 B1*/
                                    }
                                }

                            } else {
                                a = fabs((pextr4.x-pextr2.x)*(pextr2.y-pend2.y)- (pextr2.x-pend2.x)*(pextr4.y-pextr2.y));
                                if ( a< ZEROTOL) /* pend2 on 42 */ {
                                    if (countend)
                                        area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pend1,&pbeg1,NULL);  /* B B2 E2 4 E E1 B1 */
                                    else
                                        area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend1,&pbeg1,NULL,NULL);  /* B B2 E2 4 E1 B1 */
                                } else { /*pend1 and pend2 close to 34 or on beging sphere for psc*/
                                    if (line2 == 2) {
                                        if (line1 == 2) {
                                            if (countend)
                                                area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pextr3,&pend1,&pbeg1); /* B B2 E2 4 E 3 E1 B1 */
                                            else
                                                area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 4 3 E1 B1 */
                                        } else {
                                            if (countend)
                                                area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend,&pend1,&pbeg1,NULL); /* B B2 E2 4 E E1 B1 */
                                            else
                                                area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr4,&pend1,&pbeg1,NULL,NULL); /* B B2 E2 4 E1 B1 */
                                        }


                                    } else {
                                        if (line1 == 2) {
                                            if (countend)
                                                area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pend,&pextr3,&pend1,&pbeg1,NULL); /* B B2 E2 E 3 E1 B1 */
                                            else
                                                area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pextr3,&pend1,&pbeg1,NULL,NULL); /* B B2 E2 3 E1 B1 */
                                        } else {
                                            if (countend)
                                                area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pend,&pend1,&pbeg1,NULL,NULL); /* B B2 E2 E E1 B1 */
                                            else
                                                area =  areaEightPoints(&pbeg,&pbeg2,&pend2,&pend1,&pbeg1,NULL,NULL,NULL); /* B B2 E2 E1 B1 */
                                        }
                                    }
                                }
                            }
                        }
                    }
                } /*extra != 0*/
                if ((param->geotype[0] == PSC)||(param->geotype[0] == CHPSC)) {
                    if (line1==2) {
                        /* add circle segment*/
                        a = AVER(pextr1.x,pend1.x); /*end to cutoff - pextr1 ,pend1 */
                        b = AVER(pextr1.y,pend1.y);
                        c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
                        r2 = ( (partbeg.x-pend1.x)*(partbeg.x-pend1.x) + (partbeg.y-pend1.y)*(partbeg.y-pend1.y)); /*radius squared*/
                        area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
                        a = AVER(pbeg.x,pbeg1.x); /* between beginings - pbeg ,pbeg1 */
                        b = AVER(pbeg.y,pbeg1.y);
                        c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
                        area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
                    } else {
                        if (line1==0) {
                            /* add circle segment*/
                            a = AVER(pextr1.x,pbeg.x); /* begining to cutoff*/
                            b = AVER(pextr1.y,pbeg.y);
                            c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
                            r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); /*radius squared*/
                            area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
                        }
                    }
                        if (line2==2) {
                        /* add circle segment*/
                        a = AVER(pextr3.x,pend2.x); /*end to cutoff - pextr3 ,pend2 */
                        b = AVER(pextr3.y,pend2.y);
                        c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
                        r2 = ( (partbeg.x-pend2.x)*(partbeg.x-pend2.x) + (partbeg.y-pend2.y)*(partbeg.y-pend2.y)); /*radius squared*/
                        area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
                        a = AVER(pbeg.x,pbeg2.x); /* between beginings - pbeg ,pbeg2 */
                        b = AVER(pbeg.y,pbeg2.y);
                        c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); /*height of triangle to segment*/
                        area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
                    } else {
                        if (line2==0) {
                            // add circle segment
                            a = AVER(pextr3.x,pbeg.x); // begining to cutoff
                            b = AVER(pextr3.y,pbeg.y);
                            c = (partbeg.x-a)*(partbeg.x-a) + (partbeg.y-b)*(partbeg.y-b); //height of triangle to segment
                            r2 = ( (partbeg.x-pbeg.x)*(partbeg.x-pbeg.x) + (partbeg.y-pbeg.y)*(partbeg.y-pbeg.y)); //radius squared
                            area +=  r2*acos(sqrt(c/r2)) - sqrt(r2*c-c*c);
                        }
                    }
                }
            }
            //area finished

            //cm2 by average begining and end
            cm2.x = AVER(pbeg.x,pend.x);
            cm2.y = AVER(pbeg.y,pend.y);
            cm2.z = 0.0;
            //length by size of end-benining
            //length2 = sqrt( (pend.x-pbeg.x)*(pend.x-pbeg.x)+(pend.y-pbeg.y)*(pend.y-pbeg.y) );
            inters.x = cm2.x - cm1.x;
            inters.y = cm2.y - cm1.y;
            inters.z = cm2.z - cm1.z;
            //printf("cm2 %f %f %f inters %f %f %f \n",cm2.x,cm2.y,cm2.z,inters.x,inters.y,inters.z);
            *ndist = sqrt(DOT(inters,inters));
            if (*ndist < param->pdis) {
                atrenergy = -param->epsilon;
            }
            else {
                atrenergy= cos(PIH*(*ndist-param->pdis)/param->pswitch);
                atrenergy *= -atrenergy*param->epsilon;
            }
            //  scaling function1: dependence on the length of intersetions plus SCALING WITH AREA
            f0=(length1 + area / param->sigma)*0.5;
            atrenergy *= f0;
            //printf(" %f   %f    %f %f %f %f %f %d %d %d \n",conf->particleStore[target].pos.z*conf->box.z,atrenergy, area, length1, length2,f0,ndist,extra,line1,line2);
            //printf("%f %f %f %f\n",pbeg.x,pbeg.y,pend.x,pend.y);
            //printf("%f %f %f %f %f %f\n",pbeg2.x,pend2.y,pextr2.x,pextr2.y,pextr1.x,pextr1.y);
        } else {
            atrenergy = 0.0;
        }

    } else {
        if (*ndist < param->pdis)
            atrenergy = -param->epsilon;
        else  {
            atrenergy= cos(PIH*(*ndist-param->pdis)/param->pswitch);
            atrenergy *= -atrenergy*param->epsilon;
        }
        //add wall scaling wall area/ particle arear.. to reflect that we have a wall not sphere
        atrenergy *= (param->rcut*param->rcut - (*ndist)*(*ndist))/(param->sigma*param->sigma) ;
    }


    //printf("%f   %f \n",conf->particleStore[target].pos.z*conf->box.z,atrenergy);
    return atrenergy;
}

double TotalEnergyCalculator::areaEightPoints(Vector *p1, Vector *p2, Vector *p3, Vector *p4,
                                              Vector *p5, Vector *p6, Vector *p7, Vector *p8) {
    double area =0.0;
    Vector vec1,vec2;

    /*area by half vector2 cross product
     |(pbegining-pbegining)x(pend-pbegining)|/2  */
    vec1.x = p2->x - p1->x;
    vec1.y = p2->y - p1->y;
    vec2.x = p3->x - p2->x;
    vec2.y = p3->y - p2->y;
    area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
    //printf("3");

    if (p4 != NULL) {
        vec1.x = p3->x - p1->x;
        vec1.y = p3->y - p1->y;
        vec2.x = p4->x - p3->x;
        vec2.y = p4->y - p3->y;
        area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
        //printf("4");

        if (p5 != NULL) {
            vec1.x = p4->x - p1->x;
            vec1.y = p4->y - p1->y;
            vec2.x = p5->x - p4->x;
            vec2.y = p5->y - p4->y;
            area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
            //printf("5");

            if (p6 != NULL) {
                vec1.x = p5->x - p1->x;
                vec1.y = p5->y - p1->y;
                vec2.x = p6->x - p5->x;
                vec2.y = p6->y - p5->y;
                area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
                //printf("6");

                if (p7 != NULL) {
                    vec1.x = p6->x - p1->x;
                    vec1.y = p6->y - p1->y;
                    vec2.x = p7->x - p6->x;
                    vec2.y = p7->y - p6->y;
                    area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
                    //printf("7");

                    if (p8 != NULL) {
                        vec1.x = p7->x - p1->x;
                        vec1.y = p7->y - p1->y;
                        vec2.x = p8->x - p7->x;
                        vec2.y = p8->y - p7->y;
                        area += fabs(vec1.x*vec2.y -vec1.y*vec2.x)*0.5;
                        //printf("8");
                    }
                }
            }
        }
    }

    return area;
}

int TotalEnergyCalculator::pscWall(Vector *pbeg, Vector *pend, Vector *projectdir, Vector *partdir,
                                   double *cutdist, Vector *partbeg, Vector *partend) {
    Vector vec1;
    double k,x1,x2,y1,y2,a,b,c,e,d;

    if (( (positive)&& (projectdir->z > 0) ) || ( (!(positive))&& (projectdir->z < 0) ))
        return 0;
    if ( fabs(partbeg->z) > (*cutdist) )
        return 0;
    /* we might have interacting segment*/
    x2 = 0.0;
    y2 = 0.0;
    /*begining point*/
    /*if begining projected along particle direction is within cutoff */
    if (fabs(partdir->z) > ZEROTOL2) {
        projectinZ(partbeg,partdir,pbeg);
        a=0;
    }
    else {
        /*we need some starting point*/
        vec1.x = 2.0*partbeg->x - partend->x;
        vec1.y = 2.0*partbeg->y - partend->y;
        vec1.z = 2.0*partbeg->z - partend->z;
        projectinZ(&vec1,projectdir,pbeg);
        a=1;
    }
    if (partdir->z != 0) {
        b = fabs(partbeg->z / partdir->z);
    } else {
        b = (*cutdist)+1.0;
    }


    if ( (b > (*cutdist)) || (a==1)) {
        /*else beginig is at sphere, find intersections with sphere of cutoff radius*/
        if ( fabs(projectdir->z) > ZEROTOL2) {
            projectinZ(partbeg,projectdir,pend);
        } else {
            pend->x = pbeg->x + projectdir->x;
            pend->y = pbeg->y + projectdir->y;
        }
        if (pend->y == pbeg->y) {
            y1=pbeg->y;
            y2=pbeg->y;
            a=sqrt( (*cutdist)*(*cutdist) - partbeg->z*partbeg->z - (pbeg->y-partbeg->y)*(pbeg->y-partbeg->y) );
            x1 = partbeg->x + a;
            x2 = partbeg->x - a;
            if (pend->x > pbeg->x) {/*select the right intersection*/
                pbeg->x = x2;
                x2 = x1;
            } else {
                pbeg->x = x1;
            }
            pbeg->y = y1;
        } else {
            k = (pend->x - pbeg->x)/ (pend->y - pbeg->y);
            a = k*k +1;
            b = partbeg->y + k*k*pbeg->y - k*pbeg->x + k*partbeg->x;
            c = partbeg->y*partbeg->y + partbeg->z*partbeg->z - (*cutdist)*(*cutdist) + (k*pbeg->y - pbeg->x + partbeg->x)*(k*pbeg->y - pbeg->x + partbeg->x);
            e = b*b-a*c;
            if (e < 0) {
                return 0; /*tehre might be no intersection with sphere*/
            }
            d = sqrt(e);
            if (pend->y > pbeg->y) {/*select the right intersection*/
                y1 = (b - d ) /a;
                y2 = (b + d ) /a;
            }
            else {
                y1 = (b + d ) /a;
                y2 = (b - d ) /a;
            }
            x1 = k * (y1 - pbeg->y) + pbeg->x;
            x2 = k * (y2 - pbeg->y) + pbeg->x;
            pbeg->x = x1;
            pbeg->y = y1;
            pbeg->z = 0.0;
        }
    }
    //printf("pscwall beg %f %f \n",pbeg->x,pbeg->y);

    /*end point*/
    a = -(*cutdist) * projectdir->z; /*z coordinate of point where projection is in cut distance*/
    //printf("sphere end %f %f   ",a,partend->z);
    if ( ((partend->z < a)&&(positive)) || ((a < partend->z)&&(!(positive))) ){
        /*end is within cut off - second sphere*/
        /*if this is the case vec1 is end of pherocylinder and pend is its projection*/
        if (projectdir->z != 0) {
            projectinZ(partend,projectdir,pend);
        } else {
            pend->x = pbeg->x + projectdir->x;
            pend->y = pbeg->y + projectdir->y;
        }
        if (pend->y == pbeg->y) {
            y1=pend->y;
            y2=pend->y;
            a=sqrt( (*cutdist)*(*cutdist) - partend->z*partend->z - (pend->y-partend->y)*(pend->y-partend->y) );
            x1 = partend->x + a;
            x2 = partend->x - a;
            if (pbeg->x > pend->x) {/*select the right intersection*/
                pend->x = x2;
            } else {
                pend->x = x1;
            }
            pend->y = y1;
        } else {
            k = (pbeg->x - pend->x)/ (pbeg->y - pend->y);
            a = k*k +1;
            b = partend->y + k*k*pend->y - k*pend->x + k*partend->x;
            c = partend->y*partend->y + partend->z*partend->z - (*cutdist)*(*cutdist) + (k*pend->y - pend->x + partend->x)*(k*pend->y - pend->x + partend->x);
            e = b*b-a*c;
            if (e < 0) {
                return 0; /*there might be no intersection with sphere*/
            }
            d = sqrt(e);
            if (pbeg->y > pend->y) {/*select the right intersection*/
                y1 = (b - d ) /a;
                y2 = (b + d ) /a;
            }
            else {
                y1 = (b + d ) /a;
                y2 = (b - d ) /a;
            }
            x1 = k * (y1 - pend->y) + pend->x;
            x2 = k * (y2 - pend->y) + pend->x;
            pend->x = x1;
            pend->y = y1;
            pend->z = 0.0;
        }
    } else {
        if ( ((partbeg->z < a)&&(positive)) || ((a < partbeg->z)&&(!(positive))) ) {
            /*end is at cutoff going through cylindrical part*/
            //printf("cylinder ");
            b = (a - partbeg->z)/ partdir->z;
            vec1.x = partbeg->x + b * partdir->x;
            vec1.y = partbeg->y + b * partdir->y;
            vec1.z = a;
            projectinZ(&vec1,projectdir,pend);
        } else {
            /* also projected end is within the same sphere as begining- no contribution from cylinder*/
            if (x2 == 0.0 ) {
                //printf("sphere beg ");
                if (projectdir->z != 0) {
                    projectinZ(partbeg,projectdir,pend);
                } else {
                    pend->x = pbeg->x + projectdir->x;
                    pend->y = pbeg->y + projectdir->y;
                }

                if (pend->y == pbeg->y) {
                    y1=pbeg->y;
                    y2=pbeg->y;
                    a=sqrt( (*cutdist)*(*cutdist) - partbeg->z*partbeg->z - (pbeg->y-partbeg->y)*(pbeg->y-partbeg->y) );
                    x1 = partbeg->x + a;
                    x2 = partbeg->x - a;
                    if (pend->x > pbeg->x) {/*select the right intersection*/
                        pend->x = x1;
                    } else {
                        pend->x = x2;
                    }
                    pend->y = y1;
                } else {
                    k = (pend->x - pbeg->x)/ (pend->y - pbeg->y);
                    a = k*k +1;
                    b = partbeg->y + k*k*pbeg->y - k*pbeg->x + k*partbeg->x;
                    c = partbeg->y*partbeg->y + partbeg->z*partbeg->z - (*cutdist)*(*cutdist) + (k*pbeg->y - pbeg->x + partbeg->x)*(k*pbeg->y - pbeg->x + partbeg->x);
                    e = b*b-a*c;
                    if (e < 0) {
                        return 0; /*tehre might be no intersection with sphere*/
                    }
                    d = sqrt(e);
                    if (pend->y > pbeg->y) {/*select the right intersection*/
                        y1 = (b - d ) /a;
                        y2 = (b + d ) /a;
                    }
                    else {
                        y1 = (b + d ) /a;
                        y2 = (b - d ) /a;
                    }
                    x1 = k * (y1 - pbeg->y) + pbeg->x;
                    x2 = k * (y2 - pbeg->y) + pbeg->x;
                    pend->x = x1;
                    pend->y = y1;
                    pend->z = 0.0;
                }
            } else {
                pend->x = x2;
                pend->y = y2;
                pend->z = 0.0;
            }
            return 2; /*line end is on sphere of particle begining = no cylindrical cutoff*/
        }
    }


    return 1;
}


int TotalEnergyCalculator::cpscWall(Vector *pbeg, Vector *pend, Vector *projectdir, Vector *partdir,
                                    double *halfl, double *cutdist, Vector *partbeg, Vector *partend) {
    Vector vec1;
    double a;

    if (( (positive)&& (projectdir->z >= 0) ) || ( (!(positive))&& (projectdir->z <= 0) ))
        return 0;
    /*if projected closer point beoynd cutoff no interaction*/
    /*project begining of spherocylinder*/
    vec1.x = partbeg->x;
    vec1.y = partbeg->y;
    vec1.z = partbeg->z;

    if (-vec1.z/projectdir->z < (*cutdist) ) {
        projectinZ(&vec1,projectdir,pbeg);
    } else {
        return 0;
    }

    /* we have interacting segment*/
    if (-partend->z/projectdir->z < (*cutdist) ) {
        /*whole segment interacts*/
        vec1.z = partend->z;
    } else {
        vec1.z = -(*cutdist)*projectdir->z;
    }
    if (partdir->z != 0.0)
        a = (vec1.z - (rcmz)) / partdir->z;
    else {
        if (orientin)
        a = -(*halfl);
        else
        a = (*halfl);
    }
    vec1.x = partdir->x * a;
    vec1.y = partdir->y * a;
    projectinZ(&vec1,projectdir,pend);

    return 1;
}


int TotalEnergyCalculator::cutProjectAtWall(Vector *pextr1, Vector *pextr2, Vector *pextr3, Vector *pextr4,
                                            Vector *projectdir, Vector *partdir, double *cutdist, Vector *partbeg,
                                            Vector *partend, Vector *pend, double *cuttoproject) {
    double y1,y2,O2z,det,a,b,dirydirz,dir2x,dir2y,dir2z,dirzldiry;

    dirydirz = partdir->y * partdir->z;
    dir2x = partdir->x * partdir->x;
    dir2y = partdir->y * partdir->y;
    dir2z = partdir->z * partdir->z;
    a = 1/(dir2x+dir2y);

    if (partdir->x != 0) {
        O2z = partbeg->z * partbeg->z;
        b=dir2y*dir2z*O2z - (dir2x+dir2y) * (O2z*(dir2x+dir2z)- (*cutdist)*(*cutdist)*dir2x);
        if (b < 0 ) {
            /*no cutoff from cylindrical part*/
            return 0;
        }
        det = sqrt(b);
        y1 = partbeg->y + (dirydirz*partbeg->z + det )*a;
        y2 = partbeg->y + (dirydirz*partbeg->z - det )*a;

        if (( (partdir->x > 0)&&(!(orientin)) ) || ( (partdir->x < 0)&&(orientin) ))  {
            pextr1->y = y1;
            pextr2->y = y2;
        } else {
            pextr1->y = y2;
            pextr2->y = y1;
        }
        pextr1->x = partbeg->x  + (partbeg->z*partdir->z -  (pextr1->y - partbeg->y)*partdir->y) / partdir->x;
        pextr2->x = partbeg->x  + (partbeg->z*partdir->z -  (pextr2->y - partbeg->y)*partdir->y) / partdir->x;

        O2z = partend->z * partend->z;
        b= dir2y*dir2z*O2z - (dir2x+dir2y) * (O2z*(dir2x+dir2z)- (*cutdist)*(*cutdist)*dir2x);
        if (b >= 0) { /*we have intersections from end*/
            det = sqrt(b);
            y1 = partend->y + (dirydirz * partend->z + det )*a;
            y2 = partend->y + (dirydirz * partend->z - det )*a;
            //printf("det %f y1 %f y2 %f \n", det,y1,y2);
            if (( (partdir->x > 0)&&(!(orientin)) ) || ( (partdir->x < 0)&&(orientin) ))  {
                pextr3->y = y1;
                pextr4->y = y2;
            } else {
                pextr3->y = y2;
                pextr4->y = y1;
            }
            pextr3->x = partend->x  + (partend->z*partdir->z -  (pextr3->y - partend->y)*partdir->y) / partdir->x;
            pextr4->x = partend->x  + (partend->z*partdir->z -  (pextr4->y - partend->y)*partdir->y) / partdir->x;
        } else {
            /*no intersection at the end the cutoff intersects the plane
             in the perpendicular projection of line segemnt, so we have to use that point */
            if (partdir->z == 0) {
                fprintf (stderr, "\nERROR: Something went wrong in calculation of projection.\n\n");
                exit (1);
            } else {
                a = ((*cuttoproject) - partbeg->z)/ partdir->z;
                //if ( projectdir->y * partdir->x < 0  )
                pextr3->x = partbeg->x + a * partdir->x;
                pextr3->y = partbeg->y + a * partdir->y;
                pextr3->z = (*cuttoproject);
                //printf("before proj %f %f dir %f %f %f ",pextr3->x,pextr3->y,projectdir->x,projectdir->y,projectdir->z);
                projectinZ(pextr3,projectdir,pextr4);
                pextr3->x = pextr4->x;
                pextr3->y = pextr4->y;
                pextr3->z = 0.0;
                //printf("after proj %f %f \n",pextr3->x,pextr3->y);
                return 2;
            }
        }
    } else {
        if (partdir->y != 0) {
            dirzldiry = partdir->z/partdir->y;
            y1 = partbeg->y + partbeg->z * dirzldiry;
            det = sqrt( (*cutdist)*(*cutdist) - partbeg->z * partbeg->z * (1+dirzldiry*dirzldiry) );
            if (( (partdir->y > 0)&&(!(orientin)) ) || ( (partdir->y < 0)&&(orientin) ))  {
                pextr1->x = partbeg->x  + det;
                pextr2->x = partbeg->x  - det;
            } else {
                pextr1->x = partbeg->x  - det;
                pextr2->x = partbeg->x  + det;
            }
            pextr1->y = y1;
            pextr2->y = y1;

            y1 = partend->y + partend->z * dirzldiry;
            b = (*cutdist)*(*cutdist) - partend->z * partend->z * (1+dirzldiry*dirzldiry);
            if (b >= 0) { /*we have intersections from end*/
                det = sqrt(b);
                if (( (partdir->y > 0)&&(!(orientin)) ) || ( (partdir->y < 0)&&(orientin) )) {
                    pextr3->x = partend->x  + det;
                    pextr4->x = partend->x  - det;
                } else {
                    pextr3->x = partend->x  - det;
                    pextr4->x = partend->x  + det;
                }
                pextr3->y = y1;
                pextr4->y = y1;
            } else {
                /*no intersection at the end the cutoff intersects the plane
                 in the perpendicular projection of line segemnt, so we have to use that point */
                if (partdir->z == 0) {
                    fprintf (stderr, "\nERROR: Something went wrong in calculation of projection.\n\n");
                    exit (1);
                } else {
                    a = ((*cutdist) - partbeg->z)/ partdir->z;
                    y1 = a * partdir->y + partbeg->y;
                    if ( projectdir->x * partdir->y > 0  ) {
                        pextr3->x = a * partdir->x  + partbeg->x;
                        pextr3->y = y1;
                        pextr4->x = pend->x;
                        pextr4->y = pend->y;
                    }else {
                        pextr3->x = pend->x;
                        pextr3->y = pend->y;
                        pextr4->x = a * partdir->x  + partbeg->x;
                        pextr4->y = y1;
                    }
                }
            }


        } else {
            return 0; /* if perpendicular to plane we don't have any intersections*/
        }


    }

    return 1;
}
