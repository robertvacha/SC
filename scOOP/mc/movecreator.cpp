/** @file movecreator.cpp*/

#include "movecreator.h"

#ifdef ENABLE_MPI
# include <mpi.h>
extern MPI_Datatype MPI_vector, MPI_Particle, MPI_exchange;
#endif

double MoveCreator::particleMove() {
    double edriftchanges =0.0;
    long target;

    /*=== This is a particle move step ===*/
    target = ran2() * (long)conf->particleStore.size();

    if ( !( ((sim->wlm[0] == 3) || (sim->wlm[1] == 3) ) && (target == 0) ) && \
    ((ran2() < 0.5) || (topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type].geotype[0] >= SP)) ) { /* no rotation for spheres */
        //target = 1;
        //printf ("displacement\n\n");
        edriftchanges = partDisplace(target);

    } else {
        /*=== Rotation step ===*/
        edriftchanges = partRotate(target);
    }
    /*=== End particle move step ===*/

    return edriftchanges;
}

double MoveCreator::partDisplace(long target) {
    double edriftchanges,energy,enermove,wlener;
    Vector orig, dr, origsyscm;
    int reject=0,wli;
    double radiusholemax_orig=0;

    /*=== Displacement step ===*/
    edriftchanges =0.0;
    origsyscm.x = 0;
    origsyscm.y = 0;
    origsyscm.z = 0;
    energy = (*calcEnergy)(target, 1, 0);

    orig = conf->particleStore[target].pos;
    dr.random();
    //ran = sqrt(ran2(&seed));
    dr.x *= sim->trans[conf->particleStore[target].type].mx/conf->box.x;
    dr.y *= sim->trans[conf->particleStore[target].type].mx/conf->box.y;
    dr.z *= sim->trans[conf->particleStore[target].type].mx/conf->box.z;
    if ( ((sim->wlm[0] == 3)||(sim->wlm[1] == 3)) && (target == 0) ) {
        dr.z = 0;
        dr.y = 0;
        dr.x = 0;
    }
    conf->particleStore[target].pos.x += dr.x;
    conf->particleStore[target].pos.y += dr.y;
    conf->particleStore[target].pos.z += dr.z;
    //} while (conf->particleStore[target].pos.x < 0.25 || conf->particleStore[target].pos.x > 0.50);

    reject = 0;
    wlener = 0.0;
    if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wlm[wli]) {
                case 1: origsyscm = conf->syscm;
                    conf->syscm.x += dr.x * topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type].volume / conf->sysvolume;
                    conf->syscm.y += dr.y * topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type].volume / conf->sysvolume;
                    conf->syscm.z += dr.z * topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type].volume / conf->sysvolume;
                    sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                    break;
                case 2: sim->wl.origmesh = sim->wl.mesh;
                   sim->wl.neworder[wli] = meshOrderMoveOne(orig, conf->particleStore[target].pos, &sim->wl.mesh, (long)conf->particleStore.size(), target, wli);
                    break;
                case 4:
                    sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                    break;
                case 5:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    origsyscm = conf->syscm;
                    conf->syscm.x += dr.x * topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type].volume / conf->sysvolume;
                    conf->syscm.y += dr.y * topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type].volume / conf->sysvolume;
                    conf->syscm.z += dr.z * topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type].volume / conf->sysvolume;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    if ( target == 0 )
                      sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->particleStore[0].pos));
                    else
                      sim->wl.neworder[wli] = radiusholeOrderMoveOne(&orig, target,wli,&(conf->particleStore[0].pos));
                    break;
                case 7:
                    sim->wl.partincontactold = sim->wl.partincontact;
                    if ( target == 0 )
                        sim->wl.neworder[wli] = contParticlesAll(wli);
                    else
                        sim->wl.neworder[wli] = contParticlesMoveOne(&orig,target,wli);
                    break;
                default:
                    sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    break;
            }
            if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
            energy += wlener;
        }

    }

    if (!reject) {  /* wang-landaou ok, try move - calcualte energy */
        enermove =  (*calcEnergy)(target, 1, 0);
    }
    if ( reject || moveTry(energy, enermove, sim->temper) ) {  /* probability acceptance */
        conf->particleStore[target].pos = orig;
        sim->trans[conf->particleStore[target].type].rej++;
        if ( (sim->wlm[0] == 1) || (sim->wlm[0] == 5) || (sim->wlm[1] == 1) || (sim->wlm[1] == 5) )
            conf->syscm = origsyscm;
        sim->wl.reject(radiusholemax_orig, sim->wlm);
    } else { /* move was accepted */
        sim->trans[conf->particleStore[target].type].acc++;
        sim->wl.accept(sim->wlm[0]);

        edriftchanges = enermove - energy + wlener;

        //printf("%f\t%f\n", conf->particleStore[0].pos.z * conf->box.z , enermove);
        //printf("%.12f\t%.12f\t%.12f\n", energy , enermove,edriftchanges);
    }

    return edriftchanges;
}

double MoveCreator::partRotate(long target) {
    double edriftchanges,energy,enermove,wlener;
    Particle origpart;
    int reject=0,wli;

    /*=== Rotation step ===*/
    //printf ("rotation %ld npart %ld\n\n",target,npart);
    energy = (*calcEnergy)(target, 1, 0);

    origpart = conf->particleStore[target];

    pscRotate(&conf->particleStore[target],sim->rot[conf->particleStore[target].type].angle, topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type].geotype[0]);
    /*should be normalised and ortogonal but we do for safety*/
    conf->particleStore[target].dir.normalise();
    conf->particleStore[target].patchdir[0].ortogonalise(conf->particleStore[target].dir);

    reject = 0;
    edriftchanges =0.0;
    wlener = 0.0;
    if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wlm[wli]) {
                case 3:
                    if (target == 0)  sim->wl.neworder[wli] = (long) floor( (conf->particleStore[0].dir.z - sim->wl.minorder[wli])/ sim->wl.dorder[wli] );
                    else sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    /* only rotation change direction */
                    break;
                default:
                    sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    break;
            }
            if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
            energy += wlener;
        }
    }

    if (!reject) {  /* wang-landaou ok, try move - calcualte energy */
        enermove =  (*calcEnergy)(target, 1, 0);
    }
    if ( reject || moveTry(energy,enermove,sim->temper) ) {  /* probability acceptance */
        conf->particleStore[target] = origpart;
        sim->rot[conf->particleStore[target].type].rej++;
        sim->wl.reject(sim->wl.radiusholemax, sim->wlm);
    } else { /* move was accepted */
        // DEBUG
        //fprintf(fenergy, "%f\t%f\n", conf->particle[1].pos.x * conf->box.x , enermove);
        sim->rot[conf->particleStore[target].type].acc++;
        sim->wl.accept(sim->wlm[0]);
        edriftchanges = enermove - energy + wlener;
        //printf("%f\t%f\n", conf->particleStore[0].patchdir[0].z, enermove);
    }

    return edriftchanges;
}

void MoveCreator::pscRotate(Particle *psc, double max_angle, int geotype) {
    double vc, vs, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    double d1, d2, d3, d4, d5, d6, d7, d8, d9 , newx, newy, newz;
    int k,m;
    Vector newaxis;

    /* generate quaternion for rotation*/
    newaxis.random(); /*random axes for rotation*/
    //    maxcos = cos(maxorient/2/180*PI);
    // vc = maxcos + ran2(&seed)*(1-maxcos); /*cos of angle must be bigger than maxcos and smaller than one*/

    vc = cos(max_angle * ran2() );

    if (ran2() <0.5) vs = sqrt(1.0 - vc*vc);
    else vs = -sqrt(1.0 - vc*vc); /*randomly choose orientation of direction of rotation clockwise or counterclockwise*/

    Quat newquat(vc, newaxis.x*vs, newaxis.y*vs, newaxis.z*vs);

    /* do quaternion rotation*/
    t2 =  newquat.w * newquat.x;
    t3 =  newquat.w * newquat.y;
    t4 =  newquat.w * newquat.z;
    t5 = -newquat.x * newquat.x;
    t6 =  newquat.x * newquat.y;
    t7 =  newquat.x * newquat.z;
    t8 = -newquat.y * newquat.y;
    t9 =  newquat.y * newquat.z;
    t10 = -newquat.z * newquat.z;

    d1 = t8 + t10;
    d2 = t6 - t4;
    d3 = t3 + t7;
    d4 = t4 + t6;
    d5 = t5 + t10;
    d6 = t9 - t2;
    d7 = t7 - t3;
    d8 = t2 + t9;
    d9 = t5 + t8;

    /*rotate spherocylinder direction vector2*/
    newx = 2.0 * ( d1*psc->dir.x + d2*psc->dir.y + d3*psc->dir.z ) + psc->dir.x;
    newy = 2.0 * ( d4*psc->dir.x + d5*psc->dir.y + d6*psc->dir.z ) + psc->dir.y;
    newz = 2.0 * ( d7*psc->dir.x + d8*psc->dir.y + d9*psc->dir.z ) + psc->dir.z;
    psc->dir.x = newx;
    psc->dir.y = newy;
    psc->dir.z = newz;

    m=1;
    if ( (geotype != SCN) && (geotype != SCA) ) {
        if ( (geotype == TPSC) || (geotype == TCPSC) || (geotype == TCHPSC) || (geotype == TCHCPSC) )
        m=2;
        for (k=0;k<m;k++) {
        /*rotate patch direction vector2*/
        newx = 2.0 * ( d1*psc->patchdir[k].x + d2*psc->patchdir[k].y + d3*psc->patchdir[k].z ) + psc->patchdir[k].x;
        newy = 2.0 * ( d4*psc->patchdir[k].x + d5*psc->patchdir[k].y + d6*psc->patchdir[k].z ) + psc->patchdir[k].y;
        newz = 2.0 * ( d7*psc->patchdir[k].x + d8*psc->patchdir[k].y + d9*psc->patchdir[k].z ) + psc->patchdir[k].z;
        psc->patchdir[k].x = newx;
        psc->patchdir[k].y = newy;
        psc->patchdir[k].z = newz;

        /*rotate patch sides vector2s*/
        newx = 2.0 * ( d1*psc->patchsides[0+2*k].x + d2*psc->patchsides[0+2*k].y + d3*psc->patchsides[0+2*k].z ) + psc->patchsides[0+2*k].x;
        newy = 2.0 * ( d4*psc->patchsides[0+2*k].x + d5*psc->patchsides[0+2*k].y + d6*psc->patchsides[0+2*k].z ) + psc->patchsides[0+2*k].y;
        newz = 2.0 * ( d7*psc->patchsides[0+2*k].x + d8*psc->patchsides[0+2*k].y + d9*psc->patchsides[0+2*k].z ) + psc->patchsides[0+2*k].z;
        psc->patchsides[0+2*k].x = newx;
        psc->patchsides[0+2*k].y = newy;
        psc->patchsides[0+2*k].z = newz;
        newx = 2.0 * ( d1*psc->patchsides[1+2*k].x + d2*psc->patchsides[1+2*k].y + d3*psc->patchsides[1+2*k].z ) + psc->patchsides[1+2*k].x;
        newy = 2.0 * ( d4*psc->patchsides[1+2*k].x + d5*psc->patchsides[1+2*k].y + d6*psc->patchsides[1+2*k].z ) + psc->patchsides[1+2*k].y;
        newz = 2.0 * ( d7*psc->patchsides[1+2*k].x + d8*psc->patchsides[1+2*k].y + d9*psc->patchsides[1+2*k].z ) + psc->patchsides[1+2*k].z;
        psc->patchsides[1+2*k].x = newx;
        psc->patchsides[1+2*k].y = newy;
        psc->patchsides[1+2*k].z = newz;
        }
    }

    m=1;
    if ( (geotype == CHPSC) || (geotype == CHCPSC) || (geotype == TCHPSC) || (geotype == TCHCPSC) ) {
        if ( (geotype == TCHPSC) || (geotype == TCHCPSC) )
        m=2;
        for (k=0;k<m;k++) {
        /*rotate chiral direction vector2*/
        newx = 2.0 * ( d1*psc->chdir[k].x + d2*psc->chdir[k].y + d3*psc->chdir[k].z ) + psc->chdir[k].x;
        newy = 2.0 * ( d4*psc->chdir[k].x + d5*psc->chdir[k].y + d6*psc->chdir[k].z ) + psc->chdir[k].y;
        newz = 2.0 * ( d7*psc->chdir[k].x + d8*psc->chdir[k].y + d9*psc->chdir[k].z ) + psc->chdir[k].z;
        psc->chdir[k].x = newx;
        psc->chdir[k].y = newy;
        psc->chdir[k].z = newz;
        }
    }
}

double MoveCreator::switchTypeMove() {
    double edriftchanges,energy,enermove,wlener;
    int reject=0,wli;
    long target;
    double radiusholemax_orig=0;

    /*=== This is an attempt to switch a type ===*/
    edriftchanges =0.0;
    wlener = 0.0;
    target = ran2() * topo->n_switch_part;
    target = topo->switchlist[target];
    DEBUG_SIM("Switching the particle type");
    DEBUG_SIM("PARTICLE: %ld", target);
    energy = (*calcEnergy)(target, 1, 0);
    // Start switching the type
    int switched = conf->particleStore[target].switched;
    int pmone = PMONE(switched);
    DEBUG_SIM("switched = %d", switched);
    DEBUG_SIM("pmone = %d", pmone);
    int tmp_type = conf->particleStore[target].type;
    conf->particleStore[target].type = conf->particleStore[target].switchtype;
    conf->particleStore[target].switchtype = tmp_type;
    conf->particleStore[target].switched += pmone;
    conf->particleStore[target].init(&(topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type]));
    DEBUG_SIM("Particle %ld is %d switched", target, switched);
    //DEBUG
#ifdef DEBUGGING_SIM
    if ((abs(pmone) != 1) || (conf->particleStore[target].type == conf->particleStore[target].switchtype)){
        fprintf(stderr, "ERROR: Something went wrong, when switching the type of particle %ld\n", target);
        exit(1);
    }
#endif
    if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wlm[wli]) {
                /*case 1: sim->wl.neworder = z_order(&sim->wl, conf,wli);
                     break;*/
                case 2: sim->wl.origmesh = sim->wl.mesh;
                    sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                          (long)conf->particleStore.size(),
                                                                          sim->wl.wlmtype,
                                                                          conf->box,
                                                                          &conf->particleStore) - sim->wl.minorder[wli]);
                    break;
                /*case 4:
                    sim->wl.neworder = twopartdist(&sim->wl,conf,wli);
                    break;*/
                case 5:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->particleStore[0].pos));
                    break;
                case 7:
                    sim->wl.partincontactold = sim->wl.partincontact;
                    sim->wl.neworder[wli] = contParticlesAll(wli);
                    break;
                default:
                    sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    break;
            }
        if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
            energy += wlener;
        }
    }

    if (!reject) {
        enermove = conf->particleStore[target].delta_mu * pmone;
        // DEBUG
        //double dmu = enermove;
        //particleStore[target].switched += pmone;
        enermove += (*calcEnergy)( target, 1, 0);
        //printf("energy: %f \t %f\t%f\n",particleStore[target].delta_mu, dmu, enermove);
    }

    // If not accepted: switch back
    if ( reject || moveTry(energy,enermove,sim->temper) ) {  /* probability acceptance */
        DEBUG_SIM("Did NOT switch it\n");
        conf->particleStore[target].switchtype = conf->particleStore[target].type;
        conf->particleStore[target].type = tmp_type;
        conf->particleStore[target].switched -= pmone;
        conf->particleStore[target].init(&(topo->ia_params[conf->particleStore[target].type][conf->particleStore[target].type]));
        sim->wl.reject(radiusholemax_orig, sim->wlm);
    } else { /* move was accepted */
        sim->wl.accept(sim->wlm[0]);
        edriftchanges = enermove - energy + wlener;
    }

    return edriftchanges;
}

double MoveCreator::chainMove() {
    double edriftchanges =0.0;
    long target;

    /*=== This is a chain move step ===*/
    target = ran2() * conf->chainCount;
    if (ran2() < 0.5) {
        /*=== Displacement step of cluster/chain ===*/
        edriftchanges = chainDisplace(target);
    } else {
        /*=== Rotation step of cluster/chain ===*/
        edriftchanges = chainRotate(target);

    } /* ==== END OF CHAIN MOVES ===== */
    return edriftchanges;
}

double MoveCreator::chainDisplace(long target) {
    double edriftchanges,energy,enermove,wlener;
    Vector dr, origsyscm;
    int reject=0,wli;
    Vector cluscm;
    long current,i;
    Particle chorig[MAXCHL];
    double radiusholemax_orig=0;

    /*=== Displacement step of cluster/chain ===*/
    //printf ("move chain\n\n");
    energy =0.0;
    wlener = 0.0;
    edriftchanges=0.0;
    i=0;
    current = conf->chainlist[target][0];
    cluscm.x = 0;
    cluscm.y = 0;
    cluscm.z = 0;
    origsyscm.x = 0;
    origsyscm.y = 0;
    origsyscm.z = 0;
    while (current >=0 ) {   /* store old configuration calculate energy*/
        chorig[i].pos = conf->particleStore[current].pos;
        energy += (*calcEnergy)(current, 2, target);
        i++;
        current = conf->chainlist[target][i];
    }
    dr.random();
    dr.x *= sim->chainm[conf->particleStore[target].molType].mx/conf->box.x;
    dr.y *= sim->chainm[conf->particleStore[target].molType].mx/conf->box.y;
    dr.z *= sim->chainm[conf->particleStore[target].molType].mx/conf->box.z;
    i=0;
    if ( ((sim->wlm[0] == 3)||(sim->wlm[1] == 3)) && (target == 0) ) {
        dr.z = 0;
        dr.y = 0;
        dr.x = 0;
    }
    current = conf->chainlist[target][0];
    while (current >=0 ) { /* move chaine to new position  */
        if ( (sim->wlm[0] == 1) || (sim->wlm[0] == 5) || (sim->wlm[1] == 1) || (sim->wlm[1] == 5) ) { /* calculate move of center of mass  */
            cluscm.x += dr.x*topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
            cluscm.y += dr.y*topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
            cluscm.z += dr.z*topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
        }
        conf->particleStore[current].pos.x += dr.x;
        conf->particleStore[current].pos.y += dr.y;
        conf->particleStore[current].pos.z += dr.z;
        i++;
        current = conf->chainlist[target][i];
    }
    enermove = 0.0;

    reject = 0;
    if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wlm[wli]) {
                case 1: origsyscm = conf->syscm;
                    conf->syscm.x += cluscm.x / conf->sysvolume;
                    conf->syscm.y += cluscm.y / conf->sysvolume;
                    conf->syscm.z += cluscm.z / conf->sysvolume;
                    sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                    break;
                case 2: sim->wl.origmesh = sim->wl.mesh;
                    sim->wl.neworder[wli] = meshOrderMoveChain(conf->chainlist[target], &sim->wl.mesh, conf->particleStore.size(), chorig,wli);
                    break;
                case 4:
                    sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                    break;
                case 5:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    origsyscm = conf->syscm;
                    conf->syscm.x += cluscm.x / conf->sysvolume;
                    conf->syscm.y += cluscm.y / conf->sysvolume;
                    conf->syscm.z += cluscm.z / conf->sysvolume;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    if ( target == 0 )
                      sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->particleStore[0].pos));
                    else
                      sim->wl.neworder[wli] = radiusholeOrderMoveChain(conf->chainlist[target], chorig,wli,&(conf->particleStore[0].pos));
                    break;
                case 7:
                    sim->wl.partincontactold = sim->wl.partincontact;
                    if ( target == 0 )
                        sim->wl.neworder[wli] = contParticlesAll(wli);
                    else
                        sim->wl.neworder[wli] = contParticlesMoveChain(conf->chainlist[target],chorig,wli);
                    break;
                default:
                    sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    break;
            }
            if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
            energy += wlener;
        }
    }
    if (!reject) { /* wang-landaou ok, try move - calcualte energy */
        i=0;
        current = conf->chainlist[target][0];
        while (current >=0 ) {
            enermove += (*calcEnergy)(current, 2, target);
            i++;
            current = conf->chainlist[target][i];
        }
    }
    if ( reject || moveTry(energy, enermove, sim->temper) ) {  /* probability acceptance */
        i=0;
        current = conf->chainlist[target][0];
        while (current >=0 ) {
            conf->particleStore[current].pos = chorig[i].pos;
            i++;
            current = conf->chainlist[target][i];
        }
        sim->chainm[conf->particleStore[target].molType].rej++;
        if ( (sim->wlm[0] == 1) || (sim->wlm[0] == 5) || (sim->wlm[1] == 1) || (sim->wlm[1] == 5) )
            conf->syscm = origsyscm;
        sim->wl.reject(radiusholemax_orig, sim->wlm);
    } else { /* move was accepted */
        sim->chainm[conf->particleStore[target].molType].acc++;
        sim->wl.accept(sim->wlm[0]);
        edriftchanges = enermove - energy + wlener;
    }

    return edriftchanges;
}

double MoveCreator::chainRotate(long target) {
    double edriftchanges,energy,enermove,wlener;
    int reject=0,wli;
    Vector cluscm;
    double chainvolume;
    long current, i;
    Particle chorig[MAXCHL];
    double radiusholemax_orig=0;

    /*=== Rotation step of cluster/chain ===*/
    //printf ("rotation of chain\n\n");
    energy=0.0; /* set values to zero*/
    edriftchanges=0.0;
    wlener = 0.0;
    current = conf->chainlist[target][0];
    cluscm.x = conf->particleStore[current].pos.x*topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
    cluscm.y = conf->particleStore[current].pos.y*topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
    cluscm.z = conf->particleStore[current].pos.z*topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
    chorig[0] = conf->particleStore[current];
    chainvolume = topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
    energy += (*calcEnergy)(current, 2, target);
    i=1;
    current = conf->chainlist[target][i];
    while (current >=0 ) {   /* store old configuration calculate energy*/
        chorig[i] = conf->particleStore[current];
        /*We have chains whole! don't have to do PBC*/
        /*r_cm.x = conf->particleStore[current].pos.x - conf->particle[first].pos.x;
         r_cm.y = conf->particleStore[current].pos.y - conf->particle[first].pos.y;
         r_cm.z = conf->particleStore[current].pos.z - conf->particle[first].pos.z;
         if ( r_cm.x < 0  )
         r_cm.x -= (double)( (long)(r_cm.x-0.5) );
         else
         r_cm.x -= (double)( (long)(r_cm.x+0.5) );
         if ( r_cm.y < 0  )
         r_cm.y -= (double)( (long)(r_cm.y-0.5) );
         else
         r_cm.y -= (double)( (long)(r_cm.y+0.5) );
         if ( r_cm.z < 0  )
         r_cm.z -= (double)( (long)(r_cm.z-0.5) );
         else
         r_cm.z -= (double)( (long)(r_cm.z+0.5) );
         */
        cluscm.x += conf->particleStore[current].pos.x*topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
        cluscm.y += conf->particleStore[current].pos.y*topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
        cluscm.z += conf->particleStore[current].pos.z*topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
        chainvolume += topo->ia_params[conf->particleStore[current].type][conf->particleStore[current].type].volume;
        energy += (*calcEnergy)(current, 2, target);
        i++;
        current = conf->chainlist[target][i];
    }
    cluscm.x = cluscm.x/chainvolume;
    cluscm.y = cluscm.y/chainvolume;
    cluscm.z = cluscm.z/chainvolume;
    /*do actual rotations around geometrical center*/
    clusterRotate(target, cluscm, sim->chainr[conf->particleStore[target].molType].angle);
    enermove=0.0;

    reject = 0;
    if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wlm[wli]) {
                case 1:
                    if (target == 0) sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                    else sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    /* if we rotated cluster it is around its CM so no change*/
                    break;
                case 2:
                    sim->wl.origmesh = sim->wl.mesh;
                   sim->wl.neworder[wli] = meshOrderMoveChain(conf->chainlist[target], &sim->wl.mesh, conf->particleStore.size(), chorig,wli);
                    break;
                case 3:
                    if (target == 0)  sim->wl.neworder[wli] = (long) floor( (conf->particleStore[0].dir.z - sim->wl.minorder[wli])/ sim->wl.dorder[wli] );
                    else sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    /* only rotation change direction */
                    break;
                case 4:
                    sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                    break;
                case 5:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    if ( target == 0 )
                        sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->particleStore[0].pos));
                    else
                        sim->wl.neworder[wli] = radiusholeOrderMoveChain(conf->chainlist[target], chorig,wli,&(conf->particleStore[0].pos));
                    break;
                case 7:
                    sim->wl.partincontactold = sim->wl.partincontact;
                    if ( target == 0 )
                        sim->wl.neworder[wli] = contParticlesAll(wli);
                    else
                        sim->wl.neworder[wli] = contParticlesMoveChain(conf->chainlist[target],chorig,wli);
                    break;
                default:
                    sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    break;
            }
            if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
            energy += wlener;
        }
    }
    if (!reject) { /* wang-landaou ok, try move - calcualte energy */
        i=0;
        current = conf->chainlist[target][0];
        while (current >=0 ) {
            enermove +=  (*calcEnergy)(current, 2, target);
            i++;
            current = conf->chainlist[target][i];
        }
    }
    if ( reject || moveTry(energy, enermove, sim->temper) ) { /* probability acceptance */
        i=0;
        current = conf->chainlist[target][0];
        while (current >=0 ) {
            conf->particleStore[current] = chorig[i];
            i++;
            current = conf->chainlist[target][i];
        }
        sim->chainr[conf->particleStore[target].molType].rej++;
        sim->wl.reject(radiusholemax_orig, sim->wlm);
    } else { /* move was accepted */
        sim->chainr[conf->particleStore[target].molType].acc++;
        sim->wl.accept(sim->wlm[0]);
        edriftchanges = enermove - energy + wlener;
    }

    return edriftchanges;
}

double MoveCreator::pressureMove() {
    double edriftchanges,energy,enermove,wlener;
    int reject=0,wli;
    double old_side;   /* Box length before attempted change */
    double *side;      /* Box dimension to try changing */
    double psch;       /* Size of a box change during pressure */
    double pvol;       /* Size of a volume during pressure */
    double pvoln;      /* Size of a new volume during pressure */
    double rsave;      /* Saved random number */
    double area;
    double radiusholemax_orig=0;

    /*=== This is a volume change step ===*/
    /*calculate energy*/
    edriftchanges=0.0;
    wlener = 0.0;
    energy = (*calcEnergy)(0, 0, 0);
    /* Choose an edge */
    switch (sim->ptype) {
        case 0:
            /* Anisotropic pressure coupling */
            rsave = ran2();
            if (rsave < 1.0/3.0) {
                side = &(conf->box.x);
                area = conf->box.y * conf->box.z;
            } else if (rsave < 2.0/3.0) {
                side = &(conf->box.y);
                area = conf->box.x * conf->box.z;
            } else {
                side = &(conf->box.z);
                area = conf->box.x * conf->box.y;
            }
            old_side = *side;
            *side += sim->edge.mx * (ran2() - 0.5);

            reject = 0;
            if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
                for (wli=0;wli<sim->wl.wlmdim;wli++) {
                    switch (sim->wlm[wli]) {
                        case 1:
                            sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                            break;
                        case 2:
                            sim->wl.origmesh = sim->wl.mesh;
                            sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                                  conf->particleStore.size(),
                                                                                  sim->wl.wlmtype,
                                                                                  conf->box,
                                                                                  &conf->particleStore) - sim->wl.minorder[wli]);
                            break;
                        case 4:
                            sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                            break;
                        case 5:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->syscm));
                            break;
                        case 6:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->particleStore[0].pos));
                            break;
                        case 7:
                            sim->wl.partincontactold = sim->wl.partincontact;
                            sim->wl.neworder[wli] = contParticlesAll(wli);
                            break;
                        default:
                            sim->wl.neworder[wli] = sim->wl.currorder[wli];
                            break;
                    }
                    if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
                }
                if (!reject) {
                    wlener = sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
                    energy += wlener;
                }
            }
            if (!reject) { /* wang-landaou ok, try move - calculate energy */
                enermove = sim->press * area * (*side - old_side) - (double)conf->particleStore.size() * log(*side/old_side) / sim->temper;
                enermove += (*calcEnergy)(0, 0, 0);
            }
            if ( reject || *side <= 0.0 || ( moveTry(energy,enermove,sim->temper) ) ) { /* probability acceptance */
                *side = old_side;
                sim->edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wlm);
            } else {  /* move was accepted */
                sim->edge.acc++;
                sim->wl.accept(sim->wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;
        case 1:
            /* Isotropic pressure coupling */
            psch = sim->edge.mx * (ran2() - 0.5);
            pvol = conf->box.x * conf->box.y * conf->box.z;
            conf->box.x += psch;
            conf->box.y += psch;
            conf->box.z += psch;
            pvoln = conf->box.x * conf->box.y * conf->box.z;

            reject = 0;
            if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
                for (wli=0;wli<sim->wl.wlmdim;wli++) {
                    switch (sim->wlm[wli]) {
                        case 1: sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                            break;
                        case 2: sim->wl.origmesh = sim->wl.mesh;
                            sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                                  conf->particleStore.size(),
                                                                                  sim->wl.wlmtype,
                                                                                  conf->box,
                                                                                  &conf->particleStore) - sim->wl.minorder[wli]);
                            break;
                        case 4:
                            sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                            break;
                        case 5:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->syscm));
                            break;
                        case 6:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->particleStore[0].pos));
                            break;
                        case 7:
                            sim->wl.partincontactold = sim->wl.partincontact;
                            sim->wl.neworder[wli] = contParticlesAll(wli);
                            break;
                        default:
                            sim->wl.neworder[wli] = sim->wl.currorder[wli];
                            break;
                    }
                    if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
                }
                if (!reject) {
                    wlener = sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
                    energy += wlener;
                }
            }
            if (!reject) { /* wang-landaou ok, try move - calcualte energy */
                enermove = sim->press * (pvoln - pvol) - (double)conf->particleStore.size() * log(pvoln/pvol) / sim->temper;
                enermove += (*calcEnergy)(0, 0, 0);
            }
            if ( reject || moveTry(energy,enermove,sim->temper) )  { /* probability acceptance */
                conf->box.x -= psch;
                conf->box.y -= psch;
                conf->box.z -= psch;
                sim->edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wlm);
            } else { /* move was accepted */
                sim->edge.acc++;
                sim->wl.accept(sim->wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;
        case 2:
            /* Isotropic pressure coupling in xy, z constant */
            psch = sim->edge.mx * (ran2() - 0.5);
            pvol = conf->box.x * conf->box.y;
            conf->box.x += psch;
            conf->box.y += psch;
            pvoln = conf->box.x * conf->box.y;

            reject = 0;
            if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
                for (wli=0;wli<sim->wl.wlmdim;wli++) {
                    switch (sim->wlm[wli]) {
                        /*no change in case 1, it does not change box.z*/
                        case 2: sim->wl.origmesh = sim->wl.mesh;
                            sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                                  conf->particleStore.size(),
                                                                                  sim->wl.wlmtype,
                                                                                  conf->box,
                                                                                  &conf->particleStore) - sim->wl.minorder[wli]);
                            break;
                        case 4:
                            sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                            break;
                        case 5:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->syscm));
                            break;
                        case 6:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->particleStore[0].pos));
                            break;
                        case 7:
                            sim->wl.partincontactold = sim->wl.partincontact;
                            sim->wl.neworder[wli] = contParticlesAll(wli);
                            break;
                        default:
                            sim->wl.neworder[wli] = sim->wl.currorder[wli];
                            break;
                    }
                    if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
                }
                if (!reject) {
                    wlener = sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
                    energy += wlener;
                }
            }
            if (!reject) { /* wang-landaou ok, try move - calculate energy */
                enermove = sim->press * conf->box.z * (pvoln - pvol) - (double)conf->particleStore.size() * log(pvoln/pvol) / sim->temper;
                enermove += (*calcEnergy)(0, 0, 0);
            }
            if ( reject || moveTry(energy,enermove,sim->temper) )  { /* probability acceptance */
                conf->box.x -= psch;
                conf->box.y -= psch;
                sim->edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wlm);
            } else { /* move was accepted */
                sim->edge.acc++;
                sim->wl.accept(sim->wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;
        case 3:
            /* Isotropic pressure coupling in xy, z coupled to have fixed volume */
            psch = sim->edge.mx * (ran2() - 0.5);
            pvol = conf->box.x * conf->box.y * conf->box.z;
            conf->box.x += psch;
            conf->box.y += psch;
            conf->box.z = pvol / conf->box.x / conf->box.y;

            reject = 0;
            if (sim->wlm[0] > 0) {  /* get new neworder for wang-landau */
                for (wli=0;wli<sim->wl.wlmdim;wli++) {
                    switch (sim->wlm[wli]) {
                        case 1: sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                            break;
                        case 2: sim->wl.origmesh = sim->wl.mesh;
                            sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                                  conf->particleStore.size(),
                                                                                  sim->wl.wlmtype,
                                                                                  conf->box,
                                                                                  &conf->particleStore) - sim->wl.minorder[wli]);
                            break;
                        case 4:
                            sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                            break;
                        case 5:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->syscm));
                            break;
                        case 6:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = radiusholeAll(wli,&(conf->particleStore[0].pos));
                            break;
                        case 7:
                            sim->wl.partincontactold = sim->wl.partincontact;
                            sim->wl.neworder[wli] = contParticlesAll(wli);
                            break;
                        default:
                            sim->wl.neworder[wli] = sim->wl.currorder[wli];
                            break;
                    }
                    if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
                }
                if (!reject) {
                    wlener = sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]] - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
                    energy += wlener;
                }
            }
            if (!reject) { /* wang-landaou ok, try move - calculate energy */
                enermove = (*calcEnergy)(0, 0, 0);
            }
            if ( reject || moveTry(energy,enermove,sim->temper) )  { /* probability acceptance */
                conf->box.x -= psch;
                conf->box.y -= psch;
                conf->box.z = pvol / conf->box.x / conf->box.y;
                sim->edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wlm);
            } else { /* move was accepted */
                sim->edge.acc++;
                sim->wl.accept(sim->wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;

        default:
            fprintf (stderr, "ERROR: unknown type of pressure coupling %d",sim->ptype);
            exit(1);
    }

    /*=== End volume change step ===*/
    return edriftchanges;
}

void MoveCreator::clusterRotate(long target, Vector gc, double max_angle) {
    long current,i;
    double vc,vs;
    //double quatsize;

    Vector newaxis;

    // create rotation quaternion
    newaxis.random(); /*random axes for rotation*/
    //    maxcos = cos(maxorient/2/180*PI);
    //vc = maxcos + ran2(&seed)*(1-maxcos); /*cos of angle must be bigger than maxcos and smaller than one*/
    vc = cos(max_angle * ran2() );
    if (ran2() <0.5) vs = sqrt(1.0 - vc*vc);
    else vs = -sqrt(1.0 - vc*vc); /*randomly choose orientation of direction of rotation clockwise or counterclockwise*/

    Quat newquat(vc, newaxis.x*vs, newaxis.y*vs, newaxis.z*vs);

    //quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);

    //shift position to geometrical center
    i=0;
    current = conf->chainlist[target][0];
    while (current >=0 ) {
        //shift position to geometrical center
        conf->particleStore[current].pos.x -= gc.x;
        conf->particleStore[current].pos.y -= gc.y;
        conf->particleStore[current].pos.z -= gc.z;
        //scale things by box not to have them distorted
        conf->particleStore[current].pos.x *= conf->box.x;
        conf->particleStore[current].pos.y *= conf->box.y;
        conf->particleStore[current].pos.z *= conf->box.z;
        //do rotation
        conf->particleStore[current].pos.rotate(newquat);
        conf->particleStore[current].dir.rotate(newquat);
        conf->particleStore[current].patchdir[0].rotate(newquat);
        conf->particleStore[current].patchdir[1].rotate(newquat);
        conf->particleStore[current].chdir[0].rotate(newquat);
        conf->particleStore[current].chdir[1].rotate(newquat);
        conf->particleStore[current].patchsides[0].rotate(newquat);
        conf->particleStore[current].patchsides[1].rotate(newquat);
        conf->particleStore[current].patchsides[2].rotate(newquat);
        conf->particleStore[current].patchsides[3].rotate(newquat);
        //sclae back
        conf->particleStore[current].pos.x /= conf->box.x;
        conf->particleStore[current].pos.y /= conf->box.y;
        conf->particleStore[current].pos.z /= conf->box.z;
        //shift positions back
        conf->particleStore[current].pos.x += gc.x;
        conf->particleStore[current].pos.y += gc.y;
        conf->particleStore[current].pos.z += gc.z;
        i++;
        current = conf->chainlist[target][i];
    }
}

double MoveCreator::replicaExchangeMove(long sweep) {
    double edriftchanges=0.0;
#ifdef ENABLE_MPI
        cout << "replica move" << endl;
        double change, *recwlweights;
        MPI_Status status;
        int oddoreven,count,wli,sizewl = 0;
        MpiExchangeData localmpi,receivedmpi;
        bool reject;
        long localwl,receivedwl;

        //int mpi_newdatatypes();

        //mpi_newdatatypes();
        int i;
        Vector vec;
        Particle part;
        MpiExchangeData exch;
        MPI_Aint     dispstart;

        MPI_Datatype MPI_vector2;
        MPI_Datatype type[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
        int          blocklen[3] = {1, 1, 1};
        MPI_Aint     disp[3];
        MPI_Address( &vec, &dispstart);
        MPI_Address( &(vec.x), &disp[0]);
        MPI_Address( &(vec.y), &disp[1]);
        MPI_Address( &(vec.z), &disp[2]);
        for (i=0; i <3; i++) disp[i] -= dispstart;
        MPI_Type_struct( 3, blocklen, disp, type, &MPI_vector2);
        MPI_Type_commit( &MPI_vector2);

        MPI_Datatype MPI_Particle;
        MPI_Datatype type2[11] = {MPI_vector2,MPI_vector2,MPI_vector2,MPI_vector2,MPI_vector2, MPI_LONG, MPI_LONG, MPI_INT,MPI_INT,MPI_DOUBLE, MPI_INT};
        int          blocklen2[11] = {1, 1, 2,4,2,1,1,1,1,1,1,};
        MPI_Aint     disp2[11];
        MPI_Address( &part, &dispstart);
        MPI_Address( &(part.pos), &disp2[0]);
        MPI_Address( &(part.dir), &disp2[1]);
        MPI_Address( &(part.patchdir), &disp2[2]);
        MPI_Address( &(part.patchsides), &disp2[3]);
        MPI_Address( &(part.chdir), &disp2[4]);
        MPI_Address( &(part.molType), &disp2[5]);
        MPI_Address( &(part.chainIndex), &disp2[6]);
        MPI_Address( &(part.type), &disp2[7]);
        MPI_Address( &(part.switchtype), &disp2[8]);
        MPI_Address( &(part.delta_mu), &disp2[9]);
        MPI_Address( &(part.switched), &disp2[10]);
        for (i=0; i <11; i++) disp2[i] -= dispstart;
        MPI_Type_struct( 11, blocklen2, disp2, type2, &MPI_Particle);
        MPI_Type_commit( &MPI_Particle);

        if (sim->wl.length[1] > 0) {
            sizewl = sim->wl.length[1] * sim->wl.length[0];
        } else {
            sizewl = sim->wl.length[0];
        }
        MPI_Datatype MPI_exchange;
        MPI_Datatype type3[7] = {MPI_vector2, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_vector2, MPI_LONG, MPI_LONG};
        int          blocklen3[7] = {1, 1, 1, 1, 1, 1, 2};
        MPI_Aint     disp3[7];
        MPI_Address( &exch, &dispstart);
        MPI_Address( &(exch.box), &disp3[0]);
        MPI_Address( &(exch.energy), &disp3[1]);
        MPI_Address( &(exch.volume), &disp3[2]);
        MPI_Address( &(exch.accepted), &disp3[3]);
        MPI_Address( &(exch.syscm), &disp3[4]);
        MPI_Address( &(exch.radiusholemax), &disp3[5]);
        MPI_Address( &(exch.wl_order), &disp3[6]);
        for (i=0; i <7; i++) disp3[i] -= dispstart;
        MPI_Type_struct(7, blocklen3, disp3, type3, &MPI_exchange);
        MPI_Type_commit( &MPI_exchange);
        //=== This is an attempt to switch replicas ===

        localmpi.box = conf->box;
        localmpi.energy = (*calcEnergy)(0, 0, 0);
        localmpi.volume = conf->box.x * conf->box.y * conf->box.z;
        localmpi.accepted = 0;
        localmpi.syscm = conf->syscm;
        localmpi.radiusholemax = sim->wl.radiusholemax;
        recwlweights = (double*) malloc( sizeof(double) * sizewl  );
        for (wli=0;wli<2;wli++) {
            localmpi.wl_order[wli] = 0;
            receivedmpi.wl_order[wli] = 0;
        }
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            localmpi.wl_order[wli] = sim->wl.currorder[wli];
            //fprintf(stdout,"wli %d %ld  %ld\n\n", wli, localmpi.wl_order[wli], sim->wl.currorder[wli] );
        }

        if ( (sweep % (2*sim->nrepchange)) == 0)
            // exchange odd ones with even ones
            oddoreven=1;
        else
            // exchange even ones with odd ones
            oddoreven=0;
        if (sim->mpinprocs == 2)
            oddoreven=1;
        count = 1;

        if (sim->mpirank % 2 == oddoreven) {
            if (sim->mpirank > 0) {
                MPI_Send(&localmpi, 1, MPI_exchange, sim->mpirank-1, count, MPI_COMM_WORLD);
                MPI_Send(sim->wl.weights, sizewl, MPI_DOUBLE, sim->mpirank-1, count, MPI_COMM_WORLD);
                //printf("send data: rank: %d energy: %f volume: %f pressure: %f \n",sim->mpirank,localmpi.energy,localmpi.volume,localmpi.pressure);

                MPI_Recv(&receivedmpi, 1, MPI_exchange, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                //decision of accepting or rejecting the exchange was done on other process
                //here we took received configuration (if move was accepted))
                //printf("received data: rank: %d energy: %f volume: %f pressure: %f \n",sim->mpirank,receivedmpi.energy,receivedmpi.volume,receivedmpi.pressure);

                if (receivedmpi.accepted == 1) {
                    sim->mpiexch.acc++;
                    Particle *temppart;
                    temppart = (Particle*) malloc(conf->particleStore.size()*sizeof(Particle));
                    MPI_Recv(temppart, conf->particleStore.size(), MPI_Particle, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
                    //printf("received data: rank: %d\n", sim->mpirank);
                    //printf("part0  x %f y %f z %f\n",temppart[0].pos.x, temppart[0].pos.y, temppart[0].pos.z);
                    //printf("part1  x %f y %f z %f\n",temppart[1].pos.x, temppart[1].pos.y, temppart[1].pos.z);
                    //printf("part0  molType %ld chainn %ld type %d\n",temppart[0].molType,temppart[0].chainn,temppart[0].type);

                    MPI_Send(&conf->particleStore[0], conf->particleStore.size(), MPI_Particle, sim->mpirank-1, count, MPI_COMM_WORLD);
                    //printf("send data: rank: %d\n",sim->mpirank);
                    //printf("part0  x %f y %f z %f\n",conf->particleStore[0].pos.x,conf->particleStore[0].pos.y,conf->particleStore[0].pos.z);
                    //printf("part1  x %f y %f z %f\n",conf->particle[1].pos.x,conf->particle[1].pos.y,conf->particle[1].pos.z);
                    //printf("part0  molType %ld chainn %ld type %d\n",conf->particleStore[0].molType,conf->particleStore[0].chainn,conf->particleStore[0].type);

                    localmpi.accepted = receivedmpi.accepted;
                    conf->box = receivedmpi.box;
                    conf->syscm = receivedmpi.syscm;
                    memcpy(&conf->particleStore[0],temppart,conf->particleStore.size()*sizeof(Particle));
                    edriftchanges = receivedmpi.energy - localmpi.energy;
                    edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)conf->particleStore.size() * log(receivedmpi.volume / localmpi.volume) / sim->temper;
                    if ( sim->wlm[0] >0 ) {
                        for (wli=0;wli<sim->wl.wlmdim;wli++) {
                            sim->wl.neworder[wli] = receivedmpi.wl_order[wli];
                        }
                        sim->wl.accept(sim->wlm[0]);
                        //exchange wl data mesh size and radius hole s
                        for (wli=0;wli<sim->wl.wlmdim;wli++) {
                            switch (sim->wlm[wli]) {
                                case 2:
                                    //it is complicated to send because of different sizes
                                     //we would have to send sizes first and realocate corrrect mesh size and then send data
                                    // it is better to recalculate (a bit slower though)
                                    sim->wl.mesh.meshInit(sim->wl.wl_meshsize,conf->particleStore.size(),sim->wl.wlmtype,
                                                          conf->box, &conf->particleStore);
                                    break;
                                case 5:
                                    //radiushole_all(topo,conf,sim,wli,&(conf->syscm));
                                    sim->wl.radiusholeold = (long*) realloc(sim->wl.radiusholeold,sizeof(long)*receivedmpi.radiusholemax);
                                    MPI_Recv(sim->wl.radiusholeold,receivedmpi.radiusholemax, MPI_LONG, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                                    MPI_Send(sim->wl.radiushole,sim->wl.radiusholemax, MPI_LONG, sim->mpirank-1, count, MPI_COMM_WORLD);
                                    longarrayCpy(&sim->wl.radiushole,&sim->wl.radiusholeold,sim->wl.radiusholemax,receivedmpi.radiusholemax);
                                    sim->wl.radiusholemax=receivedmpi.radiusholemax;
                                    break;
                                case 6:
                                    //radiushole_all(topo,conf,sim,wli,&(conf->particleStore[0].pos));
                                    sim->wl.radiusholeold = (long*) realloc(sim->wl.radiusholeold,sizeof(long)*receivedmpi.radiusholemax);
                                    MPI_Recv(sim->wl.radiusholeold,receivedmpi.radiusholemax, MPI_LONG, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                                    MPI_Send(sim->wl.radiushole,sim->wl.radiusholemax, MPI_LONG, sim->mpirank-1, count, MPI_COMM_WORLD);
                                    longarrayCpy(&sim->wl.radiushole,&sim->wl.radiusholeold,sim->wl.radiusholemax,receivedmpi.radiusholemax);
                                    sim->wl.radiusholemax=receivedmpi.radiusholemax;
                                    break;
                                case 7:
                                    //contparticles_all(topo,conf,sim,wli);
                                    MPI_Recv(&(sim->wl.partincontactold),1, MPI_LONG, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                                    MPI_Send(&(sim->wl.partincontact),1, MPI_LONG, sim->mpirank-1, count, MPI_COMM_WORLD);
                                    sim->wl.partincontact=sim->wl.partincontactold;
                                    break;
                            }
                        }
                    }

                    free(temppart);
                } else {
                    sim->mpiexch.rej++;
                    if ( sim->wlm[0] > 0 ) {
                        sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]] -= sim->wl.alpha;
                        sim->wl.hist[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]]++;
                    }
                }

            }
        } else {
            if (sim->mpirank+1 < sim->mpinprocs) {
                //there is above process
                MPI_Recv(&receivedmpi, 1, MPI_exchange, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(recwlweights, sizewl, MPI_DOUBLE, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                //we got new configuration
                //printf("received data: rank: %d energy: %f volume: %f \n",sim->mpirank,receivedmpi.energy,receivedmpi.volume);

                //valuate if accepte or reject the configuration
                //acc = exp( (1/sim->temper - 1/(sim->temper + sim.dtemp)) * (E_here - E_received) +
                //(sim->press /sim->temper - pressure_received /(sim.temper + sim->dtemp)) * (V_here - V_received)
                //if pressure the same it it simplier
                reject = false;
                change = (1/sim->temper - 1/(sim->temper + sim->dtemp)) * (localmpi.energy - receivedmpi.energy);
                //printf("acceptance decision: change: %f localE: %f receivedE: %f tempf: %f \n",change,localmpi.energy,receivedmpi.energy,(1/sim->temper - 1/(sim->temper + sim->dtemp)));
                change += (sim->press/sim->temper - (sim->press + sim->dpress)/(sim->temper + sim->dtemp)) * (localmpi.volume - receivedmpi.volume);
                //printf("pressf: %f  \n",(sim->press/sim->temper - (sim->press + sim->dpress)/(sim->temper + sim->dtemp)));
                if (sim->wlm[0] > 0) {
                    localwl = sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0];
                    receivedwl = receivedmpi.wl_order[0] + receivedmpi.wl_order[1]*sim->wl.length[0];
                    //fprintf(stdout,"decide wl   %ld %ld %ld energychange: %f \n", receivedmpi.wl_order[0],  receivedmpi.wl_order[1], receivedwl, change );
                    //fprintf(stdout,"local weights %ld %f %ld %f \n",localwl,sim->wl.weights[localwl],receivedwl,sim->wl.weights[receivedwl]);
                    change += (-sim->wl.weights[localwl] + sim->wl.weights[receivedwl] )/sim->temper + ( -recwlweights[receivedwl] + recwlweights[localwl])/(sim->temper + sim->dtemp) ;
                    //fprintf(stdout,"wlchange %f \n\n",change);
                }
                if (  (!(reject)) && ( (change > 0) || (ran2() < exp(change))  )  ) {
                    // Exchange ACCEPTED send local stuff
                    //printf("exchange accepted \n");
                    sim->mpiexch.acc++;
                    localmpi.accepted = 1;
                    conf->box = receivedmpi.box;
                    conf->syscm = receivedmpi.syscm;
                    edriftchanges = receivedmpi.energy - localmpi.energy;
                    edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)conf->particleStore.size() * log(receivedmpi.volume / localmpi.volume) / sim->temper;
                    //printf("edrift %f\n",edriftchanges);
                    if ( sim->wlm[0] > 0 ) {
                        for (wli=0;wli<sim->wl.wlmdim;wli++) {
                            sim->wl.neworder[wli] = receivedmpi.wl_order[wli];
                        }
                        sim->wl.accept(sim->wlm[0]);
                    }
                    MPI_Send(&localmpi, 1, MPI_exchange, sim->mpirank+1, count, MPI_COMM_WORLD);
                    //printf("send data: rank: %d energy: %f volume: %f pressure: %f \n",sim->mpirank,localmpi.energy,localmpi.volume,localmpi.pressure);
                    //send and receive configuration
                    MPI_Send(&conf->particleStore[0], conf->particleStore.size(), MPI_Particle, sim->mpirank+1, count, MPI_COMM_WORLD);
                    //printf("send data: rank: %d\n",sim->mpirank);
                    //printf("part0  x %f y %f z %f\n",conf->particleStore[0].pos.x,conf->particleStore[0].pos.y,conf->particleStore[0].pos.z);
                    //printf("part1  x %f y %f z %f\n",conf->particle[1].pos.x,conf->particle[1].pos.y,conf->particle[1].pos.z);
                    //printf("part0  molType %ld chainn %ld type %d\n",conf->particleStore[0].molType,conf->particleStore[0].chainn,conf->particleStore[0].type);

                    MPI_Recv(&conf->particleStore[0], conf->particleStore.size(), MPI_Particle, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
                    //printf("recieved data: rank: %d\n",sim->mpirank);
                    //printf("part0  x %f y %f z %f\n",conf->particleStore[0].pos.x,conf->particleStore[0].pos.y,conf->particleStore[0].pos.z);
                    //printf("part1  x %f y %f z %f\n",conf->particle[1].pos.x,conf->particle[1].pos.y,conf->particle[1].pos.z);
                    //printf("part0  molType %ld chainn %ld type %d\n",conf->particleStore[0].molType,conf->particleStore[0].chainn,conf->particleStore[0].type);

                    if ( sim->wlm[0] > 0 ) {
                        //exchange wl data mesh size and radius hole s
                        for (wli=0;wli<sim->wl.wlmdim;wli++) {
                            switch (sim->wlm[wli]) {
                                case 2:
                                    //it is complicated to send because of different sizes
                                    //  we would have to send sizes first and realocate corrrect mesh size and then send data
                                    //  it is better to recalculate (a bit slower though)
                                    sim->wl.mesh.meshInit(sim->wl.wl_meshsize, conf->particleStore.size(), sim->wl.wlmtype,
                                                          conf->box, &conf->particleStore);
                                    break;
                                case 5:
                                    //radiushole_all(topo,conf,sim,wli,&(conf->syscm));
                                    sim->wl.radiusholeold = (long*) realloc(sim->wl.radiusholeold,sizeof(long)*receivedmpi.radiusholemax);
                                    MPI_Send(sim->wl.radiushole,sim->wl.radiusholemax, MPI_LONG, sim->mpirank+1, count, MPI_COMM_WORLD);
                                    MPI_Recv(sim->wl.radiusholeold,receivedmpi.radiusholemax, MPI_LONG, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                                    longarrayCpy(&sim->wl.radiushole,&sim->wl.radiusholeold,sim->wl.radiusholemax,receivedmpi.radiusholemax);
                                    sim->wl.radiusholemax=receivedmpi.radiusholemax;
                                    break;
                                case 6:
                                    //radiushole_all(topo,conf,sim,wli,&(conf->particleStore[0].pos));
                                    sim->wl.radiusholeold = (long*) realloc(sim->wl.radiusholeold,sizeof(long)*receivedmpi.radiusholemax);
                                    MPI_Send(sim->wl.radiushole,sim->wl.radiusholemax, MPI_LONG, sim->mpirank+1, count, MPI_COMM_WORLD);
                                    MPI_Recv(sim->wl.radiusholeold,receivedmpi.radiusholemax, MPI_LONG, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                                    longarrayCpy(&sim->wl.radiushole,&sim->wl.radiusholeold,sim->wl.radiusholemax,receivedmpi.radiusholemax);
                                    sim->wl.radiusholemax=receivedmpi.radiusholemax;
                                    break;
                                case 7:
                                    //contparticles_all(topo,conf,sim,wli);
                                    MPI_Send(&(sim->wl.partincontact),1, MPI_LONG, sim->mpirank+1, count, MPI_COMM_WORLD);
                                    MPI_Recv(&(sim->wl.partincontact),1, MPI_LONG, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                                    break;
                            }
                        }
                    }
                } else {
                    //if exchange rejected send back info
                    //printf("exchange rejected\n");
                    sim->mpiexch.rej++;
                    MPI_Send(&localmpi, 1, MPI_exchange, sim->mpirank+1, count, MPI_COMM_WORLD);
                    if ( sim->wlm[0] > 0 ) {
                        sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]] -= sim->wl.alpha;
                        sim->wl.hist[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]]++;
                    }
                }
            }
        }
        if(localmpi.accepted) cout << "mpi accept" << endl;
        //if ( (localmpi.accepted) && (sim->pairlist_update) ) gen pair list

        MPI_Type_free(&MPI_exchange);
        MPI_Type_free(&MPI_Particle);
        MPI_Type_free(&MPI_vector2);
        free(recwlweights);
#endif
    return edriftchanges;
}

double MoveCreator::muVTMove() {

    long target;
    double entrophy = log(conf->box.x * conf->box.y * conf->box.z)/sim->temper;
    double energy = 0.0;

    // Determine what type we will be inserting/deleting
    int molType = getRandomMuVTType();
    assert(topo->chainparam[molType].chemPot != -1 && "chempot uninitialized");

    if(ran2() > 0.5) { //  insert move

        if(topo->chainparam[molType].isAtomic()) {

            // create particle
            insert.clear();
            insert.push_back(Particle());
            insert[0].random(molType, topo->chainparam[molType].particleTypes[0]);
            insert[0].init(topo->ia_params[insert[0].type]);

            // check overlap
            if(conf->overlapAll(&insert[0], topo->ia_params))
                return 0; // overlap detected, move rejected

            // calculate energy, no neighborList used, IMPLEMENT CONLIST
            energy = calcEnergy->oneToAll(&insert[0], -1);

            //
            // accept with probability -> V/N+1 * e^(ln(a*Nav*1e-27)) + (mu - U(new))/kT)
            //
            if( ( (conf->sysvolume / (conf->molCountOfType(molType) + 1)) *
                    exp( (+topo->chainparam[molType].chemPot - energy)/sim->temper) ) > ran2()) {

                conf->addMolecule(&insert);
                insert.clear();
                return energy - entrophy;
                // entire neighborList -> generated in simulate at the end of muVTmove()
            } else { /// rejected
                insert.clear();
                return 0;
            }

        } else { // this is chain insert
            // do energy calc -> neccesity of index during calculation -> to access conectivity list and chainlist
            // chainlist? -> probably only for chain insert
            // chain insert? acceptance ratio? -> from pregenerated conf
            printf("Chain insert not yet implemented!!!\n");
            exit(1);
        }

    } else { // delete move

        // choose particle -> only of certain type -> list of certain types
        target = ran2() * conf->molCountOfType(molType);

        if(topo->chainparam[molType].isAtomic()) {
            // do energy calc
            target = conf->getStoreIndex(molType, target);
            energy = calcEnergy->oneToAll(target);

            //
            // accept with probability -> N/V * e^(3*ln(wavelenght) - mu/kT + U(del)/kT)
            //
            if( ((double)conf->molCountOfType(molType) / conf->sysvolume) *
                   exp( (energy - topo->chainparam[molType].chemPot)/sim->temper ) > ran2()) {

                conf->removeMolecule(target, 1);

                return -energy + entrophy;
            } else {
                return 0;
            }
        } else {
            /// must include inner energy in energy drift
            /// trial move without inner energy
            //energy = (*calcEnergy)(conf->particleStore[target],2,conf->particleStore[target].chainIndex);
            printf("Chain delete Not yet implemented\n");
            exit(1);
        }
    }
    return 0;
}

long MoveCreator::radiusholeAll(int wli, Vector *position) {
    long i,nr,radiusholemax;
    double rx,ry,z;

    radiusholemax = radiusholePosition(sqrt(conf->box.x*conf->box.x+conf->box.y*conf->box.y),wli);
    if ( radiusholemax > sim->wl.radiusholemax ) {
        if (sim->wl.radiushole != NULL)
            free(sim->wl.radiushole);
        sim->wl.radiushole = (long int*) malloc( sizeof(long)* (radiusholemax));
        sim->wl.radiusholemax = radiusholemax;
    }

    for (i=0;i<radiusholemax;i++) {
        sim->wl.radiushole[i] = 0;
    }

    for (i=0; i< (long)conf->particleStore.size(); i++) {
        /*calculate position of particle from z axis, and add it in array */
        if ( conf->particleStore[i].type == sim->wl.wlmtype ) {
            z=conf->particleStore[i].pos.z - (*position).z; /*if above position*/
            if (z-anInt(z) > 0) {
                rx = conf->box.x * (conf->particleStore[i].pos.x - anInt(conf->particleStore[i].pos.x));
                ry = conf->box.y * (conf->particleStore[i].pos.y - anInt(conf->particleStore[i].pos.y));
                nr = radiusholePosition(sqrt(rx*rx+ry*ry), wli);
                if (nr < 0)
                    return -100;
                sim->wl.radiushole[nr]++;
            }
        }
    }

    return radiusholeOrder();
}

long MoveCreator::radiusholeOrder() {
    for (int i=0;i<sim->wl.radiusholemax-3;i++){
        if ((sim->wl.radiushole[i] >0 ) && (sim->wl.radiushole[i+1] >0 ) && (sim->wl.radiushole[i+2] >0 ) && (sim->wl.radiushole[i+3] >0 ))
          return i-1;
    }
    return -100;
}

long MoveCreator::radiusholeOrderMoveOne(Vector *oldpos, long target, int wli, Vector *position) {
    long nr,oor; /* position in radiushole */
    double rx,ry,z;
    bool oz,nz;

    if ( conf->particleStore[target].type != sim->wl.wlmtype )
        return sim->wl.currorder[wli];

    z=conf->particleStore[target].pos.z - position->z; /*if above position*/
    if (z-anInt(z) < 0) nz = false;
    else nz=true;
    z=oldpos->z - position->z;  /*if above position*/
    if (z-anInt(z) < 0) oz = false;
    else oz=true;
    if ( !(nz) && !(oz) )
        return sim->wl.currorder[wli];

    rx = conf->box.x * (conf->particleStore[target].pos.x - anInt(conf->particleStore[target].pos.x));
    ry = conf->box.y * (conf->particleStore[target].pos.y - anInt(conf->particleStore[target].pos.y));
    nr = radiusholePosition(sqrt(rx*rx+ry*ry),wli);

    if (nr < 0)
        return -100;
    /*particle move over radius bins*/
    if (nz) {
        sim->wl.radiushole[nr]++;
    }
    if (oz) {
        rx = conf->box.x * (oldpos->x - anInt(oldpos->x));
        ry = conf->box.y * (oldpos->y - anInt(oldpos->y));
        oor = radiusholePosition(sqrt(rx*rx+ry*ry),wli);
        sim->wl.radiushole[oor]--;
        if ( sim->wl.radiushole[oor] < 0 ) {
            printf ("Error(single particle move): trying to make number of beads in radiuspore smaller than 0 at position %ld\n",oor);
            radiusholePrint(sim->wl.radiushole,sim->wl.radiusholemax);
            fflush(stdout);
        }
        if (sim->wl.radiushole[oor] ==0)
            return radiusholeOrder();
    }

    if ( (nz) && (sim->wl.radiushole[nr] ==1) )  {
        return radiusholeOrder();
    }

    return sim->wl.currorder[wli];
}

long MoveCreator::radiusholeOrderMoveChain(long chain[], Particle chorig[], int wli, Vector *position) {
    long i,current,nr;
    double rx,ry,z;
    bool change=false;

    i = 0;
    rx=0;
    current = chain[0];
    while (current >=0 ) {
        if ( conf->particleStore[current].type == sim->wl.wlmtype ) {
            z=conf->particleStore[current].pos.z - position->z; /*if above system CM*/
            if (z-anInt(z) > 0) {
                rx = conf->box.x * (conf->particleStore[current].pos.x - anInt(conf->particleStore[current].pos.x));
                ry = conf->box.y * (conf->particleStore[current].pos.y - anInt(conf->particleStore[current].pos.y));
                nr = radiusholePosition(sqrt(rx*rx+ry*ry),wli);
                if (nr < 0)
                    return -100;
                sim->wl.radiushole[nr]++;
                if ( sim->wl.radiushole[nr] == 1 ) change = true;
            }
        }
        i++;
        current = chain[i];
    }
    i = 0;
    current = chain[0];
    while (current >=0 )  {
        if ( conf->particleStore[current].type == sim->wl.wlmtype ) {
            z=chorig[i].pos.z - position->z; /*if above system CM*/
            if (z-anInt(z) > 0) {
                rx = conf->box.x * (chorig[i].pos.x - anInt(chorig[i].pos.x));
                ry = conf->box.y * (chorig[i].pos.y - anInt(chorig[i].pos.y));
                nr = radiusholePosition(sqrt(rx*rx+ry*ry),wli);
                sim->wl.radiushole[nr]--;
                if ( sim->wl.radiushole[nr] < 0 ) {
                    printf ("Error (chainmove): trying to make number of beads in radiuspore smaller than 0 at position %ld\n",nr);
                    radiusholePrint(sim->wl.radiushole,sim->wl.radiusholemax);
                    fflush(stdout);
                }
                if ( sim->wl.radiushole[nr] == 0 ) change = true;
            }
        }
        i++;
        current = chain[i];
    }

    if ( change ) {
        return radiusholeOrder();
    }
    return sim->wl.currorder[wli];
}

void MoveCreator::radiusholePrint(long *radiushole, long length) {
    printf("radiushole:\n");
    for (int i=0;i<length;i++) {
        printf("%ld ",radiushole[i]);
    }
    printf("\n");
}

bool MoveCreator::particlesInContact(Vector *vec) {
    double x,y,z;

    x = vec->x - conf->particleStore[0].pos.x;
    y = vec->y - conf->particleStore[0].pos.y;
    z = vec->z - conf->particleStore[0].pos.z;

    x = conf->box.x * (x - anInt(x));
    y = conf->box.y * (y - anInt(y));
    z = conf->box.z * (z - anInt(z));

    if ( x*x + y*y + z*z < WL_CONTACTS) {
        return true;
    }
    else {
        return false;
    }
}

long MoveCreator::contParticlesAll(int wli) {
    sim->wl.partincontact = 0;

    for (int i=1; i< (long)conf->particleStore.size(); i++) {
        /*calculate position of particle and add it if in contact */
        if ( conf->particleStore[i].type == sim->wl.wlmtype ) {
            if ( particlesInContact (&(conf->particleStore[i].pos)) )
                sim->wl.partincontact++;
        }
    }

    return contParticlesOrder(wli);
}

int MoveCreator::getRandomMuVTType() {
    int molType = 0;
    for(int i=0; i<conf->molTypeCount; i++) {
        if(topo->chainparam[i].muVTmove) molType++;
    }
    molType = (long) (ran2() * ((double)molType));
    for(int i=0; i<conf->molTypeCount; i++) {
        if(topo->chainparam[i].muVTmove) {
            if(molType == 0) {
                molType = i;
                break;
            }
            molType--;
        }
    }
    return molType;
}

long MoveCreator::contParticlesMoveOne(Vector *oldpos, long target, int wli) {
    if ( conf->particleStore[target].type != sim->wl.wlmtype )
        return sim->wl.currorder[wli];

    if ( particlesInContact (&(conf->particleStore[target].pos)) )
        sim->wl.partincontact++;
    if ( particlesInContact (oldpos) )
        sim->wl.partincontact--;

    return contParticlesOrder(wli);
}

long MoveCreator::contParticlesMoveChain(long chain[], Particle chorig[], int wli) {
    long i,current;

    i = 0;
    current = chain[0];
    while (current >=0 ) {
        if ( conf->particleStore[current].type == sim->wl.wlmtype ) {
            if ( particlesInContact (&(conf->particleStore[current].pos)) )
                sim->wl.partincontact++;
        }
        i++;
        current = chain[i];
    }
    i = 0;
    current = chain[0];
    while (current >=0 )  {
        if ( conf->particleStore[current].type == sim->wl.wlmtype ) {
            if ( particlesInContact (&(chorig[i].pos)) )
                sim->wl.partincontact--;
        }
        i++;
        current = chain[i];
    }

    return contParticlesOrder(wli);
}

int MoveCreator::moveTry(double energyold, double energynew, double temperature) {
    /*DEBUG   printf ("   Move trial:    %13.8f %13.8f %13.8f %13.8f\n",
      energynew, energyold, temperature, ran2(&seed));*/
    if (energynew <= energyold )
        return 0;
    else if (exp(-1.0*(energynew-energyold)/temperature) > ran2())
        return 0;
    else
        return 1;
}


