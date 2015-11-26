/** @file movecreator.cpp*/

#include "movecreator.h"
#include <iomanip>

#ifdef ENABLE_MPI
# include <mpi.h>
extern MPI_Datatype MPI_vector, MPI_Particle, MPI_exchange;
#endif

double MoveCreator::particleMove() {
    double edriftchanges =0.0;
    long target;

    /*=== This is a particle move step ===*/
    target = ran2() * (long)conf->pvec.size();

    if ( !( ((sim->wl.wlm[0] == 3) || (sim->wl.wlm[1] == 3) ) && (target == 0) ) && \
    ((ran2() < 0.5) || (topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].geotype[0] >= SP)) ) { /* no rotation for spheres */
        //target = 1;
        //printf ("displacement\n\n");
        edriftchanges = partDisplace(target);
    } else {
        /*=== Rotation step ===*/
        // TODO: v pripade Isotropnich kouli nema pohyb ucinost ... mozna dat vyjimku pro koule
        // BTW: partAcialRotate pro uhel 180.0 a pouziti Vector::getRandomUnitConeUniform by se mel chovat stejne jako normalni partRotate ....
        if(sim->coneAngle == 0.0){
            edriftchanges = partRotate(target);
        } else {
            edriftchanges = partAxialRotate(target);
        }

    }
    /*=== End particle move step ===*/
    return edriftchanges;
}

double MoveCreator::printClustersConf() {
    // in cluster when dist < 3
    // Breadth-first search, BFS
    double energy;
    vector<vector<unsigned int> > cluster;


    bool push = true;
    for(unsigned int w=0; w<1; w++) {
        cluster.push_back(vector<unsigned int>());
        cluster.back().push_back(w);
        for(unsigned int i=0; i<cluster.back().size(); i++) {
            for(unsigned int j=0; j<conf->pvec.size(); j++) {
                push = true;
                energy = calcEnergy->p2p(i,j);
                //cout << energy << " ";
                if(energy < -7) {
                    for(unsigned int q=0; q<cluster.back().size(); q++)
                        if(j == cluster.back()[q])
                            push=false;
                    if(push)
                        cluster.back().push_back(j);
                }
            }
        }
        cout << cluster.back().size() << " ";
    }

    for(unsigned int i=0; i< cluster.back().size(); i++) {
        cout << cluster.back()[i] << " ";
    }

    int q;
    for(unsigned int i=0; i<1; i++) {
        FILE* outfile;
        outfile = fopen("cluster", "w");

        cout << cluster[i].size() << endl;
        fprintf (outfile, "%ld\n", (long)cluster[i].size());

        fprintf (outfile, "sweep %ld; box %.10f %.10f %.10f\n", (long)0, conf->geo.box.x, conf->geo.box.y, conf->geo.box.z);

        for(unsigned int j=0; j<cluster[i].size(); j++) {
            q=cluster[i][j];
            fprintf (outfile, "%15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e %d %d\n",
                     conf->geo.box.x * ((conf->pvec[q].pos.x) - anInt(conf->pvec[q].pos.x)),
                     conf->geo.box.y * ((conf->pvec[q].pos.y) - anInt(conf->pvec[q].pos.y)),
                     conf->geo.box.z * ((conf->pvec[q].pos.z) - anInt(conf->pvec[q].pos.z)),
                     conf->pvec[q].dir.x, conf->pvec[q].dir.y, conf->pvec[q].dir.z,
                     conf->pvec[q].patchdir[0].x, conf->pvec[q].patchdir[0].y, conf->pvec[q].patchdir[0].z,
                    conf->pvec[q].switched,
                    conf->pvec[q].molType);
        }
        fclose(outfile);
    }

    return 0.0;
}

double MoveCreator::partDisplace(long target) {
    double edriftchanges,energy,enermove,wlener;
    Vector orig, dr, origsyscm;
    int reject=0,wli;
    double radiusholemax_orig=0;
    //size_t temp;

    /*=== Displacement step ===*/
    edriftchanges =0.0;
    origsyscm.x = 0;
    origsyscm.y = 0;
    origsyscm.z = 0;

    //temp = clock();
    energy = (*calcEnergy)(target, 1, 0);
    //sim->energyCalc += clock() - temp;

    orig = conf->pvec[target].pos;

    dr.randomUnitSphere();

    dr.x *= sim->stat.trans[conf->pvec[target].type].mx/conf->geo.box.x;
    dr.y *= sim->stat.trans[conf->pvec[target].type].mx/conf->geo.box.y;
    dr.z *= sim->stat.trans[conf->pvec[target].type].mx/conf->geo.box.z;
    if ( ((sim->wl.wlm[0] == 3)||(sim->wl.wlm[1] == 3)) && (target == 0) ) {
        dr.z = 0;
        dr.y = 0;
        dr.x = 0;
    }
    conf->pvec[target].pos.x += dr.x;
    conf->pvec[target].pos.y += dr.y;
    conf->pvec[target].pos.z += dr.z;
    //} while (conf->pvec[target].pos.x < 0.25 || conf->pvec[target].pos.x > 0.50);

#ifdef WEDGE
        conf->geo.usePBC(&conf->pvec[target]);
#endif

    reject = 0;
    wlener = 0.0;
    if (sim->wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wl.wlm[wli]) {
                case 1: origsyscm = conf->syscm;
                    conf->syscm.x += dr.x * topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].volume / conf->sysvolume;
                    conf->syscm.y += dr.y * topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].volume / conf->sysvolume;
                    conf->syscm.z += dr.z * topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].volume / conf->sysvolume;
                    sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                    break;
                case 2: sim->wl.origmesh = sim->wl.mesh;
                   sim->wl.neworder[wli] = meshOrderMoveOne(orig, conf->pvec[target].pos, &sim->wl.mesh, (long)conf->pvec.size(), target, wli);
                    break;
                case 4:
                    sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                    break;
                case 5:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    origsyscm = conf->syscm;
                    conf->syscm.x += dr.x * topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].volume / conf->sysvolume;
                    conf->syscm.y += dr.y * topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].volume / conf->sysvolume;
                    conf->syscm.z += dr.z * topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].volume / conf->sysvolume;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    if ( target == 0 )
                      sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->pvec[0].pos));
                    else
                      sim->wl.neworder[wli] = sim->wl.radiusholeOrderMoveOne(&orig, target,wli,&(conf->pvec[0].pos));
                    break;
                case 7:
                    sim->wl.partincontactold = sim->wl.partincontact;
                    if ( target == 0 )
                        sim->wl.neworder[wli] = sim->wl.contParticlesAll(wli);
                    else
                        sim->wl.neworder[wli] = sim->wl.contParticlesMoveOne(&orig,target,wli);
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

#ifdef WEDGE
    if(conf->geo.boundaryOverlap(&conf->pvec[target].pos)) { // reject solid overlap, soft already handled by usePBC()
        conf->pvec[target].pos = orig;
        sim->stat.trans[conf->pvec[target].type].rej++;
        if ( (sim->wl.wlm[0] == 1) || (sim->wl.wlm[0] == 5) || (sim->wl.wlm[1] == 1) || (sim->wl.wlm[1] == 5) )
            conf->syscm = origsyscm;
        sim->wl.reject(radiusholemax_orig, sim->wl.wlm);

        return 0.0;
    }
#endif

    if (!reject) {  /* wang-landaou ok, try move - calcualte energy */
        //temp = clock();
        enermove =  (*calcEnergy)(target, 1, 0);
        //sim->energyCalc += clock() - temp;
    }
    if ( reject || moveTry(energy, enermove, sim->temper) ) {  /* probability acceptance */
        conf->pvec[target].pos = orig;
        sim->stat.trans[conf->pvec[target].type].rej++;
        if ( (sim->wl.wlm[0] == 1) || (sim->wl.wlm[0] == 5) || (sim->wl.wlm[1] == 1) || (sim->wl.wlm[1] == 5) )
            conf->syscm = origsyscm;
        sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
    } else { /* move was accepted */
        sim->stat.trans[conf->pvec[target].type].acc++;
        sim->wl.accept(sim->wl.wlm[0]);

        edriftchanges = enermove - energy + wlener;

        //printf("%f\t%f\n", conf->pvec[0].pos.z * conf->geo.box.z , enermove);
        //printf("%.12f\t%.12f\t%.12f\n", energy , enermove,edriftchanges);
    }

    return edriftchanges;
}

double MoveCreator::partRotate(long target) {
    double edriftchanges,energy,enermove,wlener;
    Particle origpart;
    int reject=0,wli;
    //size_t temp;

    /*=== Rotation step ===*/
    //printf ("rotation %ld npart %ld\n\n",target,npart);
    //temp = clock();
    energy = (*calcEnergy)(target, 1, 0);
    //sim->energyCalc += clock() - temp;

    origpart = conf->pvec[target];
//    pscRotate(&conf->pvec[target], sim->stat.rot[conf->pvec[target].type].angle, topo.ia_params[origpart.type][origpart.type].geotype[0]);
    conf->pvec[target].rotateRandom(sim->stat.rot[conf->pvec[target].type].angle, topo.ia_params[origpart.type][origpart.type].geotype[0]);

    /*should be normalised and ortogonal but we do for safety*/
    conf->pvec[target].dir.normalise();
    conf->pvec[target].patchdir[0].ortogonalise(conf->pvec[target].dir);

    reject = 0;
    edriftchanges =0.0;
    wlener = 0.0;
    if (sim->wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wl.wlm[wli]) {
                case 3:
                    if (target == 0)  sim->wl.neworder[wli] = (long) floor( (conf->pvec[0].dir.z - sim->wl.minorder[wli])/ sim->wl.dorder[wli] );
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
        //temp = clock();
        enermove =  (*calcEnergy)(target, 1, 0);
        //sim->energyCalc += clock() - temp;
    }
    if ( reject || moveTry(energy,enermove,sim->temper) ) {  /* probability acceptance */
        conf->pvec[target] = origpart;
        sim->stat.rot[conf->pvec[target].type].rej++;
        sim->wl.reject(sim->wl.radiusholemax, sim->wl.wlm);
    } else { /* move was accepted */
        // DEBUG
        //fprintf(fenergy, "%f\t%f\n", conf->particle[1].pos.x * conf->geo.box.x , enermove);
        sim->stat.rot[conf->pvec[target].type].acc++;
        sim->wl.accept(sim->wl.wlm[0]);
        edriftchanges = enermove - energy + wlener;
        //printf("%f\t%f\n", conf->pvec[0].patchdir[0].z, enermove);
    }

    return edriftchanges;
}

double MoveCreator::partAxialRotate(long target){
    double   edriftchanges   =   0.0,
             energyold       =   (*calcEnergy)(target, 1, 0),
             energynew       =   0.0;

    Vector   rotaxis;

    Particle origpart        =   conf->pvec[target];

    //=============================================//
    //            Get vector from cone             //
    //=============================================//
    // Get vector which is randomly distributed in cone around patch direction. Cone is specified by angle in radians in options coneAngle
    rotaxis = Vector::getRandomUnitConeUniform( conf->pvec[target].dir,\
                                                sim->coneAngle);

    //=============================================//
    //              Rotate particle                //
    //=============================================//
    // Now rotate particle around rotaxis in specified cone around patch direction
    conf->pvec[target].pscRotate(   sim->stat.rot[conf->pvec[target].type].angle*ran2(),\
                                    topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].geotype[0],\
                                    rotaxis);

    //=============================================//
    //                MC criterium                 //
    //=============================================//
    energynew = (*calcEnergy)(target, 1, 0); // Calculate energy change of target with rest of system

    if (moveTry(energyold, energynew, sim->temper)){
        // move was rejected
        conf->pvec[target] = origpart; // return to old configuration
    } else {
        // move was accepted
        edriftchanges = energynew - energyold;
    }
    return edriftchanges;
}


double MoveCreator::switchTypeMove() {
    double edriftchanges=0.0, energy,enermove,wlener=0.0;
    int reject=0,wli;
    long target;
    double radiusholemax_orig=0;
    int switchType, sequence_num=0, delta_mu;

    /*=== This is an attempt to switch a type ===*/   
    target = ran2() * conf->pvec.switchPartCount();
    target = conf->pvec.getSwitchPart(target, sequence_num); // stores sequence number
    delta_mu = topo.moleculeParam[conf->pvec[target].molType ].deltaMu[sequence_num];
    if(conf->pvec[target].switched == 0)
        switchType = topo.moleculeParam[conf->pvec[target].molType ].switchTypes[sequence_num];
    else
        switchType = topo.moleculeParam[conf->pvec[target].molType ].particleTypes[sequence_num];

    DEBUG_SIM("Switching the particle type");
    DEBUG_SIM("PARTICLE: %ld", target);
    energy = (*calcEnergy)(target, 1, 0);
    // Start switching the type
    int switched = conf->pvec[target].switched;
    int pmone = PMONE(switched);
    DEBUG_SIM("switched = %d", switched);
    DEBUG_SIM("pmone = %d", pmone);
    int tmp_type = conf->pvec[target].type;
    conf->pvec[target].type = switchType;//conf->pvec[target].switchtype;
    /*conf->pvec[target].switchtype*/ switchType = tmp_type;
    conf->pvec[target].switched += pmone;
    conf->pvec[target].init(&(topo.ia_params[conf->pvec[target].type][conf->pvec[target].type]));
    DEBUG_SIM("Particle %ld is %d switched", target, switched);
    //DEBUG
#ifdef DEBUGGING_SIM
    if ((abs(pmone) != 1) || (conf->pvec[target].type == conf->pvec[target].switchtype)){
        fprintf(stderr, "ERROR: Something went wrong, when switching the type of particle %ld\n", target);
        exit(1);
    }
#endif
    if (sim->wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wl.wlm[wli]) {
                /*case 1: sim->wl.neworder = z_order(&sim->wl, conf,wli);
                     break;*/
                case 2: sim->wl.origmesh = sim->wl.mesh;
                    sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                          (long)conf->pvec.size(),
                                                                          sim->wl.wlmtype,
                                                                          conf->geo.box,
                                                                          &conf->pvec) - sim->wl.minorder[wli]);
                    break;
                /*case 4:
                    sim->wl.neworder = twopartdist(&sim->wl,conf,wli);
                    break;*/
                case 5:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->pvec[0].pos));
                    break;
                case 7:
                    sim->wl.partincontactold = sim->wl.partincontact;
                    sim->wl.neworder[wli] = sim->wl.contParticlesAll(wli);
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
        enermove = delta_mu * pmone;
        // DEBUG
        //double dmu = enermove;
        //pvec[target].switched += pmone;
        enermove += (*calcEnergy)( target, 1, 0);
        //printf("energy: %f \t %f\t%f\n",pvec[target].delta_mu, dmu, enermove);
    }

    // If not accepted: switch back
    if ( reject || moveTry(energy,enermove,sim->temper) ) {  /* probability acceptance */
        DEBUG_SIM("Did NOT switch it\n");
        conf->pvec[target].type = tmp_type;
        conf->pvec[target].switched -= pmone;
        conf->pvec[target].init(&(topo.ia_params[conf->pvec[target].type][conf->pvec[target].type]));
        sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
    } else { /* move was accepted */
        sim->wl.accept(sim->wl.wlm[0]);
        edriftchanges = enermove - energy + wlener;
    }
    return edriftchanges;
}

double MoveCreator::chainMove() {

    double edriftchanges =0.0;
    long target;

    if(conf->pvec.getChainCount() == 0) // no chains to displace - muVTmove deleted all
        return 0.0;

    //=== This is a chain move step ===
    target = ran2() * conf->pvec.getChainCount();

    if (ran2() < 0.5) {
        /*=== Displacement step of cluster/chain ===*/
        edriftchanges = chainDisplace(target);
    } else {
        /*=== Rotation step of cluster/chain ===*/
        edriftchanges = chainRotate(target);
    } /* ==== END OF CHAIN MOVES ===== */
    return edriftchanges;
}

double MoveCreator::chainDisplace(long target)
{
    Molecule chain = conf->pvec.getChain(target);
    double edriftchanges=0.0,energy=0.0,enermove=0.0,wlener=0.0;
    Vector dr, origsyscm(0.0, 0.0, 0.0);
    int reject=0,wli;
    Vector cluscm(0.0, 0.0, 0.0);
    Particle chorig[MAXCHL];
    double radiusholemax_orig=0.0;

    /*=== Displacement step of cluster/chain ===*/
    //printf ("move chain\n\n");
    for(unsigned int i=0; i<chain.size(); i++) // store old configuration
        chorig[i].pos = conf->pvec[chain[i]].pos;

    energy += calcEnergy->mol2others(chain); // Inner energy indiferent to displacement of chain

    dr.randomUnitSphere();
    dr.x *= sim->stat.chainm[conf->pvec[chain[0]].molType].mx/conf->geo.box.x;
    dr.y *= sim->stat.chainm[conf->pvec[chain[0]].molType].mx/conf->geo.box.y;
    dr.z *= sim->stat.chainm[conf->pvec[chain[0]].molType].mx/conf->geo.box.z;

    if ( ((sim->wl.wlm[0] == 3)||(sim->wl.wlm[1] == 3)) && (target == 0) ) {
        dr.z = 0;
        dr.y = 0;
        dr.x = 0;
    }
    for(unsigned int j=0; j<chain.size(); j++) { // move chaine to new position
        if ( (sim->wl.wlm[0] == 1) || (sim->wl.wlm[0] == 5) || (sim->wl.wlm[1] == 1) || (sim->wl.wlm[1] == 5) ) { /* calculate move of center of mass  */
            cluscm.x += dr.x*topo.ia_params[conf->pvec[chain[j]].type][conf->pvec[chain[j]].type].volume;
            cluscm.y += dr.y*topo.ia_params[conf->pvec[chain[j]].type][conf->pvec[chain[j]].type].volume;
            cluscm.z += dr.z*topo.ia_params[conf->pvec[chain[j]].type][conf->pvec[chain[j]].type].volume;
        }
        conf->pvec[chain[j]].pos.x += dr.x;
        conf->pvec[chain[j]].pos.y += dr.y;
        conf->pvec[chain[j]].pos.z += dr.z;
    }

    reject = 0;
    if (sim->wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wl.wlm[wli]) {
                case 1: origsyscm = conf->syscm;
                    conf->syscm.x += cluscm.x / conf->sysvolume;
                    conf->syscm.y += cluscm.y / conf->sysvolume;
                    conf->syscm.z += cluscm.z / conf->sysvolume;
                    sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                    break;
                case 2:
                    sim->wl.origmesh = sim->wl.mesh;
                    sim->wl.neworder[wli] = meshOrderMoveChain(chain, &sim->wl.mesh, conf->pvec.size(), chorig,wli);
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
                    sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    if ( target == 0 )
                      sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->pvec[0].pos));
                    else
                      sim->wl.neworder[wli] = sim->wl.radiusholeOrderMoveChain(chain, chorig,wli,&(conf->pvec[0].pos));
                    break;
                case 7:
                    sim->wl.partincontactold = sim->wl.partincontact;
                    if ( target == 0 )
                        sim->wl.neworder[wli] = sim->wl.contParticlesAll(wli);
                    else
                        sim->wl.neworder[wli] = sim->wl.contParticlesMoveChain(chain,chorig,wli);
                    break;
                default:
                    sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    break;
            }
            if ( (sim->wl.neworder[wli] < 0) || (sim->wl.neworder[wli] >= sim->wl.length[wli]) ) reject = 1;
        }
        if (!reject) {
            wlener += sim->wl.weights[sim->wl.neworder[0]+sim->wl.neworder[1]*sim->wl.length[0]]
                    - sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]];
            energy += wlener;
        }
    }
    if (!reject) { /* wang-landaou ok, try move - calcualte energy */
        enermove += calcEnergy->mol2others(chain);
    }
    if ( reject || moveTry(energy, enermove, sim->temper) ) {  // probability acceptance
        for(unsigned int j=0; j<chain.size(); j++)
            conf->pvec[chain[j]].pos = chorig[j].pos;

        sim->stat.chainm[conf->pvec[chain[0]].molType].rej++;
        if ( (sim->wl.wlm[0] == 1) || (sim->wl.wlm[0] == 5) || (sim->wl.wlm[1] == 1) || (sim->wl.wlm[1] == 5) )
            conf->syscm = origsyscm;
        sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
    } else { /* move was accepted */
        sim->stat.chainm[conf->pvec[chain[0]].molType].acc++;
        sim->wl.accept(sim->wl.wlm[0]);
        edriftchanges = enermove - energy + wlener;
    }

    return edriftchanges;
}

double MoveCreator::chainRotate(long target) {
    Molecule chain = conf->pvec.getChain(target);
    double edriftchanges=0.0, energy=0.0, enermove=0.0, wlener=0.0;
    int reject=0;
    Particle chorig[MAXCHL];
    double radiusholemax_orig=0;

    /*=== Rotation step of cluster/chain ===*/
    //printf ("rotation of chain\n\n");
    for(unsigned int j=0; j<chain.size(); j++) { // store old configuration calculate energy
        chorig[j] = conf->pvec[chain[j]];
        /*We have chains whole! don't have to do PBC*/
        /*r_cm.x = conf->pvec[current].pos.x - conf->particle[first].pos.x;
         r_cm.y = conf->pvec[current].pos.y - conf->particle[first].pos.y;
         r_cm.z = conf->pvec[current].pos.z - conf->particle[first].pos.z;
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
    }

    energy += calcEnergy->mol2others(chain);

    //do actual rotations around geometrical center
    clusterRotate(chain, sim->stat.chainr[conf->pvec[chain[0]].molType].angle);

    if (sim->wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
        for (int wli=0;wli<sim->wl.wlmdim;wli++) {
            switch (sim->wl.wlm[wli]) {
                case 1:
                    if (target == 0) sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                    else sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    /* if we rotated cluster it is around its CM so no change*/
                    break;
                case 2:
                    sim->wl.origmesh = sim->wl.mesh;
                   sim->wl.neworder[wli] = meshOrderMoveChain(chain, &sim->wl.mesh, conf->pvec.size(), chorig,wli);
                    break;
                case 3:
                    if (target == 0)  sim->wl.neworder[wli] = (long) floor( (conf->pvec[0].dir.z - sim->wl.minorder[wli])/ sim->wl.dorder[wli] );
                    else sim->wl.neworder[wli] = sim->wl.currorder[wli];
                    /* only rotation change direction */
                    break;
                case 4:
                    sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                    break;
                case 5:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    radiusholemax_orig = sim->wl.radiusholemax;
                    longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                    if ( target == 0 )
                        sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->pvec[0].pos));
                    else
                        sim->wl.neworder[wli] = sim->wl.radiusholeOrderMoveChain(chain, chorig,wli,&(conf->pvec[0].pos));
                    break;
                case 7:
                    sim->wl.partincontactold = sim->wl.partincontact;
                    if ( target == 0 )
                        sim->wl.neworder[wli] = sim->wl.contParticlesAll(wli);
                    else
                        sim->wl.neworder[wli] = sim->wl.contParticlesMoveChain(chain,chorig,wli);
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
    if (!reject) { // wang-landaou ok, try move - calcualte energy
        enermove += calcEnergy->mol2others(chain);
    }
    if ( reject || moveTry(energy, enermove, sim->temper) ) { // probability acceptance
        for(unsigned int j=0; j<chain.size(); j++)
            conf->pvec[chain[j]] = chorig[j];

        sim->stat.chainr[conf->pvec[chain[0]].molType].rej++;
        sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
    } else { // move was accepted
        sim->stat.chainr[conf->pvec[chain[0]].molType].acc++;
        sim->wl.accept(sim->wl.wlm[0]);
        edriftchanges = enermove - energy + wlener;
    }

    return edriftchanges;
}

double MoveCreator::pressureMove() {
    double edriftchanges,energy,enermove,wlener;
    int reject=0,wli;
    double old_side;   /* geo.box length before attempted change */
    double *side;      /* geo.box dimension to try changing */
    double psch;       /* Size of a geo.box change during pressure */
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
                side = &(conf->geo.box.x);
                area = conf->geo.box.y * conf->geo.box.z;
            } else if (rsave < 2.0/3.0) {
                side = &(conf->geo.box.y);
                area = conf->geo.box.x * conf->geo.box.z;
            } else {
                side = &(conf->geo.box.z);
                area = conf->geo.box.x * conf->geo.box.y;
            }
            old_side = *side;
            *side += sim->stat.edge.mx * (ran2() - 0.5);

            reject = 0;
            if (sim->wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
                for (wli=0;wli<sim->wl.wlmdim;wli++) {
                    switch (sim->wl.wlm[wli]) {
                        case 1:
                            sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                            break;
                        case 2:
                            sim->wl.origmesh = sim->wl.mesh;
                            sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                                  conf->pvec.size(),
                                                                                  sim->wl.wlmtype,
                                                                                  conf->geo.box,
                                                                                  &conf->pvec) - sim->wl.minorder[wli]);
                            break;
                        case 4:
                            sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                            break;
                        case 5:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->syscm));
                            break;
                        case 6:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->pvec[0].pos));
                            break;
                        case 7:
                            sim->wl.partincontactold = sim->wl.partincontact;
                            sim->wl.neworder[wli] = sim->wl.contParticlesAll(wli);
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
                enermove = sim->press * area * (*side - old_side) - (double)conf->pvec.size() * log(*side/old_side) / sim->temper;
                enermove += (*calcEnergy)(0, 0, 0);
            }
            if ( reject || *side <= 0.0 || ( moveTry(energy,enermove,sim->temper) ) ) { /* probability acceptance */
                *side = old_side;
                sim->stat.edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
            } else {  /* move was accepted */
                sim->stat.edge.acc++;
                sim->wl.accept(sim->wl.wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;
        case 1:
            /* Isotropic pressure coupling */
            psch = sim->stat.edge.mx * (ran2() - 0.5);
            pvol = conf->geo.box.x * conf->geo.box.y * conf->geo.box.z;
            conf->geo.box.x += psch;
            conf->geo.box.y += psch;
            conf->geo.box.z += psch;
            pvoln = conf->geo.box.x * conf->geo.box.y * conf->geo.box.z;

            reject = 0;
            if (sim->wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
                for (wli=0;wli<sim->wl.wlmdim;wli++) {
                    switch (sim->wl.wlm[wli]) {
                        case 1: sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                            break;
                        case 2: sim->wl.origmesh = sim->wl.mesh;
                            sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                                  conf->pvec.size(),
                                                                                  sim->wl.wlmtype,
                                                                                  conf->geo.box,
                                                                                  &conf->pvec) - sim->wl.minorder[wli]);
                            break;
                        case 4:
                            sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                            break;
                        case 5:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->syscm));
                            break;
                        case 6:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->pvec[0].pos));
                            break;
                        case 7:
                            sim->wl.partincontactold = sim->wl.partincontact;
                            sim->wl.neworder[wli] = sim->wl.contParticlesAll(wli);
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
                enermove = sim->press * (pvoln - pvol) - (double)conf->pvec.size() * log(pvoln/pvol) / sim->temper;
                enermove += (*calcEnergy)(0, 0, 0);
            }
            if ( reject || moveTry(energy,enermove,sim->temper) )  { /* probability acceptance */
                conf->geo.box.x -= psch;
                conf->geo.box.y -= psch;
                conf->geo.box.z -= psch;
                sim->stat.edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
            } else { /* move was accepted */
                sim->stat.edge.acc++;
                sim->wl.accept(sim->wl.wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;
        case 2:
            /* Isotropic pressure coupling in xy, z constant */
            psch = sim->stat.edge.mx * (ran2() - 0.5);
            pvol = conf->geo.box.x * conf->geo.box.y;
            conf->geo.box.x += psch;
            conf->geo.box.y += psch;
            pvoln = conf->geo.box.x * conf->geo.box.y;

            reject = 0;
            if (sim->wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
                for (wli=0;wli<sim->wl.wlmdim;wli++) {
                    switch (sim->wl.wlm[wli]) {
                        /*no change in case 1, it does not change geo.box.z*/
                        case 2: sim->wl.origmesh = sim->wl.mesh;
                            sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                                  conf->pvec.size(),
                                                                                  sim->wl.wlmtype,
                                                                                  conf->geo.box,
                                                                                  &conf->pvec) - sim->wl.minorder[wli]);
                            break;
                        case 4:
                            sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                            break;
                        case 5:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->syscm));
                            break;
                        case 6:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->pvec[0].pos));
                            break;
                        case 7:
                            sim->wl.partincontactold = sim->wl.partincontact;
                            sim->wl.neworder[wli] = sim->wl.contParticlesAll(wli);
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
                enermove = sim->press * conf->geo.box.z * (pvoln - pvol) - (double)conf->pvec.size() * log(pvoln/pvol) / sim->temper;
                enermove += (*calcEnergy)(0, 0, 0);
            }
            if ( reject || moveTry(energy,enermove,sim->temper) )  { /* probability acceptance */
                conf->geo.box.x -= psch;
                conf->geo.box.y -= psch;
                sim->stat.edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
            } else { /* move was accepted */
                sim->stat.edge.acc++;
                sim->wl.accept(sim->wl.wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;
        case 3:
            /* Isotropic pressure coupling in xy, z coupled to have fixed volume */
            psch = sim->stat.edge.mx * (ran2() - 0.5);
            pvol = conf->geo.box.x * conf->geo.box.y * conf->geo.box.z;
            conf->geo.box.x += psch;
            conf->geo.box.y += psch;
            conf->geo.box.z = pvol / conf->geo.box.x / conf->geo.box.y;

            reject = 0;
            if (sim->wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
                for (wli=0;wli<sim->wl.wlmdim;wli++) {
                    switch (sim->wl.wlm[wli]) {
                        case 1: sim->wl.neworder[wli] = sim->wl.zOrder(wli);
                            break;
                        case 2: sim->wl.origmesh = sim->wl.mesh;
                            sim->wl.neworder[wli] = (long) (sim->wl.mesh.meshInit(sim->wl.wl_meshsize,
                                                                                  conf->pvec.size(),
                                                                                  sim->wl.wlmtype,
                                                                                  conf->geo.box,
                                                                                  &conf->pvec) - sim->wl.minorder[wli]);
                            break;
                        case 4:
                            sim->wl.neworder[wli] = sim->wl.twoPartDist(wli);
                            break;
                        case 5:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->syscm));
                            break;
                        case 6:
                            radiusholemax_orig = sim->wl.radiusholemax;
                            longarrayCpy(&sim->wl.radiusholeold,&sim->wl.radiushole,sim->wl.radiusholemax,sim->wl.radiusholemax);
                            sim->wl.neworder[wli] = sim->wl.radiusholeAll(wli,&(conf->pvec[0].pos));
                            break;
                        case 7:
                            sim->wl.partincontactold = sim->wl.partincontact;
                            sim->wl.neworder[wli] = sim->wl.contParticlesAll(wli);
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
                conf->geo.box.x -= psch;
                conf->geo.box.y -= psch;
                conf->geo.box.z = pvol / conf->geo.box.x / conf->geo.box.y;
                sim->stat.edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
            } else { /* move was accepted */
                sim->stat.edge.acc++;
                sim->wl.accept(sim->wl.wlm[0]);
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



void MoveCreator::clusterRotate(vector<int> &cluster, double max_angle) {
    Vector cluscm;
    double vc,vs;
    Vector newaxis;

    cluscm = clusterCM(cluster);

    // create rotation quaternion
    newaxis.randomUnitSphere(); /*random axes for rotation*/
    vc = cos(max_angle * ran2() );
    if (ran2() <0.5) vs = sqrt(1.0 - vc*vc);
    else vs = -sqrt(1.0 - vc*vc); /*randomly choose orientation of direction of rotation clockwise or counterclockwise*/

    Quat newquat(vc, newaxis.x*vs, newaxis.y*vs, newaxis.z*vs);

    //quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);

    //shift position to geometrical center
    for(unsigned int i=0; i<cluster.size(); i++) {
        //shift position to geometrical center
        conf->pvec[cluster[i]].pos.x -= cluscm.x;
        conf->pvec[cluster[i]].pos.y -= cluscm.y;
        conf->pvec[cluster[i]].pos.z -= cluscm.z;
        //scale things by geo.box not to have them distorted
        conf->pvec[cluster[i]].pos.x *= conf->geo.box.x;
        conf->pvec[cluster[i]].pos.y *= conf->geo.box.y;
        conf->pvec[cluster[i]].pos.z *= conf->geo.box.z;
        //do rotation
        conf->pvec[cluster[i]].pos.rotate(newquat);
        conf->pvec[cluster[i]].dir.rotate(newquat);
        conf->pvec[cluster[i]].patchdir[0].rotate(newquat);
        conf->pvec[cluster[i]].patchdir[1].rotate(newquat);
        conf->pvec[cluster[i]].chdir[0].rotate(newquat);
        conf->pvec[cluster[i]].chdir[1].rotate(newquat);
        conf->pvec[cluster[i]].patchsides[0].rotate(newquat);
        conf->pvec[cluster[i]].patchsides[1].rotate(newquat);
        conf->pvec[cluster[i]].patchsides[2].rotate(newquat);
        conf->pvec[cluster[i]].patchsides[3].rotate(newquat);
        //sclae back
        conf->pvec[cluster[i]].pos.x /= conf->geo.box.x;
        conf->pvec[cluster[i]].pos.y /= conf->geo.box.y;
        conf->pvec[cluster[i]].pos.z /= conf->geo.box.z;
        //shift positions back
        conf->pvec[cluster[i]].pos.x += cluscm.x;
        conf->pvec[cluster[i]].pos.y += cluscm.y;
        conf->pvec[cluster[i]].pos.z += cluscm.z;
    }
}



void MoveCreator::clusterRotate(vector<Particle> &cluster, double max_angle) {
    Vector cluscm;
    double vc,vs;
    Vector newaxis;

    cluscm = clusterCM(cluster);

    // create rotation quaternion
    newaxis.randomUnitSphere(); /*random axes for rotation*/
    vc = cos(max_angle * ran2() );
    if (ran2() <0.5) vs = sqrt(1.0 - vc*vc);
    else vs = -sqrt(1.0 - vc*vc); /*randomly choose orientation of direction of rotation clockwise or counterclockwise*/

    Quat newquat(vc, newaxis.x*vs, newaxis.y*vs, newaxis.z*vs);

    //quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);

    //shift position to geometrical center
    for(unsigned int i=0; i<cluster.size(); i++) {
        //shift position to geometrical center
        cluster[i].pos.x -= cluscm.x;
        cluster[i].pos.y -= cluscm.y;
        cluster[i].pos.z -= cluscm.z;
        //scale things by geo.box not to have them distorted
        cluster[i].pos.x *= conf->geo.box.x;
        cluster[i].pos.y *= conf->geo.box.y;
        cluster[i].pos.z *= conf->geo.box.z;
        //do rotation
        cluster[i].pos.rotate(newquat);
        cluster[i].dir.rotate(newquat);
        cluster[i].patchdir[0].rotate(newquat);
        cluster[i].patchdir[1].rotate(newquat);
        cluster[i].chdir[0].rotate(newquat);
        cluster[i].chdir[1].rotate(newquat);
        cluster[i].patchsides[0].rotate(newquat);
        cluster[i].patchsides[1].rotate(newquat);
        cluster[i].patchsides[2].rotate(newquat);
        cluster[i].patchsides[3].rotate(newquat);
        //sclae back
        cluster[i].pos.x /= conf->geo.box.x;
        cluster[i].pos.y /= conf->geo.box.y;
        cluster[i].pos.z /= conf->geo.box.z;
        //shift positions back
        cluster[i].pos.x += cluscm.x;
        cluster[i].pos.y += cluscm.y;
        cluster[i].pos.z += cluscm.z;
    }
}



double MoveCreator::replicaExchangeMove(long sweep) {
    double edriftchanges=0.0;
#ifdef ENABLE_MPI
    double change; // energy
    double *recwlweights;
    double volume = conf->geo.volume();
    double entrophy = sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper;

    int corr=0;  // correction for sended and receiving processes when clearing all messages
    int sizewl = 0, receiverRank = -1, receivedRank = -1;

    //
    // TAGS for mpi communications
    //
    int tagExchangeMPI = 1001, tagDouble = 2022, tagInt = 3333, tagStat = 654654;

    long localwl,receivedwl;
    bool reject=true;

    MpiExchangeData localmpi, receivedmpi, exch;

    Statistics localStat, recStat;
    localStat = sim->stat;

    if (sim->wl.length[1] > 0) {
        sizewl = sim->wl.length[1] * sim->wl.length[0];
    } else {
        sizewl = sim->wl.length[0];
    }

    recwlweights = (double*) malloc( sizeof(double) * sizewl  );

    MPI_Status status; // int count, int cancelled, int MPI_SOURCE, int MPI_TAG, int MPI_ERROR
    MPI_Request myRequest[sim->mpinprocs];
    MPI_Datatype MPI_exchange;
    MPI_Datatype MPI_vector2;
    MPI_Datatype MPI_stat;

    exch.defDataType(&MPI_exchange, &MPI_vector2);
    localStat.defDataType(&MPI_stat);

    //
    // Init local mpi data
    //
    localmpi.box = conf->geo.box;
    localmpi.energy = calcEnergy->allToAll();
    localmpi.volume = conf->geo.box.x * conf->geo.box.y * conf->geo.box.z;
    localmpi.accepted = 0;
    localmpi.syscm = conf->syscm;
    localmpi.radiusholemax = sim->wl.radiusholemax;
    localmpi.mpiRank = sim->mpirank;
    for(int i=0; i<conf->pvec.molTypeCount; i++) {
        localmpi.partNum[i] = conf->pvec.molCountOfType(i);
    }
    localmpi.pseudoMpiRank = sim->pseudoRank;
    localmpi.temperature = sim->temper;

    for (int wli=0;wli<sim->wl.wlmdim;wli++) {
        localmpi.wl_order[wli] = sim->wl.currorder[wli];
        //fprintf(stdout,"wli %d %ld  %ld\n\n", wli, localmpi.wl_order[wli], sim->wl.currorder[wli] );
    }

    //
    //=== This is an attempt to switch replicas ===
    //

    int oddoreven;

    if ( (sweep % (2*sim->nrepchange)) == 0)
        // exchange odd ones with even ones
        oddoreven=1;
    else
        // exchange even ones with odd ones
        oddoreven=0;
    if (sim->mpinprocs == 2)
        oddoreven=1;


    MPI_Barrier(MPI_COMM_WORLD);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &(rank) );
    if(rank != sim->mpirank) {
        cout << "THIS CAN NEVER HAPPEN!!!" << endl;
        exit(0);
    }

    //
    // SEND, OTHER_PROCESS_EVALUATES, WAIT, RECEIVE, END
    //
    // Processes are sending on rank mod 2 == 0 or mod 2 == 1, which is periodically switched based on sweep value, sending to higher temp
    //
    if (sim->pseudoRank % 2 == oddoreven) {
        if(sim->pseudoRank > 0 )  { // all except for 0

            localmpi.wantedTemp = sim->pTemp[sim->pseudoRank-1];

            // sending NON-blocking, to all processes, only target process will respond based on tag, others just receive but do nothing
            for(int j=0; j<sim->mpinprocs; j++) {
                if(j != sim->mpirank) {
                    MPI_Isend(&localmpi, 1, MPI_exchange, j, sim->pseudoRank-1, MPI_COMM_WORLD, &myRequest[j]);
                }
            }
            corr=1;

            // WAITING for response
            MPI_Recv(&receiverRank, 1, MPI_INT, MPI_ANY_SOURCE, sim->pseudoRank-1+tagInt, MPI_COMM_WORLD, &status); // receive from all processes, ONLY ONE RESPONDS
            MPI_Send(sim->wl.weights, sizewl, MPI_DOUBLE, receiverRank, tagDouble, MPI_COMM_WORLD);
            MPI_Recv(&receivedmpi, 1, MPI_exchange, receiverRank, tagExchangeMPI, MPI_COMM_WORLD, &status);

            if (receivedmpi.accepted == 1) { //decision of accepting or rejecting the exchange was done on other process
                //
                // exchange statistics
                //
                MPI_Recv(&recStat, 1, MPI_stat, receiverRank, tagStat, MPI_COMM_WORLD, &status);
                MPI_Send(&localStat, 1, MPI_stat, receiverRank, tagStat, MPI_COMM_WORLD);

                sim->stat = recStat;

                sim->stat.mpiexch.acc++;

                localmpi.accepted = receivedmpi.accepted;
                conf->geo.box = receivedmpi.box;
                conf->syscm = receivedmpi.syscm;

                edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)conf->pvec.size() * log(receivedmpi.volume / localmpi.volume) / sim->temper;

                sim->temper = receivedmpi.temperature;
                sim->pseudoRank = receivedmpi.pseudoMpiRank;

                volume = conf->geo.volume();
                edriftchanges += (sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper) - entrophy;

            } else {
                sim->stat.mpiexch.rej++;
                if ( sim->wl.wlm[0] > 0 ) {
                    sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]] -= sim->wl.alpha;
                    sim->wl.hist[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]]++;
                }
            }
        }
    } else {

        //
        // WAIT, RECEIVE, EVALUATE, SEND_BACK, END
        //
        if (sim->pseudoRank+1 < sim->mpinprocs ) { // all except MAX
            //there is above process

            MPI_Recv(&receivedmpi, 1, MPI_exchange, MPI_ANY_SOURCE, sim->pseudoRank, MPI_COMM_WORLD, &status);
            corr = 1;

            assert(sim->temper == receivedmpi.wantedTemp && "Wrong message received in MPI");

            receivedRank = receivedmpi.mpiRank;

            MPI_Send(&sim->mpirank, 1, MPI_INT, receivedRank, sim->pseudoRank+tagInt, MPI_COMM_WORLD); // respond to correct process so that it knows who is the correct receiver
            MPI_Recv(recwlweights, sizewl, MPI_DOUBLE, receivedRank, tagDouble, MPI_COMM_WORLD, &status);

            //valuate if accepte or reject the configuration
            //
            // i = here, j = received
            //
            // Canonical: P(1, exp( (E_i - E_j) * (1/kT_i - 1/kT_j) ))
            //
            // Isobaric-Isotermal:
            // acc = exp( (1/T_here - 1/T_received) * (E_here - E_received) + (press /T_here - press_received /(T_received)) * (V_here - V_received) )
            //
            // GrandCanonical
            // P(1, exp( (1/T_here - 1/T_received) * (mu * (N_here - N_received) + (E_here - E_received)) )
            //
            reject = false;
            double temp = (1/sim->temper - 1/(sim->temper + sim->dtemp));

            // Canonical
            change = temp * ( (localmpi.energy - receivedmpi.energy) );
            //printf("acceptance decision: change: %f localE: %f receivedE: %f tempf: %f \n",change,localmpi.energy,receivedmpi.energy,(1/sim->temper - 1/(sim->temper + sim->dtemp)));

            // ISOBARIC-ISOTERMAL
            change += (sim->press/sim->temper - (sim->press + sim->dpress)/(sim->temper + sim->dtemp)) * (localmpi.volume - receivedmpi.volume);

            // GrandCanonical, chempot stored as mu/kT
            for(int i=0; i< conf->pvec.molTypeCount; i++) {
                if(topo.moleculeParam[i].activity != -1)
                    change += temp * topo.moleculeParam[i].chemPot * sim->temper * (localmpi.partNum[i] - receivedmpi.partNum[i]);
            }

            if (sim->wl.wlm[0] > 0) {
                localwl = sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0];
                receivedwl = receivedmpi.wl_order[0] + receivedmpi.wl_order[1]*sim->wl.length[0];
                //fprintf(stdout,"decide wl   %ld %ld %ld energychange: %f \n", receivedmpi.wl_order[0],  receivedmpi.wl_order[1], receivedwl, change );
                //fprintf(stdout,"local weights %ld %f %ld %f \n",localwl,sim->wl.weights[localwl],receivedwl,sim->wl.weights[receivedwl]);
                change += (-sim->wl.weights[localwl] + sim->wl.weights[receivedwl] )/sim->temper + ( -recwlweights[receivedwl] + recwlweights[localwl])/(sim->temper + sim->dtemp) ;
                //fprintf(stdout,"wlchange %f \n\n",change);
            }

            //
            // CRITERION FOR REPLICA EXCHANGE
            //
            if ( (!(reject)) && ( (change > 0) || (ran2() < exp(change))  ) ) {
                // Exchange ACCEPTED send local stuff
                //printf("exchange accepted \n");

                localmpi.accepted = 1;
                conf->geo.box = receivedmpi.box;
                conf->syscm = receivedmpi.syscm;

                edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)conf->pvec.size() * log(receivedmpi.volume / localmpi.volume) / sim->temper;

                // change temperature and pseudorank
                sim->temper = receivedmpi.temperature;
                sim->pseudoRank = receivedmpi.pseudoMpiRank;

                volume = conf->geo.volume();
                edriftchanges += (sim->press * volume - (double)conf->pvec.size() * log(volume) / sim->temper) - entrophy;

                if ( sim->wl.wlm[0] > 0 ) {
                    for (int wli=0;wli<sim->wl.wlmdim;wli++) {
                        sim->wl.neworder[wli] = receivedmpi.wl_order[wli];
                    }
                    sim->wl.accept(sim->wl.wlm[0]);
                }
                MPI_Send(&localmpi, 1, MPI_exchange, receivedRank, tagExchangeMPI, MPI_COMM_WORLD);

                // exchange statistics
                MPI_Send(&localStat, 1, MPI_stat, receivedRank, tagStat, MPI_COMM_WORLD);
                MPI_Recv(&recStat, 1, MPI_stat, receivedRank, tagStat, MPI_COMM_WORLD, &status);

                sim->stat = recStat;

                sim->stat.mpiexch.acc++;

            } else {
                //if exchange rejected send back info
                //printf("exchange rejected\n");
                sim->stat.mpiexch.rej++;
                MPI_Send(&localmpi, 1, MPI_exchange, receivedRank, tagExchangeMPI, MPI_COMM_WORLD);
                if ( sim->wl.wlm[0] > 0 ) {
                    sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]] -= sim->wl.alpha;
                    sim->wl.hist[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]]++;
                }
            }
        }
    }

    //
    // CLEAN MESSAGES SENDED NON-BLOCKING WAY
    //

    // Receive all generated messages, MPI_cancel doesnt guarantee succesful cancel... (even with MPI_WAIT, MPI_request_free)
    MPI_Barrier(MPI_COMM_WORLD); // to ensure we dont read a message that was neccesary for replica exchange

    int size = 0;
    for(int i=0; i<sim->mpinprocs; i++) {
        if (i % 2 == oddoreven) {
            if( i > 0)  { // all except for MAX
                size++;
            }
        }
    }  // size is now how many processors sended messages

    // determine if this process sended messages or received them;
    if(size > 0)
        size = size - corr; // we already received 1 message for odd/even processes and even/odd processes send 1 message
    for(int i=0; i<size; i++) {
        MPI_Recv(&receivedmpi, 1, MPI_exchange, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // i for tag, message from
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Type_free(&MPI_exchange);
    MPI_Type_free(&MPI_vector2);

    free(recwlweights);
#endif

    return edriftchanges;
}

double MoveCreator::muVTMove() {

#ifndef NDEBUG // For tests of energy
    double e = calcEnergy->allToAll();
#endif

    Molecule target;
    double volume = conf->geo.volume();
    double entrophy = log(volume)/sim->temper;
    double energy = 0.0;
    unsigned int molSize=0;
    Vector displace;

    // Determine what type we will be inserting/deleting
    int molType = getRandomMuVTType();
    molSize = topo.moleculeParam[molType].molSize();

    assert(molType == 0 && "delete this assert, only of one atomic type simulation");
    assert(molSize == 1 && "delete this assert, only of one atomic type simulation");

    assert(conf->pvec.molCountOfType(molType) == (int)conf->pvec.size() && "should be true for one atom type simulation");
    assert(insert.empty() && "Insert vector must be empty at the begining of grand canonical move");

    sim->stat.grand[molType].muVtSteps++;

    //////////////////////////////////////////////////////////////
    //                      INSERT MOVE                         //
    //////////////////////////////////////////////////////////////
    if(ran2() > 0.5) {
        if(topo.moleculeParam[molType].isAtomic()) { // RANDOM PARTICLE
            // create particle
            insert.push_back(Particle(conf->geo.randomPos(), Vector::getRandomUnitSphere(), Vector::getRandomUnitSphere()
                                      , molType, topo.moleculeParam[molType].particleTypes[0]));
            insert[0].init(&(topo.ia_params[insert[0].type][insert[0].type]));

            assert(insert[0].testInit( topo.ia_params[insert[0].type][insert[0].type].geotype[0] ) && "GrandCanonical, insertion, Particle initialized incorectly");
            assert(insert.size() == 1);
        } else { // RANDOM CHAIN FROM POOL + RANDOMIZE POSITION AND ROTATION
            displace.randomUnitCube();

            // get configuration
            insert = conf->getRandomPoolConf(molType);

            // randomize position
            for(unsigned int i=0; i<insert.size(); i++)
                insert[i].pos += displace;

            // randomize - rotate chain
            clusterRotate(insert, (double)PIH);

            for(unsigned int i=0; i<insert.size(); i++)
                insert[i].init(&(topo.ia_params[insert[i].type][insert[i].type]));
        }

        assert(!insert.empty());

        // calc energ
        energy = calcEnergy->mol2others(insert);

        // accept with probability -> V/N+1 * e^(ln(a*Nav*1e-27))  -U(new)/kT), NOTE: faunus uses exp(log(V/N+1) * ln(a*Nav*1e-27))  -U(new)/kT)
        if( ( (volume / (conf->pvec.molCountOfType(molType) + 1.0)) *
              (exp( topo.moleculeParam[molType].chemPot - (energy/sim->temper) ) ) ) > ran2() ) {

            if(!topo.moleculeParam[molType].isAtomic())
                energy += calcEnergy->chainInner(insert);

            conf->insertMolecule(insert);

            for(unsigned int i=0; i<insert.size(); i++)
                conf->sysvolume += topo.ia_params[insert[i].type][insert[i].type].volume;

            insert.clear();

            sim->stat.grand[molType].insAcc++;
            sim->stat.grand[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);

            assert((e + energy) > calcEnergy->allToAll()-0.0000001 && (e + energy) < calcEnergy->allToAll()+0.0000001 && "Energy calculated incorectly in grandcanonical insertion");

            return energy - molSize*entrophy;
        } else { // rejected
            insert.clear();
            sim->stat.grand[molType].insRej++;
            sim->stat.grand[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);

            assert(e == calcEnergy->allToAll() && "GrandCanonical, insertion rejected but energy of system changed");

            return 0;
        }
        //////////////////////////////////////////////////////////////
        //                      DELETE MOVE                         //
        //////////////////////////////////////////////////////////////
    } else {
        if(conf->pvec.molCountOfType(molType) == 0) { // check if there are molecules of certain type
            sim->stat.grand[molType].delRej++;
            sim->stat.grand[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);
            return 0;
        }

        target = conf->pvec.getMolecule(ran2() * conf->pvec.molCountOfType(molType), molType, topo.moleculeParam[molType].molSize()); // get random molecule of molType

        energy = calcEnergy->mol2others(target);

        // accept with probability -> N/V * e^(3*ln(wavelenght) - mu/kT + U(del)/kT)

        if( ( ((double)conf->pvec.molCountOfType(molType)/volume) * exp( (energy/sim->temper) - topo.moleculeParam[molType].chemPot) ) > ran2()) {
            for(unsigned int i=0; i<molSize; i++)
                conf->sysvolume -= topo.ia_params[conf->pvec[target[0]+i].type][conf->pvec[target[0]+i].type].volume;

            if(!topo.moleculeParam[molType].isAtomic())
                energy += calcEnergy->chainInner(target);

            conf->removeMolecule(target);

            sim->stat.grand[molType].delAcc++;
            sim->stat.grand[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);

            assert((e - energy) > calcEnergy->allToAll()-0.0000001 && (e - energy) < calcEnergy->allToAll()+0.0000001 && "Energy calculated incorectly in grandcanonical deletion");

            return -energy + molSize*entrophy;
        } else {
            sim->stat.grand[molType].delRej++;
            sim->stat.grand[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);

            assert(e == calcEnergy->allToAll() && "GrandCanonical, deletion rejected but energy of system changed");

            return 0;
        }
    }
    assert(false && "IMPOSIBRU!!!");
    return 0;
}

double MoveCreator::muVTMove2() {

#ifndef NDEBUG // For tests of energy
    double e = calcEnergy->allToAll();
#endif

    int target=-1;
    double volume = conf->geo.volume();
    double entrophy = log(volume)/sim->temper;
    double energy = 0.0;
    unsigned int molSize=1;
    Vector pos, dir, patchdir;
    Particle part;
    int molType = 0;
    int type = 1;

    //////////////////////////////////////////////////////////////
    //                      INSERT MOVE                         //
    //////////////////////////////////////////////////////////////
    if(ran2() > 0.5) {
        //////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////// INITIALIZATION OF PARTICLE ///////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        pos.randomUnitCube();
        dir.randomUnitSphere();
        patchdir.randomUnitSphere();
        part = Particle(pos, dir, patchdir, molType, type);
        part.init(&(topo.ia_params[type][type]));

        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////// TRIAL MOVE //////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        energy = calcEnergy->oneToAll(&part);

        // accept with probability -> V/N+1 * e^(ln(a*Nav*1e-27))  -U(new)/kT), NOTE: faunus uses exp(log(V/N+1) + ln(a*Nav*1e-27))  -U(new)/kT)
        if( ( (volume / (conf->pvec.size() + 1.0)) *
              (exp( topo.moleculeParam[molType].chemPot - (energy/sim->temper) ) ) ) > ran2() ) {


            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////// ACCEPT MOVE /////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            conf->pvec.push_back(part);

            conf->sysvolume += topo.ia_params[type][type].volume;

            assert((e + energy) > calcEnergy->allToAll()-0.0000001 && (e + energy) < calcEnergy->allToAll()+0.0000001 && "Energy calculated incorectly in grandcanonical insertion");

            return energy - molSize*entrophy;
        } else { // rejected
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////// REJECT MOVE /////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////

            assert(e == calcEnergy->allToAll() && "GrandCanonical, insertion rejected but energy of system changed");

            return 0.0;
        }

    //////////////////////////////////////////////////////////////
    //                      DELETE MOVE                         //
    //////////////////////////////////////////////////////////////
    } else {
        if(conf->pvec.size() == 0) { // check if there are molecules of certain type
            return 0.0;
        }

        target = ran2() * conf->pvec.size();

        energy = calcEnergy->oneToAll(target);

        // accept with probability -> N/V * e^(3*ln(wavelenght) - mu/kT + U(del)/kT)
        if( ( ((double)conf->pvec.size()/volume) * exp( (energy/sim->temper) - topo.moleculeParam[molType].chemPot)) > ran2()) {
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////// ACCEPT MOVE /////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////

            conf->sysvolume -= topo.ia_params[type][type].volume;

            conf->pvec.erase(conf->pvec.begin()+target);

            assert((e - energy) > calcEnergy->allToAll()-0.0000001 && (e - energy) < calcEnergy->allToAll()+0.0000001 && "Energy calculated incorectly in grandcanonical deletion");

            return -energy + molSize*entrophy;
        } else {
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////// REJECT MOVE /////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////

            assert(e == calcEnergy->allToAll() && "GrandCanonical, deletion rejected but energy of system changed");

            return 0.0;
        }
    }
    assert(false && "IMPOSIBRU!!!");
    return 0.0;
}


int MoveCreator::getRandomMuVTType() {
    int molType = 0;
    molType = ran2() * topo.gcSpecies;

    for(int i=0; i<conf->pvec.molTypeCount; i++) {
        if(topo.moleculeParam[i].activity != -1.0) {
            if(molType == 0) {
                molType = i;
                break;
            }
            molType--;
        }
    }

    assert(topo.moleculeParam[molType].chemPot != -1.0);
    assert(molType >= 0);
    assert(topo.gcSpecies >= 1 && "GrandCanonical with no defined activity, This should never happen");

    return molType;
}

double MoveCreator::clusterMove() {

    double edriftchanges =0.0;
    long target;

    target = ran2() * (long)conf->pvec.size();// Select random particle from config
    edriftchanges = clusterMoveGeom(target);// Call geometric cluster move
    return edriftchanges;
}

int isInCluster(double *list, int size, double value){
    for(int i=0; i< size; i++){
        if(list[i] == value){
            return 1;
        }
    }
    return 0;
}

double MoveCreator::clusterMoveGeom(long target) {
    /*
     * For reference to this move see:
     * Liu, Jiwen, and Erik Luijten. "Rejection-free geometric cluster algorithm for complex fluids." Physical review letters 92.3 (2004): 035504.
     * DOI: 10.1103/PhysRevLett.92.035504
    */

    double edriftchanges = calcEnergy->allToAll(), cluster[MAXN];
    Vector r_center;

    /*=============================================*/
    /*            Set reflection center            */
    /*=============================================*/
    /*
     * Both ways of selecting reflection center should be equal where in case of Local selection maximal displacement must be set.
     * From test simulations on rather small systems it seems Global relection have faster convergence
    */

    /*____________Global____________*/
    r_center.x=ran2()*conf->geo.box.x;
    r_center.y=ran2()*conf->geo.box.y;
    r_center.z=ran2()*conf->geo.box.z;

    /*____________Local (displacement like)____________*/
//    double max_displacement= 1.5;
//    r_center.randomUnitSphere();// create unit random vector
//    r_center *= ran2() * max_displacement;// set displacement from range [0:max_displacement]
//    r_center += conf->pvec[target].pos;// set center of reflection to be shifted by length of DISPLACEMENT in random direction from target

    Particle reflection;
    int counter= 0, num_particles=1;
    double energy_old, energy_new;

    cluster[0]= target;// first particle allways move so its set to be first in cluster

    /*=============================================*/
    /*            Cluster Creation Loop            */
    /*=============================================*/
    do{
        reflection = conf->pvec[cluster[counter]];// copy old particle into reflected particle
        //Reflect particle cluster[counter] by point reflection by center r_center point
        reflection.pos           = 2.0*r_center - reflection.pos;// reflect center of particle around r_center
        reflection.dir          *=-1.0;// reflect orientation of particle
        reflection.patchdir[0]  *=-1.0;// reflect orientation of patch1
        reflection.patchdir[1]  *=-1.0;// reflect orientation of patch2
        reflection.patchsides[0]*=-1.0;// reflect all sides of patch
        reflection.patchsides[1]*=-1.0;
        reflection.patchsides[2]*=-1.0;
        reflection.patchsides[3]*=-1.0;
        conf->geo.usePBC(&reflection); // bring reflected particle into box (if not particles could start to spread too far and numerical errors acumulate!)

        //Iterate through reflection "Neighbours"
        for (unsigned int i = 0; i < conf->pvec.size(); i++){
            if (!isInCluster(cluster, num_particles, i)){
                energy_old = calcEnergy->p2p(cluster[counter], i);
                energy_new = calcEnergy->p2p(&reflection, i);
                if (ran2() < (1-exp((energy_old-energy_new)/sim->temper))){//ran2() < (1-exp(-1.0*((energy_new-energy_old)/sim->temper))) acceptance criteria vis. Reference
                    cluster[num_particles] = i;
                    num_particles++;
                }
            }
        }
        conf->pvec[cluster[counter]] = reflection;// here old particle is chnged for its reflection
        counter++;
    }while(counter < num_particles);
    return calcEnergy->allToAll()-edriftchanges;
}
