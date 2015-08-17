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
    target = ran2() * (long)conf->pvec.size();

    if ( !( ((sim->wl.wlm[0] == 3) || (sim->wl.wlm[1] == 3) ) && (target == 0) ) && \
    ((ran2() < 0.5) || (topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].geotype[0] >= SP)) ) { /* no rotation for spheres */
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

    /*=== Displacement step ===*/
    edriftchanges =0.0;
    origsyscm.x = 0;
    origsyscm.y = 0;
    origsyscm.z = 0;

    energy = (*calcEnergy)(target, 1, 0);

    orig = conf->pvec[target].pos;

    dr.randomUnitSphere();

    dr.x *= sim->trans[conf->pvec[target].type].mx/conf->geo.box.x;
    dr.y *= sim->trans[conf->pvec[target].type].mx/conf->geo.box.y;
    dr.z *= sim->trans[conf->pvec[target].type].mx/conf->geo.box.z;
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
        sim->trans[conf->pvec[target].type].rej++;
        if ( (sim->wl.wlm[0] == 1) || (sim->wl.wlm[0] == 5) || (sim->wl.wlm[1] == 1) || (sim->wl.wlm[1] == 5) )
            conf->syscm = origsyscm;
        sim->wl.reject(radiusholemax_orig, sim->wl.wlm);

        return 0.0;
    }
#endif

    if (!reject) {  /* wang-landaou ok, try move - calcualte energy */
        enermove =  (*calcEnergy)(target, 1, 0);
    }
    if ( reject || moveTry(energy, enermove, sim->temper) ) {  /* probability acceptance */
        conf->pvec[target].pos = orig;
        sim->trans[conf->pvec[target].type].rej++;
        if ( (sim->wl.wlm[0] == 1) || (sim->wl.wlm[0] == 5) || (sim->wl.wlm[1] == 1) || (sim->wl.wlm[1] == 5) )
            conf->syscm = origsyscm;
        sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
    } else { /* move was accepted */
        sim->trans[conf->pvec[target].type].acc++;
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

    /*=== Rotation step ===*/
    //printf ("rotation %ld npart %ld\n\n",target,npart);
    energy = (*calcEnergy)(target, 1, 0);

    origpart = conf->pvec[target];

//    pscRotate(&conf->pvec[target], sim->rot[conf->pvec[target].type].angle, topo.ia_params[origpart.type][origpart.type].geotype[0]);
    conf->pvec[target].rotateRandom(sim->rot[conf->pvec[target].type].angle, topo.ia_params[origpart.type][origpart.type].geotype[0]);

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
        enermove =  (*calcEnergy)(target, 1, 0);
    }
    if ( reject || moveTry(energy,enermove,sim->temper) ) {  /* probability acceptance */
        conf->pvec[target] = origpart;
        sim->rot[conf->pvec[target].type].rej++;
        sim->wl.reject(sim->wl.radiusholemax, sim->wl.wlm);
    } else { /* move was accepted */
        // DEBUG
        //fprintf(fenergy, "%f\t%f\n", conf->particle[1].pos.x * conf->geo.box.x , enermove);
        sim->rot[conf->pvec[target].type].acc++;
        sim->wl.accept(sim->wl.wlm[0]);
        edriftchanges = enermove - energy + wlener;
        //printf("%f\t%f\n", conf->pvec[0].patchdir[0].z, enermove);
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
    dr.x *= sim->chainm[conf->pvec[chain[0]].molType].mx/conf->geo.box.x;
    dr.y *= sim->chainm[conf->pvec[chain[0]].molType].mx/conf->geo.box.y;
    dr.z *= sim->chainm[conf->pvec[chain[0]].molType].mx/conf->geo.box.z;

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

        sim->chainm[conf->pvec[chain[0]].molType].rej++;
        if ( (sim->wl.wlm[0] == 1) || (sim->wl.wlm[0] == 5) || (sim->wl.wlm[1] == 1) || (sim->wl.wlm[1] == 5) )
            conf->syscm = origsyscm;
        sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
    } else { /* move was accepted */
        sim->chainm[conf->pvec[chain[0]].molType].acc++;
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
    clusterRotate(chain, sim->chainr[conf->pvec[chain[0]].molType].angle);

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

        sim->chainr[conf->pvec[chain[0]].molType].rej++;
        sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
    } else { // move was accepted
        sim->chainr[conf->pvec[chain[0]].molType].acc++;
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
            *side += sim->edge.mx * (ran2() - 0.5);

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
                sim->edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
            } else {  /* move was accepted */
                sim->edge.acc++;
                sim->wl.accept(sim->wl.wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;
        case 1:
            /* Isotropic pressure coupling */
            psch = sim->edge.mx * (ran2() - 0.5);
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
                sim->edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
            } else { /* move was accepted */
                sim->edge.acc++;
                sim->wl.accept(sim->wl.wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;
        case 2:
            /* Isotropic pressure coupling in xy, z constant */
            psch = sim->edge.mx * (ran2() - 0.5);
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
                sim->edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
            } else { /* move was accepted */
                sim->edge.acc++;
                sim->wl.accept(sim->wl.wlm[0]);
                edriftchanges = enermove - energy + wlener;
            }
            break;
        case 3:
            /* Isotropic pressure coupling in xy, z coupled to have fixed volume */
            psch = sim->edge.mx * (ran2() - 0.5);
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
                sim->edge.rej++;
                sim->wl.reject(radiusholemax_orig, sim->wl.wlm);
            } else { /* move was accepted */
                sim->edge.acc++;
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

        //MPI_Datatype type2[11] = {MPI_vector2,MPI_vector2,MPI_vector2,MPI_vector2,MPI_vector2,
        //MPI_LONG, MPI_LONG, MPI_INT,MPI_INT,MPI_DOUBLE, MPI_INT};

        MPI_Datatype type2[10] = {MPI_vector2,MPI_vector2,MPI_vector2,MPI_vector2,MPI_vector2,
                                  MPI_LONG, MPI_INT,MPI_INT,MPI_DOUBLE, MPI_INT};

        //int          blocklen2[11] = {1, 1, 2,4,2,1,1,1,1,1,1,};
        int          blocklen2[10] = {1,1,2,4,2,  1,1,1,1,1,};

        //MPI_Aint     disp2[11];
        MPI_Aint     disp2[10];

        MPI_Address( &part, &dispstart);
        MPI_Address( &(part.pos), &disp2[0]);
        MPI_Address( &(part.dir), &disp2[1]);
        MPI_Address( &(part.patchdir), &disp2[2]);
        MPI_Address( &(part.patchsides), &disp2[3]);
        MPI_Address( &(part.chdir), &disp2[4]);
        MPI_Address( &(part.molType), &disp2[5]);
        //MPI_Address( &(part.chainIndex), &disp2[6]);
        MPI_Address( &(part.type), &disp2[6]);
        MPI_Address( &(part.switchtype), &disp2[7]);
        MPI_Address( &(part.delta_mu), &disp2[8]);
        MPI_Address( &(part.switched), &disp2[9]);

        //for (i=0; i <11; i++) disp2[i] -= dispstart;
        for (i=0; i <10; i++) disp2[i] -= dispstart;

        //MPI_Type_struct( 11, blocklen2, disp2, type2, &MPI_Particle);
        MPI_Type_struct( 10, blocklen2, disp2, type2, &MPI_Particle);

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

        localmpi.box = conf->geo.box;
        localmpi.energy = (*calcEnergy)(0, 0, 0);
        localmpi.volume = conf->geo.box.x * conf->geo.box.y * conf->geo.box.z;
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
                    temppart = (Particle*) malloc(conf->pvec.size()*sizeof(Particle));
                    MPI_Recv(temppart, conf->pvec.size(), MPI_Particle, sim->mpirank-1, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
                    //printf("received data: rank: %d\n", sim->mpirank);
                    //printf("part0  x %f y %f z %f\n",temppart[0].pos.x, temppart[0].pos.y, temppart[0].pos.z);
                    //printf("part1  x %f y %f z %f\n",temppart[1].pos.x, temppart[1].pos.y, temppart[1].pos.z);
                    //printf("part0  molType %ld chainn %ld type %d\n",temppart[0].molType,temppart[0].chainn,temppart[0].type);

                    MPI_Send(&conf->pvec[0], conf->pvec.size(), MPI_Particle, sim->mpirank-1, count, MPI_COMM_WORLD);
                    //printf("send data: rank: %d\n",sim->mpirank);
                    //printf("part0  x %f y %f z %f\n",conf->pvec[0].pos.x,conf->pvec[0].pos.y,conf->pvec[0].pos.z);
                    //printf("part1  x %f y %f z %f\n",conf->particle[1].pos.x,conf->particle[1].pos.y,conf->particle[1].pos.z);
                    //printf("part0  molType %ld chainn %ld type %d\n",conf->pvec[0].molType,conf->pvec[0].chainn,conf->pvec[0].type);

                    localmpi.accepted = receivedmpi.accepted;
                    conf->geo.box = receivedmpi.box;
                    conf->syscm = receivedmpi.syscm;
                    memcpy(&conf->pvec[0],temppart,conf->pvec.size()*sizeof(Particle));
                    edriftchanges = receivedmpi.energy - localmpi.energy;
                    edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)conf->pvec.size() * log(receivedmpi.volume / localmpi.volume) / sim->temper;
                    if ( sim->wl.wlm[0] >0 ) {
                        for (wli=0;wli<sim->wl.wlmdim;wli++) {
                            sim->wl.neworder[wli] = receivedmpi.wl_order[wli];
                        }
                        sim->wl.accept(sim->wl.wlm[0]);
                        //exchange wl data mesh size and radius hole s
                        for (wli=0;wli<sim->wl.wlmdim;wli++) {
                            switch (sim->wl.wlm[wli]) {
                                case 2:
                                    //it is complicated to send because of different sizes
                                     //we would have to send sizes first and realocate corrrect mesh size and then send data
                                    // it is better to recalculate (a bit slower though)
                                    sim->wl.mesh.meshInit(sim->wl.wl_meshsize,conf->pvec.size(),sim->wl.wlmtype,
                                                          conf->geo.box, &conf->pvec);
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
                                    //radiushole_all(topo,conf,sim,wli,&(conf->pvec[0].pos));
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
                    if ( sim->wl.wlm[0] > 0 ) {
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
                if (sim->wl.wlm[0] > 0) {
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
                    conf->geo.box = receivedmpi.box;
                    conf->syscm = receivedmpi.syscm;
                    edriftchanges = receivedmpi.energy - localmpi.energy;
                    edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)conf->pvec.size() * log(receivedmpi.volume / localmpi.volume) / sim->temper;
                    //printf("edrift %f\n",edriftchanges);
                    if ( sim->wl.wlm[0] > 0 ) {
                        for (wli=0;wli<sim->wl.wlmdim;wli++) {
                            sim->wl.neworder[wli] = receivedmpi.wl_order[wli];
                        }
                        sim->wl.accept(sim->wl.wlm[0]);
                    }
                    MPI_Send(&localmpi, 1, MPI_exchange, sim->mpirank+1, count, MPI_COMM_WORLD);
                    //printf("send data: rank: %d energy: %f volume: %f pressure: %f \n",sim->mpirank,localmpi.energy,localmpi.volume,localmpi.pressure);
                    //send and receive configuration
                    MPI_Send(&conf->pvec[0], conf->pvec.size(), MPI_Particle, sim->mpirank+1, count, MPI_COMM_WORLD);
                    //printf("send data: rank: %d\n",sim->mpirank);
                    //printf("part0  x %f y %f z %f\n",conf->pvec[0].pos.x,conf->pvec[0].pos.y,conf->pvec[0].pos.z);
                    //printf("part1  x %f y %f z %f\n",conf->particle[1].pos.x,conf->particle[1].pos.y,conf->particle[1].pos.z);
                    //printf("part0  molType %ld chainn %ld type %d\n",conf->pvec[0].molType,conf->pvec[0].chainn,conf->pvec[0].type);

                    MPI_Recv(&conf->pvec[0], conf->pvec.size(), MPI_Particle, sim->mpirank+1, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
                    //printf("recieved data: rank: %d\n",sim->mpirank);
                    //printf("part0  x %f y %f z %f\n",conf->pvec[0].pos.x,conf->pvec[0].pos.y,conf->pvec[0].pos.z);
                    //printf("part1  x %f y %f z %f\n",conf->particle[1].pos.x,conf->particle[1].pos.y,conf->particle[1].pos.z);
                    //printf("part0  molType %ld chainn %ld type %d\n",conf->pvec[0].molType,conf->pvec[0].chainn,conf->pvec[0].type);

                    if ( sim->wl.wlm[0] > 0 ) {
                        //exchange wl data mesh size and radius hole s
                        for (wli=0;wli<sim->wl.wlmdim;wli++) {
                            switch (sim->wl.wlm[wli]) {
                                case 2:
                                    //it is complicated to send because of different sizes
                                    //  we would have to send sizes first and realocate corrrect mesh size and then send data
                                    //  it is better to recalculate (a bit slower though)
                                    sim->wl.mesh.meshInit(sim->wl.wl_meshsize, conf->pvec.size(), sim->wl.wlmtype,
                                                          conf->geo.box, &conf->pvec);
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
                                    //radiushole_all(topo,conf,sim,wli,&(conf->pvec[0].pos));
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
                    if ( sim->wl.wlm[0] > 0 ) {
                        sim->wl.weights[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]] -= sim->wl.alpha;
                        sim->wl.hist[sim->wl.currorder[0]+sim->wl.currorder[1]*sim->wl.length[0]]++;
                    }
                }
            }
        }
        //if ( (localmpi.accepted) && (sim->pairlist_update) ) gen pair list

        MPI_Type_free(&MPI_exchange);
        MPI_Type_free(&MPI_Particle);
        MPI_Type_free(&MPI_vector2);
        free(recwlweights);
#endif
    return edriftchanges;
}

double MoveCreator::muVTMove() {

#ifndef NDEBUG
    double e = calcEnergy->allToAll();
    vector<Particle> pTemp;
    for(unsigned int i=0; i<conf->pvec.size(); i++) {
        pTemp.push_back(conf->pvec[i]);
    }
#endif

    Molecule target;
    double volume = conf->geo.volume();
    double entrophy = log(volume)/sim->temper;
    double energy = 0.0;
    double factor = 1.0;
    unsigned int molSize=0;
    Vector displace;

    // Determine what type we will be inserting/deleting
    int molType = getRandomMuVTType();
    molSize = topo.moleculeParam[molType].molSize();

    assert(conf->pvec.molCountOfType(molType) == (int)conf->pvec.size() && "should be true for one atom type simulation");

    topo.moleculeParam[molType].muVtSteps++;
    if(ran2() > 0.5) { //  insert move
        if(topo.moleculeParam[molType].isAtomic()) {
            // create particle
            insert.push_back(Particle(conf->geo.randomPos(), Vector::getRandomUnitSphere(), Vector::getRandomUnitSphere()
                                      , molType, topo.moleculeParam[molType].particleTypes[0]));
            insert[0].init(&topo.ia_params[insert[0].type][insert[0].type]);
            assert(insert[0].testInit());

            energy = calcEnergy->oneToAll(&insert[0], NULL, NULL);

            // accept with probability -> V/N+1 * e^(ln(a*Nav*1e-27))  -U(new)/kT)
            if( ( (volume / (conf->pvec.molCountOfType(molType) + 1.0)) *
                  exp( topo.moleculeParam[molType].chemPot - (energy/sim->temper) ) ) > ran2()) {

                conf->insertMolecule(insert);

                insert.clear();
                conf->sysvolume += topo.ia_params[insert[0].type][insert[0].type].volume;
                topo.moleculeParam[molType].insAcc++;
                topo.moleculeParam[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);

                assert((e + energy) > calcEnergy->allToAll()-0.0000001 && (e + energy) < calcEnergy->allToAll()+0.0000001 );
#ifndef NDEBUG
                for(unsigned int i=0; i<pTemp.size(); i++) {
                    assert(pTemp[i] == conf->pvec[i]);
                }
                assert(conf->pvec.size() == pTemp.size()+1);
#endif

                return energy - entrophy;
            } else { // rejected
                insert.clear();
                topo.moleculeParam[molType].insRej++;
                topo.moleculeParam[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);

                assert(e == calcEnergy->allToAll());
#ifndef NDEBUG
                for(unsigned int i=0; i<conf->pvec.size(); i++) {
                    assert(pTemp[i] == conf->pvec[i]);
                }
#endif

                return 0;
            }

        } else { // this is chain insert
            displace.randomUnitCube();

            // get configuration
            insert = conf->getRandomPoolConf(molType);

            // randomize position
            displace -= insert[0].pos;
            for(unsigned int i=0; i<insert.size(); i++)
                insert[i].pos += displace;

            // randomize - rotate chain
            clusterRotate(insert, (double)PIH);

            // check overlap
            for(unsigned int i=0; i<insert.size(); i++) {
                if(conf->overlapAll(&insert[i], topo.ia_params)) {
                    insert.clear();
                    topo.moleculeParam[molType].insRej++;
                    topo.moleculeParam[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);
                    return 0; // overlap detected, move rejected
                }
            }

            //generate conlists
            vector<ConList> con;
            for(unsigned int j=0; j<insert.size(); j++) {
                con.push_back(ConList() );

                if (j > 0) //if this is not first particle fill tail bond
                    con[j].conlist[0] = &insert[j-1];

                if ( j+1 < insert.size() ) //if there is a next particle fill it to head bond
                    con[j].conlist[1] = &insert[j+1];

                if (j > 1) //if this is not second or first particle fill second tail bond
                    con[j].conlist[2] = &insert[j-2];

                if ( j+2 < insert.size() ) //if there is a second next particle fill it second neighbour
                    con[j].conlist[3] = &insert[j+2];
            }

            // calc energ
            energy += calcEnergy->mol2others(insert);

            factor *= volume / (conf->pvec.molCountOfType(molType) + 1.0);

            // accept with probability -> V/N+1 * e^(ln(a*Nav*1e-24))  -U(new)/kT)
            if( ( factor * exp( topo.moleculeParam[molType].chemPot - (energy/sim->temper) ) ) > ran2()) {
                // add internal energy(with external)
                energy += calcEnergy->chainInner(insert, con);

                conf->insertMolecule(insert);

                for(unsigned int i=0; i<insert.size(); i++)
                    conf->sysvolume += topo.ia_params[insert[i].type][insert[i].type].volume;

                insert.clear();         

                topo.moleculeParam[molType].insAcc++;
                topo.moleculeParam[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);

                return energy - molSize*entrophy;
            } else {
                topo.moleculeParam[molType].insRej++;
                topo.moleculeParam[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);

                assert(e == calcEnergy->allToAll());

                insert.clear();
                return 0;
            }
        }

    } else { // delete move
        if(conf->pvec.molCountOfType(molType) == 0) // check if there are molecules of certain type
            return 0;

        target = conf->pvec.getMolecule(ran2() * conf->pvec.molCountOfType(molType), molType, topo.moleculeParam[molType].molSize()); // get random molecule of molType

#ifndef NDEBUG
        Particle temp = conf->pvec[target[0]];
#endif

        energy = calcEnergy->mol2others(target);

        // accept with probability -> N/V * e^(3*ln(wavelenght) - mu/kT + U(del)/kT)
        if( ( (conf->pvec.molCountOfType(molType)/volume) * exp( (energy/sim->temper) - topo.moleculeParam[molType].chemPot)) > ran2()) {
            for(unsigned int i=0; i<molSize; i++)
                conf->sysvolume -= topo.ia_params[conf->pvec[target[0]+i].type][conf->pvec[target[0]+i].type].volume;

            energy += calcEnergy->chainInner(target);

            conf->removeMolecule(target);
            topo.moleculeParam[molType].delAcc++;
            topo.moleculeParam[molType].muVtAverageParticles += conf->pvec.molCountOfType(molType);

            assert((e - energy) > calcEnergy->allToAll()-0.0000001 && (e - energy) < calcEnergy->allToAll()+0.0000001 );

#ifndef NDEBUG
            for(unsigned int i=0; i< conf->pvec.size(); i++)
                assert(temp != conf->pvec[i] && "Wrong particle deleted");
#endif

            return -energy + molSize*entrophy;
        } else {
            topo.moleculeParam[molType].delRej++;
            topo.moleculeParam[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);

            assert(e == calcEnergy->allToAll());

            return 0;
        }
    }
    return 0;
}


int MoveCreator::getRandomMuVTType() {
    int molType = 0;
    for(int i=0; i<conf->pvec.molTypeCount; i++) {
        if(topo.moleculeParam[i].activity != -1.0) molType++;
    }
    molType = (long) (ran2() * ((double)molType));
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

    return molType;
}

double MoveCreator::clusterMove() {

    double edriftchanges =0.0;
    long target;

    target = ran2() * (long)conf->pvec.size();// Select random particle from config
    edriftchanges = clusterMoveGeom(target);// Call geometric cluster move
    return edriftchanges;
}

int MoveCreator::moveTry(double energyold, double energynew, double temperature) {
    /*DEBUG   printf ("   Move trial:    %13.8f %13.8f %13.8f %13.8f\n",
      energynew, energyold, temperature, ran2(&seed));*/
    if (energynew <= energyold ) {
        return 0;
    } else {
        if (exp(-1.0*(energynew-energyold)/temperature) > ran2()) {
            return 0;
        } else {
            return 1;
        }
    }
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
