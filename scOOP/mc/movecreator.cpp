/** @file movecreator.cpp*/

#include "movecreator.h"
#include <iomanip>
#include <algorithm>

#ifdef ENABLE_MPI
# include <mpi.h>
extern MPI_Datatype MPI_vector, MPI_Particle, MPI_exchange;
#endif

double MoveCreator::particleMove() {
    double edriftchanges =0.0;
    long target;

    /*=== This is a particle move step ===*/
    target = ran2() * (long)conf->pvec.size();

    if ( /*!( ((wl.wlm[0] == 3) || (wl.wlm[1] == 3) ) && (target == 0) ) && chceme at se hybe*/ \
         ((ran2() < 0.5) || (topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].geotype[0] >= SP)) ) { // no rotation for spheres
        edriftchanges = partDisplace(target);
    } else {
        //=== Rotation step ===//
        // BTW: partAcialRotate pro uhel 180.0 a pouziti Vector::getRandomUnitConeUniform by se mel chovat stejne jako normalni partRotate ....
        if(sim->coneAngle == 0.0){
            edriftchanges = partRotate(target);
        } else {
            edriftchanges = partAxialRotate(target);
        }

    }
    //=== End particle move step ===
    return edriftchanges;
}

double MoveCreator::clusterMove() {

    double edriftchanges =0.0;
    long target;

    target = ran2() * (long)conf->pvec.size();// Select random particle from config
    edriftchanges = clusterMoveGeom(target);// Call geometric cluster move
    return edriftchanges;
}

int MoveCreator::isInCluster(double *list, int size, double value){
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

    if(conf->pvec.empty())
        return 0.0;

    double edriftchanges = calcEnergy->allToAll(), cluster[MAXN];
    Vector r_center;

    /*=============================================*/
    /*            Set reflection center            */
    /*=============================================*/
    /*There exists two ways how to select reflection center, Global and Local.
     * Global ---- select random point in simulation box as reflection center
     * Local  ---- select vector from selected particle of random length in interval (0; max_len>
     *
     * Both ways of selecting reflection center should be equal where in case of Local selection maximal displacement (max_len) must be set.
     *
     * From test simulations on rather small systems it seems Global relection have faster convergence
    */

    /*____________Global____________*/
    r_center = conf->geo.randomPos();

    /*____________Local (displacement like)____________*/
    //    double max_displacement= 1.5;
    //    r_center.randomUnitSphere();// create unit random vector
    //    r_center *= ran2() * max_displacement;// set displacement from range [0:max_displacement]
    //    r_center += conf->pvec[target].pos;// set center of reflection to be shifted by length of DISPLACEMENT in random direction from target

    Particle reflection;
    Molecule selected_chain;
    std::vector<Molecule> chainsToFix; // thse are all cahins that might be disturbed by application of PBC
    int counter= 0, num_particles=0;
    double energy_old, energy_new;

    /*=============================================*/
    /*     Addition of particles into cluster      */
    /*=============================================*/
    /*
     * Here we chose to add whole chain into cluster if one particle of chain is includedtarget in cluster
     * This make move slightly faster but
     * !!!!SIMULATION OF CHAINS WITHOUT SINGLE PARTICLE MOVES CANT CONVERGE!!!!
     *
     * TODO:    change it in way that intra chain energy is used to determine if other particles from chain should be added
     *          into cluster
    */

    double molecule_size;
    //topo.moleculeParam[conf->pvec[target].molType].particleTypes.size() == number of particles in chain ... special case is single particle of length 1
    molecule_size = topo.moleculeParam[conf->pvec[target].molType].particleTypes.size();
    if ( molecule_size == 1 ){
        cluster[num_particles] = target;
        num_particles++;
    }else{
        selected_chain = conf->pvec.getMolOfPart(target);
        chainsToFix.push_back(selected_chain);
        for(unsigned int i=0; i < selected_chain.size(); i++){
            cluster[num_particles] = selected_chain[i];
            num_particles++;
        }
    }

    /*=============================================*/
    /*            Cluster Creation Loop            */
    /*=============================================*/
    do{
        reflection = conf->pvec[cluster[counter]];// copy old particle into reflected particle
        //Reflect particle cluster[counter] by point reflection by center r_center point
        reflection.pos           = 2.0*r_center - reflection.pos;// reflect center of particle around r_center
        reflection.dir          *=1.0;// reflect orientation of particle
        reflection.patchdir[0]  *=-1.0;// reflect orientation of patch1
        reflection.patchdir[1]  *=-1.0;// reflect orientation of patch2
        reflection.patchsides[0]*=-1.0;// reflect all sides of patch
        reflection.patchsides[1]*=-1.0;
        reflection.patchsides[2]*=-1.0;
        reflection.patchsides[3]*=-1.0;
        reflection.chdir[0]     *=-1.0;
        reflection.chdir[1]     *=-1.0;
        conf->geo.usePBC(&reflection);

        // bring reflected particle into box (if not particles could start to spread too far and numerical errors acumulate!)

        //Iterate through reflection "Neighbours"
        for (unsigned int i = 0; i < conf->pvec.size(); i++){
            if (!isInCluster(cluster, num_particles, i)){
                energy_old = calcEnergy->p2p(cluster[counter], i);
                energy_new = calcEnergy->p2p(&reflection, i);
                if (ran2() < (1-exp((energy_old-energy_new)/sim->temper))){//ran2() < (1-exp(-1.0*((energy_new-energy_old)/sim->temper))) acceptance criteria vis. Reference
                    //Addition of chain into cluster
                    //-----------------------------------------------------
                    molecule_size = topo.moleculeParam[conf->pvec[i].molType].particleTypes.size();
                    if(molecule_size == 1){
                        cluster[num_particles] = i;
                        num_particles++;
                    }else{
                        selected_chain = conf->pvec.getMolOfPart(i);
                        chainsToFix.push_back(selected_chain);
                        for(unsigned int t=0; t < selected_chain.size(); t++){
                            cluster[num_particles] = selected_chain[t];
                            num_particles++;
                        }
                    }
                    //-----------------------------------------------------
                }
            }
        }
        conf->pvec[cluster[counter]] = reflection;
        counter++;
    }while(counter < num_particles);
    for ( std::vector<Molecule>::iterator it = chainsToFix.begin(); it != chainsToFix.end(); ++it ){
        conf->makeMoleculeWhole(&(*it));
    }
    return calcEnergy->allToAll()-edriftchanges;
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

double MoveCreator::switchTypeMove() {
    double edriftchanges=0.0, energy,enermove=0.0, switchE=0.0,wlener=0.0;
    int reject=0;
    long target;
    double radiusholemax_orig=0;
    int switchType, sequence_num=0, delta_mu;

    //=== This is an attempt to switch a type ===
    target = ran2() * conf->pvec.switchPartCount();
    target = conf->pvec.getSwitchPart(target, sequence_num); // stores sequence number
    delta_mu = topo.moleculeParam[conf->pvec[target].molType ].deltaMu[sequence_num];
    if(conf->pvec[target].switched == 0)
        switchType = topo.moleculeParam[conf->pvec[target].molType ].switchTypes[sequence_num];
    else
        switchType = topo.moleculeParam[conf->pvec[target].molType ].particleTypes[sequence_num];

    DEBUG_SIM("Switching the particle type");
    DEBUG_SIM("PARTICLE: %ld", target);

    energy = calcEnergy->oneToAll(target);

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

    if(wl.wlm[0] > 0)
        wlener = wl.runSwitch(reject, radiusholemax_orig);

    if (!reject) {
        switchE = delta_mu * pmone;
        // DEBUG
        //double dmu = enermove;
        //pvec[target].switched += pmone;
        enermove = calcEnergy->oneToAllTrial(target);

        //printf("energy: %f \t %f\t%f\n",pvec[target].delta_mu, dmu, enermove);
    }

    // If not accepted: switch back
    if ( reject || moveTry(energy+wlener,enermove+switchE,sim->temper) ) {  // probability acceptance
        DEBUG_SIM("Did NOT switch it\n");
        conf->pvec[target].type = tmp_type;
        conf->pvec[target].switched -= pmone;
        conf->pvec[target].init(&(topo.ia_params[conf->pvec[target].type][conf->pvec[target].type]));
        wl.reject(radiusholemax_orig, wl.wlm);
        sim->stat.switchMv[conf->pvec[target].type].rej++;
    } else { // move was accepted
        wl.accept(wl.wlm[0]);
        sim->stat.switchMv[conf->pvec[target].type].acc++;
        edriftchanges = enermove - energy;

        calcEnergy->update(target);
    }

    return edriftchanges;
}

double MoveCreator::chainMove() {

#ifndef NDEBUG
    double testE = calcEnergy->allToAllTrial();
#endif

    double edriftchanges =0.0;
    long target;

    if(conf->pvec.getChainCount() == 0) // no chains to displace - muVTmove deleted all
        return 0.0;

    //=== This is a chain move step ===
    target = ran2() * conf->pvec.getChainCount();

    if (ran2() < 0.5) { //=== Displacement step of cluster/chain ===
        edriftchanges = chainDisplace(target);
        assert(fabs( testE - calcEnergy->allToAllTrial() + edriftchanges) < 0.000001 || !(cerr << "Drift chainMove translate: " << testE - calcEnergy->allToAllTrial() + edriftchanges << endl) );
    } else {            //=== Rotation step of cluster/chain ===
        edriftchanges = chainRotate(target);
        assert(fabs( testE - calcEnergy->allToAllTrial() + edriftchanges) < 0.000001 || !(cerr << "Drift chainMove rotate: " << testE - calcEnergy->allToAllTrial() + edriftchanges << endl) );
    }
    return edriftchanges;
}

double MoveCreator::pressureMove() {
    double edriftchanges,energy,enermove=0.0,wlener;
    int reject=0;
    double old_side;   // geo.box length before attempted change
    double *side;      // geo.box dimension to try changing
    double psch;       // Size of a geo.box change during pressure
    double pvol;       // Size of a volume during pressure
    double pvoln;      // Size of a new volume during pressure
    double rsave;      // Saved random number
    double area;
    double radiusholemax_orig=0;

    //=== This is a volume change step ===
    edriftchanges=0.0;
    wlener = 0.0;

    energy = calcEnergy->allToAll();

    // Choose an edge
    switch (sim->ptype) {
    case 0:
        // Anisotropic pressure coupling
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
        if (wl.wlm[0] > 0) {  // get new neworder for wang-landau
            wlener = wl.runPress(reject, radiusholemax_orig);
        }
        if (!reject) { // wang-landaou ok, try move - calculate energy
            // enermove = (p*(V_new - V_old)                  - N * k * T * ln(V_new / V_old)
            enermove = sim->press * area * (*side - old_side) - (double)conf->pvec.size() * sim->temper * log(*side/old_side);

            enermove += calcEnergy->allToAllTrial();
        }
        //
        // probability(old -> new) = min[ 1, exp{-(E_new - E_old)/(kT) + (p*(V_new - V_old) - NkT* ln(V_new / V_old))/(kT) ]}
        //
        if ( reject || *side <= 0.0 || ( moveTry(energy+wlener,enermove,sim->temper) ) ) { // probability acceptance
            *side = old_side;
            sim->stat.edge.rej++;
            wl.reject(radiusholemax_orig, wl.wlm);
        } else {  // move was accepted
            sim->stat.edge.acc++;
            calcEnergy->update();
            wl.accept(wl.wlm[0]);
            edriftchanges = enermove - energy;
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
        if (wl.wlm[0] > 0) {  /* get new neworder for wang-landau */
            wlener = wl.runPress(reject, radiusholemax_orig);
        }
        if (!reject) { /* wang-landaou ok, try move - calcualte energy */
            // enermove = (p*(V_new - V_old)       - N * k * T * ln(V_new / V_old)
            enermove = sim->press * (pvoln - pvol) - (double)conf->pvec.size() * sim->temper * log(pvoln/pvol);

            enermove += calcEnergy->allToAllTrial();
        }
        if ( reject || moveTry(energy+wlener,enermove,sim->temper) )  { /* probability acceptance */
            conf->geo.box.x -= psch;
            conf->geo.box.y -= psch;
            conf->geo.box.z -= psch;
            sim->stat.edge.rej++;
            wl.reject(radiusholemax_orig, wl.wlm);
        } else { // move was accepted
            sim->stat.edge.acc++;
            wl.accept(wl.wlm[0]);
            calcEnergy->update();
            edriftchanges = enermove - energy;
        }
        break;
    case 2:
        // Isotropic pressure coupling in xy, z constant
        psch = sim->stat.edge.mx * (ran2() - 0.5);
        pvol = conf->geo.box.x * conf->geo.box.y;
        conf->geo.box.x += psch;
        conf->geo.box.y += psch;
        pvoln = conf->geo.box.x * conf->geo.box.y;

        reject = 0;
        if (wl.wlm[0] > 0) {  // get new neworder for wang-landau
            wlener = wl.runPress(reject, radiusholemax_orig, true);
        }
        if (!reject) { // wang-landaou ok, try move - calculate energy
            // enermove = (p*(V_new - V_old)       - N * k * T * ln(V_new / V_old)
            enermove = sim->press * conf->geo.box.z * (pvoln - pvol) - (double)conf->pvec.size() * sim->temper * log(pvoln/pvol);

            enermove += calcEnergy->allToAllTrial();
        }
        if ( reject || moveTry(energy+wlener,enermove,sim->temper) )  { // probability acceptance
            conf->geo.box.x -= psch;
            conf->geo.box.y -= psch;
            sim->stat.edge.rej++;
            wl.reject(radiusholemax_orig, wl.wlm);
        } else { // move was accepted
            sim->stat.edge.acc++;
            wl.accept(wl.wlm[0]);
            calcEnergy->update();
            edriftchanges = enermove - energy;
        }
        break;
    case 3:
        // Isotropic pressure coupling in xy, z coupled to have fixed volume
        psch = sim->stat.edge.mx * (ran2() - 0.5);
        pvol = conf->geo.box.x * conf->geo.box.y * conf->geo.box.z;
        conf->geo.box.x += psch;
        conf->geo.box.y += psch;
        conf->geo.box.z = pvol / conf->geo.box.x / conf->geo.box.y;

        reject = 0;
        if (wl.wlm[0] > 0) {  // get new neworder for wang-landau
            wlener = wl.runPress(reject, radiusholemax_orig);
        }
        if (!reject) { // wang-landaou ok, try move - calculate energy
            enermove += calcEnergy->allToAllTrial();
        }
        if ( reject || moveTry(energy+wlener,enermove,sim->temper) )  { // probability acceptance
            conf->geo.box.x -= psch;
            conf->geo.box.y -= psch;
            conf->geo.box.z = pvol / conf->geo.box.x / conf->geo.box.y;
            sim->stat.edge.rej++;
            wl.reject(radiusholemax_orig, wl.wlm);
        } else { // move was accepted
            sim->stat.edge.acc++;
            wl.accept(wl.wlm[0]);
            calcEnergy->update();
            edriftchanges = enermove - energy;
        }
        break;

    default:
        fprintf (stderr, "ERROR: unknown type of pressure coupling %d",sim->ptype);
        exit(1);
    }

    //=== End volume change step ===
    return edriftchanges;
}

double MoveCreator::replicaExchangeMove(long sweep) {
    double edriftchanges=0.0;
#ifdef ENABLE_MPI
    double change; // energy
    double *recwlweights;
    double volume = conf->geo.volume();
    double entrophy = sim->press * volume - (double)conf->pvec.size() * log(volume) * sim->temper;

    int sizewl = 0, receiverRank = -1, receivedRank = -1;
    int rank;

    //
    // TAGS for mpi communications
    //
    int tagExchangeMPI = 1001, tagDouble = 2022, tagInt = 3333, tagStat = 54321;

    long localwl,receivedwl;
    bool reject=true;

    MpiExchangeData localmpi, receivedmpi;
    MpiExchangeData* locMpi;
    MpiExchangeData* recMpi;

    Statistics localStat, recStat;
    localStat = sim->stat;

    if (wl.length[1] > 0) {
        sizewl = wl.length[1] * wl.length[0];
    } else {
        sizewl = wl.length[0];
    }

    recwlweights = (double*) malloc( sizeof(double) * sizewl  );

    MPI_Status status;
    MPI_Datatype MPI_exchange;
    MPI_Datatype MPI_vector2;
    MPI_Datatype MPI_stat;

    localmpi.defDataType(&MPI_exchange, &MPI_vector2);
    localStat.defDataType(&MPI_stat);

    //
    // Init local mpi data
    //
    localmpi.energy = calcEnergy->allToAll();
    localmpi.press = sim->press;
    localmpi.volume = conf->geo.box.x * conf->geo.box.y * conf->geo.box.z;
    localmpi.accepted = 0;
    localmpi.radiusholemax = wl.radiusholemax;
    localmpi.mpiRank = sim->mpirank;
    for(int i=0; i<conf->pvec.molTypeCount; i++) {
        localmpi.partNum[i] = conf->pvec.molCountOfType(i);
    }
    localmpi.pseudoMpiRank = sim->pseudoRank;
    localmpi.temperature = sim->temper;

    for (int wli=0;wli<wl.wlmdim;wli++) {
        localmpi.wl_order[wli] = wl.currorder[wli];
    }

    //
    //=== This is an attempt to switch replicas ===
    //
    int oddoreven;

    if ( (sweep % (2*sim->nrepchange)) == 0) // exchange odd ones with even ones
        oddoreven=1;
    else
        oddoreven=0; // exchange even ones with odd ones
    if (sim->mpinprocs == 2)
        oddoreven=1;

    if (sim->pseudoRank % 2 == oddoreven) {
        if(sim->pseudoRank > 0 ) {  // all except for 0
            localmpi.wantedTemp = sim->pTemp[sim->pseudoRank-1];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank );
    assert(rank != sim->mpirank && "THIS CAN NEVER HAPPEN!!!");

    locMpi = new MpiExchangeData[sim->mpinprocs];
    recMpi = new MpiExchangeData[sim->mpinprocs];

    for ( int i=0 ; i < sim->mpinprocs ; ++i ) {
        locMpi[i] = localmpi;
    }

    MPI_Alltoall(locMpi, 1, MPI_exchange, recMpi, 1, MPI_exchange, MPI_COMM_WORLD);

    for ( int i=0 ; i < sim->mpinprocs ; ++i ) {
        if(recMpi[i].wantedTemp == sim->temper)
            receivedmpi = recMpi[i];
    }

    delete locMpi;
    delete recMpi;


    //
    // SEND, OTHER_PROCESS_EVALUATES, WAIT, RECEIVE, END
    //
    // Processes are sending on rank mod 2 == 0 or mod 2 == 1, which is periodically switched based on sweep value, sending to higher temp
    //
    if (sim->pseudoRank % 2 == oddoreven) {
        if(sim->pseudoRank > 0 )  { // all except for 0

            MPI_Recv(&receiverRank, 1, MPI_INT, MPI_ANY_SOURCE, sim->pseudoRank-1+tagInt, MPI_COMM_WORLD, &status); // receive from all processes, ONLY ONE RESPONDS
            MPI_Send(wl.shared_weights, sizewl, MPI_DOUBLE, receiverRank, tagDouble, MPI_COMM_WORLD);
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

                edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)conf->pvec.size() * log(receivedmpi.volume / localmpi.volume) * sim->temper;

                sim->temper = receivedmpi.temperature;
                sim->press = receivedmpi.press;
                sim->pseudoRank = receivedmpi.pseudoMpiRank;

                volume = conf->geo.volume();
                edriftchanges += (sim->press * volume - (double)conf->pvec.size() * log(volume) * sim->temper) - entrophy;

            } else {
                sim->stat.mpiexch.rej++;
                if ( wl.wlm[0] > 0 ) {
                    wl.shared_weights[wl.currorder[0]+wl.currorder[1]*wl.length[0]] -= wl.alpha;
                    wl.shared_hist[wl.currorder[0]+wl.currorder[1]*wl.length[0]]++;
                }
            }
        }
    } else {

        //
        // WAIT, RECEIVE, EVALUATE, SEND_BACK, END
        //
        if (sim->pseudoRank+1 < sim->mpinprocs ) { // all except MAX
            //there is above process

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

            // ISOBARIC-ISOTERMAL
            change += (sim->press/sim->temper - (sim->press + sim->dpress)/(sim->temper + sim->dtemp)) * (localmpi.volume - receivedmpi.volume);

            // GrandCanonical, chempot stored as mu/kT
            for(int i=0; i< conf->pvec.molTypeCount; i++) {
                if(topo.moleculeParam[i].activity != -1)
                    change += temp * topo.moleculeParam[i].chemPot * sim->temper * (localmpi.partNum[i] - receivedmpi.partNum[i]);
            }

            if (wl.wlm[0] > 0) {
                localwl = wl.currorder[0]+wl.currorder[1]*wl.length[0];
                receivedwl = receivedmpi.wl_order[0] + receivedmpi.wl_order[1]*wl.length[0];
                change += (-wl.shared_weights[localwl] + wl.shared_weights[receivedwl] )/sim->temper + ( -recwlweights[receivedwl] + recwlweights[localwl])/(sim->temper + sim->dtemp) ;
            }

            //
            // CRITERION FOR REPLICA EXCHANGE
            //
            if ( (!(reject)) && ( (change > 0) || (ran2() < exp(change))  ) ) { // Exchange ACCEPTED send local stuff
                localmpi.accepted = 1;

                edriftchanges += sim->press * (receivedmpi.volume - localmpi.volume) - (double)conf->pvec.size() * log(receivedmpi.volume / localmpi.volume) * sim->temper;

                // change temperature and pseudorank
                sim->temper = receivedmpi.temperature;
                sim->press = receivedmpi.press;
                sim->pseudoRank = receivedmpi.pseudoMpiRank;

                volume = conf->geo.volume();
                edriftchanges += (sim->press * volume - (double)conf->pvec.size() * log(volume) * sim->temper) - entrophy;

                if ( wl.wlm[0] > 0 ) {
                    for (int wli=0;wli<wl.wlmdim;wli++) {
                        wl.neworder[wli] = receivedmpi.wl_order[wli];
                    }
                    wl.accept(wl.wlm[0]);
                }
                MPI_Send(&localmpi, 1, MPI_exchange, receivedRank, tagExchangeMPI, MPI_COMM_WORLD);
                // exchange statistics
                MPI_Send(&localStat, 1, MPI_stat, receivedRank, tagStat, MPI_COMM_WORLD);
                MPI_Recv(&recStat, 1, MPI_stat, receivedRank, tagStat, MPI_COMM_WORLD, &status);

                sim->stat = recStat;

                sim->stat.mpiexch.acc++;

            } else {
                //if exchange rejected send back info
                sim->stat.mpiexch.rej++;
                localmpi.accepted = 0;
                MPI_Send(&localmpi, 1, MPI_exchange, receivedRank, tagExchangeMPI, MPI_COMM_WORLD);
                if ( wl.wlm[0] > 0 ) {
                    wl.shared_weights[wl.currorder[0]+wl.currorder[1]*wl.length[0]] -= wl.alpha;
                    wl.shared_hist[wl.currorder[0]+wl.currorder[1]*wl.length[0]]++;
                }
            }         
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); // to ensure we dont read a message that was neccesary for replica exchange

    MPI_Type_free(&MPI_exchange);
    MPI_Type_free(&MPI_vector2);
    MPI_Type_free(&MPI_stat);

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
    double entrophy = log(volume)*sim->temper;
    double energy = 0.0;
    unsigned int molSize=0;
    Vector displace;

    // Determine what type we will be inserting/deleting
    int molType = getRandomMuVTType();
    molSize = topo.moleculeParam[molType].molSize();

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

            assert(insert[0].testInit( topo.ia_params[insert[0].type][insert[0].type].geotype[0], 99999999 ) && "GrandCanonical, insertion, Particle initialized incorectly");
            assert(insert.size() == 1);
        } else { // RANDOM CHAIN FROM POOL + RANDOMIZE POSITION AND ROTATION
            displace.randomUnitCube();

            // get configuration
            insert = conf->getRandomPoolConf(molType);

            // Scale to box
            for(auto &item : insert) {
                item.pos.x /= conf->geo.box.x;
                item.pos.y /= conf->geo.box.y;
                item.pos.z /= conf->geo.box.z;

                conf->geo.usePBC(&item);
            }

            // randomize position
            for(unsigned int i=0; i<insert.size(); i++)
                insert[i].pos += displace;

            // randomize - rotate chain
            clusterRotate(insert, (double)PIH);

            for(unsigned int i=0; i<insert.size(); i++)
                insert[i].init(&(topo.ia_params[insert[i].type][insert[i].type]));
        }

        assert(!insert.empty() && "Insert vector must contain a particle");

        // calc energ
        energy = calcEnergy->mol2others(insert);

        // accept with probability -> V/N+1 * e^(ln(a*Nav*1e-27))  -U(new)/kT), NOTE: faunus uses exp(log(V/N+1) * ln(a*Nav*1e-27))  -U(new)/kT)
        if(( (volume / (conf->pvec.molCountOfType(molType) + 1.0)) *
             (exp( topo.moleculeParam[molType].chemPot - (energy/sim->temper) ) ) ) > ran2() ) {

            if(!topo.moleculeParam[molType].isAtomic())
                energy += calcEnergy->chainInner(insert);

            conf->insertMolecule(insert);

            for(unsigned int i=0; i<insert.size(); i++)
                conf->sysvolume += topo.ia_params[insert[i].type][insert[i].type].volume;

            insert.clear();

            sim->stat.grand[molType].insAcc++;
            sim->stat.grand[molType].muVtAverageParticles +=  conf->pvec.molCountOfType(molType);
            calcEnergy->update(EMResize()); // just a dummy datatype for resize of energy matrix

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
            calcEnergy->update(EMResize()); // just a dummy datatype for resize of energy matrix

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

int MoveCreator::getRandomMuVTType() {
    int molType = 0;
    molType = ran2() * topo.gcSpeciesCount;

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
    assert(topo.gcSpeciesCount >= 1 && "GrandCanonical with no defined activity, This should never happen");

    return molType;
}

double MoveCreator::partDisplace(long target) {
    double edriftchanges = 0.0,energy,enermove,wlener = 0.0;
    chorig[0] = conf->pvec[target];
    Vector dr, origsyscm = conf->syscm;
    int reject=0;
    double radiusholemax_orig=0;

    energy = calcEnergy->oneToAll(target);

    dr.randomUnitSphere();

    dr.x *= sim->stat.trans[conf->pvec[target].type].mx/conf->geo.box.x;
    dr.y *= sim->stat.trans[conf->pvec[target].type].mx/conf->geo.box.y;
    dr.z *= sim->stat.trans[conf->pvec[target].type].mx/conf->geo.box.z;

    conf->pvec[target].pos.x += dr.x;
    conf->pvec[target].pos.y += dr.y;
    conf->pvec[target].pos.z += dr.z;

    if(wl.wlm[0] > 0) {
        Molecule tar;
        tar.resize(1);
        tar[0] = target;
        dr.scale(topo.ia_params[conf->pvec[target].type][conf->pvec[target].type].volume);
        wlener = wl.run(reject, chorig, dr, radiusholemax_orig, tar);
    }

    if (!reject) {  // wang-landaou ok, try move - calcualte energy
        enermove = calcEnergy->oneToAllTrial(target);
    }
    if (reject || moveTry(energy+wlener, enermove, sim->temper)) {  // probability acceptance
        conf->pvec[target].pos = chorig[0].pos;
        sim->stat.trans[conf->pvec[target].type].rej++;
        if ( (wl.wlm[0] == 1) || (wl.wlm[0] == 5) || (wl.wlm[1] == 1) || (wl.wlm[1] == 5) )
            conf->syscm = origsyscm;
        wl.reject(radiusholemax_orig, wl.wlm);

    } else { // move was accepted
        sim->stat.trans[conf->pvec[target].type].acc++;
        wl.accept(wl.wlm[0]);

        edriftchanges = enermove - energy;

        calcEnergy->update(target);
    }

    return edriftchanges;
}

double MoveCreator::partRotate(long target) {
    double edriftchanges = 0.0, energy, enermove, wlener = 0.0;
    Particle origpart = conf->pvec[target];
    int reject = 0;

    energy = calcEnergy->oneToAll(target);

    conf->pvec[target].rotateRandom(sim->stat.rot[conf->pvec[target].type].angle, topo.ia_params[origpart.type][origpart.type].geotype[0]);

    //should be normalised and ortogonal but we do for safety
    assert(isSame(conf->pvec[target].dir.size(), 1.0));
    conf->pvec[target].patchdir[0].ortogonalise(conf->pvec[target].dir);

    if(wl.wlm[0] > 0)
        wlener = wl.runRot(reject, target);

    if (!reject) {  // wang-landaou ok, try move - calcualte energy
        enermove = calcEnergy->oneToAllTrial(target);
    }
    if ( reject || moveTry(energy+wlener,enermove,sim->temper) ) {  // probability acceptance
        conf->pvec[target] = origpart;
        sim->stat.rot[conf->pvec[target].type].rej++;
        wl.reject(wl.radiusholemax, wl.wlm);
    } else { // move was accepted
        sim->stat.rot[conf->pvec[target].type].acc++;
        wl.accept(wl.wlm[0]);
        edriftchanges = enermove - energy;

        calcEnergy->update(target);
    }

    return edriftchanges;
}

double MoveCreator::partAxialRotate(long target){
    double edriftchanges = 0.0;
    double energyold;
    double energynew = 0.0;

    Vector   rotaxis;
    Particle origpart        =   conf->pvec[target];

    energyold = calcEnergy->oneToAll(target);

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
    energynew = calcEnergy->oneToAllTrial(target);

    if (moveTry(energyold, energynew, sim->temper)){
        // move was rejected
        sim->stat.rot[conf->pvec[target].type].rej++;
        conf->pvec[target] = origpart; // return to old configuration
    } else {
        // move was accepted
        sim->stat.rot[conf->pvec[target].type].acc++;
        edriftchanges = energynew - energyold;

        calcEnergy->update(target);
    }

    return edriftchanges;
}

double MoveCreator::chainDisplace(long target) {
    Molecule chain = conf->pvec.getChain(target);
#ifndef NDEBUG
    double innerE = calcEnergy->chainInner(chain);
    double allE = calcEnergy->allToAllTrial();
#endif
    assert(chain.size() > 1);
    double edriftchanges=0.0,energy=0.0,enermove=0.0,wlener=0.0;
    Vector dr, origsyscm = conf->syscm;
    int reject=0;
    Vector cluscm(0.0, 0.0, 0.0);
    double radiusholemax_orig=0.0;

    //=== Displacement step of cluster/chain ===
    for(unsigned int i=0; i<chain.size(); i++) // store old configuration
        chorig[i].pos = conf->pvec[chain[i]].pos;

    energy += calcEnergy->mol2others(chain);

    dr.randomUnitSphere();
    dr.x *= sim->stat.chainm[conf->pvec[chain[0]].molType].mx/conf->geo.box.x;
    dr.y *= sim->stat.chainm[conf->pvec[chain[0]].molType].mx/conf->geo.box.y;
    dr.z *= sim->stat.chainm[conf->pvec[chain[0]].molType].mx/conf->geo.box.z;

    if ( ((wl.wlm[0] == 3)||(wl.wlm[1] == 3)) && (target == 0) ) {
        dr.z = 0;
        dr.y = 0;
        dr.x = 0;
    }
    for(unsigned int j=0; j<chain.size(); j++) { // move chaine to new position
        if ( (wl.wlm[0] == 1) || (wl.wlm[0] == 5) || (wl.wlm[1] == 1) || (wl.wlm[1] == 5) ) { /* calculate move of center of mass  */
            cluscm.x += dr.x*topo.ia_params[conf->pvec[chain[j]].type][conf->pvec[chain[j]].type].volume;
            cluscm.y += dr.y*topo.ia_params[conf->pvec[chain[j]].type][conf->pvec[chain[j]].type].volume;
            cluscm.z += dr.z*topo.ia_params[conf->pvec[chain[j]].type][conf->pvec[chain[j]].type].volume;
        }
        conf->pvec[chain[j]].pos.x += dr.x;
        conf->pvec[chain[j]].pos.y += dr.y;
        conf->pvec[chain[j]].pos.z += dr.z;
    }

    if(wl.wlm[0] > 0)
        wlener = wl.run(reject, chorig, cluscm, radiusholemax_orig, chain);

    if (!reject) { // wang-landaou ok, try move - calcualte energy
        enermove += calcEnergy->mol2othersTrial(chain);
    }
    if ( reject || moveTry(energy+wlener, enermove, sim->temper) ) {  // probability acceptance
        for(unsigned int j=0; j<chain.size(); j++)
            conf->pvec[chain[j]].pos = chorig[j].pos;

        sim->stat.chainm[conf->pvec[chain[0]].molType].rej++;
        if ( (wl.wlm[0] == 1) || (wl.wlm[0] == 5) || (wl.wlm[1] == 1) || (wl.wlm[1] == 5) )
            conf->syscm = origsyscm;
        wl.reject(radiusholemax_orig, wl.wlm);

        assert( fabs( allE - calcEnergy->allToAllTrial()) < 0.00001
                || !(cerr << allE << " " << calcEnergy->allToAllTrial() << " allToAll energy changed on reject, chainDisplace" << endl));
    } else { // move was accepted
        sim->stat.chainm[conf->pvec[chain[0]].molType].acc++;
        wl.accept(wl.wlm[0]);

        calcEnergy->update(chain);

        edriftchanges = enermove - energy;
        assert( fabs( allE - calcEnergy->allToAllTrial() + edriftchanges ) < 0.00001
                || !(cerr << allE << " " << calcEnergy->allToAllTrial() << " " << edriftchanges << " allToAll energy not match, chainDisplace" << endl));
    }

    assert( fabs( innerE - calcEnergy->chainInner(chain) ) < 0.00001 || !(cerr << innerE - calcEnergy->chainInner(chain) << " Inner energy changed, chainDisplace" << endl));

    return edriftchanges;
}

double MoveCreator::chainRotate(long target) {
    Molecule chain = conf->pvec.getChain(target);
#ifndef NDEBUG
    double testE = calcEnergy->chainInner(chain);
#endif
    double edriftchanges=0.0, energy=0.0, enermove=0.0, wlener=0.0;
    int reject=0;
    Particle chorig[MAXCHL];
    double radiusholemax_orig=0;

    //=== Rotation step of cluster/chain ===
    for(unsigned int j=0; j<chain.size(); j++) { // store old configuration calculate energy
        chorig[j] = conf->pvec[chain[j]];
    }

    energy += calcEnergy->mol2others(chain);

    //do actual rotations around geometrical center
    clusterRotate(chain, sim->stat.chainr[conf->pvec[chain[0]].molType].angle);

    if(wl.wlm[0] > 0)
        wlener = wl.runChainRot(reject, chorig, radiusholemax_orig, chain);

    if (!reject) { // wang-landaou ok, try move - calcualte energy
        enermove += calcEnergy->mol2othersTrial(chain);
    }
    if ( reject || moveTry(energy+wlener, enermove, sim->temper) ) { // probability acceptance
        for(unsigned int j=0; j<chain.size(); j++)
            conf->pvec[chain[j]] = chorig[j];

        sim->stat.chainr[conf->pvec[chain[0]].molType].rej++;
        wl.reject(radiusholemax_orig, wl.wlm);
    } else { // move was accepted
        sim->stat.chainr[conf->pvec[chain[0]].molType].acc++;
        wl.accept(wl.wlm[0]);
        edriftchanges = enermove - energy;

        calcEnergy->update(chain);
    }

    assert( fabs( testE - calcEnergy->chainInner(chain) ) < 0.00001 || !(cerr << testE - calcEnergy->chainInner(chain) << " Inner energy change" << endl));

    return edriftchanges;
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

Vector MoveCreator::clusterCM(vector<Particle> &cluster) {
    double chainVolume=0.0;
    Vector cluscm(0.0, 0.0, 0.0);

    for(unsigned int i=0; i<cluster.size(); i++) {
        cluscm.x += cluster[i].pos.x * topo.ia_params[cluster[i].type][cluster[i].type].volume;
        cluscm.y += cluster[i].pos.y * topo.ia_params[cluster[i].type][cluster[i].type].volume;
        cluscm.z += cluster[i].pos.z * topo.ia_params[cluster[i].type][cluster[i].type].volume;

        chainVolume += topo.ia_params[cluster[i].type][cluster[i].type].volume;
    }

    cluscm.x /= chainVolume;
    cluscm.y /= chainVolume;
    cluscm.z /= chainVolume;

    return cluscm;
}

Vector MoveCreator::clusterCM(vector<int> &cluster) {
    double chainVolume=0.0;
    Vector cluscm(0.0, 0.0, 0.0);

    for(unsigned int i=0; i<cluster.size(); i++) {
        cluscm.x += conf->pvec[cluster[i]].pos.x * topo.ia_params[conf->pvec[cluster[i]].type][conf->pvec[cluster[i]].type].volume;
        cluscm.y += conf->pvec[cluster[i]].pos.y * topo.ia_params[conf->pvec[cluster[i]].type][conf->pvec[cluster[i]].type].volume;
        cluscm.z += conf->pvec[cluster[i]].pos.z * topo.ia_params[conf->pvec[cluster[i]].type][conf->pvec[cluster[i]].type].volume;

        chainVolume += topo.ia_params[conf->pvec[cluster[i]].type][conf->pvec[cluster[i]].type].volume;
    }

    cluscm.x /= chainVolume;
    cluscm.y /= chainVolume;
    cluscm.z /= chainVolume;

    return cluscm;
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
