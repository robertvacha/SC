#include "inicializer.h"

#include <iostream>

#include "simlib.h"
#include "mygetline.h"
#include "randomGenerator.h"

#ifdef ENABLE_MPI
# include <mpi.h>
extern MPI_Datatype MPI_vector, MPI_Particle, MPI_exchange;
#endif

extern Topo topo;

void Inicializer::readOptions() {

    if(sim->mpirank == 0)
        cout << "Reading options..." << endl;

    int num_options = -1;
    double transmx, rotmx, chainmmx, chainrmx, angle, chain_angle;
    long int seed;

    char *id, *value, *tokLine, *line;
    FILE *infile;

    // for new options add before the last line
    Option options[] = {
        {"write_cluster",       Long,   false, &sim->write_cluster},
        {"adjust",              Long,   false, &sim->adjust},
        {"movie",               Long,   false, &sim->movie},
        {"nequil",              Long,   false, &sim->nequil},
        {"nsweeps",             Long,   false, &sim->nsweeps},
        {"nrepchange",          Long,   false, &sim->nrepchange},
        {"nGrandCanon",         Long,   false, &sim->nGrandCanon},
        {"nClustMove",          Long,   false, &sim->nClustMove},
        {"paramfrq",            Long,   false, &sim->paramfrq},
        {"report",              Long,   false, &sim->report},
        {"seed",                Long,   false, &seed},
        {"pairlist_update",     Int,    false, &sim->pairlist_update},
        {"ptype",               Int,    false, &sim->ptype},
        {"wlm",                 Int2,   false, &sim->wl.wlm},
        {"wlmtype",             Int,    false, &sim->wl.wlmtype},
        {"press",               Double, false, &sim->press},
        {"paralpress",          Double, false, &sim->paralpress},
        {"edge_mx",             Double, false, &sim->stat.edge.mx},
        {"shave",               Double, false, &sim->shave},
        {"chainprob",           Double, false, &sim->chainprob},
        {"switchprob",          Double, false, &sim->switchprob},
        {"temper",              Double, false, &sim->temper},
        {"paraltemper",         Double, false, &sim->paraltemper},
        {"transmx",             Double, false, &transmx},
        {"rotmx",               Double, false, &rotmx},
        {"coneAngle",           Double, true, &sim->coneAngle}, // default value given in constructor of sim
        {"chainmmx",            Double, false, &chainmmx},
        {"chainrmx",            Double, false, &chainrmx},
        {"last",                Int,    false, NULL}
    };
    while(options[++num_options].var != NULL)
        ;

    //--- 1. Read in values ---
    size_t line_size = (STRLEN + 1) * sizeof(char);
    line = (char *) malloc(line_size);

    infile = fopen(files->optionsfile, "r");
    if (infile == NULL) {
        fprintf (stderr, "\nERROR: Could not open options file.\n\n");
        exit (1);
    }

    while(myGetLine(&line, &line_size, infile) != -1){

        // strip comments
        strip_comment(line);
        trim(line);
        if(strlen(line) == 0) continue;

        // tokenize
        tokLine = line;
        id = strtok(tokLine, "=");
        if(id == NULL){
            fprintf(stderr, "error parsing Configuration line (%s)", line);
            free(line);
            exit(1);
        }
        trim(id);
        tokLine = NULL;
        value = strtok(tokLine, "=");
        trim(value);
        if(value == NULL){
            fprintf(stderr, "error parsing Configuration line (%s)", line);
            free(line);
            exit(1);
        }
        //printf("id: %s; value: %s\n", id, value);
        int i = 0;
        for(i=0; i < num_options; i++){
            if(strcmp(id, options[i].id) == 0){
        if(options[i].type == Int2){
            readii2(value,*((int (*)[2]) options[i].var));
                    options[i].set = true;
                    break;
                }
                if(options[i].type == Int){
                    *((int *) options[i].var) = readi2(value);
                    options[i].set = true;
                    break;
                }
                else if(options[i].type == Long){
                    *((long *) options[i].var) = readl2(value);
                    options[i].set = true;
                    break;
                }
                else if(options[i].type == Double){
                    *((double *) options[i].var) = readd2(value);
                    options[i].set = true;
                    break;
                }
                else {
                    fprintf(stderr, "Could not determine type of %s!\n", id);
                    free(line);
                    exit(1);
                }
            }
        }
        if(i == num_options){
            fprintf(stderr, "Unknown identifier %s!\nWill procede.\n", id);
        }
    }
    fclose (infile);
    free(line);

    // Check, wheter all options have been readin
    for(int i = 0; i < num_options; i++){
        if(!options[i].set){
            fprintf(stderr, "option '%s' is not set!\n", options[i].id);
            exit(1);
        }
    }

    //--- 2. Summarize results on standard output ---
    // Density of close-packed spherocylinders
    //   rho_cp = 2.0/(sqrt(2.0) + *length * sqrt(3.0));

    if(sim->mpirank == 0 && SILENT == 1) {
        printf (" Pressure coupling type:                             %d\n", sim->ptype);
        printf (" Pressure:                                           %.8f\n", sim->press);
        printf (" Replica exchange pressure:                          %.8f\n", sim->paralpress);
        printf (" Average volume change attempts per sweep:           %.8f\n", sim->shave);
        printf (" Equilibration sweeps:                               %ld\n", sim->nequil);
        printf (" Sweeps between step size adjustments:               %ld\n", sim->adjust);
        printf (" Production sweeps:                                  %ld\n", sim->nsweeps);
        printf (" Sweeps between statistics samples:                  %ld\n", sim->paramfrq);
        printf (" Sweeps between statistics reports:                  %ld\n", sim->report);
        printf (" Average chain move attempts per sweep:              %.8f\n", sim->chainprob);
        printf (" Initial maximum displacement:                       %.8f\n", transmx);
        printf (" Inititial maximum angular change (degrees):         %.8f\n", rotmx);
        printf (" Inititial maximum angular cone angle (degrees):     %.8f\n", sim->coneAngle);
        printf (" Inititial maximum geo.box edge change:              %.8f\n", sim->stat.edge.mx);
        printf (" Initial maximum chain displacement:                 %.8f\n", chainmmx);
        printf (" Inititial maximum chain angular change (degrees):   %.8f\n", chainrmx);
        printf (" Temperature in kT/e:                                %.8f\n", sim->temper);
        printf (" Parallel tempering temperature in kT/e:             %.8f\n", sim->paraltemper);
        printf (" Sweeps between replica exchange:                    %ld\n", sim->nrepchange);
        printf (" Sweeps between Grand-Canonical move:                %ld\n", sim->nGrandCanon);
        printf (" Sweeps between Cluster moves:                       %ld\n", sim->nClustMove);
        printf (" Wang-Landau method:                                 %d %d\n", sim->wl.wlm[0],sim->wl.wlm[1]);
        printf (" Calculate the Wang-Landau method for atom type:     %d\n", sim->wl.wlmtype);
        printf (" Average type switch attempts per sweep:             %.8f\n", sim->switchprob);
        printf (" Number of Sweeps per pairlist update:               %d\n", sim->pairlist_update);
        printf (" Random number seed:                                 %ld\n", seed);
        printf (" Number of sweeps per writing out cluster info:      %ld\n", sim->write_cluster);

        if (sim->movie > 0) {
            printf (" Sweeps between movie frames:                      %ld\n", sim->movie);
        } else {
            printf (" No movie\n");
        }
        printf ("\n");

        if(sim->pairlist_update){
            printf(" A pairlist will be generated every %d steps. This is a greedy"
                   " algorithm; make sure you don't have big chains etc.!\n",
                   sim->pairlist_update);
        }
    }

    //--- 3. Validity checks ---
    if (rotmx < 0.0 || rotmx > 180) {
        fprintf (stderr, "ERROR: Maximum orientation change must be in range 0 to 180.\n\n");
        exit (1);
    }
    if (chainrmx < 0.0 || chainrmx > 180) {
        fprintf (stderr, "ERROR: Maximum orientation change for chains must be in range 0 to 180.\n\n");
        exit (1);
    }
    if ( (sim->ptype <0) || (sim->ptype>3) ) {
        fprintf (stderr, "ERROR: Unknown pressure coupling %d. Program only knows: 0 - anisotropic coupling, \
                1 - isotropic coupling, 2 - isotropic in xy z=const, 3 - isotropic xy V=const.\n\n",sim->ptype);
        exit (1);
    }
    if ( (sim->wl.wlm[0] <0) || (sim->wl.wlm[0] > 7) || (sim->wl.wlm[1] <0) || (sim->wl.wlm[1] > 7)  ) {
        fprintf (stderr, "ERROR: Unknown Wang-Landau method %d %d. Program only knows: 0 - none, \
                1 - z-direction od 1st particle, 2 - pore in membrane, 3 - zorientation of 0th particle,\
                4 - distance of fist two particles, 5 - pore around z-axis above CM,\
                6 - pore around z-axis above 0th particle, 7 - number of particles in contact \n\n",sim->wl.wlm[0],sim->wl.wlm[1]);
        exit (1);
    }
    if ( (sim->wl.wlm[0] == 0) && (sim->wl.wlm[1] > 0)  ) {
        fprintf (stderr, "ERROR: Wang-Landau method has to be set for first order parameter and then for second order parameter\n\n");
        exit (1);
    }
    if ( (sim->wl.wlm[0] == 2) || (sim->wl.wlm[0] == 5) || (sim->wl.wlm[0] == 6)  ) {
        if(sim->wl.wlmtype < 1){
            fprintf (stderr, "ERROR: Atom type for the Wang-Landau Method (%d) was false defined.\n\n",sim->wl.wlmtype);
            exit (1);
        }
        if ( (sim->wl.wlm[1] == 2) || (sim->wl.wlm[1] == 5) || (sim->wl.wlm[1] == 6) ) {
            fprintf (stderr, "ERROR: Simulaneous use of two pore order parameters has not been implemented yet.\n\n");
            exit (1);
        }
    }

    // we store maximum rotation as half angle - useful for quaterions
    angle = rotmx / 180.0 * PIH *0.5;
    rotmx = cos((rotmx)/180.0*PIH);
    chain_angle = chainrmx / 180.0 * PIH;
    chainrmx = cos((chainrmx)/180.0*PIH);
    sim->stat.edge.mx *= 2.0;   // The full range is -maxl to +maxl, i.e. spanning 2*maxl
    transmx *= 2.0;   // The full range is -maxr to +maxr, i.e. spanning 2*maxr
    chainmmx *= 2.0;   // The full range is -maxr to +maxr, i.e. spanning 2*maxr
    sim->coneAngle *= DEGTORAD; // Now transfer angle in degrees into radians

    for (int i=0;i<MAXT;i++) {
        sim->stat.trans[i].mx = transmx;
        sim->stat.rot[i].mx = rotmx;
        sim->stat.rot[i].angle = angle;
    }
    for (int i=0;i<MAXMT;i++) {
        sim->stat.chainm[i].mx = chainmmx;
        sim->stat.chainr[i].mx = chainrmx;
        sim->stat.chainr[i].angle = chain_angle;
    }

    //parallel tempering
#ifdef ENABLE_MPI
    if ( (sim->temper != sim->paraltemper) && (sim->mpinprocs <2) ) {
        printf("ERROR: Paralllel tempering at single core does not work.\n\n");
        exit(1);
    }
    sim->dtemp = (sim->paraltemper - sim->temper )/(sim->mpinprocs-1);
    for(int i=0; i<sim->mpinprocs; i++) {
        sim->pTemp.push_back(sim->temper + (sim->dtemp * i));
    }
    sim->temper += sim->dtemp * sim->mpirank;
    if ( (sim->press != sim->paralpress) && (sim->mpinprocs <2) ) {
        printf("ERROR: Pressure replica exchange at single core does not work.\n\n");
        exit(1);
    }
    sim->dpress = (sim->paralpress - sim->press )/(sim->mpinprocs-1);
    sim->press += sim->dpress * sim->mpirank;
    seed += sim->mpirank;
    sim->stat.mpiexch.mx = sim->dtemp;
    sim->stat.mpiexch.angle = sim->dpress;
#endif

    ran2.setSeed(seed);
}

void Inicializer::initTop() {

    bool exclusions[MAXT][MAXT] = {false};

    readTopoFile(exclusions); // EXCLUDE LOADED CORRECTLY 7.8. 2015

    fprintf (stdout, "\nTopology succesfully read. Generating pair interactions...\n");

    //fill ia_params combinations and topology parameters
    topo.genParamPairs(exclusions);
    topo.genTopoParams();

    setParticlesParams();
    initClusterList();
    initSwitchList();
    initGroupLists();

    conf->sysvolume = 0;
    for (unsigned int i=0; i<conf->pvec.size(); i++)
        conf->sysvolume += topo.ia_params[conf->pvec[i].type][conf->pvec[i].type].volume;

    if(sim->nGrandCanon != 0) {
        bool existGrand = false;
        int i=0;
        while(topo.moleculeParam[i].name != NULL) {
            if(topo.moleculeParam[i].activity != -1.0 )
                existGrand = true;
            i++;
        }
        if(!existGrand) {
            cout << "In options nGrandCanon != 0, but no activity set for any species in top.init" << endl;
            exit(1);
        }
    }

    DEBUG_INIT("Finished with reading the topology");

#ifdef ENABLE_MPI  // Parallel tempering check
    // probability to switch replicas = exp ( -0.5 * dT*dT * N / (1 + dT) )
    printf("Probability to switch replicas is roughly: %f\n",exp(-0.5 * conf->pvec.size() * sim->dtemp * sim->dtemp / (1.0 + sim->dtemp)) );
#endif

    topDealoc();
}

void Inicializer::initClusterList() {
    sim->clusterlist = (long int*) malloc(sizeof(long) * MAXN);
    if(sim->clusterlist == NULL){
        fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sim->clusterlist!");
        exit(1);
    }
    sim->clustersenergy = (double*) malloc(sizeof(double) * MAXN);
    if(sim->clustersenergy== NULL){
        fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sim->clustersenergy!");
        exit(1);
    }
    sim->clusters = NULL;
}

void Inicializer::initSwitchList() {
    // count switch types for all molecular types
    int count;
    bool switchPartExist = false;
    for(int i=0; i<MAXMT; i++) {
        if(topo.moleculeParam[i].particleTypes.empty())
            break;
        count =0;
        for(unsigned int j=0; j<topo.moleculeParam[i].switchTypes.size(); j++) {
            if(topo.moleculeParam[i].switchTypes[j] != -1) {
                count++;
                switchPartExist = true;
            }
        }
        topo.moleculeParam[i].switchCount = count;
    }

    if (!switchPartExist && sim->switchprob > 0){
        fprintf(stderr, "TOPOLOGY WARNING: No switchable particles found, but probability for a switch is not zero!\n");
        sim->switchprob = 0;
        fprintf(stderr, "TOPOLOGY WARNING: We changed Switch Probability to zero in this run!\n");
    }

    //  Mark particles as not switched
    for(unsigned int i = 0; i < conf->pvec.size(); i++){
        conf->pvec[i].switched = 0;
    }
}

bool Inicializer::initConfig(FILE** infile, std::vector<Particle > &pvec) {

    int err,fields,tmp_type;
    long j,current,first;
    char * line, line2[STRLEN];
    size_t line_size = (STRLEN + 1) * sizeof(char);
    line = (char *) malloc(line_size);
    Particle chorig[MAXCHL];

    double maxlength = 0.0;
    for(int i = 0; i < MAXT; i++){
        if(maxlength < topo.ia_params[i][i].len[0])
            maxlength = topo.ia_params[i][i].len[0];
    }

    if(myGetLine(&line, &line_size, *infile) == -1){
        fprintf (stderr, "ERROR: Could not read box size1.\n\n");
        return false;
    }
    strip_comment(line);
    trim(line);
#ifdef WEDGE
    double angle, innerR, outerR;
    Vector box;
    if (sscanf(line, "%le %le %le %le", &outerR, &innerR, &box.z, &angle) != 4) {
        if(myGetLine(&line, &line_size, infile) == -1){
            fprintf (stderr, "ERROR: Could not read box size.\n\n");
            return false;
        }
        aftercommand(line2,line,BOXSEP);
        strip_comment(line2);
        trim(line2);
        if (sscanf(line2, "%le %le %le %le", &box.z, &angle, &outerR, &innerR) != 4) {
            fprintf (stderr, "ERROR: Could not read box size.\n\n");
            return false;
        }
    }

    conf->geo = Wedge(box.z, angle, outerR, innerR); //(double box.z, double angle, double outerR, double innerR)
#else
    Vector box;

    if (sscanf(line, "%le %le %le", &(box.x), &(box.y), &(box.z) ) != 3) {
        if(myGetLine(&line, &line_size, *infile) == -1){
            fprintf (stderr, "ERROR: Could not read box size2.\n\n");
            return false;
        }
        aftercommand(line2,line,BOXSEP);
        strip_comment(line2);
        trim(line2);
        if (sscanf(line2, "%le %le %le", &(box.x), &(box.y), &(box.z) ) != 3) {
            fprintf (stderr, "ERROR: Could not read box size3.\n\n");
            return false;
        }
    }

    conf->geo = Cuboid(box);
#endif
    if (conf->geo.box.x < maxlength * 2.0 + 2.0) {
        printf ("WARNING: x (%f) geo.box length is less than two spherocylinders long (%f).\n\n", conf->geo.box.x, maxlength * 2.0 + 2.0);
        exit(1);
    }
    if (conf->geo.box.y < maxlength * 2.0 + 2.0) {
        printf ("WARNING: y (%f) geo.box length is less than two spherocylinders long (%f).\n\n", conf->geo.box.y, maxlength * 2.0 + 2.0);
        exit(1);
    }
    if (conf->geo.box.z < maxlength * 2.0 + 2.0) {
        printf ("WARNING: z (%f) geo.box length is less than two spherocylinders long (%f).\n\n", conf->geo.box.z, maxlength * 2.0 + 2.0);
        exit(1);
    }

    DEBUG_INIT("Position of the particle");
    for(unsigned int i=0; i < pvec.size(); i++) {
        if(myGetLine(&line, &line_size, *infile) == -1){
            break;
        }
        strip_comment(line);
        trim(line);

        fields = sscanf(line, "%le %le %le %le %le %le %le %le %le %d",
                        &pvec[i].pos.x, &pvec[i].pos.y, &pvec[i].pos.z,
                        &pvec[i].dir.x, &pvec[i].dir.y, &pvec[i].dir.z,
                        &pvec[i].patchdir[0].x, &pvec[i].patchdir[0].y, &pvec[i].patchdir[0].z,
                        &pvec[i].switched);

        pvec[i].patchdir[1].x = pvec[i].patchdir[1].y = pvec[i].patchdir[1].z =0;
        pvec[i].chdir[0].x = pvec[i].chdir[0].y = pvec[i].chdir[0].z =0;
        pvec[i].chdir[1].x = pvec[i].chdir[1].y = pvec[i].chdir[1].z =0;
        DEBUG_INIT("Line: %s\nNumber of Fields: %d", line, fields);
        if (fields == 9){
            pvec[i].switched = 0;
            fprintf(stdout, "WARNING: Particle %u is assumed to be not switched!\n", i+1);
            fields++;
        }
        if (fields != 10) {
            fprintf (stderr, "ERROR: Could not read coordinates for particle %u.\n \
                    Did you specify box size at the begining?\n\n", i+1);
            free(line);
            exit (1);
        }
        /* Scale position vector to the unit cube */
#ifdef WEDGE
        pvec[i].pos.x /= conf->geo.box.x;
        pvec[i].pos.y /= conf->geo.box.y;
        pvec[i].pos.z /= conf->geo.box.z;

        conf->geo.usePBC(&pvec[i]);
#else
        pvec[i].pos.x /= conf->geo.box.x;
        pvec[i].pos.y /= conf->geo.box.y;
        pvec[i].pos.z /= conf->geo.box.z;

        // for compatibility unfortunately
        conf->geo.usePBC(&pvec[i]);
#endif

        if ((topo.ia_params[pvec[i].type][pvec[i].type].geotype[0]<SP)&&( DOT(pvec[i].dir, pvec[i].dir) < ZEROTOL )) {
            //DEBUG_INIT("Geotype = %d < %d", conf->pvec[i].geotype,SP);
            fprintf (stderr,
                    "ERROR: Null direction vector supplied for particle %u.\n\n", i+1);
            free(line);
            return false;
        } else {
            pvec[i].dir.normalise();
        }

        if ((topo.ia_params[pvec[i].type][pvec[i].type].geotype[0]<SP)&&( DOT(pvec[i].patchdir[0], pvec[i].patchdir[0]) < ZEROTOL )) {
            fprintf (stderr,
                    "ERROR: Null patch vector supplied for particle %u.\n\n", i+1);
            free(line);
            return false;
        } else {
            ortogonalise(&pvec[i].patchdir[0],&pvec[i].dir);
            pvec[i].patchdir[0].normalise();
        }
        // Switch the type
        if(pvec[i].switched){
            if(pvec[i].switchtype == 0){
                fprintf(stderr, "ERROR: Particle %u switched even though it has no switchtype", i);
                free(line);
                exit(1);
            }
            tmp_type = pvec[i].type;
            pvec[i].type = pvec[i].switchtype;
            pvec[i].switchtype = tmp_type;
        }

        DEBUG_INIT("%ld:\t%lf\t%lf\t%lf", i, pvec[i].pos.x, pvec[i].pos.y, pvec[i].pos.z);

    }
    free(line);
    /*Make chains WHOLE*/
    for (int i=0; i<conf->pvec.getChainCount(); i++){
        j=0;
        current = conf->pvec.getChainPart(i,0);
        first = current;
        chorig[0].pos = pvec[first].pos;
        while (current >=0 ) {
            /*shift the chain particle by first one*/
            pvec[current].pos.x -= chorig[0].pos.x;
            pvec[current].pos.y -= chorig[0].pos.y;
            pvec[current].pos.z -= chorig[0].pos.z;
            /*put it in orig geo.box*/
            pvec[current].pos.x -=  anInt(pvec[current].pos.x);
            pvec[current].pos.y -=  anInt(pvec[current].pos.y);
            pvec[current].pos.z -=  anInt(pvec[current].pos.z);
            //printf("ant: %f %f %f\n",conf->pvec[current].pos.x,conf->pvec[current].pos.y,conf->pvec[current].pos.z);
            /*shot it back*/
            pvec[current].pos.x += chorig[0].pos.x;
            pvec[current].pos.y += chorig[0].pos.y;
            pvec[current].pos.z += chorig[0].pos.z;
            //printf("posstart: %f %f %f\n",conf->pvec[current].pos.x,conf->pvec[current].pos.y,conf->pvec[current].pos.z);
            j++;
            current = conf->pvec.getChainPart(i,j);
        }
    }

    err = 0;
    //for (i=0; i < topo.npart-1; i++) {
    //    for (j=i+1; j < topo.npart; j++) {
    //        if ( overlap(conf->pvec[i], conf->particle[j], conf->geo.box, topo.ia_params) ) {
    //            fprintf (stderr,
    //                    "ERROR: Overlap in initial coniguration between particles %ld and %ld.\n",
    //                    i+1, j+1);
    //            err = 1;
    //        }
    //    }
    //}
    if (err) {
        printf ("\n");
        return false;
    }
    fflush (stdout);

    return true;
}


void Inicializer::testChains() {
    if (conf->pvec.getChainCount() == 0) {    // no chain -> make the probability of moving them 0
        if (sim->chainprob > 0)
            printf ("No chains... chain move probability set to 0.\n");
        sim->chainprob = 0.0;
    } else {
        for(int i=0; i<conf->pvec.molTypeCount; i++) {
            if(topo.moleculeParam[i].isGrandCanonical() && !topo.moleculeParam[i].isAtomic()) {
                if(!poolConfig) {
                    cout << "ChainInsert with no Pool system stated! State [Pool] in top.init" << endl;
                    exit(1);
                }
            }
        }
    }
}

void Inicializer::initWriteFiles() {
    sprintf(files->configurationPool, "pool");
    sprintf(files->configurationInFile, "config.init");
    sprintf(files->configurationoutfile, "config.last");
    sprintf(files->optionsfile, "options");
    sprintf(files->topologyInFile, "top.init");
    sprintf(files->topologyOutFile, "top.last");
    sprintf(files->moviefile, "movie");
    sprintf(files->wlinfile, "wl.dat");
    sprintf(files->wloutfile, "wl-new.dat");
    sprintf(files->statfile, "stat.dat");
    sprintf(files->clusterfile, "cluster.dat");
    sprintf(files->clusterstatfile, "cluster_stat.dat");
    sprintf(files->energyfile, "energy.dat");
}

void Inicializer::initNeighborList() {
    cout << "\nAllocating memory for pairlist..." << endl;
    //sim->pairlist = (Pairs*) xMalloc(sizeof(Pairs) * topo.npart); // deprecated, see Conf->neighborlist

    // Highest guess: Every particle interacts with the others
    // TODO: Make it more sophisticated
    conf->neighborList.resize(conf->pvec.size());
    for(unsigned long i = 0; i < conf->neighborList.size(); i++){
        conf->neighborList[i].neighborID = (long int*) malloc(sizeof(long) * MAXN);
        conf->neighborList[i].neighborCount = 0;
    }
}

void Inicializer::initMPI(int argc, char** argv) {
#ifdef ENABLE_MPI
        FILE *infile;
        printf(" MPI version");
        MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD, &(sim->mpinprocs) );
        MPI_Comm_rank(MPI_COMM_WORLD, &(sim->mpirank) );
        sim->pseudoRank = sim->mpirank;

        // MPI out files
        files->initMPIRank(sim->mpirank);

        //test if there is a specific input configuration for mpi run

        infile = fopen(files->configurationInFile, "r");
        if (infile != NULL)
            fclose (infile);
        else  sprintf(files->configurationInFile, "config.init");

        //test if there is a specific input wang-landau for mpi run

        infile = fopen(files->wlinfile, "r");
        if (infile != NULL)
            fclose (infile);
        else  sprintf(files->wlinfile, "wl.dat");
#endif
}

void Inicializer::initGroupLists() {

    if(SILENT == 1)
        cout << "Generating GroupLists..." << endl;

    // setGroupList;
    int type=0;
    while(topo.moleculeParam[type].name != NULL) { // get all types
        for(unsigned int i = 0; i < conf->pvec.size(); i++) {
            if(type < conf->pvec[i].molType) { // note: we arent searching for molType of particle, could be 0 particles
                type = conf->pvec[i].molType;
                conf->pvec.first[type] = i;
                break;
            }
            // FIX for situation
            // A X>0
            // B 0
            if(type !=0 && (i+1) == conf->pvec.size())
                conf->pvec.first[type] = conf->pvec.size();
        }
        type++;
    }
    conf->pvec.molTypeCount = type;
    conf->pvec.first[type] = conf->pvec.size();
    conf->pvec.calcChainCount();

    //test grouplist consistency
    /*int size=0;
    for(int i=0; i < type; i++) {
        size = 0;
        for(unsigned int j=0; j<conf->pvec.size(); j++) {
            if(i == conf->pvec[j].molType)
                size++;
        }
        cout << size << "==" << conf->pvec.molCountOfType(i) << endl;
    }*/

    int newType = -1;
    for(unsigned int i = 0; i < conf->pool.size(); i++) {
        // set simple grouplist
        if(newType != conf->pool[i].molType) {
            newType = conf->pool[i].molType;
            conf->pool.first[newType] = i;
        }
    }
    conf->pool.molTypeCount = newType+1;
    conf->pool.first[newType+1] = conf->pool.size();
}





/************************************************************************************************
 *                                      PRIVATE METHODS                                         *
 ************************************************************************************************/



void Inicializer::setParticlesParamss(MolIO* molecules, long *sysmoln, char **sysnames, std::vector<Particle> *pvec) {
    long i=0, j=0, mol, k, maxpart=0;

    while (sysnames[i]!=NULL) {
        mol=0;
        while (strcmp(molecules[mol].name,sysnames[i]) && mol < MAXMT) {
            mol++;
            if (molecules[mol].name == NULL) {
                fprintf (stderr, "TOPOLOGY ERROR: molecules %s is not defined.\n\n",sysnames[i]);
                topDealoc();
                exit(1);
            }
        }

        for (j=0;j<sysmoln[i];j++) {
            //DEBUG	    fprintf (stdout, "molnames %s sysname %s sysnum %ld \n",molnames[mol],sysnames[i],sysmoln[i]);
            k=0;
            while (molecules[mol].type[k] != -1) {

                pvec->push_back(Particle());
                (*pvec)[maxpart].type        = molecules[mol].type[k];
                (*pvec)[maxpart].switchtype  = molecules[mol].switchtype[k];
                (*pvec)[maxpart].delta_mu    = molecules[mol].delta_mu[k];
                (*pvec)[maxpart].molType     = mol;
                //(*pvec)[maxpart].chainIndex  = maxch;

                k++;
                maxpart++;

                if (maxpart > MAXN) {
                    fprintf (stderr, "TOPOLOGY ERROR: more particles(%ld) than allowed(%d).\n",maxpart,MAXN);
                    fprintf (stderr, "Change MAXN in source and recompile the program. \n\n");
                    topDealoc();
                    exit(1);
                }
            }
        }
        i++;
    }

    assert(maxpart == (long)pvec->size());
}


void Inicializer::readTopoFile(bool exclusions[][MAXT]) {
    char *dummy=NULL;
    char line[STRLEN], keystr[STRLEN], molname[STRLEN];
    unsigned size;
    long i=0;
    FILE *infile;
    char *pline=NULL;

    if ((infile = fopen(files->topologyInFile, "r")) == NULL) {
        fprintf (stderr, "\nTOPOLOGY ERROR: Could not open top.init file.\n\n");
        exit (1);
    }

    if(SILENT == 1)
        cout << "Reading topology...\n"
             << "Species:" << endl;

    molname[0] = ' ';

    pline = (char*) malloc((size_t)STRLEN);
    while (fgets2(line,STRLEN-2,infile) != NULL) {
        strcpy(pline,line);
        if (!pline) fprintf (stderr, "\nTOPOLOGY ERROR: Empty line in topology.\n\n");

        // build one long line from several fragments
        while (continuing(line) && (fgets2(line,STRLEN-1,infile) != NULL)) {
            size=strlen(pline)+strlen(line)+1;
            free(pline);
            pline = (char*) malloc((size_t)size);
            strcat(pline,line);
        }

        strip_comment (pline);
        trim (pline);

        if ((int)strlen(pline) > 0) {
            // get the [COMMAND] key
            if (pline[0] == OPENKEY) {
                pline[0] = ' ';
                beforecommand(keystr,pline,CLOSEKEY);
                upstring (keystr);
            } else {

                //DEBUG		fprintf (stdout, "Topology read type:%s, %s \n",keystr,pline);
                if (!strcmp(keystr,"TYPES")) {
                    fflush(stdout);
                    if (!fillTypes(&pline)) {
                        DEBUG_INIT("Something went wrong with filltypes");
                        fprintf (stderr, "\nTOPOLOGY ERROR: in reading types\n\n");
                        topDealoc();
                        free(pline); pline = NULL;
                        exit (1);
                    }
                    DEBUG_INIT("back in init_top");
                    continue;
                }
                if (!strcmp(keystr,"MOLECULES")) {
                    DEBUG_INIT("Let's go to the molecules");
                    if (molname[0] == ' ') {
                        beforecommand(molname,pline,SEPARATOR);
                        i=0;
                        while (molecules[i].name != NULL)
                            i++;
                        DEBUG_INIT("in the middle of getting to fillmol");
                        molecules[i].name = (char*) malloc(strlen(molname)+1);
                        strcpy(molecules[i].name, molname);
                        fprintf (stdout, "\nTopology read for molecule: %s \n",molname);
                    }
                    if (!fillMol(molname, pline, molecules)) {
                        fprintf (stderr, "\nTOPOLOGY ERROR: in reading molecules\n\n");
                        topDealoc();
                        free(pline); pline = NULL;
                        exit (1);
                    }
                    if ((dummy = strchr (pline,CLOSEMOL)) != NULL)
                        molname[0] = ' ';
                    continue;
                }
                if (!strcmp(keystr,"SYSTEM")) {
                    char name[9] = "system: ";
                    if (!fillSystem(pline,sysnames,&sysmoln, name)) {
                        fprintf (stderr, "\nTOPOLOGY ERROR: in reading system\n\n");
                        topDealoc();
                        free(pline); pline = NULL;
                        exit (1);
                    }
                    continue;
                }
                if (!strcmp(keystr, "POOL")) {
                    poolConfig = true;
                    char name[7] = "pool: ";
                    if (!fillSystem(pline,poolNames,&poolMolNum, name)) {
                        fprintf (stderr, "\nTOPOLOGY ERROR: in reading system\n\n");
                        topDealoc();
                        free(pline); pline = NULL;
                        exit (1);
                    }
                    continue;
                }
                if (!strcmp(keystr,"EXTER")) {
                    fflush(stdout);
                    if (!fillExter(&pline)) {
                        DEBUG_INIT("Something went wrong with external potential");
                        fprintf (stderr, "\nTOPOLOGY ERROR: in reading external potential\n\n");
                        topDealoc();
                        free(pline); pline = NULL;
                        exit (1);
                    }
                    continue;
                }
                if (!strcmp(keystr,"EXCLUDE")) {
                    fflush(stdout);
                    if (!fillExclusions(&pline,exclusions)) {
                        DEBUG_INIT("Something went wrong with exclusions potential");
                        fprintf (stderr, "\nTOPOLOGY ERROR: in reading exclusions\n\n");
                        topDealoc();
                        free(pline); pline = NULL;
                        exit (1);
                    }
                    continue;
                }

                fprintf (stderr, "\nTOPOLOGY ERROR: invalid keyword:%s.\n\n", keystr);
                topDealoc();
                free(pline); pline = NULL;
                exit (1);
            }
        }
    }
    //we have sucessfully read topology
    if (pline !=NULL) free(pline);
    pline=NULL;
    fclose (infile);
    fflush (stdout);
}



void *Inicializer::xMalloc(size_t num) {
    void *neww = malloc (num);
    if (!neww){
        fprintf(stderr, "Couldn't allocate any memory!\n");
        exit(1);
    }
    return neww;
}


int Inicializer::fillExclusions(char **pline, bool exlusions[][MAXT]) {
    long num1,num2;
    char *pline1, *pline2;

    num1 = strtol(*pline, &pline2, 10);
    trim(pline2);
    if ((int)strlen(pline2) > 0) {
        num2 = strtol(pline2, &pline1, 10);
        trim(pline1);
        exlusions[num1][num2]=true;
        exlusions[num2][num1]=true;
    } else {
        fprintf(stderr, "Error in readin Topology exclusions, probably there is not even number of types \n");
        return 0;
    }
    while ((int)strlen(pline1) > 0) {
      num1 = strtol(pline1, &pline2, 10);
      trim(pline2);
      if ((int)strlen(pline2) > 0) {
        num2 = strtol(pline2, &pline1, 10);
        trim(pline1);
        exlusions[num1][num2]=true;
        exlusions[num2][num1]=true;
      } else {
        fprintf(stderr, "Error in readin Topology exclusions, probably there is not even number of types \n");
        return 0;
      }
    }
    return 1;
}


int Inicializer::fillSystem(char *pline, char *sysnames[], long **sysmoln, char* name) {
    int i,fields;
    char zz[STRLEN];

    trim(pline);
    if (!pline) {
        fprintf (stderr, "TOPOLOGY ERROR: obtained empty line in fil system.\n\n");
        return 0;
    }
    i=0;
    while (sysnames[i]!=NULL) i++;

    fields = sscanf(pline, "%s %ld", zz, &(*sysmoln)[i]);
    sysnames[i] = (char*) malloc(strlen(zz)+1);
    strcpy(sysnames[i],zz);

    if (fields != 2) {
        fprintf (stderr, "TOPOLOGY ERROR: failed reading system from (%s).\n\n", pline);
        return 0;
    }
    /*if ((*sysmoln)[i] < 1) {
        fprintf (stderr, "TOPOLOGY ERROR: cannot have %ld number of molecules.\n\n", (*sysmoln)[i]);
        return 0;
    }*/
    fprintf (stdout, "%s %s %ld\n",name, sysnames[i],(*sysmoln)[i]);
    return 1;
}


int Inicializer::fillTypes(char **pline) {
    int type;
    int geotype_i;
    int fields;
    char name[SMSTR];
    char geotype[SMSTR];

    double param[12];
    /* 0: epsilon
     * 1: sigma
     * 2: attraction dist
     * 3: sttraction switch
     * 4: patch angle
     * 5: patch switch
     * 6: length
     * 7: parallel_eps
     * 8(optional): second patche rotation
     * 9(optional): second patch angle
     * 10(optional): second patch angle switch
     * +1: chirality
     */
    char typestr[STRLEN], paramstr[STRLEN];

    beforecommand(typestr, *pline, SEPARATOR);
    aftercommand(paramstr, *pline, SEPARATOR);

    fields = sscanf(paramstr, "%s %d %s %le %le %le %le %le %le %le %le %le %le %le %le",
                    name, &type, geotype, &param[0], &param[1], &param[2], &param[3], &param[4],
                    &param[5], &param[6], &param[7], &param[8], &param[9], &param[10], &param[11]);

    cout << "Fields:" << fields << " " << param[11] << endl;

    fields -= 5; // number of parameter fields => I am too lazy to adjust everywhere below the numbers
    //DEBUG    fprintf (stdout, "Topology read geotype: %ld with parameters fields %d, str:%s and %s in pline %s\n",geotype,fields,geotypestr,paramstr,pline);

    geotype_i = convertGeotype(geotype);
    if(!geotype_i){
        fprintf(stderr, "TOPOLOGY ERROR: Unknown GEOTYPE: %s!", geotype);
        return 0;
    }
    DEBUG_INIT("geotype_i: %d; fields = %d", geotype_i, fields);
    if (( (geotype_i == SPN) ) && (fields != 0)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA" << endl;
        return 0;
    }
    if (( (geotype_i == SCN) ) && (fields != 1)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA SC_LENGTH" << endl;
        return 0;
    }
    if (( (geotype_i == SPA)) && (fields != 2)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA ATTRACT_DIST ATTRACT_SWITCH" << endl;
        return 0;
    }
    if (( (geotype_i == SCA) ) && (fields != 3)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA ATTRACT_DIST ATTRACT_SWITCH SC_LENGTH" << endl;
        return 0;
    }
    if (( (geotype_i == PSC) || (geotype_i == CPSC) ) && (fields != 6)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA ATTRACT_DIST ATTRACT_SWITCH PATCH_ANGLE PATCH_SWITCH SC_LENGTH PARALLEL_EPS" << endl;
        return 0;
    }
    if (( (geotype_i == CHCPSC) || (geotype_i == CHCPSC) )&& ( fields != 7)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA ATTRACT_DIST ATTRACT_SWITCH PATCH_ANGLE PATCH_SWITCH SC_LENGTH PARALLEL_EPS CHIRAL_ANGLE" << endl;
        return 0;
    }
    if (( (geotype_i == TPSC) || (geotype_i == TCPSC) ) && (fields != 9)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA ATTRACT_DIST ATTRACT_SWITCH PATCH_ANGLE PATCH_SWITCH SC_LENGTH PARALLEL_EPS  PATCH_ROTATION PATCH_ANGLE PATCH_SWITCH" << endl;
        return 0;
    }
    if (( (geotype_i == TCHCPSC) || (geotype_i == TCHCPSC) )&& ( fields != 10)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA ATTRACT_DIST ATTRACT_SWITCH PATCH_ANGLE PATCH_SWITCH SC_LENGTH PARALLEL_EPS  PATCH_ROTATION PATCH_ANGLE PATCH_SWITCH CHIRAL_ANGLE" << endl;
        return 0;
    }

    if ((geotype_i < 0) || (geotype_i > (MAXT + 10))) {
        fprintf (stderr, "TOPOLOGY ERROR: geotype (%s) is out of range: 0 - %d.\n\n", geotype, MAXT + 10);
        return 0;
    }

    strcpy(topo.ia_params[type][type].name, name);
    strcpy(topo.ia_params[type][type].other_name, name);
    topo.ia_params[type][type].geotype[0] = geotype_i;
    topo.ia_params[type][type].geotype[1] = geotype_i;
    topo.ia_params[type][type].epsilon = param[0];
    topo.ia_params[type][type].sigma = param[1];
    topo.ia_params[type][type].rcutwca = (topo.ia_params[type][type].sigma)*pow(2.0,1.0/6.0);

    fprintf(stdout, "Topology read of %d: %8s (geotype: %s, %d) with parameters %g %g", type, name, geotype, geotype_i, topo.ia_params[type][type].epsilon, topo.ia_params[type][type].sigma);

    if (fields > 0 && fields != 1 && fields != 3) { // all except SCN and SCA
        topo.ia_params[type][type].pdis = param[2];
        topo.ia_params[type][type].pswitch = param[3];
        topo.ia_params[type][type].rcut = topo.ia_params[type][type].pswitch+topo.ia_params[type][type].pdis;
        fprintf(stdout, " | %g %g",topo.ia_params[type][type].pdis,topo.ia_params[type][type].pswitch);
    }
    if(fields == 1) { // SCN
        for(int i = 0; i < 2; i++){
            topo.ia_params[type][type].len[i] = param[2];
            topo.ia_params[type][type].half_len[i] = param[2] / 2;
        }
    }
    if(fields == 3) { // SCA
        for(int i = 0; i < 2; i++){
            topo.ia_params[type][type].len[i] = param[4];
            topo.ia_params[type][type].half_len[i] = param[4] / 2;
        }
    }
    if (fields > 2 && fields != 3) { // except SCA
        int i;
        for(i = 0; i < 2; i++){
            topo.ia_params[type][type].len[i] = param[6];
            topo.ia_params[type][type].half_len[i] = param[6] / 2;
            topo.ia_params[type][type].pangl[i] = param[4];
            topo.ia_params[type][type].panglsw[i] = param[5];
            topo.ia_params[type][type].pcangl[i] = cos(param[4]/2.0/180*PI);                 // C1
            topo.ia_params[type][type].pcanglsw[i] = cos((param[4]/2.0+param[5])/180*PI);    // C2
            //topo.ia_params[type][type].pcangl[i] = topo.ia_params[type][type].pcangl[i];
            //topo.ia_params[type][type].pcanglsw[i] =	topo.ia_params[type][type].pcanglsw[i];
            topo.ia_params[type][type].pcoshalfi[i] = cos((param[4]/2.0+param[5])/2.0/180*PI);
            topo.ia_params[type][type].psinhalfi[i] = sqrt(1.0 - topo.ia_params[type][type].pcoshalfi[i] * topo.ia_params[type][type].pcoshalfi[i]);
            topo.ia_params[type][type].parallel = param[7];
	  
	}
        fprintf(stdout, " | %g %g | %g", topo.ia_params[type][type].pangl[0], topo.ia_params[type][type].panglsw[0], topo.ia_params[type][type].parallel);
    }
    if(fields == 7){
        int i;
        for(i = 0; i < 2; i++){
            topo.ia_params[type][type].chiral_cos[i] = cos(param[8] / 360 * PI);
            topo.ia_params[type][type].chiral_sin[i] = sqrt(1 - topo.ia_params[type][type].chiral_cos[i] * topo.ia_params[type][type].chiral_cos[i]);
            fprintf(stdout, "| chirality %g ", param[8]);
        }
    }
    if ((fields == 9)||(fields == 10)) {
        int i;
        for(i = 0; i < 2; i++){
            topo.ia_params[type][type].csecpatchrot[i] = cos(param[8] / 360 * PI);
            topo.ia_params[type][type].ssecpatchrot[i] = sqrt(1 - topo.ia_params[type][type].csecpatchrot[i] * topo.ia_params[type][type].csecpatchrot[i]);
            //fprintf(stdout, " | %g %g", topo.ia_params[type][type].csecpatchrot[0], topo.ia_params[type][type].ssecpatchrot[0]);

            topo.ia_params[type][type].pangl[i+2] = param[9];
            topo.ia_params[type][type].panglsw[i+2] = param[10];
            topo.ia_params[type][type].pcangl[i+2] = cos(param[9]/2.0/180*PI);                 // C1
            topo.ia_params[type][type].pcanglsw[i+2] = cos((param[9]/2.0+param[10])/180*PI);    // C2
            //topo.ia_params[type][type].pcangl[i] = topo.ia_params[type][type].pcangl[i];
            //topo.ia_params[type][type].pcanglsw[i] = topo.ia_params[type][type].pcanglsw[i];
            topo.ia_params[type][type].pcoshalfi[i+2] = cos((param[9]/2.0+param[10])/2.0/180*PI);
            topo.ia_params[type][type].psinhalfi[i+2] = sqrt(1.0 - topo.ia_params[type][type].pcoshalfi[i+2] * topo.ia_params[type][type].pcoshalfi[i+2]);
        }
        fprintf(stdout, " | %g  %g %g", param[8], topo.ia_params[type][type].pangl[2], topo.ia_params[type][type].panglsw[2]);
    }
    if(fields == 10){
        int i;
        for(i = 0; i < 2; i++){
            topo.ia_params[type][type].chiral_cos[i] = cos(param[11] / 360 * PI);
            topo.ia_params[type][type].chiral_sin[i] = sqrt(1 - topo.ia_params[type][type].chiral_cos[i] * topo.ia_params[type][type].chiral_cos[i]);
        }
        fprintf(stdout, " | %g ", param[11]);
    }

    // Volume
    if (geotype_i < SP)
        topo.ia_params[type][type].volume = 4.0/3.0*PI*pow((topo.ia_params[type][type].sigma)/2.0,3.0) + PI/2.0*topo.ia_params[type][type].len[0]*pow((topo.ia_params[type][type].sigma)/2.0,2.0) ;
    else
        topo.ia_params[type][type].volume = 4.0/3.0*PI*pow((topo.ia_params[type][type].sigma)/2.0,3.0);
    if ( topo.ia_params[type][type].rcutwca > topo.sqmaxcut )
        topo.sqmaxcut = topo.ia_params[type][type].rcutwca;
    if ( topo.ia_params[type][type].rcut > topo.sqmaxcut )
        topo.sqmaxcut = topo.ia_params[type][type].rcut;
    fprintf(stdout, " \n");
    DEBUG_INIT("Finished filltypes");
    return 1;
}


int Inicializer::convertGeotype(char *geotype) {
//    if (strcmp(geotype, "SC") == 0)
//        return SC;
    if (strcmp(geotype, "SCN") == 0)
        return SCN;
    if (strcmp(geotype, "SCA") == 0)
        return SCA;
    if (strcmp(geotype, "PSC") == 0)
        return PSC;
    if (strcmp(geotype, "CPSC") == 0)
        return CPSC;
    if (strcmp(geotype, "CHPSC") == 0)
        return CHPSC;
    if (strcmp(geotype, "CHCPSC") == 0)
        return CHCPSC;
    if (strcmp(geotype, "TPSC") == 0)
        return TPSC;
    if (strcmp(geotype, "TCPSC") == 0)
        return TCPSC;
    if (strcmp(geotype, "TCHPSC") == 0)
        return TCHPSC;
    if (strcmp(geotype, "TCHCPSC") == 0)
        return TCHCPSC;
    if (strcmp(geotype, "SP") == 0)
        return SP;
    if (strcmp(geotype, "SPN") == 0)
        return SPN;
    if (strcmp(geotype, "SPA") == 0)
        return SPA;
    return 0;
}


int Inicializer::fillExter(char **pline) {
    int fields;

    double param[3];
    /* 0: thickness
     * 1: epsilon
     * 2: attraction
     */
    char typestr[STRLEN], paramstr[STRLEN];

    beforecommand(typestr, *pline, SEPARATOR);
    aftercommand(paramstr, *pline, SEPARATOR);
    fields = sscanf(paramstr, "%le %le %le", &param[0], &param[1], &param[2]);
    if (fields >3) {
        fprintf (stderr, "TOPOLOGY ERROR: too many parameters for external potential. We have \
                thickness, epsilon, and attraction distance so far.\n\n");
        return 0;
    }
    if (fields >0) {
        topo.exter.exist = true;
        topo.exter.thickness = param[0];
        fprintf(stdout, "External potential with thickness: %le ",topo.exter.thickness);
        if (fields >1) {
            topo.exter.epsilon = param[1];
            fprintf(stdout, "epsilon: %le ",topo.exter.epsilon);
            if (fields >2) {
                topo.exter.attraction = param[2];
                fprintf(stdout, "and range of attraction: %le ",topo.exter.attraction);
            }
        }
    } else{
        topo.exter.exist = false;
        fprintf(stdout, "No external potential ");
    }

    fprintf(stdout, " \n");
    DEBUG_INIT("Finished filling external potential");
    return 1;
}


int Inicializer::fillMol(char *molname, char *pline, MolIO *molecules) {
    DEBUG_INIT("fillmol just has been called!");
    char str[STRLEN],str2[STRLEN],molcommand[STRLEN],molparams[STRLEN];
    int i,j,fields;
    double bondk,bonddist, activity;
    const double Nav = 6.022137e23;

    beforecommand(str2, pline, CLOSEMOL);
    aftercommand(str, str2, OPENMOL);
    trim(str);

    if (strlen(str) == 0) return 1;
    beforecommand(molcommand,str,SEPARATOR);
    aftercommand(molparams,str,SEPARATOR);
    trim(molcommand);
    trim(molparams);
    upstring(molcommand);
    DEBUG_INIT("molcommand: %s", molcommand);
    DEBUG_INIT("molparams: %s", molparams);
    i=0;
    while (strcmp(molecules[i].name, molname)) // number of molTypes already loaded
        i++;
    j=0;
    while (molecules[i].type[j] != -1) // number of particles of this molType loaded
        j++;

    if (!strcmp(molcommand,"PARTICLES")) {
        fprintf (stdout, "particle %d: \t", j + 1);
        fields =  sscanf(molparams,"%d %ld %lf",molecules[i].type + j,
                         molecules[i].switchtype + j, molecules[i].delta_mu + j);
        fprintf (stdout, "%d ",molecules[i].type[j]);

        if(j==0) {
            topo.moleculeParam[i].name = (char*) malloc(strlen(molname)+1);
            strcpy(topo.moleculeParam[i].name, molname);
        }

        topo.moleculeParam[i].particleTypes.push_back(molecules[i].type[j]);
        assert(topo.moleculeParam[i].particleTypes[j] == molecules[i].type[j]);

        if (fields == 1){
                (molecules[i].switchtype[j]) = -1;//(molecules[i].type[j]);
                (molecules[i].delta_mu[j]) = 0;
                fields = 3;
        } else{
            fprintf(stdout, "(with switchtype: %ld and delta_mu: %lf)", molecules[i].switchtype[j], molecules[i].delta_mu[j]);
            topo.moleculeParam[i].switchTypes.push_back(molecules[i].switchtype[j]);
            topo.moleculeParam[i].deltaMu.push_back(molecules[i].delta_mu[j]);
        }
        if (fields != 3) {
            fprintf (stderr, "TOPOLOGY ERROR: could not read a pacticle.\n\n");
            return 0;
        }

        if (molecules[i].type[j] < 0) {
            fprintf (stderr, "TOPOLOGY ERROR: pacticles include negative type.\n\n");
            return 0;
        }
        if (molecules[i].type[j] > MAXT) {
            fprintf (stderr, "TOPOLOGY ERROR: pacticles include type out of range 0-%ld.\n\n",(long)MAXT);
            return 0;
        }
        fprintf (stdout, "\n");
        return 1;
    }
    if (!strcmp(molcommand,"BOND1")) {
        fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
        if (fields < 2) {
            fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for bond1, should be 2.\n\n");
            return 0;
        }
        if (bonddist < 0) {
            fprintf (stderr, "TOPOLOGY ERROR: bonddist cannot be negative: %f \n\n",bonddist);
            return 0;
        }
        topo.moleculeParam[i].bond1c = bondk;
        topo.moleculeParam[i].bond1eq = bonddist;
        fprintf (stdout, "bond1: %f %f \n",topo.moleculeParam[i].bond1c,topo.moleculeParam[i].bond1eq);
        return 1;
    }
    if (!strcmp(molcommand,"BOND2")) {
        fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
        if (fields < 2) {
            fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for bond2, should be 2.\n\n");
            return 0;
        }
        if (bonddist < 0) {
            fprintf (stderr, "TOPOLOGY ERROR: bonddist cannot be negative: %f \n\n",bonddist);
            return 0;
        }
        topo.moleculeParam[i].bond2c = bondk;
        topo.moleculeParam[i].bond2eq = bonddist;
        fprintf (stdout, "bond2: %f %f \n",topo.moleculeParam[i].bond2c,topo.moleculeParam[i].bond2eq);
        return 1;
    }
    if (!strcmp(molcommand,"BONDD")) {
        fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
        if (fields < 2) {
            fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for bondd, should be 2.\n\n");
            return 0;
        }
        if (bonddist < 0) {
            fprintf (stderr, "TOPOLOGY ERROR: bonddist cannot be negative: %f \n\n",bonddist);
            return 0;
        }
        topo.moleculeParam[i].bonddc = bondk;
        topo.moleculeParam[i].bonddeq = bonddist;
        fprintf (stdout, "bondd: %f %f \n",topo.moleculeParam[i].bonddc,topo.moleculeParam[i].bonddeq);
        return 1;
    }

    if (!strcmp(molcommand,"BONDH")) {
        fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
        if (fields < 2) {
            fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for bondh, should be 2.\n\n");
            return 0;
        }
        if (bonddist < 0) {
            fprintf (stderr, "TOPOLOGY ERROR: bonddist cannot be negative: %f \n\n",bonddist);
            return 0;
        }
        topo.moleculeParam[i].bondhc = bondk;
        topo.moleculeParam[i].bondheq = bonddist;
        fprintf (stdout, "bondd: %f %f \n",topo.moleculeParam[i].bondhc,topo.moleculeParam[i].bondheq);
        return 1;
    }

    if (!strcmp(molcommand,"ANGLE1")) {
        fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
        if (fields < 2) {
            fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for angle1, should be 2.\n\n");
            return 0;
        }
        if (bonddist < 0) {
            fprintf (stderr, "TOPOLOGY ERROR: equilibrium angle cannot be negative: %f \n\n",bonddist);
            return 0;
        }
        topo.moleculeParam[i].angle1c = bondk;
        topo.moleculeParam[i].angle1eq = bonddist/180.0*PI;
        fprintf (stdout, "angle1: %f %f \n",topo.moleculeParam[i].angle1c,topo.moleculeParam[i].angle1eq);
        return 1;
    }
    if (!strcmp(molcommand,"ANGLE2")) {
        fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
        if (fields < 2) {
            fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for angle2, should be 2.\n\n");
            return 0;
        }
        if (bonddist < 0) {
            fprintf (stderr, "TOPOLOGY ERROR: equilibrium angle cannot be negative: %f \n\n",bonddist);
            return 0;
        }
        topo.moleculeParam[i].angle2c = bondk;
        topo.moleculeParam[i].angle2eq = bonddist/180.0*PI;
        fprintf (stdout, "angle2: %f %f \n",topo.moleculeParam[i].angle2c,topo.moleculeParam[i].angle2eq);
        return 1;
    }

    // INIT of muVT ensemble
    if (!strcmp(molcommand,"ACTIVITY")) {
        if(sim->nGrandCanon == 0) {
            cout << "Activity stated in top.init, But nGrandCanon=0 in options" << endl;
            exit(1);
        }
        fields = sscanf(molparams, "%le ", &activity);
        topo.moleculeParam[i].activity = activity;
        topo.moleculeParam[i].chemPot = log(activity*Nav*1e-24); // faunus log(activity*Nav*1e-27) [mol/l]
        fprintf (stdout, "activity: %f \n",topo.moleculeParam[i].activity);
        return 1;
    }

    fprintf (stderr, "TOPOLOGY ERROR: unknown parameter: %s.\n\n",molcommand);
    return 0;
}


void Inicializer::readii2(char * num, int value[2]) {
    char *end, *num2;

    value[0] = strtol(num, &num2, 10);
    trim(num2);
    if ((int)strlen(num2) > 0)
    value[1] = strtol(num2, &end, 10);
    else {
    value[1] =0;
    return;
    }
    if(*end){
        fprintf(stderr, "Could not convert %s into two integers\n", num);
        exit(1);
    }

    return;
}




double Inicializer::readd2(char * num) {
    char *end;
    double i = strtod(num, &end);
    if(*end){
        fprintf(stderr, "Could not convert %s into double\n", num);
        exit(1);
    }
    return i;
}




long Inicializer::readl2(char * num) {
    char *end;
    long i = strtol(num, &end, 10);
    if(*end){
        fprintf(stderr, "Could not convert %s into long\n", num);
        exit(1);
    }
    return i;
}




int Inicializer::readi2(char * num) {
    char *end;
    int i = strtol(num, &end, 10);
    if(*end){
        fprintf(stderr, "Could not convert %s into integer\n", num);
        exit(1);
    }
    return (int) i;
}
