#include "inicializer.h"

#ifdef MPI
extern MPI_Datatype MPI_vector, MPI_Particle, MPI_exchange;
#endif

void Inicializer::readOptions() {
    cout << "Reading options..." << endl;

    int num_options = -1;
    double transmx, rotmx, chainmmx, chainrmx, angle, chain_angle;

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
        {"paramfrq",            Long,   false, &sim->paramfrq},
        {"report",              Long,   false, &sim->report},
        {"seed",                Long,   false, &seed},
        {"pairlist_update",     Int,    false, &sim->pairlist_update},
        {"ptype",               Int,    false, &sim->ptype},
        {"wlm",                 Int2,   false, &sim->wlm},
        {"wlmtype",             Int,    false, &sim->wl.wlmtype},
        {"press",               Double, false, &sim->press},
        {"paralpress",          Double, false, &sim->paralpress},
        {"edge_mx",             Double, false, &sim->edge.mx},
        {"shave",               Double, false, &sim->shave},
        {"chainprob",           Double, false, &sim->chainprob},
        {"switchprob",          Double, false, &sim->switchprob},
        {"temper",              Double, false, &sim->temper},
        {"paraltemper",         Double, false, &sim->paraltemper},
        {"transmx",             Double, false, &transmx},
        {"rotmx",               Double, false, &rotmx},
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
        if(strlen(line) == 0){
            continue;
        }
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
    printf (" Inititial maximum box edge change:                  %.8f\n", sim->edge.mx);
    printf (" Initial maximum chain displacement:                 %.8f\n", chainmmx);
    printf (" Inititial maximum chain angular change (degrees):   %.8f\n", chainrmx);
    printf (" Temperature in kT/e:                                %.8f\n", sim->temper);
    printf (" Parallel tempering temperature in kT/e:             %.8f\n", sim->paraltemper);
    printf (" Sweeps between replica exchange:                    %ld\n", sim->nrepchange);
    printf (" Sweeps between muVT insert/delete move:             %ld\n", sim->nGrandCanon);
    printf (" Wang-Landau method:                                 %d %d\n", sim->wlm[0],sim->wlm[1]);
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
    if ( (sim->wlm[0] <0) || (sim->wlm[0] > 7) || (sim->wlm[1] <0) || (sim->wlm[1] > 7)  ) {
        fprintf (stderr, "ERROR: Unknown Wang-Landau method %d %d. Program only knows: 0 - none, \
                1 - z-direction od 1st particle, 2 - pore in membrane, 3 - zorientation of 0th particle,\
                4 - distance of fist two particles, 5 - pore around z-axis above CM,\
                6 - pore around z-axis above 0th particle, 7 - number of particles in contact \n\n",sim->wlm[0],sim->wlm[1]);
        exit (1);
    }
    if ( (sim->wlm[0] == 0) && (sim->wlm[1] > 0)  ) {
        fprintf (stderr, "ERROR: Wang-Landau method has to be set for first order parameter and then for second order parameter\n\n");
        exit (1);
    }
    if ( (sim->wlm[0] == 2) || (sim->wlm[0] == 5) || (sim->wlm[0] == 6)  ) {
        if(sim->wl.wlmtype < 1){
            fprintf (stderr, "ERROR: Atom type for the Wang-Landau Method (%d) was false defined.\n\n",sim->wl.wlmtype);
            exit (1);
        }
        if ( (sim->wlm[1] == 2) || (sim->wlm[1] == 5) || (sim->wlm[1] == 6) ) {
            fprintf (stderr, "ERROR: Simulaneous use of two pore order parameters has not been implemented yet.\n\n");
            exit (1);
        }
    }

    // we store maximum rotation as half angle - useful for quaterions
    angle = rotmx / 180.0 * PIH *0.5;
    rotmx = cos((rotmx)/180.0*PIH);
    chain_angle = chainrmx / 180.0 * PIH;
    chainrmx = cos((chainrmx)/180.0*PIH);
    sim->edge.mx *= 2.0;   // The full range is -maxl to +maxl, i.e. spanning 2*maxl
    transmx *= 2.0;   // The full range is -maxr to +maxr, i.e. spanning 2*maxr
    chainmmx *= 2.0;   // The full range is -maxr to +maxr, i.e. spanning 2*maxr

    for (int i=0;i<MAXT;i++) {
        sim->trans[i].mx = transmx;
        sim->rot[i].mx = rotmx;
        sim->rot[i].angle = angle;
    }
    for (int i=0;i<MAXMT;i++) {
        sim->chainm[i].mx = chainmmx;
        sim->chainr[i].mx = chainrmx;
        sim->chainr[i].angle = chain_angle;
    }

    //parallel tempering
#ifdef MPI
    if ( (sim->temper != sim->paraltemper) && (sim->mpinprocs <2) ) {
        printf("ERROR: Paralllel tempering at single core does not work.\n\n");
        exit(1);
    }
    sim->dtemp = (sim->paraltemper - sim->temper )/(sim->mpinprocs-1);
    sim->temper += sim->dtemp * sim->mpirank;
    if ( (sim->press != sim->paralpress) && (sim->mpinprocs <2) ) {
        printf("ERROR: Pressure replica exchange at single core does not work.\n\n");
        exit(1);
    }
    sim->dpress = (sim->paralpress - sim->press )/(sim->mpinprocs-1);
    sim->press += sim->dpress * sim->mpirank;
    seed += sim->mpirank;
    sim->mpiexch.mx = sim->dtemp;
    sim->mpiexch.angle = sim->dpress;
#endif

}

void Inicializer::initTop() {

    long i,j,k;
    char *sysnames[MAXN] = {NULL};
    long  *sysmoln /*[MAXN]*/;
    bool exclusions[MAXT][MAXT] = {false};

    Molecule molecules[MAXMT];

    sysmoln = (long int*) malloc( sizeof(long)*MAXN);
    if(sysmoln == NULL){
        fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sysmoln");
        exit(1);
    }

    readTopoFile(molecules, sysmoln, sysnames, exclusions);
    fprintf (stdout, "\nTopology succesfully read. Generating pair interactions...\n");

    //fill ia_params combinations and topology parameters
    topo->genParamPairs(&exclusions);
    topo->genTopoParams();


    //TODO fill chain list and maxch, park particle type
    fprintf (stdout, "Generating chainlist...\n");

    setParticlesParams(sysmoln, sysnames, molecules);

    // Initialize the clusterlist
    sim->clusterlist = (long int*) malloc(sizeof(long) * conf->particleStore.size());
    if(sim->clusterlist == NULL){
        fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sim->clusterlist!");
        exit(1);
    }
    sim->clustersenergy = (double*) malloc(sizeof(double) * conf->particleStore.size());
    if(sim->clustersenergy== NULL){
        fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for sim->clustersenergy!");
        exit(1);
    }
    sim->clusters = NULL;

    // get all the particles with switch type
    long switchlist[conf->particleStore.size()];
    long n_switch_part = 0;
    for(i = 0; i < (long)conf->particleStore.size(); i++){
        if(conf->particleStore[i].type != conf->particleStore[i].switchtype){
            switchlist[n_switch_part] = i;
            n_switch_part++;
        }
    }
    topo->n_switch_part = n_switch_part;
    if (n_switch_part == 0 && sim->switchprob > 0){
        fprintf(stderr, "TOPOLOGY WARNING: No switchable particles found, but probability for a switch is not zero!\n");
        sim->switchprob = 0;
        fprintf(stderr, "TOPOLOGY WARNING: We changed Switch Probability to zero in this run!\n");
    }
    topo->switchlist=NULL;
    if (n_switch_part > 0){
        topo->switchlist = (long int*) malloc(sizeof(long) * n_switch_part);
        for(i = 0; i < n_switch_part; i++){
            topo->switchlist[i] = switchlist[i];
            //DEBUG
            //printf("%ld is in switchlist\n", switchlist[i]);
        }
    }



    k=0;
    //clear connectivity and then fill it from chain list
    fprintf (stdout, "Generating connectivity...\n");
    for (i=0; i < (long)conf->particleStore.size(); i++) {
        conf->particleStore[i].conlist[0] = -1;
        conf->particleStore[i].conlist[1] = -1;
        conf->particleStore[i].conlist[2] = -1;
        conf->particleStore[i].conlist[3] = -1;
    }
    conf->sysvolume = 0;
    for (i=0; i<(long)conf->particleStore.size(); i++) {
        for (j=0; j<MAXCHL; j++) {
            if (conf->chainlist[i][j] >= 0) {
                k = conf->chainlist[i][j];
                if ((j+1 < MAXCHL)&&(conf->chainlist[i][j+1] >= 0))
                    conf->particleStore[k].conlist[1] = conf->chainlist[i][j+1]; //if there is a next particle fill it to head bond
                if (j > 0)
                    conf->particleStore[k].conlist[0] = conf->chainlist[i][j-1]; //if this is not first particle fill tail bond
                if ((j+2 < MAXCHL)&& (conf->chainlist[i][j+2] >= 0))
                    conf->particleStore[k].conlist[3] = conf->chainlist[i][j+2]; //if there is a second next particle fill it second neighbour
                if (j > 1)
                    conf->particleStore[k].conlist[2] = conf->chainlist[i][j-2]; //if this is not second or first particle fill second tail bond
            }
        }
        conf->sysvolume += topo->ia_params[conf->particleStore[i].type][conf->particleStore[i].type].volume;
    }

    /*DEBUG
      for (i=0; i<MAXN; i++) {
      for (j=0; j<MAXCHL; j++) {
      fprintf (stderr, " %d",chainlist[i][j]);
      }
      fprintf (stderr, " \n");
      }
      for (i=0; i<MAXN; i++) {
      printf (" %ld %ld %ld %ld\n",conlist[i][0],conlist[i][1],conlist[i][2],conlist[i][3]);
      }
     */

    //  Mark particles as not switched
    for(i = 0; i < (long)conf->particleStore.size(); i++){
        conf->particleStore[i].switched = 0;
    }

    topDealoc(sysnames,&sysmoln);
    DEBUG_INIT("Finished with reading the topology");

    // Parallel tempering check
    #ifdef MPI
        // probability to switch replicas = exp ( -0.5 * dT*dT * N / (1 + dT) )
        printf("Probability to switch replicas is roughly: %f\n",exp(-0.5 * maxpart * sim->dtemp * sim->dtemp / (1.0 + sim->dtemp)) );
    #endif

    return;
}

void Inicializer::initConfig() {
    cout << "\nReading configuration...\n";

    int err,fields,tmp_type;
    long i,j,current,first;
    FILE * infile;
    char * line, line2[STRLEN];
    size_t line_size = (STRLEN + 1) * sizeof(char);
    line = (char *) malloc(line_size);
    Particle chorig[MAXCHL];

    double maxlength = 0;
    for(i = 0; i < MAXT; i++){
        if(maxlength < topo->ia_params[i][i].len[0])
            maxlength = topo->ia_params[i][i].len[0];
    }


    infile = fopen(files->configurationInFile, "r");
    if (infile == NULL) {
        fprintf (stderr, "\nERROR: Could not open config.init file.\n\n");
        exit (1);
    }

    if(myGetLine(&line, &line_size, infile) == -1){
        fprintf (stderr, "ERROR: Could not read box size.\n\n");
        exit (1);
    }
    strip_comment(line);
        trim(line);
    if (sscanf(line, "%le %le %le", &(conf->box.x), &(conf->box.y), &(conf->box.z)) != 3) {
        if(myGetLine(&line, &line_size, infile) == -1){
              fprintf (stderr, "ERROR: Could not read box size.\n\n");
              exit (1);
        }
        aftercommand(line2,line,BOXSEP);
        strip_comment(line2);
        trim(line2);
        if (sscanf(line2, "%le %le %le", &(conf->box.x), &(conf->box.y), &(conf->box.z)) != 3) {
              fprintf (stderr, "ERROR: Could not read box size.\n\n");
              exit (1);
        }
    }
    if (conf->box.x < maxlength * 2.0 + 2.0) {
        printf ("WARNING: x box length is less than two spherocylinders long.\n\n");
    }
    if (conf->box.y < maxlength * 2.0 + 2.0) {
        printf ("WARNING: y box length is less than two spherocylinders long.\n\n");
    }
    if (conf->box.z < maxlength * 2.0 + 2.0) {
        printf ("WARNING: z box length is less than two spherocylinders long.\n\n");
    }

    DEBUG_INIT("Position of the particle");
    for (i=0; i < (long)conf->particleStore.size(); i++) {
        if(myGetLine(&line, &line_size, infile) == -1){
            break;
        }
        strip_comment(line);
        trim(line);
        fields = sscanf(line, "%le %le %le %le %le %le %le %le %le %d",
                &conf->particleStore[i].pos.x, &conf->particleStore[i].pos.y, &conf->particleStore[i].pos.z,
                &conf->particleStore[i].dir.x, &conf->particleStore[i].dir.y, &conf->particleStore[i].dir.z,
                &conf->particleStore[i].patchdir[0].x, &conf->particleStore[i].patchdir[0].y, &conf->particleStore[i].patchdir[0].z,
                &conf->particleStore[i].switched);
        conf->particleStore[i].patchdir[1].x = conf->particleStore[i].patchdir[1].y = conf->particleStore[i].patchdir[1].z =0;
        conf->particleStore[i].chdir[0].x = conf->particleStore[i].chdir[0].y = conf->particleStore[i].chdir[0].z =0;
        conf->particleStore[i].chdir[1].x = conf->particleStore[i].chdir[1].y = conf->particleStore[i].chdir[1].z =0;
        DEBUG_INIT("Line: %s\nNumber of Fields: %d", line, fields);
        if (fields == 9){
            conf->particleStore[i].switched = 0;
            fprintf(stdout, "WARNING: Particle %ld is assumed to be not switched!\n", i+1);
            fields++;
        }
        if (fields != 10) {
            fprintf (stderr, "ERROR: Could not read coordinates for particle %ld.\n \
                    Did you specify box size at the begining?\n\n", i+1);
            free(line);
            exit (1);
        }
        /* Scale position vector to the unit cube */
        usePBC(&conf->particleStore[i].pos, conf->box );

        conf->particleStore[i].pos.x /= conf->box.x;
        conf->particleStore[i].pos.y /= conf->box.y;
        conf->particleStore[i].pos.z /= conf->box.z;

        if ((topo->ia_params[conf->particleStore[i].type][conf->particleStore[i].type].geotype[0]<SP)&&( DOT(conf->particleStore[i].dir, conf->particleStore[i].dir) < ZEROTOL )) {
            //DEBUG_INIT("Geotype = %d < %d", conf->particleStore[i].geotype,SP);
            fprintf (stderr,
                    "ERROR: Null direction vector supplied for particle %ld.\n\n", i+1);
            free(line);
            exit (1);
        } else {
            conf->particleStore[i].dir.normalise();
        }

        if ((topo->ia_params[conf->particleStore[i].type][conf->particleStore[i].type].geotype[0]<SP)&&( DOT(conf->particleStore[i].patchdir[0], conf->particleStore[i].patchdir[0]) < ZEROTOL )) {
            fprintf (stderr,
                    "ERROR: Null patch vector supplied for particle %ld.\n\n", i+1);
            free(line);
            exit (1);
        } else {
            ortogonalise(&conf->particleStore[i].patchdir[0],&conf->particleStore[i].dir);
            conf->particleStore[i].patchdir[0].normalise();
        }
        // Switch the type
        if(conf->particleStore[i].switched){
            if(conf->particleStore[i].switchtype == 0){
                fprintf(stderr, "ERROR: Particle %ld switched even though it has no switchtype", i);
                free(line);
                exit(1);
            }
            tmp_type = conf->particleStore[i].type;
            conf->particleStore[i].type = conf->particleStore[i].switchtype;
            conf->particleStore[i].switchtype = tmp_type;
        }

        DEBUG_INIT("%ld:\t%lf\t%lf\t%lf", i, conf->particleStore[i].pos.x, conf->particleStore[i].pos.y, conf->particleStore[i].pos.z);

    }
    free(line);
    /*Make chains WHOLE*/
    for (i=0;i<conf->chainCount;i++){
        j=0;
        current = conf->chainlist[i][0];
        first = current;
        chorig[0].pos = conf->particleStore[first].pos;
        while (current >=0 ) {
            /*shift the chain particle by first one*/
            conf->particleStore[current].pos.x -= chorig[0].pos.x;
            conf->particleStore[current].pos.y -= chorig[0].pos.y;
            conf->particleStore[current].pos.z -= chorig[0].pos.z;
            /*put it in orig box*/
            conf->particleStore[current].pos.x -=  anInt(conf->particleStore[current].pos.x);
            conf->particleStore[current].pos.y -=  anInt(conf->particleStore[current].pos.y);
            conf->particleStore[current].pos.z -=  anInt(conf->particleStore[current].pos.z);
            //printf("ant: %f %f %f\n",conf->particleStore[current].pos.x,conf->particleStore[current].pos.y,conf->particleStore[current].pos.z);
            /*shot it back*/
            conf->particleStore[current].pos.x += chorig[0].pos.x;
            conf->particleStore[current].pos.y += chorig[0].pos.y;
            conf->particleStore[current].pos.z += chorig[0].pos.z;
            //printf("posstart: %f %f %f\n",conf->particleStore[current].pos.x,conf->particleStore[current].pos.y,conf->particleStore[current].pos.z);
            j++;
            current = conf->chainlist[i][j];
        }
    }

    err = 0;
    //for (i=0; i < topo->npart-1; i++) {
    //    for (j=i+1; j < topo->npart; j++) {
    //        if ( overlap(conf->particleStore[i], conf->particle[j], conf->box, topo->ia_params) ) {
    //            fprintf (stderr,
    //                    "ERROR: Overlap in initial coniguration between particles %ld and %ld.\n",
    //                    i+1, j+1);
    //            err = 1;
    //        }
    //    }
    //}
    if (err) {
        printf ("\n");
        exit (1);
    }

    fclose (infile);

    fflush (stdout);
}


void Inicializer::testChains() {
    if (conf->chainCount == 0) {    // no chain -> make the probability of moving them 0
        if (sim->chainprob > 0)
            printf ("No chains... chain move probability set to 0.\n");
        sim->chainprob = 0;
    }
}

void Inicializer::initWriteFiles() {
    cout << "\nPatchy Spherocylinders version 3.6\n";
    cout << "-------------------------------------\n";

    sprintf(files->configurationInFileMuVTChains, "configMuVT");
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

void Inicializer::initPairlist() {
    printf("\nAllocating memory for pairlist...\n");

    //sim->pairlist = (Pairs*) xMalloc(sizeof(Pairs) * topo->npart); // deprecated, solved bz allocating particle

    // Highest guess: Every particle interacts with the others
    // TODO: Make it more sophisticated
    long i;
    /*for(i = 0; i < topo->npart; i++){ // del after
        sim->pairlist[i].pairs = (long int*) malloc(sizeof(long) * topo->npart);
        sim->pairlist[i].num_pairs = 0;
    }*/

    for(i = 0; i < (long)conf->particleStore.size(); i++){
        conf->particleStore[i].neighborID = (long int*) malloc(sizeof(long) * conf->particleStore.size());
        conf->particleStore[i].neighborCount = 0;
    }
}

void Inicializer::initMPI() {
#ifdef MPI
    FILE *infile;
    printf(" MPI version");
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &(sim.mpinprocs) );
    MPI_Comm_rank(MPI_COMM_WORLD, &(sim.mpirank) );

    sprintf(files.configurationoutfile, "%dconfig.last", sim.mpirank);
    sprintf(files.moviefile, "%dmovie", sim.mpirank);
    sprintf(files.wloutfile, "%dwl-new.dat", sim.mpirank);
    sprintf(files.clusterfile, "%dcluster.dat", sim.mpirank);
    sprintf(files.clusterstatfile, "%dcluster_stat.dat", sim.mpirank);
    sprintf(files.energyfile, "%denergy.dat", sim.mpirank);
    sprintf(files.statfile, "%dstat.dat", sim.mpirank);

    /*test if there is a specific input configuration for mpi run*/
    sprintf(files.configurationinfile, "%dconfig.init", sim.mpirank);
    infile = fopen(files.configurationinfile, "r");
    if (infile != NULL)
        fclose (infile);
    else  sprintf(files.configurationinfile, "config.init");

    /*test if there is a specific input wang-landau for mpi run*/
    sprintf(files.wlinfile, "%dwl.dat", sim.mpirank);
    infile = fopen(files.wlinfile, "r");
    if (infile != NULL)
        fclose (infile);
    else  sprintf(files.wlinfile, "wl.dat");
#endif
}





/************************************************************************************************
 *                                      PRIVATE METHODS                                         *
 ************************************************************************************************/



void Inicializer::setParticlesParams(long *sysmoln, char **sysnames, Molecule *molecules) {
    long i=0, j=0, mol, k, maxpart=0, maxch=0;

    try{
        conf->particleStore.reserve(MAXN);
    } catch(std::bad_alloc& bad) {
        fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for conf->particleStore");
        exit(1);
    }

    while (sysnames[i]!=NULL) {
        mol=0;
        while (strcmp(molecules[mol].name,sysnames[i])) {
            mol++;
            if (molecules[mol].name == NULL) {
                fprintf (stderr, "TOPOLOGY ERROR: molecules %s is not defined.\n\n",sysnames[i]);
                topDealoc(sysnames,&sysmoln);
                exit(1);
            }
        }

        for (j=0;j<sysmoln[i];j++) {
            //DEBUG	    fprintf (stdout, "molnames %s sysname %s sysnum %ld \n",molnames[mol],sysnames[i],sysmoln[i]);
            k=0;
            while (molecules[mol].type[k] != -1) {

                conf->particleStore.push_back(Particle());
                conf->particleStore[maxpart].type        = molecules[mol].type[k];
                conf->particleStore[maxpart].switchtype  = molecules[mol].switchtype[k];
                conf->particleStore[maxpart].delta_mu    = molecules[mol].delta_mu[k];
                conf->particleStore[maxpart].molType     = mol;
                conf->particleStore[maxpart].chainIndex  = maxch;

                if (k > MAXCHL) {
                    fprintf (stderr, "TOPOLOGY ERROR: more particles in chan (%ld) than allowed(%d).\n",k,MAXCHL);
                    fprintf (stderr, "Change MAXCHL in source and recompile the program. \n\n");
                    topDealoc(sysnames,&sysmoln);
                    exit(1);
                }
                if (molecules[mol].type[1] != -1) {
                    conf->chainlist[maxch][k] = maxpart;
                }
                k++;
                maxpart++;
                if (maxpart > MAXN) {
                    fprintf (stderr, "TOPOLOGY ERROR: more particles(%ld) than allowed(%d).\n",maxpart,MAXN);
                    fprintf (stderr, "Change MAXN in source and recompile the program. \n\n");
                    topDealoc(sysnames,&sysmoln);
                    exit(1);
                }
            }
            if (molecules[mol].type[1] != -1) {
                maxch++;
            }
        }
        i++;
    }

    int newType = -1;
    for(i = 0; i < maxpart; i++) {
        // set simple grouplist
        if(newType != conf->particleStore[i].molType) {
            newType = conf->particleStore[i].molType;
            conf->first[newType] = i;
            conf->molSize[newType] = topo->chainparam[newType].molSize();
        }
    }
    conf->molTypeCount = newType+1;

    j = 0;
    while (conf->chainlist[j][0] >= 0) {
        j++;
    }
    conf->chainCount = j;

    if (conf->chainCount != maxch) {
        fprintf (stderr, "TOPOLOGY ERROR: Maximum number of chains(%ld) does not agree with number of chains (%ld)\n\n",maxch,conf->chainCount);
        topDealoc(sysnames,&sysmoln);
        exit (1);
    }
}



void Inicializer::readTopoFile(Molecule* molecules, long  *sysmoln, char *sysnames[MAXN], bool exclusions[][MAXT]) {
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

    fprintf (stdout, "Reading topology...\n");
    fflush(stdout);
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
        // skip trailing and leading spaces and comment text
        strip_comment (pline);
        trim (pline);
        // if there is something left...
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
                        topDealoc(sysnames,&sysmoln);
                        free(pline); pline = NULL;
                        exit (1);
                    }
                    DEBUG_INIT("back in init_top");
                } else {
                    if (!strcmp(keystr,"MOLECULES")){
                        DEBUG_INIT("Let's go to the molecules");
                        if (molname[0] == ' ') {
                            beforecommand(molname,pline,SEPARATOR);
                            i=0;
                            while (molecules[i].name != NULL)
                                i++;
                            DEBUG_INIT("in the middle of getting to fillmol");
                            molecules[i].name = (char*) malloc(strlen(molname)+1);
                            strcpy(molecules[i].name, molname);
                            fprintf (stdout, "Topology read for molecule: %s \n",molname);
                        }
                        if (!fillMol(molname, pline, molecules)) {
                            fprintf (stderr, "\nTOPOLOGY ERROR: in reading molecules\n\n");
                            topDealoc(sysnames,&sysmoln);
                            free(pline); pline = NULL;
                            exit (1);
                        }
                        if ((dummy = strchr (pline,CLOSEMOL)) != NULL)
                            molname[0] = ' ';
                    } else {
                        if (!strcmp(keystr,"SYSTEM")) {
                            if (!fillSystem(pline,sysnames,&sysmoln)) {
                                fprintf (stderr, "\nTOPOLOGY ERROR: in reading system\n\n");
                                topDealoc(sysnames,&sysmoln);
                                free(pline); pline = NULL;
                                exit (1);
                            }
                        } else {
                            if (!strcmp(keystr, "RESERVOIR")) {
                                if (!fillSystem(pline,sysnames,&sysmoln)) {
                                    fprintf (stderr, "\nTOPOLOGY ERROR: in reading system\n\n");
                                    topDealoc(sysnames,&sysmoln);
                                    free(pline); pline = NULL;
                                    exit (1);
                                }
                            }
                             else {
                                if (!strcmp(keystr,"EXTER")) {
                                    fflush(stdout);
                                    if (!fillExter(&pline)) {
                                        DEBUG_INIT("Something went wrong with external potential");
                                        fprintf (stderr, "\nTOPOLOGY ERROR: in reading external potential\n\n");
                                        topDealoc(sysnames,&sysmoln);
                                        free(pline); pline = NULL;
                                        exit (1);
                                    }
                                } else {
                                    if (!strcmp(keystr,"EXCLUDE")) {
                                        fflush(stdout);
                                        if (!fillExclusions(&pline,exclusions)) {
                                            DEBUG_INIT("Something went wrong with exclusions potential");
                                            fprintf (stderr, "\nTOPOLOGY ERROR: in reading exclusions\n\n");
                                            topDealoc(sysnames,&sysmoln);
                                            free(pline); pline = NULL;
                                            exit (1);
                                        }
                                    } else {
                                        fprintf (stderr, "\nTOPOLOGY ERROR: invalid keyword:%s.\n\n", keystr);
                                        topDealoc(sysnames,&sysmoln);
                                        free(pline); pline = NULL;
                                        exit (1);
                                    } // else exclude
                                } // else exter
                            } // else system
                        } // else reservoir
                    } // else molecules
                } // else types
            } // else get commnad
        } // if //getcommand
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




int Inicializer::topDealoc(char *sysnames[], long **sysmoln) {

    if ((*sysmoln) != NULL) free((*sysmoln));
        (*sysmoln)=NULL;

    for (int i=0;i<MAXN;i++) {
        if ((sysnames[i]) != NULL) free(sysnames[i]);
            sysnames[i]=NULL;
    }
    return 0;

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


int Inicializer::fillSystem(char *pline, char *sysnames[], long **sysmoln) {
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
    if ((*sysmoln)[i] < 1) {
        fprintf (stderr, "TOPOLOGY ERROR: cannot have %ld number of molecules.\n\n", (*sysmoln)[i]);
        return 0;
    }
    fprintf (stdout, "system: %s %ld\n",sysnames[i],(*sysmoln)[i]);
    return 1;
}


void Inicializer::usePBC(Vector *pos, Vector pbc){
    do {
        (*pos).x += pbc.x;
    } while ((*pos).x < 0.0);
    do {
        (*pos).x -= pbc.x;
    } while ((*pos).x > pbc.x);

    do {
        (*pos).y += pbc.y;
    } while ((*pos).y < 0.0);
    do {
        (*pos).y -= pbc.y;
    } while ((*pos).y > pbc.y);

    do {
        (*pos).z += pbc.z;
    } while ((*pos).z < 0.0);
    do {
        (*pos).z -= pbc.z;
    } while ((*pos).z > pbc.z);
}


int Inicializer::fillTypes(char **pline) {
    int type;
    int geotype_i;
    int fields;
    char name[SMSTR];
    char geotype[SMSTR];

    double param[11];
    /* 0: epsilon
     * 1: sigma
     * 2: attraction dist
     * 3: sttraction switch
     * 4: patch angle
     * 5: patch switch
     * 6: length
     * 7(optional): second patche rotation
     * 8(optional): second patch angle
     * 9(optional): second patch angle switch
     * +1: chirality
     */
    char typestr[STRLEN], paramstr[STRLEN];

    beforecommand(typestr, *pline, SEPARATOR);
    aftercommand(paramstr, *pline, SEPARATOR);

    fields = sscanf(paramstr, "%s %d %s %le %le %le %le %le %le %le %le %le %le %le", name, &type, geotype, &param[0], &param[1], &param[2], &param[3], &param[4], &param[5], &param[6], &param[7], &param[8], &param[9], &param[10]);
    fields -= 5; // number of parameter fields => I am too lazy to adjust everywhere below the numbers
    //DEBUG    fprintf (stdout, "Topology read geotype: %ld with parameters fields %d, str:%s and %s in pline %s\n",geotype,fields,geotypestr,paramstr,pline);

    geotype_i = convertGeotype(geotype);
    if(!geotype_i){
        fprintf(stderr, "TOPOLOGY ERROR: Unknown GEOTYPE: %s!", geotype);
        return 0;
    }
    DEBUG_INIT("geotype_i: %d; fields = %d", geotype_i, fields);
    if (( (geotype_i == SCN) || (geotype_i == SPN) ) && (fields != 0)) {
        fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 1.\n\n", geotype);
        return 0;
    }
    if (( (geotype_i == SCA) || (geotype_i == SPA)) && (fields != 2)) {
        fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 3.\n\n", geotype);
        return 0;
    }
    if (( (geotype_i == PSC) || (geotype_i == CPSC) ) && (fields != 5)) {
        fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 5.\n\n", geotype);
        return 0;
    }
    if (( (geotype_i == CHCPSC) || (geotype_i == CHCPSC) )&& ( fields != 6)) {
        fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 6.\n\n", geotype);
        return 0;
    }
    if (( (geotype_i == TPSC) || (geotype_i == TCPSC) ) && (fields != 8)) {
        fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 8, is %d.\n\n", geotype, fields);
        return 0;
    }
    if (( (geotype_i == TCHCPSC) || (geotype_i == TCHCPSC) )&& ( fields != 9)) {
        fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for %s geotype, should be 9.\n\n", geotype);
        return 0;
    }

    if ((geotype_i < 0) || (geotype_i > (MAXT + 10))) {
        fprintf (stderr, "TOPOLOGY ERROR: geotype (%s) is out of range: 0 - %d.\n\n", geotype, MAXT + 10);
        return 0;
    }

    strcpy(topo->ia_params[type][type].name, name);
    strcpy(topo->ia_params[type][type].other_name, name);
    topo->ia_params[type][type].geotype[0] = geotype_i;
    topo->ia_params[type][type].geotype[1] = geotype_i;
    topo->ia_params[type][type].epsilon = param[0];
    topo->ia_params[type][type].sigma = param[1];
    topo->ia_params[type][type].rcutwca = (topo->ia_params[type][type].sigma)*pow(2.0,1.0/6.0);
    fprintf(stdout, "Topology read of %d: %s (geotype: %s, %d) with parameters %f %f", type, name, geotype, geotype_i, topo->ia_params[type][type].epsilon, topo->ia_params[type][type].sigma);
    if (fields > 0) {
        topo->ia_params[type][type].pdis = param[2];
        topo->ia_params[type][type].pswitch = param[3];
        topo->ia_params[type][type].rcut = topo->ia_params[type][type].pswitch+topo->ia_params[type][type].pdis;
        fprintf(stdout, " %f %f",topo->ia_params[type][type].pdis,topo->ia_params[type][type].pswitch);
    }
    if (fields > 2) {
        int i;
        for(i = 0; i < 2; i++){
            topo->ia_params[type][type].len[i] = param[6];
            topo->ia_params[type][type].half_len[i] = param[6] / 2;
            topo->ia_params[type][type].pangl[i] = param[4];
            topo->ia_params[type][type].panglsw[i] = param[5];
            topo->ia_params[type][type].pcangl[i] = cos(param[4]/2.0/180*PI);                 // C1
            topo->ia_params[type][type].pcanglsw[i] = cos((param[4]/2.0+param[5])/180*PI);    // C2
            //topo->ia_params[type][type].pcangl[i] = topo->ia_params[type][type].pcangl[i];
            //topo->ia_params[type][type].pcanglsw[i] =	topo->ia_params[type][type].pcanglsw[i];
            topo->ia_params[type][type].pcoshalfi[i] = cos((param[4]/2.0+param[5])/2.0/180*PI);
            topo->ia_params[type][type].psinhalfi[i] = sqrt(1.0 - topo->ia_params[type][type].pcoshalfi[i] * topo->ia_params[type][type].pcoshalfi[i]);
        }
        fprintf(stdout, " %f %f", topo->ia_params[type][type].pangl[0], topo->ia_params[type][type].panglsw[0]);
    }
    if(fields == 6){
        int i;
        for(i = 0; i < 2; i++){
            topo->ia_params[type][type].chiral_cos[i] = cos(param[7] / 360 * PI);
            topo->ia_params[type][type].chiral_sin[i] = sqrt(1 - topo->ia_params[type][type].chiral_cos[i] * topo->ia_params[type][type].chiral_cos[i]);
            fprintf(stdout, " %f ", param[7]);
        }
    }
    if ((fields == 8)||(fields == 9)) {
        int i;
        for(i = 0; i < 2; i++){
            topo->ia_params[type][type].csecpatchrot[i] = cos(param[7] / 360 * PI);
            topo->ia_params[type][type].ssecpatchrot[i] = sqrt(1 - topo->ia_params[type][type].csecpatchrot[i] * topo->ia_params[type][type].csecpatchrot[i]);
            //fprintf(stdout, " %f %f", topo->ia_params[type][type].csecpatchrot[0], topo->ia_params[type][type].ssecpatchrot[0]);

            topo->ia_params[type][type].pangl[i+2] = param[8];
            topo->ia_params[type][type].panglsw[i+2] = param[9];
            topo->ia_params[type][type].pcangl[i+2] = cos(param[8]/2.0/180*PI);                 // C1
            topo->ia_params[type][type].pcanglsw[i+2] = cos((param[8]/2.0+param[9])/180*PI);    // C2
            //topo->ia_params[type][type].pcangl[i] = topo->ia_params[type][type].pcangl[i];
            //topo->ia_params[type][type].pcanglsw[i] = topo->ia_params[type][type].pcanglsw[i];
            topo->ia_params[type][type].pcoshalfi[i+2] = cos((param[8]/2.0+param[9])/2.0/180*PI);
            topo->ia_params[type][type].psinhalfi[i+2] = sqrt(1.0 - topo->ia_params[type][type].pcoshalfi[i+2] * topo->ia_params[type][type].pcoshalfi[i+2]);
        }
        fprintf(stdout, " %f  %f %f", param[7], topo->ia_params[type][type].pangl[2], topo->ia_params[type][type].panglsw[2]);
    }
    if(fields == 9){
        int i;
        for(i = 0; i < 2; i++){
            topo->ia_params[type][type].chiral_cos[i] = cos(param[10] / 360 * PI);
            topo->ia_params[type][type].chiral_sin[i] = sqrt(1 - topo->ia_params[type][type].chiral_cos[i] * topo->ia_params[type][type].chiral_cos[i]);
            fprintf(stdout, " %f ", param[9]);
        }
    }

    // Volume
    if (geotype_i < SP)
        topo->ia_params[type][type].volume = 4.0/3.0*PI*pow((topo->ia_params[type][type].sigma)/2.0,3.0) + PI/2.0*topo->ia_params[type][type].len[0]*pow((topo->ia_params[type][type].sigma)/2.0,2.0) ;
    else
        topo->ia_params[type][type].volume = 4.0/3.0*PI*pow((topo->ia_params[type][type].sigma)/2.0,3.0);
    if ( topo->ia_params[type][type].rcutwca > topo->sqmaxcut )
        topo->sqmaxcut = topo->ia_params[type][type].rcutwca;
    if ( topo->ia_params[type][type].rcut > topo->sqmaxcut )
        topo->sqmaxcut = topo->ia_params[type][type].rcut;
    fprintf(stdout, " \n");
    DEBUG_INIT("Finished filltypes");
    return 1;
}


int Inicializer::convertGeotype(char *geotype) {
    if (strcmp(geotype, "CPSC") == 0)
        return CPSC;
    if (strcmp(geotype, "CHCPSC") == 0)
        return CHCPSC;
    if (strcmp(geotype, "SCA") == 0)
        return SCA;
    if (strcmp(geotype, "PSC") == 0)
        return PSC;
    if (strcmp(geotype, "CHPSC") == 0)
        return CHPSC;
    if (strcmp(geotype, "TCPSC") == 0)
        return TCPSC;
    if (strcmp(geotype, "TCHCPSC") == 0)
        return TCHCPSC;
    if (strcmp(geotype, "TPSC") == 0)
        return TPSC;
    if (strcmp(geotype, "TCHPSC") == 0)
        return TCHPSC;
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
        topo->exter.exist = true;
        topo->exter.thickness = param[0];
        fprintf(stdout, "External potential with thickness: %le ",topo->exter.thickness);
        if (fields >1) {
            topo->exter.epsilon = param[1];
            fprintf(stdout, "epsilon: %le ",topo->exter.epsilon);
            if (fields >2) {
                topo->exter.attraction = param[2];
                fprintf(stdout, "and range of attraction: %le ",topo->exter.attraction);
            }
        }
    } else{
        topo->exter.exist = false;
        fprintf(stdout, "No external potential ");
    }

    fprintf(stdout, " \n");
    DEBUG_INIT("Finished filling external potential");
    return 1;
}


int Inicializer::fillMol(char *molname, char *pline, Molecule *molecules) {
    DEBUG_INIT("fillmol just has been called!");
    char str[STRLEN],str2[STRLEN],molcommand[STRLEN],molparams[STRLEN];
    int i,j,fields;
    double bondk,bonddist, mu, lnLambda;
    int muVTmove;

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
            topo->chainparam[i].name = (char*) malloc(strlen(molname)+1);
            strcpy(topo->chainparam[i].name, molname);
        }

        topo->chainparam[i].particleTypes[j] = molecules[i].type[j];

        if (fields == 1){
                (molecules[i].switchtype[j]) = (molecules[i].type[j]);
                (molecules[i].delta_mu[j]) = 0;
                fields = 3;
        } else{
            fprintf(stdout, "(with switchtype: %ld and delta_mu: %lf)", molecules[i].switchtype[j], molecules[i].delta_mu[j]);
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
        topo->chainparam[i].bond1c = bondk;
        topo->chainparam[i].bond1eq = bonddist;
        fprintf (stdout, "bond1: %f %f \n",topo->chainparam[i].bond1c,topo->chainparam[i].bond1eq);
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
        topo->chainparam[i].bond2c = bondk;
        topo->chainparam[i].bond2eq = bonddist;
        fprintf (stdout, "bond2: %f %f \n",topo->chainparam[i].bond2c,topo->chainparam[i].bond2eq);
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
        topo->chainparam[i].bonddc = bondk;
        topo->chainparam[i].bonddeq = bonddist;
        fprintf (stdout, "bondd: %f %f \n",topo->chainparam[i].bonddc,topo->chainparam[i].bonddeq);
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
        topo->chainparam[i].angle1c = bondk;
        topo->chainparam[i].angle1eq = bonddist/180.0*PI;
        fprintf (stdout, "angle1: %f %f \n",topo->chainparam[i].angle1c,topo->chainparam[i].angle1eq);
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
        topo->chainparam[i].angle2c = bondk;
        topo->chainparam[i].angle2eq = bonddist/180.0*PI;
        fprintf (stdout, "angle2: %f %f \n",topo->chainparam[i].angle2c,topo->chainparam[i].angle2eq);
        return 1;
    }

        cout << "3" << endl;

    /// INIT of muVT ensemble
    if (!strcmp(molcommand,"MU")) {
        fields = sscanf(molparams, "%le ", &mu);
        topo->chainparam[i].mu = mu;
        fprintf (stdout, "mu: %f \n",topo->chainparam[i].mu);
        return 1;
    }
    if (!strcmp(molcommand,"MUVTMOVE")) {
        fields = sscanf(molparams, "%d ", &muVTmove);
        topo->chainparam[i].muVTmove = muVTmove;
        cout << std::boolalpha << "muVTmove: " << topo->chainparam[i].muVTmove << endl;
        return 1;
    }
    if (!strcmp(molcommand,"LNLAMBDA")) {
        fields = sscanf(molparams, "%le ", &lnLambda);
        topo->chainparam[i].lnThermalWavelengh = lnLambda;
        fprintf (stdout, "ln thermal wavelenght: %f \n",topo->chainparam[i].lnThermalWavelengh);
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
