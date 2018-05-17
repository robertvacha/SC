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


bool Inicializer::initConfig(FILE** infile, std::vector<Particle > &pvec, bool scale_to_box) {

    int err,fields,tmp_type, check=0;
    long j,current;
    char * line, line2[STRLEN];
    size_t line_size = (STRLEN + 1) * sizeof(char);
    line = (char *) malloc(line_size);
    //Particle chorig[MAXCHL];

    double maxlength = 0.0;
    for(int i = 0; i < MAXT; i++){
        if(maxlength < topo.ia_params[i][i].len[0])
            maxlength = topo.ia_params[i][i].len[0];
    }

    if(myGetLine(&line, &line_size, *infile) == -1){
        fprintf (stderr, "ERROR: Could not read box size (Inicializer::initConfig)\n\n");
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
            cerr << "ERROR: Could not read box size2." << endl;
            return false;
        }
        aftercommand(line2,line,BOXSEP);
        strip_comment(line2);
        trim(line2);
        if (sscanf(line2, "%le %le %le", &(box.x), &(box.y), &(box.z) ) != 3) {
            cerr << "ERROR: Could not read box size3." << endl;
            return false;
        }
    }

    conf->geo = Cuboid(box);
#endif
    if (conf->geo.box.x < maxlength * 2.0 + 2.0) {
        mcout.get() << "WARNING: x (" << conf->geo.box.x << ") geo.box length is less than two spherocylinders long (" << maxlength * 2.0 + 2.0 << ")." << endl;
        //exit(1);
    }
    if (conf->geo.box.y < maxlength * 2.0 + 2.0) {
        mcout.get() << "WARNING: y (" << conf->geo.box.y << ") geo.box length is less than two spherocylinders long (" << maxlength * 2.0 + 2.0 << ")." << endl;
        //exit(1);
    }
    if (conf->geo.box.z < maxlength * 2.0 + 2.0) {
        mcout.get() <<"WARNING: z (" << conf->geo.box.z << ") geo.box length is less than two spherocylinders long (" << maxlength * 2.0 + 2.0 << ")." << endl;
        //exit(1);
    }

    for(unsigned int i=0; i < pvec.size(); i++) {
        if(myGetLine(&line, &line_size, *infile) == -1)
            break;

        ++check;

        strip_comment(line);
        trim(line);

        fields = sscanf(line, "%le %le %le %le %le %le %le %le %le %d",
                        &pvec[i].pos.x, &pvec[i].pos.y, &pvec[i].pos.z,
                        &pvec[i].dir.x, &pvec[i].dir.y, &pvec[i].dir.z,
                        &pvec[i].patchdir[0].x, &pvec[i].patchdir[0].y, &pvec[i].patchdir[0].z,
                        &pvec[i].switched);

        if(fields < 9){
            cerr << "Not enough particles in " << files->configurationInFile << endl;
            exit(1);
        }

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

        // For analysis of sheet
        //conf->geo.usePBC2(&pvec[i]); // range 0 - box

        if(scale_to_box) {
            pvec[i].pos.x /= conf->geo.box.x;
            pvec[i].pos.y /= conf->geo.box.y;
            pvec[i].pos.z /= conf->geo.box.z;

            // for compatibility unfortunately
            conf->geo.usePBC(&pvec[i]);
        }
#endif

        if ((topo.ia_params[pvec[i].type][pvec[i].type].geotype[0]<SP)&&( DOT(pvec[i].dir, pvec[i].dir) < ZEROTOL )) {
            //DEBUG_INIT("Geotype = %d < %d", conf->pvec[i].geotype,SP);
            fprintf (stderr, "ERROR: Null direction vector supplied for particle %u.\n\n", i+1);
            free(line);
            return false;
        } else {
            pvec[i].dir.normalise();
        }

        if ((topo.ia_params[pvec[i].type][pvec[i].type].geotype[0]<SP && topo.ia_params[pvec[i].type][pvec[i].type].geotype[0] != SCN )&&( DOT(pvec[i].patchdir[0], pvec[i].patchdir[0]) < ZEROTOL )) {
            fprintf (stderr, "ERROR: Null patch vector supplied for particle %u.\n\n", i+1);
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

    if(check+1 < conf->pvec.size()) {
        cerr << "Not enough particles in " << files->configurationInFile << endl;
        exit(1);
    }
    free(line);
    /*Make chains WHOLE*/
//    for (int i=0; i<conf->pvec.getChainCount(); i++){
//        j=0;
//        current = conf->pvec.getChainPart(i,0);
//        first = current;
//        chorig[0].pos = pvec[first].pos;
//        while (current >=0 ) {
//            /*shift the chain particle by first one*/
//            pvec[current].pos.x -= chorig[0].pos.x;
//            pvec[current].pos.y -= chorig[0].pos.y;
//            pvec[current].pos.z -= chorig[0].pos.z;
//            /*put it in orig geo.box*/
//            pvec[current].pos.x -=  anInt(pvec[current].pos.x);
//            pvec[current].pos.y -=  anInt(pvec[current].pos.y);
//            pvec[current].pos.z -=  anInt(pvec[current].pos.z);
//            //printf("ant: %f %f %f\n",conf->pvec[current].pos.x,conf->pvec[current].pos.y,conf->pvec[current].pos.z);
//            /*shot it back*/
//            pvec[current].pos.x += chorig[0].pos.x;
//            pvec[current].pos.y += chorig[0].pos.y;
//            pvec[current].pos.z += chorig[0].pos.z;
//            //printf("posstart: %f %f %f\n",conf->pvec[current].pos.x,conf->pvec[current].pos.y,conf->pvec[current].pos.z);
//            j++;
//            current = conf->pvec.getChainPart(i,j);
//        }
//    }

        for (int i=0; i<conf->pvec.getChainCount(); i++){
            j=0;
            Molecule mol;
            current = conf->pvec.getChainPart(i,0);
            while (current >=0 ) {
                mol.push_back(current);
                j++;
                current = conf->pvec.getChainPart(i,j);
            }
            conf->makeMoleculeWhole(&mol);
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
    if (conf->pvec.getChainCount() == 0 && !topo.gcSpecies) {    // no chain -> make the probability of moving them 0
        if (sim->chainprob > 0)
            mcout.get() << "No chains... chain move probability set to 0." << endl;
        sim->chainprob = 0.0;
    } else {
        for(int i=0; i<conf->pvec.molTypeCount; i++) {
            if(topo.moleculeParam[i].isGrandCanonical() && !topo.moleculeParam[i].isAtomic()) {
                if(!topo.poolConfig) {
                    mcout.get() << "ChainInsert with no Pool system stated! State [Pool] in top.init" << endl;
                    exit(1);
                }
            }
        }
    }
}

void Inicializer::initNeighborList() {
    mcout.get() << "\nAllocating memory for pairlist, " << (double) conf->neighborList.size() * sizeof(long) * MAXNEIGHBORS / (1024 *1024) << " MB" << endl;

    // Highest guess: Every particle interacts with the others
    // TODO: Make it more sophisticated
    conf->neighborList.resize(conf->pvec.size());
    for(unsigned long i = 0; i < conf->neighborList.size(); i++){
        conf->neighborList[i].neighborID = (long int*) malloc(sizeof(long) * MAXNEIGHBORS);
        conf->neighborList[i].neighborCount = 0;
    }
}

void Inicializer::initGroupLists() {

    mcout.get() << "Generating GroupLists..." << endl;

    // setGroupList;
    conf->pvec.molTypeCount = 0;
    while(!topo.moleculeParam[conf->pvec.molTypeCount].name.empty()) { // get all types
        //cout << topo.moleculeParam[conf->pvec.molTypeCount].name << " " << topo.moleculeParam[conf->pvec.molTypeCount].molType << endl;
        conf->pvec.molTypeCount++;
    }

    // Set first of each type
    int i=0;
    bool empty = true;
    int count=0;
    conf->pvec.first[0] = 0;
    for(int type=0; type < conf->pvec.molTypeCount; type++) {
        empty = true;
        count = 0;
        while(i < (int)conf->pvec.size()) {
            if(type == conf->pvec[i].molType) { // note: we arent searching for molType of particle, could be 0 particles
                conf->pvec.first[type] = i;

                // FIX all others empty after this one
                for(unsigned int j=i; j<conf->pvec.size(); j++) {
                    if(type == conf->pvec[j].molType)
                        count++;
                }
                for(int j = type+1; j < conf->pvec.molTypeCount; j++) {
                    conf->pvec.first[j] = i+count;
                }

                empty = false;
                break;
            }
            i++;
        }
        if(empty) {
            i=conf->pvec.first[type];
        }
    }

    conf->pvec.first[conf->pvec.molTypeCount] = conf->pvec.size();

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



void Inicializer::setParticlesParamss(std::vector<Particle> *pvec, System& system) {
    long allpart=0;

    for(unsigned int i = 0; i < system.moleculeName.size(); ++i) { // loop over defined types in [System]

        for(int perMol=0; perMol < system.moleculeCount[i] * topo.moleculeParam[i].particleTypes.size(); ++perMol) { // number of molecules of given type i
            pvec->push_back(Particle());

            (*pvec)[allpart].type        = topo.moleculeParam[i].particleTypes[ perMol % topo.moleculeParam[i].particleTypes.size() ];
            (*pvec)[allpart].switchtype  = (topo.moleculeParam[i].switchTypes.empty()) ? -1 :
                                           topo.moleculeParam[i].switchTypes[perMol % topo.moleculeParam[i].switchTypes.size()];
            (*pvec)[allpart].delta_mu    = (topo.moleculeParam[i].deltaMu.empty()) ? 0.0 :
                                           topo.moleculeParam[i].deltaMu[perMol % topo.moleculeParam[i].deltaMu.size()];
            (*pvec)[allpart].molType     = i;

            ++allpart;

            if (allpart > MAXN) {
                fprintf (stderr, "TOPOLOGY ERROR: more particles(%ld) than allowed(%d).\n",allpart,MAXN);
                fprintf (stderr, "Change MAXN in source and recompile the program. \n\n");
                exit(1);
            }
        }
    }

    //  Mark particles as not switched
    for(unsigned int i = 0; i < pvec->size(); ++i){
        (*pvec)[i].switched = 0;
    }

    assert(allpart == (long)pvec->size());
}




void *Inicializer::xMalloc(size_t num) {
    void *neww = malloc (num);
    if (!neww){
        fprintf(stderr, "Couldn't allocate any memory!\n");
        exit(1);
    }
    return neww;
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
