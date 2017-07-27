#include "topo.h"

#include <cstring>

void Topo::genParamPairs(bool exclusions[MAXT][MAXT]) {
    int a[2];
    int len;
    double length = 0; // The length of a PSC, currently only one is allow, ie implemented

    for (int i=0;i<MAXT;i++) {
        for (int j=0;j<MAXT;j++) {
            if (i!=j) {
                if((ia_params[j][j].geotype[0] != 0) && (ia_params[i][i].geotype[0] != 0)) {
                    a[0] = i;
                    a[1] = j;

                    for(int k = 0; k < 2; k++){
                        ia_params[i][j].geotype[k] = ia_params[a[k]][a[k]].geotype[0];
                        ia_params[i][j].len[k]     = ia_params[a[k]][a[k]].len[0];
                        if (ia_params[a[k]][a[k]].len[0] > 0){

                            if (length == 0){
                                length = ia_params[a[k]][a[k]].len[0];
                            } else {
                                if (length > 0) {
                                    if (length != ia_params[a[k]][a[k]].len[0]) {
                                        fprintf(stderr, "Error: ");
                                        fprintf(stderr, "Different lengths for spherocylinders have not been implemented yet!\n");
                                        fprintf(stderr, "\tCheck the length of type %d!\n", a[k]);
                                        exit(1);
                                    }
                                }
                            }
                        }
                        ia_params[i][j].half_len[k] = ia_params[a[k]][a[k]].half_len[0];
                        /* Handle angles only, when geotype is a patchs sphero cylinder */
                        if(ia_params[i][j].geotype[k] >= PSC && ia_params[i][j].geotype[k] < SP) {

                            ia_params[i][j].pangl[k]      = ia_params[a[k]][a[k]].pangl[0];
                            ia_params[i][j].panglsw[k]    = ia_params[a[k]][a[k]].panglsw[0];
                            ia_params[i][j].pcangl[k]     = cos(ia_params[i][j].pangl[k]/2.0/180*PI);
                            ia_params[i][j].pcanglsw[k]   = cos((ia_params[i][j].pangl[k]/2.0+ia_params[i][j].panglsw[k])/180*PI);
                            ia_params[i][j].pcoshalfi[k]  = cos((ia_params[i][j].pangl[k]/2.0+ia_params[i][j].panglsw[k])/2.0/180*PI);
                            ia_params[i][j].psinhalfi[k]  = sqrt(1.0 - ia_params[i][j].pcoshalfi[k] * ia_params[i][j].pcoshalfi[k]);

                        }

                        /* Only when the PSC is chiral */
                        if( (ia_params[i][j].geotype[k] == CHCPSC) || (ia_params[i][j].geotype[k] == CHPSC) \
                           || (ia_params[i][j].geotype[k] == TCHCPSC) || (ia_params[i][j].geotype[k] == TCHPSC) ) {
                            ia_params[i][j].chiral_cos[k] = ia_params[a[k]][a[k]].chiral_cos[0];
                            ia_params[i][j].chiral_sin[k] = ia_params[a[k]][a[k]].chiral_sin[0];
                        }
                        /* Information of two patches */
                        if( (ia_params[i][j].geotype[k] == TCPSC) ||
                            (ia_params[i][j].geotype[k] == TPSC) ||
                            (ia_params[i][j].geotype[k] == TCHCPSC) ||
                            (ia_params[i][j].geotype[k] == TCHPSC) ){

                            ia_params[i][j].csecpatchrot[k] = ia_params[a[k]][a[k]].csecpatchrot[0];
                            ia_params[i][j].ssecpatchrot[k] = ia_params[a[k]][a[k]].ssecpatchrot[0];

                            ia_params[i][j].pangl[k+2]     = ia_params[a[k]][a[k]].pangl[2];
                            ia_params[i][j].panglsw[k+2]   = ia_params[a[k]][a[k]].panglsw[2];
                            ia_params[i][j].pcangl[k+2]    = cos(ia_params[i][j].pangl[k+2]/2.0/180*PI);
                            ia_params[i][j].pcanglsw[k+2]  = cos((ia_params[i][j].pangl[k+2]/2.0+ia_params[i][j].panglsw[k+2])/180*PI);
                            ia_params[i][j].pcoshalfi[k+2] = cos((ia_params[i][j].pangl[k+2]/2.0+ia_params[i][j].panglsw[k+2])/2.0/180*PI);
                            ia_params[i][j].psinhalfi[k+2] = sqrt(1.0 - ia_params[i][j].pcoshalfi[k+2] * ia_params[i][j].pcoshalfi[k+2]);
                        }
                    }
                    len = strlen(ia_params[i][i].name);
                    strncpy(ia_params[i][j].name, ia_params[i][i].name, len + 1);
                    len = strlen(ia_params[i][i].other_name);
                    strncpy(ia_params[i][j].other_name, ia_params[i][i].other_name, len + 1);

                    ia_params[i][j].sigma   = AVER(ia_params[i][i].sigma,ia_params[j][j].sigma);
                    ia_params[i][j].sigmaSq   = ia_params[i][j].sigma * ia_params[i][j].sigma;
                    ia_params[i][j].epsilon = sqrt(ia_params[i][i].epsilon *  ia_params[j][j].epsilon);
                    ia_params[i][j].A = 4 * ia_params[i][j].epsilon * pow(ia_params[i][j].sigma, 12 );
                    ia_params[i][j].B = 4 * ia_params[i][j].epsilon * pow(ia_params[i][j].sigma, 6 );
                    ia_params[i][j].pswitch = AVER(ia_params[i][i].pswitch,ia_params[j][j].pswitch);
                    ia_params[i][j].pswitchINV = 1.0 / ia_params[i][j].pswitch;
                    ia_params[i][j].rcutwca = (ia_params[i][j].sigma)*pow(2.0,1.0/6.0);
                    ia_params[i][j].rcutwcaSq = ia_params[i][j].rcutwca * ia_params[i][j].rcutwca;
		     if ((ia_params[i][i].parallel > 0) && (ia_params[j][j].parallel > 0)) ia_params[i][j].parallel = sqrt(ia_params[i][i].parallel *  ia_params[j][j].parallel);
		     if ((ia_params[i][i].parallel < 0) && (ia_params[j][j].parallel < 0)) ia_params[i][j].parallel = -sqrt(ia_params[i][i].parallel *  ia_params[j][j].parallel);
		     
		    
		    
                    // Averaging of the flat part of attraction
                    ia_params[i][j].pdis = AVER(ia_params[i][i].pdis - ia_params[i][i].rcutwca,
                                                      ia_params[j][j].pdis - ia_params[j][j].rcutwca) + ia_params[i][j].rcutwca;
                    ia_params[i][j].pdisSq = ia_params[i][j].pdis * ia_params[i][j].pdis;
                    ia_params[i][j].rcut = ia_params[i][j].pswitch+ia_params[i][j].pdis;
                    ia_params[i][j].rcutSq = ia_params[i][j].rcut * ia_params[i][j].rcut;

                    // if not non-attractive == if attractive
                    if (!((ia_params[i][j].geotype[0] % 10 == 0) || (ia_params[i][j].geotype[1] % 10 == 0))) {

                        if (ia_params[i][j].rcutwca >  ia_params[i][j].rcut) {
                            fprintf(stderr, "Error: Repulsive cutoff is larger than the attractive cutoff!\n");
                            fprintf(stderr, "       between %d and %d: %f > %f\n", i, j, ia_params[i][j].rcutwca, ia_params[i][j].rcut);
                        }
                    }

                    if ( ia_params[i][j].rcutwca > sqmaxcut )
                        sqmaxcut = ia_params[i][j].rcutwca;
                    if ( ia_params[i][j].rcut > sqmaxcut )
                        sqmaxcut = ia_params[i][j].rcut;
                } // (ia_params[j][j].geotype[0] != 0) && (ia_params[i][i].geotype[0] != 0)
            }

        } // end of for cycle j

        /*filling interaction with external potential*/
        if( (exter.exist) && (ia_params[i][i].geotype[0] != 0)){
            /*use everything like for given particles except distance and attraction, which is generated as for other interactions*/
            exter.interactions[i] = ia_params[i][i];
            exter.interactions[i].sigma = AVER(ia_params[i][i].sigma, exter.thickness);
            exter.interactions[i].rcutwca = (exter.interactions[i].sigma)*pow(2.0,1.0/6.0);
            exter.interactions[i].epsilon = sqrt(ia_params[i][i].epsilon *  exter.epsilon);
            exter.interactions[i].pswitch = AVER(ia_params[i][i].pswitch, exter.attraction);
            exter.interactions[i].pdis = AVER(ia_params[i][i].pdis - ia_params[i][i].rcutwca, 0.0) + exter.interactions[i].rcutwca;
            exter.interactions[i].rcut = exter.interactions[i].pswitch + exter.interactions[i].pdis;
            if (exter.interactions[i].rcut > exter.sqmaxcut ) exter.sqmaxcut = exter.interactions[i].rcut;
        }
    }
    for (int i=0;i<MAXT;i++) {
        for (int j=0;j<MAXT;j++) {
            ia_params[i][j].exclude = exclusions[i][j];
        }
    }

    cout << "\n\n\n check pdis and rcutWCA validity - NOT DONE\n\n" << endl;
}

void Topo::genTopoParams() {
    double maxlength = 0;
    for(int i = 0; i < MAXT; i++){
        if(maxlength < ia_params[i][i].len[0])
            maxlength = ia_params[i][i].len[0];
    }

    sqmaxcut += maxlength;
    sqmaxcut *= 1.1;
    maxcut = sqmaxcut;
    sqmaxcut = sqmaxcut*sqmaxcut;
    exter.sqmaxcut +=  maxlength;
    exter.sqmaxcut *= exter.sqmaxcut*1.1;
}

int Topo::convertGeotype(char *geotype) {
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

int Topo::fillSystem(char *pline, char *sysnames[], long **sysmoln, char *name) {
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

int Topo::fillTypes(char **pline) {
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
    if (( (geotype_i == CHPSC) || (geotype_i == CHCPSC) )&& ( fields != 7)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA ATTRACT_DIST ATTRACT_SWITCH PATCH_ANGLE PATCH_SWITCH SC_LENGTH PARALLEL_EPS CHIRAL_ANGLE" << endl;
        return 0;
    }
    if (( (geotype_i == TPSC) || (geotype_i == TCPSC) ) && (fields != 9)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA ATTRACT_DIST ATTRACT_SWITCH PATCH_ANGLE PATCH_SWITCH SC_LENGTH PARALLEL_EPS  PATCH_ROTATION PATCH_ANGLE PATCH_SWITCH" << endl;
        return 0;
    }
    if (( (geotype_i == TCHPSC) || (geotype_i == TCHCPSC) )&& ( fields != 10)) {
        cerr << "TOPOLOGY ERROR: wrong number of parameters for " << geotype << endl;
        cerr << "Parameters are:\n" << "#NAME NUMBER GEOTYPE EPSILON SIGMA ATTRACT_DIST ATTRACT_SWITCH PATCH_ANGLE PATCH_SWITCH SC_LENGTH PARALLEL_EPS  PATCH_ROTATION PATCH_ANGLE PATCH_SWITCH CHIRAL_ANGLE" << endl;
        return 0;
    }

    if ((geotype_i < 0) || (geotype_i > (MAXT + 10))) {
        fprintf (stderr, "TOPOLOGY ERROR: geotype (%s) is out of range: 0 - %d.\n\n", geotype, MAXT + 10);
        return 0;
    }

    strcpy(ia_params[type][type].name, name);
    strcpy(ia_params[type][type].other_name, name);
    ia_params[type][type].geotype[0] = geotype_i;
    ia_params[type][type].geotype[1] = geotype_i;
    ia_params[type][type].epsilon = param[0];
    ia_params[type][type].sigma = param[1];
    ia_params[type][type].sigmaSq   = ia_params[type][type].sigma * ia_params[type][type].sigma;
    ia_params[type][type].A = 4 * ia_params[type][type].epsilon * pow(ia_params[type][type].sigma, 12 );
    ia_params[type][type].B = 4 * ia_params[type][type].epsilon * pow(ia_params[type][type].sigma, 6 );
    ia_params[type][type].rcutwca = (ia_params[type][type].sigma)*pow(2.0,1.0/6.0);
    ia_params[type][type].rcutwcaSq = ia_params[type][type].rcutwca * ia_params[type][type].rcutwca;

    fprintf(stdout, "Topology read of %d: %8s (geotype: %s, %d) with parameters %g %g", type, name, geotype, geotype_i, ia_params[type][type].epsilon, ia_params[type][type].sigma);

    if (fields > 0 && fields != 1 && fields != 3) { // all except SCN and SCA
        ia_params[type][type].pdis = param[2];
        ia_params[type][type].pdisSq = ia_params[type][type].pdis  * ia_params[type][type].pdis;
        ia_params[type][type].pswitch = param[3];
        ia_params[type][type].pswitchINV = 1.0/param[3];
        ia_params[type][type].rcut = ia_params[type][type].pswitch+ia_params[type][type].pdis;
        ia_params[type][type].rcutSq = ia_params[type][type].rcut * ia_params[type][type].rcut;
        fprintf(stdout, " | %g %g",ia_params[type][type].pdis,ia_params[type][type].pswitch);
    }
    if(fields == 1) { // SCN
        for(int i = 0; i < 2; i++){
            ia_params[type][type].len[i] = param[2];
            ia_params[type][type].half_len[i] = param[2] / 2;
        }
    }
    if(fields == 3) { // SCA
        for(int i = 0; i < 2; i++){
            ia_params[type][type].len[i] = param[4];
            ia_params[type][type].half_len[i] = param[4] / 2;
        }
    }
    if (fields > 2 && fields != 3) { // except SCA
        int i;
        for(i = 0; i < 2; i++){
            ia_params[type][type].len[i] = param[6];
            ia_params[type][type].half_len[i] = param[6] / 2;
            ia_params[type][type].pangl[i] = param[4];
            ia_params[type][type].panglsw[i] = param[5];
            ia_params[type][type].pcangl[i] = cos(param[4]/2.0/180*PI);                 // C1
            ia_params[type][type].pcanglsw[i] = cos((param[4]/2.0+param[5])/180*PI);    // C2
            //ia_params[type][type].pcangl[i] = ia_params[type][type].pcangl[i];
            //ia_params[type][type].pcanglsw[i] =	ia_params[type][type].pcanglsw[i];
            ia_params[type][type].pcoshalfi[i] = cos((param[4]/2.0+param[5])/2.0/180*PI);
            ia_params[type][type].psinhalfi[i] = sqrt(1.0 - ia_params[type][type].pcoshalfi[i] * ia_params[type][type].pcoshalfi[i]);
            ia_params[type][type].parallel = param[7];

        }
        fprintf(stdout, " | %g %g | %g", ia_params[type][type].pangl[0], ia_params[type][type].panglsw[0], ia_params[type][type].parallel);
    }
    if(fields == 7){
        int i;
        for(i = 0; i < 2; i++){
            ia_params[type][type].chiral_cos[i] = cos(param[8] / 360 * PI);
            ia_params[type][type].chiral_sin[i] = sqrt(1 - ia_params[type][type].chiral_cos[i] * ia_params[type][type].chiral_cos[i]);
            fprintf(stdout, "| chirality %g ", param[8]);
        }
    }
    if ((fields == 9)||(fields == 10)) {
        int i;
        for(i = 0; i < 2; i++){
            ia_params[type][type].csecpatchrot[i] = cos(param[8] / 360 * PI);
            ia_params[type][type].ssecpatchrot[i] = sqrt(1 - ia_params[type][type].csecpatchrot[i] * ia_params[type][type].csecpatchrot[i]);
            //fprintf(stdout, " | %g %g", ia_params[type][type].csecpatchrot[0], ia_params[type][type].ssecpatchrot[0]);

            ia_params[type][type].pangl[i+2] = param[9];
            ia_params[type][type].panglsw[i+2] = param[10];
            ia_params[type][type].pcangl[i+2] = cos(param[9]/2.0/180*PI);                 // C1
            ia_params[type][type].pcanglsw[i+2] = cos((param[9]/2.0+param[10])/180*PI);    // C2
            //ia_params[type][type].pcangl[i] = ia_params[type][type].pcangl[i];
            //ia_params[type][type].pcanglsw[i] = ia_params[type][type].pcanglsw[i];
            ia_params[type][type].pcoshalfi[i+2] = cos((param[9]/2.0+param[10])/2.0/180*PI);
            ia_params[type][type].psinhalfi[i+2] = sqrt(1.0 - ia_params[type][type].pcoshalfi[i+2] * ia_params[type][type].pcoshalfi[i+2]);
        }
        fprintf(stdout, " | %g  %g %g", param[8], ia_params[type][type].pangl[2], ia_params[type][type].panglsw[2]);
    }
    if(fields == 10){
        int i;
        for(i = 0; i < 2; i++){
            ia_params[type][type].chiral_cos[i] = cos(param[11] / 360 * PI);
            ia_params[type][type].chiral_sin[i] = sqrt(1 - ia_params[type][type].chiral_cos[i] * ia_params[type][type].chiral_cos[i]);
        }
        fprintf(stdout, " | %g ", param[11]);
    }

    // Volume
    if (geotype_i < SP)
        ia_params[type][type].volume = 4.0/3.0*PI*pow((ia_params[type][type].sigma)/2.0,3.0) + PI/2.0*ia_params[type][type].len[0]*pow((ia_params[type][type].sigma)/2.0,2.0) ;
    else
        ia_params[type][type].volume = 4.0/3.0*PI*pow((ia_params[type][type].sigma)/2.0,3.0);
    if ( ia_params[type][type].rcutwca > sqmaxcut )
        sqmaxcut = ia_params[type][type].rcutwca;
    if ( ia_params[type][type].rcut > sqmaxcut )
        sqmaxcut = ia_params[type][type].rcut;
    fprintf(stdout, " \n");
    DEBUG_INIT("Finished filltypes");
    return 1;
}

int Topo::fillExclusions(char **pline, bool exlusions[][MAXT]) {
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

int Topo::fillExter(char **pline) {
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
        exter.exist = true;
        exter.thickness = param[0];
        fprintf(stdout, "External potential with thickness: %le ",exter.thickness);
        if (fields >1) {
            exter.epsilon = param[1];
            fprintf(stdout, "epsilon: %le ",exter.epsilon);
            if (fields >2) {
                exter.attraction = param[2];
                fprintf(stdout, "and range of attraction: %le ",exter.attraction);
            }
        }
    } else{
        exter.exist = false;
        fprintf(stdout, "No external potential ");
    }

    fprintf(stdout, " \n");
    DEBUG_INIT("Finished filling external potential");
    return 1;
}

int Topo::fillMol(char *molname, char *pline, MolIO *molecules) {
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
            moleculeParam[i].name = (char*) malloc(strlen(molname)+1);
            strcpy(moleculeParam[i].name, molname);
            moleculeParam[i].molType = i;
        }

        moleculeParam[i].particleTypes.push_back(molecules[i].type[j]);
        assert(moleculeParam[i].particleTypes[j] == molecules[i].type[j]);

        if (fields == 1){
            (molecules[i].switchtype[j]) = -1;//(molecules[i].type[j]);
            (molecules[i].delta_mu[j]) = 0;
            fields = 3;
        } else{
            fprintf(stdout, "(with switchtype: %ld and delta_mu: %lf)", molecules[i].switchtype[j], molecules[i].delta_mu[j]);
            moleculeParam[i].switchTypes.push_back(molecules[i].switchtype[j]);
            moleculeParam[i].deltaMu.push_back(molecules[i].delta_mu[j]);
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
        moleculeParam[i].bond1c = bondk;
        moleculeParam[i].bond1eq = bonddist;
        fprintf (stdout, "bond1: %f %f \n",moleculeParam[i].bond1c,moleculeParam[i].bond1eq);
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
        moleculeParam[i].bond2c = bondk;
        moleculeParam[i].bond2eq = bonddist;
        fprintf (stdout, "bond2: %f %f \n",moleculeParam[i].bond2c,moleculeParam[i].bond2eq);
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
        moleculeParam[i].bonddc = bondk;
        moleculeParam[i].bonddeq = bonddist;
        fprintf (stdout, "bondd: %f %f \n",moleculeParam[i].bonddc,moleculeParam[i].bonddeq);
        return 1;
    }

    if (!strcmp(molcommand,"BONDH")) {
        fields = sscanf(molparams, "%le %le ", &bondk, &bonddist);
        if (fields < 2) {
            fprintf (stderr, "TOPOLOGY ERROR: wrong number of parameters for bondh, should be 2.\n\n");
            return 0;
        }
        if (bonddist < 0) {
            fprintf (stderr, "TOPOLOGY ERROR: bondhdist cannot be negative: %f \n\n",bonddist);
            return 0;
        }
        moleculeParam[i].bondhc = bondk;
        moleculeParam[i].bondheq = bonddist;
        fprintf (stdout, "bondh: %f %f \n",moleculeParam[i].bondhc,moleculeParam[i].bondheq);
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
        moleculeParam[i].angle1c = bondk;
        moleculeParam[i].angle1eq = bonddist*DEGTORAD;
        fprintf (stdout, "angle1: %f %f \n",moleculeParam[i].angle1c,moleculeParam[i].angle1eq);
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
        moleculeParam[i].angle2c = bondk;
        moleculeParam[i].angle2eq = bonddist*DEGTORAD;
        fprintf (stdout, "angle2: %f %f \n",moleculeParam[i].angle2c,moleculeParam[i].angle2eq);
        return 1;
    }

    // INIT of muVT ensemble
    if (!strcmp(molcommand,"ACTIVITY")) {
        fields = sscanf(molparams, "%le ", &activity);
        moleculeParam[i].activity = activity;
        moleculeParam[i].chemPot = log(activity*Nav*1e-24); // faunus log(activity*Nav*1e-27) [mol/l]
        fprintf (stdout, "activity: %f \n",moleculeParam[i].activity);
        return 1;
    }

    fprintf (stderr, "TOPOLOGY ERROR: unknown parameter: %s.\n\n",molcommand);
    return 0;
}

void Topo::topDealoc() {
    delete[] molecules;

    if (sysmoln != NULL) free(sysmoln);
    sysmoln=NULL;

    if (poolMolNum != NULL) free(poolMolNum);
    poolMolNum=NULL;

    for (int i=0;i<MAXN;i++) {
        if ((sysnames[i]) != NULL) free(sysnames[i]);
        sysnames[i]=NULL;
    }

    for (int i=0;i<MAXN;i++) {
        if ((poolNames[i]) != NULL) free(poolNames[i]);
        poolNames[i]=NULL;
    }
}
