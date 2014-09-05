#include "topo.h"


void Topo::genParamPairs(bool (*exclusions)[MAXT][MAXT]) {
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
                    ia_params[i][j].epsilon = sqrt(ia_params[i][i].epsilon *  ia_params[j][j].epsilon);
                    ia_params[i][j].pswitch = AVER(ia_params[i][i].pswitch,ia_params[j][j].pswitch);
                    ia_params[i][j].rcutwca = (ia_params[i][j].sigma)*pow(2.0,1.0/6.0);

                    // Averaging of the flat part of attraction
                    ia_params[i][j].pdis = AVER(ia_params[i][i].pdis - ia_params[i][i].rcutwca,
                                                      ia_params[j][j].pdis - ia_params[j][j].rcutwca) + ia_params[i][j].rcutwca;
                    ia_params[i][j].rcut = ia_params[i][j].pswitch+ia_params[i][j].pdis;

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

            if ( (*exclusions)[i][j] )
                ia_params[i][j].epsilon = 0.0;

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
            if ( (*exclusions)[i][j] )
                ia_params[i][j].epsilon = 0.0;
        }
    }
}

void Topo::genTopoParams() {
    double maxlength = 0;
    for(int i = 0; i < MAXT; i++){
        if(maxlength < ia_params[i][i].len[0])
            maxlength = ia_params[i][i].len[0];
    }

    sqmaxcut += maxlength+2;
    sqmaxcut *= 1.1;
    maxcut = sqmaxcut;
    sqmaxcut = sqmaxcut*sqmaxcut;
    exter.sqmaxcut +=  maxlength;
    exter.sqmaxcut *= exter.sqmaxcut*1.1;
}
