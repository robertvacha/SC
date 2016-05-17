#include "paire.h"




void PairE::initIntFCE() {
    cout << "\nInitializing energy functions...\n";
    // NB
    // Fill in the names of the functions for calculating the
    // interaction energy
    long geotype, other_geotype;
    for(int i = 0; i < MAXT; i++){
        for(int j = 0; j < MAXT; j++){
            /* Initialize them as not existing */
            intFCE[i][j] = &PairE::eNoExist;
            geotype = topo.ia_params[i][j].geotype[0];
            other_geotype = topo.ia_params[i][j].geotype[1];

            /*if ( ( (geotype == CHCPSC || geotype == CPSC || geotype == TCHCPSC || geotype == TCPSC) &&
                    (other_geotype == CHPSC || other_geotype == PSC || other_geotype == TCHPSC || other_geotype == TPSC) ) ||
                  ( (geotype == CHPSC || geotype == PSC || geotype == TCHPSC || geotype == TPSC)  &&
                    (other_geotype == CHCPSC || other_geotype == CPSC || other_geotype == TCHCPSC || other_geotype == TCPSC) ) )  {
                intFCE[i][j] = &PairE::ePscPsc;
            }
            if ( (geotype == CHCPSC || geotype == CPSC || geotype == TCHCPSC || geotype == TCPSC) &&
                    (other_geotype == CHCPSC || other_geotype == CPSC || other_geotype == TCHCPSC || other_geotype == TCPSC) ){
                intFCE[i][j] = &PairE::eCpscCpsc;
            }
            if ( (geotype == CHPSC || geotype == PSC || geotype == TCHPSC || geotype == TPSC) &&
                    (other_geotype == CHPSC || other_geotype == PSC || other_geotype == TCHPSC || other_geotype == TPSC) ){
                intFCE[i][j] = &PairE::ePscPsc;
            }

            if(geotype == SCN || geotype == SPN
                    || other_geotype == SCN || other_geotype == SPN){
                intFCE[i][j] = &PairE::eSpnOrScn;
            }*/
            if((geotype == SCA && other_geotype == SCA)
                    || (geotype == SPA && other_geotype == SPA)){
                intFCE[i][j] = &PairE::e2ScaOr2Spa;
            }
            /*if((geotype == SCA && other_geotype == SPA)
                    || (geotype == SPA && other_geotype == SCA)){
                intFCE[i][j] = &PairE::eSpaSca;
            }
            if(( (geotype == PSC || geotype == CHPSC || geotype == TCHPSC || geotype == TPSC) && other_geotype == SPA)
                    || (geotype == SPA && (other_geotype == PSC||other_geotype == CHPSC || other_geotype == TCHPSC || other_geotype == TPSC) )){
                intFCE[i][j] = &PairE::ePscSpa;
            }
            if(( (geotype == CPSC ||geotype == CHCPSC || geotype == TCHCPSC || geotype == TCPSC) && other_geotype == SPA)
                    || (geotype == SPA && (other_geotype == CPSC||other_geotype == CHCPSC || other_geotype == TCHCPSC || other_geotype == TCPSC)  )){
                intFCE[i][j] = &PairE::eCpscSpa;
            }*/

        }
    }
}
