#include "paire.h"




void PairE::initIntFCE() {
    cout << "\nInitializing energy functors...\n";
    // NB
    // Fill in the names of the functions for calculating the
    // interaction energy
    long geotype, other_geotype;
    for(int i = 0; i < MAXT; i++){
        for(int j = 0; j < MAXT; j++){
            // Initialize them as not existing
            eFce[i][j] = new EBasic();
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
            }*/

            if ( (geotype == CHPSC || geotype == TCHPSC) && (other_geotype == CHPSC || other_geotype == TCHPSC) ){
                eFce[i][j] = new Psc<WcaCos2, EBond, EAngle>();
            }

            if ( (geotype == PSC || geotype == TPSC) && (other_geotype == PSC || other_geotype == TPSC) ){ // PSC-PSC, PSC-TPSC, TPSC-PSC, TPSC-TPSC
                eFce[i][j] = new Psc<WcaCos2, EBond, EAngle>();
            }

            if(geotype == SCN || geotype == SPN
                    || other_geotype == SCN || other_geotype == SPN){
                eFce[i][j] = new Sphere<WcaTrunc, HarmonicSp, EAngle>();
            }
            if((geotype == SCA && other_geotype == SCA)
                    || (geotype == SPA && other_geotype == SPA)){
                eFce[i][j] = new Sphere<WcaCos2, HarmonicSp, EAngle>();
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




Vector SpheroCylinder::minDistSegments(const Vector &segA, const Vector &segB, double halfl1, double halfl2, const Vector &r_cm) {
    Vector u,v,w,vec;
    double a,b,c,d,e,D,sc,sN,sD,tc,tN,tD;
    bool paralel = false;

    // direction of lines
    u = (2.0*halfl1) * segA; //S1.P1 - S1.P0;
    v = (2.0*halfl2) * segB; //S2.P1 - S2.P0;

    w.x = segB.x*halfl2 - segA.x*halfl1 - r_cm.x;
    w.y = segB.y*halfl2 - segA.y*halfl1 - r_cm.y;
    w.z = segB.z*halfl2 - segA.z*halfl1 - r_cm.z; //S1.P0 - S2.P0;

    a = DOT(u,u);        // always >= 0
    b = DOT(u,v);
    c = DOT(v,v);        // always >= 0
    d = DOT(u,w);
    e = DOT(v,w);
    D = a*c - b*b;       // always >= 0
    sc = D;
    sN = D;
    sD = D;      // sc = sN / sD, default sD = D >= 0
    tc = D;
    tN = D;
    tD = D;      // tc = tN / tD, default tD = D >= 0

    // compute the linetopo.ia_paramseters of the two closest points
    if (D < 0.00000001) { // the lines are almost parallel
        paralel = true;
        sN = 0.0;        // force using point P0 on segment S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    } else {                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    if (fabs(sN) < 0.00000001) sc = 0.0 ;
    else sc = sN / sD;
    if (fabs(tN) < 0.00000001) tc = 0.0 ;
    else tc = tN / tD;

    // get the difference of the two closest points
    //Vector = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
    vec.x = u.x*sc + w.x - v.x*tc;
    vec.y = u.y*sc + w.y - v.y*tc;
    vec.z = u.z*sc + w.z - v.z*tc;

    if(paralel) { // force using point P1 on line S0 and compare the distance to point P0 on segment S1
        Vector vec2;
        // direction of lines
        // note u is v and v is u

        w.x = segA.x*halfl1 - segB.x*halfl2 + r_cm.x;
        w.y = segA.y*halfl1 - segB.y*halfl2 + r_cm.y;
        w.z = segA.z*halfl1 - segB.z*halfl2 + r_cm.z; //S1.P0 - S2.P0;

        // a is c
        // b is same
        // c is a
        d = DOT(v,w); // recalc
        e = DOT(u,w); // recalc
        D = a*c - b*b;
        sc = D;
        sN = D;
        sD = D;      // sc = sN / sD, default sD = D >= 0
        tc = D;
        tN = D;
        tD = D;      // tc = tN / tD, default tD = D >= 0

        // compute the linetopo.ia_paramseters of the two closest points
        if (D < 0.00000001) { // the lines are almost parallel
            sN = 0.0;        // force using point P0 on segment S1
            sD = 1.0;        // to prevent possible division by 0.0 later
            tN = e;
            tD = a;
        }

        if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
            tN = 0.0;
            // recompute sc for this edge
            if (-d < 0.0)
                sN = 0.0;
            else if (-d > c)
                sN = sD;
            else {
                sN = -d;
                sD = c;
            }
        } else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
            tN = tD;
            // recompute sc for this edge
            if ((-d + b) < 0.0)
                sN = 0;
            else if ((-d + b) > c)
                sN = sD;
            else {
                sN = (-d + b);
                sD = c;
            }
        }

        // finally do the division to get sc and tc
        if (fabs(sN) < 0.00000001) sc = 0.0 ;
        else sc = sN / sD;
        if (fabs(tN) < 0.00000001) tc = 0.0 ;
        else tc = tN / tD;

        // get the difference of the two closest points
        //Vector = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
        vec2.x = v.x*sc + w.x - u.x*tc;
        vec2.y = v.y*sc + w.y - u.y*tc;
        vec2.z = v.z*sc + w.z - u.z*tc;

        if(vec2.dot(vec2) < vec.dot(vec))
            return vec2;
    }

    return vec;
}


