#include "Conf.h"

extern Topo topo;

void Conf::insertMolecule(std::vector<Particle>& molecule) {
#ifndef NDEBUG
    assert(pvec.checkConsistency());
    int size = pvec.size();
    int molTypeSize = pvec.molCountOfType(molecule[0].molType);
#endif

    // add neighbors to neighborList
    if(pairlist_update) {
        for(unsigned int i = 0; i<molecule.size(); i++) {
            neighborList.push_back(Neighbors() );
            neighborList.back().neighborID = (long int*) malloc(sizeof(long) * MAXN);
        }
    }

    pvec.insertMolecule(molecule);

#ifndef NDEBUG
    assert(pvec.checkConsistency());
    assert(pvec.size() == molecule.size() + size);
    assert(pvec.molCountOfType(molecule[0].molType) == molTypeSize + 1);
#endif
}

void Conf::removeMolecule(Molecule &target) {
#ifndef NDEBUG
    assert(pvec.checkConsistency());
    assert((int)pvec.size() == pvec.vecSize());
    unsigned int pSize = pvec.size();
    int molTypeSize = pvec.molCountOfType(pvec[target[0]].molType);
    int molType = pvec[target[0]].molType;
#endif

    if(pairlist_update) {
        for(unsigned int i=0; i<target.size(); i++) {
            delete neighborList.back().neighborID;
            neighborList.pop_back();
        }
    }

    pvec.deleteMolecule(target); // delete from groupList

#ifndef NDEBUG
    assert(pvec.checkConsistency());
    assert(pvec.size() == pSize - target.size());
    assert(pvec.molCountOfType(molType) == molTypeSize - 1);
#endif
}

std::vector<Particle> Conf::getRandomPoolConf(int molType) {
    long target;
    vector<Particle> vec;

    target = ran2() * pool.molCountOfType(molType);
    target = pool.getStoreIndex(molType, target);

    for(int i=0; i < topo.moleculeParam[molType].molSize(); i++) {
        vec.push_back(pool[target+i]);
    }
    return vec;
}

void Conf::massCenter() {
    syscm.x = 0;
    syscm.y = 0;
    syscm.z = 0;
    for (unsigned long i=0; i < pvec.size(); i++) {
        //using periodic boundary conditions
        syscm.x += (pvec[i].pos.x - anInt(pvec[i].pos.x) ) *
            topo.ia_params[pvec[i].type][pvec[i].type].volume;
        syscm.y += (pvec[i].pos.y - anInt(pvec[i].pos.y) ) *
            topo.ia_params[pvec[i].type][pvec[i].type].volume;
        syscm.z += (pvec[i].pos.z - anInt(pvec[i].pos.z) ) *
            topo.ia_params[pvec[i].type][pvec[i].type].volume;
    }
    syscm.x /= sysvolume;
    syscm.y /= sysvolume;
    syscm.z /= sysvolume;
    return;
}

void Conf::partVecInit() {
    for(int i = 0; i < (long)pvec.size(); i++){
        if ( topo.ia_params[pvec[i].type][pvec[i].type].geotype[0]  < SP)
            pvec[i].init(&(topo.ia_params[pvec[i].type][pvec[i].type]));
    }
}

int Conf::overlap(Particle *part1, Particle *part2, Ia_param ia_params[][MAXT]) {

    double b, c, d, e, f;   /* Coefficients in distance quadratic */
    double boundary;     /* Half length of central boundary zone of quadratic */
    double det;
    double halfl;        /* Half length of cylinder */
    double s0, t0;       /* det times location of min separation of infinite lines */
    double ss, tt;       /* Location of min separation of line segments */
    Vector r_cm;  /* Vector between centres of mass */
    double dist;         /* Distance between particles*/
    Vector distvec; /* Distance vector between particles*/

    r_cm = geo.image(&part1->pos, &part2->pos);

    if ((part1->type >= SP) && (part2->type >= SP)) { /*we have two spheres - most common, do nothing*/
        dist=sqrt(DOT(r_cm,r_cm));
    } else {
        if ((ia_params[part1->type][part2->type].geotype[0] < SP) && (ia_params[part1->type][part2->type].geotype[1] < SP)) { /*we have two spherocylinders*/
            /*finding closes contact between them*/
            b = -DOT(part1->dir, part2->dir);
            d =  DOT(part1->dir, r_cm);
            e = -DOT(part2->dir, r_cm);
            f =  DOT(r_cm, r_cm);
            det = 1.0 - b*b;
            //halfl = length / 2.0;
            // Just take the mean
            halfl = ia_params[part1->type][part2->type].half_len[0] = ia_params[part1->type][part2->type].half_len[1];
            halfl /= 2;
            boundary = det * halfl;
            /* Location of smallest separation of the infinite lines */
            s0 = b*e - d;
            t0 = b*d - e;
            /* Location of smallest separation of line segments */
            if (s0 >= boundary) {
                if (t0 >= boundary) {
                    /* Region 2 */
                    if ( d + halfl + halfl*b < 0.0 ) {
                        ss = halfl;
                        tt = linemin( -ss*b - e, halfl );
                    } else {
                        tt = halfl;
                        ss = linemin( -tt*b - d, halfl );
                    }
                } else if (t0 >= -boundary) {
                    /* Region 1 */
                    ss = halfl;
                    tt = linemin( -ss*b - e, halfl );
                } else {
                    /* Region 8 */
                    if ( d + halfl - halfl*b < 0.0 ) {
                        ss = halfl;
                        tt = linemin( -ss*b - e, halfl );
                    } else {
                        tt = -halfl;
                        ss = linemin( -tt*b - d, halfl );
                    }
                }
            } else if (s0 >= -boundary) {
                if (t0 >= boundary) {
                    /* Region 3 */
                    tt = halfl;
                    ss = linemin( -tt*b - d, halfl );
                } else if (t0 >= -boundary) {
                    /* Region 0 */
                    ss = s0/det;
                    tt = t0/det;
                } else {
                    /* Region 7 */
                    tt = -halfl;
                    ss = linemin( -tt*b - d, halfl );
                }
            } else {
                if (t0 >= boundary) {
                    /* Region 4 */
                    if ( d - halfl + halfl*b > 0.0 ) {
                        ss = -halfl;
                        tt = linemin( -ss*b - e, halfl );
                    } else {
                        tt = halfl;
                        ss = linemin( -tt*b - d, halfl );
                    }
                } else if (t0 >= -boundary) {
                    /* Region 5 */
                    ss = -halfl;
                    tt = linemin( -ss*b - e, halfl );
                } else {
                    /* Region 6 */
                    if ( d - halfl - halfl*b > 0.0 ) {
                        ss = -halfl;
                        tt = linemin( -ss*b - e, halfl );
                    } else {
                        tt = -halfl;
                        ss = linemin( -tt*b - d, halfl );
                    }
                }
            }
            /*ss snd tt are Location of min separation of line segments */
            dist=sqrt(f + ss*ss + tt*tt + 2.0*(ss*d + tt*e + ss*tt*b));
        } else {
            if (ia_params[part1->type][part2->type].geotype[0] < SP) { /*We have one spherocylinder -it is first one*/
                //halfl=length/2;/*finding closest vector from sphyrocylinder to sphere*/
                halfl=ia_params[part1->type][part2->type].half_len[0];/*finding closest vector from sphyrocylinder to sphere*/
                c = DOT(part1->dir,r_cm);
                if (c >= halfl) d = halfl;
                else {
                    if (c > -halfl) d = c;
                    else d = -halfl;
                }
                distvec.x = - r_cm.x + part1->dir.x * d;
                distvec.y = - r_cm.y + part1->dir.y * d;
                distvec.z = - r_cm.z + part1->dir.z * d;
                dist=sqrt(DOT(distvec,distvec));
            } else { /*lst option first one is sphere second one spherocylinder*/
                //halfl=length/2; /*finding closest vector from sphyrocylinder to sphere*/
                halfl=ia_params[part1->type][part2->type].half_len[1];/*finding closest vector from sphyrocylinder to sphere*/
                c = DOT(part2->dir,r_cm);
                if (c >= halfl) d = halfl;
                else {
                    if (c > -halfl) d = c;
                    else d = -halfl;
                }
                distvec.x = r_cm.x - part2->dir.x * d;
                distvec.y = r_cm.y - part2->dir.y * d;
                distvec.z = r_cm.z - part2->dir.z * d;
                dist=sqrt(DOT(distvec,distvec));
            }
        }
    }

    /* Overlap exists if smallest separation is less than diameter of cylinder */
    if (dist < ia_params[part1->type][part2->type].sigma*0.5 ) {
        return 1;
    } else {
        return 0;
    }
}




bool Conf::overlapAll(Particle* target, Ia_param ia_params[][MAXT]) {
    for (unsigned long i=0; i<pvec.size(); i++) {
        if (&pvec[i] != target) {
            if ( overlap(target, &pvec[i], ia_params) ) {
                return true;
            }
        }
    }
    return false;
}




int Conf::checkall(Ia_param ia_params[][MAXT]) {
    unsigned long i, j;

    for (i=0; i<pvec.size()-1; i++) {
        for (j=i+1; j<pvec.size(); j++) {
            if ( overlap(&pvec[i], &pvec[j], ia_params) ) {
                return 1;
            }
        }
    }
    return 0;
}

