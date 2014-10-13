#include "Conf.h"

void Conf::addMolecule(std::vector<Particle>* molecule) {
    int insID = pvecGroupList.getInsertIndex((*molecule)[0].molType);   // store index where to insert

    // insert at end of pvec -> trivial
    if((*molecule)[0].molType == pvecGroupList.molTypeCount-1) {
        for(unsigned int i=0; i<molecule->size(); i++)
            pvec.push_back((*molecule)[i]);

    } else { // insert in middle of particleVector
        //
        // copy all particles after insert (done automatically by std::vector::insert() )
        // optimalization, when possible copy only minimal number of particles of succeeding molTypes
        //
        pvec.insert(pvec.begin()+insID, molecule->begin(), molecule->end());
    }

    // add neighbors to neighborList
    if(pairlist_update) {
        for(unsigned int i = 0; i<molecule->size(); i++) {
            neighborList.push_back(Neighbors());
            neighborList.back().neighborID = (long int*) malloc(sizeof(long) * MAXN);
        }
    }
    // change groupList
    pvecGroupList.insertMolecule((*molecule)[0].molType);
}

void Conf::removeMolecule(int target, int size) {

    //if(pvec[(*target)[0]].molType == molTypeCount-1) {


    //} else { // erasing in middle of other types
        //
        // copy all particles after delete (done automatically by std::vector::erase() )
        //
        // optimalization, when possible copy only minimal number of particles of succeeding molTypes
        //

        pvec.erase(pvec.begin()+target, pvec.begin()+target+size);

        if(pairlist_update) {
            for(int i=0; i<size; i++) {;
                delete neighborList.back().neighborID;
                neighborList.pop_back();
            }
        }
    //}
    // modify groupList
    pvecGroupList.deleteMolecule(pvec[target].molType);
}

void Conf::massCenter(Topo* topo) {
    syscm.x = 0;
    syscm.y = 0;
    syscm.z = 0;
    for (unsigned long i=0; i < pvec.size(); i++) {
        //using periodic boundary conditions
        syscm.x += (pvec[i].pos.x - anInt(pvec[i].pos.x) ) *
            topo->ia_params[pvec[i].type][pvec[i].type].volume;
        syscm.y += (pvec[i].pos.y - anInt(pvec[i].pos.y) ) *
            topo->ia_params[pvec[i].type][pvec[i].type].volume;
        syscm.z += (pvec[i].pos.z - anInt(pvec[i].pos.z) ) *
            topo->ia_params[pvec[i].type][pvec[i].type].volume;
    }
    syscm.x /= sysvolume;
    syscm.y /= sysvolume;
    syscm.z /= sysvolume;
    return;
}

void Conf::partVecInit(Topo* topo) {
    for(int i = 0; i < (long)pvec.size(); i++){
        if ( topo->ia_params[pvec[i].type][pvec[i].type].geotype[0]  < SP)
            pvec[i].init(&(topo->ia_params[pvec[i].type][pvec[i].type]));
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

    r_cm = image(&part1->pos, &part2->pos, &box);

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

