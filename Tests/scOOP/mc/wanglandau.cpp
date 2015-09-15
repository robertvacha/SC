#include "wanglandau.h"

#include "stdio.h"

extern Topo topo;

void WangLandau::init(char wlinfile[30]) {
    if ( wlm[0] >0 ) {
        if (initCalc(wlinfile) != 0)
            return;
        wlmdim = 1 ;
        if ( wlm[1] > 0 )
            wlmdim = 2 ;
        for (long wli=0; wli < wlmdim; wli++) {
            switch (wlm[wli]) {
                case 1:
                    conf->massCenter();
                    currorder[wli] = zOrder(wli);
                    break;
                case 2:
                    wl_meshsize = (topo.ia_params[wlmtype][wlmtype].sigma) / 3.0; // TODO
                    mesh.data = NULL;
                    mesh.tmp = NULL;
                    origmesh.data = NULL;
                    origmesh.tmp = NULL;
                    currorder[wli] = (long) (mesh.meshInit(wl_meshsize,
                                                                           (long)conf->pvec.size(),
                                                                           wlmtype, conf->geo.box,
                                                                           &conf->pvec) - minorder[wli]);
                    break;
                case 3:
                    currorder[wli] = (long) floor( (conf->pvec[0].dir.z - minorder[wli])/ dorder[wli] );
                    break;
                case 4:
                    currorder[wli] = twoPartDist(wli);
                    break;
                case 5:
                    conf->massCenter();
                    radiusholemax = 0;
                    radiushole = NULL;
                    radiusholeold = NULL;
                    currorder[wli] = radiusholeAll(wli,&(conf->syscm));
                    break;
                case 6:
                    radiusholemax = 0;
                    radiushole = NULL;
                    radiusholeold = NULL;
                    currorder[wli] = radiusholeAll(wli,&(conf->pvec[0].pos));
                    break;
                case 7:
                    currorder[wli] = contParticlesAll(wli);
                    break;
                default:
                    currorder[wli] = 0;
                    break;
            }
            if ( (currorder[wli] >= length[wli] ) || (currorder[wli] < 0) ) {
                printf("Error: starting Wang-Landau method with order parameter %f out of range(%f - %f)\n\n", dorder[wli]*currorder[wli] + \
                   minorder[wli], minorder[wli], minorder[wli]+dorder[wli]*length[wli]  );
                end();
                return;
            }
        }
        if (alpha < WL_ALPHATOL/100) alpha = WL_ZERO;
        fflush (stdout);
    }
}

long WangLandau::zOrder(int wli) {
    //    printf("%f %ld\n",pvec[0].pos.z * box.z,lround(pvec[0].pos.z * box.z / dorder[wli] - minorder[wli]));
    /* Because older C compilators do not know lround we can use ceil as well
       return lround(pvec[0].pos.z * box.z / dorder[wli] - minorder[wli]);*/

    /*if(cond) {
    printf("pos Z %f ",conf->pvec[0].pos.z );
    printf("%f ",conf->syscm.z);
    printf("%f ",conf->box.z);
    printf("%f ", wl->minorder[wli]);
    printf("dorder %f \n", wl->dorder[wli] );
    }*/

    return (long) ceil( ((conf->pvec[0].pos.z - conf->syscm.z) * conf->geo.box.z- minorder[wli]) / dorder[wli]  );
}

long WangLandau::twoPartDist(int wli) {
    Vector r_cm;


    r_cm = conf->geo.image(&conf->pvec[0].pos, &conf->pvec[1].pos);
/*
    r_cm.x = conf->pvec[0].pos.x - conf->pvec[1].pos.x;
    r_cm.y = conf->pvec[0].pos.y - conf->pvec[1].pos.y;
    r_cm.z = conf->pvec[0].pos.z - conf->pvec[1].pos.z;
    if ( r_cm.x < 0  )
        r_cm.x = conf->box.x * (r_cm.x - (double)( (long)(r_cm.x-0.5) ) );
    else
        r_cm.x = conf->box.x * (r_cm.x - (double)( (long)(r_cm.x+0.5) ) );
    if ( r_cm.y < 0  )
        r_cm.y = conf->box.y * (r_cm.y - (double)( (long)(r_cm.y-0.5) ) );
    else
        r_cm.y = conf->box.y * (r_cm.y - (double)( (long)(r_cm.y+0.5) ) );
    if ( r_cm.z < 0  )
        r_cm.z = conf->box.z * (r_cm.z - (double)( (long)(r_cm.z-0.5) ) );
    else
        r_cm.z = conf->box.z * (r_cm.z - (double)( (long)(r_cm.z+0.5) ) );
*/

    return (long) ceil( (( sqrt(r_cm.x*r_cm.x + r_cm.y*r_cm.y) ) - minorder[wli]) / dorder[wli]  );
}

int WangLandau::write(char filename[]) {
    long i,j;
    FILE *outfile;

    outfile = fopen(filename, "w");
    if (outfile == NULL) {
        fprintf (stderr, "\nERROR: Could not open %s file.\n\n",filename);
        return 1;
    }
    fprintf (outfile, "%15.8e \n",alpha);
    if (length[1] == 0) {
        for (i=0; i<length[0]; i++) {
            fprintf (outfile, "%15.8e %15.8e %ld \n",minorder[0] + i * dorder[0], weights[i], hist[i]);
        }
    } else {
        for (j=0; j<length[1]; j++) {
            for (i=0; i<length[0]; i++) {
                fprintf (outfile, "%15.8e %15.8e %15.8e %ld \n",minorder[0] + i * dorder[0], minorder[1]+j*dorder[1], weights[i+length[0]*j], hist[i+length[0]*j]);
            }
            fprintf (outfile, " \n");
        }
    }

    fflush(outfile);
    fclose(outfile);

    return 0;
}



void WangLandau::reject(long oldlength, int wlm[2]) {
    if ( wlm[0] > 0 ) {
        weights[currorder[0]+currorder[1]*length[0]] -= alpha;
        hist[currorder[0]+currorder[1]*length[0]]++;
        if ( (wlm[0] == 2) || (wlm[1] == 2) )
            this->mesh = this->origmesh;
        if ( (wlm[0] == 5) || (wlm[1] == 5)||(wlm[0] == 6) || (wlm[1] == 6) ) {
            longarrayCpy(&radiushole,&radiusholeold,radiusholemax,oldlength);
            radiusholemax = oldlength;
        }
        partincontact = partincontactold;
    }
}

void WangLandau::accept(int wlm) {
    if ( wlm > 0 ) {
        for (int i=0;i<2;i++)
            currorder[i] = neworder[i];
        weights[ currorder[0] + currorder[1] * length[0]] -= alpha;
        hist[ currorder[0] + currorder[1] * length[0]]++;
    }
}

int WangLandau::end() {
    free(weights);
    free(hist);
    return 0;
}


int WangLandau::initCalc(char filename[]) {
    long i,length,fields=0;
    double field[5];
    FILE *infile;
    char line[STRLEN];

    infile = fopen(filename, "r");
    if (infile == NULL) {
        fprintf (stderr, "\nERROR: Could not open %s file.\n\n",filename);
        return 1;
    }
    length=0;
    while (wlfgets(line,STRLEN-2,infile) != NULL) {
        strip_comment(line);
        trim(line);
        /* if there is something left... */
        if ((int)strlen(line) > 0) {
            length++;
        }
    }
    length--; /*there is alpha at the first line*/
    weights = (double*) malloc( sizeof(double) * length );
    hist = (long*) malloc( sizeof(long) * length );
    this->length[1] = 0;
    dorder[1] = 0;

    fseek(infile,0,SEEK_SET);
    i=0;
    while (wlfgets(line,STRLEN-2,infile) != NULL) {
        strip_comment(line);
        trim(line);
        /* if there is something left... */
        if ((int)strlen(line) > 0) {
            if (i == 0) {
                if (sscanf(line, "%le",&alpha)!= 1) {
                    fprintf (stderr, "ERROR: Could not read alpha at the begining.\n\n");
                    this->end();
                    return 1;
                } else i++;
            } else {
                fields = sscanf(line, "%le %le %le %le",&field[0],&field[1],&field[2],&field[3]);
                if ( fields == 3 ) {
                    if (i==1)
                        minorder[0] = field[0];
                    weights[i-1] = field[1];
                    hist[i-1] = field[2];
                    this->length[0]++;
                    i++;
                } else if (fields == 4 ) {
                    if (i==1) {
                        minorder[0] = field[0];
                        minorder[1] = field[1];
                    }
                    if ( minorder[1] == field[1] )
                          this->length[0]++;
                    weights[i-1] = field[2];
                    hist[i-1] = field[3];
                    i++;
                } else {
                    fprintf (stderr, "ERROR: Could not read order parameter at line %ld.\n\n", i);
                    this->end();
                    return 1;
                }
            }
        }
    }

    if (fields == 4 ) {
        this->length[1] = length / this->length[0];
        dorder[1] = (field[1] - minorder[1])/(this->length[1]-1);
    }
    dorder[0] = (field[0] - minorder[0])/(this->length[0]-1);
    if  ( ( (i-1) != this->length[0] ) && (fields==3) )  {
        fprintf (stderr, "ERROR: In reading order parameters length %ld does not fit number of lines %ld.\n\n", this->length[0],i-1);
        this->end();
        return 1;
    }
    if  ( ( (i-1) != this->length[0]*this->length[1] ) && (fields==4) )  {
        fprintf (stderr, "ERROR: In reading order parameters lengths %ld %ld does not fit number of lines %ld.\n\n", this->length[0],this->length[1],i-1);
        this->end();
        return 1;
    }
    /*DEBUG*/
    printf("Wang-Landau method init:\n");
    printf("alpha: %f\n",alpha);
    /*int j=0;
    if (length[1] == 0) {
        for (i=0; i<length[0]; i++) {
            printf ("%15.8e %15.8e %ld \n",minorder[0] + i * dorder[0], weights[i], hist[i]);
        }
    } else {
        for (j=0; j<length[1]; j++) {
            for (i=0; i<length[0]; i++) {
                printf ("%15.8e %15.8e %15.8e %ld \n",minorder[0] + i * dorder[0], minorder[1]+j*dorder[1], weights[i+length[0]*j], hist[i+length[0]*j]);
            }
            printf (" \n");
        }
    }*/
    fclose(infile);
    fflush(stdout);
    /**/
    return 0;
}

long WangLandau::radiusholeAll(int wli, Vector *position) {
    long i,nr,radiusholemax;
    double rx,ry,z;

    radiusholemax = radiusholePosition(sqrt(conf->geo.box.x * conf->geo.box.x + conf->geo.box.y * conf->geo.box.y),wli);
    if ( radiusholemax > radiusholemax ) {
        if (radiushole != NULL)
            free(radiushole);
        radiushole = (long int*) malloc( sizeof(long)* (radiusholemax));
        radiusholemax = radiusholemax;
    }

    for (i=0;i<radiusholemax;i++) {
        radiushole[i] = 0;
    }

    for (i=0; i< (long)conf->pvec.size(); i++) {
        /*calculate position of particle from z axis, and add it in array */
        if ( conf->pvec[i].type == wlmtype ) {
            z=conf->pvec[i].pos.z - (*position).z; /*if above position*/
            if (z-anInt(z) > 0) {
                rx = conf->geo.box.x * (conf->pvec[i].pos.x - anInt(conf->pvec[i].pos.x));
                ry = conf->geo.box.y * (conf->pvec[i].pos.y - anInt(conf->pvec[i].pos.y));
                nr = radiusholePosition(sqrt(rx*rx+ry*ry), wli);
                if (nr < 0)
                    return -100;
                radiushole[nr]++;
            }
        }
    }

    return radiusholeOrder();
}

long WangLandau::contParticlesAll(int wli) {
    partincontact = 0;

    for (int i=1; i< (long)conf->pvec.size(); i++) {
        /*calculate position of particle and add it if in contact */
        if ( conf->pvec[i].type == wlmtype ) {
            if ( particlesInContact (&(conf->pvec[i].pos)) )
                partincontact++;
        }
    }

    return contParticlesOrder(wli);
}

long WangLandau::contParticlesMoveOne(Vector *oldpos, long target, int wli) {
    if ( conf->pvec[target].type != wlmtype )
        return currorder[wli];

    if ( particlesInContact (&(conf->pvec[target].pos)) )
        partincontact++;
    if ( particlesInContact (oldpos) )
        partincontact--;

    return contParticlesOrder(wli);
}

long WangLandau::contParticlesMoveChain(Molecule chain, Particle chorig[], int wli) {

    for(unsigned int i=0; i<chain.size(); i++ ) {
        if ( conf->pvec[chain[i]].type == wlmtype ) {
            if ( particlesInContact (&(conf->pvec[chain[i]].pos)) )
                partincontact++;
        }
    }
    for(unsigned int i=0; i<chain.size(); i++ ) {
        if ( conf->pvec[chain[i]].type == wlmtype ) {
            if ( particlesInContact (&(chorig[i].pos)) )
                partincontact--;
        }
    }

    return contParticlesOrder(wli);
}

void WangLandau::radiusholePrint(long *radiushole, long length) {
    printf("radiushole:\n");
    for (int i=0;i<length;i++) {
        printf("%ld ",radiushole[i]);
    }
    printf("\n");
}




long WangLandau::radiusholeOrderMoveOne(Vector *oldpos, long target, int wli, Vector *position) {
    long nr,oor; /* position in radiushole */
    double rx,ry,z;
    bool oz,nz;

    if ( conf->pvec[target].type != wlmtype )
        return currorder[wli];

    z=conf->pvec[target].pos.z - position->z; /*if above position*/
    if (z-anInt(z) < 0) nz = false;
    else nz=true;
    z=oldpos->z - position->z;  /*if above position*/
    if (z-anInt(z) < 0) oz = false;
    else oz=true;
    if ( !(nz) && !(oz) )
        return currorder[wli];

    rx = conf->geo.box.x * (conf->pvec[target].pos.x - anInt(conf->pvec[target].pos.x));
    ry = conf->geo.box.y * (conf->pvec[target].pos.y - anInt(conf->pvec[target].pos.y));
    nr = radiusholePosition(sqrt(rx*rx+ry*ry),wli);

    if (nr < 0)
        return -100;
    /*particle move over radius bins*/
    if (nz) {
        radiushole[nr]++;
    }
    if (oz) {
        rx = conf->geo.box.x * (oldpos->x - anInt(oldpos->x));
        ry = conf->geo.box.y * (oldpos->y - anInt(oldpos->y));
        oor = radiusholePosition(sqrt(rx*rx+ry*ry),wli);
        radiushole[oor]--;
        if ( radiushole[oor] < 0 ) {
            printf ("Error(single particle move): trying to make number of beads in radiuspore smaller than 0 at position %ld\n",oor);
            radiusholePrint(radiushole,radiusholemax);
            fflush(stdout);
        }
        if (radiushole[oor] ==0)
            return radiusholeOrder();
    }

    if ( (nz) && (radiushole[nr] ==1) )  {
        return radiusholeOrder();
    }

    return currorder[wli];
}

long WangLandau::radiusholeOrderMoveChain(Molecule chain, Particle chorig[], int wli, Vector *position) {
    long nr;
    double rx,ry,z;
    bool change=false;

    rx=0;
    for(unsigned int i=0; i<chain.size(); i++ ) {
        if ( conf->pvec[chain[i]].type == wlmtype ) {
            z=conf->pvec[chain[i]].pos.z - position->z; /*if above system CM*/
            if (z-anInt(z) > 0) {
                rx = conf->geo.box.x * (conf->pvec[chain[i]].pos.x - anInt(conf->pvec[chain[i]].pos.x));
                ry = conf->geo.box.y * (conf->pvec[chain[i]].pos.y - anInt(conf->pvec[chain[i]].pos.y));
                nr = radiusholePosition(sqrt(rx*rx+ry*ry),wli);
                if (nr < 0)
                    return -100;
                radiushole[nr]++;
                if ( radiushole[nr] == 1 ) change = true;
            }
        }
    }
    for(unsigned int i=0; i<chain.size(); i++ ) {
        if ( conf->pvec[chain[i]].type == wlmtype ) {
            z=chorig[i].pos.z - position->z; /*if above system CM*/
            if (z-anInt(z) > 0) {
                rx = conf->geo.box.x * (chorig[i].pos.x - anInt(chorig[i].pos.x));
                ry = conf->geo.box.y * (chorig[i].pos.y - anInt(chorig[i].pos.y));
                nr = radiusholePosition(sqrt(rx*rx+ry*ry),wli);
                radiushole[nr]--;
                if ( radiushole[nr] < 0 ) {
                    printf ("Error (chainmove): trying to make number of beads in radiuspore smaller than 0 at position %ld\n",nr);
                    radiusholePrint(radiushole,radiusholemax);
                    fflush(stdout);
                }
                if ( radiushole[nr] == 0 ) change = true;
            }
        }
    }

    if ( change ) {
        return radiusholeOrder();
    }
    return currorder[wli];
}




void WangLandau::endWangLandau(char wloutfile[30]) {
    long i=0,j=0;
    if (wlm[0] > 0) {
        min = hist[0];
        for (i=0;i < length[0];i++) {
            j=0;
            if ( hist[i+j*length[0]] < min ) min = hist[i+j*length[0]];
            for (j=1;j < length[1];j++) {
                if ( hist[i+j*length[0]] < min ) min = hist[i+j*length[0]];
            }
        }
        wmin = weights[0];
        for (i=0;i < length[0];i++) {
            j=0;
            weights[i+j*length[0]] -= wmin;
            for (j=1;j < length[1];j++) {
                weights[i+j*length[0]] -= wmin;
            }
        }
        write(wloutfile);
        end();
        if ( (wlm[0] == 2)||(wlm[1] == 2) ) {
            mesh.end();
            origmesh.end();
        }
        if ( (wlm[0] == 5)||(wlm[1] == 5)||(wlm[0] == 6)||(wlm[1] == 6)  ) {
            if ( radiushole != NULL ) free(radiushole);
            if ( radiusholeold != NULL ) free(radiusholeold);
        }
    }
}





