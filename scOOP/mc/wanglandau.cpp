#include "wanglandau.h"


long WangLandau::zOrder(int wli) {
    //    printf("%f %ld\n",particleStore[0].pos.z * box.z,lround(particleStore[0].pos.z * box.z / wl.dorder[wli] - wl.minorder[wli]));
    /* Because older C compilators do not know lround we can use ceil as well
       return lround(particleStore[0].pos.z * box.z / wl.dorder[wli] - wl.minorder[wli]);*/
    /*
    printf("%f ",conf->particleStore[0].pos.z );
    printf("%f ",conf->syscm.z);
    printf("%f ",conf->box.z);
    printf("%f ", wl->minorder[wli]);
    printf("%f \n", wl->dorder[wli] );*/

    return (long) ceil( ((conf->particleStore[0].pos.z - conf->syscm.z) * conf->box.z- minorder[wli]) / dorder[wli]  );
}

long WangLandau::twoPartDist(int wli) {
    Vector r_cm;

    r_cm.x = conf->particleStore[0].pos.x - conf->particleStore[1].pos.x;
    r_cm.y = conf->particleStore[0].pos.y - conf->particleStore[1].pos.y;
    r_cm.z = conf->particleStore[0].pos.z - conf->particleStore[1].pos.z;
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
