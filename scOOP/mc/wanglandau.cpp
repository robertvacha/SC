#include "wanglandau.h"

#include "stdio.h"

extern Topo topo;

long WangLandau::meshOrderMoveMolecule(Molecule &target, Particle chorig[], Mesh *mesh, int wli) {
    assert(target.size() > 0);
    int change;

    if(target.size() == 1) { // single particle

        int nx,ny,ox,oy; // position in mesh
        double resid;

        if ( conf->pvec[target[0]].type != wlmtype )
            return currorder[wli];

        nx = (int) (INBOX(conf->pvec[target[0]].pos.x,resid) * (*mesh).dim[0]);
        ny = (int) (INBOX(conf->pvec[target[0]].pos.y,resid) * (*mesh).dim[1]);
        ox = (int) (INBOX(chorig[0].pos.x,resid) * (*mesh).dim[0]);
        oy = (int) (INBOX(chorig[0].pos.y,resid) * (*mesh).dim[1]);

        if ( (nx == ox) && (ny == oy) ) // particle stayed in the same mesh bin
            return currorder[wli];
    }

    change= 1;
    for(unsigned int i=0; (i<target.size()) && (change); i++ ) {
        if ( conf->pvec[target[i]].type == wlmtype ) {
            change = mesh->addPart(conf->pvec[target[i]].pos.x, conf->pvec[target[i]].pos.y);
        }
    }
    for(unsigned int i=0; (i<target.size()) && (change); i++ ) {
        if ( conf->pvec[target[i]].type == wlmtype ) {
            change = mesh->removePart(chorig[i].pos.x, chorig[i].pos.y);
        }
    }

    if ( !change ) {
        // fill the mesh with particles
        mesh->meshFill(conf->pvec.size(), wlmtype, &conf->pvec);
        return (long) (mesh->findHoles() - minorder[wli]);
    }

    return currorder[wli];
}

void WangLandau::init(char wlinfile[30]) {
    wl_meshsize = 0;
    radiushole = NULL;
    radiusholeold = NULL;
    radiusholemax = 0;
    partincontactold = 0;
    partincontact = 0;
    wlmdim = 0;
    wlmdim = 0;
    length[0]=0;
    length[1]=0;
    currorder[0]=0;
    currorder[1]=0;
    neworder[0]=0;
    neworder[1]=0;
    weights = NULL;
    hist = NULL;

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
                exit(0);
            }
        }
        if (alpha < WL_ALPHATOL/100) alpha = WL_ZERO;
        fflush (stdout);
    }
}


int WangLandau::write(char filename[]) {
    long i,j;
    FILE *outfile;

    outfile = fopen(filename, "w");
    if (outfile == NULL) {
        cerr << "\nERROR: Could not open " << filename << " file.\n" << endl;
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

int WangLandau::initCalc(char filename[]) {
    long i,length,fields=0;
    double field[5];
    FILE *infile;
    char line[STRLEN];

    infile = fopen(filename, "r");
    if (infile == NULL) {
        cerr << "\nERROR: Could not open " << filename << " file.\n" << endl;
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

#ifdef ENABLE_MPI
    if (mpirank == 0)
    {
       _sizeInt = length * sizeof(long int);
       _sizeDouble = length * sizeof(double);
       _sizeShared = 1 * sizeof(double);
       MPI_Win_allocate_shared(_sizeDouble, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &weights, &_winDouble);
       MPI_Win_allocate_shared(_sizeInt, sizeof(long int), MPI_INFO_NULL, MPI_COMM_WORLD, &hist, &_winInt);
       MPI_Win_allocate_shared(_sizeShared, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shared, &_winShared);
    }
    else
    {
       int disp_unitInt, disp_unitDouble, disp_shared;
       MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &weights, &_winDouble);
       MPI_Win_allocate_shared(0, sizeof(long int), MPI_INFO_NULL, MPI_COMM_WORLD, &hist, &_winInt);
       MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shared, &_winShared);

       MPI_Win_shared_query(_winDouble, 0, &_sizeDouble, &disp_unitDouble, &weights);
       MPI_Win_shared_query(_winInt, 0, &_sizeInt, &disp_unitInt, &hist);
       MPI_Win_shared_query(_winShared, 0, &_sizeShared, &disp_shared, &shared);
    }

    MPI_Barrier(MPI_COMM_WORLD);
#else
    weights = (double*) malloc( sizeof(double) * length );
    hist = (long*) malloc( sizeof(long) * length );
    shared = new double[3];
#endif

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

                    if (mpirank == 0) {
                        weights[i-1] = field[1];
                        hist[i-1] = field[2];
                    }

                    this->length[0]++;
                    i++;
                } else if (fields == 4 ) {
                    if (i==1) {
                        minorder[0] = field[0];
                        minorder[1] = field[1];
                    }
                    if ( minorder[1] == field[1] )
                        this->length[0]++;

                    if (mpirank == 0) {
                        weights[i-1] = field[2];
                        hist[i-1] = field[3];
                    }

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
    mcout.get() << "Wang-Landau method init:" << endl;
    mcout.get() << "alpha: " << alpha << endl;
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


void WangLandau::radiusholePrint(long *radiushole, long length) {
    printf("radiushole:\n");
    for (int i=0;i<length;i++) {
        printf("%ld ",radiushole[i]);
    }
    printf("\n");
}

long WangLandau::meshOrderMoveOne(Vector oldpos, Vector newpos, Mesh *mesh, long npart, long target, int wli) {
    int change;
    int nx,ny,ox,oy; /* position in mesh */
    double resid;

    if ( conf->pvec[target].type != wlmtype )
        return currorder[wli];

    nx = (int) (INBOX(newpos.x,resid) * (*mesh).dim[0]);
    ny = (int) (INBOX(newpos.y,resid) * (*mesh).dim[1]);
    ox = (int) (INBOX(oldpos.x,resid) * (*mesh).dim[0]);
    oy = (int) (INBOX(oldpos.y,resid) * (*mesh).dim[1]);
    if ( (nx == ox) && (ny == oy) ) return currorder[wli]; /* particle stayed in the same mesh bin*/

    change = mesh->addPart(newpos.x,newpos.y);
    if (change) {
        change = mesh->removePart(oldpos.x,oldpos.y);
    }
    if ( !change ) {
        /* fill the mesh with particles*/
        mesh->meshFill(npart, wlmtype, &conf->pvec);
        return (long) (mesh->findHoles() - minorder[wli]);
    }
    return currorder[wli];
}

long WangLandau::meshOrderMoveChain(Molecule chain, Mesh *mesh, long npart, Particle chorig[], int wli) {
    long current;
    int change;

    change= 1;
    for(unsigned int i=0; (i<chain.size()) && (change); i++ ) {
        current = chain[i];
        if ( conf->pvec[current].type == wlmtype ) {
            change = mesh->addPart(conf->pvec[current].pos.x, conf->pvec[current].pos.y);
        }
    }
    for(unsigned int i=0; (i<chain.size()) && (change); i++ ) {
        current = chain[i];
        if ( conf->pvec[current].type == wlmtype ) {
            change = mesh->removePart(chorig[i].pos.x, chorig[i].pos.y);
        }
    }

    if ( !change ) {
        /* fill the mesh with particles*/
        mesh->meshFill(npart, wlmtype, &conf->pvec);
        return (long) (mesh->findHoles() - minorder[wli]);
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
