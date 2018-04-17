/** @file mesh.cpp*/

#include "mesh.h"

#include "../structures/structures.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int Mesh::meshInit(double meshsize, long npart, int wlmtype, Vector box, std::vector<Particle >* pvec) {
    //  int i;
    int maxsize,length;

    dim[0] = (int)(box.x/meshsize);
    dim[1] = (int)(box.y/meshsize);
    if ( data != NULL ) free(data);
    if ( tmp != NULL ) free(tmp);
    length = dim[0] * dim[1];
    data = (int*) malloc( sizeof(int)* (length));
    tmp = (int*) malloc( sizeof(int)* (length+1));

    /* fill the mesh with particles*/
    meshFill(npart, wlmtype, pvec);
    /* perfrom hole cluster algorithm */
    maxsize = findHoles();

    /*DEBUG    printf("maxsize: %d\n",maxsize);
      printf("mesh:\n");
      for (i=0;i<mesh.dim[0]*mesh.dim[1];i++) {
      printf("%d ",mesh.data[i]);
      if ( ((i+1) % mesh.dim[0]) == 0) printf("\n");
      }*/
    return maxsize;
}

void Mesh::meshFill(long npart, int wlmtype, std::vector<Particle >* pvec) {
    for (int i=0; i<(dim[0] * dim[1]); i++) {
        data[i] = 0;
    }

    for (int i=0; i<npart; i++) {
        /*calculate position of particle on mesh and add it to all where it belongs */
        if ((*pvec)[i].type == wlmtype)
            addPart((*pvec)[i].pos.x, (*pvec)[i].pos.y);
    }
}

int Mesh::addPart(double posx, double posy) {
    int i, squares[9], onhole;
    double resid;

    onhole = 1;
    meshSquare( (int) (INBOX(posx,resid) * dim[0]), (int) (INBOX(posy,resid) * dim[1]) , squares);
    for(i=0;i<9;i++) {
        if ( (squares[i] >= dim[0]*dim[1])||(squares[i] <0) ) {
            printf ("Error: trying to write to %d\n",squares[i]);
            printf ("%d %d and  %d\n", (int) (INBOX(posx,resid) * dim[0]), (int) (INBOX(posy,resid) * dim[1]),i );
            fflush(stdout);
        }
        if ( (data[ squares[i] ]) >= 0 ) onhole =  0;
        data[ squares[i] ]--;
    }
    return onhole;
}

int Mesh::removePart(double posx, double posy) {
    int squares[9];
    double resid;

    meshSquare((int) (INBOX(posx,resid) * dim[0]), (int) (INBOX(posy,resid) * dim[1]) , squares);
    for(int i=0;i<9;i++) {
        //DEBUG	    if (square[i] >= dim[0]*dim[1]) printf ("Error: trying to write to %d\n",square[i]);
        data[ squares[i] ]++;
        if ( (data[ squares[i] ]) == 0 ) return 0;
    }

    return 1;
}

void Mesh::meshSquare(int x, int y, int square[9]) {
    int a,b;

    b=y;
    square[0] = x + dim[0]*b;
    a = x-1;
    if ( a<0 ) a = dim[0]-1;
    square[1] = a + dim[0]*b;
    a = x+1;
    if ( a==dim[0] ) a = 0;
    square[2] = a + dim[0]*b;

    b = y-1;
    if ( b<0 ) b = dim[1]-1;
    square[3] = x + dim[0]*b;
    a = x-1;
    if ( a<0 ) a = dim[0]-1;
    square[4] = a + dim[0]*b;
    a = x+1;
    if ( a==dim[0] ) a = 0;
    square[5] = a + dim[0]*b;

    b = y+1;
    if ( b==dim[1] ) b = 0;
    square[6] = x + dim[0]*b;
    a = x-1;
    if ( a<0 ) a = dim[0]-1;
    square[7] = a + dim[0]*b;
    a = x+1;
    if ( a==dim[0] ) a = 0;
    square[8] = a + dim[0]*b;
}

void Mesh::meshNeighbors(int pos, int neighbors[]) {
    int x,y,a;

    x = pos % dim[0];
    y = pos / dim[0];

    a = x-1;
    if ( a<0 ) a = dim[0]-1;
    neighbors[0] = a + dim[0]*y;
    a = x+1;
    if ( a==dim[0] ) a = 0;
    neighbors[1] = a + dim[0]*y;

    a = y-1;
    if ( a<0 ) a = dim[1]-1;
    neighbors[2] = x + dim[0]*a;
    a = y+1;
    if ( a==dim[1] ) a = 0;
    neighbors[3] = x + dim[0]*a;
}

int Mesh::findHoles() {
    return findHolesDistrib();
    /*int i,j, n, size, li, maxsize;
    int neighbors[4];

    n=0;
    maxsize = 0;
    for (i=0;i<(dim[0] * dim[1]);i++) {
        tmp[i] = 0;
        if ( data[i] > 0 ) data[i] = 0;
    }
    i=0;
    // go through all mesh points
    while ( i < (dim[0] * dim[1]) ) {
        // test if mesh point is occupied
        if ( data[i] != 0 ) {
            ++i;
        } else {
            // mesh point is free, create a new cluster
            n++;
            data[i] = n;
            // start new cluster, put mesh point as first element, and set list pointer on first element
            //DEBUG      if (n >= mesh.dim[0]*mesh.dim[1]) printf ("Error: trying to write to sizes position %d\n",n);
            size = 1;
            tmp[0] = i;
            li = 0;
            // go through all elements of the cluster
            while ( li < size ) {
                //go through all neighbors
                j =  tmp[li];
                meshNeighbors(j, neighbors);
                for ( int k=0; k<4; k++ ) {
                    // test if status is free and append it to the cluster
                    if ( data[ neighbors[k] ] == 0 ) {
                        data[ neighbors[k] ] = n;
                        // append mesh point as element in the list
                        tmp[size] = neighbors[k];
                        size++;
                    }
                    if (  data[ neighbors[k] ] > 0 &&  data[ neighbors[k] ]<n ) {
                        cerr << "Error: Mesh cluster out of range, propably going infinite through pbc." << endl;
                    }
                }
                li++;
            }
            if (size > maxsize)
                maxsize = size;
        }

    }

    return maxsize;*/
}

int Mesh::findHolesDistrib() {
    int i,j, k, n, li, maxsize;
        int neighbors[4], size[1000];
        int hist[200] = {0};

        n=0;
        maxsize = 0;

        for (i=0;i<1000;i++){
            size[i]=0;
            }
        for (i=0;i<(dim[0] * dim[1]);i++) {
            tmp[i] = 0;
            if ( data[i] > 0 ) data[i] = 0;
        }
        i=0;
        // go through all mesh points
        while ( i < (dim[0] * dim[1]) ) {
            // test if mesh point is occupied
            if ( data[i] != 0 ) { i++; }
            else {
                // mesh point is free, create a new cluster
                n++;
                data[i] = n;
                // start new cluster, put mesh point as first element, and set list pointer on first element
                //DEBUG      if (n >= mesh.dim[0]*mesh.dim[1]) printf ("Error: trying to write to sizes position %d\n",n);
                size[n] = 1;
                tmp[0] = i;
                li = 0;
                // go through all elements of the cluster
                while ( li < size[n] ) {
                    //go through all neighbors
                    j =  tmp[li];
                    meshNeighbors(j, neighbors);
                    for ( k=0; k<4; k++ ) {
                        // test if status is free and append it to the cluster
                        if ( data[ neighbors[k] ] == 0 ) {
                            data[ neighbors[k] ] = n;
                            // append mesh point as element in the list
                            tmp[size[n]] = neighbors[k];
                            size[n]++;
                        }
                        if (  data[ neighbors[k] ] > 0 &&  data[ neighbors[k] ]<n ) {
                            fprintf(stderr,"Error: Mesh cluster out of range, propably going infinite through pbc.");
                            fflush(stderr);
                        }
                    }
                    li++;
                }
                if (size[n] > maxsize) maxsize = size[n];
            }

        }

        if(alpha_init == WL_ZERO) {
            ofstream myfile;
            myfile.open ("pore_distrib");

            if( myfile.is_open() ) { // histogram of holes

                for(i=0; i<1000; i++) {
                    j=size[i];
                    if ( j > 0) {
                        hist[j]++;
                    }
                }

                bool zero = true;
                for(i=0; i<200; i++){
                    if( hist[i] != 0 ) {
                        zero = false;
                        break;
                    }
                }

                for(i=0; i<200; i++){
                    myfile << ( (i==0 && zero) ? 1 : hist[i] ) << " \t";
                }
                myfile << endl;
                myfile.close();
            }
        }

        return maxsize;
}

void Mesh::print() {
    int i;

    printf("mesh:\n");
    for (i=0;i<dim[0] * dim[1];i++) {
        printf("%d ",data[i]);
        if ( ((i+1) % dim[0]) == 0) printf("\n");
    }
    printf("hole %d:\n", findHoles() );
    printf("\n");
}

void Mesh::operator=(Mesh &source) {
    if ( this->data != NULL) {
        if ( (this->dim[0] == source.dim[0]) && (this->dim[1] == source.dim[1]) ) {
            memcpy(this->data,source.data, sizeof(int)* (this->dim[0] * this->dim[1])  );
            return;
        } else {
            free (this->data);
            if ( source.dim[0] * source.dim[1] > this->dim[0] * this->dim[1] ) {
                if (this->tmp != NULL ) free (this->tmp);
                this->tmp = (int*) malloc( sizeof(int)* (source.dim[0] * source.dim[1] + 1));
            }
        }
    }
    this->dim[0] = source.dim[0];
    this->dim[1] = source.dim[1];
    this->data = (int*) malloc( sizeof(int)* (this->dim[0] * this->dim[1]));
    if (this->tmp == NULL ) this->tmp = (int*) malloc( sizeof(int)* (source.dim[0] * source.dim[1] + 1));
    memcpy(this->data,source.data, sizeof(int)* (this->dim[0] * this->dim[1])  );

}


int Mesh::end() {
    /* free allocated memory */
    if ( data!= NULL ) free(data);
    if ( tmp!= NULL ) free(tmp);
    return 0;
}


