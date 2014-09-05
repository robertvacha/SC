#include "movecreator.h"


long MoveCreator::meshOrderMoveOne(Vector oldpos, Vector newpos, Mesh *mesh, long npart, long target, int wli) {
    int change;
    int nx,ny,ox,oy; /* position in mesh */
    double resid;

    if ( conf->particleStore[target].type != sim->wl.wlmtype )
        return sim->wl.currorder[wli];

    nx = (int) (INBOX(newpos.x,resid) * (*mesh).dim[0]);
    ny = (int) (INBOX(newpos.y,resid) * (*mesh).dim[1]);
    ox = (int) (INBOX(oldpos.x,resid) * (*mesh).dim[0]);
    oy = (int) (INBOX(oldpos.y,resid) * (*mesh).dim[1]);
    if ( (nx == ox) && (ny == oy) ) return sim->wl.currorder[wli]; /* particle stayed in the same mesh bin*/

    change = Mesh::addPart(newpos.x,newpos.y,&(*mesh).data,(*mesh).dim);
    if (change) {
        change = Mesh::removePart(oldpos.x,oldpos.y,&(*mesh).data,(*mesh).dim);
    }
    if ( !change ) {
        /* fill the mesh with particles*/
        mesh->meshFill(npart, sim->wl.wlmtype);
        return (long) (mesh->findHoles() - sim->wl.minorder[wli]);
    }
    return sim->wl.currorder[wli];
}

long MoveCreator::meshOrderMoveChain(long chain[], Mesh *mesh, long npart, Particle chorig[], int wli) {
    long i,current;
    int change;

    change= 1;
    i = 0;
    current = chain[0];
    while ( (current >=0 ) && (change) ) {
        if ( conf->particleStore[current].type == sim->wl.wlmtype )
            change = Mesh::addPart(conf->particleStore[current].pos.x, conf->particleStore[current].pos.y, &(*mesh).data, (*mesh).dim);
        i++;
        current = chain[i];
    }
    i = 0;
    current = chain[0];
    while ( (current >=0 ) && (change) ) {
        if ( conf->particleStore[current].type == sim->wl.wlmtype )
            change = Mesh::removePart(chorig[i].pos.x, chorig[i].pos.y, &(*mesh).data, (*mesh).dim);
        i++;
        current = chain[i];
    }

    if ( !change ) {
        /* fill the mesh with particles*/
        mesh->meshFill(npart, sim->wl.wlmtype);
        return (long) (mesh->findHoles() - sim->wl.minorder[wli]);
    }
    return sim->wl.currorder[wli];
}
