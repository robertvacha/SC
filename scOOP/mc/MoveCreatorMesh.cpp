#include "movecreator.h"


long MoveCreator::meshOrderMoveOne(Vector oldpos, Vector newpos, Mesh *mesh, long npart, long target, int wli) {
    int change;
    int nx,ny,ox,oy; /* position in mesh */
    double resid;

    if ( conf->pvec[target].type != sim->wl.wlmtype )
        return sim->wl.currorder[wli];

    nx = (int) (INBOX(newpos.x,resid) * (*mesh).dim[0]);
    ny = (int) (INBOX(newpos.y,resid) * (*mesh).dim[1]);
    ox = (int) (INBOX(oldpos.x,resid) * (*mesh).dim[0]);
    oy = (int) (INBOX(oldpos.y,resid) * (*mesh).dim[1]);
    if ( (nx == ox) && (ny == oy) ) return sim->wl.currorder[wli]; /* particle stayed in the same mesh bin*/

    change = mesh->addPart(newpos.x,newpos.y);
    if (change) {
        change = mesh->removePart(oldpos.x,oldpos.y);
    }
    if ( !change ) {
        /* fill the mesh with particles*/
        mesh->meshFill(npart, sim->wl.wlmtype, &conf->pvec);
        return (long) (mesh->findHoles() - sim->wl.minorder[wli]);
    }
    return sim->wl.currorder[wli];
}

long MoveCreator::meshOrderMoveChain(vector<int> chain, Mesh *mesh, long npart, Particle chorig[], int wli) {
    long current;
    int change;

    change= 1;
    current = chain[0];
    for(unsigned int i=0; (i<chain.size()) && (change); i++ ) {
        if ( conf->pvec[current].type == sim->wl.wlmtype ) {
            change = mesh->addPart(conf->pvec[current].pos.x, conf->pvec[current].pos.y);
        }
        current = chain[i];
    }
    current = chain[0];
    for(unsigned int i=0; (i<chain.size()) && (change); i++ ) {
        if ( conf->pvec[current].type == sim->wl.wlmtype )
            change = mesh->removePart(chorig[i].pos.x, chorig[i].pos.y);
        current = chain[i];
    }

    if ( !change ) {
        /* fill the mesh with particles*/
        mesh->meshFill(npart, sim->wl.wlmtype, &conf->pvec);
        return (long) (mesh->findHoles() - sim->wl.minorder[wli]);
    }
    return sim->wl.currorder[wli];
}
