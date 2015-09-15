#include "particle.h"

void Particle::init(Ia_param * ia_parami) {

    if ( (ia_parami->geotype[0] == SCA) || (ia_parami->geotype[0] == SCN) ){
        /*SCA and SCN are isotropic... nothing to initialize*/
        return;
    }

    dir.normalise();
    patchdir[0].ortogonalise(dir);
    patchdir[0].normalise();

    // calculate patch sides
    if ( (ia_parami->geotype[0] == PSC) || (ia_parami->geotype[0] == CPSC)
      || (ia_parami->geotype[0] == TPSC) || (ia_parami->geotype[0] == TCPSC)  ) {

        // rotate patch vector by half size of patch
        patchsides[0] = patchdir[0];
        patchsides[0].rotate(dir, ia_parami->pcoshalfi[0], ia_parami->psinhalfi[0]);

        // second side
        patchsides[1] = patchdir[0];
        patchsides[1].rotate(dir, ia_parami->pcoshalfi[0], -1.0*ia_parami->psinhalfi[0]);
    }

    // calculate second patchdir
    if ( (ia_parami->geotype[0] == TPSC) || (ia_parami->geotype[0] == TCPSC) ||
      (ia_parami->geotype[0] == TCHPSC) || (ia_parami->geotype[0] == TCHCPSC)) {

        patchdir[1] = patchdir[0];
        patchdir[1].rotate(dir, ia_parami->csecpatchrot[0], ia_parami->ssecpatchrot[0]);
        patchdir[1].ortogonalise(dir);
        patchdir[1].normalise();
    }

    // calculate second patch sides
    if ( (ia_parami->geotype[0] == TPSC) || (ia_parami->geotype[0] == TCPSC) ) {

        // rotate patch vector by half size of patch
        patchsides[2] = patchdir[1];
        patchsides[2].rotate(dir, ia_parami->pcoshalfi[2], ia_parami->psinhalfi[2]);

        // second side
        patchsides[3] = patchdir[1];
        patchsides[3].rotate(dir, ia_parami->pcoshalfi[2], -1.0*ia_parami->psinhalfi[2]);
    }

    // calculate chdir vector
    if ( (ia_parami->geotype[0] == CHPSC) || (ia_parami->geotype[0] == CHCPSC)
      || (ia_parami->geotype[0] == TCHPSC) || (ia_parami->geotype[0] == TCHCPSC)) {

        chdir[0] = dir;
        chdir[0].rotate(patchdir[0], ia_parami->chiral_cos[0], ia_parami->chiral_sin[0]);

        // rotate patch vector by half size of patch
        patchsides[0] = patchdir[0];
        patchsides[0].rotate(chdir[0], ia_parami->pcoshalfi[0], ia_parami->psinhalfi[0]);

        // second side
        patchsides[1] = patchdir[0];
        patchsides[1].rotate(chdir[0], ia_parami->pcoshalfi[0], -1.0*ia_parami->psinhalfi[0]);
    }

    // calculate chdir vector for seond patch
    if ( (ia_parami->geotype[0] == TCHPSC) || (ia_parami->geotype[0] == TCHCPSC) ) {

        chdir[1] = dir;
        chdir[1].rotate(patchdir[1], ia_parami->chiral_cos[0], ia_parami->chiral_sin[0]);

        // rotate patch vector by half size of patch to get sides
        patchsides[2] = patchdir[1];
        patchsides[2].rotate(chdir[1], ia_parami->pcoshalfi[2], ia_parami->psinhalfi[2]);

        // second side
        patchsides[3] = patchdir[1];
        patchsides[3].rotate(chdir[1], ia_parami->pcoshalfi[2], -1.0*ia_parami->psinhalfi[2]);
    }
}
