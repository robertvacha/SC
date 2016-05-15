/** @file statistics.h*/

#ifndef STATISTICS_H
#define STATISTICS_H

#include "macros.h"

class Statistics {
public:
    Statistics (){}

    StatsGrand grand[MAXMT];

    Disp edge;                  ///< \brief Maximum box length change and statistics -> Pressure move
    Disp rot[MAXT];             ///< \brief Maximum rotation and statistics
    Disp trans[MAXT];           ///< \brief Maximum translation  and statistics
    Disp switchMv[MAXT];        ///< \brief Switch move statistics
    Disp chainm[MAXMT];         ///< \brief Maximum translation for chain  and statistics
    Disp chainr[MAXMT];         ///< \brief Maximum rotation for chain and statistics
    Disp mpiexch;               ///< \brief MPI statistics

    int accSwitch() {
        int var=0;
        for(int i=0; i<MAXT; i++)
            var += switchMv[i].acc;
        return var;
    }

    int rejSwitch() {
        int var=0;
        for(int i=0; i<MAXT; i++)
            var += switchMv[i].rej;
        return var;
    }

    int accSTrans() {
        int var=0;
        for(int i=0; i<MAXT; i++)
            var += trans[i].acc;
        return var;
    }

    int rejSTrans() {
        int var=0;
        for(int i=0; i<MAXT; i++)
            var += trans[i].rej;
        return var;
    }

    int stepsSTrans() { return accSTrans() + rejSTrans(); }
    int stepsSwitch() { return accSwitch() + rejSwitch(); }

    int accCTrans() {
        int var=0;
        for(int i=0; i<MAXMT; i++)
            var += chainm[i].acc;
        return var;
    }

    int rejCTrans() {
        int var=0;
        for(int i=0; i<MAXMT; i++)
            var += chainm[i].rej;
        return var;
    }

    int stepsCTrans() { return accCTrans() + rejCTrans(); }

    int accSRot() {
        int var=0;
        for(int i=0; i<MAXT; i++) {
            var += rot[i].acc;
        }
        return var;
    }

    int rejSRot() {
        int var=0;
        for(int i=0; i<MAXT; i++) {
            var += rot[i].rej;
        }
        return var;
    }

    int stepsSRot() { return accSRot() + rejSRot(); }

    int accCRot() {
        int var=0;
        for(int i=0; i<MAXMT; i++) {
            var += chainr[i].acc;
        }
        return var;
    }

    int rejCRot() {
        int var=0;
        for(int i=0; i<MAXMT; i++) {
            var += chainr[i].rej;
        }
        return var;
    }

    int stepsCRot() { return accCRot() + rejCRot(); }

    void printEqStat(Disp *dat, double scale, int length) {
        for(int i=0; i<length; i++) {
            if (RATIO(dat[i]) > 0)
                printf ("   TYPE %d           %.6f  /  %.6f\n", i, dat[i].mx/scale,RATIO(dat[i]));
        }
    }

    void printEqStat() {

        printf ("   Equilibrated maximum displacement / acceptance ratio:            \n");
        printEqStat(trans,2.0,MAXT);

        printf ("   Equilibrated maximum rotation / acceptance ratio:                       \n");
        printEqStat(rot,1.0,MAXT);

        printf ("   Equilibrated maximum box length change / acceptance ratio:              \n");
        printf ("                     %.6e  /  %.6e\n", edge.mx/2.0,RATIO(edge));

        printf ("   Equilibrated maximum displacement of chain / acceptance ratio:   \n");
        printEqStat(chainm,2.0,MAXMT);

        printf ("   Equilibrated maximum rotation of chain / acceptance ratio:              \n");
        printEqStat(chainr,1.0,MAXMT);
        printf ("\n");
    }



#ifdef ENABLE_MPI
    MPI_Datatype* defDataType(MPI_Datatype* MPI_stat) {

        StatsGrand g;
        Disp d;

        MPI_Datatype MPI_DISP, MPI_GRAND;

        d.defDataType(&MPI_DISP);
        g.defDataType(&MPI_GRAND);

        MPI_Aint     dispstart;

        MPI_Datatype mpiexdataType[7] = {MPI_GRAND, MPI_DISP, MPI_DISP, MPI_DISP, MPI_DISP, MPI_DISP, MPI_DISP};
        int          mpiexdataLen[7]  = {MAXMT,     1,        MAXT,     MAXT,     MAXMT,    MAXMT,    1};

        MPI_Aint     mpiexdata[7];
        MPI_Address( this, &dispstart);
        MPI_Address( &(this->grand), &mpiexdata[0]);

        MPI_Address( &(this->edge), &mpiexdata[1]);
        MPI_Address( &(this->rot), &mpiexdata[2]);
        MPI_Address( &(this->trans), &mpiexdata[3]);
        MPI_Address( &(this->chainm), &mpiexdata[4]);
        MPI_Address( &(this->chainr), &mpiexdata[5]);
        MPI_Address( &(this->mpiexch), &mpiexdata[6]);

        for (int i=0; i <7; i++)
            mpiexdata[i] -= dispstart;

        MPI_Type_struct(7, mpiexdataLen, mpiexdata, mpiexdataType, MPI_stat);
        MPI_Type_commit( MPI_stat);

        MPI_Type_free(&MPI_DISP);
        MPI_Type_free(&MPI_GRAND);

        return MPI_stat;
    }
#endif
};


#endif // STATISTICS_H
