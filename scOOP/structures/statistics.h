/** @file statistics.h*/

#ifndef STATISTICS_H
#define STATISTICS_H

#include <iomanip>

#include "macros.h"

class StatsGrand {
public:
    StatsGrand() : delAcc(0), delRej(0), insAcc(0), insRej(0),
    muVtAverageParticles(0), muVtSteps(0) {}

    int delAcc;
    int delRej;
    int insAcc;
    int insRej;

    unsigned long long int muVtAverageParticles;
    unsigned int muVtSteps;

#ifdef ENABLE_MPI
    MPI_Datatype* defDataType(MPI_Datatype* MPI_stat) {

        MPI_Aint     dispstart;

        MPI_Datatype mpiexdataType[6] = {MPI_INT,    MPI_INT,    MPI_INT,    MPI_INT,    MPI_LONG, MPI_INT};
        int          mpiexdataLen[6]  = {1,          1,          1,          1,          1,        1};

        MPI_Aint     mpiexdata[6];
        MPI_Address( this, &dispstart);
        MPI_Address( &(this->delAcc), &mpiexdata[0]);
        MPI_Address( &(this->delRej), &mpiexdata[1]);
        MPI_Address( &(this->insAcc), &mpiexdata[2]);

        MPI_Address( &(this->insRej), &mpiexdata[3]);
        MPI_Address( &(this->muVtAverageParticles), &mpiexdata[4]);
        MPI_Address( &(this->muVtSteps), &mpiexdata[5]);

        for (int i=0; i <6; i++)
            mpiexdata[i] -= dispstart;

        MPI_Type_struct(6, mpiexdataLen, mpiexdata, mpiexdataType, MPI_stat);
        MPI_Type_commit( MPI_stat);

        return MPI_stat;
    }
#endif
};

/**
  * @brief Define step size and acceptance ratio statistics
  **/
class Disp {
public:
    double mx;          ///< \brief Maximum value displacement, cos(angle), etc.
    double angle;       ///< \brief Maximum angle, since in .mx cos(angle) is saved
    double oldrmsd;     ///< \brief Averaged mx value in previous equilibration round
    double oldmx;       ///< \brief Change in mx in last equlibrium step
    long acc;           ///< \brief Number of accepted steps
    long rej;           ///< \brief Number of rejected steps

    Disp() : mx(0.0), angle(0.0), oldrmsd(0.0), oldmx(0.0), acc(0), rej(0) {}

inline double ratio() {
    if(acc + rej > 0) {
        return ((double) acc)/(acc+rej);
    }
    else return 0.0;
}

#ifdef ENABLE_MPI
    MPI_Datatype* defDataType(MPI_Datatype* MPI_stat) {

        MPI_Aint     dispstart;

        MPI_Datatype mpiexdataType[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_LONG};
        int          mpiexdataLen[6]  = {1,          1,          1,          1,          1,        1};

        MPI_Aint     mpiexdata[6];
        MPI_Address( this, &dispstart);
        MPI_Address( &(this->mx), &mpiexdata[0]);
        MPI_Address( &(this->angle), &mpiexdata[1]);
        MPI_Address( &(this->oldrmsd), &mpiexdata[2]);

        MPI_Address( &(this->oldmx), &mpiexdata[3]);
        MPI_Address( &(this->acc), &mpiexdata[4]);
        MPI_Address( &(this->rej), &mpiexdata[5]);

        for (int i=0; i <6; i++)
            mpiexdata[i] -= dispstart;

        MPI_Type_struct(6, mpiexdataLen, mpiexdata, mpiexdataType, MPI_stat);
        MPI_Type_commit( MPI_stat);

        return MPI_stat;
    }
#endif
};

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

    void print() {
        mcout.get() << std::setprecision(2) << std::left;
        mcout.get() << "\n\n******************************************************************************" << endl;
        mcout.get() << "*                               Moves Statistics                             *" << endl;
        mcout.get() << "******************************************************************************" << endl;
        mcout.get() << setw(30) << "Move" << setw(10) << "Acc (%)" << setw(10) << "Rej (%)" << setw(10) << "Steps" << endl;

        if(stepsSTrans() > 0) {
            mcout.get() << setw(30) << "Single particle translation: "
                        << setw(10) <<  (double)accSTrans()/stepsSTrans()*100.0
                        << setw(10) << (double)rejSTrans()/stepsSTrans()*100.0
                        << setw(10) << stepsSTrans() << endl;
        }

        if(stepsSRot() > 0) {
            mcout.get() << setw(30) << "Single particle rotation: "
                        << setw(10) << 100.0*accSRot()/stepsSRot()
                        << setw(10) << (double)rejSRot()/stepsSRot()*100.0
                        << setw(10) << stepsSRot() << endl;
        }

        if(stepsCTrans() > 0) {
            mcout.get() << setw(30) << "Chain translation: "
                        << setw(10) <<  (double)accCTrans()/stepsCTrans()*100.0
                        << setw(10) << (double)rejCTrans()/stepsCTrans()*100.0
                        << setw(10) << stepsCTrans() << endl;
        }

        if(stepsCRot() > 0) {
            mcout.get() << setw(30) << "Chain rotation: "
                        << setw(10) <<  (double)accCRot()/stepsCRot()*100.0
                        << setw(10) << (double)rejCRot()/stepsCRot()*100.0
                        << setw(10) << stepsCRot() << endl;
        }

        if(edge.acc + edge.rej > 0) {
            mcout.get() << setw(30) << "Pressure move: "
                        << setw(10) << (double)edge.acc / (edge.acc + edge.rej)*100.0
                        << setw(10) << (double)edge.rej / (edge.acc + edge.rej)*100.0
                        << setw(10) << (edge.acc + edge.rej) << endl;
        }

        if(stepsSwitch() > 0) {
            mcout.get() << setw(30) << "Switch type move: "
                        << setw(10) <<  (double)accSwitch()/stepsSwitch()*100.0
                        << setw(10) << (double)rejSwitch()/stepsSwitch()*100.0
                        << setw(10) << stepsSwitch() << endl;
        }

        for(int i=0; i<MAXT; ++i) {
            if(trans[i].acc + trans[i].rej > 0) {
                mcout.get() << setw(30) << "Single particle translation:"
                            << setw(10) << (double) trans[i].acc / (trans[i].acc + trans[i].rej) * 100.0
                            << setw(10) << (double) trans[i].rej / (trans[i].acc + trans[i].rej) * 100.0
                            << setw(10) << (trans[i].acc + trans[i].rej) << " type " << i << " name.get()" << endl;
            }
        }
    }

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
