/** @file statistics.h*/

#ifndef STATISTICS_H
#define STATISTICS_H

#include <iomanip>
#include <array>
#include "macros.h"
#include "topo.h"

extern Topo topo;

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
    array<Disp, MAXT> rot;      ///< \brief Maximum rotation and statistics
    array<Disp, MAXT> trans;    ///< \brief Maximum translation  and statistics
    array<Disp, MAXT> switchMv; ///< \brief Switch move statistics
    array<Disp, MAXMT> chainm;  ///< \brief Maximum translation for chain  and statistics
    array<Disp, MAXMT> chainr;  ///< \brief Maximum rotation for chain and statistics
    Disp mpiexch;               ///< \brief MPI statistics

    void print() {
        mcout.get() << std::setprecision(2) << std::left;
        mcout.get() << "\n\n******************************************************************************" << endl;
        mcout.get() << "*                               Moves Statistics                             *" << endl;
        mcout.get() << "******************************************************************************" << endl;
        mcout.get() << setw(30) << "Move" << setw(10) << "Acc (%)" << setw(10) << "Rej (%)" << setw(10) << "Steps" << endl;

        if(steps<MAXT>(trans) > 0)
            mcout.get() << moveInfo<MAXT>(trans, "Single particle translation: ") << endl;

        if(steps<MAXT>(rot) > 0)
            mcout.get() << moveInfo<MAXT>(rot, "Single particle rotation: ") << endl;

        if(steps<MAXMT>(chainm) > 0)
            mcout.get() << moveInfo<MAXMT>(chainm, "Chain translation: ") << endl;

        if(steps<MAXMT>(chainr) > 0)
            mcout.get() << moveInfo<MAXMT>(chainr, "Chain rotation: ") << endl;

        if(edge.acc + edge.rej > 0)
            mcout.get() << moveInfo(edge, "Pressure move: ") << endl;

        if(steps<MAXT>(switchMv) > 0)
            mcout.get() << moveInfo<MAXT>(switchMv, "Switch type move: ") << endl;

        for(int i=0; i<MAXT; ++i) {
            if(trans[i].acc + trans[i].rej > 0) {
                mcout.get() << moveInfo(trans[i], "Single particle translation: ") << " type " << i << topo.ia_params[i][i].name << endl;
            }
        }
    }

    string moveInfo(Disp& disp, string name) {
        stringstream ss;

        ss << std::setprecision(3) << std::left;
        ss << setw(30) << name
           << setw(10) << (double)disp.acc / (disp.acc + disp.rej) * 100.0
           << setw(10) << (double)disp.rej / (disp.acc + disp.rej) * 100.0
           << setw(10) << (disp.acc + disp.rej);

        return ss.str();
    }

    template<std::size_t SIZE>
    string moveInfo(std::array<Disp, SIZE>& disp, string name) {
        stringstream ss;

        ss << std::setprecision(3) << std::left;
        ss << setw(30) << name
           << setw(10) <<  (double)acc<SIZE>(disp)/steps<SIZE>(disp)*100.0
           << setw(10) << (double)rej<SIZE>(disp)/steps<SIZE>(disp)*100.0
           << setw(10) << steps<SIZE>(disp);

        return ss.str();
    }

    template<std::size_t SIZE>
    int acc(std::array<Disp, SIZE>& disp) {
        int var=0;
        for(unsigned int i=0; i<SIZE; ++i)
            var += disp[i].acc;
        return var;
    }

    template<std::size_t SIZE>
    int rej(std::array<Disp, SIZE>& disp) {
        int var=0;
        for(unsigned int i=0; i<SIZE; ++i)
            var += disp[i].rej;
        return var;
    }

    template<std::size_t SIZE>
    int steps(std::array<Disp, SIZE>& disp) {
        return acc<SIZE>(disp) + rej<SIZE>(disp);
    }

    template<std::size_t SIZE>
    void printEqStat(array<Disp, SIZE> dat, double scale) {
        for(unsigned int i=0; i<SIZE; ++i) {
            if (RATIO(dat[i]) > 0)
                printf ("   TYPE %d           %.6f  /  %.6f\n", i, dat[i].mx/scale,RATIO(dat[i]));
        }
    }


    void printEqStat() {

        printf ("   Equilibrated maximum displacement / acceptance ratio:            \n");
        printEqStat<MAXT>(trans,2.0);

        printf ("   Equilibrated maximum rotation / acceptance ratio:                       \n");
        printEqStat<MAXT>(rot,1.0);

        printf ("   Equilibrated maximum box length change / acceptance ratio:              \n");
        printf ("                     %.6e  /  %.6e\n", edge.mx/2.0,RATIO(edge));

        printf ("   Equilibrated maximum displacement of chain / acceptance ratio:   \n");
        printEqStat<MAXMT>(chainm,2.0);

        printf ("   Equilibrated maximum rotation of chain / acceptance ratio:              \n");
        printEqStat<MAXMT>(chainr,1.0);
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
