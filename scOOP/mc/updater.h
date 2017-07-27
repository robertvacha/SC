/** @file updater.h*/

#ifndef UPDATER_H
#define UPDATER_H

#include "movecreator.h"
#include "../structures/Conf.h"
#include "clust.h"

class Updater
{
public:
    Updater(Sim* sim, Conf* conf, FileNames* files) :
        sim(sim), conf(conf), files(files),
        calcEnergy(sim, conf), move(sim, conf, &calcEnergy), clust(conf, sim, &calcEnergy, files) {
            if ( sim->wlm[0] > 0 ) {
            FILE* outfile = fopen(files->wlinfile, "r");
            if (outfile == NULL) {
                printf ("ERROR: Cannot open file for Wang-Landau method (%s).\n",files->wlinfile);
                exit(1);
            }
            fclose (outfile);
            }
        }

private:
    Sim* sim;            ///< \brief contains the simulation options.
    Conf* conf;          ///< \brief contains particle vector and geometry
    FileNames* files;

    TotalEnergyCalculator calcEnergy;   ///< \brief energy calculations
    MoveCreator move;                   ///< \brief move calculations

    ClusterSampler clust;

    long nsweeps;
    long adjust;
    long paramfrq;
    long report;

    long next_frame;   // Next sweep number for dumping a movie fram
    long next_adjust;  // Next sweep number for step size adjustment
    long next_calc;    // Next sweep number for order parameter calculation
    long next_dump;    // Next sweep number for reporting statistics

public:
    /**
     * @brief simulate
     * @param nsweeps
     * @param adjust
     * @param paramfrq
     * @param report
     */
    void simulate(long nsweeps, long adjust, long paramfrq, long report);

private:
    /**
     * @brief testEnergyMatrix
     * @return TRUE - Energy matrix is fine
     */
    bool testEnergyMatrix() {
        if(conf->pvec.empty())
            return true;

        for (unsigned int i = 1; i < conf->pvec.size(); ++i) {
            for (unsigned long j = 0; j < i; ++j) {
                if( !(calcEnergy.eMat.energyMatrix->operator [](i)[j] + 0.0000001 >= calcEnergy.p2p(i,j)
                        && calcEnergy.eMat.energyMatrix->operator [](i)[j] - 0.0000001 <= calcEnergy.p2p(i,j)  ) ) {
                    cout << "[" << i << "][" << j << "]= " << calcEnergy.eMat.energyMatrix->operator [](i)[j] << ", "
                         << ", calc= " << calcEnergy.p2p(i,j) << endl;
                    return false;
                }
            }
        }

        return true;
    }

    void dumpMovie(long sweep) {
        FILE* mf;
        if (sim->movie > 0) {
            mf = fopen(files->moviefile, "a");
            //fprintf (mf, "> box %.10f %.10f %.10f ; num_part %ld ; sweep %ld <\n", conf->geo.box.x, conf->geo.box.y, conf->geo.box.z, (long)conf->pvec.size(), sweep);
            fprintf (mf, "%ld\nsweep %ld; box %.10f %.10f %.10f\n",(long)conf->pvec.size(), sweep, conf->geo.box.x, conf->geo.box.y, conf->geo.box.z);
            conf->draw(mf);
            fflush (mf);
            next_frame += sim->movie;
            fclose (mf);
        }
    }

    void emptyFiles();
    void initValues();

    /**
     * @brief optimizestep Optimize the maximum displacement within the specified limits and resets the
       acceptance counters to zero.
     * @param x
     * @param hi
     * @param lo
     */
    void optimizeStep(Disp *x, double hi, double lo);

    /**
     * @brief optimizerot Optimize the maximum rotation within the specified limits and resets the
       acceptance counters to zero. Rotation is given by cos of angle
       larger rotation = smaller cos
     * @param x
     * @param hi
     * @param lo
     */
    void optimizeRot(Disp *x, double hi, double lo);

    /**
     * @brief alignment_order alignment ORDER PARAMETER
     * @return
     */
    double alignmentOrder();

    /**
     * @brief gen_pairlist Interface for the generation of the pairlist. Define other pairlist algorithms above.
     */
    inline void genPairList() {
        genSimplePairList();
    }

    /**
     * @brief gen_simple_pairlist Generates a pairlist with a very basic alogrithm
     */
    void genSimplePairList();
};

#endif // UPDATER_H
