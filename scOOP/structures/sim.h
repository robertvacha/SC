/** @file sim.h*/

#ifndef SIM_H
#define SIM_H

#include <tuple>
#include "structures.h"
#include "Conf.h"
#include "../mc/mygetline.h"
#include "../structures/statistics.h"

/**
 * @brief Should contain mostly all the simulation options and variables
 */
class Sim {
public:
    //
    //  Simulation constants (exchange in replica exchange moves)
    //
    double press;               ///< \brief Pressure
    double temper;              ///< \brief Temperature
    int pseudoRank;             ///< \brief For MPI i/o, changes depending on temp

    //
    //  Statistics (exchange in replica exchange moves)
    //
    Statistics stat;

    //
    //  Read only variables NOTE: should change to constzz
    //
    double paralpress;          ///< \brief Parallel pressure for replica exachnge
    double dpress;		        ///< \brief Pressure change for replica exchange
    double shave;               ///< \brief Average number of volume changes to attempt per sweep
    double shprob;              ///< \brief Probability of attempting a volume change
    double chainprob;           ///< \brief Average number of chain move attempt per sweep
    double switchprob;          ///< \brief Average number of type switch attempt per sweep
    int pairlist_update;        ///< \brief Number of sweep per updating the pairlist
    double paraltemper;         ///< \brief Temperature for parallel tempering
    double dtemp;               ///< \brief Temprature step
    vector<double> pTemp;       ///< \brief Exact temperatures for paralel tempering
    int ptype;                  ///< \brief Type of pressure coupling
    long adjust;                ///< \brief Number of sweeps between step size adjustments
    long movie;                 ///< \brief Number of sweeps between movie frames
    long nequil;                ///< \brief Number of equilibration sweeps
    long nsweeps;               ///< \brief Number of production sweeps
    long paramfrq;              ///< \brief Number of sweeps between order parameter samples
    long report;                ///< \brief Number of sweeps between statistics reports
    //long terms;                 ///< \brief Number of Fourier terms as smectic order parameters
    long nrepchange;            ///< \brief Number of sweeps between replica exchanges
    long nGrandCanon;           ///< \brief Number of sweeps between particle insert/delete
    long nClustMove;            ///< \brief Number of sweeps between cluster moves
    double coneAngle;             ///< \brief Prephere rotation around axis of spherocylinder in particular angle from axis
    long write_cluster;         ///< \brief Number of sweeps per writing out cluster info


    int mpirank;                ///< \brief MPI number for given process, identical to calling MPI_Comm_rank, constant during simulation
    int mpinprocs;              ///< \brief MPI number of processes

    double cell;                  ///< \brief Maximum translation of all types
    double max_dist_squared[MAXT][MAXT]; ///< \brief Stored cutoffs of all particle types for pairList

    int wlm[2];                ///< \brief Wang landau method (wl)
    int wlmtype;            ///< \brief Atom type for the Wang landau method (wl)

    size_t pairList;
    //size_t energyCalc;
    //size_t move;
    size_t all;

    Sim(Conf* conf, FileNames* files, int rank=0, int procs=1): press(0.0), temper(0.0), pseudoRank(rank), paralpress(0.0), dpress(0.0), shave(0.0), shprob(0.0), chainprob(0.0), switchprob(0.0), pairlist_update(0),
        paraltemper(0.0), dtemp(0.0), ptype(0), adjust(0), movie(0), nequil(0), nsweeps(0),paramfrq(0), report(0),
        nrepchange(0), nGrandCanon(0), nClustMove(0), coneAngle(0.0), mpirank(rank), mpinprocs(procs), cell(0.0), pairList(0), /*energyCalc(0), move(0),*/ all(0) {

        readOptions(files);
    }

    void info() {
        mcout.get() << " Pressure coupling type:                             " << ptype << endl;
        mcout.get() << " Pressure:                                           " << press << endl;
        mcout.get() << " Replica exchange pressure:                          " << paralpress << endl;
        mcout.get() << " Average volume change attempts per sweep:           " << shave << endl;
        mcout.get() << " Equilibration sweeps:                               " << nequil << endl;
        mcout.get() << " Sweeps between step size adjustments:               " << adjust << endl;
        mcout.get() << " Production sweeps:                                  " << nsweeps << endl;
        mcout.get() << " Sweeps between statistics samples:                  " << paramfrq << endl;
        mcout.get() << " Sweeps between statistics reports:                  " << report << endl;
        mcout.get() << " Average chain move attempts per sweep:              " << chainprob << endl;
        mcout.get() << " Inititial maximum angular cone angle (degrees):     " << coneAngle << endl;
        mcout.get() << " Inititial maximum geo.box edge change:              " << stat.edge.mx << endl;
        mcout.get() << " Temperature in kT/e:                                " << temper << endl;
        mcout.get() << " Parallel tempering temperature in kT/e:             " << paraltemper << endl;
        mcout.get() << " Sweeps between replica exchange:                    " << nrepchange << endl;
        mcout.get() << " Sweeps between Grand-Canonical move:                " << nGrandCanon << endl;
        mcout.get() << " Sweeps between Cluster moves:                       " << nClustMove << endl;
        mcout.get() << " Wang-Landau method:                                 " << wlm[0] << " "  << wlm[1] << endl;
        mcout.get() << " Calculate the Wang-Landau method for atom type:     " << wlmtype << endl;
        mcout.get() << " Average type switch attempts per sweep:             " << switchprob << endl;
        mcout.get() << " Number of Sweeps per pairlist update:               " << pairlist_update << endl;
        mcout.get() << " Number of sweeps per writing out cluster info:      " << write_cluster << endl;
    }

private:
    int deallocPairlist() {  // deprecated, done in memoryDealloc()
        /*long i;
        if(sim.pairlist != NULL){
            for(i = 0; i < topo.npart; i++){
                if(sim.pairlist[i].pairs != NULL){
                    free(sim.pairlist[i].pairs);
                }
            }
            free(sim.pairlist);
        }*/
        return 0;
    }


    /**
     * @brief Reads the run parameters from the external file "options".  See the end of the
       code for a template. All comments starting with '#' are stripped out.  The
       options are summarised on standard output and checked for validity of range.
     */
    void readOptions(FileNames* files) {

        mcout.get() << "Reading options..." << endl;

        int num_options = -1;
        double transmx, rotmx, chainmmx, chainrmx, angle, chain_angle;
        std::vector<tuple<int,double>> transmx_type;
        long int seed;

        char *id, *value, *num, *tokLine, *line;
        FILE *infile;

        // for new options add before the last line
        Option options[] = {
            {"write_cluster",       Long,   false, &write_cluster},
            {"adjust",              Long,   false, &adjust},
            {"movie",               Long,   false, &movie},
            {"nequil",              Long,   false, &nequil},
            {"nsweeps",             Long,   false, &nsweeps},
            {"nrepchange",          Long,   false, &nrepchange},
            {"nGrandCanon",         Long,   false, &nGrandCanon},
            {"nClustMove",          Long,   false, &nClustMove},
            {"paramfrq",            Long,   false, &paramfrq},
            {"report",              Long,   false, &report},
            {"seed",                Long,   false, &seed},
            {"pairlist_update",     Int,    false, &pairlist_update},
            {"ptype",               Int,    false, &ptype},
            {"wlm",                 Int2,   false, &wlm},
            {"wlmtype",             Int,    false, &wlmtype},
            {"press",               Double, false, &press},
            {"paralpress",          Double, false, &paralpress},
            {"edge_mx",             Double, false, &stat.edge.mx},
            {"shave",               Double, false, &shave},
            {"chainprob",           Double, false, &chainprob},
            {"switchprob",          Double, false, &switchprob},
            {"temper",              Double, false, &temper},
            {"paraltemper",         Double, false, &paraltemper},
            {"transmx",             Double, false, &transmx},
            {"transmx_type",        Tuple, true, &transmx_type}, // default se by transmx
            {"rotmx",               Double, false, &rotmx},
            {"coneAngle",           Double, true, &coneAngle}, // default value given in constructor of sim
            {"chainmmx",            Double, false, &chainmmx},
            {"chainrmx",            Double, false, &chainrmx},
            {"last",                Int,    false, NULL}
        };
        while(options[++num_options].var != NULL)
            ;

        //--- 1. Read in values ---
        size_t line_size = (STRLEN + 1) * sizeof(char);
        line = (char *) malloc(line_size);

        infile = fopen(files->optionsfile, "r");
        if (infile == NULL) {
            fprintf (stderr, "\nERROR: Could not open options file.\n\n");
            exit (1);
        }

        while(myGetLine(&line, &line_size, infile) != -1){

            // strip comments
            strip_comment(line);
            trim(line);
            if(strlen(line) == 0) continue;

            // tokenize
            tokLine = line;
            id = strtok(tokLine, "=");
            if(id == NULL){
                fprintf(stderr, "error parsing Configuration line (%s)", line);
                free(line);
                exit(1);
            }
            trim(id);
            tokLine = NULL;
            value = strtok(tokLine, "=");
            trim(value);
            if(value == NULL){
                fprintf(stderr, "error parsing Configuration line (%s)", line);
                free(line);
                exit(1);
            }
            //printf("id: %s; value: %s\n", id, value);
            int i = 0;
            for(i=0; i < num_options; i++){
                if(strcmp(id, options[i].id) == 0){
            if(options[i].type == Int2){
                readii2(value,*((int (*)[2]) options[i].var));
                        options[i].set = true;
                        break;
                    }
                    if(options[i].type == Int) {
                        *((int *) options[i].var) = readi2(value);
                        options[i].set = true;
                        break;
                    }
                    else if(options[i].type == Long) {
                        *((long *) options[i].var) = readl2(value);
                        options[i].set = true;
                        break;
                    }
                    else if(options[i].type == Double) {
                        *((double *) options[i].var) = readd2(value);
                        options[i].set = true;
                        break;
                    }
                    else if(options[i].type == Tuple) {
                        int val_int=0;
                        double val_double = 0.0;

                        num = strtok(value, " "); // we must start with value, but then only use NULL... stupid stupid stupid
                        val_int = readl2(num);

                        num = strtok(NULL, " ");
                        val_double = readd2(num);

                        transmx_type.push_back( tuple<int,double>( val_int, val_double ) );

                        while(num != NULL) {
                            num = strtok(NULL, " ");
                            if(num == NULL) // ugly, but meh...
                                break;
                            val_int = readl2(num);

                            num = strtok(NULL, " ");
                            val_double = readd2(num);

                            transmx_type.push_back( tuple<int,double>( val_int, val_double ) );
                        }

                        options[i].set = true;
                        break;
                    }
                    else {
                        fprintf(stderr, "Could not determine type of %s!\n", id);
                        free(line);
                        exit(1);
                    }
                }
            }
            if(i == num_options){
                fprintf(stderr, "Unknown identifier %s!\nWill procede.\n", id);
            }
        }
        fclose (infile);
        free(line);

        // Check, wheter all options have been readin
        for(int i = 0; i < num_options; i++){
            if(!options[i].set){
                fprintf(stderr, "option '%s' is not set!\n", options[i].id);
                exit(1);
            }
        }

        //--- 2. Summarize results on standard output ---
        // Density of close-packed spherocylinders
        //   rho_cp = 2.0/(sqrt(2.0) + *length * sqrt(3.0));

        mcout.get() << " Pressure coupling type:                             " << ptype << endl;
        mcout.get() << " Pressure:                                           " << press << endl;
        mcout.get() << " Replica exchange pressure:                          " << paralpress << endl;
        mcout.get() << " Average volume change attempts per sweep:           " << shave << endl;
        mcout.get() << " Equilibration sweeps:                               " << nequil << endl;
        mcout.get() << " Sweeps between step size adjustments:               " << adjust << endl;
        mcout.get() << " Production sweeps:                                  " << nsweeps << endl;
        mcout.get() << " Sweeps between statistics samples:                  " << paramfrq << endl;
        mcout.get() << " Sweeps between statistics reports:                  " << report << endl;
        mcout.get() << " Average chain move attempts per sweep:              " << chainprob << endl;
        mcout.get() << " Initial maximum displacement:                       " << transmx << endl;
        mcout.get() << " Inititial maximum angular change (degrees):         " << rotmx << endl;
        mcout.get() << " Inititial maximum angular cone angle (degrees):     " << coneAngle << endl;
        mcout.get() << " Inititial maximum geo.box edge change:              " << stat.edge.mx << endl;
        mcout.get() << " Initial maximum chain displacement:                 " << chainmmx << endl;
        mcout.get() << " Inititial maximum chain angular change (degrees):   " << chainrmx << endl;
        mcout.get() << " Temperature in kT/e:                                " << temper << endl;
        mcout.get() << " Parallel tempering temperature in kT/e:             " << paraltemper << endl;
        mcout.get() << " Sweeps between replica exchange:                    " << nrepchange << endl;
        mcout.get() << " Sweeps between Grand-Canonical move:                " << nGrandCanon << endl;
        mcout.get() << " Sweeps between Cluster moves:                       " << nClustMove << endl;
        mcout.get() << " Wang-Landau method:                                 " << wlm[0] << " "  << wlm[1] << endl;
        mcout.get() << " Calculate the Wang-Landau method for atom type:     " << wlmtype << endl;
        mcout.get() << " Average type switch attempts per sweep:             " << switchprob << endl;
        mcout.get() << " Number of Sweeps per pairlist update:               " << pairlist_update << endl;
        mcout.get() << " Random number seed:                                 " << seed << endl;
        mcout.get() << " Number of sweeps per writing out cluster info:      " << write_cluster << endl;

        if (movie > 0) {
            mcout.get() << " Sweeps between movie frames:                      " << movie << endl;
        } else {
            mcout.get() <<" No movie" << endl;
        }
        mcout.get() << endl;

        if(pairlist_update){
            mcout.get() << " A pairlist will be generated every " << pairlist_update << " steps. This is a greedy"
                                                                                        " algorithm; make sure you don't have big chains etc.!\n" << endl;
        }

        //--- 3. Validity checks ---
        if (rotmx < 0.0 || rotmx > 180) {
            fprintf (stderr, "ERROR: Maximum orientation change must be in range 0 to 180.\n\n");
            exit (1);
        }
        if (chainrmx < 0.0 || chainrmx > 180) {
            fprintf (stderr, "ERROR: Maximum orientation change for chains must be in range 0 to 180.\n\n");
            exit (1);
        }
        if ( (ptype <0) || (ptype>5) ) {
            fprintf (stderr, "ERROR: Unknown pressure coupling %d. Program only knows: 0 - anisotropic coupling, \
                    1 - isotropic coupling, 2 - isotropic in xy and z=const, 3 - isotropic xy and  V=const, 4 = anisotropic in xy and z=const; 5 - anisotropic in y and xz=const.\n\n",ptype);
            exit (1);
        }
        if ( (wlm[0] <0) || (wlm[0] > 9) || (wlm[1] <0) || (wlm[1] > 9)  ) {
            fprintf (stderr, "ERROR: Unknown Wang-Landau method %d %d. Program only knows: 0 - none, \
                    1 - z-direction od 1st particle, 2 - pore in membrane, 3 - zorientation of 0th particle,\
                    4 - distance of fist two particles, 5 - pore around z-axis above CM,\
                    6 - pore around z-axis above 0th particle, 7 - number of particles in contact,\
                    8 - sampling of box length in x, 9 - sampling of box length in y \n\n",wlm[0],wlm[1]);
            exit (1);
        }
        if ( (wlm[0] == 0) && (wlm[1] > 0)  ) {
            fprintf (stderr, "ERROR: Wang-Landau method has to be set for first order parameter and then for second order parameter\n\n");
            exit (1);
        }
        if ( (wlm[0] == 2) || (wlm[0] == 5) || (wlm[0] == 6)  ) {
            if(wlmtype < 1){
                fprintf (stderr, "ERROR: Atom type for the Wang-Landau Method (%d) was false defined.\n\n",wlmtype);
                exit (1);
            }
            if ( (wlm[1] == 2) || (wlm[1] == 5) || (wlm[1] == 6) ) {
                fprintf (stderr, "ERROR: Simulaneous use of two pore order parameters has not been implemented yet.\n\n");
                exit (1);
            }
        }

        // we store maximum rotation as half angle - useful for quaterions
        angle = rotmx / 180.0 * PIH *0.5;
        rotmx = cos((rotmx)/180.0*PIH);
        chain_angle = chainrmx / 180.0 * PIH;
        chainrmx = cos((chainrmx)/180.0*PIH);
        stat.edge.mx *= 2.0;   // The full range is -maxl to +maxl, i.e. spanning 2*maxl
        transmx *= 2.0;   // The full range is -maxr to +maxr, i.e. spanning 2*maxr
        chainmmx *= 2.0;   // The full range is -maxr to +maxr, i.e. spanning 2*maxr
        coneAngle *= DEGTORAD; // Now transfer angle in degrees into radians

        for (int i=0;i<MAXT;i++) {
            stat.trans[i].mx = transmx;
            stat.rot[i].mx = rotmx;
            stat.rot[i].angle = angle;
        }
        for(unsigned int i = 0; i<transmx_type.size(); ++i) {
            stat.trans[ std::get<0>(transmx_type[i]) ].mx = std::get<1>(transmx_type[i]);
        }
        for (int i=0;i<MAXMT;i++) {
            stat.chainm[i].mx = chainmmx;
            stat.chainr[i].mx = chainrmx;
            stat.chainr[i].angle = chain_angle;
        }

        //parallel tempering
#ifdef ENABLE_MPI
        if ( (temper != paraltemper) && (mpinprocs <2) ) {
            cerr << "ERROR: Paralllel tempering at single core does not work.\n" << endl;
            exit(1);
        }
        dtemp = ((1.0/temper)-(1.0/paraltemper))/(mpinprocs-1);
        for(int i=0; i<mpinprocs; i++) {
            pTemp.push_back( temper/(1.0-i*temper*dtemp) );
        }
        temper =  temper/(1.0-mpirank*temper*dtemp);
        if ( (press != paralpress) && (mpinprocs <2) ) {
            cerr << "ERROR: Pressure replica exchange at single core does not work.\n" << endl;
            exit(1);
        }
        dpress = (paralpress - press )/(mpinprocs-1);
        press += dpress * mpirank;
        seed += mpirank;
        stat.mpiexch.mx = dtemp;
        stat.mpiexch.angle = dpress;
#endif

        ran2.setSeed(seed);
    }

    /**
     * @brief convert string num into two integers
     * @param num
     * @param value
     */
    void readii2(char * num, int value[2]) {
        char *end, *num2;

        value[0] = strtol(num, &num2, 10);
        trim(num2);
        if ((int)strlen(num2) > 0)
            value[1] = strtol(num2, &end, 10);
        else {
            value[1] =0;
            return;
        }
        if(*end){
            fprintf(stderr, "Could not convert %s into two integers\n", num);
            exit(1);
        }

        return;
    }

    /**
     * @brief convert string num into double
     * @param num
     * @return
     */
    double readd2(char * num) {
        char *end;
        double i = strtod(num, &end);
        if(*end){
            fprintf(stderr, "Could not convert %s into double\n", num);
            exit(1);
        }
        return i;
    }

    /**
     * @brief convert string num into long
     * @param num
     * @return
     */
    long readl2(char * num) {
        char *end;
        long i = strtol(num, &end, 10);
        if(*end){
            fprintf(stderr, "Could not convert %s into long\n", num);
            exit(1);
        }
        return i;
    }

    /**
     * @brief convert string num into integer
     * @param num
     * @return
     */
    int readi2(char * num) {
        char *end;
        int i = strtol(num, &end, 10);
        if(*end){
            fprintf(stderr, "Could not convert %s into integer\n", num);
            exit(1);
        }
        return (int) i;
    }
};

#endif // SIM_H
