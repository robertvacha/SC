/** @file sim.h*/

#ifndef SIM_H
#define SIM_H

#include <tuple>
#include "Conf.h"
#include "../mc/mygetline.h"
#include "../structures/statistics.h"


/**
 * @brief Holds the type of a variable in struct option
 */
typedef enum {
    Int,
    Int2,
    Long,
    Double,
    Tuple
} Type;

typedef struct {        // for reading in the options
    const char *id;     // The name of the value in the option file
    Type type;          // The type (int, double or long)
    bool set;           // Whether the variable has been set
    void *var;          // The variable
} Option;

/**
 * @brief Should contain mostly all the simulation options and variables
 */
class Sim {
public:
    ///
    ///  Equilibration setting
    ///
    long nequil = 0;                ///< \brief Number of equilibration sweeps
    long adjust = 0;                ///< \brief Number of sweeps between step size adjustments

    ///
    /// Production setting
    ///
    long nsweeps = 0;               ///< \brief Number of production sweeps
    double temper = 0.0;              ///< \brief Temperature
    int pairlist_update = 0;        ///< \brief Number of sweep per updating the pairlist
    long int seed = 0;

    ///
    /// Barostat setting
    ///
    int ptype = 0;                  ///< \brief Type of pressure coupling
    double press = 0.0;               ///< \brief Pressure
    double shave = 0.0;               ///< \brief Average number of volume changes to attempt per sweep
    double shprob = 0.0;              ///< \brief Probability of attempting a volume change

    ///
    /// Wang-Landau method
    ///
    int wlm[2];                ///< \brief Wang landau method (wl)
    int wlmtype;            ///< \brief Atom type for the Wang landau method (wl)

    ///
    /// Parallel Tempering setting
    ///
    long nrepchange = 0;            ///< \brief Number of sweeps between replica exchanges
    double paraltemper = 0.0;         ///< \brief Temperature for parallel tempering
    double paralpress = 0.0;          ///< \brief Parallel pressure for replica exachnge
    double dtemp = 0.0;               ///< \brief Temprature step
    double dpress = 0.0;		      ///< \brief Pressure change for replica exchange
    vector<double> pTemp;       ///< \brief Exact temperatures for paralel tempering


    ///
    /// GrandCanonical setting
    ///
    long nGrandCanon = 0;           ///< \brief Number of sweeps between particle insert/delete

    ///
    /// Simulation settings
    ///
    long nClustMove = 0;            ///< \brief Number of sweeps between cluster moves
    double switchprob = 0.0;          ///< \brief Average number of type switch attempt per sweep
    double chainprob = 0.0;           ///< \brief Average number of chain move attempt per sweep
    double transmx = 0.0;
    std::vector<tuple<int,double>> transmx_type;
    double rotmx = 0.0;
    double coneAngle = 0.0;             ///< \brief Prephere rotation around axis of spherocylinder in particular angle from axis
    double edge_mx = 0.0;
    double chainmmx = 0.0;
    double chainrmx = 0.0;

    ///
    /// Output setting
    ///
    int thermo = 0;
    long paramfrq = 0;              ///< \brief Number of sweeps between order parameter samples
    long report = 0;                ///< \brief Number of sweeps between statistics reports
    long movie = 0;                 ///< \brief Number of sweeps between movie frames
    long write_cluster;         ///< \brief Number of sweeps per writing out cluster info

    //
    //  Statistics (exchange in replica exchange moves)
    //
    Statistics stat;

    int pseudoRank;             ///< \brief For MPI i/o, changes depending on temp
    int mpirank;                ///< \brief MPI number for given process, identical to calling MPI_Comm_rank, constant during simulation
    int mpinprocs;              ///< \brief MPI number of processes

    double cell = 0.0;                  ///< \brief Maximum translation of all types
    double max_dist_squared[MAXT][MAXT]; ///< \brief Stored cutoffs of all particle types for pairList
    double customCutoff = 0.0;

    size_t pairList = 0;
    size_t all = 0;

    Sim(string optionsfile, int rank=0, int procs=1): pseudoRank(rank), mpirank(rank), mpinprocs(procs) {
        readOptions(optionsfile);
    }

    string toString() {
        std::ostringstream o;

        o << "###  Equilibration setting" << endl;
        o << " Equilibration sweeps:                               " << nequil << endl;
        o << " Sweeps between step size adjustments:               " << adjust << endl;

        o << "### Production setting" << endl;
        o << " Production sweeps:                                  " << nsweeps << endl;
        o << " Temperature in kT/e:                                " << temper << endl;
        o << " Number of Sweeps per pairlist update:               " << pairlist_update << endl;
        o << " Random number seed:                                 " << seed << endl;

        o << "### Barostat setting" << endl;
        o << " Pressure coupling type:                             " << ptype << endl;
        o << " Pressure:                                           " << press << endl;
        o << " Average volume change attempts per sweep:           " << shave << endl;

        o << "### Wang-Landau method" << endl;
        o << " Wang-Landau method:                                 " << wlm[0] << " "  << wlm[1] << endl;
        o << " Calculate the Wang-Landau method for atom type:     " << wlmtype << endl;

        o << "### Parallel Tempering setting" << endl;
        o << " Sweeps between replica exchange:                    " << nrepchange << endl;
        o << " Parallel tempering temperature in kT/e:             " << paraltemper << endl;
        o << " Replica exchange pressure:                          " << paralpress << endl;

        o << "### GrandCanonical setting" << endl;
        o << " Sweeps between Grand-Canonical move:                " << nGrandCanon << endl;

        o << "### Simulation settings" << endl;
        o << " Sweeps between Cluster moves:                       " << nClustMove << endl;
        o << " Average type switch attempts per sweep:             " << switchprob << endl;
        o << " Average chain move attempts per sweep:              " << chainprob << endl;
        o << " Initial maximum displacement:                       " << transmx << endl;
        o << " Initial maximum displacement type specific:         ";
        for(auto& item : transmx_type)
            o << std::get<0>(item) << " " << std::get<1>(item);
        o << endl;
        o << " Initial maximum orientation change:                 " << rotmx << endl;
        o << " Inititial maximum angular cone angle (degrees):     " << coneAngle << endl;
        o << " Inititial maximum geo.box edge change:              " << stat.edge.mx << endl;
        o << " Initial maximum chain displacement:                 " << chainmmx << endl;
        o << " Initial maximum chain rotation change:              " << chainrmx << endl;
        if(customCutoff != 0.0)
            o << " Custom cutoff:                                      " << customCutoff << endl;

        o << "### Output setting" << endl;
        o << " Sweeps between between reports to stdout:           " << thermo << endl;
        o << " Sweeps between statistics samples:                  " << paramfrq << endl;
        o << " Sweeps between statistics reports:                  " << report << endl;
        o << " Number of sweeps between movie frames:              " << movie << endl;
        o << " Number of sweeps per writing out cluster info:      " << write_cluster << endl;

        return o.str();
    }

private:
    /**
     * @brief Reads the run parameters from the external file "options".  See the end of the
       code for a template. All comments starting with '#' are stripped out.  The
       options are summarised on standard output and checked for validity of range.
     */
    void readOptions(string optionsfile) {

        mcout.get() << "Reading options..." << endl;

        int num_options = -1;
        double angle, chain_angle;

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
            {"edge_mx",             Double, false, &edge_mx},
            {"shave",               Double, false, &shave},
            {"chainprob",           Double, false, &chainprob},
            {"switchprob",          Double, false, &switchprob},
            {"temper",              Double, false, &temper},
            {"paraltemper",         Double, false, &paraltemper},
            {"transmx",             Double, false, &transmx},
            {"transmx_type",        Tuple,  true,  &transmx_type}, // default set by transmx
            {"rotmx",               Double, false, &rotmx},
            {"coneAngle",           Double, true,  &coneAngle}, // default value given in constructor of sim
            {"chainmmx",            Double, false, &chainmmx},
            {"chainrmx",            Double, false, &chainrmx},
            {"thermo",              Int,    true,  &thermo}, // default 0
            {"custom_cut",          Double,    true, &customCutoff}, // default 0.0
            {"last",                Int,    false, NULL}
        };
        while(options[++num_options].var != NULL)
            ;

        //--- 1. Read in values ---
        size_t line_size = (STRLEN + 1) * sizeof(char);
        line = (char *) malloc(line_size);

        infile = fopen( optionsfile.c_str(), "r");
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
        mcout.get() << toString() << endl;

        if(pairlist_update){
            mcout.get() << " A pairlist will be generated every " << pairlist_update << " steps. This is a greedy"
                                                                                        " algorithm; make sure you don't have big chains etc.!\n" << endl;
        }

        //--- 3. Validity checks ---
        checkValidity();

        stat.edge.mx = 2.0 * edge_mx;   // The full range is -maxl to +maxl, i.e. spanning 2*maxl
        coneAngle *= DEGTORAD; // Now transfer angle in degrees into radians
        for (int i=0;i<MAXT;i++) {
            stat.trans[i].mx = 2.0*transmx; // The full range is -maxr to +maxr, i.e. spanning 2*maxr
            stat.rot[i].mx = cos( rotmx * DEGTORAD * 0.5 );
            stat.rot[i].angle = rotmx * DEGTORAD * 0.5 * 0.5;
        }
        for(unsigned int i = 0; i<transmx_type.size(); ++i) {
            stat.trans[ std::get<0>(transmx_type[i]) ].mx = std::get<1>(transmx_type[i]);
        }
        for (int i=0;i<MAXMT;i++) {
            stat.chainm[i].mx = 2.0*chainmmx; // The full range is -maxr to +maxr, i.e. spanning 2*maxr
            stat.chainr[i].mx = cos( chainrmx * DEGTORAD * 0.5 );
            stat.chainr[i].angle = chainrmx * DEGTORAD * 0.5;
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

    void checkValidity() {
        if (rotmx < 0.0 || rotmx > 180) {
            fprintf (stderr, "ERROR: Maximum orientation change must be in range 0 to 180.\n\n");
            exit (1);
        }
        if (chainrmx < 0.0 || chainrmx > 180) {
            fprintf (stderr, "ERROR: Maximum orientation change for chains must be in range 0 to 180.\n\n");
            exit (1);
        }
        if ( (ptype <0) || (ptype>3) ) {
            fprintf (stderr, "ERROR: Unknown pressure coupling %d. Program only knows: 0 - anisotropic coupling, \
                    1 - isotropic coupling, 2 - isotropic in xy z=const, 3 - isotropic xy V=const.\n\n",ptype);
            exit (1);
        }
        if ( (wlm[0] <0) || (wlm[0] > 7) || (wlm[1] <0) || (wlm[1] > 7)  ) {
            fprintf (stderr, "ERROR: Unknown Wang-Landau method %d %d. Program only knows: 0 - none, \
                    1 - z-direction od 1st particle, 2 - pore in membrane, 3 - zorientation of 0th particle,\
                    4 - distance of fist two particles, 5 - pore around z-axis above CM,\
                    6 - pore around z-axis above 0th particle, 7 - number of particles in contact \n\n",wlm[0],wlm[1]);
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
