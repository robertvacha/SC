#ifndef TOPO_H
#define TOPO_H

#include "macros.h"
#include "structures.h"
#include "../mc/simlib.h"
#include "../mc/mygetline.h"

/* It would be nice, if this struct would contain all the topo stuff in the end*/
class Topo
{
public:
    double sqmaxcut;    ///< \brief square of distance over which even spherocylinders cannot interact (distance between CM)
    double maxcut;      ///< \brief distance over which even spherocylinders cannot interact (distance between CM)
    int gcSpecies;

    MoleculeParams moleculeParam[MAXMT];   ///< \brief parameters for Molecules

    Ia_param ia_params[MAXT][MAXT];     ///< \brief parametrization of particles for all interations
    Exters exter;                       ///< \brief external potential - wall

    //
    //  METHODS
    //

    Topo():sqmaxcut(0), maxcut(0), gcSpecies(0) {}

    Topo(FileNames* files) : sqmaxcut(0), maxcut(0), gcSpecies(0)  {
        bool exclusions[MAXT][MAXT] = {false};

        readTopoFile(exclusions, files); // EXCLUDE LOADED CORRECTLY 7.8. 2015
    }

    ~Topo() {
        printf ("Deallocating Topo...\n");

        for(int i=0; i<MAXMT; i++) {
            free(moleculeParam[i].name);
        }
    }

    /**
     * @brief generate interations pairs
     * @param (*exlusions)[][]
     */
    void genParamPairs(bool exclusions[MAXT][MAXT]);

    void genTopoParams();

    void info() {
        cout << "Topology:\n";
        int i=0;
        while(moleculeParam[i].name != NULL) {
            cout << moleculeParam[i].name << ", activity:" << moleculeParam[i].activity << endl;
            i++;
        }

        for(int q=0; q<5; q++) {
            for(int w=0; w<5; w++) {
                cout << ia_params[q][w].geotype[0] << " " << ia_params[q][w].geotype[1] <<" "<< ia_params[q][w].epsilon << " | ";
            }
            cout << endl;
        }
    }

private:

    MolIO* molecules;    ///< @brief List of AtomType parameters read from init.top

    char *sysnames[MAXN];   ///< @brief List of MoleculeType names of system
    char *poolNames[MAXN];  ///< @brief List of MoleculeType names of pool

    long  *sysmoln /*[MAXN]*/;
    long  *poolMolNum /*[MAXN]*/;

    bool poolConfig;
    bool nGrandCanon;

    void readTopoFile(bool exclusions[][MAXT], FileNames* files) {
        char *dummy=NULL;
        char line[STRLEN], keystr[STRLEN], molname[STRLEN];
        unsigned size;
        long i=0;
        FILE *infile;
        char *pline=NULL;

        if ((infile = fopen(files->topologyInFile, "r")) == NULL) {
            fprintf (stderr, "\nTOPOLOGY ERROR: Could not open top.init file.\n\n");
            exit (1);
        }

        if(SILENT == 1)
            cout << "Reading topology...\n"
                 << "Species:" << endl;

        molname[0] = ' ';

        pline = (char*) malloc((size_t)STRLEN);
        while (fgets2(line,STRLEN-2,infile) != NULL) {
            strcpy(pline,line);
            if (!pline) fprintf (stderr, "\nTOPOLOGY ERROR: Empty line in topology.\n\n");

            // build one long line from several fragments
            while (continuing(line) && (fgets2(line,STRLEN-1,infile) != NULL)) {
                size=strlen(pline)+strlen(line)+1;
                free(pline);
                pline = (char*) malloc((size_t)size);
                strcat(pline,line);
            }

            strip_comment (pline);
            trim (pline);

            if ((int)strlen(pline) > 0) {
                // get the [COMMAND] key
                if (pline[0] == OPENKEY) {
                    pline[0] = ' ';
                    beforecommand(keystr,pline,CLOSEKEY);
                    upstring (keystr);
                } else {

                    //DEBUG		fprintf (stdout, "Topology read type:%s, %s \n",keystr,pline);
                    if (!strcmp(keystr,"TYPES")) {
                        fflush(stdout);
                        if (!fillTypes(&pline)) {
                            DEBUG_INIT("Something went wrong with filltypes");
                            fprintf (stderr, "\nTOPOLOGY ERROR: in reading types\n\n");
                            topDealoc();
                            free(pline); pline = NULL;
                            exit (1);
                        }
                        DEBUG_INIT("back in init_top");
                        continue;
                    }
                    if (!strcmp(keystr,"MOLECULES")) {
                        DEBUG_INIT("Let's go to the molecules");
                        if (molname[0] == ' ') {
                            beforecommand(molname,pline,SEPARATOR);
                            i=0;
                            while (molecules[i].name != NULL)
                                i++;
                            DEBUG_INIT("in the middle of getting to fillmol");
                            molecules[i].name = (char*) malloc(strlen(molname)+1);
                            strcpy(molecules[i].name, molname);
                            fprintf (stdout, "\nTopology read for molecule: %s \n",molname);
                        }
                        if (!fillMol(molname, pline, molecules)) {
                            fprintf (stderr, "\nTOPOLOGY ERROR: in reading molecules\n\n");
                            topDealoc();
                            free(pline); pline = NULL;
                            exit (1);
                        }
                        if ((dummy = strchr (pline,CLOSEMOL)) != NULL)
                            molname[0] = ' ';
                        continue;
                    }
                    if (!strcmp(keystr,"SYSTEM")) {
                        char name[9] = "system: ";
                        if (!fillSystem(pline,sysnames,&sysmoln, name)) {
                            fprintf (stderr, "\nTOPOLOGY ERROR: in reading system\n\n");
                            topDealoc();
                            free(pline); pline = NULL;
                            exit (1);
                        }
                        continue;
                    }
                    if (!strcmp(keystr, "POOL")) {
                        poolConfig = true;
                        char name[7] = "pool: ";
                        if (!fillSystem(pline,poolNames,&poolMolNum, name)) {
                            fprintf (stderr, "\nTOPOLOGY ERROR: in reading system\n\n");
                            topDealoc();
                            free(pline); pline = NULL;
                            exit (1);
                        }
                        continue;
                    }
                    if (!strcmp(keystr,"EXTER")) {
                        fflush(stdout);
                        if (!fillExter(&pline)) {
                            DEBUG_INIT("Something went wrong with external potential");
                            fprintf (stderr, "\nTOPOLOGY ERROR: in reading external potential\n\n");
                            topDealoc();
                            free(pline); pline = NULL;
                            exit (1);
                        }
                        continue;
                    }
                    if (!strcmp(keystr,"EXCLUDE")) {
                        fflush(stdout);
                        if (!fillExclusions(&pline,exclusions)) {
                            DEBUG_INIT("Something went wrong with exclusions potential");
                            fprintf (stderr, "\nTOPOLOGY ERROR: in reading exclusions\n\n");
                            topDealoc();
                            free(pline); pline = NULL;
                            exit (1);
                        }
                        continue;
                    }

                    fprintf (stderr, "\nTOPOLOGY ERROR: invalid keyword:%s.\n\n", keystr);
                    topDealoc();
                    free(pline); pline = NULL;
                    exit (1);
                }
            }
        }
        //we have sucessfully read topology
        if (pline !=NULL) free(pline);
        pline=NULL;
        fclose (infile);
        fflush (stdout);
    }



    /**
     * @brief Converts the geometrical type string into a number
     * @param geotype
     * @return
     */
    int convertGeotype(char *geotype);



    /**
     * @brief filling the system parameters
     * @param pline
     * @param sysnames
     * @param sysmoln
     * @return
     */
    int fillSystem(char *pline, char *sysnames[], long **sysmoln, char* name);



    /**
     * @brief filing the parameters for types from given strings. Returns 1 on succes.
     * @param pline
     * @return
     */
    int fillTypes(char **pline);


    /**
     * @brief filling pair for which we exlude attraction interaction. Returns 1 on succes.
     * @param pline
     * @param (*exlusions)[][]
     * @return
     */
    int fillExclusions(char **pline, bool exlusions[][MAXT]);


    /**
     * @brief filling the parameters of external potentail - wall. Returns 1 on succes.
     * @param pline
     * @return
     */
    int fillExter(char **pline);


    /**
     * @brief filling the parameters for molecules
     * @param molname
     * @param pline
     * @param molecules
     * @return
     */
    int fillMol(char *molname, char *pline, MolIO *molecules);


    void topDealoc();
};

#endif // TOPO_H
