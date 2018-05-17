#ifndef TOPO_H
#define TOPO_H

#include "macros.h"
#include "structures.h"
#include "../mc/simlib.h"
#include "../mc/mygetline.h"
#include "mpicout.h"

#include <array>

extern MpiCout mcout;

/**
 * @brief The System class
 * [SYSTEM]
 * moleculeName[i] moleculeCount[i]
 */
class System {
public:
    System(string name) : name(name) {}

    string name;
    vector<string> moleculeName;
    vector<int> moleculeCount;

    /**
     * @brief filling the system parameters
     * @param pline
     * @param sysnames
     * @param sysmoln
     * @return
     */
    int fillSystem(char *pline, char* namee) {

        long number;
        int fields;
        char zz[STRLEN];

        trim(pline);
        if (!pline) {
            fprintf (stderr, "TOPOLOGY ERROR: obtained empty line in fil system.\n\n");
            return 0;
        }

        fields = sscanf(pline, "%s %ld", zz, &number);

        if (fields != 2) {
            cerr << "TOPOLOGY ERROR: failed reading system from " << pline << ".\n" << endl;
            return 0;
        }

        cout << namee << " " << zz << " " << number << endl;

        moleculeName.push_back(zz);
        moleculeCount.push_back(number);

        return 1;
    }

    bool operator==(const System& o) const {
        return (this->moleculeCount == o.moleculeCount) && (this->moleculeName == o.moleculeName);
    }

    string toString() {
        stringstream ss;
        ss << name << endl;
        for(unsigned int i; i<moleculeName.size(); ++i) {
            ss << moleculeName[i] << " " << moleculeCount[i] << endl;
        }
        return ss.str();
    }

    void print() {
        cout << toString();
    }
};

/* It would be nice, if this struct would contain all the topo stuff in the end*/
class Topo
{
public:
    System system = System("[System]");
    System pool = System("[Pool]");

    MolIO* molecules;    ///< @brief List of AtomType parameters read from init.top

    bool poolConfig = false;
    bool switchSpecies  = false;
    bool gcSpecies      = false;
    bool existExclusion = false;
    bool exclusions[MAXT][MAXT] = {{false}};

    double sqmaxcut     = 0;    ///< \brief square of distance over which even spherocylinders cannot interact (distance between CM)
    double maxcut       = 0;      ///< \brief distance over which even spherocylinders cannot interact (distance between CM)
    int gcSpeciesCount  = 0;

    std::array<MoleculeParams, MAXMT>  moleculeParam;   ///< \brief parameters for Molecules

    // Multidimensional array not supported yet
    std::array<Ia_param, MAXT> ia_params[MAXT];     ///< \brief parametrization of particles for all interations
    Exters exter;                       ///< \brief external potential - wall

    Topo() { molecules = new MolIO[MAXMT]; }

    Topo(bool switchSpecies, bool gcSpecies, FileNames* files) : switchSpecies(switchSpecies), gcSpecies(gcSpecies) {

        molecules = new MolIO[MAXMT];

        readTopoFile(files);

        mcout.get() << "\nTopology succesfully read. Generating pair interactions..." << endl;

        genParamPairs();
        genTopoParams();

        if(switchSpecies)
            initSwitchList();

        if(gcSpecies)
            initGCList();
    }

    string toString() {
        stringstream ss;

        ss << "[Types]" << endl;
        ss << "#NAME NUM GEOTYP EPS SIGMA ATTR_DIST ATTR_SW PATCH_ANGLE PATCH_SW SC_LEN Para_EPS angle2 Patch_ANGLE2 PASW2 Chiral" << endl;
        for(int i=0; i<MAXT; ++i) {
            if( !ia_params[i][i].name.empty() )
                ss << ia_params[i][i].toString() << endl;
        }

        ss << "[Molecules]" << endl;
        for(int i=0; i<MAXMT; ++i) {
            if( !moleculeParam[i].name.empty() ) {
                ss << moleculeParam[i].toString() << endl;
            }
        }

        ss << system.toString();
        if(poolConfig)
            ss << pool.toString();

        if(existExclusion) {
            ss << "[EXCLUDE]" << endl;
            for(int i=0; i < MAXT; ++i)
                for(int j=i; j<MAXT; ++j)
                    if( exclusions[i][j] )
                        ss << i << " " << j << endl;
        }

        ss << ((exter.exist) ? exter.toString() : string());

        return ss.str();
    }

    bool operator==(Topo& o) const {
        bool same = true;
        for(int i=0; i< MAXT; ++i) {
            if( !(this->ia_params[i] == o.ia_params[i]) ) {
                same = false;
                cerr << "this->ia_params[" << i << "] != o.ia_params[" << i << "]" << endl;
            }
        }

        if( !(this->moleculeParam == o.moleculeParam) )
            cerr << "this->moleculeParam != o.moleculeParam" << endl;

        if( !(this->exter == o.exter) ) {
            cerr << "this->exter != o.exter" << endl;
        }

        if( !( (this->sqmaxcut == o.sqmaxcut) && (this->maxcut == o.maxcut) && (this->gcSpeciesCount == o.gcSpeciesCount) ) )
            cerr << "this->sqmaxcut != o.sqmaxcut || this->maxcut != o.maxcut || this->gcSpecies != o.gcSpecies" << endl;

        if(!(this->system == o.system) )
            cerr << "this->system != o.system" << endl;

        if(!(this->pool == o.pool) )
            cerr << "this->pool != o.pool" << endl;

        return same && (this->sqmaxcut == o.sqmaxcut) && (this->maxcut == o.maxcut) && (this->gcSpeciesCount == o.gcSpeciesCount)
                && (this->moleculeParam == o.moleculeParam) && (this->exter == o.exter) && (this->system == o.system) && (this->pool == o.pool);
    }


    /**
     * @brief generate interations pairs
     * @param (*exlusions)[][]
     */
    void genParamPairs();

    void genTopoParams();

    void info() {
        cout << "Topology:\n";
        int i=0;
        while(!moleculeParam[i].name.empty()) {
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

    void initSwitchList() {
        // count switch types for all molecular types
        int count;
        bool switchPartExist = false;
        for(int i=0; i<MAXMT; i++) {
            if(moleculeParam[i].particleTypes.empty())
                break;
            count =0;
            for(unsigned int j=0; j<moleculeParam[i].switchTypes.size(); j++) {
                if(moleculeParam[i].switchTypes[j] != -1) {
                    count++;
                    switchPartExist = true;
                }
            }
            moleculeParam[i].switchCount = count;
        }

        if (!switchPartExist){
            cerr << "TOPOLOGY ERROR: No switchable particles found, but probability for a switch is not zero!" << endl;
            exit(1);
        }
    }

    void initGCList() {
        bool existGrand = false;
        int i=0;
        while(!moleculeParam[i].name.empty()) {
            if(moleculeParam[i].activity != -1.0 )
                existGrand = true;
            i++;
        }
        if(!existGrand) {
            cout << "In options nGrandCanon != 0, but no activity set for any species in top.init" << endl;
            exit(1);
        }
    }

public:

      void readTopoFile(FileNames* files) {
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

        mcout.get() << "Reading topology...\n" << "Species:" << endl;

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
                        if (!system.fillSystem(pline, name)) {
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
                        if (!pool.fillSystem(pline, name)) {
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
                        existExclusion = true;
                        fflush(stdout);
                        if (!fillExclusions(&pline)) {
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
    int fillExclusions(char **pline);


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
