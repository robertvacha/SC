/** @file mcsimsystem.h*/

#ifndef MCSIMSYSTEM_H
#define MCSIMSYSTEM_H

#include <iostream>
#include <iomanip>
#include "inicializer.h"
#include "updater.h"


using namespace std;

class MCSimSystem
{
private:
    FILE *outfile,*mov;       // Handle for writing configuration

    Topo topo;                // will maybe contain all the topo stuff in future
    Sim sim;                  // Should contain the simulation options.
    Conf conf;                // Should contain fast changing particle and box(?) information
    FileNames files;

    Updater* updater;

public:
    MCSimSystem() {}
    ~MCSimSystem() {delete updater;}

    /**
     * @brief init
     */
    void init(int argc, char **argv);

    /**
     * @brief equilibrate
     */
    void equilibrate();

    /**
     * @brief productionRun
     */
    void productionRun();

    /**
     * @brief dealloc
     */
    void dealloc();

private:

    /**
     * @brief clearOutfiles Clears outfiles for movie and Wang-Landau method
     */
    void clearOutfiles();

    /**
     * @brief memorydealloc
     */
    int memoryDealloc();

    /**
     * @brief dealloc_pairlist Cleans up: deallocates the memory for the pairlist
     * @param topo
     * @param sim
     * @return
     */
    int deallocPairlist(); // deprecated, done in memoryDealloc()

};



#endif // MCSIMSYSTEM_H
