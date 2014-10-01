/** @file main.cpp*/

#include <time.h>
#include "mc/randomGenerator.h"
#include "mc/mcsimsystem.h"

using namespace std;

int main(int argc, char** argv)
{
    MCSimSystem mc;

    mc.init(argc,argv);
    mc.equilibrate();
    mc.productionRun();
    mc.dealloc();

    /*
    TODO LIST:
    GrandCanonical:
     - movie -> each snapshot requires count of types
     - conf::sysvolume -> changing after deletes/inserts
     - rewrite initTopo -> usage of IO class Molecules
     */
    return 0;
}



