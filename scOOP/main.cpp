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

    return 0;
}



