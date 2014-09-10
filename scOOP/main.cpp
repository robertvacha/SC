/** @file main.cpp*/

#include <time.h>
#include "mc/mcsimsystem.h"

using namespace std;

extern bool print = false;

int main()
{
    MCSimSystem mc;

    mc.init();
    mc.equilibrate();
    mc.productionRun();
    mc.dealloc();

    return 0;
}



