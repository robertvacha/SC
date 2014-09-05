/** @file main.cpp*/

#include <time.h>
#include "mc/mcsimsystem.h"

using namespace std;

int main()
{
    MCSimSystem mc;

    mc.init();
    mc.equilibrate();
    mc.productionRun();
    mc.dealloc();

    return 0;
}



