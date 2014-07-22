/** @file main.cpp*/

#include <time.h>
#include "mc/mcsimsystem.h"

using namespace std;



int main()
{
    clock_t time = clock();

    MCSimSystem mc;

    mc.init();
    mc.equilibrate();
    mc.productionRun();
    mc.dealloc();

    time = clock() - time;

    cout << time << endl;

    return 0;
}




