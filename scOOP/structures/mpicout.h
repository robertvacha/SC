#ifndef MPICOUT_H
#define MPICOUT_H

#include <fstream>

class MpiCout {
public:
    MpiCout(int rank) : rank(rank) {}
    int rank;
    ofstream nullstream;

    std::ostream& get() {
        if(rank == 0)
            return std::cout;
        else
            return nullstream;
    }
};

#endif // MPICOUT_H
