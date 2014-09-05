SC
==

C++ code for Monte Carlo simulations of patchy spherocylinders

Compilation:

in scOOP/ directory:

$ cmake .
$ make

executable SC will be created

Parallel tempering (enabling MPI) compile with:

$ cmake . -DENABLE_MPI=ON
$ make
$ mpirun -np x SC

, where x is the number of threads
