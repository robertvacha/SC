SC
==

C++ code for Monte Carlo simulations of patchy spherocylinders

Compiling
=========

in scOOP/ directory:

    $ cmake .
    $ make

executable SC will be created

Parallel tempering (enabling MPI) compile with:

    $ cmake . -DENABLE_MPI=ON
    $ make
    $ mpirun -np x SC

, where x is the number of threads

Currently Grand-Canonical works only for single atomic type simulation with random init

NOTE: dont use GrandCanonical with:
    Replica Exchange!!!
    
      
