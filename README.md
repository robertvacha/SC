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

Currently Grand-Canonical works for single atomic type simulation
You can mix several grandcanonical species with non-grandcanonical species

NOTE: to submit a job with grandCanonically active species use runGrand

NOTE: dont use GrandCanonical with:
    Replica Exchange!!!
    
      