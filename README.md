SC 4.0
==

C++ code for Monte Carlo simulations of patchy spherocylinders

Compiling
=========

in scOOP/ directory:

    $ module add cmake
    $ cmake .
    $ make -j number_of_cores

executable SC will be created.

Parallel tempering or Multiple Walker Wang-Landau:
Enable MPI, compile with:

    $ module add mpich
    $ cmake . -DENABLE_MPI=ON -DENABLE_OPENMP=OFF
    $ make -j number_of_cores
    $ mpirun -np number_of_threads SC

Generating documentation:

    $ doxygen

Testing:

    $ cd Tests/
    $ ./test
    $ cd Interactions_tests/
    $ ./test.sh

Comprehensive comparisons to past versions (doi 10.1063/1.4933229)

folders : SC_PSC_MEMBRANE_WANG and System_averages_tests

Use scripts/movie-scripts/python 2.x/sc35-movie.py for vmd openable movie
Use scripts/movie-scripts/python 2.x/sc35-movie.py -g 1 for grandCanonical

NOTE: to submit a job with grandCanonically active species use runGrand
