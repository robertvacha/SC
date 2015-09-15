SC
==

C++ code for Monte Carlo simulations of patchy spherocylinders

Compiling
=========

in scOOP/ directory:

    $ module add cmake
    $ cmake .       
    $ make

executable SC will be created.

Parallel tempering (enabling MPI) compile with:

    $ module add mpich
    $ cmake . -DENABLE_MPI=ON -DENABLE_OPENMP=OFF
    $ make
    $ mpirun -np x SC

, where x is the number of threads

Parallel acceleration (enabling OPENMP) compile with:

    $ cmake . -DENABLE_OPENMP=ON -DENABLE_MPI=OFF
    $ make          ($ module add openmpi)
    $ export OMP_NUM_THREADS=x
    $ ./SC

, where x is the number of threads

Generating documentation:

    $ doxygen

Testing:

    $ cd Tests/
    $ ./test

GrandCanonical finished.

Use scripts/movie-scripts/python 2.x/sc35-movie.py -g 1 for grandCanonical movie

Version tested against old *.c program. All test succesful.

NOTE: to submit a job with grandCanonically active species use runGrand

NOTE: dont use GrandCanonical with:
    Replica Exchange!!!
    
      
