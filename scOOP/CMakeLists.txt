cmake_minimum_required(VERSION 2.8.5 FATAL_ERROR)
project("SC")
SET(sampleName SC)

# -----------------
#   User options
# -----------------

option(ENABLE_OPENMP   "Try to use OpenMP parallization" off)
option(ENABLE_MPI      "Enable MPI code" off)

# ------------------
#   Compiler Flags
# ------------------

set (CMAKE_CXX_FLAGS " -Ofast -march=native -std=c++11")

# -------------------------------------
#   Mersenne Twister requires
# ----=--------------------------------
if(ENABLE_TWISTER)
    add_definitions(-DMERSENNETWISTER)
endif()

# -------------------------------------
#   OPEN_MP PARALELLIZATION
# ----=--------------------------------
if(ENABLE_OPENMP)
    add_definitions(-DOMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# -------------------------------------
#   For Testing
# ----=--------------------------------
if(TESTING)
    add_definitions(-DTESTING)
endif()

# -----------------------
#   Link with MPI
# -----------------------
if(ENABLE_MPI)
  find_package(MPI)
  include(CMakeForceCompiler)
  CMAKE_FORCE_C_COMPILER(mpicc ${CMAKE_C_COMPILER_ID})
  CMAKE_FORCE_CXX_COMPILER(mpicxx ${CMAKE_CXX_COMPILER_ID})
  if(MPI_CXX_FOUND)
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
    if (APPLE)
      set(CMAKE_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
    else()
      set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
    endif()
    set(LINKLIBS ${LINKLIBS} ${MPI_CXX_LIBRARIES})
    include_directories(${MPI_CXX_INCLUDE_PATH})
    add_definitions(-DENABLE_MPI)
  endif()
endif()

add_executable (SC main.cpp
#
# STRUCTURES
#
structures/mpicout.h     			  # mpicout
structures/macros.h                               # MACROS
structures/Vector.h                               # VECTOR
structures/quaternion.h                           # QUATERNION
structures/structures.h                           # STRUCTURES
structures/moleculeparams.h                       # MOLECULE PARAMS
structures/sim.h                                  # SIM
structures/topo.h structures/topo.cpp             # TOPO
structures/Conf.h structures/Conf.cpp             # CONF
structures/particle.h structures/particle.cpp     # PARTICLE
structures/geometry.h     			  # GEOMETRY
structures/statistics.h     			  # GEOMETRY
#
# MC
#
mc/math_calc.h                                                # MATHS
mc/wanglandau.h mc/wanglandau.cpp                             # WANG-LANDAU
mc/updater.h mc/updater.cpp                                   # UPDATER
mc/totalenergycalculator.h 			       	      # TOTAL-ENERGY CALCULATOR
mc/externalenergycalculator.h mc/externalenergycalculator.cpp # EXTERNAL-ENERGY CALCULATOR
mc/pairenergycalculator.h mc/pairenergycalculator.cpp         # PAIR-ENERGY CALCULATOR
mc/paire.h mc/paire.cpp         # PAIR-ENERGY CALCULATOR
mc/mygetline.h mc/mygetline2.cpp                              # MAC PORT FOR IO FUNCTION 
mc/movecreator.h mc/movecreator.cpp mc/MoveCreatorMesh.cpp    # MOVES - MOVECREATOR
mc/inicializer.h mc/inicializer.cpp                           # INICIALIZER
mc/mesh.h mc/mesh.cpp                                         # MESH
mc/simlib.h mc/simlib.cpp                                     # SIMLIB
mc/randomGenerator.h mc/randomGenerator.cpp                   # RANDOM GENERATOR
mc/clust.h mc/clust.cpp
mc/dSFMT-src-2.2.3/dSFMT.c
mc/dSFMT-src-2.2.3/dSFMT.h
mc/dSFMT-src-2.2.3/dSFMT-params.h
)                              # PRINTSTAT

#target_link_libraries(SC ${MPI_C_LIBRARIES}) 


