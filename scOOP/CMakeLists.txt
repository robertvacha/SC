cmake_minimum_required(VERSION 2.8.5 FATAL_ERROR)
project("SC")
SET(sampleName SC)

# -----------------
#   User options
# -----------------

#option(ENABLE_OPENMP   "Try to use OpenMP parallization" off)
option(ENABLE_MPI      "Enable MPI code" off)

# ------------------
#   Compiler Flags
# ------------------

set (CMAKE_CXX_FLAGS " -O3 -msse2 -mfpmath=sse ")

# -----------------------
#   Link with MPI
# -----------------------
if(ENABLE_MPI)
  find_package(MPI)
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
structures/macros.h                               # MACROS
structures/Vector.h                               # VECTOR
structures/quaternion.h                           # QUATERNION
structures/structures.h                           # STRUCTURES
structures/moleculeparams.h                       # MOLECULE PARAMS
structures/sim.h                                  # SIM
structures/topo.h structures/topo.cpp             # TOPO
structures/Conf.h structures/Conf.cpp             # CONF
structures/particle.h structures/particle.cpp     # PARTICLE
#
# MC
#
mc/math_calc.h                                               # MATHS
mc/wanglandau.h mc/wanglandau.cpp                            # WANG-LANDAU
mc/updater.h mc/updater.cpp                                  # UPDATER
mc/totalenergycalculator.h mc/totalenergycalculator.cpp      # TOTAL-ENERGY CALCULATOR
mc/pairenergycalculator.h mc/pairenergycalculator.cpp        # PAIR-ENERGY CALCULATOR
mc/mygetline.h mc/mygetline2.cpp                             # MAC PORT FOR IO FUNCTION 
mc/movecreator.h mc/movecreator.cpp mc/MoveCreatorMesh.cpp   # MOVES - MOVECREATOR
mc/mcsimsystem.h mc/mcsimsystem.cpp                          # BASE - MCSIMSYSTEM
mc/inicializer.h mc/inicializer.cpp                          # INICIALIZER
mc/mesh.h mc/mesh.cpp                                        # MESH
mc/simlib.h mc/simlib.cpp                                    # SIMLIB
mc/ran2.h mc/ran2.cpp                                        # RANDOM GENERATOR
mc/printStat.h mc/printStat.cpp)                             # PRINTSTAT


