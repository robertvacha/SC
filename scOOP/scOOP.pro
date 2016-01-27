TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -O3 -msse2 -mfpmath=sse  #-g -fno-inline

#QMAKE_CXXFLAGS+= -fopenmp
#QMAKE_LFLAGS +=  -fopenmp

#QMAKE_CXXFLAGS += -DENABLE_MPI

#QMAKE_CXX = mpicxx
#QMAKE_CXX_RELEASE = $$QMAKE_CXX
#QMAKE_CXX_DEBUG = $$QMAKE_CXX
#QMAKE_LINK = $$QMAKE_CXX
#QMAKE_CC = mpicc

#QMAKE_CFLAGS += $$system(mpicc --showme:compile)
#QMAKE_LFLAGS += $$system(mpicxx --showme:link)
#QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK


SOURCES += main.cpp \
    mc/wanglandau.cpp \
    mc/updater.cpp \
    mc/totalenergycalculator.cpp \
    mc/pairenergycalculator.cpp \
    mc/mygetline2.cpp \
    mc/MoveCreatorMesh.cpp \
    mc/movecreator.cpp \
    mc/inicializer.cpp \
    mc/mesh.cpp \
    mc/simlib.cpp \
    structures/Conf.cpp \
    structures/particle.cpp \
    mc/printStat.cpp \
    structures/topo.cpp \
    mc/randomGenerator.cpp \
    mc/externalenergycalculator.cpp \
    mc/dSFMT-src-2.2.3/dSFMT.c

HEADERS += \
    mc/wanglandau.h \
    structures/Vector.h \
    mc/updater.h \
    mc/totalenergycalculator.h \
    structures/sim.h \
    structures/quaternion.h \
    mc/pairenergycalculator.h \
    mc/movecreator.h \
    mc/inicializer.h \
    structures/macros.h \
    mc/math_calc.h \
    mc/mesh.h \
    mc/mygetline.h \
    mc/simlib.h \
    structures/structures.h \
    structures/particle.h \
    structures/Conf.h \
    mc/printStat.h \
    structures/moleculeparams.h \
    structures/topo.h \
    mc/randomGenerator.h \
    mc/externalenergycalculator.h \
    structures/geometry.h \
    mc/dSFMT-src-2.2.3/dSFMT.h \
    mc/dSFMT-src-2.2.3/dSFMT-params.h \
    mc/dSFMT-src-2.2.3/dSFMT-params521.h \
    mc/dSFMT-src-2.2.3/dSFMT-params1279.h \
    mc/dSFMT-src-2.2.3/dSFMT-params2203.h \
    mc/dSFMT-src-2.2.3/dSFMT-params4253.h \
    mc/dSFMT-src-2.2.3/dSFMT-params11213.h \
    mc/dSFMT-src-2.2.3/dSFMT-params19937.h \
    mc/dSFMT-src-2.2.3/dSFMT-params44497.h \
    mc/dSFMT-src-2.2.3/dSFMT-params86243.h \
    mc/dSFMT-src-2.2.3/dSFMT-params132049.h \
    mc/dSFMT-src-2.2.3/dSFMT-params216091.h \
    mc/dSFMT-src-2.2.3/dSFMT-common.h \
    unitTests/pvectester.h \
    unitTests/vectortester.h

OTHER_FILES += \
    mc/dSFMT-src-2.2.3/FILES.txt \
    mc/dSFMT-src-2.2.3/LICENSE.txt \
    mc/dSFMT-src-2.2.3/README.jp.txt \
    mc/dSFMT-src-2.2.3/README.txt

