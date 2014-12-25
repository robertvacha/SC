TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -O3 -msse2 -mfpmath=sse -g

#QMAKE_CXXFLAGS+= -fopenmp
#QMAKE_LFLAGS +=  -fopenmp

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
    mc/externalenergycalculator.cpp

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
    mc/externalenergycalculator.h

