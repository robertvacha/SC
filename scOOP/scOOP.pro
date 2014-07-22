TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

QMAKE_CXXFLAGS_RELEASE += -O3 -march=native -msse2 -mfpmath=sse

SOURCES += main.cpp \
    mc/wanglandau.cpp \
    mc/updater.cpp \
    mc/totalenergycalculator.cpp \
    mc/pairenergycalculator.cpp \
    mc/mygetline2.cpp \
    mc/MoveCreatorMesh.cpp \
    mc/movecreator.cpp \
    mc/mcsimsystem.cpp \
    mc/inicializer.cpp \
    mc/mesh.cpp \
    mc/simlib.cpp \
    mc/ran2.cpp \
    structures/Conf.cpp \
    structures/particle.cpp \
    mc/printStat.cpp \
    structures/topo.cpp

HEADERS += \
    mc/wanglandau.h \
    structures/Vector.h \
    mc/updater.h \
    mc/totalenergycalculator.h \
    structures/sim.h \
    structures/quaternion.h \
    mc/pairenergycalculator.h \
    mc/movecreator.h \
    mc/mcsimsystem.h \
    mc/inicializer.h \
    structures/macros.h \
    mc/math_calc.h \
    mc/mesh.h \
    mc/mygetline.h \
    mc/simlib.h \
    structures/structures.h \
    mc/ran2.h \
    structures/particle.h \
    structures/Conf.h \
    mc/printStat.h \
    structures/moleculeparams.h \
    structures/topo.h

