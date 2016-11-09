TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -larmadillo -llapack -lblas

QMAKE_CXXFLAGS_RELEASE += -DARMA_NO_DEBUG

SOURCES += main.cpp \
    lib.cpp \
    System.cpp \
    Potential.cpp \
    pBar.cpp \
    LennardJones.cpp

HEADERS += \
    lib.h \
    System.h \
    Potential.h \
    Structs.h \
    pBar.h \
    LennardJones.h

OTHER_FILES += \
    lib.o \
    _lib.so \
    lib.i

QMAKE_CXXFLAGS += -std=c++11
