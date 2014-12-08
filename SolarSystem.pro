TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    celestialbody.cpp \
    solarsystem.cpp \
    vec3.cpp \
    rk4.cpp \
    verlet.cpp \
    gaussiandeviate.cpp

HEADERS += \
    celestialbody.h \
    solarsystem.h \
    vec3.h \
    rk4.h \
    verlet.h \
    gaussiandeviate.h

unix: LIBS += -L/usr/local/lib -larmadillo

INCLUDEPATH += /usr/include
DEPENDPATH += /usr/include



