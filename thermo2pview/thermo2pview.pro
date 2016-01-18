TEMPLATE = app
CONFIG += qt core gui widgets warn_on release embed_manifest_exe
QT += core gui widgets
FORMS += thermo2pview.ui
SOURCES += thermo2pview.cpp b64.c
HEADERS += thermo2pview.hpp centroider.hpp b64.h
QMAKE_CFLAGS += "/D_HAS_ITERATOR_DEBUGGING=0 /D_SECURE_SCL=0"
QMAKE_CXXFLAGS += "/D_HAS_ITERATOR_DEBUGGING=0 /D_SECURE_SCL=0"
