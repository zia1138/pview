# Default pview configuration.
TEMPLATE = app
QT += core gui widgets

CONFIG += qt core gui widgets warn_on release embed_manifest_exe
CONFIG -= app_bundle # no MacOS app bundle for now
#CONFIG += static staticlib
FORMS += MSWin.ui MSMSWin.ui XICWin.ui MSConfigWin.ui 
FORMS += MassWin.ui AAWin.ui RecalWin.ui
SOURCES += MSWin.cpp MSMSWin.cpp dialogs.cpp xml.cpp lcms.cpp 
SOURCES += analysis.cpp ms2.cpp b64.c fdr.cpp 
HEADERS += MSWin.hpp MSMSWin.hpp dialogs.hpp xml.hpp msconfig.hpp 2dtree.hpp 
HEADERS += lcms.hpp analysis.hpp ms2.hpp b64.h util.hpp fft.hpp
HEADERS += fdr.hpp centroider.hpp
## Might fix issue with MacOS mountain lion compile
##QMAKE_CXXFLAGS += -fpermissive
# Fast XML parser.
unix {
   LIBS += -lexpat
}
macx {
   # Build a universal binary for each architecture. 
   CONFIG += x86_64
   LIBS += -lexpat
}
win32 {
   QMAKE_CFLAGS += "/D_HAS_ITERATOR_DEBUGGING=0 /D_SECURE_SCL=0"
   QMAKE_CXXFLAGS += "/D_HAS_ITERATOR_DEBUGGING=0 /D_SECURE_SCL=0"
# for the 32-bit version
#   INCLUDEPATH += "C:/Program Files (x86)/Expat 2.0.1/Source/lib"
#   LIBS += "C:/Program Files (x86)/Expat 2.0.1/bin/libexpat.lib"
# for the 64-bit version
   INCLUDEPATH += "C:/Program Files/Expat 2.1.0/Source/lib"
   LIBS += "C:/Program Files/Expat 2.1.0/bin/libexpat.lib"
}
