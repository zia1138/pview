@echo off

rem This batch script is used to create a 64-bit Visual Studio C++ 2008 
rem build environment for PVIEW.  See the README for details.

rem Create a command line shortcut that's like this.
rem Note you might need to update the version number.
rem %comspec /k ""C:\qt-everywhere-opensource-src-4.6.3\bin\qtvars.bat""
rem

echo Setting up a Qt environment...

set QTDIR=C:\qt-everywhere-opensource-src-4.6.3
echo -- QTDIR set to C:\qt-everywhere-opensource-src-4.6.3
set PATH=C:\qt-everywhere-opensource-src-4.6.3\bin;%PATH%
echo -- Added C:\qt-everywhere-opensource-src-4.6.3 to PATH
set QMAKESPEC=win32-msvc2008
echo -- QMAKESPEC set to "win32-msvc2008"

rem Use 64-bit environment variables.
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat" amd64


