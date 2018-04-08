@if /i "%1" == "clean" GOTO Clean

@ECHO Compiling static library and executable

copy TasmanianSparseGrid.windows.hpp TasmanianSparseGrid.hpp
copy TasmanianSparseGrid.windows.h TasmanianSparseGrid.h

del *.obj

cl -c *.cpp /DTSG_STATIC /D_USE_MATH_DEFINES /D_TASMANIAN_WINDOWS_ /Ox /EHsc /openmp

lib tsg*.obj Tasmanian*.obj /OUT:libtasmaniansparsegrid_static.lib

cl tasgrid*.obj libtasmaniansparsegrid_static.lib /Fe:tasgrid.exe

del *.obj

@ECHO Compiling shared library

cl /EHsc /Ox /openmp /DTSG_DLL /D_USE_MATH_DEFINES /D_TASMANIAN_WINDOWS_ /LD tsg*.cpp TasmanianSparseGrid.cpp /Fe:libtasmaniansparsegrid.dll

@GOTO End

:Clean

@ECHO Cleaning

del *.obj *.dll *.lib *.exe *.exp

del tasmanianConfig.hpp

@GOTO End

:End
