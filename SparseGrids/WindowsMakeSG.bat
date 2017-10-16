@if /i "%1" == "clean" GOTO Clean

@ECHO Compiling static library and executable

copy TasmanianSparseGrid.windows.hpp TasmanianSparseGrid.hpp

del *.obj

cl -c *.cpp /DTSG_STATIC /Ox /EHsc /openmp

lib tsg*.obj Tasmanian*.obj /OUT:libtasmaniansparsegrid_static.lib

cl tasgrid*.obj libtasmaniansparsegrid_static.lib /Fe:tasgrid.exe

del *.obj

cl /EHsc /Ox /openmp /DTSG_DLL /LD tsg*.cpp TasmanianSparseGrid.cpp /Fe:libtasmaniansparsegrid.dll

@GOTO End

:Clean

@ECHO Cleaning

del *.obj *.dll *.lib *.exe

@GOTO End

:End
