@if /i "%1" == "clean" GOTO Clean

@ECHO Compiling static library and executable

del *.obj

cl /I ../SparseGrids -c *.cpp /DTSG_STATIC /D_USE_MATH_DEFINES /D_TASMANIAN_WINDOWS_ /Ox /EHsc /openmp

lib tdr*.obj Tasmanian*.obj /OUT:libtasmaniandream_static.lib

cl tasdream*.obj libtasmaniandream_static.lib ..\SparseGrids\libtasmaniansparsegrid_static.lib /Fe:tasdream.exe

@ECHO Skipping shared library

@GOTO End

:Clean

@ECHO Cleaning

del *.obj *.dll *.lib *.exe

@GOTO End

:End
