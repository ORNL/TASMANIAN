@if /i "%1" == "clean" GOTO Clean
@if /i "%1" == "test" GOTO Test

@ECHO Compiling SparseGrids

copy Config\AltBuildSystems\TasmanianConfig.hpp SparseGrids\TasmanianConfig.hpp

cd SparseGrids
call WindowsMakeSG.bat
cd..

cd DREAM
call WindowsMakeDREAM.bat
cd ..

copy SparseGrids\libtasmaniansparsegrid_static.lib libtasmaniansparsegrid_static.lib
copy SparseGrids\libtasmaniansparsegrid.lib libtasmaniansparsegrid.lib
copy SparseGrids\libtasmaniansparsegrid.dll libtasmaniansparsegrid.dll
copy SparseGrids\tasgrid.exe tasgrid.exe
copy SparseGrids\gridtest.exe gridtest.exe
copy Config\AltBuildSystems\example_sparse_grids.py example_sparse_grids.py

copy SparseGrids\GaussPattersonRule.table GaussPattersonRule.table

copy DREAM\libtasmaniandream_static.lib libtasmaniandream_static.lib
copy DREAM\tasdream.exe tasdream.exe

copy Config\AltBuildSystems\TasmanianSG.windows.py TasmanianSG.py
copy Examples\example_sparse_grids.py example_sparse_grids.py
copy Config\AltBuildSystems\testTSG.py testTSG.py
copy InterfacePython\sandbox.py sandbox.py

copy Config\AltBuildSystems\tsgGetPaths.m InterfaceMATLAB\tsgGetPaths.m

@ECHO To enable MATLAB interface
@ECHO edit InterfaceMATLAB\tsgGetPaths.m
@ECHO add InterfaceMATLAB to the MATLAB path
@ECHO.
@ECHO addpath('%cd%\InterfaceMATLAB')
@ECHO.

@ECHO WARNING: the WindowsMake.bat script is deprecated and will be soon removed
@ECHO          CMake is the recommended way to use Tasmanian under MS Windows

@GOTO End

:Clean

@ECHO Cleaning SparseGrids

cd SparseGrids
call WindowsMakeSG.bat clean
cd..

cd DREAM
call WindowsMakeDREAM.bat clean
cd ..

del *.obj *.dll *.lib *.exe *.py *.pyc *.table
del testSave

@GOTO End

:Test

gridtest.exe
tasdream.exe -test

@GOTO End

:End
