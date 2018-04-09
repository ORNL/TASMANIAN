In order to compile TASMANIAN using VC++

1. Install Microsoft Visual C++

2. Extract the .zip files, for example into 
   C:\Users\Bill\Downloads\Tasmanian

3. Run "Development Command Prompt VS..." from the windows menu

4. Navigate to the folder with the source code using the cd command, 
   in the example above: cd \Users\Bill\Downloads\Tasmanian

5. Run the WindowsMake.bat script, i.e., type: 
   WindowsMake.bat

The WindowsMake scrip will compile the libraries and the 
tasgrid.exe and tasdream.exe executables.
The files will be copied in \Users\Bill\Downloads\Tasmanian

6. Run basic test:
   WindowsMake.bat test

(optional) Enabling MATLAB interface:

7. Right-click and edit
   C:\Users\Bill\Downloads\Tasmanian\InterfaceMATLAB\tsgGetPaths.m
   
   Choose a work folder, for example:
   C:\Users\Bill\Downloads\Tasmanian\InterfaceMATLAB\WorkFiles\
   (make sure to create the folder)
   
   Set the path to the executable:
   C:\Users\Bill\Downloads\Tasmanian\tasgrid.exe
   
8. Add the Tasmanian interface files to the MATLAB path
   see the final command of the WindowsMake.bat command
   
   addpath('C:\Users\Bill\Downloads\Tasmanian\InterfaceMATLAB\')


Note: to delete all files created by the compile use: 
      WindowsMake.bat clean

Note: visual studio seems to default to 32-bit even on 64-bit systems
      64-bit code will perform better, and Python scripts need matching
      dynamic libraries (e.g., 64-bit Anaconda will only work with
      64-bit Tasmanian)
      In order to enable 64-bit Tasmanian, after step 3
      cd VC
      vcvarsall.bat amd64
      (amd64 is the name of the 64-bit standard and is used by both
      Intel and AMD cpus)

Note: if you experience problems with OpenMP or other compile options
      you can manually edit the commands in:
      
      C:\Users\Bill\Downloads\Tasmanian\SparseGrids\WindowsMakeSG.bat
      C:\Users\Bill\Downloads\Tasmanian\DREAM\WindowsMakeDREAM.bat

      for example, removing /openmp will remove OpenMP

Note: cmake can be used under Windows and the projects compile from
      the command prompt using the msbuild.exe <project> command.
      The INSTALL.vcxproj works as intended, i.e., it install the
      source files as appropriate.

