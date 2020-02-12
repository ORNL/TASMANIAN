# Installation

### Document Sections
* Requirements
* Install using CMake: the preferred way
* Install with the `install` script-wrapper around CMake
* Install with (basic) GNU Make
* Install with Python Pip
* Install with Spack
* Install on MS Windows platform
* Install folder structure
* Linking to Tasmanian: CMake Package Config
* Known issues

### Requirements

Minimum requirements to use Tasmanian:
* a C/C++ compiler and either [CMake](https://cmake.org/) or [GNU Make](https://www.gnu.org/software/make/) build engines.

Recommended additional features:
* [Python](https://www.python.org/) with with [NumPy](http://www.numpy.org/) 1.10 (or newer) and [CTypes](https://docs.python.org/3/library/ctypes.html) packages
* [Basic Linear Algebra Subroutine (BLAS)](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms) implementation
* [OpenMP](https://en.wikipedia.org/wiki/OpenMP) implementation (usually included with the compiler)

Optional features:
* Acceleration using Nvidia [linear algebra libraries](https://developer.nvidia.com/cublas) and custom [CUDA kernels](https://developer.nvidia.com/cuda-zone)
* GPU accelerated linear algebra using [UTK MAGMA library](http://icl.cs.utk.edu/magma/)
* Basic [Python matplotlib](https://matplotlib.org/) support
* Fully featured [MATLAB/Octave](https://www.gnu.org/software/octave/) interface via wrappers around the command-line tool
* Fortran 90/95 interface using [gfortran](https://gcc.gnu.org/wiki/GFortran) or [ifort](https://software.intel.com/en-us/intel-compilers) or [pgf90](https://www.pgroup.com/index.htm)
* Addon templates for the [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface)
* [Doxygen](http://www.doxygen.org/) documentation

**Note:** with the exception of the Intel and PGI compilers and the MAGMA library, the rest of the software is included in the repositories of most Linux distributions, e.g., [Ubuntu](https://www.ubuntu.com/) or [Fedora](https://getfedora.org/), as well as [Mac OSX Homebrew](https://brew.sh/).

| Feature | Tested versions     | Recommended      |
|----|----|----|
| gcc     | 5 - 8               | any              |
| clang   | 4 - 8               | 4 or 5           |
| icc     | 18.0                | 18.0             |
| xl      | 16.0                | 16.0             |
| pgi     | 19.4                | 19.4             |
| cmake   | 3.10 - 3.15         | 3.10             |
| python  | 2.7, 3.5, 3.6       | 3.5 or 3.6       |
| anaconda| 5.3                 | 5.3              |
| OpenBlas| 0.2.18, 0.2.20      | 0.2.18 or 0.2.20 |
| ATLAS   | 3.10                | 3.10             |
| ESSL    | 6.2                 | 6.2              |
| CUDA    | 8.0 - 10.2          | 10.2             |
| libiomp | 5.0                 | 5.0              |
| MAGMA   | 2.5.1               | 2.5.1            |
| Doxygen | 1.8.13              | 1.8.13           |
| MPI     | 2.1, 3.1            | 3.1              |

### Install using CMake: the preferred way

The preferred way to install Tasmanian is to use the included CMake build script, which requires CMake version 3.10 or newer.

* The commands for an out-of-source CMake build:
```
  mkdir Build
  cd Build
  cmake <options> <path-to-Tasmanian-source>
  make
  make test
  make install
  make test_install
```

* Standard CMake options are accepted:
```
  -D CMAKE_INSTALL_PREFIX:PATH=<install-prefix> (install folder for the make install command)
  -D CMAKE_BUILD_TYPE:STRING=<Debug/Release>    (set debug flags or default optimization flags)
  -D BUILD_SHARED_LIBS:BOOL=<ON/OFF>            (pick shared/static libs, undefined builds both)
  -D CMAKE_CXX_COMPILER:PATH=<path>             (specify the C++ compiler)
  -D CMAKE_CUDA_COMPILER:PATH=<path>            (specify the CUDA nvcc compiler)
  -D CMAKE_CXX_FLAGS:STING=<flags>              (set additional flags)
```

* List of Tasmanian specific CMake options (all default to **OFF**):
```
  -D Tasmanian_ENABLE_OPENMP:BOOL=<ON/OFF>      (recommended)
  -D Tasmanian_ENABLE_BLAS:BOOL=<ON/OFF>        (recommended)
  -D Tasmanian_ENABLE_PYTHON:BOOL=<ON/OFF>      (recommended)
  -D Tasmanian_ENABLE_RECOMMENDED:BOOL=<ON/OFF> (enable the above and the -O3 flag)
  -D Tasmanian_ENABLE_CUDA:BOOL=<ON/OFF>        (stable)
  -D Tasmanian_ENABLE_MAGMA:BOOL=<ON/OFF>       (stable)
  -D Tasmanian_MATLAB_WORK_FOLDER:PATH=""       (stable)
  -D Tasmanian_ENABLE_DOXYGEN:BOOL=<ON/OFF>     (stable)
  -D Tasmanian_ENABLE_FORTRAN:BOOL=<ON/OFF>     (mostly stable)
  -D Tasmanian_ENABLE_MPI:BOOL=<ON/OFF>         (mostly stable)
```

* Acceleration options:
    * OpenMP allows Tasmanian to use more than one CPU core, which greatly increases the performance
    * Basic Linear Algebra Subroutines (BLAS) is a standard with many implementations,
      e.g., [https://www.openblas.net/](https://www.openblas.net/); optimized BLAS improves the
      performance when using evaluate commands on grids with many points or working with models with many outputs
    * CUDA is a C++ language extension that allows Tasmanian to leverage the computing power of Nvidia GPU devices,
      which greatly enhances the performance of `evaluateFast()` and `evaluateBatch()` calls
    * Matrix Algebra on GPU and Multicore Architectures (MAGMA) is a library for CUDA accelerated linear
      algebra developed at the University of Tennessee at Knoxville
    * MPI allows the use of distributed memory in Bayesian inference, parallel model construction, and send/receive grid through an MPI comm
* The **Tasmanian_ENABLE_RECOMMENDED** option searches for OpenMP, BLAS, and Python, enables the options (if possible) and also sets the `-O3` flag
* Additional interfaces are available, beyond the default C/C++ library and the command line tools:
    * Python and Fortran require appropriate interpreter and compiler
    * The MATLAB/Octave interface requires a work-folder with read/write permission for temporary files.
      The interface is enabled by setting **Tasmanian_MATLAB_WORK_FOLDER** to a valid read/write location.
* The Doxygen option will build the HTML documentation

* Options to adjust the testing environment: by default Tasmanian testing will use the system provided OpenMP parameters and run tests on all visible Nvidia GPU devices; specific number of threads and device can be selected (note that only the testing environment is affected here):
```
 -D Tasmanian_TESTS_OMP_NUM_THREADS=<number-of-threads-for-testing> (only used with OpenMP)
 -D Tasmanian_TESTS_GPU_ID=<cuda-device-id-for-testing>             (only used with CUDA)
```

* Additional commands to guide the CMake `find_package()` modules:
```
  -D PYTHON_EXECUTABLE:PATH         (specify the Python interpreter)
  -D CMAKE_CUDA_COMPILER:PATH       (specify the CUDA nvcc compiler)
  -D MPI_CXX_COMPILER:PATH=<path>   (specify the MPI compiler wrapper)
  -D CMAKE_Fortran_COMPILER:PATH    (specify the Fortran compiler)
  -D Tasmanian_MAGMA_ROOT_DIR:PATH  (specify the path to the MAGMA installation)
```

* Alternatives allowing to directly specify libraries and bypass `find_package()` altogether:
```
  -D BLAS_LIBRARIES
  -D Tasmanian_MAGMA_LIBRARIES
  -D Tasmanian_MAGMA_INCLUDE_DIRS
```

* Extra options are available in case CMake fails to find a required dependency, e.g., `find_package()` sometimes fails to acknowledge that the ACML implementation of BLAS depends on both `libgfortran` and `libgomp`; the manual options below should not be necessary in the majority of cases:
```
  -D Tasmanian_EXTRA_LIBRARIES:STRING   (add more libraries as dependencies)
  -D Tasmanian_EXTRA_INCLUDE_DIRS:PATH  (add more include paths to search for headers)
  -D Tasmanian_EXTRA_LINK_DIRS:PATH     (appends more link paths to search for libraries)
```

* Options helpful to Tasmanian developers:
```
  -D DOXYGEN_INTERNAL_DOCS=YES  (include the documentation of the Tasmanian internals)
```

### Install with the `install` script-wrapper around CMake

Tasmanian also includes an `install` script that wraps around CMake and automatically calls the build commands.
The script uses `Tasmanian_ENABLE_RECOMMENDED` option and allows for other options to be enabled/disabled
with command line switches.

* Basic usage of the `install` script
```
  ./install <install-path> <optional: matlab-work-folder> <extra switches>
  ./install --help  (list all switches)
```
* Example call that enables MATLAB/Octave and CUDA
```
  ./install /home/me/Tasmanian /home/me/Tasmanian/WorkFolder -cuda=/usr/local/cuda-9.2/bin/nvcc
```
* Additional notes:
    * the script must be called from the main source code folder
    * using absolute paths is strongly recommended
    * if the MATLAB work-folder is omitted, the MATLAB interface will be disabled

### Install with (basic) GNU Make

The core capabilities of Tasmanian can be build with a few simple GNU Make commands.
The basic build engine is useful for quick testing and exploring Tasmanian, or
if CMake is unavailable or unwanted.
Acceleration options other than OpenMP are not supported in the basic mode, but could be enabled
by manually editing `Config/AltBuildSystems/Makefile.in` configuration file.

* Using GNU Make with `g++` and optionally `gfortran` and `/usr/bin/env python`
```
  make
  make test     (will fail if /usr/bin/env python is missing the numpy or ctypes modules)
  make matlab   (optional: sets matlab work folder to ./tsgMatlabWorkFolder/)
  make python3  (optional: sets #!/usr/bin/env python3 in place of /usr/bin/env python)
  make fortran  (optional: compile Fortran libraries)
  make examples
  make clean
```
In the basic mode, the source folder will become the installation folder, i.e.,
the libraries, executables and Python modules will be build in the source folder
and the headers and the Fortran module will be copied to the `include` folder.

### Install with Python Pip

Tasmanian is included in the Python Pip index: [https://pypi.org/project/Tasmanian/](https://pypi.org/project/Tasmanian/)
```
  python3 -m pip install scikit-build packaging numpy --user (required dependencies)
  python3 -m pip install Tasmanian --user                    (user installation)
  python3 -m pip install Tasmanian                           (virtual env installation)
```
The Tasmanian module is not a regular Python-only project but a wrapper around C++ libraries, hence some limitations apply:
* The `scikit-build`, `packaging` and `numpy` dependencies have to be manually installed first.
* Only user installations are supported, installation for all users is possible with CMake but not Pip.
* Python virtual environments are supported, as well as Linux, Mac and Windows operating systems.

The pip installer will enable only the recommended options, if the required libraries are found automatically by CMake.
CUDA acceleration is not available through Pip and there is currently no way to manually specify the BLAS libraries.
Only the C++ and Python interfaces can be installed through Pip.

### Install with Spack

Tasmanian is also included in Spack: [https://spack.io/](https://spack.io/)
```
 spack install tasmanian@7.0+openmp+blas+cuda+magma+python+fortran
```

### Install on MS Windows platform

Tasmanian has been tested with MS Visual Studio 2017 and CMake 3.11.

* First use the CMake GUI to set the folders and options
* Then use the command prompt (`cmd.exe`) to enter the build folder
```
  cd <cmake-build-folder>
  cmake --build . --config Release
  ctest -C Release
  cmake --build . --config Release --target install
```
* Both Debug and Release are supported config modes, but do not use them simultaneously,
  pick only one Release or Debug.

### Install folder structure

Tasmanian follows standard Linux conventions, the install path could
be set to `/usr/local/`, although it is recommended to install in
a location inside the user home folder to avoid potential system-wide conflicts.

* Install folder structure:
```
  <install-path>/bin/                       (tagrid executable tool)
  <install-path>/lib/                       (shared and static libraries)
  <install-path>/lib/Tasmanian/             (cmake package-config files)
  <install-path>/lib/pythonX.Y/             (python module)
  <install-path>/include/                   (headers .h and .hpp, and Fortran .mod)
  <install-path>/share/Tasmanian            (bash env scripts, install log, table)
  <install-path>/share/Tasmanian/examples/  (reference examples)
  <install-path>/share/Tasmanian/matlab/    (matlab scripts)
  <install-path>/share/Tasmanian/python/    (copy of <install-path>/lib/pythonX.Y)
```

Additional notes:
* A summary of all compile options is stored in
```
  <install-path>/share/Tasmanian/Tasmanian.log
```
* The log file can be displayed with wither command:
```
  <install-path>/bin/tasgrid -log  (using the command line tool)
  python -m Tasmanian              (only if Python is enabled)
```
* The executable and library paths, as well as the Python path can be set in `bash` by sourcing
```
  source <install-path>/share/Tasmanian/TasmanianENVsetup.sh
```
* When using pip, the python module location may be different depending on the OS and the Python virtual environment.
* The Python module is version independent, i.e., the files work with all tested Python versions; the `share/Tasmanian/python` copy allows for easy and version independent access to the Tasmanian modules.
* Under MS Windows the shared libraries (e.g., the .dll files) are installed in `bin`.

### Linking to Tasmanian: CMake Package Config

Tasmanian will install CMake package-config files in `<install-path>/lib/Tasmanian`, the files will contain all necessary information to import the Tasmanian targets into other CMake projects using the CMake command:
```
  find_package(Tasmanian 7.0 PATHS "<install-path>")
```
See the included `CMakeLists.txt` in `<install-path>/share/Tasmanian/examples`.

Note that the `PATHS` do not have to be specified explicitly
if the `TasmanianENVsetup.sh` is sourced or if the `Tasmanian_ROOT` environment variable is set.
The correct `find_package()` command is displayed in the log (see the previous section).

The imported targets will be named:
```
  Tasmanian::shared     (link to all shared libraries, if shared libraries were build)
  Tasmanian::static     (link to all static libraries, if static libraries were build)
  Tasmanian::Tasmanian  (always available and equivalent to either static or shared)
  Tasmanian::tasgrid    (imported executable pointing to the command line tool)
  Tasmanian::Fortran::shared   (shared libraries for Fortran)
  Tasmanian::Fortran::static   (static libraries for Fortran)
  Tasmanian::Fortran           (equivalent to Tasmanian::Tasmanian but using Fortran)
```
Note that the `Tasmanian::Tasmanian` and `Tasmanian::Fortran` targets are no longer equivalent (as of 7.1).

In addition, the following variables will be set:
```
  Tasmanian_PYTHONPATH         (path to the Python module, if Python was enabled)
  Tasmanian_MATLAB_WORK_FOLDER (path to the MATLAB work folder, if set during build)
  Tasmanian_MATLABPATH         (path to the MATLAB scripts, if MATLAB was enabled)
  Tasmanian_<component>_FOUND  (set to ON for each available component)
```
The possible components are:
```
  SHARED STATIC OPENMP BLAS CUDA MAGMA MPI PYTHON MATLAB FORTRAN
```
The modules correspond to shared and static libraries and the cmake options used during build.

All available components will be included even if the component is not explicitly requested.
Requesting components can alter the behavior of `Tasmanian::Tasmanian`,
and help catch errors early in the build process,
and/or print useful log messages. For example:
```
  find_package(Tasmanian 7.0 REQUIRED SHARED PYTHON CUDA OPTIONAL_COMPONENTS OPENMP)
```
In the above example:
* an error will be generated if Tasmanian was build without shared libraries, CUDA or Python support
* the `Tasmanian::Tasmanian` target will be set to the shared libraries
* a status message will report whether Tasmanian was build with OpenMP support

The `Tasmanian::Tasmanian` target will point to the shared libraries if only shared libraries are available,
or if the `SHARED` component is available and requested without the `STATIC` component.
Otherwise, `Tasmanian::Tasmanian` will point to the `STATIC` libraries.
For example:
```
  # suppose both shared and static libraries are available, then
  find_package(... REQUIRED SHARED)                            # Tasmanian::Tasmanian is shared
  find_package(... OPTIONAL_COMPONENTS SHARED)                 # Tasmanian::Tasmanian is shared
  find_package(... REQUIRED SHARED OPTIONAL_COMPONENTS STATIC) # Tasmanian::Tasmanian is static
  find_package(... <no SHARED/STATIC component specified>)     # Tasmanian::Tasmanian is static
```

### Known Issues

Several known issues and work-around fixes:
* The addon tests sometime fail due to thread scheduling
    * The overhead associated with thread scheduling is much larger than the simple test models used,
    which leads to unrealistically large fluctuations in sample run-time, which in turn leads to
    randomness in the results, most notably on machines with few cpu cores.
    * Rerun the tests and/or installation to see if the problem is persistent
* Mixing the GCC and Clang compilers and linkers sometimes fails with an error about the architecture
    * use shared libraries only, i.e., `-D BUILD_SHARED_LIBS=ON` in CMake
* The PGI compiler fails when using optimization `-O2` with an error about an empty `free()`
    * the bug happens when `std::vector<double>::resize()` is called on an empty vector
    * use the `-O1` instead, this issue needs further investigation
* Older versions of CUDA do not work with newer versions of some compilers, e.g., `gcc`
    * consult the CUDA manual for a list of acceptable compilers
