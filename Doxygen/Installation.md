# Installation

### Requirements

The minimum requirements to use Tasmanian are a C/C++ compilers and either CMake or GNU Make build engines. Additionally, we recommend Python with NumPy and CTypes packages, a BAsic Linear Algebra Subroutine (BLAS) implementation, and an OpenMP implementation (usually included in the corresponding compiler). Optionally, Tasmanian provides acceleration using Nvidia CUDA libraries with custom kernels, as well as basic Python matplotlib support. MATLAB/Octave wrappers around the command-line tool can be enabled, and a Fortran 90/95 interfaces is available, with the gfortran and ifort compilers corresponding to gcc and icc.

| Feature | Tested versions     | Recommended      |
|----|----|----|
| gcc     | 5, 6, 7, 8          | any              |
| clang   | 4, 5, 6             | 4 or 5           |
| icc     | 18.0                | 18.0             |
| cmake   | 3.10, 3.11, 3.12    | 3.10             |
| python  | 2.7, 3.5, 3.6       | 3.5 or 3.6       |
| OpenBlas| 0.2.18, 0.2.20      | 0.2.18 or 0.2.20 |
| CUDA    | 9.0, 9.1, 9.2, 10.0 | 9.1 or 9.2       |
| libiomp | 5.0                 | 5.0              |
| MAGMA   | 2.3, 2.4            | 2.4              |

### Using CMake: the preferred way

The preferred way to install Tasmanian is to use the included CMake build script, which requires CMake 3.10
that can be obtained from [https://cmake.org/](https://cmake.org/). Note that CMake is usually included
in the native package repositories of most Linux distributions (e.g., Ubuntu 18.04) as well as Mac OSX Homebrew.

* The commands for out-of-source CMake build are:
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
  -D CMAKE_INSTALL_PREFIX:PATH=<install-path> (install folder for the make install command)
  -D CMAKE_BUILD_TYPE:STRING=<Debug/Release>  (set debug flags or default optimization flags)
  -D BUILD_SHARED_LIBS:BOOL=<ON/OFF>          (if defined, build only shared/static libraries, otherwise build both)
  -D CMAKE_CXX_COMPILER:PATH=<path>           (specify the C++ compiler)
  -D CMAKE_CXX_FLAGS:STING=<flags>            (set additional flags)
```

* List of Tasmanian specific CMake options:
```
  -D Tasmanian_ENABLE_OPENMP:BOOL=OFF       (recommended)
  -D Tasmanian_ENABLE_BLAS:BOOL=OFF         (recommended)
  -D Tasmanian_ENABLE_PYTHON:BOOL=OFF       (recommended)
  -D Tasmanian_ENABLE_RECOMMENDED:BOOL=OFF  (enable the above and the -O3 flag)
  -D Tasmanian_ENABLE_CUDA:BOOL=OFF         (stable)
  -D Tasmanian_ENABLE_MAGMA:BOOL=OFF        (stable)
  -D Tasmanian_MATLAB_WORK_FOLDER:PATH=""   (stable)
  -D Tasmanian_ENABLE_FORTRAN:BOOL=OFF      (mostly stable)
  -D Tasmanian_ENABLE_MPI:BOOL=OFF          (experimental)
  -D Tasmanian_ENABLE_DOXYGEN:BOOL=OFF      (work in progress)
```
* Acceleration options:
    * OpenMP allows Tasmanian to use more than one CPU core, which greatly increases performance
    * Basic Linear Algebra Subroutines (BLAS) is a standard with many implementations,
      i.e., [https://www.openblas.net/](https://www.openblas.net/); optimized BLAS greatly improves
      performance when evaluating grid at many points simultaneously or working with models with many outputs
    * CUDA is a C++ language extension that allows programming of Nvidia GPU devices, which greatly enhances
      performance of `evaluateFast()` and `evaluateBatch()` calls
    * Matrix Algebra on GPU and Multicore Architectures (MAGMA) is a library for GPU accelerated linear
      algebra developed at University of Tennessee at Knoxville
    * MPI allows the use of distributed memory in Bayesian inference
* The `Tasmanian_ENABLE_RECOMMENDED` option searches for OpenMP, BLAS, and Python, enables the options (if possible) and also sets `-O3` flag
* Additional interfaces are available, beyond the default C/C++ library and Command Line tools:
    * Python and Fortran require appropriate interpreter and compiler
    * The MATLAB/Octave interface requires a work folder with read/write permission for temporary files.
      The interface is enabled by setting `Tasmanian_MATLAB_WORK_FOLDER` to a valid read/write location.
* Doxygen [http://www.doxygen.org/](http://www.doxygen.org/) will build the HTML documentation

* Options to adjust the testing environment, the default behavior is to use the system provided OpenMP parameters and to run tests on all visible GPU devices; this behavior can be modified and the commands below affect only the ctest environment in `make test`, the Tasmanian library can use any number of threads or Nvidia devices regardless of the test options:
```
 -D Tasmanian_TESTS_OMP_NUM_THREADS=<number-of-threads-for-testing>
 -D Tasmanian_TESTS_GPU_ID=<cuda-device-id-for-testing>
```

* Additional commands to guide the CMake `find_package()` modules:
```
  -D PYTHON_EXECUTABLE:PATH         (specify path to Python)
  -D CUDA_TOOLKIT_ROOT_DIR:PATH     (specify path to CUDA)
  -D CMAKE_Fortran_COMPILER:PATH    (specify Fortran compiler)
  -D Tasmanian_MAGMA_ROOT_DIR:PATH  (specify path to MAGMA installation)
```

* Alternatives to directly specify libraries and bypass `find_package()` altogether:
```
  -D BLAS_LIBRARIES
  -D Tasmanian_MAGMA_LIBRARIES
  -D Tasmanian_MAGMA_INCLUDE_DIRS
  -D MPI_CXX_LIBRARIES
  -D MPI_CXX_INCLUDE_PATH
  -D MPI_COMPILE_FLAGS
  -D MPI_LINK_FLAGS
```

* Extra options are available in case CMake fails to find a required dependency, e.g., `find_package()` sometimes fails to acknowledge that the ACML implementation of BLAS depends on both `libgfortran` and `libgomp`, then manual options are available but those should not be necessary in the majority of cases:
```
  -D Tasmanian_EXTRA_LIBRARIES:STRING   (add more libraries as dependencies)
  -D Tasmanian_EXTRA_INCLUDE_DIRS:PATH  (add more include paths to search for headers)
  -D Tasmanian_EXTRA_LINK_DIRS:PATH     (appends more link paths to search for libraries)
```

### Using the `install` script-wrapper around CMake

Tasmanian also includes an `install` script that wraps around CMake and automatically calls the build commands.
The script uses `Tasmanian_ENABLE_RECOMMENDED` option by default and allows for other options to be enabled/disabled
with command line switches.

* Basic usage of the `install` script
```
  ./install <install-path> <optional: matlab-work-folder> <additional options>
  ./install --help  (lists all options)
```
* Example call that enables MATLAB and CUDA
```
  ./install /home/<user>/Tasmanian /home/<user>/Tasmanian/WorkFolder -cuda=/usr/local/cuda-9.2
```
* Additional notes:
    * the script must be called from the root of the source code
    * using absolute paths is strongly recommended
    * if the MATLAB work-folder is omitted, the MATLAB interface will be disabled

### Alternative (basic) build using GNU Make

The core capabilities of Tasmanian can be build with a few simple GNU Make commands.
The basic build engine is useful for quick testing and exploring option, or
if CMake is unavailable or unwanted.
Acceleration options other than OpenMP are not supported in the basic mode, but could be enabled
by manually editing `Config/AltBuildSystems/Makefile.in` configuration file.

* Using GNU Make with `g++` and optionally `gfortran` and `/usr/bin/env python`
```
  make
  make test     (will fail if /usr/bin/env python is missing numpy or ctypes modules)
  make matlab   (optional: sets matlab work folder to ./tsgMatlabWorkFolder/)
  make python3  (optional: sets #!/usr/bin/env python3 in place of /usr/bin/env python)
  make fortran  (optional: compile Fortran libraries)
  make examples
  make clean
```

### Using Spack

Tasmanian is also included in Spack: [https://spack.io/](https://spack.io/)
```
  spack install tasmanian@develop+blas+python+cuda
```

### MS Windows Installation

Tasmanian has been tested with MS Visual Studio 2015 and 2017 and CMake 3.11.

* First: use the CMake GUI to set the folders and options
* Second: use the command prompt (`cmd.exe`) to enter the build folder
```
  cd <cmake-build-folder>
  cmake --build . --config Release
  ctest -C Release
  cmake --build . --config Release --target install
```
* Both Debug and Release are supported config modes above

### Install folder structure

Tasmanian follows standard Linux conventions, the install path could potentially
be set to `/usr/local/`, although it is recommended to install in
a location inside the user home folder to avoid potential system-wide conflicts.
Note that the MATLAB work folder requires read and write privileges.

* Install folder structure:
```
  <install-path>/bin/                       (tagrid and tasdream executables)
  <install-path>/lib/                       (shared and static libraries)
  <install-path>/lib/Tasmanian/             (cmake package-config files)
  <install-path>/lib/pythonX.Y/             (python module)
  <install-path>/include/                   (headers .h and .hpp, and Fortran .mod)
  <install-path>/share/Tasmanian            (bash env scripts, install log, table)
  <install-path>/share/Tasmanian/examples/  (reference examples)
  <install-path>/share/Tasmanian/matlab/    (matlab scripts)
  <install-path>/share/Tasmanian/python/    (sym-link to <install-path>/lib/pythonX.Y)
```

Additional notes:
* A summary of all compile options is stored in `<install-path>/share/Tasmanian/Tasmanian.log`
* The executable and library paths, as well as Python path can be set in `bash` by sourcing `<install-path>/share/Tasmanian/TasmanianENVsetup.sh`
* The Python module is version independent, i.e., the same file works with all tested versions, the version independent sym-link `share/Tasmanian/python` allows to use the Python interface regardless of the version of Python used during the install
* Under MS Windows the shared library (e.g., the .dll files) are installed in `bin` and sym-links are not supported

### CMake Package Config

Tasmanian will install CMake package-config files in `<install-path>/lib/Tasmanian`, the files will contain all necessary information to import Tasmanian targets into another CMake project using the CMake command:
```
 find_package(Tasmanian 6.1 PATHS "<install-path>")
```
The imported targets will be called:
```
  Tasmanian_libsparsegrid_shared  Tasmanian_libsparsegrid_static
  Tasmanian_libdream_shared       Tasmanian_libdream_static
  Tasmanian_libfortran90_shared   Tasmanian_libfortran90_static
```
but depending on the install options, not all targets may be available, e.g., using `-D BUILD_SHARED_LIBS=ON` will disable the static targets. In order to simplify the user code, the package-config will also create CMake interface targets without the shared/static suffixes:
```
  Tasmanian_libsparsegrid    Tasmanian_libdream    Tasmanian_libfortran90
```
The interface targets will always depend on valid imported targets, and if both shared and static libraries are present, the static libraries will be chosen by default.

See `<install-path>/share/Tasmanian/examples/CMakeLists.txt` for an example of how to use the Tasmanian package-config.
