Changelog for version 7.7
--------------

* added new method for constructing quadrature rules with negative weight functions
    * the weight function needs to be bounded from below
    * the methods are added to the Addons module under the Exotic Quadrature name
    * credit goes to William Kong
* stabilized the DPC++ backend and improved the oneMKL compatibility
    * cleaned the build system and improved the way the dependencies are handled
* stabilized the Fortrain 2003 interface to cover all Sparse Grid functionality
    * GPU and MPI methods are now fully supported through Fortran 2003
    * added Fortran 2003 methods for complex numbers for Fourier grids
    * the old Fortran 90/95 interface is now deprecated
* improved the grid compression methods
    * now can remove points by keeping a specific number of points/basis-functions
* moved the post-install testing to CMake
    * simplifies the spack testing for the E4S compatibility
* improved the libEnsemble-Tasmanian integration, see the documentation
* improved the logic for accepting external handles and queues for all GPUs
* deprecated the simple install script, use CMake directly

Changelog for version 7.5
--------------

* added DPC++ backend using oneMKL
    * added both kernels and calls to BLAS/LAPACK methods
    * the kernels for the global grids is missing atomics (currently very slow)
* updated Rocm interface
    * using the Rocm OpenMP implementation in libiomp5.so and libomp.so
* added github CI covering a MacOSX build and basic Ubuntu build
    * credit goes to Viktor Reshniak
* extended the Fortrain 2003 methods
    * added functions, e.g., returnPoints() as members to TasmanianSparseGrid
    * added MPI methods to the Fortain interface
* updates to the GNU Make build system
    * removed Fortran from the GNU Make build system, use CMake
        * the added complexity to the SWIG build system is too much to maintain
    * GNU Make now defaults to Python 3, removed the make python3 target

Changelog for version 7.3
--------------

* added AMD HIP/ROCm backend
    * all CUDA capabilities have an AMD HIP/ROCm equivalent
    * this includes MAGMA-HIP capabilities
    * capabilities are accessed with accel_gpu_cuda and accel_gpu_cublas
* removed the XSDK names for options
    * the XSDK requirement was alleviated due to spack
* CMake now builds only one type of library shared or static
    * matlab interface works well with shared libs thanks to rpaths
    * magma search is simpler looking for either shared or static libs
    * ::shared and ::static targets are deprecated

Changelog for version 7.1
--------------

* the Python pip installer can now enable CUDA, MAGMA, MPI, etc.
    * requires setting environment variables

* the Python pip installer now handles dependencies automatically
    * requires an up-to-date version of pip

* added the API to manually set the CUDA handles
    * manual handles can be set for cuBlas, cuSparse and cuSolverDn

* new addon method for constructing a surrogate from unstructured data
    * utilizes QR factorization from cuSolver and MAGMA
    * using out-of-core implementation within MAGMA

* an LAPACK implementation is required for the BLAS option
    * BLAS and LAPACK are packaged together in all tested libraries

* added Swig generated Fortran 2003 interface
    * naming conventions and objects mimic C++, similar to python

* split the Tasmanian::Tasmanian and Tasmanian::Fortran targets
    * also added Tasmanian::Fortran::shared and Tasmanian::Fortran::static

* added mixed-precision templates for
    * `evaluateBatchGPU()` and `evaluateHierarchicalFunctionsGPU()`
    * `evaluateBatch()` when the acceleration mode is cuda or magma


Changelog for version 7.0
--------------

* improved MPI capability, see the updated documentation
    * fully automated distributed adaptive sparse grid construction

* new module "Tasmanian Addons" consisting of miscellaneous templates
    * algorithms for automated sampling from a lambda model
    * use through the master target and header, cannot be used independently

* new Tasmanian CMake master targets that handle all modules
    * now available Tasmanian::shared and Tasmanian::static
    * also Tasmanian::Tasmanian will always point to an available target
    * added Tasmanian::tasgrid imported executable target

* new Tasmanian master Python module: `import Tasmanian`
    * the python changes are backwards compatible, the old module is still there
    * the sparse grid class can be called with `grid = Tasmanian.SparseGrid()`
    * the DREAM binding is available through `Tasmanian.DREAM`

* new Tasmanian master header to add all modules: Tasmanian.hpp

* new dynamic construction algorithms have been added
    * removed the need for blocking between refinement iterations
        * Tasmanian can accept data one-point at a time
        * new candidate points can be requested at any time
    * candidate points for dynamic construction are weighted by "importance"
    * the dynamic construction process is available through C++ and Python interfaces

* modernized C++ compatibility
    * see the updated DREAM api notes
    * TasmanianSparseGrid has move and copy constructors and `operator=` overloads
    * Tasmanian C++ API no-longer returns raw-pointers, only STL containers
    * **broke backward api** affects mostly get points and weights
    * externally allocated raw-pointers are still accepted and used by C/Python/Fortran APIs

* updated the Fortran interface, the library no longer has a global state
    * replaced the integer grid-ids with a derived type(TasmanianSparseGrid)
    * grid variables now require tsgAllocateGrid(grid) and tsgDeallocateGrid(grid)
    * tsgAllocateGrid() is now thread safe (when called on different grid variables)
    * no other changes appear on the front-end, see the Fortran examples

* improved the `add_subdirectory()` capability
    * can specify the export name used by the Tasmanian install commands
        * `set(Tasmanian_export_name <name> CACHE INTERNAL "")`
        * `add_subdirectory(<path-to-Tasmanian-source> <work-folder>)`
        * `install(EXPORT ${Tasmanian_export_name} ...)`
    * the export name functionality is required to import the transitive dependencies
    * in addition, when using `add_subdirectory()`:
        * will not enable testing, tests are still set if enabled by the master project
        * will not install package-config, that is the job of the master project
        * will not install example CMakeLists.txt or post-install tests
        * (`make test_install` and examples require Tasmanian package-config)

* updated the DREAM interface:
    * excessive polymorphism is replaced by lambdas
    * sampling is done by a template
    * can specify arbitrary domain (using lambdas)
    * the new API is also available though Python
    * Note: the old API is completely obsolete and incompatible

* added Doxygen documentation (CMake option, extra pages, etc.)

* improved CUDA support
    * added support for CUDA 10
    * all grids now benefit from all acceleration modes
    * added `evaluateBatchGPU()` where both `x` and `y` sit on the GPU

* removed the deprecated MS Windows build system with batch scripts

* require cmake 3.10 or newer
    * removed the clumsy work-around for OpenMP on old cmake systems
    * removed other small legacy cmake fixes
    * added cuda as a cmake lang, added multiple hacks to avoid feature regression


Changelog for version 6.0
--------------

* as always, a new version includes numerous bug fixes and performance enhancements

* CXX standard 2011 is now required and enabled by default

* required cmake 3.5 or newer (as opposed to 2.8 in version 5.1)

* merged Tasmanian_ENABLE_CUBLAS option into Tasmanian_ENABLE_CUDA
    * the option is OFF by default
    * two distinct acceleration modes use this variable: `gpu-cublas` and `gpu-cuda`

* added new acceleration mode `gpu-magma` that uses UTK MAGMA library (requires CUDA)

* removed Tasmanian_STRICT_OPTIONS, now all options are considered strict by default

* new option Tasmanian_ENABLE_RECOMMENDED
    * searches for OpenMP, BLAS, and Python, and enables if found
    * set the `-O3` flag for Debug and Release
    * adjusted the install script, see `./install --help`

* modified install folder structure now everything sits in four places
    * `<prefix>/bin` takes the executable files
    * `<prefix>/lib` takes the libraries and Fortran `.mod` files
    * `<prefix>/lib/Tasmanian` takes the cmake package-config
    * `<prefix>/lib/PythonX.Y` takes the python module
    * `<prefix>/include` takes the headers
    * `<prefix>/share/Tasmanian` takes everything else

* added package-config file, now it is possible to use the command
```
find_package(Tasmanian 6.0 PATHS "<Tasmanian install prefix>")
# PATHS "<Tasmanian install prefix>" not needed if the install folder
# is included in CMAKE_PREFIX_PATH
```

* added new grids and rules
    * `GridFourier` that uses trigonometric basis functions (see Manual)
    * `localpb` rule to local polynomial grids that favors the boundary

* API change for error handling (not backward compatible)
    * errors in C++ now throw exceptions
    * `std::runtime_error` is thrown on wrong file format
    * `std::runtime_error` is thrown when requesting level too high
    * `std::invalid_argument` is thrown when calling function for wrong grid
    * `logstream` is no longer needed and has been removed
    * the error handling in the interface calls remains the same as before

* updated C++ API with overloaded functions
    * using `std::vector<int/double>` as opposed to arrays
    * error checking is provided based on the vector size
    * massive portion of the internal API has also changed

* acceleration with custom CUDA kernels is now available
    * Sequence grids, i.e., build with `makeSequenceGrid()`
    * Fourier grids, i.e., build with `makeFourierGrid()`

* the library is more robust when working with large data sets
    * number of points, outputs, or batch size is still `int`
    * in most cases the product (i.e., matrix size) can exceed `MAX_INT`

* improved MS Windows support using cmake
    * recommended optimization flags are properly set (Debug and Release)
    * shared symbols are exported by cmake
    * CUDA is now supported under Windows
    * `getGPUName()` under Python currently does not work
    * the old `WindowsMake.bat` script is deprecated


Changelog for version 5.1
--------------

* visit us on github.com/ORNL/TASMANIAN

* added functionality to access hierarchical functions and coefficients
  which allows constructing model approximation from an arbitrary
  cloud of samples, as opposed to requiring the use of the specific
  points in the hypercube

* added custom CUDA kernels for evaluations for LocalPolynomial grids,
  both sparse/dense basis of hierarchical functions and batch evaluations
  can be performed with the CUDA kernels for extra speed.
  evaluateFast for LocalPolynomial also benefits from acceleration

* grids can be written/read in either binary or ascii format,
  MATLAB defaults to binary format now

* CUDA accelerated evaluations can be triggered from MATLAB,
  just specify `lGrid.gpuDevice = X,` where X is a valid CUDA device
  run `tsgCoreTests()` to see a list of the available devices

* added tsgCoreTests() script to MATLAB

* added tweaks and simplifications to the build system,
  the main options remain the same, but check the manual for changes
  to the additional options

* added Fortran 90/95 interface (still very EXPERIMENTAL)

* significantly expanded the automated testing procedures and code
  coverage

* numerous bug fixes and performance enhancements


Changelog for version 5.0
--------------

* added initial implementation of the DiffeRential Evolution Adaptive
  Metropolis (DREAM) capabilities (C++ only, for now)

* batch evaluation in Global and Sequence grids can be accelerated using
  BLAS and Nvidia CUDA (cuBlas). Useful for interpolants with large
  number of outputs.

* overhaul of the folder structure and build system, consult the updated
  manual for all included build options.

- Install script: install.sh (in v5.1 renamed to install) using cmake, compile the code in TasmanianSource/Build
~~~
     cd <tasmanian source folder>
     ./install.sh -help
     ./install.sh <install prefix> <matlab work folder>
~~~
- Manual cmake: !use out-of-source build!
~~~
     cd <build folder>
     cmake <options> <tasmanian source folder>
     make
     make test
     make install
     ./test_post_install.sh (replaced by "make test_install" in v5.1)
~~~

- GNU make
~~~
     make
     make test
     make examples
     ./example_sparse_grids
     ./example_dream
     make matlab
~~~


Changelog for version 4.0
--------------

* the build scrips compile both shared (dynamic) and static libraries

* added C-interface, the library compiles C++ and C functions and
  TasmanianSparseGrids.h header defines the C interface that is based
  on a void-pointer that holds an instance of the C++ class

* added Python module interface based on ctypes and the C interface.
  The Python module is called TasmanianSG.py and it defines the TasmanianSparseGrids
  python class which mimics closely the C++ class with enumerates replaced by strings
  The Python module raises TasmanianInputError exception on incorrect input

* support for cmake, starting with version 4.0 cmake is the preffered build system
  the simple GNU-Make will continue to be supported.

  Example cmake call (Linux): (!DOES NOT WORK ON VERSION 5.x, SEE MANUAL)
~~~
  cd /home/myself/TasmanianBuild
  cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/home/myself/TSG \
    -DENABLE_PYTHON=ON \
    -DENABLE_MATLAB=ON \
    -DMATLAB_WORK_FOLDER=/home/myself/tsgWorkFiles/ \
    -DENABLE_OPENMP=ON \
    /home/myself/TasmanianSparseGrids/
  make
  make test
  make install
~~~

* in the above example, add the following to .bashrc
  in order to use tasgrid executable:
    `export PATH=$PATH:/home/myself/TSG/bin`

  in order to use tasgrid shared library:
    `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/myself/TSG/lib`

  in order to compile C or C++ code:
    `export C_INCLUDE_PATH=$C_INCLUDE_PATH:/home/myself/TSG/include`
    `export CPP_INCLUDE_PATH=$CPP_INCLUDE_PATH:/home/myself/TSG/include`

  alternatively you can add -I/home/myself/TSG/include to your build command

* in the above example, the python module will be installed under
    `/home/myself/TSG/python/TasmanianSG.py`
    in your python script use:

~~~
    import sys
    sys.path.append("/home/myself/TSG/python")
    import TasmanianSG
~~~

* in the above example, the MATLAB files will be installed under
    /home/myself/TSG/matlab/

    in MATLAB use:
    `addpath("/home/myself/TSG/matlab/")`

    or add the path using the MATLAB gui

* added piece-wise constant polynomial rule

* added functions to directly evaluate hierarhical functions for piece-wise
  polynomial, wavelet, and sequence grids. Also, the hierarchical surpluses
  associated with those functions can be overwritten. This functionality is
  useful in constructing a grid that is the result of a projection of an
  arbitrary data set, as opposed to interpolaiton at pre-defined nodes.
  However, the functionality currently intrudes into the data structures and
  is hence to fully complient with the rest of the workflow and codeflow, thus
  the functions are to be considered "experimental" and to be used with lots of
  caution.

* fixed various small bugs


Changelog for version 3.1
--------------

* fixed a major bug when using the tasgrid executable, `TasmanianSparseGrid::write`, and `TasmanianSparseGrid::read` funcitons
  the custom domain transform was not saved properly
  using grid files saved with version 3.0 will generate some warning, but it will work fine otherwise (except for the custom transfrom)

* added example CC, COMPILE_OPTIONS, and OPTC Makefile directives for clang + OpenMP as well as Intel icc compilers

* added a windows bat script that mimics the make command but is compatible with Microsoft Visual C++

* changed the signature of
~~~
   TasmanianSparseGrid::setAnisotropicRefinement
   TasmanianSparseGrid::estimateAnisotropicCoefficients
   TasmanianSparseGrid::setSurplusRefinement
   TasmanianSparseGrid::setSurplusRefinement
~~~
  All functions accept an output that specifies which output to be used
  Sequence, Local Polynomial, and Wavelet grids can now refine based on a specified output as opposed to all output
  to regain the old behavior use `output = -1`
