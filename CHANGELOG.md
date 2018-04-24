Changelog for version 6.0
--------------

* required cmake 3.5 or newer (as opposed to 2.8 in version 5.1)

* CXX standard 2011 is now enabled by default even if CUDA and MPI
  are disabled; removed the option to force-disable CXX 2011


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
