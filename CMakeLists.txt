cmake_minimum_required(VERSION 3.10)

cmake_policy(VERSION 3.10)
project(Tasmanian VERSION 7.7.0 LANGUAGES CXX)
set(Tasmanian_version_comment " (release candidate)") # e.g., " (release candidate)", " (development)", ""
set(Tasmanian_license "BSD 3-Clause with UT-Battelle disclaimer") # used in some headers and python modules (only human readable)

########################################################################
#  User specified options:
#    -D Tasmanian_ENABLE_RECOMMENDED:BOOL=OFF (includes some flags)
#    -D Tasmanian_ENABLE_OPENMP:BOOL=OFF      (recommended)
#    -D Tasmanian_ENABLE_BLAS:BOOL=OFF        (recommended)
#    -D Tasmanian_ENABLE_PYTHON:BOOL=OFF      (recommended)
#    -D Tasmanian_ENABLE_CUDA:BOOL=OFF        (stable)
#    -D Tasmanian_ENABLE_HIP:BOOL=OFF         (stable)
#    -D Tasmanian_ENABLE_DPCPP:BOOL=OFF       (experimental)
#    -D Tasmanian_ENABLE_MAGMA:BOOL=OFF       (stable)
#    -D Tasmanian_MATLAB_WORK_FOLDER:PATH=""  (stable)
#    -D Tasmanian_ENABLE_FORTRAN:BOOL=OFF     (mostly stable)
#    -D Tasmanian_ENABLE_MPI:BOOL=OFF         (stable)
#    -D Tasmanian_ENABLE_SWIG:BOOL=OFF        (mostly stable)
#
# Additional options: (if default test behavior is not desired)
#    -D Tasmanian_TESTS_OMP_NUM_THREADS:INT="sets OpenMP number of threads"
#    -D Tasmanian_TESTS_GPU_ID:INT="specifies which GPU to use for testing"
#
# Extra options: (allow for workarounds of bugs and exotic systems)
#    -D Tasmanian_EXTRA_LIBRARIES:STRING="appends more link libraries"
#    -D Tasmanian_EXTRA_INCLUDE_DIRS:LIST="appends more include paths"
#    -D Tasmanian_EXTRA_LINK_DIRS:PATH="appends more link paths"
#
# Selecting specific packages for find_package()
#     -D PYTHON_EXECUTABLE:PATH
#     -D CMAKE_CUDA_COMPILER:PATH
#     -D Tasmanian_MAGMA_ROOT:PATH (or just MAGMA_ROOT)
#     -D MPI_CXX_COMPILER
#     -D MPI_Fortran_COMPILER (must be compatible with the MPI_CXX_COMPILER)
#     -D MPIEXEC_EXECUTABLE   (must be compatible with the MPI_CXX_COMPILER)
#     -D SWIG_EXECUTABLE      (path to swigfortran)
#     (note: swig may require extra includes, use Tasmanian_EXTRA_INCLUDE_DIRS)
#
# Directly specifying libraries and bypassing find_package()
#     -D BLAS_LIBRARIES
#     -D LAPACK_LIBRARIES
#
# Note: directly specifying libraries still requires the corresponding
#       Tasmanian_ENABLE_ option
#
########################################################################

########################################################################
# The build process follows these steps:
# 1. Define all options and assign default values
# 2. Perform sanity check
#    - resolve conflicts, e.g., Python requires shared libs
#    - call find_package() for options that need packages
#    - fail or fallback if find_package() fails for some option
#    - add Tasmanian_dependencies target with all TPLs
# 3. Enable testing and set the master config file
# 4. Add all subdirectories in the correct order
#    - SparseGirds, DREAM, Addons
#    - Interfaces: Python, Fortran, MATLAB
#    - master target and examples
# 5. Install top level files
#    - add the "make test_install" target
#    - cmake export file
#    - master config file
#    - CMakeLists.txt for the examples (include OpenMP hack if needed)
#    - shell (bash) source file that defines paths
# 6. Print final messages
########################################################################

option(USE_XSDK_DEFAULTS            "Enable xSDK compatibility"              OFF)
option(BUILD_SHARED_LIBS            "Specify shared or static libraries for Tasmanian"  ON)
option(Tasmanian_ENABLE_RECOMMENDED "Enable (if found) OpenMP, BLAS, and Python, also set optimization flags" OFF)
option(Tasmanian_ENABLE_OPENMP      "Enable OpenMP support for Tasmanian (recommended)"    OFF)
option(Tasmanian_ENABLE_BLAS        "Enable CPU Blas support for Tasmanian (recommended)"  OFF)
option(Tasmanian_ENABLE_PYTHON      "Enable Python interface for Tasmanian (recommended)"  OFF)
option(Tasmanian_ENABLE_CUDA        "Enable Nvidia CUDA kernels and libraries within Tasmanian (stable)"    OFF)
option(Tasmanian_ENABLE_HIP         "Enable AMD HIP kernels and ROCm libraries within Tasmanian (stable)"   OFF)
option(Tasmanian_ENABLE_DPCPP       "Enable Intel DPC++ kernels and libraries within Tasmanian (experimental)"    OFF)
option(Tasmanian_ENABLE_MAGMA       "Enable acceleration using UTK Magma library within Tasmanian (stable)" OFF)
option(Tasmanian_ENABLE_FORTRAN     "Enable Fortran interface for Tasmanian (mostly stable)" OFF)
option(Tasmanian_ENABLE_MPI         "Enable MPI support for Tasmanian (stable)"       OFF)
option(Tasmanian_ENABLE_SWIG        "Regenerate Fortran 2003 bindings with SWIG"       OFF)
option(Tasmanian_ENABLE_DOXYGEN     "Enable Doxygen documentation for Tasmanian (complete)" OFF)

# treat those similar to options, but give values other than ON/OFF
set(Tasmanian_MATLAB_WORK_FOLDER     ""  CACHE PATH   "Enable MATLAB interface and use the path for the MATLAB work folder (stable)")

set(Tasmanian_TESTS_GPU_ID          "-1" CACHE STRING "(integer) specify GPU ID for testing, -1 means running tests on all visible devices (optional)")
set(Tasmanian_TESTS_OMP_NUM_THREADS "-1" CACHE STRING "(integer) specify OMP_NUM_THREADS for the tests, if less than 0, this option will be ignored (optional)")

# allow an external project to specify the export variable name
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(Tasmanian_export_name "Tasmanian-export")
    enable_testing() # Tasmanian is the main project, enable testing
else()
    set(Tasmanian_export_name "Tasmanian-export" CACHE STRING "Specify the export variable name when merging Tasmanian into an external project")
endif()


########################################################################
# Sanity check, xSDK compatibility, and find_package() calls are all
# done in one place for consistency
# after this call, every enabled option comes with proper found set of
# find_package() provided libraries, includes, etc.
# Also adds Tasmanian_dependencies that links to all dependencies
########################################################################
include("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/helper_macros.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/sanity_check_and_xsdk.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/add_dependencies_target.cmake")


########################################################################
# Core project configuration, passes all options to CXX
########################################################################
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianConfig.in.hpp"
               "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.hpp")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.hpp"
        DESTINATION include
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
# configure the headers for the logs
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/SparseGrids/tasgridLogs.in.hpp" "${CMAKE_CURRENT_BINARY_DIR}/configured/tasgridLogs.hpp")


########################################################################
# Setup targets
# CXX subdirs have to come first, SparseGrids must precede DREAM
# the Addons must come last
# each sub-directory creates a set of targets and links those to
# existing targets (hence the order is important)
# After CXX, the interfaces Fortran, Python, MATLAB can come in any order
# Finally, add the master target to link to all of the above
########################################################################
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/SparseGrids" "${CMAKE_CURRENT_BINARY_DIR}/SparseGrids")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/DREAM"       "${CMAKE_CURRENT_BINARY_DIR}/DREAM")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Addons"      "${CMAKE_CURRENT_BINARY_DIR}/Addons")

if (Tasmanian_ENABLE_FORTRAN)
    if (NOT Tasmanian_DISABLE_F90)
        add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfaceFortran"  "${CMAKE_CURRENT_BINARY_DIR}/Fortran")
    endif()
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfaceSwig"
      "${CMAKE_CURRENT_BINARY_DIR}/Fortran03")
endif()

if (Tasmanian_ENABLE_PYTHON)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfacePython"  "${CMAKE_CURRENT_BINARY_DIR}/Python")
endif()

if (NOT "${Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfaceMATLAB"  "${CMAKE_CURRENT_BINARY_DIR}/MATLAB")
endif()

if (Tasmanian_ENABLE_DOXYGEN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Doxygen"  "${CMAKE_CURRENT_BINARY_DIR}/Doxygen")
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/master_targets.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/all_examples.cmake")


########################################################################
# Tasmanian is merged into an external master project then skip the exports
# The master should be responsible for post-install testing and exporting
########################################################################
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    include("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/exports.cmake")
endif()


########################################################################
# Setup Message
########################################################################
include("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/setup_message.cmake")

