cmake_minimum_required(VERSION 2.8)
if (${CMAKE_MAJOR_VERSION} LESS 3)
    project(Tasmanian)
    set(Tasmanian_VERSION_MAJOR 5)
    set(Tasmanian_VERSION_MINOR 1)
    set(Tasmanian_VERSION_PATCH 0)
else()
    cmake_policy(VERSION 3.0.2)
    project(Tasmanian VERSION 5.1.0 LANGUAGES CXX)
endif()

########################################################################
#  User specified options:
#    -D Tasmanian_STRICT_OPTIONS:BOOL=OFF
#    -D Tasmanian_ENABLE_OPENMP:BOOL=ON
#    -D Tasmanian_ENABLE_BLAS:BOOL=ON
#    -D Tasmanian_ENABLE_MPI:BOOL=OFF
#    -D Tasmanian_ENABLE_CUBLAS:BOOL=ON
#    -D Tasmanian_ENABLE_CUDA:BOOL=OFF
#    -D Tasmanian_ENABLE_MATLAB:BOOL=OFF
#    -D Tasmanian_ENABLE_PYTHON:BOOL=ON
#    -D Tasmanian_SHARED_LIBRARY:BOOL=ON
#    -D Tasmanian_STATIC_LIBRARY:BOOL=ON
#    -D Tasmanian_MATLAB_WORK_FOLDER:PATH=""
#
# Additional options:
#    -D Tasmanian_TESTS_OMP_NUM_THREADS:INT="sets OpenMP number of threads"
#    -D Tasmanian_EXTRA_CXX_FLAGS:STRING="appends to the default flags"
#    -D Tasmanian_EXTRA_INCLUDE_DIRS:PATH="appends more include paths"
#    -D Tasmanian_EXTRA_LIBRARIES:STRING="appends more link libraries"
#    -D Tasmanian_EXTRA_LINK_DIRS:PATH="appends more link paths"
#    -D Tasmanian_ENABLE_CXX_2011:BOOL=OFF
#      (note: unless manually specified, CUDA, CUBLAS, and MPI set this to ON)
#
# Development only options:
#    -D Tasmanian_ENABLE_DEVELOPMENT_DEFAULTS:BOOL=OFF
#    -D Tasmanian_DEVELOPMENT_BACKWARDS:BOOL=OFF
#
# xSDK options:
#    -D USE_XSDK_DEFAULTS:BOOL=OFF
#    -D XSDK_ENABLE_OPENMP:BOOL=OFF
#    -D TPL_ENABLE_BLAS:BOOL=OFF
#    -D XSDK_ENABLE_PYTHON:BOOL=OFF
#    -D XSDK_ENABLE_FORTRAN:BOOL=OFF
#    -D XSDK_ENABLE_CUDA:BOOL=OFF
#
# Alternatives to find_package()
#     -D PYTHON_EXECUTABLE:PATH
#     -D CUDA_TOOLKIT_ROOT_DIR:PATH
#     -D OpenMP_CXX_FLAGS:STRING
#     -D BLAS_LIBRARIES:STRING
#     -D CUDA_CUBLAS_LIBRARIES:STRING
#     -D CUDA_cusparse_LIBRARY:STRING
#     -D MPI_CXX_LIBRARIES:STRING
#
# PYTHON_EXECUTABLE and CUDA_TOOLKIT_ROOT_DIR help cmake find a specific
#     implementation, but find_package() is still used
# OpenMP_CXX_FLAGS bypasses find_package(OpenMP)
#
# BLAS_LIBRARIES bypasses find_package(BLAS)
# MPI_CXX_LIBRARIES bypasses find_package(MPI)
# CUDA_CUBLAS_LIBRARIES bypasses find_package(CUDA)
#   but only if Tasmanian_ENABLE_CUDA=OFF,
#   Tasmanian_ENABLE_CUDA=ON absolutely requires find_package(CUDA)
#
# Optional to be added to the alternatives above
#     -D MPI_CXX_INCLUDE_PATH
#     -D MPI_COMPILE_FLAGS
#     -D MPI_LINK_FLAGS
#     -D CUDA_INCLUDE_DIR
#     -D CUDA_LIBRARIES
#     -D CUDA_cusparse_LIBRARY
#
########################################################################
option(USE_XSDK_DEFAULTS "Enable xSDK compatibility" OFF)
option(Tasmanian_STRICT_OPTIONS "Disable automatic adjustment of the build options for Tasmanian" OFF)
option(Tasmanian_ENABLE_OPENMP "Enable OpenMP support for Tasmanian" ON)
option(Tasmanian_ENABLE_BLAS "Enable CPU Blas support for Tasmanian" ON)
option(Tasmanian_ENABLE_MPI "Enable MPI support for Tasmanian" OFF)
option(Tasmanian_ENABLE_MATLAB "Enable MATLAB interface for Tasmanian" OFF)
option(Tasmanian_ENABLE_PYTHON "Enable Python interface for Tasmanian" ON)
option(Tasmanian_ENABLE_FORTRAN "Enable Fortran interface for Tasmanian" OFF)
option(Tasmanian_ENABLE_CUBLAS "Enable Nvidia GPU acceleration using cuBLAS for Tasmanian" ON)
option(Tasmanian_ENABLE_CUDA "Enable Nvidia GPU acceleration using CUDA kernels for Tasmanian" OFF)
option(Tasmanian_STATIC_LIBRARY "Compile static library and the executables with static linkage (Tasmanian)" ON)
option(Tasmanian_SHARED_LIBRARY "Compile dynamic library (Tasmanian)" ON)

# version and license
set(Tasmanian_version_string "${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}")
set(Tasmanian_license "BSD 3-Clause with UT-Battelle disclamer")


########################################################################
# Defines the source files
# NOTE: sanity check modifies the CUDA source lists before setup target
########################################################################
include(${CMAKE_SOURCE_DIR}/Config/CMakeIncludes/source_file_list.cmake)


########################################################################
# Sanity check, xSDK compatibility, and find_package() calls
########################################################################
include(${CMAKE_SOURCE_DIR}/Config/CMakeIncludes/sanity_check_and_xsdk.cmake)


########################################################################
# Setup targets
########################################################################
include(${CMAKE_SOURCE_DIR}/Config/CMakeIncludes/target_setup.cmake)


########################################################################
# Setup options, BLAS, OpenMP, MPI, etc.
########################################################################
include(${CMAKE_SOURCE_DIR}/Config/CMakeIncludes/option_setup.cmake)


########################################################################
# Setup compiler flags and extra libs
########################################################################
include(${CMAKE_SOURCE_DIR}/Config/CMakeIncludes/compiler_flags_extras.cmake)


########################################################################
# Configure various files
########################################################################
include(${CMAKE_SOURCE_DIR}/Config/CMakeIncludes/configure_files.cmake)


########################################################################
# Testing
########################################################################
enable_testing()
add_test(SparseGridsAcceleration tasgrid -test acceleration)
add_test(SparseGridsDomain tasgrid -test domain)
add_test(SparseGridsRefinement tasgrid -test refinement)
add_test(SparseGridsGlobal tasgrid -test global)
add_test(SparseGridsLocal tasgrid -test local)
add_test(SparseGridsWavelet tasgrid -test wavelet)
add_test(DREAM_Analytic tasdream -test analytic)
add_test(DREAM_Model tasdream -test model)
set(Tasmanian_tests_list SparseGrids_Acceleration SparseGrids_Domain SparseGrids_Refinement
                         SparseGrids_Global SparseGrids_Local SparseGrids_Wavelet
                         DREAM_Analytic DREAM_Model)
if (Tasmanian_ENABLE_PYTHON)
    add_test(NAME GenericComprehensive COMMAND ${PYTHON_EXECUTABLE} testTSG.py)
    set(Tasmanian_tests_list ${Tasmanian_tests_list} LibSG)
endif()
if (Tasmanian_ENABLE_FORTRAN)
    add_test(Fortran fortester)
    set(Tasmanian_tests_list ${Tasmanian_tests_list} Fortran)
endif()
if (Tasmanian_TESTS_OMP_NUM_THREADS)
    set_tests_properties(${Tasmanian_tests_list} PROPERTIES
        PROCESSORS ${Tasmanian_TESTS_OMP_NUM_THREADS}
        ENVIRONMENT OMP_NUM_THREADS=${Tasmanian_TESTS_OMP_NUM_THREADS})
endif()
add_custom_target(test_install COMMAND ${CMAKE_BINARY_DIR}/test_post_install.sh -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER})



########################################################################
# Setup Message
########################################################################
message("")
message("TASMANIAN: summary of all compile options")
if (NOT USE_XSDK_DEFAULTS)
    message(" -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
    foreach(Tasmanian_option Tasmanian_ENABLE_OPENMP Tasmanian_ENABLE_BLAS Tasmanian_ENABLE_MPI
                             Tasmanian_ENABLE_MATLAB Tasmanian_ENABLE_PYTHON
                             Tasmanian_ENABLE_CUBLAS Tasmanian_ENABLE_CUDA
                             Tasmanian_ENABLE_FORTRAN
                             Tasmanian_SHARED_LIBRARY Tasmanian_STATIC_LIBRARY)
        if (${Tasmanian_option})
            message(" -D ${Tasmanian_option}:BOOL=ON")
        else()
            message(" -D ${Tasmanian_option}:BOOL=OFF")
        endif()
    endforeach()
    if (Tasmanian_ENABLE_MATLAB)
        message(" -D Tasmanian_MATLAB_WORK_FOLDER:PATH=${Tasmanian_MATLAB_WORK_FOLDER}")
        message(" pre-install MATLAB folder: addpath('${CMAKE_BINARY_DIR}/matlab/')")
    endif()
endif()
message("")


########################################################################
# installation
########################################################################
include(${CMAKE_SOURCE_DIR}/Config/CMakeIncludes/install.cmake)
