cmake_minimum_required(VERSION 3.5)

cmake_policy(VERSION 3.5)
project(Tasmanian VERSION 6.0.0 LANGUAGES CXX)
message(STATUS "This is a development version of Tasmanian")
set(Tasmanian_license "BSD 3-Clause with UT-Battelle disclaimer") # used in some headers and python modules (only human readable)

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
#    -D Tasmanian_ENABLE_FORTRAN:BOOL=OFF
#    -D Tasmanian_SHARED_LIBRARY:BOOL=ON
#    -D Tasmanian_STATIC_LIBRARY:BOOL=ON
#    -D Tasmanian_MATLAB_WORK_FOLDER:PATH=""
#
# Additional options:
#    -D Tasmanian_TESTS_OMP_NUM_THREADS:INT="sets OpenMP number of threads"
#    -D Tasmanian_TESTS_GPU_ID:INT="specifies which GPU to use for testing"
#    -D Tasmanian_EXTRA_CXX_FLAGS:STRING="appends to the default flags"
#    -D Tasmanian_EXTRA_INCLUDE_DIRS:PATH="appends more include paths"
#    -D Tasmanian_EXTRA_LIBRARIES:STRING="appends more link libraries"
#    -D Tasmanian_EXTRA_LINK_DIRS:PATH="appends more link paths"
#
# Development only options:
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
#     -D MPI_CXX_LIBRARIES:STRING
#
# PYTHON_EXECUTABLE and CUDA_TOOLKIT_ROOT_DIR help cmake find a specific
#     implementation, but find_package() is still used
# OpenMP_CXX_FLAGS bypasses find_package(OpenMP)
#
# BLAS_LIBRARIES bypasses find_package(BLAS)
# MPI_CXX_LIBRARIES bypasses find_package(MPI)
# CUDA_CUBLAS_LIBRARIES bypasses find_package(CUDA),
#   but only if Tasmanian_ENABLE_CUDA=OFF,
#   Tasmanian_ENABLE_CUDA=ON absolutely requires find_package(CUDA)
#
# Optional that can be be added to the alternatives above
#     -D OpenMP_CXX_LIBRARIES:STRING
#     -D MPI_CXX_INCLUDE_PATH
#     -D MPI_COMPILE_FLAGS
#     -D MPI_LINK_FLAGS
#     -D CUDA_INCLUDE_DIRS
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

# threat those similar to options, but give values other than ON/OFF
set(Tasmanian_TESTS_GPU_ID "-1" CACHE STRING "specify GPU ID for testing, -1 means running tests on all visible devices")
set(Tasmanian_TESTS_OMP_NUM_THREADS "-1" CACHE STRING "specify OMP_NUM_THREADS for the tests, if less than 0, this option will be ignored")
set(Tasmanian_MATLAB_WORK_FOLDER "" CACHE PATH "specify the path for the MATLAB work folder")


########################################################################
# Sanity check, xSDK compatibility, and find_package() calls are all
# done in one place for consistency
# after this call, every enabled option comes with proper found set of
# find_package() provided libraries, includes, etc.
########################################################################
include(${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/sanity_check_and_xsdk.cmake)


########################################################################
# Core project configuration
########################################################################
enable_testing() # call in one place, each subfolder gives additional tests

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianConfig.in.hpp"
               "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.hpp")


########################################################################
# Setup targets
# CXX subdirs have to come first, SparseGrids must preceed DREAM
# The order of the interfaces is not important
# each sub-directory creates a set of targets and links those to
# existing targets (hence the order is important)
########################################################################
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/SparseGrids" "${CMAKE_CURRENT_BINARY_DIR}/SparseGrids")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/DREAM"       "${CMAKE_CURRENT_BINARY_DIR}/DREAM")

if (Tasmanian_ENABLE_FORTRAN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfaceFortran"  "${CMAKE_CURRENT_BINARY_DIR}/Fortran")
endif()

if (Tasmanian_ENABLE_PYTHON)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfacePython"  "${CMAKE_CURRENT_BINARY_DIR}/Python")
endif()

if (Tasmanian_ENABLE_MATLAB)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfaceMATLAB"  "${CMAKE_CURRENT_BINARY_DIR}/MATLAB")
endif()


########################################################################
# Testing: post install, make test_install
# checks if the executables can be run and if the examples compile
########################################################################
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Testing/test_post_install.in.sh" "${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh")
if (Tasmanian_ENABLE_FORTRAN)
    add_custom_target(test_install COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
else()
    add_custom_target(test_install COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
endif()


########################################################################
# Setup Message
########################################################################
if (NOT USE_XSDK_DEFAULTS)
    message("")
    message("TASMANIAN: summary of all compile options")
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
        message(" pre-install MATLAB folder: addpath('${CMAKE_CURRENT_BINARY_DIR}/matlab/')")
    endif()
    message("")
endif()


########################################################################
# Install exports, config header, examples cmake
########################################################################
install(EXPORT Tasmanian_exports DESTINATION "config" FILE "Tasmanian.cmake")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.hpp"
        DESTINATION include
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

# cmake file for the examples, to be used post-install
# see the OpenMP remark in the SparseGirds sub-directory
set(Tasmanian_openmp_hack "# ALL GOOD, NO EXTRA OPENMP FLAG NEEDED HERE")
if (Tasmanian_ENABLE_OPENMP AND (NOT (DEFINED OpenMP_CXX_LIBRARIES)))
    set(Tasmanian_openmp_hack "set(CMAKE_CXX_FLAGS \"\$\{CMAKE_CXX_FLAGS\} ${OpenMP_CXX_FLAGS}\")")
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeLists.examples.txt" "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt"
        DESTINATION "examples"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

# configure environment shell script that can be sourced to set path and lib path
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianENVsetup.in.sh" "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh"
        DESTINATION "config"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)


########################################################################
# Print final message
# The message is written in CMakeLists.txt as the subdir is added last
# this ensures that the message will appear last in the install process
########################################################################
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/")
