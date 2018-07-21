cmake_minimum_required(VERSION 3.5)

cmake_policy(VERSION 3.5)
project(Tasmanian VERSION 6.0.0 LANGUAGES CXX)
message(STATUS "This is a development version of Tasmanian")
set(Tasmanian_license "BSD 3-Clause with UT-Battelle disclaimer") # used in some headers and python modules (only human readable)

########################################################################
#  User specified options:
#    -D Tasmanian_ENABLE_RECOMMENDED:BOOL=OFF (includes some flags)
#    -D Tasmanian_ENABLE_OPENMP:BOOL=OFF      (recommended)
#    -D Tasmanian_ENABLE_BLAS:BOOL=OFF        (recommended)
#    -D Tasmanian_ENABLE_PYTHON:BOOL=OFF      (recommended)
#    -D Tasmanian_ENABLE_CUDA:BOOL=OFF        (stable)
#    -D Tasmanian_ENABLE_MAGMA:BOOL=OFF       (stable)
#    -D Tasmanian_MATLAB_WORK_FOLDER:PATH=""  (stable)
#    -D Tasmanian_ENABLE_FORTRAN:BOOL=OFF     (experimental)
#    -D Tasmanian_ENABLE_MPI:BOOL=OFF         (experimental)
#
# Additional options: (if default test behavior is not desired)
#    -D Tasmanian_TESTS_OMP_NUM_THREADS:INT="sets OpenMP number of threads"
#    -D Tasmanian_TESTS_GPU_ID:INT="specifies which GPU to use for testing"
#
# Extra options: (allow for workarounds of bugs and exotic systems)
#    -D Tasmanian_EXTRA_LIBRARIES:STRING="appends more link libraries"
#    -D Tasmanian_EXTRA_INCLUDE_DIRS:PATH="appends more include paths"
#    -D Tasmanian_EXTRA_LINK_DIRS:PATH="appends more link paths"
#
# Development only options: (only developers should use this)
#    -D Tasmanian_DEVELOPMENT_BACKWARDS:BOOL=OFF
#
# xSDK options: (USE_XSDK_DEFAULTS must be ON to use the other XSDK options)
#    -D USE_XSDK_DEFAULTS:BOOL=OFF
#    -D XSDK_ENABLE_OPENMP:BOOL=OFF
#    -D TPL_ENABLE_BLAS:BOOL=OFF
#    -D TPL_ENABLE_MAGMA:BOOL=OFF
#    -D XSDK_ENABLE_PYTHON:BOOL=OFF
#    -D XSDK_ENABLE_FORTRAN:BOOL=OFF
#    -D XSDK_ENABLE_CUDA:BOOL=OFF
#
# Selecting specific packages for find_package()
#     -D PYTHON_EXECUTABLE:PATH
#     -D CUDA_TOOLKIT_ROOT_DIR:PATH
#     -D Tasmanian_MAGMA_ROOT_DIR:PATH
#
# Directly specifying libraries and bypassing find_package()
#     -D BLAS_LIBRARIES
#     -D Tasmanian_MAGMA_LIBRARIES
#     -D Tasmanian_MAGMA_INCLUDE_DIRS
#     -D MPI_CXX_LIBRARIES
#     -D MPI_CXX_INCLUDE_PATH
#     -D MPI_COMPILE_FLAGS
#     -D MPI_LINK_FLAGS
#
# Note: directly specifying libraries still requires the corresponding
#       Tasmanian_ENABLE_ option
#
########################################################################

########################################################################
# The build process follows these steps:
# 1. Define all options and assign default values
# 2. Perform sanity check
#    - compliance with XSDK build policies
#    - resolve conflicts, e.g., Python requires shared libs
#    - call find_package() for options that need packages
#    - fail or fallback if find_package() fails for some option
# 3. Enable testing and set the master config file
# 4. Add all subdirectories in the correct order
#    - Interfaces require SparseGirds, but do not depend on each other
#    - DREAM requires SparseGrids
#    - SparseGrids does not require other sub-dirs
# 5. Add the "make test_install" command
# 6. Install top level files
#    - cmake export file
#    - master config file
#    - CMakeLists.txt for the examples (include OpenMP hack if needed)
#    - shell (bash) source file that defines paths
# 7. Print final messages
########################################################################

option(USE_XSDK_DEFAULTS            "Enable xSDK compatibility"              OFF)
option(Tasmanian_ENABLE_RECOMMENDED "Enable (if found) OpenMP, BLAS, and Python, also set optimization flags" OFF)
option(Tasmanian_ENABLE_OPENMP      "Enable OpenMP support for Tasmanian (recommended)"    OFF)
option(Tasmanian_ENABLE_BLAS        "Enable CPU Blas support for Tasmanian (recommended)"  OFF)
option(Tasmanian_ENABLE_PYTHON      "Enable Python interface for Tasmanian (recommended)"  OFF)
option(Tasmanian_ENABLE_CUDA        "Enable Nvidia CUDA kernels and libraries within Tasmanian (stable)"    OFF)
option(Tasmanian_ENABLE_MAGMA       "Enable acceleration using UTK Magma library within Tasmanian (stable)" OFF)
option(Tasmanian_ENABLE_MPI         "Enable MPI support for Tasmanian (experimental)"       OFF)
option(Tasmanian_ENABLE_FORTRAN     "Enable Fortran interface for Tasmanian (experimental)" OFF)

# treat those similar to options, but give values other than ON/OFF
set(Tasmanian_MATLAB_WORK_FOLDER     ""  CACHE PATH   "Enable MATLAB interface and use the path for the MATLAB work folder (stable)")

set(Tasmanian_TESTS_GPU_ID          "-1" CACHE STRING "(integer) specify GPU ID for testing, -1 means running tests on all visible devices")
set(Tasmanian_TESTS_OMP_NUM_THREADS "-1" CACHE STRING "(integer) specify OMP_NUM_THREADS for the tests, if less than 0, this option will be ignored")


########################################################################
# Sanity check, xSDK compatibility, and find_package() calls are all
# done in one place for consistency
# after this call, every enabled option comes with proper found set of
# find_package() provided libraries, includes, etc.
########################################################################
include("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/sanity_check_and_xsdk.cmake")


########################################################################
# Core project configuration
########################################################################
enable_testing() # call in one place, each subfolder gives additional tests

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianConfig.in.hpp"
               "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.hpp")


########################################################################
# Setup targets
# CXX subdirs have to come first, SparseGrids must precede DREAM
# each sub-directory creates a set of targets and links those to
# existing targets (hence the order is important)
# After CXX, the interfaces Fortran, Python, MATLAB can come in any order
########################################################################
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/SparseGrids" "${CMAKE_CURRENT_BINARY_DIR}/SparseGrids")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/DREAM"       "${CMAKE_CURRENT_BINARY_DIR}/DREAM")

if (Tasmanian_ENABLE_FORTRAN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfaceFortran"  "${CMAKE_CURRENT_BINARY_DIR}/Fortran")
endif()

if (Tasmanian_ENABLE_PYTHON)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfacePython"  "${CMAKE_CURRENT_BINARY_DIR}/Python")
endif()

if (NOT "${Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/InterfaceMATLAB"  "${CMAKE_CURRENT_BINARY_DIR}/MATLAB")
endif()


########################################################################
# Testing: post install, make test_install
# checks if the executables run and if the examples compile and run
########################################################################
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Testing/test_post_install.in.sh" "${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh")
if (Tasmanian_ENABLE_FORTRAN)
    add_custom_target(test_install COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh"
                                           "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
else()
    add_custom_target(test_install COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh"
                                           "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
endif()


########################################################################
# Install exports, config header, examples cmake
########################################################################
# if Tasmanian is included as a sub-dir, then CMAKE_PROJECT_NAME is the name of the parent project, PROJECT_NAME is Tasmanian
# ${CMAKE_PROJECT_NAME}-exports is generated by the Tasmanian sub-dirs, the parent can add to the -exports and install in separate locaiton
# see the comments below about exports and package config
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    install(EXPORT ${CMAKE_PROJECT_NAME}-exports DESTINATION "lib/${CMAKE_PROJECT_NAME}" FILE "${CMAKE_PROJECT_NAME}.cmake")
endif()

# Tasmanian package config files: TasmanianConfig gives the targets, TasmanianConfigVersion gives the version
# TasmanianConfig also requires "lib/${CMAKE_PROJECT_NAME}/${CMAKE_PROJECT_NAME}-<config>.cmake" which is automatically generated
#
# NOTE: if Tasmanian is included in another project with "add_subdirectory()"
# then the files below assume that the parent project calls
# install(EXPORT ${CMAKE_PROJECT_NAME}-exports DESTINATION "lib/${CMAKE_PROJECT_NAME}" FILE "${CMAKE_PROJECT_NAME}.cmake")
# if Tasmanian makes this call from the sub-dir, then "${CMAKE_PROJECT_NAME}-<config>.cmake" is not generated
# with the install command above, find_package(Tasmanian PATHS ${CMAKE_INSTALL_PREFIX}/lib/) can still be used after the parent is installed
# without the command, the CMakeLists.txt file for the Tasmanian examples will be incorrect and Tasmanian "test_install" will fail
#
# The current setup allows the parent project to append more targets to "${CMAKE_PROJECT_NAME}-exports" and install in a different spot
# The parent project can also completely ignore any Tasmanian exports and config files and make thier own files
# END NOTE
configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianConfig.in.cmake"
                              "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.cmake"
                              INSTALL_DESTINATION "lib/Tasmanian/")
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfigVersion.cmake"
                                 COMPATIBILITY AnyNewerVersion)
# not sure why it is necessary to explicitly install TasmanianConfig.cmake, INSTALL_DESTINATION above doesn't work
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.cmake"
        DESTINATION "lib/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfigVersion.cmake"
        DESTINATION "lib/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

# configured header to be included with the headers from the sub-dirs
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.hpp"
        DESTINATION include
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

# cmake file for the examples, to be used post-install
# see the OpenMP remark in the SparseGirds sub-directory
# unset(OpenMP_CXX_LIBRARIES) # use this to test on newer cmake
set(Tasmanian_openmp_hack "# ALL GOOD, the commands above are sufficient to link any external project to Tasmanian")
if (Tasmanian_ENABLE_OPENMP AND (NOT (DEFINED OpenMP_CXX_LIBRARIES)))
    set(Tasmanian_openmp_hack "set(CMAKE_CXX_FLAGS \"\$\{CMAKE_CXX_FLAGS\} ${OpenMP_CXX_FLAGS}\") # see the long comment below")
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeLists.examples.txt" "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt"
        DESTINATION "share/Tasmanian/examples/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

# configure environment shell script that can be sourced to set path and lib path
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianENVsetup.in.sh" "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh"
        DESTINATION "share/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)


########################################################################
# Setup Message
########################################################################
if (NOT USE_XSDK_DEFAULTS)
    message("")
    message("TASMANIAN: summary of all compile options")
    message(" -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
    message(" -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}")
    message(" -D CMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}")
    if (Tasmanian_ENABLE_CUDA)
        message(" -D CUDA_NVCC_FLAGS:STRING=${CUDA_NVCC_FLAGS}")
    endif()
    if (DEFINED BUILD_SHARED_LIBS)
        message(" -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}")
    else()
        message(" -D BUILD_SHARED_LIBS=Undefined")
    endif()
    foreach(Tasmanian_option Tasmanian_ENABLE_OPENMP  Tasmanian_ENABLE_BLAS
                             Tasmanian_ENABLE_MPI     Tasmanian_ENABLE_PYTHON
                             Tasmanian_ENABLE_CUDA    Tasmanian_ENABLE_MAGMA
                             Tasmanian_ENABLE_FORTRAN)
        if (${Tasmanian_option})
            message(" -D ${Tasmanian_option}:BOOL=ON")
        else()
            message(" -D ${Tasmanian_option}:BOOL=OFF")
        endif()
    endforeach()
    if (Tasmanian_MAGMA AND Tasmanian_MAGMA_ROOT)
        message(" -D Tasmanian_MAGMA_ROOT:PATH=${Tasmanian_MAGMA_ROOT}")
    endif()
    if (NOT "${Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
        message(" -D Tasmanian_MATLAB_WORK_FOLDER:PATH=${Tasmanian_MATLAB_WORK_FOLDER}")
        message(" pre-install MATLAB folder: addpath('${CMAKE_CURRENT_BINARY_DIR}/MATLAB/matlab/')")
    endif()
    message("")
endif()


########################################################################
# Print final message (a bit of a hack)
# The message is written in the CMakeLists.txt, as the subdir is added last
# this ensures that the message will appear last in the install process
# do not print if USE_XSDK or Tasmanian has been imported with addsubdir
########################################################################
if (NOT USE_XSDK_DEFAULTS AND (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME}))
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/")
endif()
