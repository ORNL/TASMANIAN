cmake_minimum_required(VERSION 3.10)

cmake_policy(VERSION 3.10)
project(Tasmanian VERSION 6.1.0 LANGUAGES CXX)
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
#    -D Tasmanian_ENABLE_FORTRAN:BOOL=OFF     (mostly stable)
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
#     -D CMAKE_CUDA_COMPILER:PATH
#     -D Tasmanian_MAGMA_ROOT_DIR:PATH (or just MAGMA_ROOT_DIR)
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
option(Tasmanian_ENABLE_FORTRAN     "Enable Fortran interface for Tasmanian (mostly stable)" OFF)
option(Tasmanian_ENABLE_MPI         "Enable MPI support for Tasmanian (experimental)"       OFF)
option(Tasmanian_ENABLE_DOXYGEN     "Enable Dogygen documentation for Tasmanian (not fully complete)" OFF)

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
########################################################################
include("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/sanity_check_and_xsdk.cmake")


########################################################################
# Core project configuration, passes all options to CXX
########################################################################
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianConfig.in.hpp"
               "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.hpp")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.hpp"
        DESTINATION include
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)


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

if (Tasmanian_ENABLE_DOXYGEN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Doxygen"  "${CMAKE_CURRENT_BINARY_DIR}/Doxygen")
endif()


########################################################################
# Skip test post-install, exports and package config,
# in the case when Tasmanian is merged into an external master project.
# The master should be responsible for such testing and exporting
########################################################################
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
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
# Install exports, package-config files, and examples cmake file
########################################################################
    install(EXPORT "${Tasmanian_export_name}" DESTINATION "lib/${CMAKE_PROJECT_NAME}" FILE "${CMAKE_PROJECT_NAME}.cmake")

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

    # cmake file for the examples, to be used post-install
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeLists.examples.txt" "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt" @ONLY)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt"
            DESTINATION "share/Tasmanian/examples/"
            PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
endif()


# configure environment shell script that can be sourced to set path and lib path
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianENVsetup.in.sh" "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh"
        DESTINATION "share/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)


########################################################################
# Setup Message
########################################################################
if (NOT USE_XSDK_DEFAULTS)
    message(STATUS "")
    message(STATUS "Tasmanian ${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR} (development version): summary of build options")
    message(STATUS " -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
    message(STATUS " -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}")
    message(STATUS " -D CMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        message(STATUS " -D CMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}") # useful for Windows debugging
        message(STATUS " -D CMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}")
    endif()
    if (Tasmanian_ENABLE_CUDA)
        message(STATUS " -D CUDA_NVCC_FLAGS:STRING=${CUDA_NVCC_FLAGS}")
    endif()
    if (DEFINED BUILD_SHARED_LIBS)
        message(STATUS " -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}")
    else()
        message(STATUS " -D BUILD_SHARED_LIBS=Undefined")
    endif()
    foreach(Tasmanian_option Tasmanian_ENABLE_OPENMP  Tasmanian_ENABLE_BLAS
                             Tasmanian_ENABLE_MPI     Tasmanian_ENABLE_PYTHON
                             Tasmanian_ENABLE_CUDA    Tasmanian_ENABLE_MAGMA
                             Tasmanian_ENABLE_FORTRAN Tasmanian_ENABLE_DOXYGEN)
        if (${Tasmanian_option})
            message(STATUS " -D ${Tasmanian_option}:BOOL=ON")
        else()
            message(STATUS " -D ${Tasmanian_option}:BOOL=OFF")
        endif()
    endforeach()
    if (Tasmanian_MAGMA AND Tasmanian_MAGMA_ROOT)
        message(STATUS " -D Tasmanian_MAGMA_ROOT:PATH=${Tasmanian_MAGMA_ROOT}")
    endif()
    if (NOT "${Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
        message(STATUS " -D Tasmanian_MATLAB_WORK_FOLDER:PATH=${Tasmanian_MATLAB_WORK_FOLDER}")
        message(STATUS " pre-install MATLAB folder: addpath('${CMAKE_CURRENT_BINARY_DIR}/MATLAB/matlab/')")
    endif()
    message(STATUS "")
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
