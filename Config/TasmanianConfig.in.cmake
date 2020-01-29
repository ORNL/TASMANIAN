cmake_minimum_required(VERSION 3.10)

@PACKAGE_INIT@

# https://cmake.org/cmake/help/v3.5/module/CMakePackageConfigHelpers.html#module:CMakePackageConfigHelpers
# seems to indicate that I need to pre-pend PACKAGE_ to the variable names after the @PACKAGE_INIT@ macro
# but this doesn't seem to work, not sure if this is a "relocatable package" (low concern)
include("@Tasmanian_final_install_path@/lib/@CMAKE_PROJECT_NAME@/@CMAKE_PROJECT_NAME@.cmake")

add_executable(Tasmanian::tasgrid IMPORTED)
set_property(TARGET Tasmanian::tasgrid PROPERTY IMPORTED_LOCATION "@Tasmanian_final_install_path@/bin/tasgrid${CMAKE_EXECUTABLE_SUFFIX_CXX}")

add_library(Tasmanian::Tasmanian INTERFACE IMPORTED GLOBAL)

add_library(Tasmanian_libsparsegrid INTERFACE) # for backwards compatibility
add_library(Tasmanian_libdream INTERFACE)

if (TARGET Tasmanian_shared)
    add_library(Tasmanian::shared INTERFACE IMPORTED GLOBAL)
    set_target_properties(Tasmanian::shared PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_shared)
    set(Tasmanian_SHARED_FOUND "ON")
endif()

# define _static/_shared independent libraries, default to static if both types are present
if (TARGET Tasmanian_static)
    target_link_libraries(Tasmanian_libsparsegrid INTERFACE Tasmanian_libsparsegrid_static)
    target_link_libraries(Tasmanian_libdream INTERFACE Tasmanian_libdream_static)

    add_library(Tasmanian::static INTERFACE IMPORTED GLOBAL)
    set_target_properties(Tasmanian::static PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_static)
    set(Tasmanian_STATIC_FOUND "ON")

    if (@Tasmanian_ENABLE_CUDA@)
    # Since Tasmanian does not transitively include <cuda.h> and since all CUDA calls are wrapped in CXX API,
    # projects do not require cuda language to link to Tasmanian; however, CMake adds the extraneous dependence.
    # If Tasmanian was build with CUDA and if the user has not explicitly enabled the CUDA language,
    # then overwrite the CMake generated extraneous CUDA requirements and link with the CXX compiler only.
    # This hack is not necessary when building shared libraries.
        if (NOT CMAKE_CUDA_COMPILER)
            set_target_properties(Tasmanian_libsparsegrid_static PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX")
            set_target_properties(Tasmanian_libsparsegrid_static PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX")
            get_target_property(TasLibs Tasmanian_libsparsegrid_static INTERFACE_LINK_LIBRARIES)
            if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
                set_target_properties(Tasmanian_libsparsegrid_static PROPERTIES INTERFACE_LINK_LIBRARIES "${TasLibs};@Tasmanian_cudaruntime@")
            else()
                set_target_properties(Tasmanian_libsparsegrid_static PROPERTIES INTERFACE_LINK_LIBRARIES "${TasLibs};@Tasmanian_cudaruntime@;dl")
            endif()
        endif()
    endif()
else()
    target_link_libraries(Tasmanian_libsparsegrid INTERFACE Tasmanian_libsparsegrid_shared)
    target_link_libraries(Tasmanian_libdream INTERFACE Tasmanian_libdream_shared)
endif()

if (@Tasmanian_ENABLE_FORTRAN@)
    add_library(Tasmanian_libfortran90 INTERFACE)
    if (TARGET Tasmanian_libfortran90_static)
        target_link_libraries(Tasmanian_libfortran90 INTERFACE Tasmanian_libfortran90_static)
    else()
        target_link_libraries(Tasmanian_libfortran90 INTERFACE Tasmanian_libfortran90_shared)
    endif()
    set(Tasmanian_FORTRAN_FOUND "ON")

    add_library(Tasmanian::Fortran INTERFACE IMPORTED GLOBAL)
    set_target_properties(Tasmanian::Fortran PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian::Tasmanian)
endif()

# export the python path so other projects can configure python scripts
if (@Tasmanian_ENABLE_PYTHON@)
    set_and_check(Tasmanian_PYTHONPATH "@Tasmanian_final_install_path@/share/Tasmanian/python/")
    set(Tasmanian_PYTHON_FOUND "ON")
endif()

# export the MATLAB paths so other projects can write files directly to MATLAB
if (NOT "@Tasmanian_MATLAB_WORK_FOLDER@" STREQUAL "")
    set_and_check(Tasmanian_MATLAB_WORK_FOLDER "@Tasmanian_MATLAB_WORK_FOLDER@")
    set_and_check(Tasmanian_MATLABPATH "@Tasmanian_final_install_path@/share/Tasmanian/matlab/")
    set(Tasmanian_MATLAB_FOUND "ON")
endif()

set(Tasmanian_OPENMP_FOUND "@Tasmanian_ENABLE_OPENMP@")
set(Tasmanian_BLAS_FOUND   "@Tasmanian_ENABLE_BLAS@")
set(Tasmanian_MPI_FOUND    "@Tasmanian_ENABLE_MPI@")
set(Tasmanian_CUDA_FOUND   "@Tasmanian_ENABLE_CUDA@")
set(Tasmanian_MAGMA_FOUND  "@Tasmanian_ENABLE_MAGMA@")

# write component info
foreach(_comp ${Tasmanian_FIND_COMPONENTS})
    if (Tasmanian_${_comp}_FOUND)
        message(STATUS "Tasmanian component ${_comp}: found")
    else()
        if (Tasmanian_FIND_REQUIRED_${_comp})
            message(WARNING "Tasmanian required component ${_comp}: missing (error)")
        else()
            message(STATUS "Tasmanian optional component ${_comp}: missing")
        endif()
    endif()
endforeach()
unset(_comp)

check_required_components(Tasmanian)

# if find_package(Tasmanian REQUIRED SHARED) is called without STATIC then default to shared libraries
if ((SHARED IN_LIST Tasmanian_FIND_COMPONENTS) AND (NOT STATIC IN_LIST Tasmanian_FIND_COMPONENTS) AND (TARGET Tasmanian_shared))
    set_target_properties(Tasmanian::Tasmanian PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_shared)
else() # otherwise use the default (static if existing, else shared)
    set_target_properties(Tasmanian::Tasmanian PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_master)
endif()
