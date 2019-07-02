cmake_minimum_required(VERSION 3.5)

@PACKAGE_INIT@

# https://cmake.org/cmake/help/v3.5/module/CMakePackageConfigHelpers.html#module:CMakePackageConfigHelpers
# seems to indicate that I need to pre-pend PACKAGE_ to the variable names after the @PACKAGE_INIT@ macro
# but this doesn't seem to work, not sure if this is a "relocatable package" (low concern)
include("@CMAKE_INSTALL_PREFIX@/lib/@CMAKE_PROJECT_NAME@/@CMAKE_PROJECT_NAME@.cmake")

check_required_components(Tasmanian)

add_library(Tasmanian_libsparsegrid INTERFACE)
add_library(Tasmanian_libdream INTERFACE)

add_library(Tasmanian::Tasmanian INTERFACE IMPORTED GLOBAL)
set_target_properties(Tasmanian::Tasmanian PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_master)

if (TARGET Tasmanian_shared)
    add_library(Tasmanian::Tasmanian_shared INTERFACE IMPORTED GLOBAL)
    set_target_properties(Tasmanian::Tasmanian_shared PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_shared)
endif()

# define _static/_shared independent libraries, default to static if both types are present
if (TARGET Tasmanian_static)
    target_link_libraries(Tasmanian_libsparsegrid INTERFACE Tasmanian_libsparsegrid_static)
    target_link_libraries(Tasmanian_libdream INTERFACE Tasmanian_libdream_static)

    add_library(Tasmanian::Tasmanian_static INTERFACE IMPORTED GLOBAL)
    set_target_properties(Tasmanian::Tasmanian_static PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_static)

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
endif()

# export the python path so other projects can configure python scripts
if (@Tasmanian_ENABLE_PYTHON@)
    set(Tasmanian_PYTHONPATH "@Tasmanian_PYTHONPATH@")
endif()
