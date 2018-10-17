cmake_minimum_required(VERSION 3.5)

@PACKAGE_INIT@

# https://cmake.org/cmake/help/v3.5/module/CMakePackageConfigHelpers.html#module:CMakePackageConfigHelpers
# seems to indicate that I need to pre-pend PACKAGE_ to the variable names after the @PACKAGE_INIT@ macro
# but this doesn't seem to work, not sure if this is a "relocatable package" (low concern)
include("@CMAKE_INSTALL_PREFIX@/lib/@CMAKE_PROJECT_NAME@/@CMAKE_PROJECT_NAME@.cmake")

check_required_components(Tasmanian)

add_library(Tasmanian_libsparsegrid INTERFACE)
add_library(Tasmanian_libdream INTERFACE)

# define _static/_shared independent libraries, default to static if both types are present
if (TARGET Tasmanian_libsparsegrid_static)
    target_link_libraries(Tasmanian_libsparsegrid INTERFACE Tasmanian_libsparsegrid_static)
    target_link_libraries(Tasmanian_libdream INTERFACE Tasmanian_libdream_static)
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
