cmake_minimum_required(VERSION 3.5)

########################################################################
# Fortran librareis and command line tools
########################################################################

macro(Tasmanian_macro_add_libdfortran Tasmanian_shared_or_static)

if (${Tasmanian_shared_or_static} STREQUAL SHARED)
    set(Tasmanian_libfor_target_name "Tasmanian_libfortran90_shared")
else()
    set(Tasmanian_libfor_target_name "Tasmanian_libfortran90_static")
endif()

add_library(${Tasmanian_libfor_target_name} ${Tasmanian_shared_or_static} ${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.f90
                                                                          ${CMAKE_CURRENT_SOURCE_DIR}/tsgC2FortranBridge.f90
                                                                          ${CMAKE_CURRENT_SOURCE_DIR}/tsgC2Fortran.cpp)

target_include_directories(${Tasmanian_libfor_target_name} PUBLIC $<INSTALL_INTERFACE:include>)

set_target_properties(${Tasmanian_libfor_target_name} PROPERTIES OUTPUT_NAME "tasmanianfortran90")

target_include_directories(${Tasmanian_libfor_target_name} PUBLIC $<INSTALL_INTERFACE:include>)
target_include_directories(${Tasmanian_libfor_target_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>)

# setup the name, in MacOSX and Linux static libraries use .a extension and shared libs use .dylib and .so
# in Windows, shared .dll libs need another file for linking which shares .lib extension with the static libs
# thus on Windows, we specifically add _static to the name of the lib
if ((${CMAKE_SYSTEM_NAME} STREQUAL "Windows") AND (${Tasmanian_shared_or_static} STREQUAL STATIC))
    set_target_properties(${Tasmanian_libtsg_target_name} PROPERTIES OUTPUT_NAME "tasmanianfortran90_static")
endif()

endmacro()

########################################################################
# add the fortran tester and examples executables
########################################################################
add_executable(Tasmanian_fortester "${CMAKE_CURRENT_SOURCE_DIR}/../Testing/fortester.f90")
add_executable(Tasmanian_example_sparse_grids_f90 "${CMAKE_CURRENT_SOURCE_DIR}/Examples/example_sparse_grids.f90")

set_target_properties(Tasmanian_fortester                PROPERTIES OUTPUT_NAME "fortester")
set_target_properties(Tasmanian_example_sparse_grids_f90 PROPERTIES OUTPUT_NAME "example_sparse_grids_f90")

########################################################################
# add the libraries
########################################################################
if (Tasmanian_SHARED_LIBRARY)
    Tasmanian_macro_add_libdfortran(SHARED)

    target_link_libraries(Tasmanian_libfortran90_shared Tasmanian_libsparsegrid_shared)

    list(APPEND Tasmanian_install_fortran_targets Tasmanian_libfortran90_shared)
endif()

if (Tasmanian_STATIC_LIBRARY)
    Tasmanian_macro_add_libdfortran(STATIC)

    target_link_libraries(Tasmanian_libfortran90_static      Tasmanian_libsparsegrid_static)
    target_link_libraries(Tasmanian_fortester                Tasmanian_libfortran90_static)
    target_link_libraries(Tasmanian_example_sparse_grids_f90 Tasmanian_libfortran90_static)

    list(APPEND Tasmanian_install_fortran_targets Tasmanian_libfortran90_static)
else()
    target_link_libraries(Tasmanian_fortester                Tasmanian_libfortran90_shared)
    target_link_libraries(Tasmanian_example_sparse_grids_f90 Tasmanian_libfortran90_shared)
endif()


########################################################################
# Testing
########################################################################
add_test(Fortran90 ${CMAKE_CURRENT_BINARY_DIR}/fortester)
set(Tasmanian_tests_list_fortran Fortran90)


########################################################################
# Installation
########################################################################
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tasmaniansg.mod"
        DESTINATION include
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/Examples/example_sparse_grids.f90"
        DESTINATION "examples"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)

install(TARGETS ${Tasmanian_install_fortran_targets}
        EXPORT Tasmanian_exports
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")
