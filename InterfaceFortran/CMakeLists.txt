########################################################################
# Fortran librareis and command line tools
########################################################################

add_library(Tasmanian_libfortran90 TasmanianSG.f90
                                   tsgC2FortranBridge.f90
                                   tsgC2Fortran.cpp)

target_include_directories(Tasmanian_libfortran90 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>)
target_link_libraries(Tasmanian_libfortran90  Tasmanian_addons)

set_target_properties(Tasmanian_libfortran90 PROPERTIES OUTPUT_NAME "tasmanianfortran90"
                                             SOVERSION ${Tasmanian_VERSION_MAJOR}
                                             VERSION   ${PROJECT_VERSION})
Tasmanian_rpath_target(TARGET Tasmanian_libfortran90 COMPONENTS SparseGrids DREAM)

install(TARGETS Tasmanian_libfortran90
        EXPORT "${Tasmanian_export_name}"
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")


########################################################################
# add the fortran tester and examples executables
########################################################################
add_executable(Tasmanian_fortester fortester.f90)

set_target_properties(Tasmanian_fortester  PROPERTIES OUTPUT_NAME "fortester" LINKER_LANGUAGE Fortran)
Tasmanian_rpath_target(TARGET Tasmanian_fortester USE_CURRENT COMPONENTS SparseGrids DREAM)

target_link_libraries(Tasmanian_fortester  Tasmanian_libfortran90)


########################################################################
# handle the MPI dependence and MPI tests
########################################################################
if (Tasmanian_ENABLE_MPI)
    target_link_libraries(Tasmanian_libfortran90 MPI::MPI_Fortran)

    add_executable(Tasmanian_mpif90    mpitester.f90)
    set_target_properties(Tasmanian_mpif90     PROPERTIES OUTPUT_NAME "mpitester" LINKER_LANGUAGE Fortran)
    target_link_libraries(Tasmanian_mpif90     Tasmanian_libfortran90)
    add_test(MPIFortranGridIO ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/mpitester ${MPIEXEC_POSTFLAGS})
    set_tests_properties(MPIFortranGridIO PROPERTIES RUN_SERIAL ON)
    Tasmanian_set_test_properties(TESTS MPIFortranGridIO)
endif()


########################################################################
# Testing
########################################################################
add_test(Fortran90 fortester)
Tasmanian_set_test_properties(TESTS Fortran90)


########################################################################
# Installation
########################################################################
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tasmaniansg.mod"
        DESTINATION include
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
