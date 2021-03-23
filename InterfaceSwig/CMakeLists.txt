if (Tasmanian_ENABLE_SWIG)
  if (CMAKE_VERSION VERSION_LESS 3.18)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/backport-cmake-318")
  endif()

  find_package(SWIG COMPONENTS fortran REQUIRED)

  # SWIG is requested and available; make sure it's the Fortran fork.
  cmake_policy(SET CMP0078 "NEW")
  cmake_policy(SET CMP0086 "NEW")
  include(UseSWIG)
endif()

set(Tasmanian_GENERATE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/generated")

# usage: Tasmanian_set_fortran_props(TARGET Tasmanian_libdream OUTPUT_NAME "tasmanianfortran")
# sets the common properties of Tasmanian Fortran targets
macro(Tasmanian_set_fortran_props)
    cmake_parse_arguments(Tasmanian_fprop "" "TARGET;OUTPUT_NAME" "" ${ARGN})
    target_link_libraries(${Tasmanian_fprop_TARGET} Tasmanian_addons)
    set_target_properties(${Tasmanian_fprop_TARGET} PROPERTIES
                                                    OUTPUT_NAME ${Tasmanian_fprop_OUTPUT_NAME}
                                                    INSTALL_RPATH "${Tasmanian_rpath}"
                                                    SOVERSION ${Tasmanian_VERSION_MAJOR}
                                                    VERSION   ${PROJECT_VERSION})
    target_include_directories(${Tasmanian_fprop_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>)

    if (Tasmanian_ENABLE_MPI)
        # see the comments in the Addons/CMakeLists.txt
        target_link_libraries(${Tasmanian_fprop_TARGET} MPI::MPI_Fortran)
    endif()

    unset(Tasmanian_fprop_TARGET)
    unset(Tasmanian_fprop_OUTPUT_NAME)
endmacro()

function(tasmanian_add_swig_module tsg_name)
#
# Generates the wrappers using SWIG
#
    # SWIG is available; actually generate the library dynamically.
    set(src_file "${CMAKE_CURRENT_SOURCE_DIR}/${tsg_name}.i")
    # We're using C++
    set_property(SOURCE "${src_file}" PROPERTY CPLUSPLUS ON)
    # We need to include the source directory
    set_property(SOURCE "${src_file}" PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON)

    # Create the library (testing the generated wrappers)
    swig_add_library(${tsg_name}
      LANGUAGE Fortran
      TYPE USE_BUILD_SHARED_LIBS
      OUTPUT_DIR "${Tasmanian_GENERATE_DIR}"
      SOURCES "${src_file}" ${ARGN}
    )

    Tasmanian_set_fortran_props(TARGET ${tsg_name} OUTPUT_NAME ${tsg_name})

endfunction()

macro(Tasmanian_macro_add_libdfortran03)
#
# Build the Fortran module from an existing generated source file
#

    file(READ "${Tasmanian_GENERATE_DIR}/tasmanian_swig.f90" Tasmanian_generated_f90)
    string(REPLACE "module tasmanian_swig" "module tasmanian" Tasmanian_generated_f90 ${Tasmanian_generated_f90})
    string(REPLACE " type, public :: TasmanianSparseGrid
  type(SwigClassWrapper), public :: swigdata
 contains"
    "type, public :: TasmanianSparseGrid
  type(SwigClassWrapper), public :: swigdata
 contains
    procedure :: returnPoints => tsgGetPoints
    procedure :: returnLoadedPoints => tsgGetLoadedPoints
    procedure :: returnNeededPoints => tsgGetNeededPoints
    procedure :: returnQuadratureWeights => tsgGetQuadratureWeights
    "
    Tasmanian_generated_f90 ${Tasmanian_generated_f90})

    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/regenerated/tasmanian.f90" "${Tasmanian_generated_f90}")

    add_library(Tasmanian_libfortran03  "${CMAKE_CURRENT_BINARY_DIR}/regenerated/tasmanian.f90"
                                        "${Tasmanian_GENERATE_DIR}/tasmanian_swigFORTRAN_wrap.cxx")

    Tasmanian_set_fortran_props(TARGET Tasmanian_libfortran03 OUTPUT_NAME "tasmanianfortran")

    install(TARGETS Tasmanian_libfortran03
            EXPORT "${Tasmanian_export_name}"
            RUNTIME DESTINATION "bin"
            LIBRARY DESTINATION "lib"
            ARCHIVE DESTINATION "lib")

    if (Tasmanian_ENABLE_MPI)
        file(READ "${Tasmanian_GENERATE_DIR}/tasmanian_mpi_swig.f90" Tasmanian_generated_f90)
        string(REPLACE "module tasmanian_mpi_swig" "module tasmanian_mpi" Tasmanian_generated_f90 ${Tasmanian_generated_f90})
        string(REPLACE "use tasmanian_swig" "use tasmanian" Tasmanian_generated_f90 ${Tasmanian_generated_f90})

        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/regenerated/tasmanian_mpi.f90" "${Tasmanian_generated_f90}")

        add_library(Tasmanian_libfortranmpi03  "${CMAKE_CURRENT_BINARY_DIR}/regenerated/tasmanian_mpi.f90"
                                               "${Tasmanian_GENERATE_DIR}/tasmanian_mpi_swigFORTRAN_wrap.cxx")

        Tasmanian_set_fortran_props(TARGET Tasmanian_libfortranmpi03 OUTPUT_NAME "tasmanianfortranmpi")
        target_link_libraries(Tasmanian_libfortranmpi03  Tasmanian_libfortran03)

        install(TARGETS Tasmanian_libfortranmpi03
                EXPORT "${Tasmanian_export_name}"
                RUNTIME DESTINATION "bin"
                LIBRARY DESTINATION "lib"
                ARCHIVE DESTINATION "lib")
        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tasmanian_mpi.mod"
                DESTINATION include
                PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    endif()

endmacro()

########################################################################
# Fortran librareis and command line tools
########################################################################

add_executable(Tasmanian_fortran_swig_tester ${CMAKE_CURRENT_SOURCE_DIR}/FortranTests/fortranswigtester.f90
                                             ${CMAKE_CURRENT_SOURCE_DIR}/FortranTests/test_common.f90
                                             ${CMAKE_CURRENT_SOURCE_DIR}/FortranTests/test_make_grid.f90
                                             ${CMAKE_CURRENT_SOURCE_DIR}/FortranTests/test_domain.f90
                                             ${CMAKE_CURRENT_SOURCE_DIR}/FortranTests/test_evaluate.f90
                                             ${CMAKE_CURRENT_SOURCE_DIR}/FortranTests/test_update_grid.f90
                                             ${CMAKE_CURRENT_SOURCE_DIR}/FortranTests/test_refinement.f90
                                             )
set_target_properties(Tasmanian_fortran_swig_tester PROPERTIES OUTPUT_NAME "forswigtester")
Tasmanian_set_fortran_linker(TARGET Tasmanian_fortran_swig_tester)

if (Tasmanian_ENABLE_SWIG)
    # generate the interface
    tasmanian_add_swig_module(tasmanian_swig)
    if (Tasmanian_ENABLE_MPI)
        tasmanian_add_swig_module(tasmanian_mpi_swig)
        target_link_libraries(tasmanian_mpi_swig  tasmanian_swig)
    endif()
endif()

Tasmanian_macro_add_libdfortran03()

target_link_libraries(Tasmanian_fortran_swig_tester  Tasmanian_libfortran03)

add_test(Fortran2003 forswigtester)
Tasmanian_set_test_properties(TESTS Fortran2003)

if (Tasmanian_ENABLE_MPI)
    add_executable(Tasmanian_fortran_mpiswig_tester ${CMAKE_CURRENT_SOURCE_DIR}/FortranTests/mpitester.f90
                                                    )
    set_target_properties(Tasmanian_fortran_mpiswig_tester PROPERTIES OUTPUT_NAME "mpiforswigtester")
    Tasmanian_set_fortran_linker(TARGET Tasmanian_fortran_mpiswig_tester)

    target_link_libraries(Tasmanian_fortran_mpiswig_tester  Tasmanian_libfortranmpi03)

    add_test(MPIFortran2003 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/mpiforswigtester ${MPIEXEC_POSTFLAGS})
    set_tests_properties(MPIFortran2003 PROPERTIES RUN_SERIAL ON)
    Tasmanian_set_test_properties(TESTS MPIFortran2003)
endif()

########################################################################
# Installation
########################################################################

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tasmanian.mod"
        DESTINATION include
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/FortranExamples/"
        DESTINATION "share/Tasmanian/examples/"
        FILES_MATCHING PATTERN "*.f90"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
