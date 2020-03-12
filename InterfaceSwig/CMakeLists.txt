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

    target_link_libraries(${tsg_name} Tasmanian_addons Tasmanian_libdream_${Tasmanian_lib_default})
    set_target_properties(${tsg_name} PROPERTIES INSTALL_RPATH "${Tasmanian_final_install_path}/lib"
                                                 SOVERSION ${Tasmanian_VERSION_MAJOR}
                                                 VERSION   ${Tasmanian_version_string})

    target_include_directories(${tsg_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>)

    if (Tasmanian_ENABLE_MPI)
        # see the comments in the Addons/CMakeLists.txt
        target_link_libraries(${tsg_name} MPI::MPI_Fortran)
    endif()

    install(TARGETS ${tsg_name}
            EXPORT "${Tasmanian_export_name}"
            RUNTIME DESTINATION "bin"
            LIBRARY DESTINATION "lib"
            ARCHIVE DESTINATION "lib")
endfunction()

macro(Tasmanian_macro_add_libdfortran03 Tasmanian_shared_or_static)
#
# Build the Fortran module from an existing generated source file
#

    if (${Tasmanian_shared_or_static} STREQUAL SHARED)
        set(Tasmanian_libfor_target_name "Tasmanian_libfortran03_shared")
    else()
        set(Tasmanian_libfor_target_name "Tasmanian_libfortran03_static")
    endif()

    add_library(${Tasmanian_libfor_target_name} ${Tasmanian_shared_or_static}
                                                    "${Tasmanian_GENERATE_DIR}/tasmanian.f90"
                                                    "${Tasmanian_GENERATE_DIR}/tasmanianFORTRAN_wrap.cxx")

    target_link_libraries(${Tasmanian_libfor_target_name} Tasmanian_addons)
    set_target_properties(${Tasmanian_libfor_target_name} PROPERTIES OUTPUT_NAME "tasmanianfortran"
                                                                    INSTALL_RPATH "${Tasmanian_final_install_path}/lib"
                                                                    SOVERSION ${Tasmanian_VERSION_MAJOR}
                                                                    VERSION   ${Tasmanian_version_string}
                                                                    LINKER_LANGUAGE Fortran)

    target_include_directories(${Tasmanian_libfor_target_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>)

    if (Tasmanian_ENABLE_MPI)
        # see the comments in the Addons/CMakeLists.txt
        target_link_libraries(${Tasmanian_libfor_target_name} MPI::MPI_Fortran)
    endif()

    install(TARGETS ${Tasmanian_libfor_target_name}
            EXPORT "${Tasmanian_export_name}"
            RUNTIME DESTINATION "bin"
            LIBRARY DESTINATION "lib"
            ARCHIVE DESTINATION "lib")

    unset(Tasmanian_libfor_target_name)
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
                                             )
set_target_properties(Tasmanian_fortran_swig_tester PROPERTIES OUTPUT_NAME "forswigtester"
                                                               LINKER_LANGUAGE Fortran)

if (Tasmanian_ENABLE_SWIG)
    # generate the interface
    tasmanian_add_swig_module(tasmanian)

    target_link_libraries(Tasmanian_fortran_swig_tester  tasmanian)

    add_library(Tasmanian_libfortran03_${Tasmanian_lib_default} ALIAS tasmanian)
else()
    # build the currently generated interface
    foreach(_tsglibtype ${Tasmanian_libs_type})
        string(TOUPPER ${_tsglibtype} Tasmanian_ltype)
        Tasmanian_macro_add_libdfortran03(${Tasmanian_ltype})

        set_source_files_properties("${Tasmanian_GENERATE_DIR}/tasmanian.f90" PROPERTIES LANGUAGE Fortran)
        target_link_libraries(Tasmanian_libfortran03_${_tsglibtype} Tasmanian_libdream_${_tsglibtype})
    endforeach()
    unset(_tsglibtype)
    unset(Tasmanian_ltype)

    target_link_libraries(Tasmanian_fortran_swig_tester  Tasmanian_libfortran03_${Tasmanian_lib_default})
endif()

add_test(Fortran2003 forswigtester)
if (Tasmanian_TESTS_OMP_NUM_THREADS GREATER 0)
    set_tests_properties(Fortran2003
        PROPERTIES
        PROCESSORS "${Tasmanian_TESTS_OMP_NUM_THREADS}"
        ENVIRONMENT "OMP_NUM_THREADS=${Tasmanian_TESTS_OMP_NUM_THREADS}")
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
