########################################################################
# Various CMake macros to be used by Tasmanian
########################################################################

# CMake native includes used in Tasmanian CMake scripts
include(FindPackageHandleStandardArgs)
include(CMakePackageConfigHelpers)

# usage: Tasmanian_find_libraries(REQUIRED foo1 foo2 OPTIONAL foo3 PREFIX bar LIST saloon NO_DEFAULT_PATH)
# this will search for foo1/2/3 in bar/<lib/lib64/arch> with NO_DEFAULT_PATH (skip if defaults are to be used)
# the found libraries will be added to a list Tasmanian_saloon
# the search results will be also added to cached variables Tasmanian_foo1 and Tasmanian_foo2
# that can be used to call find_package_handle_standard_args(... DEFAULT_MSG Tasmanian_foo1 Tasmanian_foo2)
# the macro will not call find_package_handle_standard_args()
# the optional libraries will not create cached entries and mission optional will not be added to the LIST
macro(Tasmanian_find_libraries)
    cmake_parse_arguments(Tasmanian_findlibs "NO_DEFAULT_PATH" "PREFIX;LIST" "REQUIRED;OPTIONAL" ${ARGN})

    foreach(_tsg_lib ${Tasmanian_findlibs_REQUIRED})
        list(APPEND Tasmanian_findlibs_required "Tasmanian_${_tsg_lib}")
    endforeach()

    set(Tasmanian_findlibs_default "")
    if (${Tasmanian_findlibs_NO_DEFAULT_PATH})
        set(Tasmanian_findlibs_default "NO_DEFAULT_PATH")
    endif()

    foreach(_tsg_lib ${Tasmanian_findlibs_REQUIRED} ${Tasmanian_findlibs_OPTIONAL})

        find_library(Tasmanian_${_tsg_lib} ${_tsg_lib}
                     HINTS "${Tasmanian_findlibs_PREFIX}"
                     HINTS "${Tasmanian_findlibs_PREFIX}/lib/"
                     HINTS "${Tasmanian_findlibs_PREFIX}/lib/intel64"
                     HINTS "${Tasmanian_findlibs_PREFIX}/${CMAKE_LIBRARY_ARCHITECTURE}/lib/"
                     HINTS "${Tasmanian_findlibs_PREFIX}/lib/${CMAKE_LIBRARY_ARCHITECTURE}"
                     HINTS "${Tasmanian_findlibs_PREFIX}/lib64/"
                     HINTS "${Tasmanian_findlibs_PREFIX}/${CMAKE_LIBRARY_ARCHITECTURE}/lib64/"
                     HINTS "${Tasmanian_findlibs_PREFIX}/lib64/${CMAKE_LIBRARY_ARCHITECTURE}"
                     ${Tasmanian_findlibs_default})

        if (CMAKE_FIND_DEBUG_MODE)
            message(STATUS "Tasmanian searching library: ${_tsg_lib} => ${Tasmanian_${_tsg_lib}}")
        endif()

        list(FIND Tasmanian_findlibs_required "Tasmanian_${_tsg_lib}" Tasmanian_findlibs_is_required)
        if (${Tasmanian_findlibs_is_required} EQUAL -1) # not a required library
            if (Tasmanian_${_tsg_lib})
                list(APPEND Tasmanian_${Tasmanian_findlibs_LIST} ${Tasmanian_${_tsg_lib}})
            endif()
            unset(Tasmanian_${_tsg_lib} CACHE)
        else()
            list(APPEND Tasmanian_${Tasmanian_findlibs_LIST} ${Tasmanian_${_tsg_lib}})
        endif()

    endforeach()

    foreach(_tsg_lib default required NO_DEFAULT_PATH PREFIX LIST REQUIRED OPTIONAL) # cleanup
        unset(Tasmanian_${_tsg_lib})
    endforeach()
    unset(_tsg_lib)
endmacro()

# usage: Tasmanian_find_header(FILE foo.h RESULT foo_h ROOT bar HINT saloon NO_DEFAULT_PATH)
# will search for file named foo.h and will add the folder in variable Tasmanian_foo_h
# NO_DEFAULT_PATH defines whether to use NO_DEFAULT_PATH in the find_path() command
# in addition to the default paths, the method will search for folders bar and bar/include
# as well as saloon_dir, saloon_dir/include, saloon_dir/../include and saloon_dir/../../include
# where saloon_dir is the directory of the file saloon
macro(Tasmanian_find_header)
    cmake_parse_arguments(Tasmanian_findh "NO_DEFAULT_PATH" "FILE;ROOT;HINT;RESULT" "" ${ARGN})

    if (Tasmanian_findh_HINT)
        get_filename_component(Tasmanian_findh_dir ${Tasmanian_findh_HINT} DIRECTORY)
    endif()

    set(Tasmanian_findh_default "")
    if (${Tasmanian_findh_NO_DEFAULT_PATH})
        set(Tasmanian_findh_default "NO_DEFAULT_PATH")
    endif()

    find_path(Tasmanian_${Tasmanian_findh_RESULT} "${Tasmanian_findh_FILE}"
              HINTS "${Tasmanian_findh_ROOT}"
              HINTS "${Tasmanian_findh_ROOT}/include"
              HINTS "${Tasmanian_findh_dir}"
              HINTS "${Tasmanian_findh_dir}/include"
              HINTS "${Tasmanian_findh_dir}/../include"
              HINTS "${Tasmanian_findh_dir}/../../include"
              ${Tasmanian_findh_default})

    if (CMAKE_FIND_DEBUG_MODE)
        message(STATUS "Tasmanian searching header: ${Tasmanian_findh_FILE} => ${Tasmanian_${Tasmanian_findh_RESULT}}")
    endif()

    foreach(_tsg_opts dir default NO_DEFAULT_PATH FILE ROOT HINT RESULT) # cleanup
        unset(Tasmanian_${_tsg_opts})
    endforeach()
    unset(_tsg_opts)
endmacro()

# usage: Tasmanian_find_rpath(LIBRARIES ${BLAS_LIBRARIES} LIST rpath)
# will find the rpaths for each library in the BLAS_LIBRARIES list
# and will append the rpath to a list called Tasmanian_rpath
macro(Tasmanian_find_rpath)
    cmake_parse_arguments(Tasmanian_findrpath "" "LIST" "LIBRARIES" ${ARGN})

    foreach(_tsg_lib ${Tasmanian_findrpath_LIBRARIES})
        get_filename_component(_tsg_libpath ${_tsg_lib} DIRECTORY)
        list(APPEND Tasmanian_${Tasmanian_findrpath_LIST} ${_tsg_libpath})
        #message(STATUS "rpath for ${_tsg_lib} => ${_tsg_libpath}")
    endforeach()

    foreach(_tsg_lib LIST LIBRARIES) # cleanup
        unset(Tasmanian_${_tsg_lib})
    endforeach()
    unset(_tsg_lib)
    unset(_tsg_libpath)
endmacro()

# usage: Tasmanian_rpath_target(TARGET Tasmanian_libdream USE_CURRENT COMPONENTS SparseGrid)
# will set the rpath for the given target and will handle the work-around for the missing libomp.so under HIP
# the USE_CURRENT adds the current folder, use for executable targets
# COMPONENTS lists the additional folders, e.g., DREAM needs SparseGrid and Addons needs both DREAM and SparseGrid
macro(Tasmanian_rpath_target)
    cmake_parse_arguments(Tasmanian_rpt "USE_CURRENT" "TARGET" "COMPONENTS" ${ARGN})
    if (Tasmanian_ENABLE_HIP)
        set(_rpt_rpath "${Tasmanian_rpath}")
        if (Tasmanian_rpt_USE_CURRENT)
            list(APPEND _rpt_rpath "${CMAKE_CURRENT_BINARY_DIR}")
        endif()
        foreach(_rpt_comp ${Tasmanian_rpt_COMPONENTS})
            list(APPEND _rpt_rpath "${CMAKE_CURRENT_BINARY_DIR}/../${_rpt_comp}")
        endforeach()
        set_target_properties(${Tasmanian_rpt_TARGET} PROPERTIES
                              BUILD_WITH_INSTALL_RPATH "ON"
                              INSTALL_RPATH "${_rpt_rpath}")
        unset(_rpt_rpath)
        unset(_rpt_comp)
    else()
        set_target_properties(${Tasmanian_rpt_TARGET} PROPERTIES
                              INSTALL_RPATH "${Tasmanian_rpath}")
    endif()
    unset(Tasmanian_rpt_TARGET)
    unset(Tasmanian_rpt_USE_CURRENT)
    unset(Tasmanian_rpt_COMPONENTS)
endmacro()

# usage: Tasmanian_rpath_target(TARGET Tasmanian_libdream)
# must be called for every Fortran executable target to account for the missing main() when using ifort
macro(Tasmanian_set_fortran_linker)
    cmake_parse_arguments(Tasmanian_flink "" "TARGET" "" ${ARGN})
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set_target_properties(${Tasmanian_flink_TARGET} PROPERTIES LINKER_LANGUAGE CXX)
    else()
        set_target_properties(${Tasmanian_flink_TARGET} PROPERTIES LINKER_LANGUAGE Fortran)
    endif()
    unset(Tasmanian_flink_TARGET)
endmacro()

# usage: Tasmanian_compiler_type(COMPILER ${CMAKE_CXX_COMPILER} TYPE "DPC++" RESULT Tasmanian_has_dpcpp)
# checks whether the compiler called with "--version" returns a string that contains the TYPE
# the result is written to the boolean variable
macro(Tasmanian_compiler_type)
    cmake_parse_arguments(Tasmanian_ccheck "" "COMPILER;TYPE;RESULT;SWITCH" "" ${ARGN})
    if (NOT Tasmanian_ccheck_SWITCH)
        set(Tasmanian_ccheck_SWITCH "--version")
    endif()
    execute_process(COMMAND ${Tasmanian_ccheck_COMPILER} ${Tasmanian_ccheck_SWITCH} OUTPUT_VARIABLE ${Tasmanian_ccheck_RESULT})
    string(FIND "${${Tasmanian_ccheck_RESULT}}" "${Tasmanian_ccheck_TYPE}" ${Tasmanian_ccheck_RESULT})
    if (NOT ${${Tasmanian_ccheck_RESULT}} LESS 0)
        set(${Tasmanian_ccheck_RESULT} ON)
    else()
        set(${Tasmanian_ccheck_RESULT} OFF)
    endif()
    unset(Tasmanian_ccheck_RESULT)
    unset(Tasmanian_ccheck_TYPE)
    unset(Tasmanian_ccheck_COMPILER)
endmacro()

# usage: Tasmanian_set_test_properties(TESTS Test1 Test2 Test3)
# sets common properties for the Tasmanian tests, e.g., the OMP_NUM_THREADS variable
macro(Tasmanian_set_test_properties)
    cmake_parse_arguments(Tasmanian_tprops "" "" "TESTS" ${ARGN})

    if (Tasmanian_TESTS_OMP_NUM_THREADS GREATER 0)
        set_tests_properties(${Tasmanian_tprops_TESTS}
                             PROPERTIES
                             PROCESSORS "${Tasmanian_TESTS_OMP_NUM_THREADS}"
                             ENVIRONMENT "OMP_NUM_THREADS=${Tasmanian_TESTS_OMP_NUM_THREADS}")
    endif()

    unset(Tasmanian_tprops_TESTS)
endmacro()
