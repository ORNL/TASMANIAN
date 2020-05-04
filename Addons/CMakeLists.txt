########################################################################
# Addon templates and related tests and interface options
########################################################################

macro(Tasmanian_addon_sources)
    cmake_parse_arguments(_tsg_sources "" "PATH" "FILES" ${ARGN} )
    target_include_directories(Tasmanian_addons INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_tsg_sources_PATH}>)
    foreach(_tsg_file ${_tsg_sources_FILES})
        target_sources(Tasmanian_addons INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_tsg_sources_PATH}/${_tsg_file}>
                                                  $<INSTALL_INTERFACE:${Tasmanian_final_install_path}/include/${_tsg_file}>)
    endforeach()
    unset(_tsg_file)
    unset(_tsg_sources_FILES)
    unset(_tsg_sources_PATH)
endmacro(Tasmanian_addon_sources)

# Tasmanian_addons target is not sufficient to use as a stand-alone as it requires dream_static or dream_shared
add_library(Tasmanian_addons INTERFACE)

Tasmanian_addon_sources(PATH .
                        FILES tsgAddonsCommon.hpp
                              tsgMPIScatterGrid.hpp
                              tsgMPIScatterDream.hpp
                              tsgMPISampleDream.hpp
                              tsgLoadNeededValues.hpp
                              tsgLoadUnstructuredPoints.hpp
                              tsgCandidateManager.hpp
                              tsgConstructSurrogate.hpp
                              tsgMPIConstructGrid.hpp
                              TasmanianAddons.hpp)

if (Tasmanian_ENABLE_PYTHON) # the C addons are used by Python
    add_library(Tasmanian_caddons SHARED tsgCLoadNeededValues.cpp tsgCConstructSurrogate.cpp tsgCLoadUnstructuredPoints.cpp)

    set_target_properties(Tasmanian_caddons PROPERTIES OUTPUT_NAME "tasmaniancaddons"
                                                       INSTALL_RPATH "${Tasmanian_final_install_path}/lib"
                                                       LINK_LIBRARIES "Tasmanian_addons;Tasmanian_libdream_shared"
                                                       SOVERSION ${Tasmanian_VERSION_MAJOR}
                                                       VERSION   ${PROJECT_VERSION})

    install(TARGETS Tasmanian_caddons
            EXPORT  "${Tasmanian_export_name}"
            RUNTIME DESTINATION "bin"
            LIBRARY DESTINATION "lib"
            ARCHIVE DESTINATION "lib")
endif()

if (NOT Tasmanian_ENABLE_OPENMP) # using threading without OpenMP
    target_link_libraries(Tasmanian_addons INTERFACE ${CMAKE_THREAD_LIBS_INIT})
endif()

# The Tasmanian MPI capabilities are templated into the Addons
if (Tasmanian_ENABLE_MPI)
    target_link_libraries(Tasmanian_addons INTERFACE MPI::MPI_CXX)

    add_executable(Tasmanian_mpitester testMPI.cpp testMPI.hpp testMPIDream.hpp)
    add_test(MPISparseGridsIO ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/mpitester ${MPIEXEC_POSTFLAGS})

    list(APPEND Tasmanian_addons_executables "mpitester")
    list(APPEND Tasmanian_addons_tests "MPISparseGridsIO")
endif()

# test for non-MPI capabilities
add_executable(Tasmanian_addontester testAddons.cpp testConstructSurrogate.hpp testLoadUnstructured.hpp)
add_test(AddonTests addontester -gpuid ${Tasmanian_TESTS_GPU_ID})

list(APPEND Tasmanian_addons_executables "addontester")
list(APPEND Tasmanian_addons_tests "AddonTests")

# set properties for all executable targets
foreach(_tsg_exec ${Tasmanian_addons_executables})
    set_target_properties(Tasmanian_${_tsg_exec} PROPERTIES OUTPUT_NAME "${_tsg_exec}" CXX_EXTENSIONS OFF)
    target_link_libraries(Tasmanian_${_tsg_exec} Tasmanian_addons Tasmanian_libdream)
endforeach()
unset(_tsg_exec)

# cheap models have huge random variation in time-per-sample which can cause tests to fail
# make sure other tests are not holding resources and making things even worse
set_tests_properties(${Tasmanian_addons_tests} PROPERTIES RUN_SERIAL ON)
if (Tasmanian_TESTS_OMP_NUM_THREADS GREATER 0)
    set_tests_properties(${Tasmanian_addons_tests}
        PROPERTIES
        PROCESSORS "${Tasmanian_TESTS_OMP_NUM_THREADS}"
        ENVIRONMENT "OMP_NUM_THREADS=${Tasmanian_TESTS_OMP_NUM_THREADS}")
endif()

# Windows specific support (DLL export/import directives and names)
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    if (TARGET Tasmanian_caddons)
        set_target_properties(Tasmanian_caddons PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    endif()
    if (NOT "static" IN_LIST Tasmanian_libs_type)
        set(Tasmanian_MSVC_PATH_STRING "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Release;${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Debug")
        set(Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING};${CMAKE_CURRENT_BINARY_DIR}/../DREAM/Release;${CMAKE_CURRENT_BINARY_DIR}/../DREAM/Debug;$ENV{PATH}")
        string(REPLACE ";" "\\;" Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING}")
        set_tests_properties(${Tasmanian_addons_tests} PROPERTIES ENVIRONMENT "PATH=${Tasmanian_MSVC_PATH_STRING}")
    endif()
endif()

# install directives
target_include_directories(Tasmanian_addons INTERFACE $<INSTALL_INTERFACE:include>)
get_target_property(_tsg_sources Tasmanian_addons INTERFACE_SOURCES)
install(FILES ${_tsg_sources} DESTINATION include)
install(TARGETS Tasmanian_addons EXPORT "${Tasmanian_export_name}" PUBLIC_HEADER DESTINATION include)

# cleanup
unset(_tsg_sources)
unset(Tasmanian_addons_executables)
unset(Tasmanian_addons_tests)
