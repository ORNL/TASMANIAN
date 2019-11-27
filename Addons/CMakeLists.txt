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
                              tsgCandidateManager.hpp
                              tsgConstructSurrogate.hpp
                              tsgMPIConstructGrid.hpp
                              TasmanianAddons.hpp)

if (Tasmanian_ENABLE_PYTHON) # the C addons are used by Python
    get_target_property(Tasmanian_addon_sources Tasmanian_addons INTERFACE_SOURCES)
    add_library(Tasmanian_caddons SHARED tsgCLoadNeededValues.cpp tsgCConstructSurrogate.cpp)
    target_link_libraries(Tasmanian_caddons Tasmanian_addons Tasmanian_libdream_shared)
    set_target_properties(Tasmanian_caddons PROPERTIES OUTPUT_NAME "tasmaniancaddons"
                                                       INSTALL_RPATH "${Tasmanian_final_install_path}/lib"
                                                       CXX_EXTENSIONS OFF)

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
    # the find_package(MPI) cmake scripts creates an MPI::MPI_CXX target but that cannot be exported
    # this creates a potential issue on a system with multiple MPI versions and communicating the version used on build
    # adding the flags and libraries directly works, but the flags are strings while the target options accept lists
    # hence converting strings to lists and appending one entry at a time
    target_link_libraries(Tasmanian_addons INTERFACE ${MPI_CXX_LIBRARIES})

    if (NOT "${MPI_CXX_INCLUDE_PATH}" STREQUAL "")
        target_include_directories(Tasmanian_addons INTERFACE ${MPI_CXX_INCLUDE_PATH})
    endif()

    if (NOT "${MPI_CXX_COMPILE_FLAGS}" STREQUAL "")
        string(REPLACE " " ";" _tsg_temp_list ${MPI_CXX_COMPILE_FLAGS})
        foreach(_flag ${_tsg_temp_list})
            target_compile_options(Tasmanian_addons INTERFACE ${_flag})
        endforeach()
    endif()

    if (NOT "${MPI_CXX_LINK_FLAGS}" STREQUAL "")
        string(REPLACE " " ";" _tsg_temp_list ${MPI_CXX_LINK_FLAGS})
        foreach(_flag ${_tsg_temp_list})
            set_target_properties(Tasmanian_addons PROPERTIES INTERFACE_LINK_OPTIONS ${_flag})
        endforeach()
    endif()
    unset(_tsg_temp_list)
    unset(_flag)

    add_executable(Tasmanian_mpitester testMPI.cpp testMPI.hpp testMPIDream.hpp)
    set_target_properties(Tasmanian_mpitester PROPERTIES OUTPUT_NAME "mpitester" CXX_EXTENSIONS OFF)
    target_link_libraries(Tasmanian_mpitester Tasmanian_addons Tasmanian_libdream)
    add_test(MPISparseGridsIO ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} mpitester ${MPIEXEC_POSTFLAGS})
    set_tests_properties(MPISparseGridsIO PROPERTIES RUN_SERIAL ON)
    # cheap models have huge random variation in time-per-sample which can cause tests to fail
    # make sure other tests are not holding resources and making things even worse
endif()

# test for non-MPI capabilities
add_executable(Tasmanian_addontester testAddons.cpp testConstructSurrogate.hpp)
set_target_properties(Tasmanian_addontester PROPERTIES OUTPUT_NAME "addontester" CXX_EXTENSIONS OFF)
target_link_libraries(Tasmanian_addontester Tasmanian_addons Tasmanian_libdream)
add_test(AddonTests addontester)
set_tests_properties(AddonTests PROPERTIES RUN_SERIAL ON) # see MPISparseGridsIO

# Windows specific support (DLL export/import directives and names)
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    if (TARGET Tasmanian_caddons)
        set_target_properties(Tasmanian_caddons PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    endif()
    if (NOT "static" IN_LIST Tasmanian_libs_type)
        set(Tasmanian_MSVC_PATH_STRING "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Release;${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Debug")
        set(Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING};${CMAKE_CURRENT_BINARY_DIR}/../DREAM/Release;${CMAKE_CURRENT_BINARY_DIR}/../DREAM/Debug;$ENV{PATH}")
        string(REPLACE ";" "\\;" Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING}")
        set_tests_properties(AddonTests PROPERTIES ENVIRONMENT "PATH=${Tasmanian_MSVC_PATH_STRING}")
    endif()
endif()

# install directives
target_include_directories(Tasmanian_addons INTERFACE $<INSTALL_INTERFACE:include>)
get_target_property(_tsg_sources Tasmanian_addons INTERFACE_SOURCES)
install(FILES ${_tsg_sources} DESTINATION include)
install(TARGETS Tasmanian_addons EXPORT "${Tasmanian_export_name}" PUBLIC_HEADER DESTINATION include)
unset(_tsg_sources)
