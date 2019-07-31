########################################################################
# Addon templates and related tests and interface options
########################################################################

macro(Tasmanian_addon_sources)
    cmake_parse_arguments(_tsg_sources "" "PATH" "FILES" ${ARGN} )
    target_include_directories(Tasmanian_addons INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_tsg_sources_PATH}>)
    foreach(_tsg_file ${_tsg_sources_FILES})
        target_sources(Tasmanian_addons INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_tsg_sources_PATH}/${_tsg_file}>
                                                  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include/${_tsg_file}>)
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
                              tsgLoadNeededValues.hpp
                              tsgCandidateManager.hpp
                              tsgConstructSurrogate.hpp
                              tsgMPIConstructGrid.hpp
                              TasmanianAddons.hpp)

if (NOT Tasmanian_ENABLE_OPENMP) # using threading without OpenMP
    target_link_libraries(Tasmanian_addons INTERFACE ${CMAKE_THREAD_LIBS_INIT})
endif()

# The Tasmanian MPI capabilities are templated into the Addons
if (Tasmanian_ENABLE_MPI)
    target_link_libraries(Tasmanian_addons INTERFACE ${MPI_CXX_LIBRARIES})

    if (DEFINED MPI_CXX_INCLUDE_PATH)
        target_include_directories(Tasmanian_addons INTERFACE "${MPI_CXX_INCLUDE_PATH}")
    endif()

    if(DEFINED MPI_CXX_COMPILE_FLAGS)
        target_compile_options(Tasmanian_addons INTERFACE "${MPI_CXX_COMPILE_FLAGS}")
    endif()

    if(DEFINED MPI_CXX_LINK_FLAGS)
        set_target_properties(Tasmanian_addons PROPERTIES INTERFACE_LINK_OPTIONS "${MPI_CXX_LINK_FLAGS}")
    endif()

    add_executable(Tasmanian_mpitester testMPI.cpp testMPI.hpp)
    set_target_properties(Tasmanian_mpitester PROPERTIES OUTPUT_NAME "mpitester")
    target_link_libraries(Tasmanian_mpitester Tasmanian_addons Tasmanian_libdream)
    add_test(MPISparseGridsIO ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} mpitester ${MPIEXEC_POSTFLAGS})
endif()

# test for non-MPI capabilities
add_executable(Tasmanian_addontester testAddons.cpp testConstructSurrogate.hpp)
set_target_properties(Tasmanian_addontester PROPERTIES OUTPUT_NAME "addontester")
target_link_libraries(Tasmanian_addontester Tasmanian_addons Tasmanian_libdream)
add_test(AddonTests addontester)

# install directives
target_include_directories(Tasmanian_addons INTERFACE $<INSTALL_INTERFACE:include>)
get_target_property(tsg_sources Tasmanian_addons INTERFACE_SOURCES)
install(FILES ${tsg_sources} DESTINATION include)
install(TARGETS Tasmanian_addons EXPORT "${Tasmanian_export_name}" PUBLIC_HEADER DESTINATION include)
