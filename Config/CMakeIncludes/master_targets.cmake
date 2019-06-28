########################################################################
# adds a single mater target for all Tasmanian modules
########################################################################
macro(Tasmanian_add_mpi Tasmanian_target)
    if (Tasmanian_ENABLE_MPI)
        target_link_libraries(${Tasmanian_target} INTERFACE ${MPI_CXX_LIBRARIES})

        if (DEFINED MPI_CXX_INCLUDE_PATH)
            target_include_directories(${Tasmanian_target} INTERFACE "${MPI_CXX_INCLUDE_PATH}")
        endif()

        if(DEFINED MPI_CXX_COMPILE_FLAGS)
            target_compile_options(${Tasmanian_target} INTERFACE "${MPI_CXX_COMPILE_FLAGS}")
        endif()

        if(DEFINED MPI_CXX_LINK_FLAGS)
            set_target_properties(${Tasmanian_target} PROPERTIES INTERFACE_LINK_OPTIONS "${MPI_CXX_LINK_FLAGS}")
        endif()
    endif()
endmacro(Tasmanian_add_mpi)

# add master target
add_library(Tasmanian_master INTERFACE)
target_link_libraries(Tasmanian_master INTERFACE Tasmanian_addons)

# add :: interface, useful when using with add_subdirectory()
add_library(Tasmanian::Tasmanian INTERFACE IMPORTED)
set_target_properties(Tasmanian::Tasmanian PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_master)

# add shared variant
if (NOT "${Tasmanian_libs_type}" STREQUAL "STATIC_ONLY")
    add_library(Tasmanian_shared INTERFACE)
    target_link_libraries(Tasmanian_shared INTERFACE Tasmanian_libdream_shared)
    Tasmanian_add_mpi(Tasmanian_shared)
    install(TARGETS Tasmanian_shared EXPORT "${Tasmanian_export_name}")

    add_library(Tasmanian::Tasmanian_shared INTERFACE IMPORTED)
    set_target_properties(Tasmanian::Tasmanian_shared PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_shared)
endif()

# add static variant and link the master to the first available from {static, shared}
if (NOT "${Tasmanian_libs_type}" STREQUAL "SHARED_ONLY")
    add_library(Tasmanian_static INTERFACE)
    target_link_libraries(Tasmanian_static INTERFACE Tasmanian_libdream_static)
    Tasmanian_add_mpi(Tasmanian_static)
    install(TARGETS Tasmanian_static EXPORT "${Tasmanian_export_name}")

    add_library(Tasmanian::Tasmanian_static INTERFACE IMPORTED)
    set_target_properties(Tasmanian::Tasmanian_static PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_static)

    target_link_libraries(Tasmanian_master INTERFACE Tasmanian_static)
else()
    target_link_libraries(Tasmanian_master INTERFACE Tasmanian_shared)
endif()

install(TARGETS Tasmanian_master EXPORT "${Tasmanian_export_name}")

# add executable that has the sole purpose of testing the master target
add_executable(tasmanian_version "${CMAKE_CURRENT_SOURCE_DIR}/Testing/tasmanian_version.cpp")
target_link_libraries(tasmanian_version Tasmanian::Tasmanian)
