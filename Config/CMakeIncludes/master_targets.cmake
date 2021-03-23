########################################################################
# adds a single mater target for all Tasmanian modules
########################################################################

# add master target
add_library(Tasmanian_master INTERFACE)
target_link_libraries(Tasmanian_master INTERFACE Tasmanian_addons)
install(TARGETS Tasmanian_master EXPORT "${Tasmanian_export_name}")

# add :: interface, useful when using with add_subdirectory()
add_library(Tasmanian::Tasmanian INTERFACE IMPORTED GLOBAL)
target_link_libraries(Tasmanian::Tasmanian INTERFACE Tasmanian_master)

if (BUILD_SHARED_LIBS) # for backwards compatibility
    add_library(Tasmanian::shared INTERFACE IMPORTED GLOBAL)
    target_link_libraries(Tasmanian::shared INTERFACE Tasmanian_master)
else()
    add_library(Tasmanian::static INTERFACE IMPORTED GLOBAL)
    target_link_libraries(Tasmanian::static INTERFACE Tasmanian_master)
endif()

if (Tasmanian_ENABLE_FORTRAN)
    add_library(Tasmanian::Fortran INTERFACE IMPORTED GLOBAL)
    target_link_libraries(Tasmanian::Fortran INTERFACE Tasmanian_libfortran90)
    target_link_libraries(Tasmanian::Fortran INTERFACE Tasmanian_libfortran03)
    if (Tasmanian_ENABLE_MPI)
        target_link_libraries(Tasmanian::Fortran INTERFACE Tasmanian_libfortranmpi03)
    endif()
endif()

# add executable that has the sole purpose of testing the master target
add_executable(tasmanian_version "${CMAKE_CURRENT_SOURCE_DIR}/Testing/tasmanian_version.cpp")
target_link_libraries(tasmanian_version Tasmanian::Tasmanian)
