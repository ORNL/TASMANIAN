########################################################################
# adds a single mater target for all Tasmanian modules
########################################################################

# add master target
add_library(Tasmanian_master INTERFACE)

# add :: interface, useful when using with add_subdirectory()
add_library(Tasmanian::Tasmanian INTERFACE IMPORTED GLOBAL)
set_target_properties(Tasmanian::Tasmanian PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_master)

# add shared variant
if (NOT "${Tasmanian_libs_type}" STREQUAL "STATIC_ONLY")
    add_library(Tasmanian_shared INTERFACE)
    target_link_libraries(Tasmanian_shared INTERFACE "Tasmanian_libdream_shared;Tasmanian_addons")
    install(TARGETS Tasmanian_shared EXPORT "${Tasmanian_export_name}")

    add_library(Tasmanian::Tasmanian_shared INTERFACE IMPORTED GLOBAL)
    set_target_properties(Tasmanian::Tasmanian_shared PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_shared)
endif()

# add static variant and link the master to the first available from {static, shared}
if (NOT "${Tasmanian_libs_type}" STREQUAL "SHARED_ONLY")
    add_library(Tasmanian_static INTERFACE)
    target_link_libraries(Tasmanian_static INTERFACE "Tasmanian_libdream_static;Tasmanian_addons")
    install(TARGETS Tasmanian_static EXPORT "${Tasmanian_export_name}")

    add_library(Tasmanian::Tasmanian_static INTERFACE IMPORTED GLOBAL)
    set_target_properties(Tasmanian::Tasmanian_static PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_static)

    target_link_libraries(Tasmanian_master INTERFACE Tasmanian_static)
else()
    target_link_libraries(Tasmanian_master INTERFACE Tasmanian_shared)
endif()

if (TARGET Tasmanian_libfortran90_shared)
    target_link_libraries(Tasmanian_shared INTERFACE Tasmanian_libfortran90_shared)
endif()
if (TARGET Tasmanian_libfortran90_static)
    target_link_libraries(Tasmanian_static INTERFACE Tasmanian_libfortran90_static)
endif()

install(TARGETS Tasmanian_master EXPORT "${Tasmanian_export_name}")

# add executable that has the sole purpose of testing the master target
add_executable(tasmanian_version "${CMAKE_CURRENT_SOURCE_DIR}/Testing/tasmanian_version.cpp")
target_link_libraries(tasmanian_version Tasmanian::Tasmanian)
