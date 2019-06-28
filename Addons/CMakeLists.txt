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
                              TasmanianAddons.hpp)

# install directives
target_include_directories(Tasmanian_addons INTERFACE $<INSTALL_INTERFACE:include>)
get_target_property(tsg_sources Tasmanian_addons INTERFACE_SOURCES)
install(FILES ${tsg_sources} DESTINATION include)
install(TARGETS Tasmanian_addons EXPORT "${Tasmanian_export_name}" PUBLIC_HEADER DESTINATION include)
