########################################################################
# adds a single mater target for all Tasmanian modules
########################################################################

# add master target
add_library(Tasmanian_master INTERFACE)

# add :: interface, useful when using with add_subdirectory()
add_library(Tasmanian::Tasmanian INTERFACE IMPORTED GLOBAL)
set_target_properties(Tasmanian::Tasmanian PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_master)

# add Tasmanian_shared and Tasmanian_static, associate with Tasmanian::shared and Tasmanian::static
foreach(_tsglibtype ${Tasmanian_libs_type})
    if (TARGET Tasmanian_libdream_${_tsglibtype})
        add_library(Tasmanian_${_tsglibtype} INTERFACE)
        target_link_libraries(Tasmanian_${_tsglibtype} INTERFACE "Tasmanian_libdream_${_tsglibtype};Tasmanian_addons")
        install(TARGETS Tasmanian_${_tsglibtype} EXPORT "${Tasmanian_export_name}")

        add_library(Tasmanian::${_tsglibtype} INTERFACE IMPORTED GLOBAL)
        set_target_properties(Tasmanian::${_tsglibtype} PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_${_tsglibtype})
    endif()
    if (TARGET Tasmanian_libfortran90_${_tsglibtype})
        add_library(Tasmanian::Fortran::${_tsglibtype} INTERFACE IMPORTED GLOBAL)
        target_link_libraries(Tasmanian::Fortran::${_tsglibtype} INTERFACE Tasmanian_libfortran90_${_tsglibtype})
        target_link_libraries(Tasmanian::Fortran::${_tsglibtype} INTERFACE Tasmanian_libfortran03_${_tsglibtype})
    endif()
endforeach()
unset(_tsglibtype)

target_link_libraries(Tasmanian_master INTERFACE Tasmanian_${Tasmanian_lib_default})

install(TARGETS Tasmanian_master EXPORT "${Tasmanian_export_name}")

if (Tasmanian_ENABLE_FORTRAN)
    add_library(Tasmanian::Fortran INTERFACE IMPORTED GLOBAL)
    set_target_properties(Tasmanian::Fortran PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_libfortran90_${Tasmanian_lib_default})

    set_target_properties(Tasmanian::Fortran PROPERTIES INTERFACE_LINK_LIBRARIES Tasmanian_libfortran03_${Tasmanian_lib_default})
endif()

# add executable that has the sole purpose of testing the master target
add_executable(tasmanian_version "${CMAKE_CURRENT_SOURCE_DIR}/Testing/tasmanian_version.cpp")
target_link_libraries(tasmanian_version Tasmanian::Tasmanian)
