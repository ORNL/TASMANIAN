########################################################################
# Testing: post install, make test_install
# checks if the executables run and if the examples compile and run
########################################################################
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Testing/test_post_install.in.sh" "${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh")
if (Tasmanian_ENABLE_FORTRAN)
    add_custom_target(test_install COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh"
                                           "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
else()
    add_custom_target(test_install COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh"
                                           "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
endif()

########################################################################
# Install exports, package-config files, and examples cmake file
########################################################################
install(EXPORT "${Tasmanian_export_name}" DESTINATION "lib/${CMAKE_PROJECT_NAME}" FILE "${CMAKE_PROJECT_NAME}.cmake")

configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianConfig.in.cmake"
                              "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.cmake"
                              INSTALL_DESTINATION "lib/Tasmanian/")
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfigVersion.cmake"
                                 COMPATIBILITY AnyNewerVersion)
# not sure why it is necessary to explicitly install TasmanianConfig.cmake, INSTALL_DESTINATION above doesn't work
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.cmake"
        DESTINATION "lib/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfigVersion.cmake"
        DESTINATION "lib/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

# TasmanianMake.in for GNU Make include/link
get_target_property(Tasmanian_list Tasmanian_libsparsegrid INTERFACE_LINK_LIBRARIES)
foreach(Tasmanian_lib_ ${Tasmanian_list})
    set(Tasmanian_libs "${Tasmanian_libs} ${Tasmanian_lib_}")
endforeach()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianMakefile.in" "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianMakefile.in" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianMakefile.in"
        DESTINATION "share/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

# cmake file for the examples, to be used post-install
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeLists.examples.txt" "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt"
        DESTINATION "share/Tasmanian/examples/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)


# configure environment shell script that can be sourced to set path and lib path
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianENVsetup.in.sh" "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh"
        DESTINATION "share/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
