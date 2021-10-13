########################################################################
# Testing: post install, make test_install
# checks if the executables run and if the examples compile and run
########################################################################
set(Tasmanian_langs "CXX")
set(Tasmanian_compilers  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
if (Tasmanian_ENABLE_FORTRAN)
    set(Tasmanian_langs "${Tasmanian_langs} Fortran")
    set(Tasmanian_compilers  "${Tasmanian_compilers} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Testing/test_post_install.in.sh" "${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Testing/CMakeLists.test.cmake" "${CMAKE_CURRENT_BINARY_DIR}/configured/testing/CMakeLists.txt" @ONLY)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tasmanian_test_install)
add_custom_target(test_install COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_post_install.sh"
                               WORKING_DIRECTORY tasmanian_test_install)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/testing/CMakeLists.txt"
        DESTINATION "share/Tasmanian/testing/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

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

# master header Tasmanian.hpp
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/Config/Tasmanian.hpp"
        DESTINATION "include"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/Config/Tasmanian.h"
        DESTINATION "include"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)


# TasmanianMake.in for GNU Make include/link
get_target_property(Tasmanian_list Tasmanian_dependencies INTERFACE_LINK_LIBRARIES)
foreach(Tasmanian_lib_ ${Tasmanian_list})
    if (NOT TARGET ${Tasmanian_lib_})
        set(Tasmanian_libs "${Tasmanian_libs} ${Tasmanian_lib_}")
    endif()
endforeach()
unset(Tasmanian_lib_)
unset(Tasmanian_list)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianMakefile.in" "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianMakefile.in" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianMakefile.in"
        DESTINATION "share/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)


# configure the cmake file for the examples, to be used post-install
if (BUILD_SHARED_LIBS)
    set(Tasmanian_components " SHARED")
else()
    set(Tasmanian_components " STATIC")
endif()

foreach(_comp OPENMP BLAS PYTHON CUDA HIP DPCPP MAGMA FORTRAN MPI)
    if (Tasmanian_ENABLE_${_comp})
        set(Tasmanian_components "${Tasmanian_components} ${_comp}")
    endif()
endforeach()
unset(_comp)
if (NOT "${Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
    set(Tasmanian_components "${Tasmanian_components} MATLAB")
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeLists.examples.cmake" "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/CMakeLists.txt"
        DESTINATION "share/Tasmanian/examples/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)


# configure environment shell script that can be sourced to set path and lib path
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config/TasmanianENVsetup.in.sh" "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianENVsetup.sh"
        DESTINATION "share/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
