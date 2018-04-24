########################################################################
# Configure the .in. files
# using three stage configuration
#   stage 0: (optional) configure files with default values for simple GNU-Make
#            by default, GNU-Make doesn't support any advanced options
#   stage 1: everything (build, test, etc) should work in the build folder
#   stage 2: everything should work after make install
########################################################################

# stage 0: for support of simple GNU-Make, ignore for release and intermediate builds
if (Tasmanian_DEVELOPMENT_BACKWARDS)
    set(Tasmanian_libsparsegrid_path "libtasmaniansparsegrid.dll")
    configure_file("${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py" "${PROJECT_SOURCE_DIR}/Config/AltBuildSystems/TasmanianSG.windows.py")

    set(Tasmanian_libsparsegrid_path "./libtasmaniansparsegrid.so")
    configure_file("${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py" "${PROJECT_SOURCE_DIR}/Config/AltBuildSystems/TasmanianSG.py")

    set(Tasmanian_string_python_hashbang "/usr/bin/env python")
    set(Tasmanian_cmake_synctest_enable "False") # disable tests that check if the library reports the same options as given by the cmake variables
    configure_file("${PROJECT_SOURCE_DIR}/Testing/testTSG.in.py" "${PROJECT_SOURCE_DIR}/Config/AltBuildSystems/testTSG.py") # also using Tasmanian_libsparsegrid_path

    set(Tasmanian_python_example_import "#")
    configure_file("${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.in.py" "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang
endif()


# stage 1: build folder paths
configure_file("${PROJECT_SOURCE_DIR}/Config/TasmanianConfig.in.hpp"  "${CMAKE_BINARY_DIR}/configured/TasmanianConfig.hpp")

configure_file("${PROJECT_SOURCE_DIR}/SparseGrids/GaussPattersonRule.table"  "${CMAKE_BINARY_DIR}/GaussPattersonRule.table" COPYONLY)

if (Tasmanian_ENABLE_PYTHON)
    set(Tasmanian_libsparsegrid_path "${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    # windows puts the temp .dll files in ${CMAKE_BUILD_TYPE} subfolder, as opposed to directly in ${CMAKE_BINARY_DIR}
        set(Tasmanian_libsparsegrid_path "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/${CMAKE_SHARED_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
    configure_file("${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py" "${CMAKE_BINARY_DIR}/TasmanianSG.py")

    set(Tasmanian_string_python_hashbang "${PYTHON_EXECUTABLE}")
    set(Tasmanian_cmake_synctest_enable "True") # enable tests that check if the library reports the same options as given by the cmake variables
    configure_file("${PROJECT_SOURCE_DIR}/Testing/testTSG.in.py" "${CMAKE_BINARY_DIR}/testTSG.py") # also uses Tasmanian_libsparsegrid_path

    set(Tasmanian_python_example_import "sys.path.append(\"${CMAKE_BINARY_DIR}\")\n")
    configure_file("${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.in.py" "${CMAKE_BINARY_DIR}/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang

    configure_file("${PROJECT_SOURCE_DIR}/Testing/sandbox.py" "${CMAKE_BINARY_DIR}/sandbox.py") # only uses Tasmanian_string_python_hashbang
endif()

if (Tasmanian_ENABLE_MATLAB)
    set(Tasmanian_tasgrid_path "${CMAKE_BINARY_DIR}/${Tasmanian_name_tasgrid}")
    set(Tasmanian_string_matlab_work ${CMAKE_BINARY_DIR}/tempMATLAB/)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tempMATLAB/)
    configure_file("${PROJECT_SOURCE_DIR}/InterfaceMATLAB/tsgGetPaths.in.m" "${CMAKE_BINARY_DIR}/matlab/tsgGetPaths.m")

    foreach(Tasmanian_matlab_file ${Tasmanian_source_matlabsg})
        configure_file("${PROJECT_SOURCE_DIR}/InterfaceMATLAB/${Tasmanian_matlab_file}" "${CMAKE_BINARY_DIR}/matlab/${Tasmanian_matlab_file}" COPYONLY)
    endforeach()
endif()


# stage 2: install folder paths
if (Tasmanian_ENABLE_PYTHON)
    set(Tasmanian_libsparsegrid_path "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    # windows puts the .dll files in bin, as opposed to lib
        set(Tasmanian_libsparsegrid_path "${CMAKE_INSTALL_PREFIX}/bin/${CMAKE_SHARED_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
    configure_file("${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py" "${CMAKE_BINARY_DIR}/install/python/TasmanianSG.py")

    set(Tasmanian_python_example_import "sys.path.append(\"${Tasmanian_python_install_path}\")\n")
    configure_file("${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.in.py" "${CMAKE_BINARY_DIR}/install/examples/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang
endif()

if (Tasmanian_ENABLE_MATLAB)
    set(Tasmanian_tasgrid_path "${CMAKE_INSTALL_PREFIX}/bin/${Tasmanian_name_tasgrid}")
    configure_file("${PROJECT_SOURCE_DIR}/InterfaceMATLAB/tsgGetPaths.in.m" "${CMAKE_BINARY_DIR}/install/matlab/tsgGetPaths.m")
endif()

# configure post-install tests, i.e., make test_install
configure_file("${PROJECT_SOURCE_DIR}/Testing/test_post_install.in.sh" "${CMAKE_BINARY_DIR}/test_post_install.sh")

# cmake file for the examples, to be used post-install
set(Tasmanian_cmake_export "${CMAKE_INSTALL_PREFIX}/config/Tasmanian.cmake")
configure_file("${PROJECT_SOURCE_DIR}/Examples/CMakeLists.in.txt" "${CMAKE_BINARY_DIR}/install/examples/CMakeLists.txt")

# environment setup
set(Tasmanian_config_line1 "export PATH=$PATH:${CMAKE_INSTALL_PREFIX}/bin/")
set(Tasmanian_config_line2 "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib/")
configure_file("${PROJECT_SOURCE_DIR}/Config/TasmanianENVsetup.in.sh" "${CMAKE_BINARY_DIR}/install/config/TasmanianENVsetup.sh")

set(Tasmanian_config_line1 "export C_INCLUDE_PATH=$C_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/")
set(Tasmanian_config_line2 "export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/")
set(Tasmanian_config_line3 "export LIBRARY_PATH=$LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib/")
configure_file("${PROJECT_SOURCE_DIR}/Config/TasmanianDEVsetup.in.sh" "${CMAKE_BINARY_DIR}/install/config/TasmanianDEVsetup.sh")
