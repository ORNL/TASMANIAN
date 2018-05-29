cmake_minimum_required(VERSION 3.5)

########################################################################
# Python module, examples and tests
#
# Configure the .in. files using three stage configuration
#   stage 0: (optional) keeps files up-to-date for the simple GNU-Make
#   stage 1: everything (build, test, etc) should work in the build folder
#   stage 2: everything should work after make install
#
########################################################################

set(Tasmanian_python_install_path "${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")

########################################################################
# Stage 0: for support of simple GNU-Make
#          to be used only for development
########################################################################
if (Tasmanian_DEVELOPMENT_BACKWARDS)
    set(Tasmanian_libsparsegrid_path "libtasmaniansparsegrid.dll")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.in.py" "${CMAKE_CURRENT_SOURCE_DIR}/../Config/AltBuildSystems/TasmanianSG.windows.py")

    set(Tasmanian_libsparsegrid_path "./libtasmaniansparsegrid.so")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.in.py" "${CMAKE_CURRENT_SOURCE_DIR}/../Config/AltBuildSystems/TasmanianSG.py")

    set(Tasmanian_string_python_hashbang "/usr/bin/env python")
    set(Tasmanian_cmake_synctest_enable "False") # disable tests that check if the library reports the same options as given by the cmake variables
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/testTSG.in.py" "${CMAKE_CURRENT_SOURCE_DIR}/../Config/AltBuildSystems/testTSG.py") # also using Tasmanian_libsparsegrid_path

    set(Tasmanian_python_example_import "#")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/example_sparse_grids.in.py" "${CMAKE_CURRENT_SOURCE_DIR}/../Config/AltBuildSystems/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang
endif()


########################################################################
# Stage 1: Build folder paths
########################################################################
set(Tasmanian_libsparsegrid_path "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
# windows puts the temp .dll files in ${CMAKE_BUILD_TYPE} subfolder, as opposed to directly in ${CMAKE_BINARY_DIR}
    set(Tasmanian_libsparsegrid_path "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/${CMAKE_BUILD_TYPE}/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.in.py" "${CMAKE_CURRENT_BINARY_DIR}/TasmanianSG.py")

set(Tasmanian_string_python_hashbang "${PYTHON_EXECUTABLE}")
set(Tasmanian_cmake_synctest_enable "True") # enable tests that check if the library reports the same options as given by the cmake variables
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/testTSG.in.py" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py") # also uses Tasmanian_libsparsegrid_path

set(Tasmanian_python_example_import "sys.path.append(\"${CMAKE_CURRENT_BINARY_DIR}\")\n")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/example_sparse_grids.in.py" "${CMAKE_CURRENT_BINARY_DIR}/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/sandbox.py" "${CMAKE_CURRENT_BINARY_DIR}/sandbox.py") # only uses Tasmanian_string_python_hashbang

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/GaussPattersonRule.table"  "${CMAKE_CURRENT_BINARY_DIR}/GaussPattersonRule.table" COPYONLY) # needed for testing


########################################################################
# Stage 2: Install folder paths
########################################################################
set(Tasmanian_libsparsegrid_path "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
# windows puts the .dll files in bin, as opposed to lib
    set(Tasmanian_libsparsegrid_path "${CMAKE_INSTALL_PREFIX}/bin/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.in.py" "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianSG.py")

set(Tasmanian_python_example_import "sys.path.append(\"${Tasmanian_python_install_path}\")\n")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/example_sparse_grids.in.py" "${CMAKE_CURRENT_BINARY_DIR}/configured/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang


########################################################################
# Testing
########################################################################
add_test(NAME PythonIO           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/testTSG.py TestTasmanian.testBasicIO)
add_test(NAME PythonAcceleration COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/testTSG.py TestTasmanian.testAcceleratedEvaluate)
add_test(NAME PythonExceptions   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/testTSG.py TestTasmanian.testBasicException)
add_test(NAME PythonMakeUpdate   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/testTSG.py TestTasmanian.testFullCoverageA)
add_test(NAME PythonRefine       COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/testTSG.py TestTasmanian.testFullCoverageB)
add_test(NAME PythonLearning     COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/testTSG.py TestTasmanian.testFullCoverageC)
add_test(NAME PythonMisc         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/testTSG.py TestTasmanian.testFullCoverageZ)
if (Tasmanian_TESTS_OMP_NUM_THREADS GREATER 0)
    set_tests_properties(PythonIO PythonAcceleration PythonExceptions PythonMakeUpdate PythonRefine PythonLearning PythonMisc
        PROPERTIES
        PROCESSORS ${Tasmanian_TESTS_OMP_NUM_THREADS}
        ENVIRONMENT OMP_NUM_THREADS=${Tasmanian_TESTS_OMP_NUM_THREADS})
endif()

########################################################################
# Installation
########################################################################
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianSG.py"
        DESTINATION "${Tasmanian_python_install_path}"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)

# Create symlink for backward compatibility
install(CODE "execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${Tasmanian_python_install_path} ${CMAKE_INSTALL_PREFIX}/python )" )
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/example_sparse_grids.py"
        DESTINATION "examples"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
