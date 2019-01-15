########################################################################
# Python module, examples and tests
#
# Configure the .in. files using three stage configuration
#   stage 0: (optional) keeps files up-to-date for the simple GNU-Make
#   stage 1: everything (build, test, etc) should work in the build folder
#   stage 2: everything should work after make install
#
########################################################################

if (NOT DEFINED Tasmanian_python_install_path)
    # calling "make install" invalidates PYTHON_VERSION_MAJOR and PYTHON_VERSION_MINOR
    # thus, we cache the string on first call of cmake and the user could edit this variable
    set(Tasmanian_python_install_path "${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages" CACHE PATH "install path for the Python module")
endif()


########################################################################
# Stage 0: for support of simple GNU-Make
#          to be used only for development
########################################################################
if (Tasmanian_DEVELOPMENT_BACKWARDS)
    set(Tasmanian_libsparsegrid_path "./libtasmaniansparsegrid.so")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.in.py" "${CMAKE_CURRENT_SOURCE_DIR}/../Config/AltBuildSystems/TasmanianSG.py")

    set(Tasmanian_string_python_hashbang "/usr/bin/env python")
    set(Tasmanian_python_example_import "#")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/example_sparse_grids.in.py" "${CMAKE_CURRENT_SOURCE_DIR}/../Config/AltBuildSystems/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang
endif()


########################################################################
# Stage 1: Build folder paths
########################################################################
set(Tasmanian_libsparsegrid_path "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.in.py" "${CMAKE_CURRENT_BINARY_DIR}/TasmanianSG.py")

set(Tasmanian_string_python_hashbang "${PYTHON_EXECUTABLE}")

if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
# windows puts the temp .dll files in Release or Debug subfolder, as opposed to directly in ${CMAKE_CURRENT_BINARY_DIR}
# add different configure file for each build type and copy the appropriate one with a custom command
    set(Tasmanian_libsparsegrid_path "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Release/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.in.py" "${CMAKE_CURRENT_BINARY_DIR}/TasmanianSG_Release.py")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/testConfigureData.in.py" "${CMAKE_CURRENT_BINARY_DIR}/testConfigureData_Release.py")
    set(Tasmanian_libsparsegrid_path "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Debug/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.in.py" "${CMAKE_CURRENT_BINARY_DIR}/TasmanianSG_Debug.py")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/testConfigureData.in.py" "${CMAKE_CURRENT_BINARY_DIR}/testConfigureData_Debug.py")
    add_custom_target(Tasmanian_python_interface ALL DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.in.py" "${CMAKE_CURRENT_BINARY_DIR}/testConfigureData.py")
    add_custom_command(TARGET Tasmanian_python_interface PRE_BUILD
                       COMMAND "${CMAKE_COMMAND}"
                       ARGS -E copy ${CMAKE_CURRENT_BINARY_DIR}/TasmanianSG_$<CONFIG>.py ${CMAKE_CURRENT_BINARY_DIR}/TasmanianSG.py
                       COMMENT "Copying Python module for config $<CONFIG>")
    add_custom_command(TARGET Tasmanian_python_interface PRE_BUILD
                       COMMAND "${CMAKE_COMMAND}"
                       ARGS -E copy ${CMAKE_CURRENT_BINARY_DIR}/testConfigureData_$<CONFIG>.py ${CMAKE_CURRENT_BINARY_DIR}/testConfigureData.py
                       COMMENT "Copying Python module for config $<CONFIG>")
endif()

set(Tasmanian_python_example_import "sys.path.append(\"${CMAKE_CURRENT_BINARY_DIR}\")\n")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/example_sparse_grids.in.py" "${CMAKE_CURRENT_BINARY_DIR}/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/sandbox.py" "${CMAKE_CURRENT_BINARY_DIR}/sandbox.py") # only uses Tasmanian_string_python_hashbang

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/GaussPattersonRule.table"  "${CMAKE_CURRENT_BINARY_DIR}/GaussPattersonRule.table" COPYONLY) # needed for testing

# See the matlab CMakeText.txt on how to copy multiple scripts
set(Tasmanian_python_test_files testTSG.py testCommon.py testBasicIO.py testAcceleration.py testExceptions.py testMakeUpdate.py testRefinement.py testUnstructuredData.py testMisc.py)

foreach(Tasmanian_python_testing_file ${Tasmanian_python_test_files})
    add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${Tasmanian_python_testing_file}"
                       COMMAND "${CMAKE_COMMAND}"
                       ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_python_testing_file} ${CMAKE_CURRENT_BINARY_DIR}//${Tasmanian_python_testing_file}
                       DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_python_testing_file}"
                       COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_python_testing_file}")
    list(APPEND Tasmanian_python_testing_files_stage1 "${CMAKE_CURRENT_BINARY_DIR}/${Tasmanian_python_testing_file}")
endforeach()
add_custom_target(Tasmanian_python_testing ALL DEPENDS "${Tasmanian_python_testing_files_stage1}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/testConfigureData.in.py"  "${CMAKE_CURRENT_BINARY_DIR}/testConfigureData.py")


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
add_test(NAME PythonIO           COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py" TestTasmanian.testBasicIO)
add_test(NAME PythonAcceleration COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py" TestTasmanian.testAcceleratedEvaluate)
add_test(NAME PythonExceptions   COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py" TestTasmanian.testBasicException)
add_test(NAME PythonMakeUpdate   COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py" TestTasmanian.testAMakeUpdate)
add_test(NAME PythonRefine       COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py" TestTasmanian.testBRefinement)
add_test(NAME PythonLearning     COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py" TestTasmanian.testCUnsructuredData)
add_test(NAME PythonMisc         COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py" TestTasmanian.testZMisc)
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
install(CODE "execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${Tasmanian_python_install_path} ${CMAKE_INSTALL_PREFIX}/share/Tasmanian/python )" )
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/example_sparse_grids.py"
        DESTINATION "share/Tasmanian/examples/"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
