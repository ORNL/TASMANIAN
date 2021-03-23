########################################################################
# Python module, examples and tests
#
# Configure the .in. files using three stage configuration
#   stage 1: everything (build, test, etc) should work in the build folder
#   stage 2: everything should work after make install
#
########################################################################

if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(Tasmanian_python_default_path "${CMAKE_INSTALL_PREFIX}/python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}/site-packages")
else()
    set(Tasmanian_python_default_path "${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")
endif()

if (NOT "${Tasmanian_python_install_path}" STREQUAL "${Tasmanian_python_default_path}")
    set(Tasmanian_python_install_path "${Tasmanian_python_default_path}" CACHE INTERNAL "install path for the Python module (version specific)" FORCE)
    set(Tasmanian_PYTHONPATH "${CMAKE_INSTALL_PREFIX}/share/Tasmanian/python/" CACHE PATH "install path for the Python module (version independent)" FORCE)
endif()

# CMake has permissive implementation of boolean true, which includes
# "On", "True" and "Yes" with different cases;
# Here we convert boolen true to "ON" so it can
# be safely passed to Python as a string and interpreted via
# python_boolean = ("@cmake_boolean@" == "ON")
foreach(_tsg_opt Tasmanian_ENABLE_BLAS Tasmanian_ENABLE_CUDA Tasmanian_ENABLE_HIP Tasmanian_ENABLE_DPCPP)
    if (${_tsg_opt})
        set(${_tsg_opt} "ON")
    endif()
endforeach()
unset(_tsg_opt)


########################################################################
# Stage 1: Build folder paths
########################################################################
set(Tasmanian_string_python_hashbang "${PYTHON_EXECUTABLE}")

if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
# windows puts the temp .dll files in Release or Debug subfolder, as opposed to directly in ${CMAKE_CURRENT_BINARY_DIR}
# add different configure file for each build type and copy the appropriate one with a custom command
    foreach(_tsg_configtype Release Debug)
        set(Tasmanian_libsparsegrid_path "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/${_tsg_configtype}/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
        set(Tasmanian_libdream_path      "${CMAKE_CURRENT_BINARY_DIR}/../DREAM/${_tsg_configtype}/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniandream${CMAKE_SHARED_LIBRARY_SUFFIX}")
        set(Tasmanian_libcaddons_path    "${CMAKE_CURRENT_BINARY_DIR}/../Addons/${_tsg_configtype}/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniancaddons${CMAKE_SHARED_LIBRARY_SUFFIX}")
        # all libraries that are needed transitively must be included here, currently on the sparse grid library is needed
        set(Tasmanian_libpaths           "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/${_tsg_configtype}/")

        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/testConfigureData.in.py"  "${CMAKE_CURRENT_BINARY_DIR}/testConfigureData${_tsg_configtype}.py")
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianConfig.in.py"    "${CMAKE_CURRENT_BINARY_DIR}/TasmanianConfig${_tsg_configtype}.py")

        list(APPEND Tasmanian_python_files  "${CMAKE_CURRENT_BINARY_DIR}/testConfigureData${_tsg_configtype}.py"
                                            "${CMAKE_CURRENT_BINARY_DIR}/TasmanianConfig${_tsg_configtype}.py")
    endforeach()
    unset(_tsg_configtype)

    add_custom_target(Tasmanian_python_interface ALL DEPENDS ${Tasmanian_python_files})
    unset(Tasmanian_python_files)

    add_custom_command(TARGET Tasmanian_python_interface PRE_BUILD
                       COMMAND "${CMAKE_COMMAND}"
                       ARGS -E copy ${CMAKE_CURRENT_BINARY_DIR}/TasmanianConfig$<CONFIG>.py ${CMAKE_CURRENT_BINARY_DIR}/TasmanianConfig.py
                       COMMENT "Copying Python module for config $<CONFIG>")
    add_custom_command(TARGET Tasmanian_python_interface PRE_BUILD
                       COMMAND "${CMAKE_COMMAND}"
                       ARGS -E copy ${CMAKE_CURRENT_BINARY_DIR}/testConfigureData$<CONFIG>.py ${CMAKE_CURRENT_BINARY_DIR}/testConfigureData.py
                       COMMENT "Copying Python module for config $<CONFIG>")
else()
    set(Tasmanian_libsparsegrid_path "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(Tasmanian_libdream_path      "${CMAKE_CURRENT_BINARY_DIR}/../DREAM/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniandream${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(Tasmanian_libcaddons_path    "${CMAKE_CURRENT_BINARY_DIR}/../Addons/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniancaddons${CMAKE_SHARED_LIBRARY_SUFFIX}")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/testConfigureData.in.py"  "${CMAKE_CURRENT_BINARY_DIR}/testConfigureData.py")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianConfig.in.py"    "${CMAKE_CURRENT_BINARY_DIR}/TasmanianConfig.py")
endif()

set(Tasmanian_python_example_import "sys.path.append(\"${CMAKE_CURRENT_BINARY_DIR}\")\n")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/example_sparse_grids.in.py" "${CMAKE_CURRENT_BINARY_DIR}/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/example_dream.in.py" "${CMAKE_CURRENT_BINARY_DIR}/example_dream.py") # also uses Tasmanian_string_python_hashbang

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/sandbox.py" "${CMAKE_CURRENT_BINARY_DIR}/sandbox.py") # only uses Tasmanian_string_python_hashbang

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/GaussPattersonRule.table"  "${CMAKE_CURRENT_BINARY_DIR}/GaussPattersonRule.table" COPYONLY) # needed for testing

# See the matlab CMakeText.txt on how to copy multiple scripts
set(Tasmanian_python_test_files TasmanianSG.py TasmanianAddons.py Tasmanian.py
                                TasmanianDREAM.py TasmanianDreamLikely.py TasmanianDreamState.py TasmanianDreamSampler.py
                                testTSG.py testCommon.py testBasicIO.py testAcceleration.py testExceptions.py testMakeUpdate.py
                                testRefinement.py testUnstructuredData.py testMisc.py testAddons.py testDream.py
                                example_sparse_grids_01.py example_sparse_grids_02.py example_sparse_grids_03.py
                                example_sparse_grids_04.py example_sparse_grids_05.py example_sparse_grids_06.py
                                example_sparse_grids_07.py example_sparse_grids_08.py example_sparse_grids_09.py
                                example_sparse_grids_10.py
                                example_dream_01.py example_dream_02.py example_dream_03.py example_dream_04.py example_dream_05.py)

foreach(Tasmanian_python_testing_file ${Tasmanian_python_test_files})
    add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${Tasmanian_python_testing_file}"
                       COMMAND "${CMAKE_COMMAND}"
                       ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_python_testing_file} ${CMAKE_CURRENT_BINARY_DIR}//${Tasmanian_python_testing_file}
                       DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_python_testing_file}"
                       COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_python_testing_file}")
    list(APPEND Tasmanian_python_testing_files_stage1 "${CMAKE_CURRENT_BINARY_DIR}/${Tasmanian_python_testing_file}")
endforeach()
unset(Tasmanian_python_testing_file)
add_custom_target(Tasmanian_python_testing ALL DEPENDS "${Tasmanian_python_testing_files_stage1}")


########################################################################
# Stage 2: Install folder paths
########################################################################
set(Tasmanian_shared_lib_destination "lib")
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows") # windows puts the .dll files in bin, as opposed to lib
    set(Tasmanian_shared_lib_destination "bin")
endif()

set(Tasmanian_libsparsegrid_path "${Tasmanian_final_install_path}/${Tasmanian_shared_lib_destination}/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniansparsegrid${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(Tasmanian_libdream_path      "${Tasmanian_final_install_path}/${Tasmanian_shared_lib_destination}/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniandream${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(Tasmanian_libcaddons_path    "${Tasmanian_final_install_path}/${Tasmanian_shared_lib_destination}/${CMAKE_SHARED_LIBRARY_PREFIX}tasmaniancaddons${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(Tasmanian_libpaths           "${Tasmanian_final_install_path}/bin/")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianConfig.in.py"  "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.py")


set(Tasmanian_python_example_import "sys.path.append(\"${Tasmanian_python_install_path}\")\n")
if (SKBUILD)
    set(Tasmanian_python_example_import "") # if installed through pip, no path is needed
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/example_sparse_grids.in.py" "${CMAKE_CURRENT_BINARY_DIR}/configured/example_sparse_grids.py") # also uses Tasmanian_string_python_hashbang
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/example_dream.in.py" "${CMAKE_CURRENT_BINARY_DIR}/configured/example_dream.py") # also uses Tasmanian_string_python_hashbang


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
add_test(NAME PythonAddons       COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py" TestTasmanian.testAddons)
add_test(NAME PythonDream        COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/testTSG.py" TestTasmanian.testDream)
set_tests_properties(PythonIO PROPERTIES RUN_SERIAL ON) # IO and Refine read/write to potentially the same files
set_tests_properties(PythonRefine PROPERTIES RUN_SERIAL ON)
Tasmanian_set_test_properties(TESTS PythonIO PythonAcceleration PythonExceptions PythonMakeUpdate PythonRefine PythonLearning PythonMisc PythonAddons PythonDream)


########################################################################
# Windows system paths
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(Tasmanian_MSVC_PATH_STRING "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Release;${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Debug")
    set(Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING};${CMAKE_CURRENT_BINARY_DIR}/../DREAM/Release;${CMAKE_CURRENT_BINARY_DIR}/../DREAM/Debug")
    set(Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING};${CMAKE_CURRENT_BINARY_DIR}/../Addons/Release;${CMAKE_CURRENT_BINARY_DIR}/../Addons/Debug;$ENV{PATH}")
    string(REPLACE ";" "\\;" Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING}")
    set_tests_properties(PythonIO PythonAcceleration PythonExceptions PythonMakeUpdate PythonRefine PythonLearning PythonMisc PythonAddons PythonDream
                         PROPERTIES ENVIRONMENT "PATH=${Tasmanian_MSVC_PATH_STRING}")
endif()

########################################################################
# Installation
########################################################################
macro(Tasmanian_python_install)
    cmake_parse_arguments(_tsg_py "" "" "FILES" ${ARGN})
    foreach(_tsg_py_file ${_tsg_py_FILES})
        get_filename_component(_tsg_py_name ${_tsg_py_file} NAME)

        # install in the version specific site-package location
        install(FILES ${_tsg_py_file} DESTINATION "${Tasmanian_python_install_path}"
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

        # install in the version independent location (note that pip-install doesn't handle sym-links)
        install(FILES ${_tsg_py_file} DESTINATION "${Tasmanian_PYTHONPATH}"
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

        if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows" AND Tasmanian_windows_virtual)
            # windows virtual environments use a different path and naming convention, annoying ...
            install(FILES ${_tsg_py_file} DESTINATION "${CMAKE_INSTALL_PREFIX}/Lib/site-packages/"
                    PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
        elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin" AND Tasmanian_osx_framework)
            install(FILES ${_tsg_py_file} DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python/site-packages/"
                    PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
        endif()
    endforeach()
endmacro(Tasmanian_python_install)

Tasmanian_python_install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSG.py"
                               "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianAddons.py"
                               "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianDREAM.py"
                               "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianDreamLikely.py"
                               "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianDreamState.py"
                               "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianDreamSampler.py"
                               "${CMAKE_CURRENT_BINARY_DIR}/configured/TasmanianConfig.py"
                               "${CMAKE_CURRENT_SOURCE_DIR}/Tasmanian.py")

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
        DESTINATION "share/Tasmanian/examples/"
        FILES_MATCHING PATTERN "example_*.py"
        PATTERN "*.in.*" EXCLUDE
        PATTERN "*PipInstaller*" EXCLUDE
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/example_sparse_grids.py"
              "${CMAKE_CURRENT_BINARY_DIR}/configured/example_dream.py"
        DESTINATION "share/Tasmanian/examples/"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
