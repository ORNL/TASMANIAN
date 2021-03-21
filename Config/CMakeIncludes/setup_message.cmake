########################################################################
# Print the configuration options
########################################################################
message(STATUS "")
message(STATUS "Tasmanian ${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}${Tasmanian_version_comment}: summary of build options")
message(STATUS " -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
message(STATUS " -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}")
message(STATUS " -D CMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(STATUS " -D CMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}") # useful for Windows debugging
    message(STATUS " -D CMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}")
endif()
if (Tasmanian_ENABLE_CUDA)
    message(STATUS " -D CMAKE_CUDA_FLAGS:STRING=${CMAKE_CUDA_FLAGS}")
endif()
message(STATUS " -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}")

foreach(Tasmanian_option Tasmanian_ENABLE_OPENMP  Tasmanian_ENABLE_BLAS
                         Tasmanian_ENABLE_MPI     Tasmanian_ENABLE_PYTHON
                         Tasmanian_ENABLE_CUDA    Tasmanian_ENABLE_MAGMA
                         Tasmanian_ENABLE_HIP     Tasmanian_ENABLE_DPCPP
                         Tasmanian_ENABLE_SWIG
                         Tasmanian_ENABLE_FORTRAN Tasmanian_ENABLE_DOXYGEN)

    message(STATUS " -D ${Tasmanian_option}:BOOL=${${Tasmanian_option}}")
endforeach()

if (Tasmanian_MAGMA AND Tasmanian_MAGMA_ROOT)
    message(STATUS " -D Tasmanian_MAGMA_ROOT:PATH=${Tasmanian_MAGMA_ROOT}")
endif()
if (NOT "${Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
    message(STATUS " -D Tasmanian_MATLAB_WORK_FOLDER:PATH=${Tasmanian_MATLAB_WORK_FOLDER}")
    message(STATUS " pre-install MATLAB folder: addpath('${CMAKE_CURRENT_BINARY_DIR}/MATLAB/matlab/')")
endif()
message(STATUS "")


########################################################################
# Print final message (a bit of a hack)
# The message is written in the CMakeLists.txt, as the subdir is added last
# this ensures that the message will appear last in the install process
# do not print if USE_XSDK or Tasmanian has been imported with addsubdir
########################################################################
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/")
endif()
