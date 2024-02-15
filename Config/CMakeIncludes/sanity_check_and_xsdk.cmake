########################################################################
# sanity check and xSDK compatibility
########################################################################
if ((NOT DEFINED CMAKE_BUILD_TYPE) OR (NOT CMAKE_BUILD_TYPE))
    set(CMAKE_BUILD_TYPE Debug)
endif()

# get extra options from the ENV variables (pip-installer)
if (SKBUILD)
    if ("$ENV{Tasmanian_ENABLE_RECOMMENDED}" STREQUAL "OFF")
        set(Tasmanian_ENABLE_RECOMMENDED OFF)
    endif()
    if (NOT "$ENV{Tasmanian_ENABLE_BLAS}" STREQUAL "")
        set(Tasmanian_ENABLE_BLAS ON)
        set(BLAS_LIBRARIES   "$ENV{Tasmanian_ENABLE_BLAS}")
        set(LAPACK_LIBRARIES "$ENV{Tasmanian_ENABLE_BLAS}")
    endif()
    if (NOT "$ENV{Tasmanian_ENABLE_CUDA}" STREQUAL "")
        set(Tasmanian_ENABLE_CUDA ON)
        set(CMAKE_CUDA_COMPILER "$ENV{Tasmanian_ENABLE_CUDA}")
    endif()
    if (NOT "$ENV{Tasmanian_ENABLE_HIP}" STREQUAL "" OR NOT "$ENV{Tasmanian_ENABLE_ROCM}" STREQUAL "")
        set(Tasmanian_ENABLE_HIP ON)
        set(Tasmanian_ENABLE_ROCM "$ENV{Tasmanian_ENABLE_ROCM}")
    endif()
    if (NOT "$ENV{Tasmanian_ENABLE_DPCPP}" STREQUAL "")
        set(Tasmanian_ENABLE_DPCPP ON)
        set(CMAKE_CXX_COMPILER "$ENV{Tasmanian_ENABLE_DPCPP}")
    endif()
    if (NOT "$ENV{Tasmanian_ENABLE_MAGMA}" STREQUAL "")
        set(Tasmanian_ENABLE_MAGMA ON)
        set(MAGMA_ROOT "$ENV{Tasmanian_ENABLE_MAGMA}")
    endif()
    if (NOT "$ENV{Tasmanian_MAGMA_DOWNLOAD}" STREQUAL "")
        set(Tasmanian_MAGMA_DOWNLOAD ON)
    endif()
    if (NOT "$ENV{Tasmanian_ENABLE_MPI}" STREQUAL "")
        set(Tasmanian_ENABLE_MPI ON)
        set(MPI_CXX_COMPILER "$ENV{Tasmanian_ENABLE_MPI}")
    endif()
    if (NOT "$ENV{Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
        set(Tasmanian_MATLAB_WORK_FOLDER "$ENV{Tasmanian_MATLAB_WORK_FOLDER}")
    endif()
endif()

# check for Fortran, note that enable_language always gives FATAL_ERROR if the compiler is missing
if (Tasmanian_ENABLE_FORTRAN)
    enable_language(Fortran)
    Tasmanian_compiler_type(COMPILER ${CMAKE_Fortran_COMPILER} TYPE "ifort" RESULT Tasmanian_ifort_compiler)
endif()

# swig requires Fortran and cannot handle both types of libs
if (Tasmanian_ENABLE_SWIG AND NOT Tasmanian_ENABLE_FORTRAN)
    message(FATAL_ERROR "Tasmanian_ENABLE_SWIG=ON requires Tasmanian_ENABLE_FORTRAN=ON")
endif()
if (Tasmanian_ENABLE_FORTRAN90 AND NOT Tasmanian_ENABLE_FORTRAN)
    message(FATAL_ERROR "Tasmanian_ENABLE_FORTRAN90=ON requires Tasmanian_ENABLE_FORTRAN=ON")
endif()

# OpenMP setup
if ((Tasmanian_ENABLE_OPENMP OR Tasmanian_ENABLE_RECOMMENDED) AND NOT Tasmanian_ENABLE_DPCPP)
    if (Tasmanian_ENABLE_OPENMP)
        find_package(OpenMP REQUIRED) # OpenMP requested explicitly, require regardless of ENABLE_RECOMMENDED
    else()
        find_package(OpenMP)
        set(Tasmanian_ENABLE_OPENMP ${OPENMP_CXX_FOUND}) # set ENABLE_OPENMP if OpenMP_CXX_FOUND
    endif()
endif()

# multi-threading is required by the Addons module even if OpenMP is not enabled
if (NOT Tasmanian_ENABLE_OPENMP)
    find_package(Threads REQUIRED)
endif()

# check for BLAS
if (Tasmanian_ENABLE_BLAS OR Tasmanian_ENABLE_RECOMMENDED)
    include(CheckLibraryExists)
    if (NOT DEFINED BLAS_LIBRARIES) # user defined BLAS libraries are an XSDK requirement
        if (Tasmanian_ENABLE_BLAS)
            find_package(BLAS REQUIRED) # if BLAS enabled explicitly, require
            find_package(LAPACK REQUIRED)
        else()
            find_package(BLAS)
            find_package(LAPACK)
            if (BLAS_FOUND AND LAPACK_FOUND)
                set(Tasmanian_ENABLE_BLAS ON)
            else()
                set(Tasmanian_ENABLE_BLAS OFF)
            endif()
        endif()
    endif()
    set(CMAKE_REQUIRED_QUIET ON)
    foreach(_tsglib ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        check_library_exists(${_tsglib} zgelq_ "" _tsgfound)
        if (_tsgfound)
            set(Tasmanian_BLAS_HAS_ZGELQ ON)
        endif()
        unset(_tsgfound CACHE)
    endforeach()
    unset(_tsglib)
    unset(_tsgfound)
    if (Tasmanian_BLAS_HAS_ZGELQ)
        message(STATUS "Checking for LAPACK LQ factorization: found")
    else()
        message(STATUS "Checking for LAPACK LQ factorization: not found (using QR instead)")
    endif()
endif()

# Python module requires a shared library
if (Tasmanian_ENABLE_PYTHON AND NOT BUILD_SHARED_LIBS)
    message(FATAL_ERROR "BUILD_SHARED_LIBS is OFF, but shared libraries are required by the Tasmanian Python module")
endif()

# Python setup, look for python
if (Tasmanian_ENABLE_PYTHON OR (Tasmanian_ENABLE_RECOMMENDED AND BUILD_SHARED_LIBS))
    find_package(Python 3.0 COMPONENTS Interpreter)

    if (Python_FOUND)
        # manually test for NumPy, the cmake native lookup for NumPy is not reliable when using non-native cmake in Ubuntu
        # the search itself is designed to find the NumPy headers for compiling with C and not just the import numpy
        execute_process(COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/PythonImportTest.py"
                        RESULT_VARIABLE Tasmanian_python_has_numpy)
        if (NOT "${Tasmanian_python_has_numpy}" STREQUAL "0") # exit code 0 means all OK
            if (Tasmanian_ENABLE_PYTHON)
                message(FATAL_ERROR "CMake found Python, but simple 'import numpy' failed. Probably missing numpy module.")
            endif()
        else() # found python and numpy
            set(Tasmanian_ENABLE_PYTHON ON)
        endif()
    else()
        if (Tasmanian_ENABLE_PYTHON)
            message(FATAL_ERROR "CMake could not find Python")
        endif()
    endif()
endif()

# CUDA and HIP cannot be used simultaneously
if ((Tasmanian_ENABLE_CUDA AND Tasmanian_ENABLE_HIP) OR
    (Tasmanian_ENABLE_CUDA AND Tasmanian_ENABLE_DPCPP) OR
    (Tasmanian_ENABLE_DPCPP AND Tasmanian_ENABLE_HIP))
    message(FATAL_ERROR "Tasmanian can only use one GPU backend at a time, pick either CUDA, HIP or DPCPP")
endif()

# using the Tasmanian find modules requires the path
if (Tasmanian_ENABLE_DPCPP OR Tasmanian_ENABLE_MAGMA)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/")
endif()

# Tasmanian_ENABLE_CUDA support
if (Tasmanian_ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

# AMD HIP support
if (Tasmanian_ENABLE_HIP)
    if (CMAKE_VERSION VERSION_LESS 3.21)
        message(FATAL_ERROR "Tasmanian HIP/ROCm GPU backend requires CMake version 3.21 or newer")
    endif()
    if (ROCM_PATH)
        list (APPEND CMAKE_PREFIX_PATH ${ROCM_PATH}/hip ${ROCM_PATH})
    elseif (DEFINED ENV{ROCM_PATH})
        list (APPEND CMAKE_PREFIX_PATH ENV{ROCM_PATH}/hip ENV{ROCM_PATH})
    elseif (IS_DIRECTORY /opt/rocm)
        list (APPEND CMAKE_PREFIX_PATH /opt/rocm/hip /opt/rocm)
    endif()
    if (IS_DIRECTORY Tasmanian_ENABLE_ROCM)
        list (APPEND CMAKE_PREFIX_PATH ${Tasmanian_ENABLE_ROCM})
    endif()
    enable_language(HIP)
    set(Tasmanian_rocm_dependencies "rocblas;rocsparse;rocsolver")
    if (Tasmanian_ENABLE_MAGMA)
        list(APPEND Tasmanian_rocm_dependencies "hipblas;hipsparse")
    endif()
    foreach(_tsg_roclib ${Tasmanian_rocm_dependencies})
        find_package(${_tsg_roclib} REQUIRED)
    endforeach()
endif()

# Intel DPC++ support
if (Tasmanian_ENABLE_DPCPP)
    find_package(TasmanianDpcpp REQUIRED)

    if (NOT Tasmanian_ENABLE_BLAS)
        message(FATAL_ERROR "Tasmanian DPC++ capabilities require Tasmanian_ENABLE_BLAS which must be set to MKL")
    endif()
endif()

# check for MAGMA
if (Tasmanian_ENABLE_MAGMA)
    if (NOT Tasmanian_ENABLE_CUDA AND NOT Tasmanian_ENABLE_HIP)
        message(FATAL_ERROR "Tasmanian uses the GPU capabilities of MAGMA thus either CUDA or HIP must be enabled too.")
    endif()

    find_package(TasmanianMagma REQUIRED)
endif()

# check for MPI
if (Tasmanian_ENABLE_MPI)
    find_package(MPI REQUIRED)
endif()

# check if building with Python scikit-build, i.e., pip install
if (SKBUILD)
    # scikit build compiles and install in one place, then moves the files to a new location
    # Tasmanian needs the final install path so scripts can find the libraries and
    # the libraries can find each-other with rpath without LD_LIBRARY_PATH
    set(Tasmanian_final_install_path "${Tasmanian_python_pip_final}")
else()
    set(Tasmanian_final_install_path "${CMAKE_INSTALL_PREFIX}")
endif()


########################################################################
# Recommended compiler flags: Intel hasn't been tested in a while
########################################################################
if (Tasmanian_ENABLE_RECOMMENDED)
    if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel"))
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3") # there is no point in making a slow debug
        if (Tasmanian_ENABLE_FORTRAN)
            set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O3")
        endif()
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        set(CMAKE_CXX_FLAGS_RELEASE "/MD /Ox /DNDEBUG ${OpenMP_CXX_FLAGS}") # cmake overwrites flags a lot, careful how you set those
        set(CMAKE_CXX_FLAGS_DEBUG "/MDd /Ox ${OpenMP_CXX_FLAGS}")
    endif()
endif()


########################################################################
# Check for the git commit hash, if using a git repo
########################################################################
if ("${Tasmanian_version_comment}" STREQUAL "")
    set(Tasmanian_git_hash "Release ${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}")
elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git") # this is a git repo
    find_package(Git)
    # do not set the hash if git is missing or
    # if we are generating files for simple GNU Make compatibility
    if (Git_FOUND)
        execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:%H -n 1
                        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                        OUTPUT_VARIABLE   Tasmanian_git_hash)
    endif()
endif()
if (NOT Tasmanian_git_hash)
   # if not in a git repo, or there is no git executable, or if generating GNU Make files
   set(Tasmanian_git_hash "Tasmanian git hash is not available here")
endif()


########################################################################
# Extra directories:
# rarely but sometimes find_package() fails to recognize all dependencies
# the extra variables here allow the user to circumvent the problem by
# including additional directories, also check SparseGrids/CMakeLists.txt
# comment about Tasmanian_EXTRA_LIBRARIES and Tasmanian_EXTRA_INCLUDE_DIRS
########################################################################
if (DEFINED Tasmanian_EXTRA_LINK_DIRS)
    link_directories(${Tasmanian_EXTRA_LINK_DIRS}) # cannot be done per-target
endif()


########################################################################
# Report build flags based on the compiler and options
########################################################################
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(Tasmanian_cxx_flags "${CMAKE_CXX_FLAGS}")
else()
    set(Tasmanian_cxx_flags "${CMAKE_BUILD_TYPE}, ${CMAKE_CXX_FLAGS}")
endif()
