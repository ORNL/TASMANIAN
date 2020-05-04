########################################################################
# sanity check and xSDK compatibility
########################################################################
if ((NOT DEFINED CMAKE_BUILD_TYPE) OR (NOT CMAKE_BUILD_TYPE))
    set(CMAKE_BUILD_TYPE Debug)
endif()

# XSDK mode:
#   - Never overwrite user preferences, Tasmanian_ENABLE_RECOMMENDED OFF
#   - All Tasmanian_ENABLE options are disabled (by default)
#   - Options are enabled with XSDK switches, e.g., XSDK_ENABLE_FORTRAN
if (USE_XSDK_DEFAULTS)
    set(Tasmanian_ENABLE_RECOMMENDED OFF)
    set(Tasmanian_ENABLE_OPENMP      OFF)
    set(Tasmanian_ENABLE_BLAS        OFF)
    set(Tasmanian_ENABLE_MPI         OFF)
    set(Tasmanian_ENABLE_PYTHON      OFF)
    set(Tasmanian_ENABLE_FORTRAN     OFF)
    set(Tasmanian_ENABLE_SWIG        OFF)
    set(Tasmanian_ENABLE_CUDA        OFF)
    set(Tasmanian_ENABLE_MAGMA       OFF)
    if (DEFINED XSDK_ENABLE_OPENMP)
        set(Tasmanian_ENABLE_OPENMP ${XSDK_ENABLE_OPENMP})
    endif()
    if (DEFINED TPL_ENABLE_MPI)
        set(Tasmanian_ENABLE_MPI ${TPL_ENABLE_MPI})
    endif()
    if (DEFINED TPL_ENABLE_BLAS)
        set(Tasmanian_ENABLE_BLAS ${TPL_ENABLE_BLAS})
    endif()
    if (DEFINED TPL_ENABLE_MAGMA)
        set(Tasmanian_ENABLE_MAGMA ${TPL_ENABLE_MAGMA})
        # don't really like the xSDK "TPL" convention, prefer to work
        # with Tasmanian_ variables to avoid confusion on who defined what
        if (TPL_MAGMA_LIBRARIES)
            set(Tasmanian_MAGMA_LIBRARIES ${TPL_MAGMA_LIBRARIES})
        endif()
        if (TPL_MAGMA_INCLUDE_DIRS)
            set(Tasmanian_MAGMA_INCLUDE_DIRS ${TPL_MAGMA_INCLUDE_DIRS})
        endif()
    endif()
    if (DEFINED XSDK_ENABLE_PYTHON)
        set(Tasmanian_ENABLE_PYTHON ${XSDK_ENABLE_PYTHON})
    endif()
    if (DEFINED XSDK_ENABLE_FORTRAN)
        set(Tasmanian_ENABLE_FORTRAN ${XSDK_ENABLE_FORTRAN})
    endif()
    if (DEFINED XSDK_ENABLE_CUDA)
        set(Tasmanian_ENABLE_CUDA ${XSDK_ENABLE_CUDA})
    endif()
endif()

if (SKBUILD) # get extra options from the ENV variables (pip-installer)
    if (NOT "$ENV{Tasmanian_ENABLE_BLAS}" STREQUAL "")
        set(Tasmanian_ENABLE_BLAS ON)
        set(BLAS_LIBRARIES   "$ENV{Tasmanian_ENABLE_BLAS}")
        set(LAPACK_LIBRARIES "$ENV{Tasmanian_ENABLE_BLAS}")
    endif()
    if (NOT "$ENV{Tasmanian_ENABLE_CUDA}" STREQUAL "")
        set(Tasmanian_ENABLE_CUDA ON)
        set(CMAKE_CUDA_COMPILER "$ENV{Tasmanian_ENABLE_CUDA}")
    endif()
    if (NOT "$ENV{Tasmanian_ENABLE_MAGMA}" STREQUAL "")
        set(Tasmanian_ENABLE_MAGMA ON)
        set(MAGMA_ROOT_DIR "$ENV{Tasmanian_ENABLE_MAGMA}")
    endif()
    if (NOT "$ENV{Tasmanian_ENABLE_MPI}" STREQUAL "")
        set(Tasmanian_ENABLE_MPI ON)
        set(MPI_CXX_COMPILER "$ENV{Tasmanian_ENABLE_MPI}")
    endif()
    if (NOT "$ENV{Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
        set(Tasmanian_MATLAB_WORK_FOLDER "$ENV{Tasmanian_MATLAB_WORK_FOLDER}")
    endif()
endif()

# when choosing shared/static libraries, pick the first mode that applies
# - BUILD_SHARED_LIBS=OFF: build only static libs regardless of USE_XSDK_DEFAULTS
# - BUILD_SHARED_LIBS=ON or USE_XSDK_DEFAULTS=ON: build only shared libs
# - BUILD_SHARED_LIBS=Undefined and USE_XSDK_DEFAULTS=OFF: build both types
if ((NOT "${BUILD_SHARED_LIBS}" STREQUAL "") AND (NOT BUILD_SHARED_LIBS)) # BUILD_SHARED_LIBS is defined and not an empty string
    list(APPEND Tasmanian_libs_type "static")
    set(Tasmanian_lib_default "static") # build static libs and default to static
elseif (BUILD_SHARED_LIBS OR USE_XSDK_DEFAULTS)
    list(APPEND Tasmanian_libs_type "shared")
    set(Tasmanian_lib_default "shared") # build shared libs and default to shared
else()
    list(APPEND Tasmanian_libs_type "static" "shared")
    set(Tasmanian_lib_default "static") # build both types of libs and default to static
endif()

# check for Fortran, note that enable_language always gives FATAL_ERROR if the compiler is missing
if (Tasmanian_ENABLE_FORTRAN)
    enable_language(Fortran)
endif()

# swig requires Fortran and cannot handle both types of libs
if (Tasmanian_ENABLE_SWIG)
    if (NOT Tasmanian_ENABLE_FORTRAN)
        message(FATAL_ERROR "Tasmanian_ENABLE_SWIG=ON requires Tasmanian_ENABLE_FORTRAN=ON")
    endif()
    if ("${BUILD_SHARED_LIBS}" STREQUAL "")
        message(FATAL_ERROR "Tasmanian_ENABLE_SWIG=ON requires BUILD_SHARED_LIBS to be defined (ON or OFF)")
    endif()
endif()

# OpenMP setup
if (Tasmanian_ENABLE_OPENMP OR Tasmanian_ENABLE_RECOMMENDED)
    if (Tasmanian_ENABLE_OPENMP)
        find_package(OpenMP REQUIRED) # OpenMP requested explicitly, require regardless of ENABLE_RECOMMENDED
    else()
        find_package(OpenMP)
        set(Tasmanian_ENABLE_OPENMP ${OPENMP_CXX_FOUND}) # set ENABLE_OPENMP if OpenMP_CXX_FOUND
    endif()
endif()

# fallback threads library if OpenMP is disabled, needed for Addons
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
    foreach(_tsglib ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        check_library_exists(${_tsglib} zgelq_ "" _tsgfound)
        if (_tsgfound)
            set(Tasmanian_BLAS_HAS_ZGELQ ON)
        endif()
    endforeach()
    unset(_tsglib)
    unset(_tsgfound)
endif()

# Python module requires a shared library
if (Tasmanian_ENABLE_PYTHON AND (NOT "shared" IN_LIST Tasmanian_libs_type))
    message(FATAL_ERROR "BUILD_SHARED_LIBS is OFF, but shared libraries are required by the Tasmanian Python module")
endif()

# Python setup, look for python
if (Tasmanian_ENABLE_PYTHON OR (Tasmanian_ENABLE_RECOMMENDED AND ("shared" IN_LIST Tasmanian_libs_type)))
    find_package(PythonInterp)

    if (PYTHONINTERP_FOUND)
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/PythonImportTest.py"
                        RESULT_VARIABLE Tasmanian_python_has_numpy OUTPUT_QUIET)
    else()
        set(Tasmanian_python_has_numpy "1") # 1 is error code if the script above fails
    endif()

    if (NOT "${Tasmanian_python_has_numpy}" STREQUAL "0")
        if (NOT Tasmanian_ENABLE_PYTHON) # means we are using RECOMMENDED only
            message(STATUS "Tasmanian could not find Python with numpy and ctypes modules, the python interface will not be installed and some tests will be omitted")
        else()
            message(FATAL_ERROR "-D Tasmanian_ENABLE_PYTHON is ON, but either find_package(PythonInterp) failed python executable could not 'import numpy, ctypes'\nuse -D PYTHON_EXECUTABLE:PATH to specify suitable python interpreter")
        endif()
    else()
        set(Tasmanian_ENABLE_PYTHON ON) # just in case we are using RECOMMENDED only
    endif()
endif()

# CUDA and HIP cannot be used simultaneously
if (Tasmanian_ENABLE_CUDA AND Tasmanian_ENABLE_HIP)
    message(FATAL_ERROR "Tasmanian can only use one GPU backend at a time, pick either HIP or CUDA, not both")
endif()

# Tasmanian_ENABLE_CUDA support
if (Tasmanian_ENABLE_CUDA)
    enable_language(CUDA)

    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/")
    find_package(TasmanianCudaMathLibs REQUIRED)
endif()

# AMD HIP support
if (Tasmanian_ENABLE_HIP)
    message(WARNING "The AMD HIP backend is currently a placeholder, all operations will fallback to the CPU, do NOT use in production")
endif()

# check for MAGMA
if (Tasmanian_ENABLE_MAGMA)
    if (NOT Tasmanian_ENABLE_CUDA)
        message(FATAL_ERROR "Currently Tasmanian can use only CUDA related capability from MAGMA, hence Tasmanian_ENABLE_CUDA must be set ON")
    endif()

    find_package(TasmanianMAGMA)

    if (Tasmanian_MAGMA_FOUND)
        message(STATUS "Tasmanian will use UTK MAGMA libraries (static link): ${Tasmanian_MAGMA_LIBRARIES}")
        if ("shared" IN_LIST Tasmanian_libs_type) # requesting shared libraries for Tasmanian
            message(STATUS "Tasmanian will use UTK MAGMA libraries (shared link): ${Tasmanian_MAGMA_SHARED_LIBRARIES}")
            if (NOT Tasmanian_MAGMA_SHARED_FOUND)
                message(WARNING "Setting up build with shared libraries for Tasmanian but the UTK MAGMA appears to provide static libraries only \n attempting to link anyway, but this is likely to fail\nif encountering a problem call cmake again with -D BUILD_SHARED_LIBS=OFF")
            endif()
        endif()
        message(STATUS "Tasmanian will use UTK MAGMA include: ${Tasmanian_MAGMA_INCLUDE_DIRS}")
    else()
        message(FATAL_ERROR "Tasmanian_ENABLE_MAGMA is ON, but find_package(TasmanianMAGMA) failed\n please provide valid Tasmanian_MAGMA_ROOT:PATH or Tasmanian_MAGMA_LIBRARIES with Tasmanian_MAGMA_INCLUDE_DIRS:PATH")
    endif()
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
# Needed for TasmanianConfig.cmake and TasmanianConfigVersion.cmake
# enables the use of "find_package(Tasmanian <version>)"
########################################################################
include(CMakePackageConfigHelpers)


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
