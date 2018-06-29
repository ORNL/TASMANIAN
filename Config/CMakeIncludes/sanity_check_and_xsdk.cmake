########################################################################
# sanity check and xSDK compatibility
########################################################################
if ((NOT DEFINED CMAKE_BUILD_TYPE) OR (NOT CMAKE_BUILD_TYPE))
    set(CMAKE_BUILD_TYPE Debug)
endif()

# XSDK mode:
#   - Assume Tasmanian_STRICT_OPTIONS (never overwrite user preferences)
#   - All Tasmanian_ENABLE options are disbaled (by default)
#   - Options are enabled with XSDK switches, e.g., XSDK_ENABLE_FORTRAN
if (USE_XSDK_DEFAULTS)
    set(Tasmanian_STRICT_OPTIONS ON)
    set(Tasmanian_ENABLE_OPENMP  OFF)
    set(Tasmanian_ENABLE_BLAS    OFF)
    set(Tasmanian_ENABLE_MPI     ON)
    set(Tasmanian_ENABLE_MATLAB  OFF)
    set(Tasmanian_ENABLE_PYTHON  OFF)
    set(Tasmanian_ENABLE_FORTRAN OFF)
    set(Tasmanian_ENABLE_CUBLAS  OFF)
    set(Tasmanian_ENABLE_CUDA    OFF)
    set(Tasmanian_SHARED_LIBRARY ON) # xSDK defaults to shared
    set(Tasmanian_STATIC_LIBRARY OFF)
    if (DEFINED XSDK_ENABLE_OPENMP)
        set(Tasmanian_ENABLE_OPENMP ${XSDK_ENABLE_OPENMP})
    endif()
    if (DEFINED TPL_ENABLE_BLAS)
        set(Tasmanian_ENABLE_BLAS ${TPL_ENABLE_BLAS})
    endif()
    if (DEFINED TPL_ENABLE_MAGMA)
        set(Tasmanian_ENABLE_MAGMA TPL_ENABLE_MAGMA)
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
    if (DEFINED XSDK_ENABLE_MATLAB)
        set(Tasmanian_ENABLE_MATLAB ${XSDK_ENABLE_MATLAB})
    endif()
    if (DEFINED XSDK_ENABLE_CUDA)
        set(Tasmanian_ENABLE_CUBLAS ${XSDK_ENABLE_CUDA})
        set(Tasmanian_ENABLE_CUDA   ${XSDK_ENABLE_CUDA})
    else()
        set(Tasmanian_ENABLE_CUBLAS OFF)
        set(Tasmanian_ENABLE_CUDA   OFF)
    endif()
endif()

if (DEFINED BUILD_SHARED_LIBS) # respect BUILD_SHARED_LIBS
    if (BUILD_SHARED_LIBS)
        set(Tasmanian_SHARED_LIBRARY ON)
        set(Tasmanian_STATIC_LIBRARY OFF)
    else()
        set(Tasmanian_SHARED_LIBRARY OFF)
        set(Tasmanian_STATIC_LIBRARY ON)
    endif()
endif()

# OpenMP setup
if (Tasmanian_ENABLE_OPENMP)
    find_package(OpenMP)

    if ((NOT OPENMP_FOUND) AND (NOT OPENMP_CXX_FOUND)) # OPENMP_FOUND is used prior to cmake 3.9
        if (Tasmanian_STRICT_OPTIONS)
            message(FATAL_ERROR "-D Tasmanian_ENABLE_OPENMP is ON, but find_package(OpenMP) failed")
        else()
            message("-D Tasmanian_ENABLE_OPENMP is ON, but find_package(OpenMP) failed\noverwritting option: -D Tasmanian_ENABLE_OPENMP:BOOL=OFF")
            set(Tasmanian_ENABLE_OPENMP OFF)
        endif()
    else()
        if (NOT DEFINED OpenMP_CXX_LIBRARIES)
            # on older versions of cmake use global flags, see comment in SparseGrids/CMakeLists.txt
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif()
    endif()
endif()

# Python setup, look for python
if (Tasmanian_ENABLE_PYTHON)
    find_package(PythonInterp)
    if (NOT PYTHONINTERP_FOUND)
        if (Tasmanian_STRICT_OPTIONS)
            message(FATAL_ERROR "-D Tasmanian_ENABLE_PYTHON is ON, but find_package(PythonInterp) failed\nuse -D PYTHON_EXECUTABLE:PATH to specify valid python interpreter")
        else()
            message("-D Tasmanian_ENABLE_PYTHON is ON, but find_package(PythonInterp) failed\nuse -D PYTHON_EXECUTABLE:PATH to specify python interpreter\noverwritting option: -D Tasmanian_ENABLE_PYTHON:BOOL=OFF")
            set(Tasmanian_ENABLE_PYTHON OFF)
        endif()
    endif()
endif()

# Safeguard against silly option selection
if ((NOT Tasmanian_STATIC_LIBRARY) AND (NOT Tasmanian_SHARED_LIBRARY))
    if (Tasmanian_STRICT_OPTIONS)
        message(FATAL_ERROR "Must specifiy at least one of Tasmanian_STATIC_LIBRARY or Tasmanian_SHARED_LIBRARY")
    else()
        message(WARNING "must specifiy at least one of Tasmanian_STATIC_LIBRARY or Tasmanian_SHARED_LIBRARY\noverwritting option: -DTasmanian_SHARED_LIBRARY:BOOL=ON")
        set(Tasmanian_SHARED_LIBRARY ON)
    endif()
endif()

# Python module requires a shared library
if (Tasmanian_ENABLE_PYTHON AND (NOT Tasmanian_SHARED_LIBRARY))
    if (Tasmanian_STRICT_OPTIONS)
        message(FATAL_ERROR "Tasmanian_SHARED_LIBRARY (or BUILD_SHARED_LIBS) has to be ON to use the Tasmanian Python interface")
    else()
        message(WARNING "Tasmanian_SHARED_LIBRARY is needed by the python interface\nto disable the shared library, you must disable python as well\noverwritting option: -DTasmanian_SHARED_LIBRARY:BOOL=ON")
        set(Tasmanian_SHARED_LIBRARY ON)
    endif()
endif()

# MATLAB interface requires Tasmanian_MATLAB_WORK_FOLDER path
if (Tasmanian_ENABLE_MATLAB AND (Tasmanian_MATLAB_WORK_FOLDER STREQUAL ""))
    message(FATAL_ERROR "Tasmanian_ENABLE_MATLAB is ON, but -DTasmanian_MATLAB_WORK_FOLDER:PATH=<path> is missing\nIn order to enable the MATLAB interface, you must specify -DTasmanian_MATLAB_WORK_FOLDER:PATH=")
endif()

# ENABLE_CUDA and ENABLE_CUBLAS support for the add_library vs cuda_add_library
# note that ENABLE_CUBLAS indicates both cublas and cusparse libraries
if (Tasmanian_ENABLE_CUDA OR (Tasmanian_ENABLE_CUBLAS AND (NOT DEFINED CUDA_CUBLAS_LIBRARIES)))
    find_package(CUDA)

    if (CUDA_FOUND)
        if(${CMAKE_VERSION} VERSION_LESS "3.7.0")
            # cmake versions prior to 3.7.0 do not respect c++11 flags and include directories defined per target
            # there is no other way to pass the "c++11" flag to CUDA or to include the condifured folder
            # cmake versions 3.7.0 and above do not need the statements below
            list(APPEND CUDA_NVCC_FLAGS "-std=c++11")
            include_directories("${CMAKE_CURRENT_BINARY_DIR}/configured/")
        endif()
    else()
        if (Tasmanian_STRICT_OPTIONS)
            message(FATAL_ERROR "Tasmanian_ENABLE_CUBLAS or Tasmanian_ENABLE_CUDA is on, but find_package(CUDA) failed")
        else()
            message(WARNING "find_package(CUDA) failed, could not find cuda")
            if (Tasmanian_ENABLE_CUBLAS)
                message(WARNING "overwritting option: -D Tasmanian_ENABLE_CUBLAS:BOOL=OFF")
                set(Tasmanian_ENABLE_CUBLAS OFF)
            endif()
            if (Tasmanian_ENABLE_CUDA)
                message(WARNING "overwritting option: -D Tasmanian_ENABLE_CUDA:BOOL=OFF")
                set(Tasmanian_ENABLE_CUDA OFF)
            endif()
        endif()
    endif()
endif()

# check for MAGMA
if (Tasmanian_ENABLE_MAGMA)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Config/CMakeIncludes/")
    find_package(TasmanianMAGMA)

    if (Tasmanian_MAGMA_FOUND)
        message(STATUS "Tasmanian will use UTK MAGMA libraries (static link): ${Tasmanian_MAGMA_LIBRARIES}")
        if (Tasmanian_SHARED_LIBRARY)
            message(STATUS "Tasmanian will use UTK MAGMA libraries (shared link): ${Tasmanian_MAGMA_SHARED_LIBRARIES}")
            if (NOT Tasmanian_MAGMA_SHARED_FOUND)
                message(WARNING "Tasmanian_SHARED_LIBRARY=ON but the UTK MAGMA appears to provide static libraries \n attempting to link anyway, which is likely to fail.")
            endif()
        endif()
        message(STATUS "Tasmanian will use UTK MAGMA include: ${Tasmanian_MAGMA_INCLUDE_DIRS}")
    else()
        if (Tasmanian_STRICT_OPTIONS)
            message(FATAL_ERROR "Tasmanian_ENABLE_MAGMA is ON, but find_package(TasmanianMAGMA) failed\n please provide valid Tasmanian_MAGMA_ROOT:PATH or Tasmanian_MAGMA_LIBRARIES with Tasmanian_MAGMA_INCLUDE_DIRS:PATH")
        else()
            message(WARNING "find_package(TasmanianMAGMA) failed, could not find UTK MAGMA\n overwritting option: -D Tasmanian_ENABLE_MAGMA:BOOL=OFF")
            set(Tasmanian_ENABLE_MAGMA OFF)
        endif()
    endif()
endif()

# check for BLAS
if (Tasmanian_ENABLE_BLAS)
    if (NOT DEFINED BLAS_LIBRARIES)
        find_package(BLAS)
        if (NOT BLAS_FOUND)
            if (Tasmanian_STRICT_OPTIONS)
                message(FATAL_ERROR "Tasmanian_ENABLE_BLAS is on, but find_package(BLAS) failed")
            else()
                message(WARNING "find_package(BLAS) failed, could not find BLAS\noverwritting option: -DTasmanian_ENABLE_BLAS:BOOL=OFF")
                set(Tasmanian_ENABLE_BLAS OFF)
            endif()
        endif()
    endif()
endif()

# check for MPI
if (Tasmanian_ENABLE_MPI)
    if (NOT DEFINED MPI_CXX_LIBRARIES)
        find_package(MPI)

        if (NOT MPI_CXX_FOUND)
            if (Tasmanian_STRICT_OPTIONS)
                message(FATAL_ERROR "Tasmanian_ENABLE_MPI is on, but find_package(MPI) failed")
            else()
                message(WARNING "find_package(MPI) failed, could not find MPI\noverwritting option: -D Tasmanian_ENABLE_MPI:BOOL=OFF")
                set(Tasmanian_ENABLE_MPI OFF)
            endif()
        endif()
    endif()
endif()

# check for Fortran, note that enable_language always gives FATAL_ERROR if the compiler is missing (no way to recover or autodisable)
if (Tasmanian_ENABLE_FORTRAN)
    enable_language(Fortran)
endif()

########################################################################
# Needed for TasmanianConfig.cmake and TasmanianConfigVersion.cmake
# enables the use of "find_package(Tasmanian <version>)"
########################################################################
include(CMakePackageConfigHelpers)


########################################################################
# Compiler specific flags: Intel hasn't been tested in a while
# Tasmanian_STRICT_OPTIONS=ON will prevent Tasmanian from setting
# compiler flags, only the user provided flags will be used
########################################################################
if (NOT Tasmanian_STRICT_OPTIONS)
    if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3") # there is no point in making a slow debug
        if (Tasmanian_ENABLE_FORTRAN)
            set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O3 -fno-f2c")
        endif()
#    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native -diag-disable 11074 -diag-disable 11076 -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic")
    elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Windows") # TODO: VC++ seems to ignore the cmake CXX flags
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox /EHsc")
    endif()
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
