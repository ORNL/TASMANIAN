########################################################################
# FindTasmanianDpcpp.cmake module
########################################################################
#
# DPC++ and SYCL compatibility
# 1. Check if the compiler --version and --help return words like "sycl" and "DPC++", i.e., didn't accidentally set the compiler to g++
# 2. Use -DMKLROOT or environment MKLROOT to find both MKL and OneMKL
# 3. Manually set the OpenMP flags, since CMake doesn't know the correct flags for DPC++
# 4. Print status message
#
# Defines Tasmanian_mklsycl variables with all the relevant variables for the dependencies
# Note: the -fsycl flag is applied privately to lib-sparsegrid
#

Tasmanian_compiler_type(COMPILER ${CMAKE_CXX_COMPILER} SWITCH "--version" TYPE "DPC++" RESULT Tasmanian_dpcpp_compiler)
Tasmanian_compiler_type(COMPILER ${CMAKE_CXX_COMPILER} SWITCH "--help"    TYPE "sycl"  RESULT Tasmanian_sycl_compiler)

get_filename_component(Tasmanian_dpcpproot ${CMAKE_CXX_COMPILER} DIRECTORY)  # convert <path>/bin/dpcpp to <path>/bin
get_filename_component(Tasmanian_dpcpproot ${Tasmanian_dpcpproot} DIRECTORY)  # convert <path>/bin to <path>

if (NOT Tasmanian_dpcpp_compiler AND NOT Tasmanian_sycl_compiler)
    message(FATAL_ERROR "Tasmanian_ENABLE_DPCPP requires that the CMAKE_CXX_COMPILER is set to the Intel dpcpp compiler or a compatible sycl compiler.")
endif()

if (MKLROOT)
    set(Tasmanian_MKL_SYCL_ROOT "${MKLROOT}"    CACHE PATH "The root folder for the Intel oneMKL installation")
else()
    set(Tasmanian_MKL_SYCL_ROOT "$ENV{MKLROOT}" CACHE PATH "The root folder for the Intel oneMKL installation")
endif()

Tasmanian_find_libraries(REQUIRED mkl_sycl
                         OPTIONAL mkl_intel_lp64 mkl_intel_thread mkl_core
                         PREFIX ${Tasmanian_MKL_SYCL_ROOT}
                         LIST mklsycl)

if (Tasmanian_ENABLE_OPENMP OR Tasmanian_ENABLE_RECOMMENDED)
    set(OpenMP_CXX_FLAGS "-qopenmp")
    set(OpenMP_CXX_LIBRARY "${Tasmanian_syclomp};-pthread")
    set(OpenMP_CXX_LIB_NAMES "CXX")
    set(Tasmanian_ENABLE_OPENMP ON)
endif()

foreach(_tsg_intellib ${Tasmanian_mklsycl})
    if (NOT _tsg_intellib)
        message(FATAL_ERROR "Missing a oneAPI/oneMKL library: ${_tsg_intellib}, try setting MKLROOT to the oneMKL installation path.")
    endif()
    get_filename_component(_tsg_mklsycl_rpath ${_tsg_intellib} DIRECTORY)
    list(APPEND Tasmanian_mklsycl_rpath "${_tsg_mklsycl_rpath}")
endforeach()
unset(_tsg_mklsycl_rpath)
unset(_tsg_intellib)

find_package_handle_standard_args(TasmanianDpcpp DEFAULT_MSG Tasmanian_mklsycl)
