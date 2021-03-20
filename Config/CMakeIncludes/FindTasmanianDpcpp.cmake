########################################################################
# FindTasmanianDpcpp.cmake module
########################################################################
#
#
#

Tasmanian_compiler_type(COMPILER ${CMAKE_CXX_COMPILER} TYPE "DPC++" RESULT Tasmanian_dpcpp_compiler)

get_filename_component(Tasmanian_dpcpproot ${CMAKE_CXX_COMPILER} DIRECTORY)  # convert <path>/bin/dpcpp to <path>/bin
get_filename_component(Tasmanian_dpcpproot ${Tasmanian_dpcpproot} DIRECTORY)  # convert <path>/bin to <path>

if (NOT Tasmanian_dpcpp_compiler)
    message(FATAL_ERROR "Tasmanian_ENABLE_DPCPP requires that the CMAKE_CXX_COMPILER is set to the Intel dpcpp compiler.")
endif()

set(Tasmanian_MKL_SYCL_ROOT "${MKL_SYCL_ROOT}" CACHE PATH "The root folder for the Intel OneAPI framework installation")

Tasmanian_find_libraries(REQUIRED libmkl_sycl.a
                         OPTIONAL libmkl_intel_lp64.a libmkl_intel_thread.a libmkl_core.a
                         PREFIX ${Tasmanian_MKL_SYCL_ROOT}
                         LIST mklsycl)
Tasmanian_find_libraries(REQUIRED OpenCL
                         PREFIX ${Tasmanian_dpcpproot}/lib/
                         LIST mklsycl)
Tasmanian_find_libraries(REQUIRED iomp5
                         PREFIX ${Tasmanian_dpcpproot}/compiler/lib/intel64/
                         LIST syclomp)
list(APPEND Tasmanian_mklsycl ${Tasmanian_syclomp})

if (Tasmanian_ENABLE_OPENMP OR Tasmanian_ENABLE_RECOMMENDED)
    set(OpenMP_CXX_FLAGS "-openmp;-fopenmp=libiomp5")
    set(OpenMP_CXX_LIBRARY "${Tasmanian_syclomp};-pthread")
    set(OpenMP_CXX_LIB_NAMES "CXX")
    set(Tasmanian_ENABLE_OPENMP ON)
endif()

find_package_handle_standard_args(TasmanianDpcpp DEFAULT_MSG Tasmanian_mklsycl)

foreach(_tsg_intellib ${Tasmanian_mklsycl})
    get_filename_component(_tsg_mklsycl_rpath ${_tsg_intellib} DIRECTORY)
    list(APPEND Tasmanian_mklsycl_rpath "${_tsg_mklsycl_rpath}")
endforeach()
unset(_tsg_mklsycl_rpath)
unset(_tsg_intellib)
