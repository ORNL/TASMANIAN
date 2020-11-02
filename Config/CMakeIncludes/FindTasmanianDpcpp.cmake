########################################################################
# FindTasmanianDpcpp.cmake module
########################################################################
#
#
#

execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE Tasmanian_dpcpp_compiler)
string(FIND "${Tasmanian_dpcpp_compiler}" "DPC++" Tasmanian_dpcpp_compiler)

if (Tasmanian_dpcpp_compiler LESS 0)
    message(FATAL_ERROR "Tasmanian_ENABLE_DPCPP requires that the CMAKE_CXX_COMPILER is set to the Intel dpcpp compiler.")
endif()

set(Tasmanian_MKL_SYCL_ROOT "${MKL_SYCL_ROOT}" CACHE PATH "The root folder for the Intel OneAPI framework installation")

Tasmanian_find_libraries(REQUIRED libmkl_sycl.a
                         OPTIONAL libmkl_intel_lp64.a libmkl_intel_thread.a libmkl_core.a OpenCL iomp5
                         PREFIX ${Tasmanian_MKL_SYCL_ROOT}
                         LIST mklsycl)

find_package_handle_standard_args(TasmanianDpcpp DEFAULT_MSG Tasmanian_mklsycl)

foreach(_tsg_intellib ${Tasmanian_mklsycl})
    get_filename_component(_tsg_mklsycl_rpath ${_tsg_intellib} DIRECTORY)
    list(APPEND Tasmanian_mklsycl_rpath "${_tsg_mklsycl_rpath}")
endforeach()
unset(_tsg_mklsycl_rpath)
unset(_tsg_intellib)
