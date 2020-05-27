########################################################################
# FindTasmanianCudaMathLibs.cmake module
########################################################################
#
#
#

set(Tasmanian_ROCM_ROOT "${ROCM_ROOT}" CACHE PATH "The root folder for the Rocm framework installation")
list(APPEND CMAKE_PREFIX_PATH "${Tasmanian_ROCM_ROOT}")

foreach(_tsg_roclib hip rocblas rocsparse)
    find_package(${_tsg_roclib} REQUIRED)
endforeach()

get_filename_component(Tasmanian_hiproot ${HIP_HIPCC_EXECUTABLE} DIRECTORY) # convert <path>/bin/hipcc to <path>/bin
get_filename_component(Tasmanian_hiproot ${Tasmanian_hiproot} DIRECTORY)   # convert <path>/bin to <path>

Tasmanian_find_libraries(REQUIRED hip_hcc
                         OPTIONAL
                         PREFIX ${Tasmanian_hiproot}
                         LIST hiplibs)
Tasmanian_find_rpath(LIBRARIES ${Tasmanian_hiplibs} LIST hip_rpath)

foreach(_tsg_roclib rocblas rocsparse)
    get_filename_component(Tasmanian_roclib_root ${${_tsg_roclib}_INCLUDE_DIR} DIRECTORY)
    list(APPEND Tasmanian_hiplibs roc::${_tsg_roclib})
    list(APPEND Tasmanian_hip_rpath "${Tasmanian_roclib_root}/lib")
endforeach()
unset(_tsg_roclib)

find_package_handle_standard_args(TasmanianRocm DEFAULT_MSG Tasmanian_hiplibs)

if (Tasmanian_ENABLE_OPENMP)
    set(OpenMP_CXX_FLAGS "-fopenmp=libiomp5")
    Tasmanian_find_libraries(REQUIRED iomp5
                             OPTIONAL
                             PREFIX /usr/lib/x86_64-linux-gnu
                             LIST hipomp)
    if (Tasmanian_hipomp)
        set(OpenMP_CXX_LIBRARIES "${Tasmanian_hipomp}")
        set(Tasmanian_ENABLE_OPENMP ON)
    else()
        if (Tasmanian_ENABLE_OPENMP)
            message(FATAL_ERROR "Cannot find libiomp5 which is needed by HIP/Clang to enable OpenMP")
        endif()
    endif()
endif()
