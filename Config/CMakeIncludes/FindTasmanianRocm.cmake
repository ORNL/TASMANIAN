########################################################################
# FindTasmanianCudaMathLibs.cmake module
########################################################################
#
#
#

execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE Tasmanian_hip_compiler_version)
string(REGEX MATCH "^HIP" Tasmanian_haship "${Tasmanian_hip_compiler_version}")

if (NOT Tasmanian_haship)
    message(WARNING "Tasmanian_ENABLE_HIP requires that the CMAKE_CXX_COMPILER is set to the Rocm hipcc compiler.")
endif()

get_filename_component(Tasmanian_hipccroot ${CMAKE_CXX_COMPILER} DIRECTORY)  # convert <path>/bin/hipcc to <path>/bin
get_filename_component(Tasmanian_hipccroot ${Tasmanian_hipccroot} DIRECTORY)  # convert <path>/bin to <path>

set(Tasmanian_ROCM_ROOT "${ROCM_ROOT}" CACHE PATH "The root folder for the Rocm framework installation")
list(APPEND CMAKE_PREFIX_PATH "${Tasmanian_ROCM_ROOT}")
list(APPEND CMAKE_PREFIX_PATH "${Tasmanian_hipccroot}")

if (Tasmanian_ENABLE_OPENMP)
    set(Tasmanian_HIP_IOMP5_PATH "${Tasmanian_hipccroot}/llvm" CACHE PATH "The search path for libiomp5")
endif()

foreach(_tsg_roclib hip rocblas rocsparse rocsolver)
    find_package(${_tsg_roclib} REQUIRED)
endforeach()

get_filename_component(Tasmanian_hiproot ${HIP_HIPCC_EXECUTABLE} DIRECTORY) # convert <path>/bin/hipcc to <path>/bin
get_filename_component(Tasmanian_hiproot ${Tasmanian_hiproot} DIRECTORY)   # convert <path>/bin to <path>

foreach(_tsg_roclib rocblas rocsparse rocsolver)
    get_filename_component(Tasmanian_roclib_root ${${_tsg_roclib}_INCLUDE_DIR} DIRECTORY)
    list(APPEND Tasmanian_hiplibs roc::${_tsg_roclib})
    list(APPEND Tasmanian_hip_rpath "${Tasmanian_roclib_root}/lib")
endforeach()
unset(_tsg_roclib)

find_package_handle_standard_args(TasmanianRocm DEFAULT_MSG Tasmanian_hiplibs)

if (Tasmanian_ENABLE_OPENMP OR Tasmanian_ENABLE_RECOMMENDED)
    set(OpenMP_CXX_FLAGS "-fopenmp=libiomp5")
    Tasmanian_find_libraries(REQUIRED iomp5 omp
                             OPTIONAL
                             PREFIX ${Tasmanian_HIP_IOMP5_PATH}
                             LIST hipomp)
    if (Tasmanian_hipomp)
        set(OpenMP_CXX_LIBRARIES "${Tasmanian_hipomp}")
        set(Tasmanian_ENABLE_OPENMP ON)
        foreach(_tsg_roclib ${Tasmanian_hipomp})
            get_filename_component(Tasmanian_omp_root ${_tsg_roclib} DIRECTORY)
            list(APPEND Tasmanian_hip_rpath "${Tasmanian_omp_root}")
            list(APPEND Tasmanian_hipomp_rpath "${Tasmanian_omp_root}")
        endforeach()
        list(REMOVE_DUPLICATES Tasmanian_hipomp_rpath)
        unset(_tsg_roclib)
        unset(Tasmanian_omp_root)
    else()
        if (Tasmanian_ENABLE_OPENMP)
            message(FATAL_ERROR "Cannot find libiomp5 which is needed by HIP/Clang to enable OpenMP")
        endif()
    endif()
endif()
