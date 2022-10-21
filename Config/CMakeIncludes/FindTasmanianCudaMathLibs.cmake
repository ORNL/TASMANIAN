########################################################################
# FindTasmanianCudaMathLibs.cmake module
########################################################################
#
# requires: must be called AFTER enable_language(CUDA)
#
# Uses the CUDA compiler location to deduce the location of
# the math libraries: cuBlas cuSparse cuSolver
# The CUDA Run-time libraries are also located for fallback purposes.
#
# Defines variables:
#   Tasmanian_cudamathlibs: list of cublas, cusparse, cusolver, and helper libraries
#   Tasmanian_cudaruntime: list of runtime libraries
#
#   Tasmanian_<cublas, cusparse, cusolver>: path to individual libraries

get_filename_component(Tasmanian_nvccroot ${CMAKE_CUDA_COMPILER} DIRECTORY) # convert <path>/bin/nvcc to <path>/bin
get_filename_component(Tasmanian_nvccroot ${Tasmanian_nvccroot} DIRECTORY)  # convert <path>/bin to <path>

if (NOT CMAKE_LIBRARY_ARCHITECTURE)
    set(CMAKE_LIBRARY_ARCHITECTURE "x64") # sometimes missing under Windows
endif()

if (CMAKE_VERSION VERSION_LESS 3.17)
    # different versions of CUDA require different combination of cublas_device and cublas_lt libraries
    # assume that if something is missing then it is not needed, hence the OPTIONAL qualifier
    # the alternative is to build a database with CUDA versions and libraries, which is hard to maintain
    Tasmanian_find_libraries(REQUIRED cublas cusparse cusolver
                             OPTIONAL culibos cublas_device cublasLt
                             PREFIX ${Tasmanian_nvccroot}
                             LIST cudamathlibs)
    find_package_handle_standard_args(TasmanianCudaMathLibs DEFAULT_MSG Tasmanian_cublas Tasmanian_cusparse Tasmanian_cusolver)
else()
    # using the native CUDA finder
    find_package(CUDAToolkit REQUIRED)
    Tasmanian_find_libraries(REQUIRED OPTIONAL culibos cublas_device cublasLt
                             PREFIX ${Tasmanian_nvccroot}
                             LIST cudamathlibs)
    list(APPEND Tasmanian_cudamathlibs CUDA::cublas CUDA::cusparse CUDA::cusolver)
endif()

# runtime libraries are needed if the user wants to link to Tasmanian without enabling the CUDA language
Tasmanian_find_libraries(OPTIONAL cuda cudart
                         PREFIX ${Tasmanian_nvccroot}
                         LIST cudaruntime)
