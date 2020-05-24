########################################################################
# FindTasmanianCudaMathLibs.cmake module
########################################################################
#
#
#

set(Tasmanian_ROCM_ROOT "${ROCM_ROOT}" CACHE PATH "The root folder for the Rocm framework installation")
list(APPEND CMAKE_PREFIX_PATH "${Tasmanian_ROCM_ROOT}")

find_package(hip REQUIRED)
find_package(rocblas REQUIRED)

get_filename_component(Tasmanian_hiproot ${HIP_HIPCC_EXECUTABLE} DIRECTORY) # convert <path>/bin/hipcc to <path>/bin
get_filename_component(Tasmanian_hiproot ${Tasmanian_hiproot} DIRECTORY)   # convert <path>/bin to <path>

get_filename_component(Tasmanian_rocblas_root ${rocblas_INCLUDE_DIR} DIRECTORY)

list(APPEND Tasmanian_hip_rpath "${Tasmanian_rocblas_root}/lib")
