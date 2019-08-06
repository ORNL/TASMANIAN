########################################################################
# FindTasmanianCudaMathLibs.cmake module
########################################################################
#
# requires: must be called AFTER enable_language(CUDA)
#
# Uses the CUDA compiler location to deduce the location of
# the cuBlas and cuSparse Math libraries.
# The CUDA Run-time libraries are also located for fallback purposes.

# the macro is called as Tasmanian_find_cuda_libraries(NAMES foo1 foo2 PREFIX bar LIST saloon)
# which will search for foo1 and foo2 using bar as hints, then put the search result in
# variables called Tasmanian_foo1 and Tasmanian_foo2
# finally, all results will be appended into a list called Tasmanian_saloon
macro(Tasmanian_find_cuda_libraries)
# it is tempting to say NO_DEFAULT_PATH here, but newer versions of CUDA
# install cuBlas in /usr/lib while nvcc and the other libs sit in /usr/local/cuda
    cmake_parse_arguments(Tasmanian_cudalibs "" "PREFIX;LIST" "NAMES" ${ARGN} )
    foreach(_lib ${Tasmanian_cudalibs_NAMES})
        find_library(Tasmanian_${_lib} ${_lib}
                     HINTS "${Tasmanian_cudalibs_PREFIX}"
                     HINTS "${Tasmanian_cudalibs_PREFIX}/lib/"
                     HINTS "${Tasmanian_cudalibs_PREFIX}/${CMAKE_LIBRARY_ARCHITECTURE}/lib/"
                     HINTS "${Tasmanian_cudalibs_PREFIX}/lib/${CMAKE_LIBRARY_ARCHITECTURE}"
                     HINTS "${Tasmanian_cudalibs_PREFIX}/lib64/"
                     HINTS "${Tasmanian_cudalibs_PREFIX}/${CMAKE_LIBRARY_ARCHITECTURE}/lib64/"
                     HINTS "${Tasmanian_cudalibs_PREFIX}/lib64/${CMAKE_LIBRARY_ARCHITECTURE}")

        list(APPEND Tasmanian_${Tasmanian_cudalibs_LIST} ${Tasmanian_${_lib}})
    endforeach()
endmacro()

include(FindPackageHandleStandardArgs)

get_filename_component(Tasmanian_nvccroot ${CMAKE_CUDA_COMPILER} DIRECTORY) # convert <path>/bin/nvcc to <path>/bin
get_filename_component(Tasmanian_nvccroot ${Tasmanian_nvccroot} DIRECTORY)  # convert <path>/bin to <path>

if (NOT CMAKE_LIBRARY_ARCHITECTURE)
    set(CMAKE_LIBRARY_ARCHITECTURE "x64") # sometimes missing under Windows
endif()

Tasmanian_find_cuda_libraries(NAMES cublas_static cublas_device cublasLt cublasLt_static cublas
                                    cusparse_static cusparse culibos
                              PREFIX ${Tasmanian_nvccroot}
                              LIST cudamathlibs)

find_package_handle_standard_args(TasmanianCudaMathLibs DEFAULT_MSG Tasmanian_cublas Tasmanian_cusparse)

# different versions of CUDA require different combination of cublas_device and cublas_lt libraries
# in place of building a database with detailed library list for each version of cuda,
# assume the correct libraries are present and ignore any missing ones
list(FILTER Tasmanian_cudamathlibs EXCLUDE REGEX "-NOTFOUND$")

Tasmanian_find_cuda_libraries(NAMES cuda cudart
                              PREFIX ${Tasmanian_nvccroot}
                              LIST cudaruntime)

list(FILTER Tasmanian_cudaruntime EXCLUDE REGEX "-NOTFOUND$")

#message(STATUS "Tasmanian found CUDA math libraries: ${Tasmanian_cudamathlibs}")
#message(STATUS "Tasmanian found CUDA runtime libraries: ${Tasmanian_cudaruntime}")
