########################################################################
# FindTasmanianCudaMathLibs.cmake module
########################################################################
#
# requires: must be called AFTER enable_language(CUDA)
#
# Uses the CUDA compiler location to deduce the location of
# the math libraries: cuBlas cuSparse cuSolver
# The CUDA Run-time libraries are also located for fallback purposes.

# the macro is called as Tasmanian_find_cuda_libraries(NAMES foo1 foo2 REQUIRED foo1 PREFIX bar LIST saloon)
# which will search for foo1 and foo2 using bar as hints, then put the search result in
# variables called Tasmanian_foo1 and Tasmanian_foo2
# finally, all results will be appended into a list called Tasmanian_saloon
# the REQUIRED libraries will be kept as variables, i.e., Tasmanian_foo1 to be used
# with find_package_handle_standard_args()
macro(Tasmanian_find_cuda_libraries)
# it is tempting to say NO_DEFAULT_PATH here, but newer versions of CUDA
# install cuBlas in /usr/lib while nvcc and the other libs sit in /usr/local/cuda
    cmake_parse_arguments(Tasmanian_cudalibs "" "PREFIX;LIST" "REQUIRED;NAMES" ${ARGN} )
    foreach(_lib ${Tasmanian_cudalibs_REQUIRED})
        list(APPEND Tasmanian_required_cudalibs "Tasmanian_${_lib}")
    endforeach()
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
        list(FIND Tasmanian_required_cudalibs "Tasmanian_${_lib}" Tasmanian_is_required)
        if (${Tasmanian_is_required} EQUAL -1) # is_required corresponds to the position in the list, -1 means missing
            unset(Tasmanian_${_lib} CACHE)
        endif()
    endforeach()
    if (Tasmanian_cudalibs_REQUIRED)
        find_package_handle_standard_args(TasmanianCudaMathLibs DEFAULT_MSG ${Tasmanian_required_cudalibs})
        foreach(_lib ${Tasmanian_required_cudalibs})
            get_filename_component(Tasmanian_libdir ${${_lib}} DIRECTORY)
            list(APPEND Tasmanian_cuda_rpath "${Tasmanian_libdir}")
        endforeach()
    endif()
    unset(_lib)
    unset(Tasmanian_is_required)
    unset(Tasmanian_required_cudalibs)
endmacro()

include(FindPackageHandleStandardArgs)

get_filename_component(Tasmanian_nvccroot ${CMAKE_CUDA_COMPILER} DIRECTORY) # convert <path>/bin/nvcc to <path>/bin
get_filename_component(Tasmanian_nvccroot ${Tasmanian_nvccroot} DIRECTORY)  # convert <path>/bin to <path>

if (NOT CMAKE_LIBRARY_ARCHITECTURE)
    set(CMAKE_LIBRARY_ARCHITECTURE "x64") # sometimes missing under Windows
endif()

Tasmanian_find_cuda_libraries(NAMES cusolver cusparse culibos cublas cublas_device cublasLt
                              REQUIRED cublas cusparse cusolver
                              PREFIX ${Tasmanian_nvccroot}
                              LIST cudamathlibs)

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
