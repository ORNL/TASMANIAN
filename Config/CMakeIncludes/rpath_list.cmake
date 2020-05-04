#
# included by sanity_check_and_xsdk.cmake
# for each TPL option, adds the correct RPATH
#

list(APPEND Tasmanian_rpath "${Tasmanian_final_install_path}/lib")

if (Tasmanian_ENABLE_BLAS)
    foreach(_tsg_lib ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        get_filename_component(_tsg_libpath ${_tsg_lib} DIRECTORY)
        list(APPEND Tasmanian_rpath ${_tsg_libpath})
    endforeach()
    unset(_tsg_lib)
    unset(_tsg_libpath)
endif()

if (Tasmanian_ENABLE_CUDA)
    list(APPEND Tasmanian_rpath ${Tasmanian_cuda_rpath})
    unset(Tasmanian_cuda_rpath)
endif()

if (Tasmanian_ENABLE_MAGMA)
    list(APPEND Tasmanian_rpath "${Tasmanian_libmagma_dir}")
    unset(Tasmanian_libmagma_dir)
endif()

# message(STATUS "Tasmanian RPATH: ${Tasmanian_rpath}")
