########################################################################
# Adds a target that includes all Tasmanian dependencies
########################################################################

add_library(Tasmanian_dependencies INTERFACE)

if (Tasmanian_ENABLE_OPENMP)
    target_link_libraries(Tasmanian_dependencies INTERFACE ${OpenMP_CXX_LIBRARIES})
else()
    target_link_libraries(Tasmanian_dependencies INTERFACE ${CMAKE_THREAD_LIBS_INIT})
endif()

if (Tasmanian_ENABLE_BLAS)
    target_link_libraries(Tasmanian_dependencies INTERFACE ${BLAS_LIBRARIES})
    target_link_libraries(Tasmanian_dependencies INTERFACE ${LAPACK_LIBRARIES})
endif()

if (Tasmanian_ENABLE_CUDA)
    target_link_libraries(Tasmanian_dependencies INTERFACE ${Tasmanian_cudamathlibs})
endif()

if (Tasmanian_ENABLE_MAGMA)
    if (BUILD_SHARED_LIBS)
        target_link_libraries(Tasmanian_dependencies INTERFACE ${Tasmanian_MAGMA_SHARED_LIBRARIES})
    else()
        target_link_libraries(Tasmanian_dependencies INTERFACE ${Tasmanian_MAGMA_LIBRARIES})
    endif()

    target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${Tasmanian_MAGMA_INCLUDE_DIRS}/>)
    target_include_directories(Tasmanian_dependencies INTERFACE $<INSTALL_INTERFACE:${Tasmanian_MAGMA_INCLUDE_DIRS}/>)
endif()

target_include_directories(Tasmanian_dependencies INTERFACE $<INSTALL_INTERFACE:${Tasmanian_final_install_path}/include>)

target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Config/>)
target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/InterfaceTPL/>)
target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/configured/>)

# Tasmanian_EXTRA_LIBRARIES gives the user an option to force extra dependencies,
# for example, on some systems (e.g., OLCF) find_package(BLAS) fails to
# recognize that libacml_mp requires libgomp, so the build fails with either clang or ENABLE_OPENMP=OFF
# -D Tasmanian_EXTRA_LIBRARIES=/path/to/libgomp.so circumvents the issue
# NOTE: adding Tasmanian_EXTRA_LIBRARIES to SparseGrids will propagate to all other targets
# same holds for Tasmanian_EXTRA_INCLUDE_DIRS
target_link_libraries(Tasmanian_dependencies ${Tasmanian_EXTRA_LIBRARIES})

foreach(_tsg_include ${Tasmanian_EXTRA_INCLUDE_DIRS})
    target_include_directories(Tasmanian_dependencies PUBLIC $<BUILD_INTERFACE:${_tsg_include}>)
    target_include_directories(Tasmanian_dependencies PUBLIC $<INSTALL_INTERFACE:${_tsg_include}>)
endforeach()

install(TARGETS Tasmanian_dependencies EXPORT "${Tasmanian_export_name}")

########################################################################
# Create an RPATH list for all dependencies
########################################################################
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
