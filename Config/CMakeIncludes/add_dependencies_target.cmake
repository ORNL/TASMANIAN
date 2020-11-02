########################################################################
# Adds a target that includes all Tasmanian dependencies
########################################################################

add_library(Tasmanian_dependencies INTERFACE)
target_compile_features(Tasmanian_dependencies INTERFACE cxx_std_11)
list(APPEND Tasmanian_rpath "${Tasmanian_final_install_path}/lib")

if (Tasmanian_ENABLE_OPENMP)
    target_link_libraries(Tasmanian_dependencies INTERFACE ${OpenMP_CXX_LIBRARIES})
    Tasmanian_find_rpath(LIBRARIES ${OpenMP_CXX_LIBRARIES} LIST rpath)
else()
    target_link_libraries(Tasmanian_dependencies INTERFACE ${CMAKE_THREAD_LIBS_INIT})
    Tasmanian_find_rpath(LIBRARIES ${CMAKE_THREAD_LIBS_INIT} LIST rpath)
endif()

if (Tasmanian_ENABLE_DPCPP)  # must come before BLAS to pick MKL libraries
    target_link_libraries(Tasmanian_dependencies INTERFACE ${Tasmanian_mklsycl})
    list(APPEND Tasmanian_rpath ${Tasmanian_mklsycl_rpath})
endif()

if (Tasmanian_ENABLE_BLAS)
    target_link_libraries(Tasmanian_dependencies INTERFACE ${BLAS_LIBRARIES})
    target_link_libraries(Tasmanian_dependencies INTERFACE ${LAPACK_LIBRARIES})
    Tasmanian_find_rpath(LIBRARIES ${BLAS_LIBRARIES}   LIST rpath)
    Tasmanian_find_rpath(LIBRARIES ${LAPACK_LIBRARIES} LIST rpath)
endif()

if (Tasmanian_ENABLE_CUDA)
    target_link_libraries(Tasmanian_dependencies INTERFACE ${Tasmanian_cudamathlibs})
    Tasmanian_find_rpath(LIBRARIES ${Tasmanian_cublas} ${Tasmanian_cusparse} ${Tasmaniana_cusolver} LIST rpath)
endif()

if (Tasmanian_ENABLE_HIP)
    target_link_libraries(Tasmanian_dependencies INTERFACE ${Tasmanian_hiplibs})
    list(APPEND Tasmanian_rpath ${Tasmanian_hip_rpath})

    target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${Tasmanian_hiproot}/include/>)
endif()

if (Tasmanian_ENABLE_MAGMA)
    target_link_libraries(Tasmanian_dependencies INTERFACE ${Tasmanian_magmalibs})
    Tasmanian_find_rpath(LIBRARIES ${Tasmanian_magma} LIST rpath)

    target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${Tasmanian_magma_h}/>)
endif()

target_include_directories(Tasmanian_dependencies INTERFACE $<INSTALL_INTERFACE:${Tasmanian_final_install_path}/include>)

target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Config/>)
target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/InterfaceTPL/>)
target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/configured/>)

# Tasmanian_EXTRA_LIBRARIES gives the user an option to force extra dependencies,
# for example, on some systems (e.g., OLCF) find_package(BLAS) fails to
# recognize that libacml_mp requires libgomp, so the build fails with either clang or ENABLE_OPENMP=OFF
# -D Tasmanian_EXTRA_LIBRARIES=/path/to/libgomp.so circumvents the issue
target_link_libraries(Tasmanian_dependencies ${Tasmanian_EXTRA_LIBRARIES})

foreach(_tsg_include ${Tasmanian_EXTRA_INCLUDE_DIRS})
    target_include_directories(Tasmanian_dependencies INTERFACE $<BUILD_INTERFACE:${_tsg_include}>)
    target_include_directories(Tasmanian_dependencies INTERFACE $<INSTALL_INTERFACE:${_tsg_include}>)
endforeach()

install(TARGETS Tasmanian_dependencies EXPORT "${Tasmanian_export_name}")

list(REMOVE_DUPLICATES Tasmanian_rpath)
#message(STATUS "Tasmanian RPATH: ${Tasmanian_rpath}")
