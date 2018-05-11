########################################################################
# Setup the specific options
########################################################################

########################################################################
# BLAS setup
########################################################################
if (Tasmanian_ENABLE_BLAS)
    foreach(Tasmanian_loop_target ${Tasmanian_target_list})
        target_link_libraries(${Tasmanian_loop_target} ${BLAS_LIBRARIES})
    endforeach()
endif()

########################################################################
# MPI setup (experimental, DREAM distributed posterior only)
########################################################################
if (Tasmanian_ENABLE_MPI)
    target_link_libraries(tasdream ${MPI_CXX_LIBRARIES})
    if (Tasmanian_SHARED_LIBRARY)
        target_link_libraries(Tasmanian_libdream_shared ${MPI_CXX_LIBRARIES})
    endif()
    if (Tasmanian_STATIC_LIBRARY)
        target_link_libraries(Tasmanian_libdream_static ${MPI_CXX_LIBRARIES})
    endif()

    if (DEFINED MPI_CXX_INCLUDE_PATH)
        include_directories(${MPI_CXX_INCLUDE_PATH})
        if (Tasmanian_SHARED_LIBRARY)
            set_property(TARGET Tasmanian_libdream_shared APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_PATH}")
        endif()
        if (Tasmanian_STATIC_LIBRARY)
            set_property(TARGET Tasmanian_libdream_static APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_PATH}")
        endif()
    endif()

    if(DEFINED MPI_COMPILE_FLAGS)
        set_target_properties(Tasmanian_tasdream PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
        if (Tasmanian_SHARED_LIBRARY)
            set_target_properties(Tasmanian_libdream_shared PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
        endif()
        if (Tasmanian_STATIC_LIBRARY)
            set_target_properties(Tasmanian_libdream_static PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
        endif()
    endif()

    if(DEFINED MPI_LINK_FLAGS)
        set_target_properties(tasdream PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
        if (Tasmanian_SHARED_LIBRARY)
            set_target_properties(Tasmanian_libdream_shared PROPERTIES  LINK_FLAGS "${MPI_LINK_FLAGS}")
        endif()
        if (Tasmanian_STATIC_LIBRARY)
            set_target_properties(Tasmanian_libdream_static PROPERTIES  LINK_FLAGS "${MPI_LINK_FLAGS}")
        endif()
    endif()
endif()

########################################################################
# CUDA setup (see also the cuda pre-process in the sanity check section)
########################################################################
if (Tasmanian_ENABLE_CUBLAS)
    foreach(Tasmanian_loop_target ${Tasmanian_target_list})
        target_link_libraries(${Tasmanian_loop_target} ${CUDA_cusparse_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES})
    endforeach()

    if (DEFINED CUDA_INCLUDE_DIRS)
        include_directories(${CUDA_INCLUDE_DIRS})
    endif()
endif()

########################################################################
# OpenMP setup
########################################################################
if (Tasmanian_ENABLE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

########################################################################
# Python install setup
########################################################################
if (Tasmanian_ENABLE_PYTHON)
    set(Tasmanian_python_install_path "${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")
endif()

########################################################################
# Cholmod support (experimental, wavelets only) (deosn't work better)
########################################################################
#if (TPL_ENABLE_CHOLMOD OR (DEFINED TPL_CHOLMOD_INCLUDE_DIRS) OR (DEFINED TPL_CHOLMOD_LIBRARIES))
#    add_definitions(-DTASMANIAN_CHOLMOD)
#    if ((NOT DEFINED TPL_CHOLMOD_INCLUDE_DIRS) AND (NOT DEFINED TPL_CHOLMOD_LIBRARIES))
#        set(TPL_CHOLMOD_LIBRARIES cholmod)
#        set(TPL_CHOLMOD_INCLUDE_DIRS "")
#    endif()
#
#    include_directories(${TPL_CHOLMOD_INCLUDE_DIRS})
#
#    foreach(tas_target ${tas_target_list})
#        target_link_libraries(${tas_target} ${TPL_CHOLMOD_LIBRARIES})
#    endforeach()
#endif()

