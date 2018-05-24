########################################################################
# Setup the specific options
########################################################################

macro(Tasmanian_target_link_libraries_sparsegrid Tasmanian_target_libraries)
    if (Tasmanian_SHARED_LIBRARY)
        target_link_libraries(Tasmanian_libsparsegrid_shared ${Tasmanian_target_libraries})
    endif()
    if (Tasmanian_STATIC_LIBRARY)
        target_link_libraries(Tasmanian_libsparsegrid_static ${Tasmanian_target_libraries})
    endif()
endmacro()

macro(Tasmanian_target_link_libraries_dream Tasmanian_target_libraries)
    if (Tasmanian_SHARED_LIBRARY)
        target_link_libraries(Tasmanian_libdream_shared ${Tasmanian_target_libraries})
    endif()
    if (Tasmanian_STATIC_LIBRARY)
        target_link_libraries(Tasmanian_libdream_static ${Tasmanian_target_libraries})
    endif()
endmacro()

macro(Tasmanian_target_include_directories_sparsegrid Tasmanian_option Tasmanian_target_dirs) # option PUBLIC/PRIVATE
    if (Tasmanian_SHARED_LIBRARY)
        target_include_directories(Tasmanian_libsparsegrid_shared ${Tasmanian_option} ${Tasmanian_target_dirs})
    endif()
    if (Tasmanian_STATIC_LIBRARY)
        target_include_directories(Tasmanian_libsparsegrid_static ${Tasmanian_option} ${Tasmanian_target_dirs})
    endif()
endmacro()

macro(Tasmanian_target_include_directories_dream Tasmanian_option Tasmanian_target_dir_options) # option: PUBLIC/PRIVATE
    if (Tasmanian_SHARED_LIBRARY)
        target_include_directories(Tasmanian_libdream_shared ${Tasmanian_option} ${Tasmanian_target_dir_options})
    endif()
    if (Tasmanian_STATIC_LIBRARY)
        target_include_directories(Tasmanian_libdream_static ${Tasmanian_option} ${Tasmanian_target_dir_options})
    endif()
endmacro()

macro(Tasmanian_target_compile_options_sparsegrid Tasmanian_option Tasmanian_compile_options) # option PUBLIC/PRIVATE
    if (Tasmanian_SHARED_LIBRARY)
        target_compile_options(Tasmanian_libsparsegrid_shared ${Tasmanian_option} ${Tasmanian_compile_options})
    endif()
    if (Tasmanian_STATIC_LIBRARY)
        target_compile_options(Tasmanian_libsparsegrid_static ${Tasmanian_option} ${Tasmanian_compile_options})
    endif()
endmacro()

macro(Tasmanian_target_compile_options_dream Tasmanian_option Tasmanian_compile_options) # option: PUBLIC/PRIVATE
    if (Tasmanian_SHARED_LIBRARY)
        target_compile_options(Tasmanian_libdream_shared ${Tasmanian_option} ${Tasmanian_compile_options})
    endif()
    if (Tasmanian_STATIC_LIBRARY)
        target_compile_options(Tasmanian_libdream_static ${Tasmanian_option} ${Tasmanian_compile_options})
    endif()
endmacro()

########################################################################
# BLAS setup
########################################################################
if (Tasmanian_ENABLE_BLAS)
    Tasmanian_target_link_libraries_sparsegrid(${BLAS_LIBRARIES})
    Tasmanian_target_link_libraries_dream(${BLAS_LIBRARIES})
endif()

########################################################################
# MPI setup (experimental, DREAM distributed posterior only)
########################################################################
if (Tasmanian_ENABLE_MPI)
    Tasmanian_target_link_libraries_dream(${MPI_CXX_LIBRARIES})

    if (DEFINED MPI_CXX_INCLUDE_PATH)
        Tasmanian_target_include_directories_dream("PUBLIC ${MPI_CXX_INCLUDE_PATH}")
    endif()

    if(DEFINED MPI_CXX_COMPILE_FLAGS)
        target_compile_options(Tasmanian_tasdream PUBLIC ${MPI_CXX_COMPILE_FLAGS})
        Tasmanian_target_compile_options_sparsegrid(PUBLIC ${MPI_CXX_COMPILE_FLAGS})
    endif()

    if(DEFINED MPI_CXX_LINK_FLAGS)
        set_target_properties(Tasmanian_tasdream PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
        if (Tasmanian_SHARED_LIBRARY)
            set_target_properties(Tasmanian_libdream_shared PROPERTIES  LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
        endif()
        if (Tasmanian_STATIC_LIBRARY)
            set_target_properties(Tasmanian_libdream_static PROPERTIES  LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
        endif()
    endif()
endif()

########################################################################
# CUDA setup (see also the cuda pre-process in the sanity check section)
########################################################################
if (Tasmanian_ENABLE_CUBLAS)
    Tasmanian_target_link_libraries_sparsegrid("${CUDA_cusparse_LIBRARY};${CUDA_CUBLAS_LIBRARIES};${CUDA_LIBRARIES}")

    if (DEFINED CUDA_INCLUDE_DIRS)
        Tasmanian_target_include_directories_sparsegrid(PUBLIC ${CUDA_INCLUDE_DIRS})
    endif()
endif()

########################################################################
# OpenMP setup
########################################################################
if (Tasmanian_ENABLE_OPENMP)
    if (OpenMP_CXX_LIBRARIES)
        # using the OpenMP target leads to a problem with the exports
        # the OpenMP target cannot be exported, which means that a
        # project importing an already installed Tasmanian would
        # have to "know" whether Tasmanian was build with OpenMP and
        # call find_package(OpenMP) manually
        # Furthermore, using find_package(OpenMP) from a different
        # compiler can generate a wrong target, e.g., building Tasmanian
        # with gcc links to libgomp, but calling find_package(OpenMP)
        # from clang will create an OpenMP target that uses libiomp
        Tasmanian_target_link_libraries_sparsegrid(${OpenMP_CXX_LIBRARIES})
        target_compile_options(Tasmanian_tasgrid PRIVATE ${OpenMP_CXX_FLAGS})
        Tasmanian_target_compile_options_sparsegrid(PRIVATE ${OpenMP_CXX_FLAGS})
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
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

