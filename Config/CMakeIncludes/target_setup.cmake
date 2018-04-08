########################################################################
# setup targets
########################################################################

########################################################################
# MaxOSX support, must come first before "add_library()"
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(CMAKE_MACOSX_RPATH ON)
    # code copied from https://cmake.org/Wiki/CMake_RPATH_handling
    # put RPATH in the build three (stage 1, see below)
	set(CMAKE_SKIP_BUILD_RPATH OFF)
	# use the build as opposed install RPATH
	set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
	set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

	set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON) # needed by TPL
endif()

# names of the compiled files
set(Tasmanian_name_libsparsegrid "tasmaniansparsegrid")
set(Tasmanian_name_libdream "tasmaniandream")
set(Tasmanian_name_libsparsegrid_fortran "tasmanianfortran")
set(Tasmanian_name_tasgrid "tasgrid")
set(Tasmanian_name_tasdream "tasdream")
set(Tasmanian_name_fortester "fortester")

include_directories("${PROJECT_SOURCE_DIR}/SparseGrids/") # compile time header files (needed by DREAM)
include_directories("${CMAKE_BINARY_DIR}/configured/") # holds the additional headers configured by cmake

if (Tasmanian_CUDA) # add executables
    cuda_add_executable(Tasmanian_tasgrid ${Tasmanian_source_tasgrid})
    cuda_add_executable(Tasmanian_tasdream ${Tasmanian_source_tasdream})
else()
    add_executable(Tasmanian_tasgrid ${Tasmanian_source_tasgrid})
    add_executable(Tasmanian_tasdream ${Tasmanian_source_tasdream})
endif()

set_target_properties(Tasmanian_tasgrid PROPERTIES OUTPUT_NAME "${Tasmanian_name_tasgrid}")
set_target_properties(Tasmanian_tasdream PROPERTIES OUTPUT_NAME "${Tasmanian_name_tasdream}")

set(Tasmanian_target_list Tasmanian_tasgrid Tasmanian_tasdream) # keep a list of all targets
set(Tasmanian_install_targets Tasmanian_tasgrid Tasmanian_tasdream) # install targets

if (Tasmanian_ENABLE_FORTRAN) # fortran tester
    add_executable(Tasmanian_fortester ${Tasmanian_source_fortester})
    set(Tasmanian_target_list ${Tasmanian_target_list} Tasmanian_fortester)
    set_target_properties(Tasmanian_fortester PROPERTIES OUTPUT_NAME "${Tasmanian_name_fortester}")
endif()

# setup the shared liberaries
if (Tasmanian_SHARED_LIBRARY)
    if (Tasmanian_ENABLE_CUDA)
        cuda_add_library(Tasmanian_libsparsegrid_shared SHARED ${Tasmanian_source_libsparsegrid})
        cuda_add_library(Tasmanian_libdream_shared SHARED ${Tasmanian_source_libdream})
    else()
        add_library(Tasmanian_libsparsegrid_shared SHARED ${Tasmanian_source_libsparsegrid})
        add_library(Tasmanian_libdream_shared SHARED ${Tasmanian_source_libdream})
    endif()

    target_link_libraries(Tasmanian_libdream_shared Tasmanian_libsparsegrid_shared) # DREAM uses Sparse Grids

    set_target_properties(Tasmanian_libsparsegrid_shared PROPERTIES OUTPUT_NAME "${Tasmanian_name_libsparsegrid}")
    set_target_properties(Tasmanian_libdream_shared PROPERTIES OUTPUT_NAME "${Tasmanian_name_libdream}")

    set(Tasmanian_target_list ${Tasmanian_target_list} Tasmanian_libsparsegrid_shared Tasmanian_libdream_shared) # keep track of all targets
    set(Tasmanian_install_targets ${Tasmanian_install_targets} Tasmanian_libsparsegrid_shared Tasmanian_libdream_shared)

    set_property(TARGET Tasmanian_libsparsegrid_shared APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)
    set_property(TARGET Tasmanian_libdream_shared APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)

    if (Tasmanian_ENABLE_FORTRAN)
        add_library(Tasmanian_libsparsegrid_fortran_shared SHARED ${Tasmanian_source_libsparsegrid_fortran})
        target_link_libraries(Tasmanian_libsparsegrid_fortran_shared Tasmanian_libsparsegrid_shared)

        set(Tasmanian_install_targets ${Tasmanian_install_targets} Tasmanian_libsparsegrid_fortran_shared)
        set(Tasmanian_target_list ${Tasmanian_target_list} Tasmanian_libsparsegrid_fortran_shared)

        set_property(TARGET Tasmanian_libsparsegrid_fortran_shared APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)

        set_target_properties(Tasmanian_libsparsegrid_fortran_shared PROPERTIES OUTPUT_NAME "${Tasmanian_name_libsparsegrid_fortran}")
    endif()
endif()

if (Tasmanian_STATIC_LIBRARY)
    if (Tasmanian_ENABLE_CUDA)
        cuda_add_library(Tasmanian_libsparsegrid_static STATIC ${Tasmanian_source_libsparsegrid})
        cuda_add_library(Tasmanian_libdream_static STATIC ${Tasmanian_source_libdream})
    else()
        add_library(Tasmanian_libsparsegrid_static STATIC ${Tasmanian_source_libsparsegrid})
        add_library(Tasmanian_libdream_static STATIC ${Tasmanian_source_libdream})
    endif()

    set(Tasmanian_target_list ${Tasmanian_target_list} Tasmanian_libsparsegrid_static Tasmanian_libdream_static)
    target_link_libraries(Tasmanian_libdream_static Tasmanian_libsparsegrid_static)
    if (Tasmanian_SHARED_LIBRARY)
        # without this, if both static and shared libs are enabled, make -j tries to compile shared cuda kernels twice
        # which creates a race condition and the build randomly fails
        add_dependencies(Tasmanian_libsparsegrid_shared Tasmanian_libsparsegrid_static)
    endif()

    set_target_properties(Tasmanian_libsparsegrid_static PROPERTIES OUTPUT_NAME "${Tasmanian_name_libsparsegrid}")
    set_target_properties(Tasmanian_libdream_static PROPERTIES OUTPUT_NAME "${Tasmanian_name_libdream}")

    set(Tasmanian_install_targets ${Tasmanian_install_targets} Tasmanian_libsparsegrid_static Tasmanian_libdream_static)

    # executables prefer linking to static libraries
    target_link_libraries(Tasmanian_tasgrid Tasmanian_libsparsegrid_static)
    target_link_libraries(Tasmanian_tasdream Tasmanian_libdream_static Tasmanian_libsparsegrid_static)

    set_property(TARGET Tasmanian_libsparsegrid_static APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)
    set_property(TARGET Tasmanian_libdream_static APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)

    if (Tasmanian_ENABLE_FORTRAN)
        add_library(Tasmanian_libsparsegrid_fortran_static STATIC ${Tasmanian_source_libsparsegrid_fortran})
        target_link_libraries(Tasmanian_libsparsegrid_fortran_static Tasmanian_libsparsegrid_static)

        set(Tasmanian_target_list ${Tasmanian_target_list} Tasmanian_libsparsegrid_fortran_static)
        set(Tasmanian_install_targets ${Tasmanian_install_targets} Tasmanian_libsparsegrid_fortran_static)

        set_property(TARGET Tasmanian_libsparsegrid_fortran_static APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)

        set_target_properties(Tasmanian_libsparsegrid_fortran_static PROPERTIES OUTPUT_NAME "${Tasmanian_name_libsparsegrid_fortran}")
        target_link_libraries(Tasmanian_fortester Tasmanian_libsparsegrid_fortran_static)
    endif()
else()
    target_link_libraries(Tasmanian_tasgrid Tasmanian_libsparsegrid_shared)
    target_link_libraries(Tasmanian_tasdream Tasmanian_libdream_shared Tasmanian_libsparsegrid_shared)
    if (Tasmanian_ENABLE_FORTRAN)
        target_link_libraries(Tasmanian_fortester Tasmanian_libsparsegrid_fortran_shared)
    endif()
endif()

# for development purposes, compile the exmaples together with everything else
if (Tasmanian_ENABLE_DEVELOPMENT_DEFAULTS)
    if (Tasmanian_ENABLE_CUDA)
        cuda_add_executable(Tasmanian_example_sparse_grids "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.cpp")
        cuda_add_executable(Tasmanian_example_dream "${PROJECT_SOURCE_DIR}/Examples/example_dream.cpp")
    else()
        add_executable(Tasmanian_example_sparse_grids "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.cpp")
        add_executable(Tasmanian_example_dream "${PROJECT_SOURCE_DIR}/Examples/example_dream.cpp")
    endif()
    set(Tasmanian_target_list ${Tasmanian_target_list} Tasmanian_example_sparse_grids Tasmanian_example_dream)
    include_directories("${PROJECT_SOURCE_DIR}/DREAM/") # examples need the DREAM include folder
    if (Tasmanian_STATIC_LIBRARY)
        target_link_libraries(Tasmanian_example_sparse_grids Tasmanian_libsparsegrid_static)
        target_link_libraries(Tasmanian_example_dream Tasmanian_libdream_static)
    else()
        target_link_libraries(Tasmanian_example_sparse_grids Tasmanian_libsparsegrid_shared)
        target_link_libraries(Tasmanian_example_dream Tasmanian_libdream_shared)
    endif()
    if (Tasmanian_ENABLE_FORTRAN)
        add_executable(Tasmanian_example_sparse_grids_f90 "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.f90")
        set_target_properties(Tasmanian_example_sparse_grids_f90 PROPERTIES OUTPUT_NAME "example_sparse_grids_fortran")
        set(Tasmanian_target_list ${Tasmanian_target_list} Tasmanian_example_sparse_grids_f90)
        if (Tasmanian_STATIC_LIBRARY)
            target_link_libraries(Tasmanian_example_sparse_grids_f90 Tasmanian_libsparsegrid_fortran_static Tasmanian_libdream_static)
        else()
            target_link_libraries(Tasmanian_example_sparse_grids_f90 Tasmanian_libsparsegrid_fortran_shared Tasmanian_libsparsegrid_shared)
        endif()
    endif()
endif()


########################################################################
# Windows specific support (DLL export/import directives and names)
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    configure_file("${PROJECT_SOURCE_DIR}/SparseGrids/TasmanianSparseGrid.windows.hpp"  "${PROJECT_SOURCE_DIR}/SparseGrids/TasmanianSparseGrid.hpp" COPYONLY)

    if (Tasmanian_SHARED_LIBRARY)
        target_compile_definitions(Tasmanian_libsparsegrid_shared PUBLIC -DTSG_DLL)
        target_compile_definitions(Tasmanian_libdream_shared PUBLIC -DTSG_DLL)
    else()
        target_compile_definitions(Tasmanian_tasgrid PUBLIC -DTSG_DYNAMIC)
        target_compile_definitions(Tasmanian_tasdream PUBLIC -DTSG_DYNAMIC)
    endif()

    if (Tasmanian_STATIC_LIBRARY)
        target_compile_definitions(Tasmanian_libsparsegrid_static PUBLIC -DTSG_STATIC)
        target_compile_definitions(Tasmanian_libdream_static PUBLIC -DTSG_STATIC)

        # Unix uses .a and .so to distinguish static and shared
        # Windows uses .lib in both cases, append _static to differentiate
        set_target_properties(Tasmanian_libsparsegrid_static PROPERTIES OUTPUT_NAME "${Tasmanian_name_libsparsegrid}_static")
        set_target_properties(Tasmanian_libdream_static PROPERTIES OUTPUT_NAME "${Tasmanian_name_libdream}_static")

        if (Tasmanian_ENABLE_FORTRAN)
            set_target_properties(Tasmanian_libsparsegrid_fortran_static PROPERTIES OUTPUT_NAME "${Tasmanian_name_libsparsegrid_fortran}_static")
        endif()
    endif()

    add_definitions(-D_SCL_SECURE_NO_WARNINGS) # suppresses warnings regarding pointers to the middle of an array
    add_definitions(-D_TASMANIAN_WINDOWS_) # overwrittes gettime() and sets a seed
    add_definitions(-D_USE_MATH_DEFINES) # needed to include M_PI constant
endif()

