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

########################################################################
# first define all targets including two alias targets
# both _shared or _static libraries can be build,
# in which case executables link statically, hence the alias targets
########################################################################
if (Tasmanian_CUDA) # add executables for testing and cli interface
    cuda_add_executable(Tasmanian_tasgrid ${Tasmanian_source_tasgrid})
    cuda_add_executable(Tasmanian_tasdream ${Tasmanian_source_tasdream})
    cuda_add_executable(Tasmanian_example_sparse_grids "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.cpp")
    cuda_add_executable(Tasmanian_example_dream "${PROJECT_SOURCE_DIR}/Examples/example_dream.cpp")
else()
    add_executable(Tasmanian_tasgrid ${Tasmanian_source_tasgrid})
    add_executable(Tasmanian_tasdream ${Tasmanian_source_tasdream})
    add_executable(Tasmanian_example_sparse_grids "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.cpp")
    add_executable(Tasmanian_example_dream "${PROJECT_SOURCE_DIR}/Examples/example_dream.cpp")
endif()

if (Tasmanian_ENABLE_FORTRAN) # fortran tests and examples
    add_executable(Tasmanian_fortester "${PROJECT_SOURCE_DIR}/Testing/fortester.f90")
    add_executable(Tasmanian_example_sparse_grids_f90 "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.f90")
endif()

# form the lists of targets
set(Tasmanian_install_targets Tasmanian_tasgrid Tasmanian_tasdream) # install targets

if (Tasmanian_SHARED_LIBRARY) # add the shared libraries
    if (Tasmanian_ENABLE_CUDA)
        cuda_add_library(Tasmanian_libsparsegrid_shared SHARED ${Tasmanian_source_libsparsegrid})
        cuda_add_library(Tasmanian_libdream_shared SHARED ${Tasmanian_source_libdream})
    else()
        add_library(Tasmanian_libsparsegrid_shared SHARED ${Tasmanian_source_libsparsegrid})
        add_library(Tasmanian_libdream_shared SHARED ${Tasmanian_source_libdream})
    endif()

    set(Tasmanian_install_targets ${Tasmanian_install_targets} Tasmanian_libsparsegrid_shared Tasmanian_libdream_shared)
    target_link_libraries(Tasmanian_libdream_shared Tasmanian_libsparsegrid_shared)

    if (Tasmanian_ENABLE_FORTRAN)
        add_library(Tasmanian_libsparsegrid_fortran_shared SHARED ${Tasmanian_source_libsparsegrid_fortran})

        set(Tasmanian_install_targets ${Tasmanian_install_targets} Tasmanian_libsparsegrid_fortran_shared)
        target_link_libraries(Tasmanian_libsparsegrid_fortran_shared Tasmanian_libsparsegrid_shared)
    endif()
endif()

if (Tasmanian_STATIC_LIBRARY) # add the static libraries
    if (Tasmanian_ENABLE_CUDA)
        cuda_add_library(Tasmanian_libsparsegrid_static STATIC ${Tasmanian_source_libsparsegrid})
        cuda_add_library(Tasmanian_libdream_static STATIC ${Tasmanian_source_libdream})
    else()
        add_library(Tasmanian_libsparsegrid_static STATIC ${Tasmanian_source_libsparsegrid})
        add_library(Tasmanian_libdream_static STATIC ${Tasmanian_source_libdream})
    endif()

    set(Tasmanian_install_targets ${Tasmanian_install_targets} Tasmanian_libsparsegrid_static Tasmanian_libdream_static)
    target_link_libraries(Tasmanian_libdream_static Tasmanian_libsparsegrid_static)

    add_library(Tasmanian_libsparsegrid ALIAS Tasmanian_libsparsegrid_static) # link executables statically
    add_library(Tasmanian_libdream ALIAS Tasmanian_libdream_static)

    if (Tasmanian_ENABLE_FORTRAN)
        add_library(Tasmanian_libsparsegrid_fortran_static STATIC ${Tasmanian_source_libsparsegrid_fortran})

        set(Tasmanian_install_targets ${Tasmanian_install_targets} Tasmanian_libsparsegrid_fortran_static)
        target_link_libraries(Tasmanian_libsparsegrid_fortran_static Tasmanian_libsparsegrid_static)

        add_library(Tasmanian_libfortran ALIAS Tasmanian_libsparsegrid_fortran_static)
    endif()
else()
    add_library(Tasmanian_libsparsegrid ALIAS Tasmanian_libsparsegrid_shared) # link executables dynamically (no static library)
    add_library(Tasmanian_libdream ALIAS Tasmanian_libdream_shared)

    if (Tasmanian_ENABLE_FORTRAN)
        add_library(Tasmanian_libfortran ALIAS Tasmanian_libsparsegrid_fortran_shared)
    endif()
endif()

# hack a dependency problem in parallel make
# without this, if both static and shared libs are enabled, make -j tries to compile shared cuda kernels twice
# which creates a race condition and the build randomly fails
if (Tasmanian_SHARED_LIBRARY AND Tasmanian_STATIC_LIBRARY)
    add_dependencies(Tasmanian_libsparsegrid_shared Tasmanian_libsparsegrid_static)
endif()

########################################################################
# give names for the executables and libraries
########################################################################
set_target_properties(Tasmanian_tasgrid PROPERTIES OUTPUT_NAME "tasgrid")
set_target_properties(Tasmanian_tasdream PROPERTIES OUTPUT_NAME "tasdream")
set_target_properties(Tasmanian_example_sparse_grids PROPERTIES OUTPUT_NAME "example_sparse_grids")
set_target_properties(Tasmanian_example_dream PROPERTIES OUTPUT_NAME "example_dream")
if (Tasmanian_ENABLE_FORTRAN)
    set_target_properties(Tasmanian_fortester PROPERTIES OUTPUT_NAME "fortester")
    set_target_properties(Tasmanian_example_sparse_grids_f90 PROPERTIES OUTPUT_NAME "example_sparse_grids_f90")
endif()
if (Tasmanian_SHARED_LIBRARY)
    set_target_properties(Tasmanian_libsparsegrid_shared PROPERTIES OUTPUT_NAME "tasmaniansparsegrid")
    set_target_properties(Tasmanian_libdream_shared PROPERTIES OUTPUT_NAME "tasmaniandream")
    if (Tasmanian_ENABLE_FORTRAN)
        set_target_properties(Tasmanian_libsparsegrid_fortran_shared PROPERTIES OUTPUT_NAME "tasmanianfortran")
    endif()
endif()
if (Tasmanian_STATIC_LIBRARY)
    set_target_properties(Tasmanian_libsparsegrid_static PROPERTIES OUTPUT_NAME "tasmaniansparsegrid")
    set_target_properties(Tasmanian_libdream_static PROPERTIES OUTPUT_NAME "tasmaniandream")
    if (Tasmanian_ENABLE_FORTRAN)
        set_target_properties(Tasmanian_libsparsegrid_fortran_static PROPERTIES OUTPUT_NAME "tasmanianfortran")
    endif()
endif()

########################################################################
# link executable to libraries (libraries are linked to libraries
# when the targets are created
########################################################################
target_link_libraries(Tasmanian_tasgrid Tasmanian_libsparsegrid)
target_link_libraries(Tasmanian_tasdream Tasmanian_libdream)
target_link_libraries(Tasmanian_example_sparse_grids Tasmanian_libsparsegrid)
target_link_libraries(Tasmanian_example_dream Tasmanian_libdream)
if (Tasmanian_ENABLE_FORTRAN)
    target_link_libraries(Tasmanian_fortester Tasmanian_libfortran)
    target_link_libraries(Tasmanian_example_sparse_grids_f90 Tasmanian_libfortran)
endif()


########################################################################
# add include directories for BUILD and INSTALL interfaces
########################################################################
if (Tasmanian_SHARED_LIBRARY)
    target_include_directories(Tasmanian_libsparsegrid_shared PUBLIC $<INSTALL_INTERFACE:include>)
    target_include_directories(Tasmanian_libdream_shared PUBLIC $<INSTALL_INTERFACE:include>)

    target_include_directories(Tasmanian_libsparsegrid_shared PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/SparseGrids/>)
    target_include_directories(Tasmanian_libsparsegrid_shared PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/configured/>)
    target_include_directories(Tasmanian_libdream_shared PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/DREAM/>)
    target_include_directories(Tasmanian_libdream_shared PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/configured/>)

    if (Tasmanian_ENABLE_FORTRAN)
        target_include_directories(Tasmanian_libsparsegrid_fortran_shared PUBLIC $<INSTALL_INTERFACE:include>)
        target_include_directories(Tasmanian_libdream_shared PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>)
    endif()
endif()

if (Tasmanian_STATIC_LIBRARY)
    target_include_directories(Tasmanian_libsparsegrid_static PUBLIC $<INSTALL_INTERFACE:include>)
    target_include_directories(Tasmanian_libdream_static PUBLIC $<INSTALL_INTERFACE:include>)

    target_include_directories(Tasmanian_libsparsegrid_static PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/SparseGrids/>)
    target_include_directories(Tasmanian_libsparsegrid_static PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/configured/>)
    target_include_directories(Tasmanian_libdream_static PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/DREAM/>)
    target_include_directories(Tasmanian_libdream_static PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/configured/>)

    if (Tasmanian_ENABLE_FORTRAN)
        target_include_directories(Tasmanian_libsparsegrid_fortran_static PUBLIC $<INSTALL_INTERFACE:include>)
        target_include_directories(Tasmanian_libdream_static PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>)
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
        set_target_properties(Tasmanian_libsparsegrid_static PROPERTIES OUTPUT_NAME "tasmaniansparsegrid_static")
        set_target_properties(Tasmanian_libdream_static PROPERTIES OUTPUT_NAME "tasmaniandream_static")

        if (Tasmanian_ENABLE_FORTRAN)
            set_target_properties(Tasmanian_libsparsegrid_fortran_static PROPERTIES OUTPUT_NAME "tasmanianfortran_static")
        endif()
    endif()

    target_compile_definitions(Tasmanian_tasgrid -D_TASMANIAN_WINDOWS_) # overwrittes gettime()
    add_definitions(-D_SCL_SECURE_NO_WARNINGS) # suppresses warnings regarding pointers to the middle of an array
    add_definitions(-D_USE_MATH_DEFINES) # needed to include M_PI constant (lots of targets need this, will figure it out)
endif()

