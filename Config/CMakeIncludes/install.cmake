########################################################################
# Installing
########################################################################
install(TARGETS ${Tasmanian_install_targets}
        EXPORT Tasmanian_all_export_targets
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")
install(EXPORT Tasmanian_all_export_targets DESTINATION "config" FILE "Tasmanian.cmake")

install(DIRECTORY "${CMAKE_SOURCE_DIR}/SparseGrids/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "*.windows.*" EXCLUDE
        PATTERN "Example" EXCLUDE
        PATTERN "tsgHiddenExternals.hpp" EXCLUDE)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/DREAM/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "*.windows.*" EXCLUDE
        PATTERN "*.in.*" EXCLUDE
        PATTERN "Example" EXCLUDE)
install(DIRECTORY "${CMAKE_BINARY_DIR}/configured/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")
install(FILES "${CMAKE_SOURCE_DIR}/SparseGrids/TasmanianSparseGrid.h"
        DESTINATION "include"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(FILES "${CMAKE_BINARY_DIR}/install/examples/CMakeLists.txt"
        DESTINATION "examples"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/Examples/"
        DESTINATION "examples"
        FILES_MATCHING PATTERN "*.cpp"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(DIRECTORY "${CMAKE_BINARY_DIR}/install/config/"
        DESTINATION "config"
        FILES_MATCHING PATTERN "Tasmanian*setup.sh"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

if (Tasmanian_ENABLE_PYTHON)
    install(FILES "${CMAKE_BINARY_DIR}/install/python/TasmanianSG.py"
            DESTINATION "${Tasmanian_python_install_path}"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    # Create symlink for backward compatibility
    install(CODE "execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${Tasmanian_python_install_path} ${CMAKE_INSTALL_PREFIX}/python )" )
    install(FILES "${CMAKE_BINARY_DIR}/install/examples/example_sparse_grids.py"
            DESTINATION "examples"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
endif()

if (Tasmanian_ENABLE_MATLAB)
    install(DIRECTORY "${CMAKE_SOURCE_DIR}/InterfaceMATLAB/"
            DESTINATION "matlab"
            FILES_MATCHING PATTERN "*.m"
            PATTERN "tsgGetPaths.*" EXCLUDE)
    install(FILES "${CMAKE_BINARY_DIR}/install/matlab/tsgGetPaths.m"
            DESTINATION "matlab"
            PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
endif()

if (Tasmanian_ENABLE_FORTRAN)
    install(FILES "${CMAKE_BINARY_DIR}/tasmaniansg.mod"
            DESTINATION "include"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    install(FILES "${CMAKE_SOURCE_DIR}/Examples/example_sparse_grids.f90"
            DESTINATION "examples"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
endif()


########################################################################
# Install Message
########################################################################
set(Tasmanian_post_install "--------------------------------------------------------------------------------\n")
set(Tasmanian_post_install "${Tasmanian_post_install}Tasmanian Sparse Grids install complete\n\n")
set(Tasmanian_post_install "${Tasmanian_post_install}executables:      ${Tasmanian_name_tasgrid}\n")
set(Tasmanian_post_install "${Tasmanian_post_install}                  ${Tasmanian_name_tasdream}\n")
if (Tasmanian_SHARED_LIBRARY)
    set(Tasmanian_post_install "${Tasmanian_post_install}shared libraries: ${CMAKE_SHARED_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid}${CMAKE_SHARED_LIBRARY_SUFFIX}\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}                  ${CMAKE_SHARED_LIBRARY_PREFIX}${Tasmanian_name_libdream}${CMAKE_SHARED_LIBRARY_SUFFIX}\n")
    if (Tasmanian_ENABLE_FORTRAN)
        set(Tasmanian_post_install "${Tasmanian_post_install}                  ${CMAKE_SHARED_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid_fortran}${CMAKE_SHARED_LIBRARY_SUFFIX}\n")
    endif()
endif()
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    if (Tasmanian_STATIC_LIBRARY)
        set(Tasmanian_post_install "${Tasmanian_post_install}static libraries: ${CMAKE_STATIC_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid}_static${CMAKE_STATIC_LIBRARY_SUFFIX}\n")
        set(Tasmanian_post_install "${Tasmanian_post_install}                  ${CMAKE_STATIC_LIBRARY_PREFIX}${Tasmanian_name_libdream}_static${CMAKE_STATIC_LIBRARY_SUFFIX}\n")
        if (Tasmanian_ENABLE_FORTRAN)
            set(Tasmanian_post_install "${Tasmanian_post_install}                  ${CMAKE_STATIC_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid_fortran}_static${CMAKE_STATIC_LIBRARY_SUFFIX}\n")
        endif()
    endif()
else()
    if (Tasmanian_STATIC_LIBRARY)
        set(Tasmanian_post_install "${Tasmanian_post_install}static libraries: ${CMAKE_STATIC_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid}${CMAKE_STATIC_LIBRARY_SUFFIX}\n")
        set(Tasmanian_post_install "${Tasmanian_post_install}                  ${CMAKE_STATIC_LIBRARY_PREFIX}${Tasmanian_name_libdream}${CMAKE_STATIC_LIBRARY_SUFFIX}\n")
        if (Tasmanian_ENABLE_FORTRAN)
            set(Tasmanian_post_install "${Tasmanian_post_install}                  ${CMAKE_STATIC_LIBRARY_PREFIX}${Tasmanian_name_libsparsegrid_fortran}${CMAKE_STATIC_LIBRARY_SUFFIX}\n")
        endif()
    endif()
endif()
set(Tasmanian_post_install "${Tasmanian_post_install}C++ examples:     ${CMAKE_INSTALL_PREFIX}/examples/example_sparse_grids.cpp\n")
set(Tasmanian_post_install "${Tasmanian_post_install}                  ${CMAKE_INSTALL_PREFIX}/examples/example_dream.cpp\n")
if (Tasmanian_ENABLE_FORTRAN)
    set(Tasmanian_post_install "${Tasmanian_post_install}                  ${CMAKE_INSTALL_PREFIX}/examples/example_sparse_grids.f90\n")
endif()
set(Tasmanian_post_install "${Tasmanian_post_install}                  see CMakeLists.txt in ${CMAKE_INSTALL_PREFIX}/examples\n\n")

set(Tasmanian_post_install "${Tasmanian_post_install}bash export commands (add to your environment, i.e., ~/.bashrc):\n\n")
set(Tasmanian_post_install "${Tasmanian_post_install}export PATH=$PATH:${CMAKE_INSTALL_PREFIX}/bin/\n")
set(Tasmanian_post_install "${Tasmanian_post_install}export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib/\n")
set(Tasmanian_post_install "${Tasmanian_post_install}export C_INCLUDE_PATH=$C_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/\n")
set(Tasmanian_post_install "${Tasmanian_post_install}export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/\n")
set(Tasmanian_post_install "${Tasmanian_post_install}export LIBRARY_PATH=$LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib/\n\n")

set(Tasmanian_post_install "${Tasmanian_post_install}above commands are writte in the following\n")
set(Tasmanian_post_install "${Tasmanian_post_install}bash source files: ${CMAKE_INSTALL_PREFIX}/config/TasmanianENVsetup.sh\n")
set(Tasmanian_post_install "${Tasmanian_post_install}                   ${CMAKE_INSTALL_PREFIX}/config/TasmanianDEVsetup.sh\n\n")
set(Tasmanian_post_install "${Tasmanian_post_install}cmake import:             ${CMAKE_INSTALL_PREFIX}/config/Tasmanian.cmake\n")
set(Tasmanian_post_install "${Tasmanian_post_install}cmake directive:          include(\\\"${CMAKE_INSTALL_PREFIX}/config/Tasmanian.cmake\\\")\n")
set(Tasmanian_post_install "${Tasmanian_post_install}cmake executable targets: Tasmanian_tasgrid\n")
set(Tasmanian_post_install "${Tasmanian_post_install}                          Tasmanian_tasdream\n")
if (Tasmanian_ENABLE_FORTRAN)
    if (Tasmanian_STATIC_LIBRARY)
        set(Tasmanian_post_install "${Tasmanian_post_install}cmake static targets:     Tasmanian_libsparsegrid_static\n")
        set(Tasmanian_post_install "${Tasmanian_post_install}                          Tasmanian_libdream_static\n")
        set(Tasmanian_post_install "${Tasmanian_post_install}                          Tasmanian_libsparsegrid_fortran_static\n")
    endif()
    if (Tasmanian_SHARED_LIBRARY)
        set(Tasmanian_post_install "${Tasmanian_post_install}cmake shared targets:     Tasmanian_libsparsegrid_shared\n")
        set(Tasmanian_post_install "${Tasmanian_post_install}                          Tasmanian_libdream_shared\n")
        set(Tasmanian_post_install "${Tasmanian_post_install}                          Tasmanian_libsparsegrid_fortran_shared\n\n")
    endif()
else()
    if (Tasmanian_STATIC_LIBRARY)
        set(Tasmanian_post_install "${Tasmanian_post_install}cmake static targets:     Tasmanian_libsparsegrid_static\n")
        set(Tasmanian_post_install "${Tasmanian_post_install}                          Tasmanian_libdream_static\n")
    endif()
    if (Tasmanian_SHARED_LIBRARY)
        set(Tasmanian_post_install "${Tasmanian_post_install}cmake shared targets:     Tasmanian_libsparsegrid_shared\n")
        set(Tasmanian_post_install "${Tasmanian_post_install}                          Tasmanian_libdream_shared\n\n")
    endif()
endif()

if (Tasmanian_ENABLE_PYTHON)
    set(Tasmanian_post_install "${Tasmanian_post_install}python module:  ${Tasmanian_python_install_path}/TasmanianSG.py\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}python example: ${CMAKE_INSTALL_PREFIX}/examples/example_sparse_grids.py\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}                to call the python module\n\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}                import sys\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}                sys.path.append(\\\"${CMAKE_INSTALL_PREFIX}/python\\\")\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}                import TasmanianSG\n\n")
endif()

if (Tasmanian_ENABLE_MATLAB)
    set(Tasmanian_post_install "${Tasmanian_post_install}matlab flles in: ${CMAKE_INSTALL_PREFIX}/matlab/\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}                 NOTE: add ${CMAKE_INSTALL_PREFIX}/matlab/ to the MATLAB path\n\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}                 addpath('${CMAKE_INSTALL_PREFIX}/matlab/')\n")
endif()
set(Tasmanian_post_install "${Tasmanian_post_install}--------------------------------------------------------------------------------\n")

install(CODE "message(\"${Tasmanian_post_install}\")")
install(CODE "file(WRITE ${CMAKE_INSTALL_PREFIX}/config/Tasmanian.log \"${Tasmanian_post_install}\")")
install(CODE "message(\"information stored in: ${CMAKE_INSTALL_PREFIX}/config/Tasmanian.log\n\")")
