########################################################################
# Install Message
# This is written in this CMakeLists.txt and add subdir is called last
# thus, the message will appear at the very end of the install process
#
# Do not put any other code here, this is for human readability only
########################################################################
set(Tasmanian_post_install "--------------------------------------------------------------------------------\n\n")
set(Tasmanian_post_install "${Tasmanian_post_install}Tasmanian Install Complete: Installed Version ${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}${Tasmanian_version_comment}\n\n")
set(Tasmanian_post_install "${Tasmanian_post_install}executable:\n    tasgrid${CMAKE_EXECUTABLE_SUFFIX_CXX}\n")

if (BUILD_SHARED_LIBS)
    set(TSGLTYPE "SHARED")
else()
    set(TSGLTYPE "STATIC")
endif()

set(Tasmanian_post_install "${Tasmanian_post_install}libraries:\n")
foreach(_tsgtarget Tasmanian_libsparsegrid Tasmanian_libdream Tasmanian_libfortran90 Tasmanian_libfortran03)
    if (TARGET ${_tsgtarget})
        get_target_property(_tsgname ${_tsgtarget} OUTPUT_NAME)
        set(Tasmanian_post_install "${Tasmanian_post_install}    ${CMAKE_${TSGLTYPE}_LIBRARY_PREFIX}${_tsgname}${CMAKE_${TSGLTYPE}_LIBRARY_SUFFIX}\n")
    endif()
    unset(_tsgname)
endforeach()
unset(TSGLTYPE)
unset(_tsgtarget)
unset(Tasmanian_fortran_target)

if (Tasmanian_ENABLE_FORTRAN)
    set(Tasmanian_post_install "${Tasmanian_post_install}Fortran module:\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    tasmanian.mod\n")
    if (Tasmanian_ENABLE_MPI)
        set(Tasmanian_post_install "${Tasmanian_post_install}    tasmanian_mpi.mod\n")
    endif()
    set(Tasmanian_post_install "${Tasmanian_post_install}    tasmaniansg.mod\n")
endif()

set(Tasmanian_post_install "${Tasmanian_post_install}\nsee the examples:\n    ${Tasmanian_final_install_path}/share/Tasmanian/examples/\n\n")

set(Tasmanian_post_install "${Tasmanian_post_install}bash environment setup:\n    source ${Tasmanian_final_install_path}/share/Tasmanian/TasmanianENVsetup.sh\n\n")

set(Tasmanian_post_install "${Tasmanian_post_install}cmake package config:\n    find_package(Tasmanian ${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}.${Tasmanian_VERSION_PATCH} PATHS \\\"${Tasmanian_final_install_path}/lib/\\\")\n\n")
set(Tasmanian_post_install "${Tasmanian_post_install}    CMake components:\n")
set(Tasmanian_post_install "${Tasmanian_post_install}             ${Tasmanian_components}\n")
set(Tasmanian_post_install "${Tasmanian_post_install}    CXX targets:\n")
set(Tasmanian_post_install "${Tasmanian_post_install}             Tasmanian::tasgrid      (executable)\n")
set(Tasmanian_post_install "${Tasmanian_post_install}             Tasmanian::Tasmanian    (CXX libraries)\n\n")

if (Tasmanian_ENABLE_FORTRAN)
    set(Tasmanian_post_install "${Tasmanian_post_install}    Fortran targets:\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}             Tasmanian::Fortran      (Fortran libraries)\n\n")
endif()

set(Tasmanian_post_install "${Tasmanian_post_install}    see also:\n        ${Tasmanian_final_install_path}/share/Tasmanian/examples/CMakeLists.txt\n")
set(Tasmanian_post_install "${Tasmanian_post_install}\n")

if (Tasmanian_ENABLE_PYTHON AND ("${Tasmanian_final_install_path}" STREQUAL "${CMAKE_INSTALL_PREFIX}"))
    set(Tasmanian_post_install "${Tasmanian_post_install}python module path and alternatives:\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    ${Tasmanian_python_install_path}/Tasmanian.py\n    ${CMAKE_INSTALL_PREFIX}/share/Tasmanian/python/Tasmanian.py\n\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    import sys\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    sys.path.append(\\\"${CMAKE_INSTALL_PREFIX}/share/Tasmanian/python\\\")\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    import Tasmanian\n\n")
endif()

if (NOT "${Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
    set(Tasmanian_post_install "${Tasmanian_post_install}matlab scripts location:\n    ${Tasmanian_final_install_path}/share/Tasmanian/matlab/\n\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    addpath('${Tasmanian_final_install_path}/share/Tasmanian/matlab/')\n")
endif()
set(Tasmanian_post_install "${Tasmanian_post_install}--------------------------------------------------------------------------------\n")

string(REPLACE "\\\"" "\"" _tsg_message "${Tasmanian_post_install}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/Tasmanian.log" "${_tsg_message}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Tasmanian.log"
        DESTINATION "share/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
unset(_tsg_message)

install(CODE "message(\"${Tasmanian_post_install}\")")
install(CODE "message(\"information stored in: ${CMAKE_INSTALL_PREFIX}/share/Tasmanian/Tasmanian.log\n\")")
