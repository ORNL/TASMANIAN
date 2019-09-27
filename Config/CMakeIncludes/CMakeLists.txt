########################################################################
# Install Message
# This is written in this CMakeLists.txt and add subdir is called last
# thus, the message will appear at the very end of the install process
#
# Do not put any other code here, this is for human readability only
########################################################################
set(Tasmanian_post_install "--------------------------------------------------------------------------------\n")
set(Tasmanian_post_install "${Tasmanian_post_install}Tasmanian Install Complete: Installed Version ${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}\n\n")
set(Tasmanian_post_install "${Tasmanian_post_install}executable:\n    tasgrid${CMAKE_EXECUTABLE_SUFFIX_CXX}\n")

foreach(_tsglibtype ${Tasmanian_libs_type})
    string(TOUPPER ${_tsglibtype} TSGLTYPE)
    set(Tasmanian_post_install "${Tasmanian_post_install}${_tsglibtype} libraries:\n")
    foreach(_tsgtarget Tasmanian_libsparsegrid_${_tsglibtype} Tasmanian_libdream_${_tsglibtype} Tasmanian_libfortran90_${_tsglibtype})
        if (TARGET ${_tsgtarget})
            get_target_property(_tsgname ${_tsgtarget} OUTPUT_NAME)
            set(Tasmanian_post_install "${Tasmanian_post_install}    ${CMAKE_${TSGLTYPE}_LIBRARY_PREFIX}${_tsgname}${CMAKE_${TSGLTYPE}_LIBRARY_SUFFIX}\n")
        endif()
        unset(_tsgname)
    endforeach()
    unset(_tsgtarget)
endforeach()
unset(TSGLTYPE)
unset(_tsglibtype)

set(Tasmanian_post_install "${Tasmanian_post_install}\nsee the examples:\n    ${CMAKE_INSTALL_PREFIX}/share/Tasmanian/examples/\n\n")

set(Tasmanian_post_install "${Tasmanian_post_install}bash environment setup:\n    source ${CMAKE_INSTALL_PREFIX}/config/TasmanianENVsetup.sh\n\n")

set(Tasmanian_post_install "${Tasmanian_post_install}cmake package config:\n    find_package(Tasmanian ${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}.${Tasmanian_VERSION_PATCH} PATHS \\\"${CMAKE_INSTALL_PREFIX}/lib/\\\")\n\n")
set(Tasmanian_post_install "${Tasmanian_post_install}    targets: Tasmanian::tasgrid    (executable)\n")

foreach(_tsglibtype ${Tasmanian_libs_type})
    set(Tasmanian_post_install "${Tasmanian_post_install}             Tasmanian::${_tsglibtype}     (${_tsglibtype} libraries)\n")
endforeach()
unset(_tsglibtype)

set(Tasmanian_post_install "${Tasmanian_post_install}             Tasmanian::Tasmanian  (alias to static/shared)\n")
set(Tasmanian_post_install "${Tasmanian_post_install}\n")

if (Tasmanian_ENABLE_PYTHON)
    set(Tasmanian_post_install "${Tasmanian_post_install}python module and the sym-link alternative:\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    ${Tasmanian_python_install_path}/Tasmanian.py\n    ${CMAKE_INSTALL_PREFIX}/share/Tasmanian/python/Tasmanian.py\n\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    import sys\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    sys.path.append(\\\"${CMAKE_INSTALL_PREFIX}/share/Tasmanian/python\\\")\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    import Tasmanian\n\n")
endif()

if (NOT "${Tasmanian_MATLAB_WORK_FOLDER}" STREQUAL "")
    set(Tasmanian_post_install "${Tasmanian_post_install}matlab scripts location:\n    ${CMAKE_INSTALL_PREFIX}/share/Tasmanian/matlab/\n\n")
    set(Tasmanian_post_install "${Tasmanian_post_install}    addpath('${CMAKE_INSTALL_PREFIX}/share/Tasmanian/matlab/')\n")
endif()
set(Tasmanian_post_install "${Tasmanian_post_install}--------------------------------------------------------------------------------\n")

install(CODE "message(\"${Tasmanian_post_install}\")")
install(CODE "file(WRITE \"${CMAKE_INSTALL_PREFIX}/share/Tasmanian/Tasmanian.log\" \"${Tasmanian_post_install}\")")
install(CODE "message(\"information stored in: ${CMAKE_INSTALL_PREFIX}/share/Tasmanian/Tasmanian.log\n\")")
