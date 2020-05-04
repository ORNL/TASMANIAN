########################################################################
# DREAM librareis and command line tool
########################################################################

macro(Tasmanian_macro_add_libdream Tasmanian_shared_or_static)

if (${Tasmanian_shared_or_static} STREQUAL SHARED)
    set(Tasmanian_libtdr_target_name "Tasmanian_libdream_shared")
else()
    set(Tasmanian_libtdr_target_name "Tasmanian_libdream_static")
endif()

add_library(${Tasmanian_libtdr_target_name} ${Tasmanian_shared_or_static} TasmanianDREAM.hpp
                                                                          tsgDreamState.hpp
                                                                          tsgDreamState.cpp
                                                                          tsgDreamSample.hpp
                                                                          tsgDreamSampleWrapC.cpp
                                                                          tsgDreamLikelihoodCore.hpp
                                                                          tsgDreamLikelyGaussian.hpp
                                                                          tsgDreamLikelyGaussian.cpp
                                                                          tsgDreamCoreRandom.hpp
                                                                          tsgDreamCorePDF.hpp
                                                                          tsgDreamEnumerates.hpp)

target_include_directories(${Tasmanian_libtdr_target_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>)

set_target_properties(${Tasmanian_libtdr_target_name} PROPERTIES OUTPUT_NAME "tasmaniandream"
                                                                 INSTALL_RPATH "${Tasmanian_rpath}"
                                                                 SOVERSION ${Tasmanian_VERSION_MAJOR}
                                                                 VERSION   ${PROJECT_VERSION})

install(TARGETS "${Tasmanian_libtdr_target_name}"
        EXPORT  "${Tasmanian_export_name}"
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")

unset(Tasmanian_libtdr_target_name)
endmacro()

########################################################################
# add dreamtest and examples executables
########################################################################
add_executable(Tasmanian_dreamtest dreamtest_main.cpp
                                   tasdreamExternalTests.hpp
                                   tasdreamExternalTests.cpp)

set_target_properties(Tasmanian_dreamtest  PROPERTIES OUTPUT_NAME "dreamtest" CXX_EXTENSIONS OFF)

########################################################################
# add the libraries
########################################################################
foreach(_tsglibtype ${Tasmanian_libs_type})
    string(TOUPPER ${_tsglibtype} Tasmanian_ltype)
    Tasmanian_macro_add_libdream(${Tasmanian_ltype})

    target_link_libraries(Tasmanian_libdream_${_tsglibtype} Tasmanian_libsparsegrid_${_tsglibtype})
endforeach()
unset(_tsglibtype)
unset(Tasmanian_ltype)

add_library(Tasmanian_libdream ALIAS Tasmanian_libdream_${Tasmanian_lib_default})

target_link_libraries(Tasmanian_dreamtest  Tasmanian_libdream)


########################################################################
# Testing
########################################################################
add_test(DreamAnalytic dreamtest analytic)
add_test(DreamModel dreamtest posterior)
if (Tasmanian_TESTS_OMP_NUM_THREADS GREATER 0)
    set_tests_properties(DreamAnalytic DreamModel
        PROPERTIES
        PROCESSORS "${Tasmanian_TESTS_OMP_NUM_THREADS}"
        ENVIRONMENT "OMP_NUM_THREADS=${Tasmanian_TESTS_OMP_NUM_THREADS}")
endif()

########################################################################
# Windows specific support (DLL export/import directives and names)
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    # setup the name, in MacOSX and Linux static libraries use .a extension and shared libs use .dylib and .so
    # in Windows, shared .dll libs need another file for linking which shares .lib extension with the static libs
    # thus on Windows, we specifically add _static to the name of the lib
    # also, Windows requires export/import directives (double-check those with cmake defaults)
    macro(Tasmanian_macro_dream_windows_defines Tasmanian_dream_target)
        target_compile_definitions(${Tasmanian_dream_target} PRIVATE -D_SCL_SECURE_NO_WARNINGS) # suppresses warnings regarding pointers to the middle of an array
    endmacro()

    if (TARGET Tasmanian_libdream_shared)
        set_target_properties(Tasmanian_libdream_shared PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
        Tasmanian_macro_dream_windows_defines(Tasmanian_libdream_shared)
    endif()

    if (TARGET Tasmanian_libdream_static)
        Tasmanian_macro_dream_windows_defines(Tasmanian_libdream_static)
        set_target_properties(Tasmanian_libdream_static PROPERTIES OUTPUT_NAME "tasmaniandream_static")
    endif()

    if (NOT "static" IN_LIST Tasmanian_libs_type)
        set(Tasmanian_MSVC_PATH_STRING "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Release;${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Debug;$ENV{PATH}")
        string(REPLACE ";" "\\;" Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING}")
        set_tests_properties(DreamAnalytic DreamModel PROPERTIES ENVIRONMENT "PATH=${Tasmanian_MSVC_PATH_STRING}")
    endif()

    Tasmanian_macro_dream_windows_defines(Tasmanian_dreamtest)
endif()

########################################################################
# Install headers and config files
########################################################################
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "*.windows.*" EXCLUDE
        PATTERN "*.in.*" EXCLUDE
        PATTERN "Examples" EXCLUDE)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Examples/"
        DESTINATION "share/Tasmanian/examples/"
        FILES_MATCHING PATTERN "*.cpp"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
