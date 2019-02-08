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
                                                                          tsgDreamSampleGrid.hpp
                                                                          tsgDreamSamplePosterior.hpp
                                                                          tsgDreamSamplePosteriorGrid.hpp
                                                                          tsgDreamLikelihoodCore.hpp
                                                                          tsgDreamLikelyGaussian.hpp
                                                                          tsgDreamLikelyGaussian.cpp
                                                                          tsgDreamInternalBlas.hpp
                                                                          tsgDreamCoreRandom.hpp
                                                                          tsgDreamCorePDF.hpp
                                                                          tsgDreamEnumerates.hpp)

target_include_directories(${Tasmanian_libtdr_target_name} PUBLIC $<INSTALL_INTERFACE:include>)

target_include_directories(${Tasmanian_libtdr_target_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>)
target_include_directories(${Tasmanian_libtdr_target_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../configured/>)

set_target_properties(${Tasmanian_libtdr_target_name} PROPERTIES OUTPUT_NAME "tasmaniandream")

set_property(TARGET ${Tasmanian_libtdr_target_name} PROPERTY CXX_STANDARD 11)

########################################################################
# Option setup
########################################################################
if (Tasmanian_ENABLE_BLAS)
    target_link_libraries(${Tasmanian_libtdr_target_name} ${BLAS_LIBRARIES})
endif()

########################################################################
# MPI setup (experimental, DREAM distributed posterior only)
########################################################################
if (Tasmanian_ENABLE_MPI)
    target_link_libraries(${Tasmanian_libtdr_target_name} ${MPI_CXX_LIBRARIES})

    if (DEFINED MPI_CXX_INCLUDE_PATH)
        target_include_directories(${Tasmanian_libtdr_target_name} PUBLIC "${MPI_CXX_INCLUDE_PATH}")
    endif()

    if(DEFINED MPI_CXX_COMPILE_FLAGS)
        target_compile_options(${Tasmanian_libtdr_target_name} PUBLIC "${MPI_CXX_COMPILE_FLAGS}")
    endif()

    if(DEFINED MPI_CXX_LINK_FLAGS)
        set_target_properties(${Tasmanian_libtdr_target_name} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif()
endif()

install(TARGETS "${Tasmanian_libtdr_target_name}"
        EXPORT  "${Tasmanian_export_name}"
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")

endmacro()

########################################################################
# add dreamtest and examples executables
########################################################################
add_executable(Tasmanian_dreamtest dreamtest_main.cpp
                                   tasdreamExternalTests.hpp
                                   tasdreamExternalTests.cpp)

add_executable(Tasmanian_example_dream Examples/example_dream_01.cpp
                                       Examples/example_dream_02.cpp
                                       Examples/example_dream.cpp)

set_property(TARGET Tasmanian_dreamtest     PROPERTY CXX_STANDARD 11)
set_property(TARGET Tasmanian_example_dream PROPERTY CXX_STANDARD 11)

set_target_properties(Tasmanian_dreamtest     PROPERTIES OUTPUT_NAME "dreamtest")
set_target_properties(Tasmanian_example_dream PROPERTIES OUTPUT_NAME "example_dream")

########################################################################
# add the libraries
########################################################################
if (NOT "${Tasmanian_libs_type}" STREQUAL "STATIC_ONLY")
    Tasmanian_macro_add_libdream(SHARED)

    target_link_libraries(Tasmanian_libdream_shared Tasmanian_libsparsegrid_shared)
endif()

if (NOT "${Tasmanian_libs_type}" STREQUAL "SHARED_ONLY")
    Tasmanian_macro_add_libdream(STATIC)

    # if static libraries are available, executables link statically
    target_link_libraries(Tasmanian_libdream_static  Tasmanian_libsparsegrid_static)
    target_link_libraries(Tasmanian_dreamtest        Tasmanian_libdream_static)
    target_link_libraries(Tasmanian_example_dream    Tasmanian_libdream_static)
else() # if there are not static libs, use shared libs for the executables
    target_link_libraries(Tasmanian_dreamtest      Tasmanian_libdream_shared)
    target_link_libraries(Tasmanian_example_dream  Tasmanian_libdream_shared)
endif()

########################################################################
# Testing
########################################################################
add_test(DREAM_Analytic dreamtest analytic)
add_test(DREAM_Model dreamtest posterior)
if (Tasmanian_TESTS_OMP_NUM_THREADS GREATER 0)
    set_tests_properties(DREAM_Analytic DREAM_Model
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
        target_compile_definitions(${Tasmanian_dream_target} PRIVATE -D_TASMANIAN_WINDOWS_) # overwrittes gettime()
        target_compile_definitions(${Tasmanian_dream_target} PRIVATE -D_SCL_SECURE_NO_WARNINGS) # suppresses warnings regarding pointers to the middle of an array
        target_compile_definitions(${Tasmanian_dream_target} PUBLIC  -D_USE_MATH_DEFINES) # needed to include M_PI constant
    endmacro()

    if (TARGET Tasmanian_libdream_shared)
        set_target_properties(Tasmanian_libdream_shared PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
        Tasmanian_macro_dream_windows_defines(Tasmanian_libdream_shared)
    endif()

    if (TARGET Tasmanian_libdream_static)
        Tasmanian_macro_dream_windows_defines(Tasmanian_libdream_static)
        set_target_properties(Tasmanian_libdream_static PROPERTIES OUTPUT_NAME "tasmaniandream_static")
    endif()

    if ("${Tasmanian_libs_type}" STREQUAL "SHARED_ONLY")
        set(Tasmanian_MSVC_PATH_STRING "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Release;${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Debug;$ENV{PATH}")
        string(REPLACE ";" "\\;" Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING}")
        set_tests_properties(DREAM_Analytic DREAM_Model PROPERTIES ENVIRONMENT "PATH=${Tasmanian_MSVC_PATH_STRING}")
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

#install(TARGETS Tasmanian_tasdream
#        EXPORT "${CMAKE_PROJECT_NAME}-exports"
#        RUNTIME DESTINATION "bin"
#        LIBRARY DESTINATION "lib"
#        ARCHIVE DESTINATION "lib")
