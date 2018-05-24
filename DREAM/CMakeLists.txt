cmake_minimum_required(VERSION 3.5)

########################################################################
# DREAM librareis and command line tool
########################################################################

macro(Tasmanian_macro_add_libdream Tasmanian_shared_or_static)

if (${Tasmanian_shared_or_static} STREQUAL SHARED)
    set(Tasmanian_libtdr_target_name "Tasmanian_libdream_shared")
else()
    set(Tasmanian_libtdr_target_name "Tasmanian_libdream_static")
endif()

add_library(${Tasmanian_libtdr_target_name} ${Tasmanian_shared_or_static} ${CMAKE_CURRENT_SOURCE_DIR}/TasmanianDREAM.hpp
                                                                          ${CMAKE_CURRENT_SOURCE_DIR}/TasmanianDREAM.cpp
                                                                          ${CMAKE_CURRENT_SOURCE_DIR}/tdrEnumerates.hpp
                                                                          ${CMAKE_CURRENT_SOURCE_DIR}/tdrCorePDF.hpp
                                                                          ${CMAKE_CURRENT_SOURCE_DIR}/tdrCorePDF.cpp)

target_include_directories(${Tasmanian_libtdr_target_name} PUBLIC $<INSTALL_INTERFACE:include>)

target_include_directories(${Tasmanian_libtdr_target_name} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/DREAM/>)
target_include_directories(${Tasmanian_libtdr_target_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/configured/>)

set_target_properties(${Tasmanian_libtdr_target_name} PROPERTIES OUTPUT_NAME "tasmaniandream")

# setup the name, in MacOSX and Linux static libraries use .a extension and shared libs use .dylib and .so
# in Windows, shared .dll libs need another file for linking which shares .lib extension with the static libs
# thus on Windows, we specifically add _static to the name of the lib
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    target_compile_definitions(${Tasmanian_libtdr_target_name} -D_TASMANIAN_WINDOWS_)
    target_compile_definitions(${Tasmanian_libtdr_target_name} -D_SCL_SECURE_NO_WARNINGS)
    target_compile_definitions(${Tasmanian_libtdr_target_name} -D_USE_MATH_DEFINES)
    if (${Tasmanian_shared_or_static} STREQUAL STATIC)
        set_target_properties(${Tasmanian_libtdr_target_name} PROPERTIES OUTPUT_NAME "tasmaniandream_static")
    endif()
endif()

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
        target_include_directories(${Tasmanian_libtdr_target_name} PUBLIC ${MPI_CXX_INCLUDE_PATH})
    endif()

    if(DEFINED MPI_CXX_COMPILE_FLAGS)
        target_compile_options(${Tasmanian_libtdr_target_name} PUBLIC ${MPI_CXX_COMPILE_FLAGS})
    endif()

    if(DEFINED MPI_CXX_LINK_FLAGS)
        set_target_properties(${Tasmanian_libtdr_target_name} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif()
endif()

endmacro()

########################################################################
# add tasdream and examples executables
########################################################################
add_executable(Tasmanian_tasdream ${CMAKE_CURRENT_SOURCE_DIR}/tasdream_main.cpp
                                  ${CMAKE_CURRENT_SOURCE_DIR}/tasdreamExternalTests.hpp
                                  ${CMAKE_CURRENT_SOURCE_DIR}/tasdreamExternalTests.cpp
                                  ${CMAKE_CURRENT_SOURCE_DIR}/tasdreamTestPDFs.hpp
                                  ${CMAKE_CURRENT_SOURCE_DIR}/tasdreamTestPDFs.cpp
                                  ${CMAKE_CURRENT_SOURCE_DIR}/tasdreamBenchmark.hpp
                                  ${CMAKE_CURRENT_SOURCE_DIR}/tasdreamBenchmark.cpp)
add_executable(Tasmanian_example_dream "${CMAKE_CURRENT_SOURCE_DIR}/Examples/example_dream.cpp")

list(APPEND Tasmanian_install_dream_targets Tasmanian_tasdream)

set_target_properties(Tasmanian_tasdream PROPERTIES OUTPUT_NAME "tasdream")
set_target_properties(Tasmanian_example_dream PROPERTIES OUTPUT_NAME "example_dream")

########################################################################
# add the libraries
########################################################################
if (Tasmanian_SHARED_LIBRARY)
    Tasmanian_macro_add_libdream(SHARED)

    target_link_libraries(Tasmanian_libdream_shared Tasmanian_libsparsegrid_shared)

    list(APPEND Tasmanian_install_dream_targets Tasmanian_libdream_shared)
endif()

if (Tasmanian_STATIC_LIBRARY)
    Tasmanian_macro_add_libdream(STATIC)

    target_link_libraries(Tasmanian_libdream_static Tasmanian_libsparsegrid_static)
    target_link_libraries(Tasmanian_tasdream        Tasmanian_libdream_static)
    target_link_libraries(Tasmanian_example_dream   Tasmanian_libdream_static)

    list(APPEND Tasmanian_install_dream_targets Tasmanian_libdream_static)
else()
    target_link_libraries(Tasmanian_tasdream      Tasmanian_libdream_shared)
    target_link_libraries(Tasmanian_example_dream Tasmanian_libdream_shared)
endif()

########################################################################
# Windows specific support (DLL export/import directives and names)
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    if (Tasmanian_SHARED_LIBRARY)
        target_compile_definitions(Tasmanian_libdeam_shared PUBLIC -DTSG_DLL)
    else()
        target_compile_definitions(Tasmanian_tasdream PUBLIC -DTSG_DYNAMIC)
    endif()

    if (Tasmanian_STATIC_LIBRARY)
        target_compile_definitions(Tasmanian_libdream_static PUBLIC -DTSG_STATIC)
    endif()

    target_compile_definitions(Tasmanian_tasdream -D_TASMANIAN_WINDOWS_) # overwrittes gettime()
    target_compile_definitions(Tasmanian_tasdream -D_SCL_SECURE_NO_WARNINGS) # suppresses warnings regarding pointers to the middle of an array
    target_compile_definitions(Tasmanian_tasdream -D_USE_MATH_DEFINES) # needed to include M_PI constant (lots of targets need this, will figure it out)
endif()

########################################################################
# Testing
########################################################################
add_test(DREAM_Analytic ${CMAKE_CURRENT_BINARY_DIR}/tasdream -test analytic)
add_test(DREAM_Model    ${CMAKE_CURRENT_BINARY_DIR}/tasdream -test model)
set(Tasmanian_tests_list_dream DREAM_Analytic DREAM_Model)

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
        DESTINATION "examples"
        FILES_MATCHING PATTERN "*.cpp"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(TARGETS ${Tasmanian_install_dream_targets}
        EXPORT Tasmanian_exports
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")
