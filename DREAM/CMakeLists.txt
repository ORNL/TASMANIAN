########################################################################
# DREAM librareis and command line tool
########################################################################

add_library(Tasmanian_libdream  TasmanianDREAM.hpp
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

target_include_directories(Tasmanian_libdream PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>)
target_link_libraries(Tasmanian_libdream  Tasmanian_libsparsegrid)

set_target_properties(Tasmanian_libdream PROPERTIES OUTPUT_NAME "tasmaniandream"
                                                    SOVERSION ${Tasmanian_VERSION_MAJOR}
                                                    VERSION   ${PROJECT_VERSION})
Tasmanian_rpath_target(TARGET Tasmanian_libdream COMPONENTS SparseGrids)

install(TARGETS Tasmanian_libdream
        EXPORT  "${Tasmanian_export_name}"
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")


########################################################################
# add dreamtest and examples executables
########################################################################
add_executable(Tasmanian_dreamtest dreamtest_main.cpp
                                   tasdreamExternalTests.hpp
                                   tasdreamExternalTests.cpp)

set_target_properties(Tasmanian_dreamtest  PROPERTIES OUTPUT_NAME "dreamtest"
                                                      CXX_EXTENSIONS OFF)
Tasmanian_rpath_target(TARGET Tasmanian_dreamtest USE_CURRENT COMPONENTS SparseGrids)
target_link_libraries(Tasmanian_dreamtest  Tasmanian_libdream)


########################################################################
# Testing
########################################################################
add_test(DreamAnalytic dreamtest analytic)
add_test(DreamModel dreamtest posterior)
Tasmanian_set_test_properties(TESTS DreamAnalytic DreamModel)


########################################################################
# Windows specific support (DLL export/import directives and names)
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")

    foreach(_tsg_target libdream dreamtest)
        target_compile_definitions(Tasmanian_${_tsg_target} PRIVATE -D_SCL_SECURE_NO_WARNINGS)
    endforeach()
    unset(_tsg_target)

    if (BUILD_SHARED_LIBS)
        set_target_properties(Tasmanian_libdream PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)

        set(Tasmanian_MSVC_PATH_STRING "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Release;${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/Debug;$ENV{PATH}")
        string(REPLACE ";" "\\;" Tasmanian_MSVC_PATH_STRING "${Tasmanian_MSVC_PATH_STRING}")
        set_tests_properties(DreamAnalytic DreamModel PROPERTIES ENVIRONMENT "PATH=${Tasmanian_MSVC_PATH_STRING}")
    endif()

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
