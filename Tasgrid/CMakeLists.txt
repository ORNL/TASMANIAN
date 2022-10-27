########################################################################
# tasgrid command line tool
########################################################################

########################################################################
# add the tasgrid executable
########################################################################
add_executable(Tasmanian_tasgrid tasgrid_main.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/gridtestCLICommon.hpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/gridtestExternalTests.hpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/gridtestExternalTests.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/gridtestTestHelpers.hpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/gridtestTestFunctions.hpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/gridtestTestFunctions.cpp
                                 tasgridWrapper.hpp
                                 tasgridWrapper.cpp)

target_link_libraries(Tasmanian_tasgrid Tasmanian_addons)

set_target_properties(Tasmanian_tasgrid PROPERTIES OUTPUT_NAME "tasgrid" CXX_EXTENSIONS OFF)
Tasmanian_rpath_target(TARGET Tasmanian_tasgrid COMPONENTS SparseGrids DREAM)
if (Tasmanian_ENABLE_OPENMP)
    # the OpenMP libraries are carried transitively from sparse grids library
    target_compile_options(Tasmanian_tasgrid PRIVATE ${OpenMP_CXX_FLAGS})
endif()
if (Tasmanian_ENABLE_HIP)
    target_link_libraries(Tasmanian_tasgrid hip::host)
endif()

########################################################################
# Windows specific support (DLL export/import directives and names)
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")

    # suppresses warnings regarding pointers to the middle of an array
    target_compile_definitions(Tasmanian_tasgrid PRIVATE -D_SCL_SECURE_NO_WARNINGS)
    # needed to prevent crash on using STL vector iterators
    target_compile_definitions(Tasmanian_tasgrid PUBLIC  -D_HAS_ITERATOR_DEBUGGING=0)

endif()

########################################################################
# Install binaries and shared objects
########################################################################

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/../configured/tasgridLogs.hpp"
        DESTINATION include
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(TARGETS Tasmanian_tasgrid
        EXPORT "${Tasmanian_export_name}"
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")
