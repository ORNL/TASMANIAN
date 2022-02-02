########################################################################
# Optimization Libraries
########################################################################

add_library(Tasmanian_optimization  TasmanianOptimization.hpp
                                    tsgParticleSwarm.hpp
                                    tsgParticleSwarm.cpp)

target_include_directories(Tasmanian_optimization PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>)
target_link_libraries(Tasmanian_optimization  Tasmanian_libsparsegrid)

set_target_properties(Tasmanian_optimization PROPERTIES OUTPUT_NAME "tasmanianoptimization"
                                                    SOVERSION ${Tasmanian_VERSION_MAJOR}
                                                    VERSION   ${PROJECT_VERSION})
Tasmanian_rpath_target(TARGET Tasmanian_optimization COMPONENTS SparseGrids)

install(TARGETS Tasmanian_optimization
        EXPORT  "${Tasmanian_export_name}"
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")

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
