########################################################################
# Optimization Libraries
########################################################################

add_library(Tasmanian_liboptimization  TasmanianOptimization.hpp
                                       tsgParticleSwarmSolver.hpp
                                       tsgParticleSwarmSolver.cpp
                                       tsgSolver.hpp
                                       tsgSolver.cpp
                                       tsgOptimizationEnumerates.hpp)

target_include_directories(Tasmanian_liboptimization PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>)
target_link_libraries(Tasmanian_liboptimization  Tasmanian_libsparsegrid)

set_target_properties(Tasmanian_liboptimization PROPERTIES OUTPUT_NAME "tasmanianoptimization"
                                                    SOVERSION ${Tasmanian_VERSION_MAJOR}
                                                    VERSION   ${PROJECT_VERSION})
Tasmanian_rpath_target(TARGET Tasmanian_liboptimization COMPONENTS SparseGrids)

install(TARGETS Tasmanian_liboptimization
        EXPORT  "${Tasmanian_export_name}"
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")

# Tester executable
add_executable(Tasmanian_optimizationtest testOptimization.cpp)
set_target_properties(Tasmanian_optimizationtest  PROPERTIES OUTPUT_NAME "optimizationtest" CXX_EXTENSIONS OFF)
Tasmanian_rpath_target(TARGET Tasmanian_optimizationtest USE_CURRENT COMPONENTS SparseGrids)
target_link_libraries(Tasmanian_optimizationtest  Tasmanian_liboptimization)


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
