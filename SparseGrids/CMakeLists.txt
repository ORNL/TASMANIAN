########################################################################
# Sparse Grid librareis and command line tool
########################################################################

########################################################################
# Define source files
# in order to avoid GLOB, list all source files for add_library()
########################################################################
set(Tasmanian_source_libsparsegrid
        TasmanianSparseGrid.hpp
        TasmanianSparseGrid.cpp
        TasmanianSparseGridWrapC.cpp
        tsgAcceleratedHandles.hpp
        tsgAcceleratedDataStructures.hpp
        tsgAcceleratedDataStructures.cpp
        tsgCacheLagrange.hpp
        tsgCoreOneDimensional.hpp
        tsgCoreOneDimensional.cpp
        tsgDConstructGridGlobal.hpp
        tsgDConstructGridGlobal.cpp
        tsgCudaLoadStructures.hpp
        tsgEnumerates.hpp
        tsgGridCore.hpp
        tsgGridGlobal.hpp
        tsgGridGlobal.cpp
        tsgGridWavelet.hpp
        tsgGridWavelet.cpp
        tsgHardCodedTabulatedRules.hpp
        tsgHardCodedTabulatedRules.cpp
        tsgGridLocalPolynomial.hpp
        tsgGridLocalPolynomial.cpp
        tsgGridSequence.hpp
        tsgGridSequence.cpp
        tsgGridFourier.hpp
        tsgGridFourier.cpp
        tsgIndexManipulator.hpp
        tsgIndexManipulator.cpp
        tsgHierarchyManipulator.hpp
        tsgHierarchyManipulator.cpp
        tsgIOHelpers.hpp
        tsgIndexSets.hpp
        tsgIndexSets.cpp
        tsgLinearSolvers.hpp
        tsgLinearSolvers.cpp
        tsgOneDimensionalWrapper.hpp
        tsgRuleLocalPolynomial.hpp
        tsgRuleWavelet.hpp
        tsgRuleWavelet.cpp
        tsgSequenceOptimizer.hpp
        tsgSequenceOptimizer.cpp
        tsgMathUtils.hpp
        tsgUtils.hpp)

# source files specific to cuda, used for both static and shared libs
set(Tasmanian_source_libsparsegrid_cuda
       ${CMAKE_CURRENT_SOURCE_DIR}/../InterfaceTPL/tsgCudaWrappers.cpp
       tsgCudaKernels.cu)


########################################################################
# add the tasgrid and examples executables
########################################################################
add_executable(Tasmanian_tasgrid tasgrid_main.cpp
                                 TasmanianSparseGrid.hpp
                                 tasgridCLICommon.hpp
                                 tasgridExternalTests.hpp
                                 tasgridExternalTests.cpp
                                 tasgridTestHelpers.hpp
                                 tasgridTestFunctions.hpp
                                 tasgridTestFunctions.cpp
                                 tasgridWrapper.hpp
                                 tasgridWrapper.cpp)
add_executable(Tasmanian_gridtest gridtest_main.cpp
                                  TasmanianSparseGrid.hpp
                                  tasgridCLICommon.hpp
                                  tasgridExternalTests.hpp
                                  tasgridExternalTests.cpp
                                  tasgridTestHelpers.hpp
                                  tasgridTestFunctions.hpp
                                  tasgridTestFunctions.cpp
                                  tasgridUnitTests.hpp
                                  tasgridUnitTests.cpp
                                  tasgridTestInterfaceC.cpp)
add_executable(Tasmanian_benchmarksgrid  Benchmarks/benchCommon.hpp
                                         tasgridCLICommon.hpp
                                         Benchmarks/benchMakeGrid.hpp
                                         Benchmarks/benchLoadNeeded.hpp
                                         Benchmarks/benchEvaluate.hpp
                                         Benchmarks/bench_main.cpp)

if (Tasmanian_ENABLE_CUDA)
    # CMake does not automatically add the CUDA include path, any source file that adds <cuda.h> directly or transitively
    # must be labeled as "CUDA" source files and compiled with the nvcc compiler.
    # Beware that the nvcc compiler seems to struggle with some STL algorithms, transitive inclusion of <cuda.h> must be minimized.
    set_source_files_properties(${Tasmanian_source_libsparsegrid_cuda} PROPERTIES LANGUAGE CUDA)
    add_library(Tasmanian_libsparsegrid ${Tasmanian_source_libsparsegrid} ${Tasmanian_source_libsparsegrid_cuda})
    set_property(TARGET Tasmanian_libsparsegrid PROPERTY CUDA_STANDARD 11)
elseif(Tasmanian_ENABLE_HIP)
    add_library(Tasmanian_libsparsegrid ${Tasmanian_source_libsparsegrid}
                ${CMAKE_CURRENT_SOURCE_DIR}/../InterfaceTPL/tsgHipWrappers.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/tsgHipKernels.hip.cpp)
elseif(Tasmanian_ENABLE_DPCPP)
    add_library(Tasmanian_libsparsegrid ${Tasmanian_source_libsparsegrid}
                ${CMAKE_CURRENT_SOURCE_DIR}/../InterfaceTPL/tsgDpcppWrappers.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/tsgDpcppKernels.cpp)
else()
    add_library(Tasmanian_libsparsegrid ${Tasmanian_source_libsparsegrid}
                ${CMAKE_CURRENT_SOURCE_DIR}/../InterfaceTPL/tsgGpuNull.cpp)
endif()

target_link_libraries(Tasmanian_libsparsegrid  Tasmanian_dependencies)

target_include_directories(Tasmanian_libsparsegrid PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>)
# the following are needed for Exotic quadrature only
target_include_directories(Tasmanian_libsparsegrid PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../Addons>)
target_include_directories(Tasmanian_libsparsegrid PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../DREAM>)

set_target_properties(Tasmanian_libsparsegrid PROPERTIES OUTPUT_NAME "tasmaniansparsegrid"
                                                         SOVERSION ${Tasmanian_VERSION_MAJOR}
                                                         VERSION   ${PROJECT_VERSION})
Tasmanian_rpath_target(TARGET Tasmanian_libsparsegrid)

if (Tasmanian_ENABLE_OPENMP)
    # the nvcc compiler does nor recognize OpenMP, add the flag only to non-CUDA source files
    target_compile_options(Tasmanian_libsparsegrid PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>)
endif()

install(TARGETS Tasmanian_libsparsegrid
        EXPORT  "${Tasmanian_export_name}"
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")


foreach(_tsgtarget tasgrid gridtest benchmarksgrid)
    target_link_libraries(Tasmanian_${_tsgtarget} Tasmanian_libsparsegrid)
    set_target_properties(Tasmanian_${_tsgtarget} PROPERTIES OUTPUT_NAME "${_tsgtarget}" CXX_EXTENSIONS OFF)
    Tasmanian_rpath_target(TARGET Tasmanian_${_tsgtarget} USE_CURRENT)
    if (Tasmanian_ENABLE_OPENMP)
        # the OpenMP libraries are carried transitively from sparse grids library
        target_compile_options(Tasmanian_${_tsgtarget} PRIVATE ${OpenMP_CXX_FLAGS})
    endif()
endforeach()
unset(_tsgtarget)
#set_target_properties(Tasmanian_tasgrid PROPERTIES INSTALL_RPATH "${Tasmanian_rpath}")

# data file, needed for testing and reference about custom rule definitions
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaussPattersonRule.table"  "${CMAKE_CURRENT_BINARY_DIR}/GaussPattersonRule.table" COPYONLY)


########################################################################
# Windows specific support (DLL export/import directives and names)
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")

    foreach(_tsg_target libsparsegrid tasgrid gridtest)
        # suppresses warnings regarding pointers to the middle of an array
        target_compile_definitions(Tasmanian_${_tsg_target} PRIVATE -D_SCL_SECURE_NO_WARNINGS)
        # needed to prevent crash on using STL vector iterators
        target_compile_definitions(Tasmanian_${_tsg_target} PUBLIC  -D_HAS_ITERATOR_DEBUGGING=0)
    endforeach()
    unset(_tsg_target)

    if (BUILD_SHARED_LIBS)
        set_target_properties(Tasmanian_libsparsegrid PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    endif()

endif()


########################################################################
# Testing
########################################################################
add_test(SparseGridsAcceleration gridtest acceleration -gpuid ${Tasmanian_TESTS_GPU_ID})
add_test(SparseGridsDomain       gridtest domain)
add_test(SparseGridsRefinement   gridtest refinement)
add_test(SparseGridsGlobal       gridtest global)
add_test(SparseGridsLocal        gridtest local)
add_test(SparseGridsWavelet      gridtest wavelet)
add_test(SparseGridsFourier      gridtest fourier)
add_test(SparseGridsExceptions   gridtest errors)
add_test(SparseGridsAPI          gridtest api)
add_test(SparseGridsC            gridtest c)
Tasmanian_set_test_properties(TESTS SparseGridsAcceleration SparseGridsDomain SparseGridsRefinement SparseGridsGlobal SparseGridsLocal SparseGridsWavelet SparseGridsFourier SparseGridsExceptions SparseGridsAPI SparseGridsC)

if (Tasmanian_ENABLE_BLAS)
    add_test(SparseGridsLAPACK gridtest lapack)
    Tasmanian_set_test_properties(TESTS SparseGridsLAPACK)
endif()

########################################################################
# Install headers and config files
########################################################################
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "*.windows.*" EXCLUDE
        PATTERN "Examples" EXCLUDE
        PATTERN "Benchmarks" EXCLUDE
        PATTERN "tsgHiddenExternals.hpp" EXCLUDE)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSparseGrid.h"
        DESTINATION include
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/../configured/tasgridLogs.hpp"
        DESTINATION include
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Examples/"
        DESTINATION "share/Tasmanian/examples/"
        FILES_MATCHING PATTERN "*.cpp"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/GaussPattersonRule.table"
        DESTINATION "share/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(TARGETS Tasmanian_tasgrid
        EXPORT "${Tasmanian_export_name}"
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")
