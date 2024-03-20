########################################################################
# Sparse Grid libraries and gridtest command line tool
########################################################################

########################################################################
# Define source files
# in order to avoid GLOB, list all source files for add_library()
########################################################################
# source files specific to cuda/hip, used for both static and shared libs
set(Tasmanian_source_libsparsegrid_cuda
       ${CMAKE_CURRENT_SOURCE_DIR}/../InterfaceTPL/tsgCudaWrappers.cpp
       tsgCudaKernels.cu)
set(Tasmanian_source_libsparsegrid_hip
       ${CMAKE_CURRENT_SOURCE_DIR}/../InterfaceTPL/tsgHipWrappers.cpp
       tsgHipKernels.hip.cpp)

set(Tasmanian_sparse_grid_null_backend ON)
if (Tasmanian_ENABLE_CUDA OR Tasmanian_ENABLE_HIP OR Tasmanian_ENABLE_DPCPP)
    set(Tasmanian_sparse_grid_null_backend OFF)
endif()

add_library(Tasmanian_libsparsegrid
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
        tsgUtils.hpp
        $<$<BOOL:${Tasmanian_sparse_grid_null_backend}>:${CMAKE_CURRENT_SOURCE_DIR}/../InterfaceTPL/tsgGpuNull.cpp>
        $<$<BOOL:${Tasmanian_ENABLE_CUDA}>:${Tasmanian_source_libsparsegrid_cuda}>
        $<$<BOOL:${Tasmanian_ENABLE_HIP}>:${Tasmanian_source_libsparsegrid_hip}>
        $<$<BOOL:${Tasmanian_ENABLE_DPCPP}>:${CMAKE_CURRENT_SOURCE_DIR}/../InterfaceTPL/tsgDpcppWrappers.cpp>
        $<$<BOOL:${Tasmanian_ENABLE_DPCPP}>:${CMAKE_CURRENT_SOURCE_DIR}/tsgDpcppKernels.cpp>
        )

########################################################################
# add the gridtest and examples executables
########################################################################
add_executable(Tasmanian_gridtest gridtest_main.cpp
                                  TasmanianSparseGrid.hpp
                                  gridtestCLICommon.hpp
                                  gridtestExternalTests.hpp
                                  gridtestExternalTests.cpp
                                  gridtestTestHelpers.hpp
                                  gridtestTestFunctions.hpp
                                  gridtestTestFunctions.cpp
                                  gridtestUnitTests.hpp
                                  gridtestUnitTests.cpp
                                  gridtestTestInterfaceC.cpp)

add_executable(Tasmanian_benchmarksgrid  Benchmarks/benchCommon.hpp
                                         gridtestCLICommon.hpp
                                         Benchmarks/benchMakeGrid.hpp
                                         Benchmarks/benchLoadNeeded.hpp
                                         Benchmarks/benchDifferentiate.hpp
                                         Benchmarks/benchEvaluate.hpp
                                         Benchmarks/benchRefine.hpp
                                         Benchmarks/bench_main.cpp)

if (Tasmanian_ENABLE_CUDA)
    # CMake does not automatically add the CUDA include path, any source file that adds <cuda.h> directly or transitively
    # must be labeled as "CUDA" source files and compiled with the nvcc compiler.
    # Beware that the nvcc compiler seems to struggle with some STL algorithms, transitive inclusion of <cuda.h> must be minimized.
    set_source_files_properties(${Tasmanian_source_libsparsegrid_cuda} PROPERTIES LANGUAGE CUDA)
    set_property(TARGET Tasmanian_libsparsegrid PROPERTY CUDA_STANDARD 11)
elseif(Tasmanian_ENABLE_HIP)
    set_source_files_properties(${Tasmanian_source_libsparsegrid_hip} PROPERTIES LANGUAGE HIP)
elseif(Tasmanian_ENABLE_DPCPP)
    target_compile_options(Tasmanian_libsparsegrid PRIVATE -fsycl)
    target_compile_options(Tasmanian_gridtest      PRIVATE -fsycl)
    target_link_options(Tasmanian_libsparsegrid    PRIVATE -fsycl)
    target_link_options(Tasmanian_gridtest         PRIVATE -fsycl)
endif()

target_link_libraries(Tasmanian_libsparsegrid  Tasmanian_dependencies)

target_include_directories(Tasmanian_libsparsegrid PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>)

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

foreach(_tsgtarget gridtest benchmarksgrid)
    target_link_libraries(Tasmanian_${_tsgtarget} Tasmanian_libsparsegrid)
    set_target_properties(Tasmanian_${_tsgtarget} PROPERTIES OUTPUT_NAME "${_tsgtarget}" CXX_EXTENSIONS OFF)
    Tasmanian_rpath_target(TARGET Tasmanian_${_tsgtarget} USE_CURRENT)
    if (Tasmanian_ENABLE_OPENMP)
        # the OpenMP libraries are carried transitively from sparse grids library
        target_compile_options(Tasmanian_${_tsgtarget} PRIVATE ${OpenMP_CXX_FLAGS})
    endif()
    if (Tasmanian_ENABLE_HIP)
        target_link_libraries(Tasmanian_${_tsgtarget} hip::host)
    endif()
endforeach()
unset(_tsgtarget)

# data file, needed for testing and reference about custom rule definitions
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaussPattersonRule.table"  "${CMAKE_CURRENT_BINARY_DIR}/GaussPattersonRule.table" COPYONLY)


########################################################################
# Windows specific support (DLL export/import directives and names)
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")

    foreach(_tsg_target libsparsegrid gridtest)
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
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Examples/"
        DESTINATION "share/Tasmanian/examples/"
        FILES_MATCHING PATTERN "*.cpp"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/GaussPattersonRule.table"
        DESTINATION "share/Tasmanian/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
