cmake_minimum_required(VERSION 3.5)

########################################################################
# Sparse Grid librareis and command line tool
########################################################################

macro(Tasmanian_macro_add_libsparsegrid Tasmanian_shared_or_static)

if (${Tasmanian_shared_or_static} STREQUAL SHARED)
    set(Tasmanian_libtsg_target_name "Tasmanian_libsparsegrid_shared")
else()
    set(Tasmanian_libtsg_target_name "Tasmanian_libsparsegrid_static")
endif()

if (Tasmanian_ENABLE_CUDA)
    cuda_add_library(${Tasmanian_libtsg_target_name} ${Tasmanian_shared_or_static} ${Tasmanian_source_libsparsegrid} ${Tasmanian_source_libsparsegrid_cuda})
else()
    add_library(${Tasmanian_libtsg_target_name} ${Tasmanian_shared_or_static} ${Tasmanian_source_libsparsegrid})
endif()

target_include_directories(${Tasmanian_libtsg_target_name} PUBLIC $<INSTALL_INTERFACE:include>)

target_include_directories(${Tasmanian_libtsg_target_name} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/SparseGrids/>)
target_include_directories(${Tasmanian_libtsg_target_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/configured/>)

set_target_properties(${Tasmanian_libtsg_target_name} PROPERTIES OUTPUT_NAME "tasmaniansparsegrid")

# setup the name, in MacOSX and Linux static libraries use .a extension and shared libs use .dylib and .so
# in Windows, shared .dll libs need another file for linking which shares .lib extension with the static libs
# thus on Windows, we specifically add _static to the name of the lib
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    target_compile_definitions(${Tasmanian_libtsg_target_name} -D_TASMANIAN_WINDOWS_)
    target_compile_definitions(${Tasmanian_libtsg_target_name} -D_SCL_SECURE_NO_WARNINGS)
    target_compile_definitions(${Tasmanian_libtsg_target_name} -D_USE_MATH_DEFINES)
    if (${Tasmanian_shared_or_static} STREQUAL STATIC)
        set_target_properties(${Tasmanian_libtsg_target_name} PROPERTIES OUTPUT_NAME "tasmaniansparsegrid_static")
    endif()
endif()

if (Tasmanian_ENABLE_BLAS)
    target_link_libraries(${Tasmanian_libtsg_target_name} ${BLAS_LIBRARIES})
endif()

if (Tasmanian_ENABLE_CUBLAS)
    target_link_libraries(${Tasmanian_libtsg_target_name} ${CUDA_cusparse_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES})

    if (DEFINED CUDA_INCLUDE_DIRS)
        target_include_directories(${Tasmanian_libtsg_target_name} PUBLIC ${CUDA_INCLUDE_DIRS})
    endif()
endif()

if (Tasmanian_ENABLE_OPENMP)
    if (OpenMP_CXX_LIBRARIES)
        # using the OpenMP target leads to a problem with the exports
        # the OpenMP target cannot be exported, which means that a
        # project importing an already installed Tasmanian would
        # have to "know" whether Tasmanian was build with OpenMP and
        # call find_package(OpenMP) manually
        # Furthermore, using find_package(OpenMP) from a different
        # compiler can generate a wrong target, e.g., building Tasmanian
        # with gcc links to libgomp, but calling find_package(OpenMP)
        # from clang will create an OpenMP target that uses libiomp
        target_link_libraries(${Tasmanian_libtsg_target_name} ${OpenMP_CXX_LIBRARIES})
        #target_compile_options(Tasmanian_tasgrid PRIVATE ${OpenMP_CXX_FLAGS})
        target_compile_options(${Tasmanian_libtsg_target_name} PRIVATE ${OpenMP_CXX_FLAGS})
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

endmacro()

########################################################################
# Define source files
# in order to avoid GLOB, list all source files so they can be used for
# add_library and cuda_add_library in both shared and static cases
########################################################################
set(Tasmanian_source_libsparsegrid
               ${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSparseGrid.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSparseGrid.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgAcceleratedDataStructures.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgAcceleratedDataStructures.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgCacheLagrange.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgCoreOneDimensional.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgCoreOneDimensional.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgCudaMacros.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgEnumerates.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridCore.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridCore.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridGlobal.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridGlobal.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridWavelet.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridWavelet.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgHardCodedTabulatedRules.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgHardCodedTabulatedRules.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgHiddenExternals.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridLocalPolynomial.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridLocalPolynomial.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridSequence.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgGridSequence.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgIndexManipulator.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgIndexManipulator.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgIndexSets.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgIndexSets.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgLinearSolvers.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgLinearSolvers.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgOneDimensionalWrapper.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgOneDimensionalWrapper.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgRuleLocalPolynomial.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgRuleLocalPolynomial.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgRuleWavelet.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgRuleWavelet.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgSequenceOptimizer.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgSequenceOptimizer.cpp)

# source files specific to cuda, used for both static and shared libs
set(Tasmanian_source_libsparsegrid_cuda
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgCudaKernels.cu
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgCudaBasisEvaluations.hpp
               ${CMAKE_CURRENT_SOURCE_DIR}/tsgCudaLinearAlgebra.hpp)

set(Tasmanian_source_tasgrid
                ${CMAKE_CURRENT_SOURCE_DIR}/tasgrid_main.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSparseGrid.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/tasgridExternalTests.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/tasgridExternalTests.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/tasgridTestFunctions.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/tasgridTestFunctions.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/tasgridWrapper.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/tasgridWrapper.cpp)

########################################################################
# add the tasgrid and examples executables
########################################################################
add_executable(Tasmanian_tasgrid ${CMAKE_CURRENT_SOURCE_DIR}/tasgrid_main.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSparseGrid.hpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/tasgridExternalTests.hpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/tasgridExternalTests.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/tasgridTestFunctions.hpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/tasgridTestFunctions.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/tasgridWrapper.hpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/tasgridWrapper.cpp)
add_executable(Tasmanian_example_sparse_grids "${CMAKE_CURRENT_SOURCE_DIR}/Examples/example_sparse_grids.cpp")
list(APPEND Tasmanian_install_sparsegrid_targets Tasmanian_tasgrid)

set_target_properties(Tasmanian_tasgrid PROPERTIES OUTPUT_NAME "tasgrid")
set_target_properties(Tasmanian_example_sparse_grids PROPERTIES OUTPUT_NAME "example_sparse_grids")

########################################################################
# add the libraries
########################################################################
if (Tasmanian_SHARED_LIBRARY)
    Tasmanian_macro_add_libsparsegrid(SHARED)
    list(APPEND Tasmanian_install_sparsegrid_targets Tasmanian_libsparsegrid_shared)
endif()

if (Tasmanian_STATIC_LIBRARY)
    Tasmanian_macro_add_libsparsegrid(STATIC)

    target_link_libraries(Tasmanian_tasgrid              Tasmanian_libsparsegrid_static)
    target_link_libraries(Tasmanian_example_sparse_grids Tasmanian_libsparsegrid_static)

    list(APPEND Tasmanian_install_sparsegrid_targets Tasmanian_libsparsegrid_static)
else()
    target_link_libraries(Tasmanian_tasgrid              Tasmanian_libsparsegrid_shared)
    target_link_libraries(Tasmanian_example_sparse_grids Tasmanian_libsparsegrid_shared)
endif()

# hack a dependency problem in parallel make
# without this, if both static and shared libs are enabled, make -j tries to compile shared cuda kernels twice
# which creates a race condition and the build randomly fails
if (Tasmanian_SHARED_LIBRARY AND Tasmanian_STATIC_LIBRARY)
    add_dependencies(Tasmanian_libsparsegrid_shared Tasmanian_libsparsegrid_static)
endif()

# data file, needed for testing and reference about custom rule definitions
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaussPattersonRule.table"  "${CMAKE_BINARY_DIR}/GaussPattersonRule.table" COPYONLY)

########################################################################
# Windows specific support (DLL export/import directives and names)
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSparseGrid.windows.hpp"  "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSparseGrid.hpp" COPYONLY)

    if (Tasmanian_SHARED_LIBRARY)
        target_compile_definitions(Tasmanian_libsparsegrid_shared PUBLIC -DTSG_DLL)
    else()
        target_compile_definitions(Tasmanian_tasgrid PUBLIC -DTSG_DYNAMIC)
    endif()

    if (Tasmanian_STATIC_LIBRARY)
        target_compile_definitions(Tasmanian_libsparsegrid_static PUBLIC -DTSG_STATIC)
    endif()

    target_compile_definitions(Tasmanian_tasgrid -D_TASMANIAN_WINDOWS_) # overwrittes gettime()
    target_compile_definitions(Tasmanian_tasgrid -D_SCL_SECURE_NO_WARNINGS) # suppresses warnings regarding pointers to the middle of an array
    target_compile_definitions(Tasmanian_tasgrid -D_USE_MATH_DEFINES) # needed to include M_PI constant (lots of targets need this, will figure it out)
endif()


########################################################################
# Testing
########################################################################
add_test(SparseGridsAcceleration ${CMAKE_CURRENT_BINARY_DIR}/tasgrid -test acceleration -gpuid ${Tasmanian_TESTS_GPU_ID})
add_test(SparseGridsDomain       ${CMAKE_CURRENT_BINARY_DIR}/tasgrid -test domain)
add_test(SparseGridsRefinement   ${CMAKE_CURRENT_BINARY_DIR}/tasgrid -test refinement)
add_test(SparseGridsGlobal       ${CMAKE_CURRENT_BINARY_DIR}/tasgrid -test global)
add_test(SparseGridsLocal        ${CMAKE_CURRENT_BINARY_DIR}/tasgrid -test local)
add_test(SparseGridsWavelet      ${CMAKE_CURRENT_BINARY_DIR}/tasgrid -test wavelet)
set(Tasmanian_tests_list_sparse_grids SparseGridsAcceleration SparseGridsDomain SparseGridsRefinement SparseGridsGlobal SparseGridsLocal SparseGridsWavelet)


########################################################################
# Install headers and config files
########################################################################
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "*.windows.*" EXCLUDE
        PATTERN "Examples" EXCLUDE
        PATTERN "tsgHiddenExternals.hpp" EXCLUDE)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TasmanianSparseGrid.h"
        DESTINATION include
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Examples/"
        DESTINATION "examples"
        FILES_MATCHING PATTERN "*.cpp"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(TARGETS ${Tasmanian_install_sparsegrid_targets}
        EXPORT Tasmanian_exports
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")
