cmake_minimum_required(VERSION 2.8)
project(Tasmanian LANGUAGES CXX)

# Test push
########################################################################
#  User specified options:
#    -D Tasmanian_ENABLE_OPENMP:BOOL=ON
#    -D Tasmanian_ENABLE_BLAS:BOOL=ON
#    -D Tasmanian_ENABLE_MPI:BOOL=OFF
#    -D Tasmanian_ENABLE_CUBLAS:BOOL=ON
#    -D Tasmanian_ENABLE_CUDA:BOOL=OFF
#    -D Tasmanian_ENABLE_MATLAB:BOOL=OFF
#    -D Tasmanian_ENABLE_PYTHON:BOOL=ON
#    -D Tasmanian_SHARED_LIBRARY:BOOL=ON
#    -D Tasmanian_STATIC_LIBRARY:BOOL=ON
#    -D Tasmanian_MATLAB_WORK_FOLDER:PATH="${CMAKE_INSTALL_PREFIX}/tsgMatlabWork/"
#
# Additional options:
#    -D TPL_CUBLAS_INCLUDE_DIRS:PATH="do not use FindCUDA"
#    -D TPL_CUBLAS_LIBRARIES:PATH="do not use FindCUDA"
#    -D TPL_BLAS_INCLUDE_DIRS:PATH="do not use FindBLAS"
#    -D TPL_BLAS_LIBRARIES:PATH="do not use FindBLAS"
#    -D Tasmanian_FORCE_C11:BOOL=OFF
#    -D Tasmanian_CXX_FLAGS:STRING="overwrittes default flags"
#    -D Tasmanian_EXTRA_CXX_FLAGS:STRING="appends to the default flags"
#    -D Tasmanian_EXTRA_INCLUDE_DIRS:PATH="appends more include paths"
#    -D Tasmanian_EXTRA_LIBRARIES:STRING="appends more link libraries"
#    -D Tasmanian_EXTRA_LINK_DIRS:PATH="appends more link paths"
#
# Development only options:
#    -D Tasmanian_DEVELOPMENT_EXAMPLES:BOOL=OFF
#    -D Tasmanian_DEVELOPMENT_BACKWARDS:BOOL=OFF
#
# xSDK options:
#    -D USE_XSDK_DEFAULTS:BOOL=OFF
########################################################################
option(USE_XSDK_DEFAULTS "Enable xSDK compatibility" OFF)
option(Tasmanian_ENABLE_OPENMP "Enable OpenMP support for Tasmanian" ON)
option(Tasmanian_ENABLE_BLAS "Enable CPU Blas support for Tasmanian" ON)
option(Tasmanian_ENABLE_MPI "Enable MPI support for Tasmanian" OFF)
option(Tasmanian_ENABLE_MATLAB "Enable MATLAB interface for Tasmanian" OFF)
option(Tasmanian_ENABLE_PYTHON "Enable Python interface for Tasmanian" ON)
option(Tasmanian_ENABLE_FORTRAN "Enable Fortran interface for Tasmanian" OFF)
option(Tasmanian_ENABLE_CUBLAS "Enable Nvidia GPU acceleration using cuBLAS for Tasmanian" ON)
option(Tasmanian_ENABLE_CUDA "Enable Nvidia GPU acceleration using CUDA kernels for Tasmanian" OFF)
option(Tasmanian_STATIC_LIBRARY "Compile static library and the executables with static linkage" ON)
option(Tasmanian_SHARED_LIBRARY "Compile dynamic library" ON)
set(Tasmanian_PREFER_C11 OFF)

if ((NOT DEFINED CMAKE_BUILD_TYPE) OR (NOT CMAKE_BUILD_TYPE))
    set(CMAKE_BUILD_TYPE Debug)
endif()

# version and license
set(Tasmanian_VersionMajor 5)
set(Tasmanian_VersionMinor 1)
set(Tasmanian_VersionString "${Tasmanian_VersionMajor}.${Tasmanian_VersionMinor}")
set(Tasmanian_License "BSD 3-Clause with UT-Battelle disclamer")

########################################################################
# list all source files
########################################################################
set(tasgrid_src SparseGrids/tasgrid_main.cpp
                SparseGrids/TasmanianSparseGrid.hpp
                SparseGrids/tasgridExternalTests.hpp
                SparseGrids/tasgridExternalTests.cpp
                SparseGrids/tasgridTestFunctions.hpp
                SparseGrids/tasgridTestFunctions.cpp
                SparseGrids/tasgridWrapper.hpp
                SparseGrids/tasgridWrapper.cpp)

set(libtsg_src SparseGrids/TasmanianSparseGrid.hpp
               SparseGrids/TasmanianSparseGrid.cpp
               SparseGrids/tsgAcceleratedDataStructures.hpp
               SparseGrids/tsgAcceleratedDataStructures.cpp
               SparseGrids/tsgCacheLagrange.hpp
               SparseGrids/tsgCacheLagrange.cpp
               SparseGrids/tsgCoreOneDimensional.hpp
               SparseGrids/tsgCoreOneDimensional.cpp
               SparseGrids/tsgEnumerates.hpp
               SparseGrids/tsgGridCore.hpp
               SparseGrids/tsgGridCore.cpp
               SparseGrids/tsgGridGlobal.hpp
               SparseGrids/tsgGridGlobal.cpp
               SparseGrids/tsgGridWavelet.hpp
               SparseGrids/tsgGridWavelet.cpp
               SparseGrids/tsgHardCodedTabulatedRules.hpp
               SparseGrids/tsgHardCodedTabulatedRules.cpp
               SparseGrids/tsgHiddenExternals.hpp
               SparseGrids/tsgGridLocalPolynomial.hpp
               SparseGrids/tsgGridLocalPolynomial.cpp
               SparseGrids/tsgGridSequence.hpp
               SparseGrids/tsgGridSequence.cpp
               SparseGrids/tsgIndexManipulator.hpp
               SparseGrids/tsgIndexManipulator.cpp
               SparseGrids/tsgIndexSets.hpp
               SparseGrids/tsgIndexSets.cpp
               SparseGrids/tsgLinearSolvers.hpp
               SparseGrids/tsgLinearSolvers.cpp
               SparseGrids/tsgOneDimensionalWrapper.hpp
               SparseGrids/tsgOneDimensionalWrapper.cpp
               SparseGrids/tsgRuleLocalPolynomial.hpp
               SparseGrids/tsgRuleLocalPolynomial.cpp
               SparseGrids/tsgRuleWavelet.hpp
               SparseGrids/tsgRuleWavelet.cpp
               SparseGrids/tsgSequenceOptimizer.hpp
               SparseGrids/tsgSequenceOptimizer.cpp)

set(libtdr_src DREAM/TasmanianDREAM.hpp
               DREAM/TasmanianDREAM.cpp
               DREAM/tdrEnumerates.hpp
               DREAM/tdrCorePDF.hpp
               DREAM/tdrCorePDF.cpp
               SparseGrids/TasmanianSparseGrid.hpp)

set(tasdre_src DREAM/tasdream_main.cpp
               DREAM/tasdreamExternalTests.hpp
               DREAM/tasdreamExternalTests.cpp
               DREAM/tasdreamTestPDFs.hpp
               DREAM/tasdreamTestPDFs.cpp
               DREAM/tasdreamBenchmark.hpp
               DREAM/tasdreamBenchmark.cpp)

set(tasmat_src tsgCancelRefine.m
               tsgCleanTempFiles.m
               tsgCopyGrid.m
               tsgDeleteGridByName.m
               tsgDeleteGrid.m
               tsgEstimateAnisotropicCoefficients.m
               tsgEvaluateHierarchy.m
               tsgEvaluate.m
               tsgExample.m
               tsgGetHCoefficients.m
               tsgGetInterpolationWeights.m
               tsgGetNeededIndexes.m
               tsgGetNeededPoints.m
               tsgGetPointsIndexes.m
               tsgGetPoints.m
               tsgGetPolynomialSpaces.m
               tsgGetQuadrature.m
               tsgIntegrate.m
               tsgListGridsByName.m
               tsgLoadHCoefficients.m
               tsgLoadValues.m
               tsgMakeFilenames.m
               tsgMakeGlobal.m
               tsgMakeLocalPolynomial.m
               tsgMakeQuadrature.m
               tsgMakeSequence.m
               tsgMakeWavelet.m
               tsgMergeRefine.m
               tsgPlotPoints2D.m
               tsgReadMatrix.m
               tsgRefineAnisotropic.m
               tsgRefineSurplus.m
               tsgReloadGrid.m
               tsgSummary.m
               tsgCoreTests.m
               tsgWriteCustomRuleFile.m
               tsgWriteMatrix.m)

set(libfor_src InterfaceFortran/TasmanianSG.f90
               InterfaceFortran/tsgC2FortranBridge.f90
               InterfaceFortran/tsgC2Fortran.cpp)

set(tasfor_src Testing/fortester.f90)


########################################################################
# sanity check and xSDK compatibility
########################################################################
if (USE_XSDK_DEFAULTS)
    add_definitions(-DTASMANIAN_XSDK)
    if (DEFINED BUILD_SHARED_LIBS)
        if (BUILD_SHARED_LIBS)
            set(Tasmanian_SHARED_LIBRARY ON)
        else()
            set(Tasmanian_SHARED_LIBRARY OFF)
            set(Tasmanian_STATIC_LIBRARY ON)
        endif()
    endif()
    if (DEFINED XSDK_ENABLE_PYTHON)
        set(Tasmanian_ENABLE_PYTHON XSDK_ENABLE_PYTHON)
    else()
        set(Tasmanian_ENABLE_PYTHON OFF)
    endif()
    if (DEFINED XSDK_ENABLE_CUDA)
        set(Tasmanian_ENABLE_CUBLAS XSDK_ENABLE_CUDA)
    else()
        set(Tasmanian_ENABLE_CUBLAS OFF)
    endif()
else()
    if (DEFINED BUILD_SHARED_LIBS)
        if (NOT ((BUILD_SHARED_LIBS AND Tasmanian_SHARED_LIBRARY) OR ((NOT BUILD_SHARED_LIBS) AND (NOT Tasmanian_SHARED_LIBRARY))))
            message(WARNING "-D BUILD_SHARED_LIBS set, but conflicting with Tasmanian_SHARED_LIBRARY, defaulting to Tasmanian_SHARED_LIBRARY, did you forget to set -D USE_XSDK_DEFAULTS:BOOL=ON")
        endif()
    endif()
    if (DEFINED XSDK_ENABLE_PYTHON)
        if (NOT ((XSDK_ENABLE_PYTHON AND Tasmanian_ENABLE_PYTHON) OR ((NOT XSDK_ENABLE_PYTHON) AND (NOT Tasmanian_ENABLE_PYTHON))))
            message(WARNING "-D XSDK_ENABLE_PYTHON set, but conflicting with Tasmanian_ENABLE_PYTHON, defaulting to Tasmanian_ENABLE_PYTHON, did you forget to set -D USE_XSDK_DEFAULTS:BOOL=ON")
        endif()
    endif()
    if (DEFINED XSDK_ENABLE_CUDA)
        if (NOT ((XSDK_ENABLE_CUDA AND Tasmanian_ENABLE_CUBLAS) OR ((NOT XSDK_ENABLE_CUDA) AND (NOT Tasmanian_ENABLE_CUBLAS))))
            message(WARNING "-D XSDK_ENABLE_CUDA set, but conflicting with Tasmanian_ENABLE_CUBLAS, defaulting to Tasmanian_ENABLE_CUBLAS, did you forget to set -D USE_XSDK_DEFAULTS:BOOL=ON")
        endif()
    endif()
endif()

if (Tasmanian_ENABLE_PYTHON)
    find_package(PythonInterp)
    if (NOT PYTHONINTERP_FOUND)
        if (USE_XSDK_DEFAULTS)
            message(FATAL_ERROR "-D XSDK_ENABLE_PYTHON set, but find_package(PythonInterp) failed\nuse -D PYTHON_EXECUTABLE:PATH to specify python interpreter")
        else()
            message("-D Tasmanian_ENABLE_PYTHON is ON, but find_package(PythonInterp) failed\nuse -D PYTHON_EXECUTABLE:PATH to specify python interpreter\noverwritting option: -D Tasmanian_ENABLE_PYTHON:BOOL=ON")
            set(Tasmanian_ENABLE_PYTHON OFF)
        endif()
    endif()
endif()

if ((NOT Tasmanian_STATIC_LIBRARY) AND (NOT Tasmanian_SHARED_LIBRARY))
    if (Tasmanian_ENABLE_PYTHON)
        message(WARNING "must specifiy at least one of Tasmanian_STATIC_LIBRARY or Tasmanian_SHARED_LIBRARY\noverwritting option: -D Tasmanian_SHARED_LIBRARY:BOOL=ON")
        set(Tasmanian_SHARED_LIBRARY ON)
    else()
        message(WARNING "must specifiy at least one of Tasmanian_STATIC_LIBRARY or Tasmanian_SHARED_LIBRARY\noverwritting option: -D Tasmanian_STATIC_LIBRARY:BOOL=ON")
        set(Tasmanian_STATIC_LIBRARY ON)
    endif()
endif()

if (Tasmanian_ENABLE_PYTHON AND (NOT Tasmanian_SHARED_LIBRARY))
    if (USE_XSDK_DEFAULTS)
        message(FATAL_ERROR "BUILD_SHARED_LIBS has to be ON for the Tasmanian Python interface")
    else()
        message(WARNING "Tasmanian_SHARED_LIBRARY is needed by the python interface\nto disable the shared library, you must disable python as well\noverwritting option: -D Tasmanian_ENABLE_PYTHON:BOOL=ON")
    endif()
    set(Tasmanian_SHARED_LIBRARY ON)
endif()

if (NOT DEFINED Tasmanian_MATLAB_WORK_FOLDER)
    if (Tasmanian_ENABLE_MATLAB)
        message(FATAL_ERROR "Tasmanian_ENABLE_MATLAB is ON, but Tasmanian_MATLAB_WORK_FOLDER is missing\nIn order to enable the MATLAB interface, you must specify -D Tasmanian_MATLAB_WORK_FOLDER:PATH=")
    endif()
else()
    set(TSGMatlabWork ${Tasmanian_MATLAB_WORK_FOLDER})
endif()

# Pre-process CUDA for the add_library vs cuda_add_library
if (Tasmanian_ENABLE_CUBLAS OR Tasmanian_ENABLE_CUDA)
    if (((NOT DEFINED TPL_CUBLAS_INCLUDE_DIRS) AND (NOT DEFINED TPL_CUBLAS_LIBRARIES)) OR Tasmanian_ENABLE_CUDA)
        find_package(CUDA)

        if (CUDA_FOUND)
            set(Tasmanian_USING_FINDCUDA ON)
            set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Wno-deprecated-gpu-targets")
            set(TPL_CUBLAS_INCLUDE_DIRS ${CUDA_INCLUDE_DIR})
            set(TPL_CUBLAS_LIBRARIES ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY})
            if (Tasmanian_ENABLE_CUDA)
                set(libtsg_src SparseGrids/tsgCudaKernels.cu ${libtsg_src})
                add_definitions(-DTASMANIAN_CUDA)
            endif()
        else()
            if (USE_XSDK_DEFAULTS)
                message(FATAL_ERROR "Tasmanian_ENABLE_CUBLAS or Tasmanian_ENABLE_CUDA is on, but find_package(CUDA) failed\nTPL_CUBLAS_INCLUDE_DIRS and TPL_CUBLAS_LIBRARIES were not manually set either")
            else()
                message(WARNING "find_package(CUDA) failed, could not find cuda")
                if (Tasmanian_ENABLE_CUBLAS)
                    message(WARNING "overwritting option: -D Tasmanian_ENABLE_CUBLAS:BOOL=OFF")
                    set(Tasmanian_ENABLE_CUBLAS OFF)
                endif()
                if (Tasmanian_ENABLE_CUDA)
                    message(WARNING "overwritting option: -D Tasmanian_ENABLE_CUDA:BOOL=OFF")
                    set(Tasmanian_ENABLE_CUDA OFF)
                endif()
                set(Tasmanian_USING_FINDCUDA OFF)
            endif()
        endif()
    else()
        set(Tasmanian_USING_FINDCUDA OFF)
    endif()
else()
    set(Tasmanian_USING_FINDCUDA OFF)
endif()

if (Tasmanian_ENABLE_CUBLAS)
    include_directories(${TPL_CUBLAS_INCLUDE_DIRS})
    add_definitions(-DTASMANIAN_CUBLAS)
endif()

if (Tasmanian_ENABLE_FORTRAN)
    enable_language(Fortran)
    if(NOT CMAKE_Fortran_COMPILER)
        if (USE_XSDK_DEFAULTS)
            message(FATAL_ERROR "Cannot find Fortran!")
        else()
            message(WARNING "Cannot find Fortran\noverwritting option: -D Tasmanian_ENABLE_FORTRAN:BOOL=OFF")
            set(Tasmanian_ENABLE_FORTRAN OFF)
        endif()
    endif()
endif()

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    add_definitions(-D_TASMANIAN_DEBUG_)
    if (Tasmanian_ENABLE_FORTRAN)
        set(CMAKE_Fortran_FLAGS_DEBUG "-g")
    endif()
else()
    # this is a bug in some versions of CUDA
    if (Tasmanian_USING_FINDCUDA AND ("${CUDA_VERSION_MAJOR}" LESS 8))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES")
    endif()
endif()

########################################################################
# setup targets
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(CMAKE_MACOSX_RPATH ON)
    # code copied from https://cmake.org/Wiki/CMake_RPATH_handling
    # put RPATH in the build three (stage 1, see below)
	set(CMAKE_SKIP_BUILD_RPATH OFF)
	# use the build as opposed install RPATH
	set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
	set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

	set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON) # needed by TPL
endif()

include_directories("${PROJECT_SOURCE_DIR}/SparseGrids/")

set(TSGBaseLibName "tasmaniansparsegrid")
set(TDRBaseLibName "tasmaniandream")
set(FortranBaseLibName "tasmanianfortran")

if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(TSGStaticLibName "${TSGBaseLibName}_static") # Unix uses .a and .so to distinguish static and shared
    set(TDRStaticLibName "${TDRBaseLibName}_static") # Windows uses .lib in both cases, append _static
    set(FortranStaticLibName "${FortranBaseLibName}_static")
else()
    set(TSGStaticLibName ${TSGBaseLibName})
    set(TDRStaticLibName ${TDRBaseLibName})
    set(FortranStaticLibName ${FortranBaseLibName})
endif()

if (Tasmanian_USING_FINDCUDA)
    cuda_add_executable(tasgrid ${tasgrid_src})
    cuda_add_executable(tasdream ${tasdre_src})
else()
    add_executable(tasgrid ${tasgrid_src})
    add_executable(tasdream ${tasdre_src})
endif()

set(tas_target_list tasgrid tasdream)
set(install_targets tasgrid tasdream)

if (Tasmanian_ENABLE_FORTRAN)
    add_executable(fortester ${tasfor_src})
    set(tas_target_list ${tas_target_list} fortester)
endif()

if (Tasmanian_SHARED_LIBRARY)
    if (Tasmanian_USING_FINDCUDA)
        cuda_add_library(libtsg_shared SHARED ${libtsg_src})
        cuda_add_library(libtdr_shared SHARED ${libtdr_src})
    else()
        add_library(libtsg_shared SHARED ${libtsg_src})
        add_library(libtdr_shared SHARED ${libtdr_src})
    endif()

    set(tas_target_list ${tas_target_list} libtsg_shared libtdr_shared)
    target_link_libraries(libtdr_shared libtsg_shared)

    set_target_properties(libtsg_shared PROPERTIES OUTPUT_NAME "${TSGBaseLibName}")
    set_target_properties(libtdr_shared PROPERTIES OUTPUT_NAME "${TDRBaseLibName}")

    set(install_targets ${install_targets} libtsg_shared libtdr_shared)

    set_property(TARGET libtsg_shared APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)
    set_property(TARGET libtdr_shared APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)

    if (Tasmanian_ENABLE_FORTRAN)
        add_library(libtasfor_shared SHARED ${libfor_src})
        target_link_libraries(libtasfor_shared libtsg_shared)
        set(install_targets ${install_targets} libtasfor_shared)
        set_property(TARGET libtasfor_shared APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)
        set(tas_target_list ${tas_target_list} libtasfor_shared)
        set_target_properties(libtasfor_shared PROPERTIES OUTPUT_NAME "${FortranBaseLibName}")
    endif()
endif()

if (Tasmanian_STATIC_LIBRARY)
    if (Tasmanian_USING_FINDCUDA)
        cuda_add_library(libtsg_static STATIC ${libtsg_src})
        cuda_add_library(libtdr_static STATIC ${libtdr_src})
    else()
        add_library(libtsg_static STATIC ${libtsg_src})
        add_library(libtdr_static STATIC ${libtdr_src})
    endif()

    set(tas_target_list ${tas_target_list} libtsg_static libtdr_static)
    target_link_libraries(libtdr_static libtsg_static)

    set_target_properties(libtsg_static PROPERTIES OUTPUT_NAME "${TSGStaticLibName}")
    set_target_properties(libtdr_static PROPERTIES OUTPUT_NAME "${TDRStaticLibName}")

    set(install_targets ${install_targets} libtsg_static libtdr_static)

    target_link_libraries(tasgrid libtsg_static)
    target_link_libraries(tasdream libtdr_static libtsg_static)

    set_property(TARGET libtsg_static APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)
    set_property(TARGET libtdr_static APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)

    if (Tasmanian_ENABLE_FORTRAN)
        add_library(libtasfor_static STATIC ${libfor_src})
        target_link_libraries(libtasfor_static libtsg_static)
        set(tas_target_list ${tas_target_list} libtasfor_static)
        set(install_targets ${install_targets} libtasfor_static)
        set_property(TARGET libtasfor_static APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>)
        set_target_properties(libtasfor_static PROPERTIES OUTPUT_NAME "${FortranStaticLibName}")
        target_link_libraries(fortester libtasfor_static)
    endif()
else()
    target_link_libraries(tasgrid libtsg_shared)
    target_link_libraries(tasdream libtdr_shared libtsg_shared)
    if (Tasmanian_ENABLE_FORTRAN)
        target_link_libraries(fortester libtasfor_shared)
    endif()
endif()

# for development purposes, compile the exmaples together with everything else
if (DEFINED Tasmanian_DEVELOPMENT_EXAMPLES AND Tasmanian_DEVELOPMENT_EXAMPLES)
    if (Tasmanian_USING_FINDCUDA)
        cuda_add_executable(example_sparse_grids "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.cpp")
        cuda_add_executable(example_dream "${PROJECT_SOURCE_DIR}/Examples/example_dream.cpp")
    else()
        add_executable(example_sparse_grids "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.cpp")
        add_executable(example_dream "${PROJECT_SOURCE_DIR}/Examples/example_dream.cpp")
    endif()
    set(tas_target_list ${tas_target_list} example_sparse_grids example_dream)
    include_directories("${PROJECT_SOURCE_DIR}/DREAM/")
    if (Tasmanian_STATIC_LIBRARY)
        target_link_libraries(example_sparse_grids libtsg_static)
        target_link_libraries(example_dream libtdr_static)
    else()
        target_link_libraries(example_sparse_grids libtsg_shared)
        target_link_libraries(example_dream libtdr_shared)
    endif()
    if (Tasmanian_ENABLE_FORTRAN)
        add_executable(example_sparse_grids_f90 "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.f90")
        set_target_properties(example_sparse_grids_f90 PROPERTIES OUTPUT_NAME "example_sparse_grids_fortran")
        set(tas_target_list ${tas_target_list} example_sparse_grids_f90)
        if (Tasmanian_STATIC_LIBRARY)
            target_link_libraries(example_sparse_grids_f90 libtasfor_static libtsg_static)
        else()
            target_link_libraries(example_sparse_grids_f90 libtasfor_shared libtsg_shared)
        endif()
    endif()
endif()


########################################################################
# Windows specific support (DLL export/import directives)
########################################################################
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    configure_file("${PROJECT_SOURCE_DIR}/SparseGrids/TasmanianSparseGrid.windows.hpp"  "${PROJECT_SOURCE_DIR}/SparseGrids/TasmanianSparseGrid.hpp" COPYONLY)

    if (Tasmanian_SHARED_LIBRARY)
        target_compile_definitions(libtsg_shared PUBLIC -DTSG_DLL)
        target_compile_definitions(libtdr_shared PUBLIC -DTSG_DLL)
    else()
        target_compile_definitions(tasgrid PUBLIC -DTSG_DYNAMIC)
        target_compile_definitions(tasdream PUBLIC -DTSG_DYNAMIC)
    endif()

    if (Tasmanian_STATIC_LIBRARY)
        target_compile_definitions(libtsg_static PUBLIC -DTSG_STATIC)
        target_compile_definitions(libtdr_static PUBLIC -DTSG_STATIC)
    endif()

    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif()


########################################################################
# BLAS setup
########################################################################
if (USE_XSDK_DEFAULTS)
    if (DEFINED TPL_ENABLE_BLAS)
        if (TPL_ENABLE_BLAS)
            set(Tasmanian_ENABLE_BLAS ON)
        else()
            set(Tasmanian_ENABLE_BLAS OFF)
        endif()
    else()
        set(Tasmanian_ENABLE_BLAS OFF)
    endif()
else()
    if (DEFINED TPL_ENABLE_BLAS)
        if (NOT ((TPL_ENABLE_BLAS AND Tasmanian_ENABLE_BLAS) OR ((NOT TPL_ENABLE_BLAS) AND (NOT Tasmanian_ENABLE_BLAS))))
            message(WARNING "-D TPL_ENABLE_BLAS defined to ON, but it conflicts with Tasmanian_ENABLE_BLAS. Defaulting to Tasmanian_ENABLE_BLAS. Did you forget to set -D USE_XSDK_DEFAULTS:BOOL=ON")
        endif()
    endif()
endif()

if (Tasmanian_ENABLE_BLAS)
    if ((NOT DEFINED TPL_BLAS_INCLUDE_DIRS) AND (NOT DEFINED TPL_BLAS_LIBRARIES))
        find_package(BLAS)

        if (BLAS_FOUND)
            message(STATUS "BLAS found: ${BLAS_LIBRARIES}")

            add_definitions(-DTASMANIAN_CPU_BLAS)
        else()
            if (USE_XSDK_DEFAULTS)
                message(FATAL_ERROR "BLAS requested, but TPL_BLAS_LIBRARIES and TPL_BLAS_INCLUDE_DIRS are not defined, and FindBLAS failed")
            else()
                message(STATUS "Cannont find BLAS, disabling -D Tasmanian_ENABLE_BLAS:BOOL=OFF")
                set(Tasmanian_ENABLE_BLAS OFF)
            endif()
        endif()
    else()
        add_definitions(-DTASMANIAN_CPU_BLAS)

        if (NOT DEFINED TPL_BLAS_INCLUDE_DIRS)
            set(TPL_BLAS_INCLUDE_DIRS "")
        endif()
        if (NOT DEFINED TPL_BLAS_LIBRARIES)
            set(TPL_BLAS_INCLUDE_DIRS "")
        endif()
        set(BLAS_LIBRARIES ${TPL_BLAS_LIBRARIES})
        set(BLAS_INCLUDE_DIR ${TPL_BLAS_INCLUDE_DIRS})
    endif()
endif()

if (Tasmanian_ENABLE_BLAS)
    include_directories(${BLAS_INCLUDE_DIR})

    foreach(tas_target ${tas_target_list})
        target_link_libraries(${tas_target} m ${BLAS_LIBRARIES})
    endforeach()
endif()

########################################################################
# MPI setup (experimental, DREAM distributed posterior only)
########################################################################
if (USE_XSDK_DEFAULTS)
    if (DEFINED TPL_ENABLE_MPI)
        if (TPL_ENABLE_MPI)
            set(Tasmanian_ENABLE_MPI ON)
        else()
            set(Tasmanian_ENABLE_MPI OFF)
        endif()
    else()
        set(Tasmanian_ENABLE_MPI ON) # for XSDK enable MPI by default
    endif()
else()
    if (DEFINED TPL_ENABLE_MPI)
        if (NOT ((TPL_ENABLE_MPI AND Tasmanian_ENABLE_MPI) OR ((NOT TPL_ENABLE_MPI) AND (NOT Tasmanian_ENABLE_MPI))))
            message(WARNING "-D TPL_ENABLE_MPI defined to ON, but it conflicts with Tasmanian_ENABLE_MPI. Defaulting to Tasmanian_ENABLE_MPI. Did you forget to set -D USE_XSDK_DEFAULTS:BOOL=ON")
        endif()
    endif()
endif()

if (Tasmanian_ENABLE_MPI)
    find_package(MPI REQUIRED)

    set(Tasmanian_PREFER_C11 ON)

    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(tasdream ${MPI_LIBRARIES})
    if (Tasmanian_SHARED_LIBRARY)
        target_link_libraries(libtdr_shared ${MPI_LIBRARIES})
        set_property(TARGET libtdr_shared APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_INCLUDE_PATH}")
    endif()
    if (Tasmanian_STATIC_LIBRARY)
        target_link_libraries(libtdr_static ${MPI_LIBRARIES})
        set_property(TARGET libtdr_static APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_INCLUDE_PATH}")
    endif()

    set(tdrIncludeMPI "#include <mpi.h>")
    #configure_file("${PROJECT_SOURCE_DIR}/DREAM/tdrEnableMPI.in.hpp" "${PROJECT_SOURCE_DIR}/DREAM/tdrEnableMPI.hpp")

    if(MPI_COMPILE_FLAGS)
        set_target_properties(tasdream PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
        if (Tasmanian_SHARED_LIBRARY)
            set_target_properties(libtdr_shared PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
        endif()
        if (Tasmanian_STATIC_LIBRARY)
            set_target_properties(libtdr_static PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
        endif()
    endif()
    if(MPI_LINK_FLAGS)
        set_target_properties(tasdream PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
        if (Tasmanian_SHARED_LIBRARY)
            set_target_properties(libtdr_shared PROPERTIES  LINK_FLAGS "${MPI_LINK_FLAGS}")
        endif()
        if (Tasmanian_STATIC_LIBRARY)
            set_target_properties(libtdr_static PROPERTIES  LINK_FLAGS "${MPI_LINK_FLAGS}")
        endif()
    endif()
else()
    set(tdrIncludeMPI "//#include <mpi.h>")
    #configure_file("${PROJECT_SOURCE_DIR}/DREAM/tdrEnableMPI.in.hpp"  "${PROJECT_SOURCE_DIR}/DREAM/tdrEnableMPI.hpp")
endif()

########################################################################
# CUDA setup (see also the cuda pre-process in the sanity check section)
########################################################################
if (Tasmanian_ENABLE_CUBLAS)
    set(Tasmanian_PREFER_C11 ON)
    foreach(tas_target ${tas_target_list})
        target_link_libraries(${tas_target} ${TPL_CUBLAS_LIBRARIES})
    endforeach()
endif()

########################################################################
# OpenMP setup
########################################################################
set(TAS_SUPPRESS_OMP_WARN OFF)
if (Tasmanian_ENABLE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else()
        if (USE_XSDK_DEFAULTS)
            message(FATAL_ERROR "could not find OpenMP")
        else()
            message(WARNING "could not find OpenMP, defaulting to single core\noverwritting option: -D Tasmanian_ENABLE_OPENMP:BOOL=OFF")
            set(Tasmanian_ENABLE_OPENMP OFF)
            set(TAS_SUPPRESS_OMP_WARN ON)
        endif()
    endif()
else()
    set(TAS_SUPPRESS_OMP_WARN ON)
endif()

########################################################################
# Cholmod support (experimental, wavelets only)
########################################################################
if (TPL_ENABLE_CHOLMOD OR (DEFINED TPL_CHOLMOD_INCLUDE_DIRS) OR (DEFINED TPL_CHOLMOD_LIBRARIES))
    add_definitions(-DTASMANIAN_CHOLMOD)
    if ((NOT DEFINED TPL_CHOLMOD_INCLUDE_DIRS) AND (NOT DEFINED TPL_CHOLMOD_LIBRARIES))
        set(TPL_CHOLMOD_LIBRARIES cholmod)
        set(TPL_CHOLMOD_INCLUDE_DIRS "")
    endif()

    include_directories(${TPL_CHOLMOD_INCLUDE_DIRS})

    foreach(tas_target ${tas_target_list})
        target_link_libraries(${tas_target} ${TPL_CHOLMOD_LIBRARIES})
    endforeach()
endif()

########################################################################
# C++ 2011 support
########################################################################
if (NOT DEFINED Tasmanian_FORCE_C11)
    set(Tasmanian_FORCE_C11 Tasmanian_PREFER_C11)
endif()

if (Tasmanian_FORCE_C11)
    foreach(tas_target ${tas_target_list})
        set_property(TARGET ${tas_target} PROPERTY CXX_STANDARD 11)
    endforeach()
endif()

########################################################################
# Compiler specific flags: Intel hasn't been tested in a while
########################################################################
if (NOT DEFINED Tasmanian_CXX_FLAGS)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -Wshadow -Wuninitialized -Wstrict-aliasing -pedantic")
        #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -mtune=native -Wall -Wextra -Wshadow -pedantic") # -Wno-unused-parameter
        #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage") # ./tasgrid -test; gcov CMakeFiles/libtsg_static.dir/SparseGrids/tsgIndexSets.cpp.; geany tsgIndexSets.cpp.gcov
        if (TAS_SUPPRESS_OMP_WARN)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
        endif()
    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native -diag-disable 11074 -diag-disable 11076 -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic")
    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -mtune=native -Wall -Wextra -Wshadow -pedantic")
        if (TAS_SUPPRESS_OMP_WARN)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
        endif()
    elseif (MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox /EHsc -D_SCL_SECURE_NO_WARNINGS")
    endif()
else()
    set(CMAKE_CXX_FLAGS ${Tasmanian_CXX_FLAGS})
endif()

if (Tasmanian_ENABLE_FORTRAN)
    if (DEFINED Tasmanian_Fortran_FLAGS)
        set(CMAKE_Fortran_FLAGS ${Tasmanian_Fortran_FLAGS})
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fno-f2c -mtune=native -Wall -Wextra -Wshadow -pedantic")
    endif()
endif()

if (DEFINED Tasmanian_EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Tasmanian_EXTRA_CXX_FLAGS}")
endif()

if (DEFINED Tasmanian_EXTRA_LIBRARIES)
    foreach(tas_target ${tas_target_list})
        target_link_libraries(${tas_target} ${Tasmanian_EXTRA_LIBRARIES})
    endforeach()
endif()

if (DEFINED Tasmanian_EXTRA_INCLUDE_DIRS)
    include_directories(${Tasmanian_EXTRA_INCLUDE_DIRS})
endif()

if (DEFINED Tasmanian_EXTRA_LINK_DIRS)
    link_directories(${Tasmanian_EXTRA_LINK_DIRS})
endif()

########################################################################
# Configure the .in. files
# using two stage configuration forcing two files for each .in. file
#   - everything should work in the build folder
#   - everything should work after make install
########################################################################

# stage 0: for support of simple GNU-Make, ignore for release and intermediate builds
if (DEFINED Tasmanian_DEVELOPMENT_BACKWARDS AND Tasmanian_DEVELOPMENT_BACKWARDS)
    configure_file("${PROJECT_SOURCE_DIR}/Config/tsgVersion.in.hpp"  "${PROJECT_SOURCE_DIR}/Config/tsgVersion.hpp")
    set(tdrIncludeMPIBackup ${tdrIncludeMPI})
    set(tdrIncludeMPI "//#include <mpi.h>")
    configure_file("${PROJECT_SOURCE_DIR}/Config/tdrEnableMPI.in.hpp"  "${PROJECT_SOURCE_DIR}/Config/tdrEnableMPI.hpp")
    set(tdrIncludeMPI ${tdrIncludeMPIBackup})
    set(tsgPython "/usr/bin/env python")
    set(TSGLibPyWithBlas "False")
    set(TSGLibPyWithCublas "False")
    set(TSGLibPyWithCuda "False")
    set(TSGEnableCmakeSyncTest "False")
    set(TSGLibShared "./libtasmaniansparsegrid.so")
    configure_file("${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py" "${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.py")
    configure_file("${PROJECT_SOURCE_DIR}/Testing/testTSG.in.py" "${PROJECT_SOURCE_DIR}/Testing/testTSG.py")
    set(TSGLibShared "libtasmaniansparsegrid.dll")
    configure_file("${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py" "${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.windows.py")
    set(TSGExampleImport "#")
    configure_file("${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.in.py" "${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.py")
    set(TSGExecutable "ENTER THE PATH TO tasgrid EXECUTABLE")
    set(TSGMatlabWork "ENTER THE PATH TO MATLAB WORK FOLDER")
    configure_file("${PROJECT_SOURCE_DIR}/InterfaceMATLAB/tsgGetPaths.in.m" "${PROJECT_SOURCE_DIR}/InterfaceMATLAB/tsgGetPaths.m")
endif()

# stage 1: build folder paths
set(TSGLibPath "${CMAKE_BINARY_DIR}")
set(TSGLibShared "${TSGLibPath}/${CMAKE_SHARED_LIBRARY_PREFIX}${TSGBaseLibName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(TSGLibStatic "${TSGLibPath}/${CMAKE_STATIC_LIBRARY_PREFIX}${TSGStaticLibName}${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(TDRLibShared "${TSGLibPath}/${CMAKE_SHARED_LIBRARY_PREFIX}${TDRBaseLibName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(TDRLibStatic "${TSGLibPath}/${CMAKE_STATIC_LIBRARY_PREFIX}${TDRStaticLibName}${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(TSGExecutable "${CMAKE_BINARY_DIR}/tasgrid")
set(TDRExecutable "${CMAKE_BINARY_DIR}/tasdream")
set(tsgPython "${PYTHON_EXECUTABLE}")
set(TSGLibPyWithBlas "False")
if (Tasmanian_ENABLE_BLAS)
    set(TSGLibPyWithBlas "True")
endif()
set(TSGLibPyWithCublas "False")
if (Tasmanian_ENABLE_CUBLAS)
    set(TSGLibPyWithCublas "True")
endif()
set(TSGLibPyWithCuda "False")
if (Tasmanian_ENABLE_CUDA)
    set(TSGLibPyWithCuda "True")
endif()
set(TSGPythonPath "${CMAKE_BINARY_DIR}")
set(TSGMatlabWork "${CMAKE_BINARY_DIR}/tempMATLAB/")
set(TSGLibOPTC "${CMAKE_CXX_FLAGS}")
set(TSGExampleImport "sys.path.append(\"${TSGPythonPath}\")\n")
set(OMPBasicTest "")

include_directories(${CMAKE_BINARY_DIR}/configured/)
configure_file("${PROJECT_SOURCE_DIR}/Config/tsgVersion.in.hpp"  "${CMAKE_BINARY_DIR}/configured/tsgVersion.hpp")
configure_file("${PROJECT_SOURCE_DIR}/Config/tdrEnableMPI.in.hpp"  "${CMAKE_BINARY_DIR}/configured/tdrEnableMPI.hpp")
configure_file("${PROJECT_SOURCE_DIR}/SparseGrids/GaussPattersonRule.table"  "${CMAKE_BINARY_DIR}/GaussPattersonRule.table" COPYONLY)
if (Tasmanian_ENABLE_PYTHON)
    set(TSGEnableCmakeSyncTest "True")
    configure_file("${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py" "${CMAKE_BINARY_DIR}/TasmanianSG.py")
    configure_file("${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.in.py" "${CMAKE_BINARY_DIR}/example_sparse_grids.py")
    configure_file("${PROJECT_SOURCE_DIR}/Testing/sandbox.py" "${CMAKE_BINARY_DIR}/sandbox.py" COPYONLY)
    configure_file("${PROJECT_SOURCE_DIR}/Testing/testTSG.in.py" "${CMAKE_BINARY_DIR}/testTSG.py")
endif()
if (Tasmanian_ENABLE_MATLAB)
    configure_file("${PROJECT_SOURCE_DIR}/InterfaceMATLAB/tsgGetPaths.in.m" "${CMAKE_BINARY_DIR}/matlab/tsgGetPaths.m")
    file(MAKE_DIRECTORY ${TSGMatlabWork})
    foreach(TasMatlabFile ${tasmat_src})
        configure_file("${PROJECT_SOURCE_DIR}/InterfaceMATLAB/${TasMatlabFile}" "${CMAKE_BINARY_DIR}/matlab/${TasMatlabFile}")
    endforeach()
    message("   pre-install MATLAB folder: addpath('${CMAKE_BINARY_DIR}/matlab/')")
endif()

# stage 2: install folder paths
set(TSGLibPath "${CMAKE_INSTALL_PREFIX}/lib")
set(TSGLibShared "${TSGLibPath}/${CMAKE_SHARED_LIBRARY_PREFIX}${TSGBaseLibName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(TSGLibStatic "${TSGLibPath}/${CMAKE_STATIC_LIBRARY_PREFIX}${TSGStaticLibName}${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(TDRLibShared "${TSGLibPath}/${CMAKE_SHARED_LIBRARY_PREFIX}${TDRBaseLibName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(TDRLibStatic "${TSGLibPath}/${CMAKE_STATIC_LIBRARY_PREFIX}${TDRStaticLibName}${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(FortranLibShared "${TSGLibPath}/${CMAKE_STATIC_LIBRARY_PREFIX}${FortranBaseLibName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(FortranLibStatic "${TSGLibPath}/${CMAKE_STATIC_LIBRARY_PREFIX}${FortranStaticLibName}${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(TSGExecutable "${CMAKE_INSTALL_PREFIX}/bin/tasgrid")
set(TDRExecutable "${CMAKE_INSTALL_PREFIX}/bin/tasdream")
set(TSGPythonPath "${TSGLibPath}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")
set(TSGMatlabPath "${CMAKE_INSTALL_PREFIX}/matlab")
set(TSGMatlabWork "${Tasmanian_MATLAB_WORK_FOLDER}/")
set(TSGLibOPTC "${CMAKE_CXX_FLAGS}")
set(TSGExampleImport "import sys\nsys.path.append(\"${TSGPythonPath}\")\n")
set(OMPBasicTest "")

if (Tasmanian_ENABLE_PYTHON)
    if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    # windows puts the .dll files in bin, as opposed to lib
        set(TSGLibShared "${CMAKE_INSTALL_PREFIX}/bin/${CMAKE_SHARED_LIBRARY_PREFIX}${TSGBaseLibName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
    configure_file("${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py" "${CMAKE_BINARY_DIR}/install/python/TasmanianSG.py")
    configure_file("${PROJECT_SOURCE_DIR}/Examples/example_sparse_grids.in.py" "${CMAKE_BINARY_DIR}/install/examples/example_sparse_grids.py")
else()
    set(TSGPythonPath "PythonNotIncluded")
endif()
if (Tasmanian_ENABLE_MATLAB)
    configure_file("${PROJECT_SOURCE_DIR}/InterfaceMATLAB/tsgGetPaths.in.m" "${CMAKE_BINARY_DIR}/install/matlab/tsgGetPaths.m")
endif()
configure_file("${PROJECT_SOURCE_DIR}/Testing/test_post_install.in.sh" "${CMAKE_BINARY_DIR}/test_post_install.sh")

# machine readable configuration data
set(TasLine1 "export PATH=$PATH:${CMAKE_INSTALL_PREFIX}/bin/")
set(TasLine2 "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib/")
configure_file("${PROJECT_SOURCE_DIR}/Config/TasmanianENVsetup.in.sh" "${CMAKE_BINARY_DIR}/install/config/TasmanianENVsetup.sh")
set(TasLine1 "export C_INCLUDE_PATH=$C_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/")
set(TasLine2 "export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/")
set(TasLine3 "export LIBRARY_PATH=$LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib/")
configure_file("${PROJECT_SOURCE_DIR}/Config/TasmanianDEVsetup.in.sh" "${CMAKE_BINARY_DIR}/install/config/TasmanianDEVsetup.sh")
set(TasExportCMake "${CMAKE_INSTALL_PREFIX}/config/Tasmanian.cmake")
configure_file("${PROJECT_SOURCE_DIR}/Examples/CMakeLists.in.txt" "${CMAKE_BINARY_DIR}/install/examples/CMakeLists.txt")


########################################################################
# Testing
########################################################################
enable_testing()
add_test(Tasgrid tasgrid -test)
add_test(Tasdream tasdream -test)
if (Tasmanian_ENABLE_PYTHON)
    add_test(NAME LibSG COMMAND ${PYTHON_EXECUTABLE} testTSG.py)
endif()
if (Tasmanian_ENABLE_FORTRAN)
    add_test(Fortran fortester)
endif()


########################################################################
# Installing
########################################################################
install(TARGETS ${install_targets}
        EXPORT tas_targets
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")
install(EXPORT tas_targets DESTINATION "config" FILE "Tasmanian.cmake")

install(DIRECTORY "${CMAKE_SOURCE_DIR}/SparseGrids/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "*.windows.*" EXCLUDE
        PATTERN "Example" EXCLUDE
        PATTERN "tsgHiddenExternals.hpp" EXCLUDE)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/DREAM/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "*.windows.*" EXCLUDE
        PATTERN "*.in.*" EXCLUDE
        PATTERN "Example" EXCLUDE)
install(DIRECTORY "${CMAKE_BINARY_DIR}/configured/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")
install(FILES "${CMAKE_SOURCE_DIR}/SparseGrids/TasmanianSparseGrid.h"
        DESTINATION "include"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(FILES "${CMAKE_BINARY_DIR}/install/examples/CMakeLists.txt"
        DESTINATION "examples"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/Examples/"
        DESTINATION "examples"
        FILES_MATCHING PATTERN "*.cpp"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(DIRECTORY "${CMAKE_BINARY_DIR}/install/config/"
        DESTINATION "config"
        FILES_MATCHING PATTERN "Tasmanian*setup.sh"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

if (Tasmanian_ENABLE_PYTHON)
    install(FILES "${CMAKE_BINARY_DIR}/install/python/TasmanianSG.py"
            DESTINATION "${TSGPythonPath}"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    # Create symlink for backward compatibility
    install(CODE "execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${TSGPythonPath} ${CMAKE_INSTALL_PREFIX}/python )" )
    install(FILES "${CMAKE_BINARY_DIR}/install/examples/example_sparse_grids.py"
            DESTINATION "examples"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
endif()

if (Tasmanian_ENABLE_MATLAB)
    install(DIRECTORY "${CMAKE_SOURCE_DIR}/InterfaceMATLAB/"
            DESTINATION "matlab"
            FILES_MATCHING PATTERN "*.m"
            PATTERN "tsgGetPaths.*" EXCLUDE)
    install(FILES "${CMAKE_BINARY_DIR}/install/matlab/tsgGetPaths.m"
            DESTINATION "matlab"
            PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
endif()

if (Tasmanian_ENABLE_FORTRAN)
    install(FILES "${CMAKE_BINARY_DIR}/tasmaniansg.mod"
            DESTINATION "include"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    install(FILES "${CMAKE_SOURCE_DIR}/Examples/example_sparse_grids.f90"
            DESTINATION "examples"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
endif()

########################################################################
# Setup Message
########################################################################
message("")
message("TASMANIAN: summary of all compile options")
if (USE_XSDK_DEFAULTS)
else()
    message(" -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
    foreach(TAS_OPTIONS Tasmanian_ENABLE_OPENMP Tasmanian_ENABLE_BLAS Tasmanian_ENABLE_MPI
                        Tasmanian_ENABLE_MATLAB Tasmanian_ENABLE_PYTHON
                        Tasmanian_ENABLE_CUBLAS Tasmanian_ENABLE_CUDA
                        Tasmanian_ENABLE_FORTRAN
                        Tasmanian_SHARED_LIBRARY Tasmanian_STATIC_LIBRARY)
        if (${TAS_OPTIONS})
            message(" -D ${TAS_OPTIONS}:BOOL=ON")
        else()
            message(" -D ${TAS_OPTIONS}:BOOL=OFF")
        endif()
    endforeach()
    if (Tasmanian_ENABLE_MATLAB)
        message(" -D Tasmanian_MATLAB_WORK_FOLDER:PATH=${Tasmanian_MATLAB_WORK_FOLDER}")
    endif()
endif()
message("")


########################################################################
# Install Message
########################################################################
set(TAS_POST_INSTALL "--------------------------------------------------------------------------------\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}Tasmanian Sparse Grids install complete\n\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}executables:      ${TSGExecutable}\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                  ${TDRExecutable}\n")
if (Tasmanian_SHARED_LIBRARY)
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}shared libraries: ${TSGLibShared}\n")
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                  ${TDRLibShared}\n")
    if (Tasmanian_ENABLE_FORTRAN)
        set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                  ${FortranLibShared}\n")
    endif()
endif()
if (Tasmanian_STATIC_LIBRARY)
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}static libraries: ${TSGLibStatic}\n")
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                  ${TDRLibStatic}\n")
    if (Tasmanian_ENABLE_FORTRAN)
        set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                  ${FortranLibStatic}\n")
    endif()
endif()
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}C++ examples:     ${CMAKE_INSTALL_PREFIX}/examples/example_sparse_grids.cpp\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                  ${CMAKE_INSTALL_PREFIX}/examples/example_dream.cpp\n")
if (Tasmanian_ENABLE_FORTRAN)
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                  ${CMAKE_INSTALL_PREFIX}/examples/example_sparse_grids.f90\n")
endif()
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                  see CMakeLists.txt in ${CMAKE_INSTALL_PREFIX}/examples\n\n")

set(TAS_POST_INSTALL "${TAS_POST_INSTALL}bash export commands (add to your environment, i.e., ~/.bashrc):\n\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}export PATH=$PATH:${CMAKE_INSTALL_PREFIX}/bin/\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib/\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}export C_INCLUDE_PATH=$C_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}export LIBRARY_PATH=$LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib/\n\n")

set(TAS_POST_INSTALL "${TAS_POST_INSTALL}above commands are writte in the following\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}bash source files: ${CMAKE_INSTALL_PREFIX}/config/TasmanianENVsetup.sh\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                   ${CMAKE_INSTALL_PREFIX}/config/TasmanianDEVsetup.sh\n\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}cmake import:             ${CMAKE_INSTALL_PREFIX}/config/Tasmanian.cmake\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}cmake directive:          include(\\\"${CMAKE_INSTALL_PREFIX}/config/Tasmanian.cmake\\\")\n")
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}cmake executable targets: tasgrid          tasdream\n")
if (Tasmanian_ENABLE_FORTRAN)
    if (Tasmanian_STATIC_LIBRARY)
        set(TAS_POST_INSTALL "${TAS_POST_INSTALL}cmake static targets:     libtsg_static    libtdr_static    libtasfor_static\n")
    endif()
    if (Tasmanian_SHARED_LIBRARY)
        set(TAS_POST_INSTALL "${TAS_POST_INSTALL}cmake shared targets:     libtsg_shared    libtdr_shared    libtasfor_shared\n\n")
    endif()
else()
    if (Tasmanian_STATIC_LIBRARY)
        set(TAS_POST_INSTALL "${TAS_POST_INSTALL}cmake static targets:     libtsg_static    libtdr_static\n")
    endif()
    if (Tasmanian_SHARED_LIBRARY)
        set(TAS_POST_INSTALL "${TAS_POST_INSTALL}cmake shared targets:     libtsg_shared    libtdr_shared\n\n")
    endif()
endif()

if (Tasmanian_ENABLE_PYTHON)
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}python module:  ${TSGPythonPath}/TasmanianSG.py\n")
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}python example: ${CMAKE_INSTALL_PREFIX}/examples/example_sparse_grids.py\n")
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                to call the python module\n\n")
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                import sys\n")
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                sys.path.append(\\\"${CMAKE_INSTALL_PREFIX}/python\\\")\n")
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                import TasmanianSG\n\n")
endif()

if (Tasmanian_ENABLE_MATLAB)
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}matlab flles in: ${CMAKE_INSTALL_PREFIX}/matlab/\n")
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                 NOTE: add ${CMAKE_INSTALL_PREFIX}/matlab/ to the MATLAB path\n\n")
    set(TAS_POST_INSTALL "${TAS_POST_INSTALL}                        addpath('${CMAKE_INSTALL_PREFIX}/matlab/')\n")
endif()
set(TAS_POST_INSTALL "${TAS_POST_INSTALL}--------------------------------------------------------------------------------\n")

install(CODE "message(\"${TAS_POST_INSTALL}\")")
install(CODE "file(WRITE ${CMAKE_INSTALL_PREFIX}/config/Tasmanian.log \"${TAS_POST_INSTALL}\")")
install(CODE "message(\"information stored in: ${CMAKE_INSTALL_PREFIX}/config/Tasmanian.log\n\")")
