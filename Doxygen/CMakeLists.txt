find_package(Doxygen REQUIRED OPTIONAL_COMPONENTS dot)

set(DOXYGEN_GENERATE_HTML      "YES")
set(DOXYGEN_DISABLE_INDEX       "NO")
set(DOXYGEN_GENERATE_TREEVIEW   "YES") # left-hand index
set(DOXYGEN_SORT_GROUP_NAMES    "NO") # disable automatic sorting
set(DOXYGEN_SORT_MEMBER_DOCS    "NO") # use the order in the source files (logical order)
set(DOXYGEN_SORT_BRIEF_DOCS     "NO") # which is more logical than the alphabetical one

# skip sections marked with __TASMANIAN_DOXYGEN_SKIP to avoid extraneous warnings
# document all options regardless of the build, mostly affects the internal API
set(DOXYGEN_PREDEFINED "__TASMANIAN_DOXYGEN_SKIP;Tasmanian_ENABLE_BLAS;Tasmanian_ENABLE_CUDA;Tasmanian_ENABLE_MPI")

set(DOXYGEN_HTML_EXTRA_STYLESHEET  ${CMAKE_CURRENT_SOURCE_DIR}/tasmanian.css)
set(DOXYGEN_HTML_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/header.html")
set(DOXYGEN_HTML_FOOTER "${CMAKE_CURRENT_SOURCE_DIR}/footer.html")
set(DOXYGEN_HTML_COLORSTYLE_HUE    "175") # green-ish pages
set(DOXYGEN_HTML_COLORSTYLE_SAT    "50")
set(DOXYGEN_HTML_COLORSTYLE_GAMMA  "90")


set(DOXYGEN_FORMULA_FONTSIZE "14")
set(DOXYGEN_USE_MATHJAX "YES")

set(DOXYGEN_PROJECT_NAME   "Toolkit for Adaptive Stochastic Modeling and Non-Intrusive ApproximatioN: ${PROJECT_NAME} v${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}${Tasmanian_version_comment}")
set(DOXYGEN_PROJECT_NUMBER "") # if project number is given, the rendered font is too small

set(DOXYGEN_EXAMPLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../") # finding the examples for all modules requires the use of the root dir

if (NOT DOXYGEN_INTERNAL_DOCS)
# the \internal command works for functions and namespaces, but not classes
# have to exclude the classes manually
    set(DOXYGEN_EXCLUDE_SYMBOLS TasGrid::Data2D TasGrid::MultiIndexSet TasGrid::StorageSet TasGrid::CustomTabulated TasGrid::OneDimensionalWrapper
        TasGrid::TasSparse::WaveletBasisMatrix TasGrid::GpuEngine TasGrid::GpuVector TasGrid::AccelerationDomainTransform
        TasGrid::AccelerationContext TasGrid::TableGaussPatterson TasGrid::InternalSyclQueue
        TasBLAS TasGrid::TasmanianDenseSolver TasGrid::TasmanianTridiagonalSolver TasGrid::Utils
        TasGrid::NodeData TasGrid::TensorData TasGrid::CacheLagrange TasGrid::MultiIndexManipulations::ProperWeights
        TasGrid::HierarchyManipulations::SplitDirections
        TasGrid::Utils::Wrapper2D TasGrid::SimpleConstructData TasGrid::DynamicConstructorDataGlobal
        TasGrid::Optimizer::CurrentNodes TasGrid::Optimizer::HasDerivative TasGrid::Optimizer::OptimizerResult
        TasGrid::VectorToStreamBuffer TasGrid::CandidateManager TasGrid::CompleteStorage
        TasGrid::IO::mode_ascii_type TasGrid::IO::mode_binary_type TasGrid::HandleDeleter
        TasGrid::AccHandle::Cublas TasGrid::AccHandle::Cusparse TasGrid::AccHandle::Cusolver
        TasGrid::AccHandle::Rocblas TasGrid::AccHandle::Rocsparse TasGrid::AccHandle::Syclqueue
        TasGrid::tsg_gpu_selector TasGrid::int_gpu_lapack
        TasDREAM::IO TasDREAM::Utils)
    set(DOXYGEN_PREDEFINED "__TASMANIAN_DOXYGEN_SKIP;__TASMANIAN_DOXYGEN_SKIP_INTERNAL;Tasmanian_ENABLE_BLAS;Tasmanian_ENABLE_CUDA;Tasmanian_ENABLE_DPCPP;Tasmanian_ENABLE_HIP;Tasmanian_ENABLE_MPI")
endif()

### Add the git-hash to the main page, suppress at release (the stable page needs no hash)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/../README.md" Tasmanian_readme)
string(REPLACE "# Tasmanian" "# Tasmanian\n    * *Generated from git commit ${Tasmanian_git_hash}*" Tasmanian_readme ${Tasmanian_readme})
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/README.md" "${Tasmanian_readme}")
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_CURRENT_BINARY_DIR}/README.md")

doxygen_add_docs(Tasmanian_doxygen
                 ${CMAKE_CURRENT_BINARY_DIR}/README.md
                 Doxygen/Installation.md
                 Doxygen/InterfacePython.md
                 Doxygen/InterfaceLibEnsemble.md
                 Doxygen/InterfaceCLI.md
                 Doxygen/InterfaceMATLAB.md
                 Doxygen/InterfaceFortran90.md
                 Doxygen/InterfaceFortran2003.md
                 Doxygen/CustomRuleFileFormat.md
                 Config/Tasmanian.hpp
                 SparseGrids/TasmanianSparseGrid.hpp
                 DREAM/TasmanianDREAM.hpp
                 Addons/TasmanianAddons.hpp
                 SparseGrids/tsgEnumerates.hpp
                 SparseGrids/tsgUtils.hpp
                 SparseGrids/tsgMathUtils.hpp
                 SparseGrids/tsgIndexSets.hpp
                 SparseGrids/tsgIndexManipulator.hpp
                 SparseGrids/tsgIndexManipulator.cpp
                 SparseGrids/tsgHierarchyManipulator.hpp
                 SparseGrids/tsgIOHelpers.hpp
                 SparseGrids/tsgAcceleratedHandles.hpp
                 SparseGrids/tsgAcceleratedDataStructures.hpp
                 SparseGrids/tsgCudaLoadStructures.hpp
                 SparseGrids/tsgCacheLagrange.hpp
                 SparseGrids/tsgCoreOneDimensional.hpp
                 SparseGrids/tsgOneDimensionalWrapper.hpp
                 SparseGrids/tsgHardCodedTabulatedRules.hpp
                 SparseGrids/tsgLinearSolvers.hpp
                 SparseGrids/tsgDConstructGridGlobal.hpp
                 SparseGrids/tsgSequenceOptimizer.hpp
                 SparseGrids/tsgSequenceOptimizer.cpp
                 SparseGrids/tsgRuleWavelet.hpp
                 SparseGrids/tsgRuleLocalPolynomial.hpp
                 SparseGrids/tsgGridCore.hpp
                 SparseGrids/tsgGridGlobal.hpp
                 SparseGrids/tsgGridSequence.hpp
                 SparseGrids/tsgGridLocalPolynomial.hpp
                 SparseGrids/tsgGridWavelet.hpp
                 SparseGrids/tsgGridFourier.hpp
                 DREAM/tsgDreamState.hpp
                 DREAM/tsgDreamSample.hpp
                 DREAM/tsgDreamEnumerates.hpp
                 DREAM/tsgDreamCoreRandom.hpp
                 DREAM/tsgDreamCorePDF.hpp
                 DREAM/tsgDreamLikelihoodCore.hpp
                 DREAM/tsgDreamLikelyGaussian.hpp
                 Addons/tsgExoticQuadrature.hpp
                 Addons/tsgConstructSurrogate.hpp
                 Addons/tsgCandidateManager.hpp
                 Addons/tsgLoadUnstructuredPoints.hpp
                 Addons/tsgMPIConstructGrid.hpp
                 Addons/tsgMPISampleDream.hpp
                 Addons/tsgMPIScatterDream.hpp
                 Addons/tsgMPIScatterGrid.hpp
                 Addons/tsgLoadNeededValues.hpp
                 Addons/tsgAddonsCommon.hpp
                 InterfaceTPL/tsgTPLWrappers.hpp
                 InterfaceTPL/tsgBlasWrappers.hpp
                 InterfaceTPL/tsgGpuWrappers.hpp
                 SparseGrids/Examples/example_sparse_grids.cpp
                 SparseGrids/Examples/example_sparse_grids_01.cpp
                 SparseGrids/Examples/example_sparse_grids_02.cpp
                 SparseGrids/Examples/example_sparse_grids_03.cpp
                 SparseGrids/Examples/example_sparse_grids_04.cpp
                 SparseGrids/Examples/example_sparse_grids_05.cpp
                 SparseGrids/Examples/example_sparse_grids_06.cpp
                 SparseGrids/Examples/example_sparse_grids_07.cpp
                 SparseGrids/Examples/example_sparse_grids_08.cpp
                 SparseGrids/Examples/example_sparse_grids_09.cpp
                 SparseGrids/Examples/example_sparse_grids_10.cpp
                 DREAM/Examples/example_dream.cpp
                 DREAM/Examples/example_dream_01.cpp
                 DREAM/Examples/example_dream_02.cpp
                 DREAM/Examples/example_dream_03.cpp
                 DREAM/Examples/example_dream_04.cpp
                 DREAM/Examples/example_dream_05.cpp
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../
                 COMMENT "Tasmanian Doxygen documentation")

# doxygen_add_docs() doesn't add the target to "all"
# need to manually add dependence to ensure the docs are build before "install"
# adding Tasmanian_doxygen to both sparse grids libraries forces the docs to be build first
# even if parallel build is used (e.g., make -j) this will keep the c++ and doxygen outputs separate
add_dependencies(Tasmanian_libsparsegrid  Tasmanian_doxygen)


install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION "share/Tasmanian/docs")
