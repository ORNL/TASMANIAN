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
set(DOXYGEN_HTML_COLORSTYLE_GAMMA  "20")


set(IMAGE_PATH "images")

set(DOXYGEN_FORMULA_FONTSIZE "14")

set(DOXYGEN_PROJECT_NAME   "Toolkit for Adaptive Stochastic Modeling and Non-Intrusive ApproximatioN: ${PROJECT_NAME} v${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}")
set(DOXYGEN_PROJECT_NUMBER "") # if project number is given, the rendered font is too small

set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "README.md")

set(DOXYGEN_EXAMPLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../") # finding the examples for all modules requires the use of the root dir

if (NOT DOXYGEN_INTERNAL_DOCS)
# the \internal command works for functions and namespaces, but not classes
# have to exclude the classes manually
    set(DOXYGEN_EXCLUDE_SYMBOLS TasGrid::Data2D TasGrid::MultiIndexSet TasGrid::StorageSet TasGrid::CustomTabulated TasGrid::OneDimensionalWrapper
        TasGrid::TasSparse::SparseMatrix TasGrid::CudaEngine TasGrid::AccelerationDomainTransform TasGrid::TableGaussPatterson
        TasGrid::NodeData TasGrid::TensorData TasGrid::CacheLagrange TasGrid::MultiIndexManipulations::ProperWeights
        TasGrid::Utils::Wrapper2D TasGrid::SimpleConstructData TasGrid::DynamicConstructorDataGlobal
        TasGrid::Optimizer::CurrentNodes TasGrid::Optimizer::HasDerivative TasGrid::Optimizer::OptimizerResult
        TasGrid::VectorToStreamBuffer)
endif()

doxygen_add_docs(Tasmanian_doxygen
                 README.md
                 Doxygen/Installation.md
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
                 SparseGrids/tsgHiddenExternals.hpp
                 SparseGrids/tsgAcceleratedDataStructures.hpp
                 SparseGrids/tsgCacheLagrange.hpp
                 SparseGrids/tsgCoreOneDimensional.hpp
                 SparseGrids/tsgCoreOneDimensional.hpp
                 SparseGrids/tsgOneDimensionalWrapper.hpp
                 SparseGrids/tsgLinearSolvers.hpp
                 SparseGrids/tsgDConstructGridGlobal.hpp
                 SparseGrids/tsgSequenceOptimizer.hpp
                 SparseGrids/tsgSequenceOptimizer.cpp
                 SparseGrids/tsgGridCore.hpp
                 DREAM/tsgDreamEnumerates.hpp
                 DREAM/tsgDreamState.hpp
                 DREAM/tsgDreamSample.hpp
                 DREAM/tsgDreamSampleGrid.hpp
                 DREAM/tsgDreamSamplePosterior.hpp
                 DREAM/tsgDreamSamplePosteriorGrid.hpp
                 DREAM/tsgDreamCoreRandom.hpp
                 DREAM/tsgDreamCorePDF.hpp
                 DREAM/tsgDreamLikelihoodCore.hpp
                 DREAM/tsgDreamLikelyGaussian.hpp
                 DREAM/tsgDreamInternalBlas.hpp
                 Addons/tsgConstructSurrogate.hpp
                 Addons/tsgCandidateManager.hpp
                 Addons/tsgMPIScatterGrid.hpp
                 Addons/tsgAddonsCommon.hpp
                 SparseGrids/Examples/example_sparse_grids.cpp
                 DREAM/Examples/example_dream.cpp
                 DREAM/Examples/example_dream_01.cpp
                 DREAM/Examples/example_dream_02.cpp
                 DREAM/Examples/example_dream_03.cpp
                 DREAM/Examples/example_dream_04.cpp
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../
                 COMMENT "Tasmanian Doxygen documentation")

# doxygen_add_docs() doesn't add the target to "all"
# need to manually add dependence to ensure the docs are build before "install"
# adding Tasmanian_doxygen to both sparse grids libraries forces the docs to be build first
# even if parallel build is used (e.g., make -j) this will keep the c++ and doxygen outputs separate
if (TARGET Tasmanian_libsparsegrid_shared)
    add_dependencies(Tasmanian_libsparsegrid_shared Tasmanian_doxygen)
endif()
if (TARGET Tasmanian_libsparsegrid_static)
    add_dependencies(Tasmanian_libsparsegrid_static Tasmanian_doxygen)
endif()

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION "share/Tasmanian/docs")
