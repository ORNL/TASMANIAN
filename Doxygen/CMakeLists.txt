find_package(Doxygen REQUIRED OPTIONAL_COMPONENTS dot)

set(DOXYGEN_GENERATE_HTML      "YES")
set(DOXYGEN_DISABLE_INDEX       "NO")
set(DOXYGEN_GENERATE_TREEVIEW  "YES") # left-hand index
set(DOXYGEN_SORT_GROUP_NAMES    "NO") # use the order in the source files (logical order)
set(DOXYGEN_SORT_MEMBER_DOCS    "NO") # use the order in the source files (logical order)

set(DOXYGEN_PREDEFINED "__TASMANIAN_DOXYGEN_SKIP") # indicate sections to skip from documentation

set(DOXYGEN_HTML_EXTRA_STYLESHEET  ${CMAKE_CURRENT_SOURCE_DIR}/tasmanian.css)
set(DOXYGEN_HTML_COLORSTYLE_HUE    "110") # green-ish pages
set(DOXYGEN_HTML_COLORSTYLE_SAT    "100")
set(DOXYGEN_HTML_COLORSTYLE_GAMMA  "240")

set(DOXYGEN_PROJECT_NAME   "Toolkit for Adaptive Stochastic Modeling and Non-Intrusive ApproximatioN: ${PROJECT_NAME} v${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}")
set(DOXYGEN_PROJECT_NUMBER "") # if project number if given, the rendered font is too small

set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_CURRENT_SOURCE_DIR}/../README.md")

if (NOT DOXYGEN_INTERNAL_DOCS)
# the \internal command works for functions and namespaces, but not classes
# have to exclude the classes manually
    set(DOXYGEN_EXCLUDE_SYMBOLS TasGrid::Data2D TasGrid::MultiIndexSet TasGrid::StorageSet TasGrid::CustomTabulated TasGrid::OneDimensionalWrapper)
endif()

doxygen_add_docs(Tasmanian_doxygen
                 ${CMAKE_CURRENT_SOURCE_DIR}/../README.md
                 ${CMAKE_CURRENT_SOURCE_DIR}/Installation.md
                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/tsgIndexManipulator.hpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/tsgIndexSets.hpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/tsgCoreOneDimensional.hpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/tsgOneDimensionalWrapper.hpp
                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                 COMMENT "Tasmanian Doxygen documentation")
