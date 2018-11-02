find_package(Doxygen REQUIRED OPTIONAL_COMPONENTS dot)

set(DOXYGEN_GENERATE_HTML      "YES")
set(DOXYGEN_DISABLE_INDEX       "NO")
set(DOXYGEN_GENERATE_TREEVIEW  "YES") # left-hand index
set(DOXYGEN_SORT_GROUP_NAMES    "NO") # use the order in the source files (logical order)
set(DOXYGEN_SORT_MEMBER_DOCS    "NO") # use the order in the source files (logical order)
set(DOXYGEN_INTERNAL_DOCS      "YES") # silence eventually, especially for a release

set(DOXYGEN_PREDEFINED "__TASMANIAN_DOXYGEN_SKIP") # indicate sections to skip from documentation

set(DOXYGEN_HTML_EXTRA_STYLESHEET  ${CMAKE_CURRENT_SOURCE_DIR}/tasmanian.css)
set(DOXYGEN_HTML_COLORSTYLE_HUE    "110") # green-ish pages
set(DOXYGEN_HTML_COLORSTYLE_SAT    "100")
set(DOXYGEN_HTML_COLORSTYLE_GAMMA  "240")

set(DOXYGEN_PROJECT_NAME   "Toolkit for Adaptive Stochastic Modeling and Non-Intrusive ApproximatioN: ${PROJECT_NAME} v${Tasmanian_VERSION_MAJOR}.${Tasmanian_VERSION_MINOR}")
set(DOXYGEN_PROJECT_NUMBER "") # if project number if given, the rendered font is too small

set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_CURRENT_SOURCE_DIR}/../README.md")

doxygen_add_docs(Tasmanian_docs
                 ${CMAKE_CURRENT_SOURCE_DIR}/../README.md
                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/tsgIndexManipulator.hpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/../SparseGrids/tsgIndexSets.hpp
                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                 COMMENT "Tasmanian Doxygen documentation")

add_custom_target(Tasmanian_doxygen_documentation ALL DEPENDS Tasmanian_docs)
