cmake_minimum_required(VERSION @CMAKE_MAJOR_VERSION@.@CMAKE_MINOR_VERSION@)
cmake_policy(VERSION @CMAKE_MAJOR_VERSION@.@CMAKE_MINOR_VERSION@)
project(TasmanianTesting VERSION @Tasmanian_VERSION_MAJOR@.@Tasmanian_VERSION_MINOR@.@Tasmanian_VERSION_PATCH@ LANGUAGES @Tasmanian_langs@)
enable_testing()

message(STATUS "Tasmanian post-installation testing")

# the following find_package() command will help us locate an existing Tasmanian installation.
find_package(Tasmanian @Tasmanian_VERSION_MAJOR@.@Tasmanian_VERSION_MINOR@.@Tasmanian_VERSION_PATCH@ PATHS "@Tasmanian_final_install_path@"
             REQUIRED @Tasmanian_components@)

add_subdirectory("@Tasmanian_final_install_path@/share/Tasmanian/examples" examples_cxx)

add_test(SparseGrids   "${CMAKE_CURRENT_BINARY_DIR}/examples_cxx/example_sparse_grids"  -fast)
add_test(DREAM         "${CMAKE_CURRENT_BINARY_DIR}/examples_cxx/example_dream"         -fast)
if (Tasmanian_FORTRAN_FOUND)
    add_test(Fortran      "${CMAKE_CURRENT_BINARY_DIR}/examples_cxx/example_sparse_grids_fortran"     -fast)
endif()

if (Tasmanian_PYTHON_FOUND)
    add_test(Python::SparseGrids  @PYTHON_EXECUTABLE@ "@Tasmanian_final_install_path@/share/Tasmanian/examples/example_sparse_grids.py"  -fast)
    add_test(Python::DREAM        @PYTHON_EXECUTABLE@ "@Tasmanian_final_install_path@/share/Tasmanian/examples/example_dream.py"         -fast)
endif()
