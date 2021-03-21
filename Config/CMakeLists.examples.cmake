cmake_minimum_required(VERSION @CMAKE_MAJOR_VERSION@.@CMAKE_MINOR_VERSION@)
cmake_policy(VERSION @CMAKE_MAJOR_VERSION@.@CMAKE_MINOR_VERSION@)
project(Tasmanian_Examples VERSION @Tasmanian_VERSION_MAJOR@.@Tasmanian_VERSION_MINOR@.@Tasmanian_VERSION_PATCH@ LANGUAGES @Tasmanian_langs@)

# the following find_package() command will help us locate an existing Tasmanian installation.
find_package(Tasmanian @Tasmanian_VERSION_MAJOR@.@Tasmanian_VERSION_MINOR@.@Tasmanian_VERSION_PATCH@ PATHS "@Tasmanian_final_install_path@"
             REQUIRED @Tasmanian_components@)

# The REQUIRED clause above gives a list of all available modules
# those correspond to the Tasmanian CMake options.
# The SHARED and STATIC modules correspond to the shared and static variants of the libraries.
# Each module also comes with a Boolean variable: Tasmanian_<module-name>_FOUND

# A project can skip PATHS, and use either:
#  CMAKE_PREFIX_PATH variable in CMake that includes "@CMAKE_INSTALL_PREFIX@/lib/"
#                    for example using -DCMAKE_PREFIX_PATH="@CMAKE_INSTALL_PREFIX@/lib/"
#  Tasmanian_DIR or Tasmanian_ROOT environment variables (depending on the version of CMake)

# Tasmanian::Tasmanian will be the name of the IMPORTED C++ target.
# Tasmanian::Fortran will be the name of the IMPORTED Fortran target.

# Additional targets:
#  Tasmanian::tasgrid will point to the executable ./tasgrid@CMAKE_EXECUTABLE_SUFFIX_CXX@

# Additional variables (if the corresponding options have been enabled):
#  Tasmanian_PYTHONPATH is the path to the python scripts
#  Tasmanian_MATLABPATH is the path to the MATLAB scripts
#  Tasmanian_MATLAB_WORK_FOLDER


add_executable(example_sparse_grids  example_sparse_grids_01.cpp
                                     example_sparse_grids_02.cpp
                                     example_sparse_grids_03.cpp
                                     example_sparse_grids_04.cpp
                                     example_sparse_grids_05.cpp
                                     example_sparse_grids_06.cpp
                                     example_sparse_grids_07.cpp
                                     example_sparse_grids_08.cpp
                                     example_sparse_grids_09.cpp
                                     example_sparse_grids_10.cpp
                                     example_sparse_grids.cpp)

add_executable(example_dream         example_dream_01.cpp
                                     example_dream_02.cpp
                                     example_dream_03.cpp
                                     example_dream_04.cpp
                                     example_dream_05.cpp
                                     example_dream.cpp)

target_link_libraries(example_sparse_grids  Tasmanian::Tasmanian)
target_link_libraries(example_dream         Tasmanian::Tasmanian)

# if the Fortran component was found
# can also use "if (TARGET Tasmanian::Fortran)"
if (Tasmanian_FORTRAN_FOUND)
    add_executable(example_sparse_grids_fortran  example_sparse_grids_01.f90
                                                 example_sparse_grids_02.f90
                                                 example_sparse_grids_03.f90
                                                 example_sparse_grids_04.f90
                                                 example_sparse_grids.f90)
    target_link_libraries(example_sparse_grids_fortran  Tasmanian::Fortran)
    # note that as of 7.1 Tasmanian::Fortran is not equivalent to Tasmanian::Tasmanian

    if (@Tasmanian_ifort_compiler@)
        set_target_properties(example_sparse_grids_fortran PROPERTIES LINKER_LANGUAGE Fortran)
    else()
        set_target_properties(example_sparse_grids_fortran PROPERTIES LINKER_LANGUAGE CXX)
    endif()
    # the LINKER_LANGUAGE property is required by some compilers, e.g., Intel ifort
endif()

# Tasmanian also includes the executable target Tasmanian::tasgrid
add_custom_command(TARGET example_sparse_grids PRE_BUILD COMMAND Tasmanian::tasgrid -v)
