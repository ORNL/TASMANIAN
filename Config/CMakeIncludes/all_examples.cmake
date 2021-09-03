add_executable(Tasmanian_example_sparse_grids SparseGrids/Examples/example_sparse_grids_01.cpp
                                              SparseGrids/Examples/example_sparse_grids_02.cpp
                                              SparseGrids/Examples/example_sparse_grids_03.cpp
                                              SparseGrids/Examples/example_sparse_grids_04.cpp
                                              SparseGrids/Examples/example_sparse_grids_05.cpp
                                              SparseGrids/Examples/example_sparse_grids_06.cpp
                                              SparseGrids/Examples/example_sparse_grids_07.cpp
                                              SparseGrids/Examples/example_sparse_grids_08.cpp
                                              SparseGrids/Examples/example_sparse_grids_09.cpp
                                              SparseGrids/Examples/example_sparse_grids_10.cpp
                                              SparseGrids/Examples/example_sparse_grids.cpp)

set_target_properties(Tasmanian_example_sparse_grids PROPERTIES OUTPUT_NAME "example_sparse_grids")
target_link_libraries(Tasmanian_example_sparse_grids Tasmanian_master)


add_executable(Tasmanian_example_dream DREAM/Examples/example_dream_01.cpp
                                       DREAM/Examples/example_dream_02.cpp
                                       DREAM/Examples/example_dream_03.cpp
                                       DREAM/Examples/example_dream_04.cpp
                                       DREAM/Examples/example_dream_05.cpp
                                       DREAM/Examples/example_dream.cpp)

set_target_properties(Tasmanian_example_dream PROPERTIES OUTPUT_NAME "example_dream")
target_link_libraries(Tasmanian_example_dream  Tasmanian_master)


if (Tasmanian_ENABLE_FORTRAN)
    if (NOT Tasmanian_DISABLE_F90)
        add_executable(Tasmanian_example_sparse_grids_f90 InterfaceFortran/Examples/example_sparse_grids.f90)

        set_target_properties(Tasmanian_example_sparse_grids_f90 PROPERTIES OUTPUT_NAME "example_sparse_grids_f90")
        target_link_libraries(Tasmanian_example_sparse_grids_f90 Tasmanian_libfortran90)

        Tasmanian_set_fortran_linker(TARGET Tasmanian_example_sparse_grids_f90)
    endif()

    # adding the Fortran 2003 examples
    add_executable(Tasmanian_example_sparse_grids_fortran InterfaceSwig/FortranExamples/example_sparse_grids.f90
                                                          InterfaceSwig/FortranExamples/example_sparse_grids_01.f90
                                                          InterfaceSwig/FortranExamples/example_sparse_grids_02.f90
                                                          InterfaceSwig/FortranExamples/example_sparse_grids_03.f90
                                                          InterfaceSwig/FortranExamples/example_sparse_grids_04.f90
                                                          )
    set_target_properties(Tasmanian_example_sparse_grids_fortran PROPERTIES OUTPUT_NAME "example_sparse_grids_fortran")
    target_link_libraries(Tasmanian_example_sparse_grids_fortran Tasmanian_libfortran03)

    Tasmanian_set_fortran_linker(TARGET Tasmanian_example_sparse_grids_fortran)
endif()

get_filename_component(_tsg_cpath "${CMAKE_CURRENT_BINARY_DIR}" NAME)
foreach(_tsg_example Tasmanian_example_sparse_grids Tasmanian_example_dream Tasmanian_example_sparse_grids_f90 Tasmanian_example_sparse_grids_fortran)
    if (TARGET ${_tsg_example})
        Tasmanian_rpath_target(TARGET ${_tsg_example} USE_CURRENT
                               COMPONENTS ${_tsg_cpath}/SparseGrids ${_tsg_cpath}/DREAM
                                          ${_tsg_cpath}/Fortran ${_tsg_cpath}/Fortran03)
    endif()
endforeach()
unset(_tsg_cpath)
unset(_tsg_example)
