cmake_minimum_required(VERSION 3.10)

if (TARGET Tasmanian::Tasmanian)
    return() # exit if Tasmanian has already been found and defined
endif()

@PACKAGE_INIT@

# https://cmake.org/cmake/help/v3.5/module/CMakePackageConfigHelpers.html#module:CMakePackageConfigHelpers
# seems to indicate that I need to pre-pend PACKAGE_ to the variable names after the @PACKAGE_INIT@ macro
# but this doesn't seem to work, not sure if this is a "relocatable package" (low concern)
include("@Tasmanian_final_install_path@/lib/@CMAKE_PROJECT_NAME@/@CMAKE_PROJECT_NAME@.cmake")

if ("@Tasmanian_ENABLE_MPI@" AND NOT TARGET MPI::MPI_CXX)
    if (NOT MPI_CXX_COMPILER)
        set(MPI_CXX_COMPILER "@MPI_CXX_COMPILER@")
    endif()
    if ("@Tasmanian_ENABLE_FORTRAN@" AND NOT MPI_Fortran_COMPILER)
        set(MPI_Fortran_COMPILER "@MPI_Fortran_COMPILER@")
    endif()
    find_package(MPI REQUIRED)
endif()

add_executable(Tasmanian::tasgrid IMPORTED)
set_property(TARGET Tasmanian::tasgrid PROPERTY IMPORTED_LOCATION "@Tasmanian_final_install_path@/bin/tasgrid${CMAKE_EXECUTABLE_SUFFIX_CXX}")

add_library(Tasmanian::Tasmanian INTERFACE IMPORTED GLOBAL) # master target
target_link_libraries(Tasmanian::Tasmanian INTERFACE Tasmanian_master)

if (@BUILD_SHARED_LIBS@)
    set(Tasmanian_SHARED_FOUND "ON")
    add_library(Tasmanian::shared INTERFACE IMPORTED GLOBAL) # master target
    target_link_libraries(Tasmanian::shared INTERFACE Tasmanian_master)
else()
    set(Tasmanian_STATIC_FOUND "ON")
    add_library(Tasmanian::static INTERFACE IMPORTED GLOBAL) # master target
    target_link_libraries(Tasmanian::static INTERFACE Tasmanian_master)

    if (@Tasmanian_ENABLE_CUDA@)
    # Since Tasmanian does not transitively include <cuda.h> and since all CUDA calls are wrapped in CXX API,
    # projects do not require cuda language to link to Tasmanian; however, CMake adds the extraneous dependence.
    # If Tasmanian was build with CUDA and if the user has not explicitly enabled the CUDA language,
    # then overwrite the CMake generated extraneous CUDA requirements and link with the CXX compiler only.
    # This hack is not necessary when building shared libraries.
        if (NOT CMAKE_CUDA_COMPILER)
            set_target_properties(Tasmanian_libsparsegrid PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX")
            set_target_properties(Tasmanian_libsparsegrid PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX")
            get_target_property(_TasLibs Tasmanian_libsparsegrid INTERFACE_LINK_LIBRARIES)
            if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
                set_target_properties(Tasmanian_libsparsegrid PROPERTIES INTERFACE_LINK_LIBRARIES "${_TasLibs};@Tasmanian_cudaruntime@")
            else()
                set_target_properties(Tasmanian_libsparsegrid PROPERTIES INTERFACE_LINK_LIBRARIES "${_TasLibs};@Tasmanian_cudaruntime@;dl")
            endif()
            unset(_TasLibs)
        endif()
    endif()
endif()

if (@Tasmanian_ENABLE_HIP@)
    if (NOT "@Tasmanian_ROCM_ROOT@" STREQUAL "")
        list(APPEND CMAKE_PREFIX_PATH "@Tasmanian_ROCM_ROOT@")
    endif()
    list(APPEND CMAKE_PREFIX_PATH "@Tasmanian_hiproot@")
    foreach(_tsg_package hip rocblas rocsparse rocsolver)
        find_package(${_tsg_package} REQUIRED)
    endforeach()
endif()

if (TARGET tasmanian) # swig Tasmanian Fortran
    add_library(Tasmanian_libfortran03 INTERFACE IMPORTED GLOBAL)
    target_link_libraries(Tasmanian_libfortran03 INTERFACE tasmanian)
endif()

if (TARGET Tasmanian_libfortran03)
    add_library(Tasmanian::Fortran INTERFACE IMPORTED GLOBAL)
    target_link_libraries(Tasmanian::Fortran INTERFACE Tasmanian_libfortran90 Tasmanian_libfortran03)

    if (TARGET Tasmanian_libfortranmpi03)
        target_link_libraries(Tasmanian::Fortran INTERFACE Tasmanian_libfortranmpi03)
    endif()

    set(Tasmanian_FORTRAN_FOUND "ON")
endif()

# export the python path so other projects can configure python scripts
if (@Tasmanian_ENABLE_PYTHON@)
    set_and_check(Tasmanian_PYTHONPATH "@Tasmanian_final_install_path@/share/Tasmanian/python/")
    set(Tasmanian_PYTHON_FOUND "ON")
endif()

# export the MATLAB paths so other projects can write files directly to MATLAB
if (NOT "@Tasmanian_MATLAB_WORK_FOLDER@" STREQUAL "")
    set_and_check(Tasmanian_MATLAB_WORK_FOLDER "@Tasmanian_MATLAB_WORK_FOLDER@")
    set_and_check(Tasmanian_MATLABPATH "@Tasmanian_final_install_path@/share/Tasmanian/matlab/")
    set(Tasmanian_MATLAB_FOUND "ON")
endif()

set(Tasmanian_OPENMP_FOUND "@Tasmanian_ENABLE_OPENMP@")
set(Tasmanian_BLAS_FOUND   "@Tasmanian_ENABLE_BLAS@")
set(Tasmanian_MPI_FOUND    "@Tasmanian_ENABLE_MPI@")
set(Tasmanian_CUDA_FOUND   "@Tasmanian_ENABLE_CUDA@")
set(Tasmanian_HIP_FOUND    "@Tasmanian_ENABLE_HIP@")
set(Tasmanian_DPCPP_FOUND  "@Tasmanian_ENABLE_DPCPP@")
set(Tasmanian_MAGMA_FOUND  "@Tasmanian_ENABLE_MAGMA@")

# component GPU uses either GPU backend (CUDA or HIP)
if (Tasmanian_CUDA_FOUND OR Tasmanian_HIP_FOUND OR Tasmanian_DPCPP_FOUND)
    set(Tasmanian_GPU_FOUND "ON")
endif()

# write component info
foreach(_tsg_comp ${Tasmanian_FIND_COMPONENTS})
    if (Tasmanian_${_tsg_comp}_FOUND)
        message(STATUS "Tasmanian component ${_tsg_comp}: found")
    else()
        if (Tasmanian_FIND_REQUIRED_${_tsg_comp})
            message(WARNING "Tasmanian required component ${_tsg_comp}: missing (error)")
        else()
            message(STATUS "Tasmanian optional component ${_tsg_comp}: missing")
        endif()
    endif()
endforeach()
unset(_tsg_comp)

check_required_components(Tasmanian)
