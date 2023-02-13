########################################################################
# FindTasmanianMAGMA.cmake module
########################################################################
#
# searches for libmagma and libmagma_sparse
# and the associated include folder
# uses paths defined by MAGMA_ROOT and Tasmanian_MAGMA_ROOT

# Defines variables:
#    Tasmanian_magmalibs: list of libmagma and libmagma_sparse
#    Tasmanian_magma_h: path to the folder containing magma.h
#    Tasmanian_magma: the libmagma library

if (NOT MAGMA_ROOT)
    set(MAGMA_ROOT "$ENV{MAGMA_ROOT}")
endif()

if (Tasmanian_MAGMA_DOWNLOAD)
    include (FetchContent)

    if (Tasmanian_ENABLE_CUDA)
        set(MAGMA_ENABLE_CUDA ON CACHE BOOL "pass into MAGMA")
        set(USE_FORTRAN OFF CACHE BOOL "pass into MAGMA")
    endif()

    FetchContent_Declare(TasmanianMAGMA
                         URL "http://icl.utk.edu/projectsfiles/magma/downloads/magma-2.7.0.tar.gz"
                         URL_HASH "SHA256=fda1cbc4607e77cacd8feb1c0f633c5826ba200a018f647f1c5436975b39fd18"
                         DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/_magma_download/
                         SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/_magma_download/magma
                         PATCH_COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/Config/patch_magma.sh
                         )
    FetchContent_MakeAvailable(TasmanianMAGMA)

    add_library(Tasmanian::MAGMA INTERFACE IMPORTED)
    target_link_libraries(Tasmanian::MAGMA INTERFACE magma magma_sparse)
    target_include_directories(Tasmanian::MAGMA INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/_magma_download/magma/include")
    target_include_directories(Tasmanian::MAGMA INTERFACE "${tasmanianmagma_BINARY_DIR}/include")
    if (BUILD_SHARED_LIBS)
        set(Tasmanian_magma "${CMAKE_INSTALL_PREFIX}/lib/libmagma.so;${CMAKE_INSTALL_PREFIX}/lib/libmagma_sparse.so;")
    else()
        set(Tasmanian_magma "${CMAKE_INSTALL_PREFIX}/lib/libmagma.a;${CMAKE_INSTALL_PREFIX}/lib/libmagma_sparse.a;")
    endif()
    set(Tasmanian_magma_h "${CMAKE_INSTALL_PREFIX}/include/")

else()

    set(Tasmanian_MAGMA_ROOT "${MAGMA_ROOT}" CACHE PATH "The root folder for the MAGMA installation, e.g., containing lib and include folders")

    if (Tasmanian_MAGMA_ROOT)
        set(Tasmanian_magma_nodefault "NO_DEFAULT_PATH")
    endif()

    Tasmanian_find_libraries(REQUIRED magma
                             OPTIONAL magma_sparse
                             PREFIX ${Tasmanian_MAGMA_ROOT}
                             LIST magmalibs
                             ${Tasmanian_magma_nodefault})

    if (Tasmanian_magma) # if the library is missing, there's no point in searching for the header
        Tasmanian_find_header(FILE "magma.h"
                              RESULT magma_h
                              ROOT ${Tasmanian_MAGMA_ROOT}
                              HINT ${Tasmanian_magma}
                              NO_DEFAULT_PATH)
    endif()

    add_library(Tasmanian::MAGMA INTERFACE IMPORTED)
    target_link_libraries(Tasmanian::MAGMA INTERFACE ${Tasmanian_magma})
    target_include_directories(Tasmanian::MAGMA INTERFACE ${Tasmanian_magma_h})

endif()

find_package_handle_standard_args(TasmanianMagma DEFAULT_MSG Tasmanian_magma Tasmanian_magma_h)

unset(Tasmanian_magma_nodefault)
unset(MAGMA_ROOT)
