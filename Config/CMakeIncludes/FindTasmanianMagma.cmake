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

find_package_handle_standard_args(TasmanianMagma DEFAULT_MSG Tasmanian_magma Tasmanian_magma_h)

unset(Tasmanian_magma_nodefault)
unset(MAGMA_ROOT)
