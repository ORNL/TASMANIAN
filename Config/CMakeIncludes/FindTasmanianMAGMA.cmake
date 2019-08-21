########################################################################
# FindTasmanianMAGMA.cmake module
########################################################################
#
# searches for libmagma.{a,so,dylib} and libmagma_sparse.{a,so,dylib}
# and the associated include folder
#
# using (optional) variables:
#     Tasmanian_MAGMA_ROOT_DIR
#     MAGMA_ROOT_DIR
#     Tasmanian_MAGMA_LIBRARIES
#     Tasmanian_MAGMA_SHARED_LIBRARIES
#     Tasmanian_MAGMA_INCLUDE_DIRS
#
# defines variables:
#     Tasmanian_MAGMA_FOUND: ON/OFF
#         indicates success/fail
#     Tasmanian_MAGMA_SHARED_FOUND: ON/OFF
#         ON if libraries with not ".a" suffix are found
#
#     Tasmanian_MAGMA_LIBRARIES: list
#         list of libraries to use for static linking
#         (the libraries could be shared libraries w/o warning)
#
#     Tasmanian_MAGMA_SHARED_LIBRARIES: list
#         list of libraries to use for shared linking
#         (if the libraries appear static with ".a" extension
#          then Tasmanian_MAGMA_SHARED_FOUND will be OFF)
#
#     Tasmanian_MAGMA_INCLUDE_DIRS: path
#         path to MAGMA headers
#
# recommended use:
#     set(Tasmanian_MAGMA_ROOT_DIR "<path-to-magma>")
#     find_package(TasmanianMagma)
#
# alternatively, use -D Tasmanian_MAGMA_ROOT_DIR=<path-to-magma>
#
#
# detailed use:
#
# MODE 1: User provides a list of libraries Tasmanian_MAGMA_LIBRARIES
#    (required by sXDK cmake policies)
#
#    the module will check if files exist for all libraries given in
#    Tasmanian_MAGMA_LIBRARIES and Tasmanian_MAGMA_SHARED_LIBRARIES
#    if successful, Tasmanian_MAGMA_FOUND is set ON
#    otherwise, Tasmanian_MAGMA_FOUND is set OFF
#
#    Note: this MODE will not create variable Tasmanian_MAGMA_INCLUDE_DIRS
#          in this MODE the user must specify this manually
#
#    if Tasmanian_MAGMA_SHARED_LIBRARIES is not provided
#        the value in Tasmanian_MAGMA_LIBRARIES is used instead
#
#    if no library in Tasmanian_MAGMA_SHARED_LIBRARIES has ".a" extension
#        Tasmanian_MAGMA_SHARED_FOUND is also set ON
#
# MODE 2: User provides a path Tasmanian_MAGMA_ROOT_DIR
#         (assuming Tasmanian_MAGMA_LIBRARIES is not set)
#
#    use find_library() and look for
#        libmagma.{so,dylib,a}
#        libmagma_sparse.{so,dylib,a}
#    using Tasmanian_MAGMA_ROOT_DIR as HINTS
#
#    if Tasmanian_MAGMA_ROOT_DIR is given
#        standard cmake search paths will be omitted from find_library()
#        (xSDK policy: do not silently replace a package if the user has
#         explicitly specified a version)
#
#    if all magma or _sparse libraries are missing
#        set Tasmanian_MAGMA_FOUND to OFF
#
#    if all variants of either magma or magma_sparse end on ".a" suffix,
#        then Tasmanian_MAGMA_SHARED_FOUND is set OFF
#
#    if only shared libraries are found
#        shared libraries will be used in Tasmanian_MAGMA_LIBRARIES
#
#    use find_path() to search for magma.h using HINTS
#        ${Tasmanian_MAGMA_ROOT_DIR}/include
#        ${Tasmanian_MAGMA_ROOT_DIR}
#        <path>/../../include
#        <path>/../include
#        <path>/include
#    where <path> is the folder containing libmagma.so (or .a), usually
#    <magma root>/lib/ or <magma root>/lib/<arch>/, thus we search two
#    levels up from <path>=<magma root>/lib/<arch>/
#    standard cmake search paths are also considered, since if libmagma
#    was found in the Tasmanian_MAGMA_ROOT_DIR, then the headers will be
#    found there first
#
# NOTE: if Tasmanian_MAGMA_ROOT_DIR is missing,
#           find_library() will be called using standard search,
#           if CMAKE_PREFIX_PATH is set correctly, this will find MAGMA
#
# NOTE: the complexity of the logic is needed to ensure
#       - both static and shared libraries are found (if available)
#       - the user can explicitly specify libraries (if desired)
#       - compliance with xSDK policies (if using xSDK defaults)
# https://figshare.com/articles/xSDK_Community_Installation_Policies_GNU_Autoconf_and_CMake_Options/4495133
#
#    Additional confusion can be created by projects that compile static
#    libraries with position independent code, effectively creating
#    shared libraries with .a extension. Thus, the check for the shared
#    library extension cannot be considered as FATAL_ERROR

macro(Tasmanian_magma_checkexists Tasmanian_magma_list)
    # check if the libs in the list exist, if not, set found-off
    foreach(Tasmanian_magma_user_lib ${Tasmanian_magma_list})
        if (NOT EXISTS ${Tasmanian_magma_user_lib})
            message(ERROR " Tasmanian cannon find library ${Tasmanian_magma_user_lib}.")
            set(Tasmanian_MAGMA_FOUND OFF) # missing user provided libraries
            set(Tasmanian_MAGMA_SHARED_FOUND OFF) # assume all is missing
        endif()
    endforeach()
    unset(Tasmanian_magma_user_lib)
endmacro(Tasmanian_magma_checkexists)

macro(Tasmanian_magma_checkestatic Tasmanian_magma_list)
    # checks if any of the libs in the list end on ".a", if so, set shared-found to OFF
    # on systems that don't use ".a" extension, this test will simply pass even if libs are static
    # the build will still fail, but without a meaningful Tasmanian warning (very low concern)
    foreach(Tasmanian_magma_user_lib ${Tasmanian_magma_list})
        get_filename_component(Tasmanian_magma_extension ${Tasmanian_magma_user_lib} EXT)
        if ("${Tasmanian_magma_extension}" STREQUAL ".a") # library appears static
            set(Tasmanian_MAGMA_SHARED_FOUND OFF) # this will only generate a warning
        endif()
    endforeach()
    unset(Tasmanian_magma_user_lib)
    unset(Tasmanian_magma_extension)
endmacro(Tasmanian_magma_checkestatic)

macro(Tasmanian_magma_searchlib Tasmanian_magma_var Tasmanian_magma_libname)
    # usage: Tasmanian_magma_searchlib(Tasmanian_magma_var, Tasmanian_magma_libname)
    # Tasmanian_magma_var -> VAR in the call to find_library(), returns the result
    # Tasmanian_magma_libname -> NAMES in the call to find_lbrary(), e.g., magma or magma_sparse
    # if Tasmanian_MAGMA_ROOT is given, then use as HINTS and skip the cmake and system paths
    if (Tasmanian_MAGMA_ROOT_DIR) # xSDK mode is ON, do no search default paths
        find_library(${Tasmanian_magma_var} ${Tasmanian_magma_libname} HINTS "${Tasmanian_MAGMA_ROOT_DIR}"
                                                                       HINTS "${Tasmanian_MAGMA_ROOT_DIR}/lib/"
                                                                       HINTS "${Tasmanian_MAGMA_ROOT_DIR}/lib/${CMAKE_CXX_LIBRARY_ARCHITECTURE}/"
                                                                       HINTS "${Tasmanian_MAGMA_ROOT_DIR}/${CMAKE_CXX_LIBRARY_ARCHITECTURE}/lib/"
                                                                       NO_DEFAULT_PATH NO_CMAKE_PATH)
    else()
        find_library(${Tasmanian_magma_var} ${Tasmanian_magma_libname})
    endif()
endmacro(Tasmanian_magma_searchlib)

macro(Tasmanian_magma_searchinclude)
    # search for the include path for magma.h using given path as "hint"
    # assumes that Tasmanian_libmagma_dir is the path to libmagma and is set outside of the macro (could be empty)
    find_path(Tasmanian_MAGMA_INCLUDE_DIRS "magma.h" HINTS "${Tasmanian_MAGMA_ROOT_DIR}"
                                                     HINTS "${Tasmanian_MAGMA_ROOT_DIR}/include/"
                                                     HINTS "${Tasmanian_libmagma_dir}"
                                                     HINTS "${Tasmanian_libmagma_dir}/include/"
                                                     HINTS "${Tasmanian_libmagma_dir}/../include/"
                                                     HINTS "${Tasmanian_libmagma_dir}/../../include/")
    if (NOT Tasmanian_MAGMA_INCLUDE_DIRS)
        message(WARNING " Tasmanian cannot find \"magma.h\" header, please provide valid Tasmanian_MAGMA_INCLUDE_DIRS or Tasmanian_MAGMA_ROOT_DIR.")
        set(Tasmanian_MAGMA_FOUND OFF) # cannot work without magma.h header
    endif()
endmacro(Tasmanian_magma_searchinclude)

set(Tasmanian_MAGMA_ROOT_DIR "" CACHE PATH "The root folder for the MAGMA installation, e.g., containing lib and include folders")

set(Tasmanian_MAGMA_FOUND ON) # assume the user has provided everything good, overwrite if check fails in the macros
set(Tasmanian_MAGMA_SHARED_FOUND ON) # overwritten in macro, if all libs appear static

if (NOT Tasmanian_MAGMA_LIBRARIES) # user has not provided specific libs, search for the libs
    # use MAGMA_ROOT_DIR if Tasmanian_MAGMA_ROOT_DIR is not defined
    if ((DEFINED MAGMA_ROOT_DIR) AND (NOT Tasmanian_MAGMA_ROOT_DIR))
        set(Tasmanian_MAGMA_ROOT_DIR ${MAGMA_ROOT_DIR})
    endif()

    Tasmanian_magma_searchlib(Tasmanian_libmagma               magma)
    Tasmanian_magma_searchlib(Tasmanian_libmagma_sparse        magma_sparse)
    Tasmanian_magma_searchlib(Tasmanian_libmagma_static        libmagma.a)
    Tasmanian_magma_searchlib(Tasmanian_libmagma_sparse_static libmagma_sparse.a)

    unset(Tasmanian_MAGMA_LIBRARIES)
    unset(Tasmanian_MAGMA_SHARED_LIBRARIES)

    list(APPEND Tasmanian_MAGMA_SHARED_LIBRARIES ${Tasmanian_libmagma} ${Tasmanian_libmagma_sparse}) # check if exist below

    if (Tasmanian_libmagma_static AND Tasmanian_libmagma_sparse_static)
        list(APPEND Tasmanian_MAGMA_LIBRARIES ${Tasmanian_libmagma_static} ${Tasmanian_libmagma_sparse_static})
        get_filename_component(Tasmanian_libmagma_dir ${Tasmanian_libmagma} DIRECTORY) # used as HINTS to search for "magma.h"
    else()
        set(Tasmanian_MAGMA_LIBRARIES ${Tasmanian_MAGMA_SHARED_LIBRARIES}) # it is OK if shared libs don't exist, will test below
    endif()

    if (NOT Tasmanian_MAGMA_INCLUDE_DIRS) # if not user provided includes, search for "magma.h"
        Tasmanian_magma_searchinclude() # if fail, sets found to OFF
    endif()

    unset(Tasmanian_libmagma_dir)
endif()

# at this point, Tasmanian_MAGMA_LIBRARIES is set, either by the module or the user
# check if the libraries exist and if they appear static

# check if the libraries exist
Tasmanian_magma_checkexists(${Tasmanian_MAGMA_LIBRARIES})

if (Tasmanian_MAGMA_SHARED_LIBRARIES) # check the selected shared libs
    Tasmanian_magma_checkexists(${Tasmanian_MAGMA_SHARED_LIBRARIES})
else()
    set(Tasmanian_MAGMA_SHARED_LIBRARIES ${Tasmanian_MAGMA_LIBRARIES}) # use the given libs
endif()

Tasmanian_magma_checkestatic(${Tasmanian_MAGMA_SHARED_LIBRARIES}) # check if static libs (unix only)
