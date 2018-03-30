########################################################################
# C++ 2011 support
########################################################################
if (NOT DEFINED Tasmanian_FORCE_C11)
    set(Tasmanian_FORCE_C11 Tasmanian_PREFER_C11)
endif()

if (Tasmanian_FORCE_C11)
    foreach(tas_target ${tas_target_list})
        set_property(TARGET ${tas_target} PROPERTY CXX_STANDARD 11)
    endforeach()
endif()

########################################################################
# Compiler specific flags: Intel hasn't been tested in a while
########################################################################
if (NOT DEFINED Tasmanian_CXX_FLAGS)
    if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
        #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage") # ./tasgrid -test; gcov CMakeFiles/libtsg_static.dir/SparseGrids/tsgIndexSets.cpp.; geany tsgIndexSets.cpp.gcov
        if (Tasmanian_ENABLE_DEVELOPMENT_DEFAULTS)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wuninitialized -Wstrict-aliasing -pedantic")
            if (TAS_SUPPRESS_OMP_WARN)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
            endif()
        endif()
    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native -diag-disable 11074 -diag-disable 11076 -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic")
    elseif (MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox /EHsc -D_SCL_SECURE_NO_WARNINGS")
    endif()
else()
    set(CMAKE_CXX_FLAGS ${Tasmanian_CXX_FLAGS})
endif()

if (Tasmanian_ENABLE_FORTRAN)
    if (DEFINED Tasmanian_Fortran_FLAGS)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Tasmanian_Fortran_FLAGS}")
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O3 -fno-f2c")
        if (Tasmanian_ENABLE_DEVELOPMENT_DEFAULTS)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wshadow -pedantic")
        endif()
    endif()
endif()

if (DEFINED Tasmanian_EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Tasmanian_EXTRA_CXX_FLAGS}")
endif()

########################################################################
# Extra libraries and directories
########################################################################
if (DEFINED Tasmanian_EXTRA_LIBRARIES)
    foreach(tas_target ${tas_target_list})
        target_link_libraries(${tas_target} ${Tasmanian_EXTRA_LIBRARIES})
    endforeach()
endif()

if (DEFINED Tasmanian_EXTRA_INCLUDE_DIRS)
    include_directories(${Tasmanian_EXTRA_INCLUDE_DIRS})
endif()

if (DEFINED Tasmanian_EXTRA_LINK_DIRS)
    link_directories(${Tasmanian_EXTRA_LINK_DIRS})
endif()
