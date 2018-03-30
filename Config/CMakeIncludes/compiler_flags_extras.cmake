########################################################################
# C++ 2011 support
########################################################################
if (DEFINED Tasmanian_ENABLE_CXX_2011)
    if (Tasmanian_ENABLE_CXX_2011)
        foreach(Tasmanian_loop_target ${Tasmanian_target_list})
            set_property(TARGET ${Tasmanian_loop_target} PROPERTY CXX_STANDARD 11)
        endforeach()
    endif()
else()
    if (Tasmanian_ENABLE_CUDA OR Tasmanian_ENABLE_CUBLAS OR Tasmanian_ENABLE_MPI)
        foreach(Tasmanian_loop_target ${Tasmanian_target_list})
            set_property(TARGET ${Tasmanian_loop_target} PROPERTY CXX_STANDARD 11)
        endforeach()
    endif()
endif()

########################################################################
# Compiler specific flags: Intel hasn't been tested in a while
########################################################################
if (NOT Tasmanian_STRICT_OPTIONS)
    if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
        if (Tasmanian_ENABLE_FORTRAN)
            set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O3 -fno-f2c")
        endif()
        if (Tasmanian_ENABLE_DEVELOPMENT_DEFAULTS)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wuninitialized -Wstrict-aliasing -pedantic")
            if (NOT Tasmanian_ENABLE_OPENMP)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
            endif()
            if (Tasmanian_ENABLE_FORTRAN)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wshadow -pedantic")
            endif()
        endif()
#    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native -diag-disable 11074 -diag-disable 11076 -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic")
    elseif (MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox /EHsc -D_SCL_SECURE_NO_WARNINGS")
    endif()
endif()


########################################################################
# Extra flags, libraries, and directories
########################################################################
if (DEFINED Tasmanian_EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Tasmanian_EXTRA_CXX_FLAGS}")
endif()

if (DEFINED Tasmanian_EXTRA_LIBRARIES)
    foreach(Tasmanian_loop_target ${Tasmanian_target_list})
        target_link_libraries(${Tasmanian_loop_target} ${Tasmanian_EXTRA_LIBRARIES})
    endforeach()
endif()

if (DEFINED Tasmanian_EXTRA_INCLUDE_DIRS)
    include_directories(${Tasmanian_EXTRA_INCLUDE_DIRS})
endif()

if (DEFINED Tasmanian_EXTRA_LINK_DIRS)
    link_directories(${Tasmanian_EXTRA_LINK_DIRS})
endif()
