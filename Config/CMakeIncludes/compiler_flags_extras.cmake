########################################################################
# Compiler specific flags: Intel hasn't been tested in a while
########################################################################
if (NOT Tasmanian_STRICT_OPTIONS)
    if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
        if (Tasmanian_ENABLE_FORTRAN)
            set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O3 -fno-f2c")
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
