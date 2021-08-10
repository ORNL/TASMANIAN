########################################################################
# Add tpltest
########################################################################
add_executable(Tasmanian_tpltest tpltest_main.cpp
                                 tplUnitTests.hpp
                                 tplUnitTests.cpp)

set_target_properties(Tasmanian_tpltest  PROPERTIES OUTPUT_NAME "tpltest"
                                                    CXX_EXTENSIONS OFF)
target_link_libraries(Tasmanian_tpltest  Tasmanian_dependencies)
                                                    
########################################################################
# Testing
########################################################################
add_test(LAPACK tpltest all)

