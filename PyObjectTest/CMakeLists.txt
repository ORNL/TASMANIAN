cmake_minimum_required(VERSION 3.10)

cmake_policy(VERSION 3.10)
project(CppLib VERSION 0.0.1 LANGUAGES CXX)

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Threads REQUIRED)


add_library(cpplib MODULE cpplib.cpp)
target_include_directories(cpplib PUBLIC $<BUILD_INTERFACE:${PYTHON_INCLUDE_DIRS}>)
target_link_libraries(cpplib ${PYTHON_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})


configure_file("${CMAKE_CURRENT_SOURCE_DIR}/run.sh" "${CMAKE_CURRENT_BINARY_DIR}/run.sh")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/main.py" "${CMAKE_CURRENT_BINARY_DIR}/main.py")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pytestmod.py" "${CMAKE_CURRENT_BINARY_DIR}/pytestmod.py")

message(STATUS "NOTE: run the example with: ${PYTHON_EXECUTABLE} main.py")
