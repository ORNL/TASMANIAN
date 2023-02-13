#!/bin/bash

pwd
sed -i -e 's|{CMAKE_SOURCE_DIR|{CMAKE_CURRENT_SOURCE_DIR|g' CMakeLists.txt
sed -i -e 's|{CMAKE_BINARY_DIR|{CMAKE_CURRENT_BINARY_DIR|g' CMakeLists.txt
sed -i -e 's|set(CUDA_ARCHITECTURES|set(tsgdummy|g' CMakeLists.txt
