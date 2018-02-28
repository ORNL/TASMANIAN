#!/usr/bin/env bash

set -e

[ -d build ] && rm -rf build
mkdir build && cd build
cmake \
  -D CMAKE_CXX_FLAGS="-Wall -Wextra" \
  -D CMAKE_CXX_COMPILER=g++ \
  -D PYTHON_EXECUTABLE=/usr/bin/python3 \
  -D Tasmanian_ENABLE_MATLAB=OFF \
  ..
make -j${NPROC}
ctest -j${NPROC} --no-compress-output -T Test
