#!/usr/bin/env bash

set -e

[ -d build ] && rm -rf build
mkdir build && cd build
if [ "${BUILD_TYPE}" == "clang50-python3" ]
then
    cmake \
      -D CMAKE_CXX_FLAGS="-Wall -Wextra" \
      -D CMAKE_CXX_COMPILER=clang++ \
      -D PYTHON_EXECUTABLE=/usr/bin/python3 \
      ..
elif [ "${BUILD_TYPE}" == "gcc73-python2" ]
then
    cmake \
      -D CMAKE_CXX_FLAGS="-Wall -Wextra" \
      -D CMAKE_CXX_COMPILER=g++ \
      -D PYTHON_EXECUTABLE=/usr/bin/python2 \
      ..
else
    echo "Unknown BUILD_TYPE"
    exit 1
fi
make -j${NPROC}
ctest -j${NPROC} --no-compress-output -T Test
