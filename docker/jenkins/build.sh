#!/usr/bin/env bash

set -e

[ -d build ] && rm -rf build
mkdir build && cd build
if [ "${BUILD_TYPE}" == "clang50-python3" ]
then
    cmake \
      -D CMAKE_INSTALL_PREFIX=./TasmanianInstall
      -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow -pedantic" \
      -D CMAKE_CXX_COMPILER=clang++ \
      -D PYTHON_EXECUTABLE=/usr/bin/python3 \
      -D Tasmanian_TESTS_OMP_NUM_THREADS=4 \
      ..
elif [ "${BUILD_TYPE}" == "gcc73-python2" ]
then
    cmake \
      -D CMAKE_INSTALL_PREFIX=./TasmanianInstall
      -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow -pedantic" \
      -D CMAKE_CXX_COMPILER=g++ \
      -D PYTHON_EXECUTABLE=/usr/bin/python2 \
      -D Tasmanian_TESTS_OMP_NUM_THREADS=4 \
      ..
else
    echo "Unknown BUILD_TYPE"
    exit 1
fi
make -j${NPROC}
ctest -j${NPROC} --no-compress-output -T Test
make install
make test_install

# NPROC is currently set to 4
# combined with Tasmanian_TESTS_OMP_NUM_THREADS=4
# this leads to 4x4 = 16 threads for the ctest
# test_install runs 3 very short tests using the
# system provided OMP_NUM_THREADS
