#!/usr/bin/env bash

set -e

[ -d build ] && rm -rf build
mkdir build && cd build
if [ "${BUILD_TYPE}" == "clang50-python3" ]
then
# Enables only core libraries and PYTHON
    cmake \
      -D CMAKE_INSTALL_PREFIX=./TasmanianInstall \
      -D CMAKE_CXX_FLAGS="-O3 -Wall -Wextra -Wshadow -pedantic" \
      -D CMAKE_CXX_COMPILER=clang++ \
      -D Tasmanian_ENABLE_PYTHON=ON \
      -D PYTHON_EXECUTABLE=/usr/bin/python3 \
      -D Tasmanian_TESTS_OMP_NUM_THREADS=4 \
      ..
elif [ "${BUILD_TYPE}" == "gcc73-python2" ]
then
# Enables core libraries, PYTHON and OpenMP
# Attempt to enable BLAS, BLAS is missing but this should still work
    cmake \
      -D CMAKE_INSTALL_PREFIX=./TasmanianInstall \
      -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow -pedantic" \
      -D CMAKE_CXX_COMPILER=g++ \
      -D Tasmanian_ENABLE_RECOMMENDED=ON \
      -D Tasmanian_ENABLE_FORTRAN=ON \
      -D PYTHON_EXECUTABLE=/usr/bin/python2 \
      -D Tasmanian_TESTS_OMP_NUM_THREADS=4 \
      ..
elif [ "${BUILD_TYPE}" == "gcc54-cuda90" ]
then
# Explicitly enable and test all options considered "stable" and CI worthy
# MATLAB is "stable" but tested separately
# MPI is "experimental", FORTRAN is "very experimental"
#
# NOTE: CUDA and -pedantic results in many warning messages
#       as far as I can tell, this has nothing to do with Tasmanian
    cmake \
      -D CMAKE_INSTALL_PREFIX=./TasmanianInstall \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow" \
      -D Tasmanian_ENABLE_OPENMP=ON \
      -D Tasmanian_ENABLE_BLAS=ON \
      -D Tasmanian_ENABLE_CUDA=ON \
      -D Tasmanian_ENABLE_MAGMA=OFF \
      -D Tasmanian_ENABLE_PYTHON=ON \
      -D Tasmanian_ENABLE_MPI=OFF \
      -D Tasmanian_ENABLE_FORTRAN=ON \
      -D Tasmanian_TESTS_OMP_NUM_THREADS=4 \
      ..
else
    echo "Unknown BUILD_TYPE"
    exit 1
fi
make -j${NPROC}
ctest -j1 --no-compress-output -V -T Test
make install
make test_install

# NPROC is currently set to 4
# combined with Tasmanian_TESTS_OMP_NUM_THREADS=4
# this leads to 4x4 = 16 threads for the ctest
# test_install runs 3 very short tests using the
# system provided OMP_NUM_THREADS
