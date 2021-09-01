#!/usr/bin/env bash

if (( @Tasmanian_TESTS_OMP_NUM_THREADS@ != -1 )); then
    export OMP_NUM_THREADS=@Tasmanian_TESTS_OMP_NUM_THREADS@
fi

source "@Tasmanian_final_install_path@"/share/Tasmanian/TasmanianENVsetup.sh || { echo "ERROR: Could not source <install_prefix>/share/Tasmanian/TasmanianENVsetup.sh"; exit 1; }

if [[ $(basename "$(pwd)") != "tasmanian_test_install" ]]; then
    if [ -d tasmanian_test_install ]; then
        cd tasmanian_test_install
    else
        echo "ERROR: must run this script from the tasmanian_test_install folder or the CMake build root folder"
        exit 1
    fi
fi

rm -fr ../tasmanian_test_install/*

@CMAKE_COMMAND@ @Tasmanian_compilers@ "@Tasmanian_final_install_path@/share/Tasmanian/testing"
make -j
make test
