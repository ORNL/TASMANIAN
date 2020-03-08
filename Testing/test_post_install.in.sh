#!/usr/bin/env bash

#TasmanianPostInstallTest

if (( @Tasmanian_TESTS_OMP_NUM_THREADS@ != -1 )); then
    export OMP_NUM_THREADS=@Tasmanian_TESTS_OMP_NUM_THREADS@
fi

if [ ! -d "TasmanianPostInstallTest" ]; then
    mkdir TasmanianPostInstallTest
fi
cd TasmanianPostInstallTest || { echo "ERROR: Could not cd into TasmanianPostInstallTest, terminating"; exit 1; }

rm -fr ../TasmanianPostInstallTest/*

echo ""
echo "--------------------------------------------------------------------------------"
echo " Test 1: source the PATH and the LD_LIBRARY_PATH and run the executable"
echo "--------------------------------------------------------------------------------"

source "@Tasmanian_final_install_path@"/share/Tasmanian/TasmanianENVsetup.sh || { echo "ERROR: Could not source <install_prefix>/share/Tasmanian/TasmanianENVsetup.sh"; exit 1; }
tasgrid -v  || { echo "ERROR: Could not execute ./tasgrid -v"; exit 1; }


echo "--------------------------------------------------------------------------------"
echo " Test 2: compile and run the C++ examples"
echo "--------------------------------------------------------------------------------"
echo 'Building  "cmake @Tasmanian_final_install_path@/share/Tasmanian/examples"'
if [[ -z $1 ]]; then
    @CMAKE_COMMAND@ @Tasmanian_compilers@ "@Tasmanian_final_install_path@/share/Tasmanian/examples" || { echo "ERROR: Could not cmake the C++ examples"; exit 1; }
else
    @CMAKE_COMMAND@ $1 "@Tasmanian_final_install_path@/share/Tasmanian/examples" || { echo "ERROR: Could not cmake the C++ examples"; exit 1; }
fi
echo 'Compiling "make"'
make -j3 || { echo "ERROR: Could not compile the C++ examples"; exit 1; }
echo 'Executing "./example_sparse_grids"'
./example_sparse_grids -fast >/dev/null || { echo "ERROR: Could not run the C++ Sparse Grid example"; exit 1; }
if [ -f "@Tasmanian_final_install_path@"/share/Tasmanian/examples/example_sparse_grids.f90 ]; then
    echo 'Executing "./example_sparse_grids_fortran"'
    ./example_sparse_grids_fortran -fast >/dev/null 2>&1 || { echo "ERROR: Could not run the Fortran Sparse Grid example"; exit 1; }
fi
echo 'Executing "./example_dream"'
./example_dream -fast >/dev/null || { echo "ERROR: Could not run the C++ DREAM example"; exit 1; }


echo ""
echo "--------------------------------------------------------------------------------"
echo " Test 3: run a basic python test"
echo "--------------------------------------------------------------------------------"
sPSuccess=1
if [[ "@Tasmanian_ENABLE_PYTHON@" == "ON" ]]; then
    echo 'Executing "@Tasmanian_final_install_path@/share/examples/example_sparse_grids.py"'
    "@PYTHON_EXECUTABLE@" "@Tasmanian_final_install_path@"/share/Tasmanian/examples/example_sparse_grids.py -fast > /dev/null || { echo "ERROR: could not run the python example post install!"; sPSuccess=0; }
    echo 'import Tasmanian' >> hello_world.py
    echo 'print("Tasmanian Python module version: {0:1s}".format(Tasmanian.__version__))' >> hello_world.py
    "@PYTHON_EXECUTABLE@" hello_world.py || { echo "ERROR: Could not run the dummy python test"; echo "      This is a problem either with Python install or the Tasmanian library."; sPSuccess=0; }
else
    echo "Python not enabled, skipping"
fi


if (( $sPSuccess == 0 )); then
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "   SOME TESTS FAILED, but the install may be OK"
    echo "--------------------------------------------------------------------------------"
    echo ""
    exit 1;
else
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "   ALL POST INSTALL TESTS COMPLETED SUCCESSFULLY"
    echo "--------------------------------------------------------------------------------"
    echo ""
fi

exit 0;
