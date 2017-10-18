#!/bin/bash

sPWD=`pwd`
if [ $sPWD != @CMAKE_BINARY_DIR@ ]; then
    echo "NOPE: you must run this inside @CMAKE_BINARY_DIR@"
    exit 1;
fi

if [ ! -d "FinalTest" ]; then
    mkdir FinalTest
fi
cd FinalTest || { echo "ERROR: Could not cd into FinalTest, aborting"; exit 1; }

sPWD=`pwd`
if [ $sPWD != "@CMAKE_BINARY_DIR@/FinalTest" ]; then
    echo "ERROR: somehow we failed to cd into FinalTest"
    exit 1;
fi

rm -fr *

echo "--------------------------------------------------------------------------------"
echo " Test 1: source the PATH and LD_LIBRARY_PATH and run the executable"
echo "--------------------------------------------------------------------------------"

source @CMAKE_INSTALL_PREFIX@/config/TasmanianENVsetup.sh || { echo "ERROR: Could not source <install_prefix>/config/TasmanianENVsetup.sh"; exit 1; }
tasgrid -v  || { echo "ERROR: Could not execute ./tasgrid -v"; exit 1; }


echo "--------------------------------------------------------------------------------"
echo " Test 2: compile and run the C++ examples"
echo "--------------------------------------------------------------------------------"
echo 'Building  "cmake @CMAKE_INSTALL_PREFIX@/examples"'
cmake @CMAKE_INSTALL_PREFIX@/examples > /dev/null || { echo "ERROR: Could not cmake the C++ examples"; exit 1; }
echo 'Compiling "make"'
make > /dev/null || { echo "ERROR: Could not compile the C++ examples"; exit 1; }
echo 'Executing "./example_sparse_grids"'
./example_sparse_grids -fast >/dev/null || { echo "ERROR: Could not run the C++ Sparse Grid example"; exit 1; }
echo 'Executing "./example_sparse_grids_fortran"'
./example_sparse_grids_fortran -fast >/dev/null || { echo "ERROR: Could not run the Fortran Sparse Grid example"; exit 1; }
echo 'Executing "./example_dream"'
./example_dream -fast >/dev/null || { echo "ERROR: Could not run the C++ DREAM example"; exit 1; }


echo "--------------------------------------------------------------------------------"
echo " Test 3: run a basic python test"
echo "--------------------------------------------------------------------------------"
sPSuccess=1
if [ "@TSGPythonPath@" == "@CMAKE_INSTALL_PREFIX@/python" ]; then
    echo 'Executing "@CMAKE_INSTALL_PREFIX@/examples/example_sparse_grids.py"'
    @CMAKE_INSTALL_PREFIX@/examples/example_sparse_grids.py -fast > /dev/null || { echo "ERROR: could not run the python example post install!"; sPSuccess=0; }
    head -n 1 @CMAKE_INSTALL_PREFIX@/python/TasmanianSG.py > dummy.py
    echo 'import sys' >> dummy.py
    echo 'sys.path.append("@CMAKE_INSTALL_PREFIX@/python")' >> dummy.py
    echo 'import TasmanianSG' >> dummy.py
    echo 'print("TasmanianSG Python module version: {0:1s}".format(TasmanianSG.__version__))' >> dummy.py
    chmod 755 dummy.py
    ./dummy.py || { echo "ERROR: Could not run the dummy python test"; echo "      This is a problem either with Python install or the Tasmanian library."; sPSuccess=0; }
else
    echo "Python not enabled, skipping"
fi


echo "--------------------------------------------------------------------------------"
echo " Test 4: source the CPLUS_INCLUDE_PATH and LIBRARY_PATH and compile a dummy test"
echo "--------------------------------------------------------------------------------"

sSuccess=1
if [[ ! -z `which g++` ]]; then

    source @CMAKE_INSTALL_PREFIX@/config/TasmanianDEVsetup.sh || { echo "ERROR: Could not source <install_prefix>/config/TasmanianDEVsetup.sh"; exit 1; }

    echo '#include <iostream>' > dummy.cpp
    echo '#include "TasmanianSparseGrid.hpp"' >> dummy.cpp
    echo 'using namespace std;' >> dummy.cpp
    echo 'int main(int argc, const char ** argv){' >> dummy.cpp
    echo 'cout << "Tasmanian Sparse Grids  version: " << TasGrid::TasmanianSparseGrid::getVersion() << endl;' >> dummy.cpp
    echo 'return 0;' >> dummy.cpp
    echo '}' >> dummy.cpp

    # check is running on mac
    sUname="$(uname -s)"
    
    if [ $sUname == Linux ]; then
        g++ -fopenmp dummy.cpp -o dummy_test -ltasmaniansparsegrid || { sSuccess=0; }
        ./dummy_test || { sSuccess=0; }
    else
        g++ dummy.cpp -o dummy_test -ltasmaniansparsegrid || { sSuccess=0; }
    fi
    if (( $sSuccess == 0 )); then
        echo "ERROR: could not compile simple g++ test, but cmake exmaples worked."
        echo "       This is probably an issue with the compiler or the simple compile command."
        echo "       If you use cmake exported targets to link to the libraries, you don't need this anyway."
    fi
else
    echo "Could not find g++ command, cannot test this. If you are using cmake exported targets, you don't need this anyway."
fi

if (( $sPSuccess == 0 )) || (( $sSuccess == 0 )); then
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "   SOME TESTS FAILED, but the install is probably OK"
    echo "--------------------------------------------------------------------------------"
    echo ""
    exit 1;
else
    echo ""
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "   ALL POST INSTALL TESTS COMPLETED SUCCESSFULLY"
    echo "--------------------------------------------------------------------------------"
    echo ""
fi
exit 0;
