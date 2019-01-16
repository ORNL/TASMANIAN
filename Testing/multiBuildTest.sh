#!/bin/bash

# analyze current system, cuda, OpenMP and BLAS should be detected automatically
# look for fortran, octave, and various compilers gcc and clang
bGfortran=1
bOctave=1
bGCC=1
bMacOS=0
bPython=1
bPython3=1
bNVCC=1
bMAGMA230=0
/usr/bin/env g++ --version | grep g++
if ((! $? == 0 )); then
    echo "NOTICE: No g++"
    bGCC=0
fi
/usr/bin/env gfortran --version | grep -i fortran
if ((! $? == 0 )); then
    echo "NOTICE: No GNU Fortran"
    bGfortran=0
fi
/usr/bin/env octave --version | grep version
if ((! $? == 0 )); then
    echo "NOTICE: No Octave"
    bOctave=0
fi
/usr/bin/env python --version
if ((! $? == 0 )); then
    echo "NOTICE: No Python"
    bPython=0
fi
/usr/bin/env python3 --version
if ((! $? == 0 )); then
echo "NOTICE: No Python3"
    bPython3=0
fi
/usr/bin/env nvcc --version
if ((! $? == 0 )); then
echo "NOTICE: No NVCC"
    bNVCC=0
fi
if [[ `uname -s` == 'Darwin' ]]; then
    bMacOS=1
fi
if [ -d ~/.magma230/ ]; then
    bMAGMA230=1
    echo "Using MAGMA 2.3.0 in ~/.magma230/"
fi

sPWD=`pwd`
sMultibuildLogFile="$sPWD\multiBuildLog.log"

bShowHelp=1

if (( ${#@} < 1 )); then
    bShowHelp=0
    echo "ERROR: not enough input parameters"
    echo ""
fi

if [[ $1 == *"help" ]]; then
    bShowHelp=0
fi

if (( bShowHelp == 1 && ${#@} < 2 )); then
    bShowHelp=0
    echo "ERROR: not enough input parameters"
    echo ""
fi

if (( bShowHelp == 0 )); then
    echo "Usage: ./multiBuildTest.sh <source folder> <test_dir> <log file>"
    echo "Note: everything in the <test_dir> will be erased"
    echo "      recommend running the script from safe (erasable) location"
    echo ""
    exit 0
fi

echo "--------------------------------------------------------------"
echo "  Performing multiple Build/Tests for Tasmanian"
echo "--------------------------------------------------------------"

echo "Source folder: " $1
echo "Build Folder:  " $2
if [[ ! -z $3 ]]; then
    sMultibuildLogFile=$3
fi
echo "Log file:      " $sMultibuildLogFile
echo "WARNING: everything in $2 will be deleted!"
echo ""

read -p "Press ENTER to begin." sBegin

if [ ! -d "$2" ]; then
    mkdir -p $2
    cd $2
else
    cd $2
    sPWD=`pwd`
    if [[ "$sPWD" == "$2" ]] || [[ "$sPWD/" == "$2" ]]; then
        rm -fr *
    else
        echo "$sPWD"
        echo "$2"
        echo "ERROR: WHERE THE HELL AM I?"
        exit 1
    fi
fi

#cp $1 . || { exit 1; }
set -x

########################################################################
# Here we go!
########################################################################
echo "" > $sMultibuildLogFile
echo "Tasmanian Internal Testing Suite" >> $sMultibuildLogFile
if (( $bMacOS == 1 )); then
    echo "Testing environment: MacOSX" >> $sMultibuildLogFile
else
    echo "Testing environment: Linux" >> $sMultibuildLogFile
fi
if (( $bGCC == 1 )); then
    echo "g++        --- yes" >> $sMultibuildLogFile
else
    echo "g++        ---  no" >> $sMultibuildLogFile
fi
if (( $bGfortran == 1 )); then
    echo "gfortran   --- yes" >> $sMultibuildLogFile
else
    echo "gfortran   ---  no" >> $sMultibuildLogFile
fi
if (( $bOctave == 1 )); then
    echo "octave     --- yes" >> $sMultibuildLogFile
else
    echo "octave     ---  no" >> $sMultibuildLogFile
fi
if (( $bPython == 1 )); then
    echo "python     --- yes" >> $sMultibuildLogFile
else
    echo "python     ---  no" >> $sMultibuildLogFile
fi
if (( $bPython3 == 1 )); then
    echo "python3    --- yes" >> $sMultibuildLogFile
else
    echo "python3    ---  no" >> $sMultibuildLogFile
fi
if (( $bNVCC == 1 )); then
    echo "nvcc       --- yes" >> $sMultibuildLogFile
else
    echo "nvcc       ---  no" >> $sMultibuildLogFile
fi
if (( $bMAGMA230 == 1 )); then
    echo "magma 2.3  --- yes" >> $sMultibuildLogFile
else
    echo "magma 2.3  ---  no" >> $sMultibuildLogFile
fi

sDashFort=""
if (( $bGfortran == 1 )); then
    sDashFort="-fortran"
fi

sTestRoot=$2
cd $sTestRoot
cp -r $1 ./TempTasmanian
cd TempTasmanian
sTempSource=`pwd`
cd ..
mkdir Run
cd Run
sTempBuild=`pwd`
cd $sTestRoot


########################################################################
# Test: simpel GNU Make build system
########################################################################
if (( $bGCC == 1 )); then
    #unzip $1 || { exit 1; }
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    make -j || { echo "Legacy build failed!"; exit 1; }
    make test || { echo "Legacy test failed!"; exit 1; }
    make examples || { echo "Legacy build examples failed!"; exit 1; }
    ./example_sparse_grids -fast || { echo "Legacy SG examples failed!"; exit 1; }
    ./example_dream -fast || { echo "Legacy dream examples failed!"; exit 1; }
    if (( $bGfortran == 1 )) && (( $bMacOS == 0 )); then
        make fortran  || { exit 1; }
        ./example_sparse_grids_f90 -fast || { echo "Legacy fortran examples failed!"; exit 1; }
    else
        echo "Either using OSX or no gfortran, skipping GNU Make fortran test" >> $sMultibuildLogFile
    fi
    if (( $bOctave == 1 )); then
        make matlab  || { exit 1; }
        octave --eval "addpath('$sTempBuild/Tasmanian/InterfaceMATLAB/'); tsgCoreTests()" || { exit 1; }
    else
        echo "No octave executable, skipping GNU Make MATLAB test" >> $sMultibuildLogFile
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "===========================================================================================" >> $sMultibuildLogFile
    echo "======= PASSED: simpel GNU Make build system" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


########################################################################
# Test: default install, install script
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder  -make-j -verbose -nobashrc || { exit 1; }
if [[ ! -z `./TasInstall/bin/tasgrid -v | grep gpu-cuda` ]]; then
    echo "CUDA not disabled by default"
    exit 1;
fi
if (( $bOctave == 1 )); then
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/share/Tasmanian/matlab/'); tsgCoreTests()" || { exit 1; }
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: default install, using install script" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


#########################################################################
## Test: default install, pure cmake
#########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/ || { exit 1; }
mkdir Build || { exit 1; }
cd Build || { exit 1; }
cmake -D CMAKE_BUILD_TYPE=Release -D Tasmanian_ENABLE_RECOMMENDED=ON -D CMAKE_INSTALL_PREFIX="$sTempBuild/Install" $sTempBuild/Tasmanian || { exit 1; }
make -j || { exit 1; }
make test || { exit 1; }
make install || { exit 1; }
make test_install || { exit 1; }
if [ ! -f $sTempBuild/Install/lib/libtasmaniansparsegrid.so ]; then
    echo "cmake did not install correctly"
    exit 1;
fi
cd $sTempBuild
rm -fr Tasmanian/
rm -fr Build/
rm -fr Install/
cd $sTestRoot
echo "======= PASSED: default install, using pure cmake" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


#########################################################################
## Test: pure cmake with lots of warning flags
#########################################################################
if (( $bNVCC == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/ || { exit 1; }
    mkdir Build || { exit 1; }
    cd Build || { exit 1; }
    cmake -D CMAKE_INSTALL_PREFIX="$sTempBuild/Install" \
          -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow" \
          -D BUILD_SHARED_LIBS=ON \
          -D Tasmanian_ENABLE_OPENMP=ON \
          -D Tasmanian_ENABLE_BLAS=ON \
          -D Tasmanian_ENABLE_CUDA=ON \
          -D Tasmanian_ENABLE_PYTHON=ON \
          -D Tasmanian_ENABLE_MPI=OFF \
          -D Tasmanian_ENABLE_FORTRAN=OFF \
          $sTempBuild/Tasmanian || { exit 1; }
    make -j || { exit 1; }
    make test || { exit 1; }
    make install || { exit 1; }
    make test_install || { exit 1; }
    if [ ! -f $sTempBuild/Install/lib/libtasmaniansparsegrid.so ]; then
        echo "cmake did not install correctly"
        exit 1;
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    rm -fr Build/
    rm -fr Install/
    cd $sTestRoot
    echo "======= PASSED: pure cmake with lots of warning flags, all options" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


#########################################################################
## Test: pure cmake with lots of warning flags, no cuda
#########################################################################
if (( $bNVCC == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/ || { exit 1; }
    mkdir Build || { exit 1; }
    cd Build || { exit 1; }
    cmake -D CMAKE_INSTALL_PREFIX="$sTempBuild/Install" \
          -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow" \
          -D BUILD_SHARED_LIBS=ON \
          -D Tasmanian_ENABLE_OPENMP=ON \
          -D Tasmanian_ENABLE_BLAS=ON \
          -D Tasmanian_ENABLE_CUDA=OFF \
          -D Tasmanian_ENABLE_PYTHON=ON \
          -D Tasmanian_ENABLE_MPI=OFF \
          -D Tasmanian_ENABLE_FORTRAN=OFF \
          $sTempBuild/Tasmanian || { exit 1; }
    make -j || { exit 1; }
    make test || { exit 1; }
    make install || { exit 1; }
    make test_install || { exit 1; }
    if [ ! -f $sTempBuild/Install/lib/libtasmaniansparsegrid.so ]; then
        echo "cmake did not install correctly"
        exit 1;
    fi
    if [ -f $sTempBuild/Install/lib/libtasmaniansparsegrid.a ]; then
        echo "cmake did not disable static libs"
        exit 1;
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    rm -fr Build/
    rm -fr Install/
    cd $sTestRoot
    echo "======= PASSED: pure cmake with lots of warning flags, no cuda" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


#########################################################################
## Test: pure cmake build with nvcc and magma 2.3.0
#########################################################################
if (( $bNVCC == 1)) && (( $bMAGMA230 == 1)); then
    OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/.magma230/lib/
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/ || { exit 1; }
    mkdir Build || { exit 1; }
    cd Build || { exit 1; }
    cmake -D CMAKE_INSTALL_PREFIX="$sTempBuild/Install" \
          -D CMAKE_BUILD_TYPE=Release \
          -D BUILD_SHARED_LIBS=ON \
          -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow" \
          -D Tasmanian_ENABLE_OPENMP=ON \
          -D Tasmanian_ENABLE_BLAS=ON \
          -D Tasmanian_ENABLE_CUDA=ON \
          -D Tasmanian_ENABLE_MAGMA=ON \
          -D Tasmanian_MAGMA_ROOT_DIR=~/.magma230 \
          -D Tasmanian_ENABLE_PYTHON=ON \
          -D Tasmanian_ENABLE_MPI=OFF \
          -D Tasmanian_ENABLE_FORTRAN=OFF \
          $sTempBuild/Tasmanian || { exit 1; }
    make -j || { exit 1; }
    make test || { exit 1; }
    make install || { exit 1; }
    make test_install || { exit 1; }
    if [ ! -f $sTempBuild/Install/lib/libtasmaniansparsegrid.so ]; then
        echo "cmake did not install correctly"
        exit 1;
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    rm -fr Build/
    rm -fr Install/
    cd $sTestRoot
    echo "======= PASSED: pure cmake build with nvcc and magma 2.3.0" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
    export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
fi


########################################################################
# Test: -cuda
########################################################################
if (( $bNVCC == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -nostatic -make-j -cuda -verbose -nobashrc || { exit 1; }
    if [[ -z `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$sTempBuild/Tasmanian/TasInstall/lib ./TasInstall/bin/tasgrid -v | grep gpu-cuda` ]]; then
        echo "Failed to enable Nvidia CUDA"
        exit 1;
    fi
    if (( $bOctave == 1 )); then
        LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$sTempBuild/Tasmanian/TasInstall/lib octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/share/Tasmanian/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    if [ -f $sTempBuild/Install/lib/libtasmaniansparsegrid.a ]; then
        echo "cmake did not disable static libs"
        exit 1;
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: with CUDA install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: with CUDA install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


########################################################################
# Test: -python=/usr/bin/python3 -fortran
########################################################################
sUseFortHere=""
if (( $bGfortran == 1 )); then
    sUseFortHere="-fortran"
fi
if (( $bPython3 == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder $sUseFortHere -make-j -python=/usr/bin/python3 -verbose -nobashrc || { exit 1; }
    ./TasInstall/share/Tasmanian/examples/example_sparse_grids.py -fast || { echo "Could not run python3 version of examples"; exit 1; }
    if [[ -z `head -n 1 ./TasInstall/share/Tasmanian/examples/example_sparse_grids.py | grep python3` ]]; then
        echo "Failed to set python3 in the hash-bang of the examples"
        exit 1;
    fi
    if (( $bGfortran == 1 )); then
        if [ ! -f ./TasInstall/lib/libtasmanianfortran90.so ] && [ ! -f ./TasInstall/lib/libtasmanianfortran90.dylib ]; then
            echo "Failed to enable Fortran"
            exit 1;
        fi
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: python3 install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: python3 install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


########################################################################
# Mix different compilers
########################################################################
if [ -f /usr/bin/clang++-5.0 ]; then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    ./install ./TasInstall $sDashFort -cuda -make-j -verbose -nobashrc || { echo "Failed to make a mixed compiler release, gcc to clang"; exit 1; }
    mkdir -p TasExamples || { echo "Failed to build mixed compiler examples"; exit 1; }
    cd TasExamples
    cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++-5.0 ../TasInstall/share/Tasmanian/examples/ || { echo "Failed to cmake mixed compiler examples"; exit 1; }
    make -j || { echo "Failed to make mixed compiler examples"; exit 1; }
    ./example_sparse_grids -fast || { echo "Failed to run mixed compiler examples, sparse grid"; exit 1; }
    ./example_dream -fast || { echo "Failed to run mixed compiler examples, dream"; exit 1; }
    if [ -f example_sparse_grids_f90 ]; then
        ./example_sparse_grids_f90 -fast || { echo "Failed to run mixed compiler examples, fortran"; exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot

    # something is broken here, g++ cannot link static libraries build with clang, shared works but static is somehow broken
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    ./install ./TasInstall $sDashFort -cuda -make-j -nostatic -verbose -nobashrc -compiler=/usr/bin/clang++-5.0 || { echo "Failed to make a mixed compiler release, clang to gcc"; exit 1; }
    mkdir -p TasExamples || { echo "Failed to build mixed compiler examples"; exit 1; }
    cd TasExamples
    source ../TasInstall/share/Tasmanian/TasmanianENVsetup.sh
    cmake -DCMAKE_CXX_COMPILER=/usr/bin/g++ ../TasInstall/share/Tasmanian/examples/ || { echo "Failed to cmake mixed compiler examples"; exit 1; }
    make -j || { echo "Failed to make mixed compiler examples"; exit 1; }
    ./example_sparse_grids -fast || { echo "Failed to run mixed compiler examples, sparse grid"; exit 1; }
    ./example_dream -fast || { echo "Failed to run mixed compiler examples, dream"; exit 1; }
    if [ -f example_sparse_grids_f90 ]; then
       ./example_sparse_grids_f90 -fast || { echo "Failed to run mixed compiler examples, fortran"; exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
fi
echo "======= PASSED: mixed compilers" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Clean up
########################################################################
cd $sTestRoot
sPWD=`pwd`
if [[ "$sPWD" == "$sTestRoot" ]] || [[ "$sPWD/" == "$sTestRoot" ]]; then
    rm -fr TempTasmanian
    rm -fr Run
else
    echo "$sPWD"
    echo "$sTestRoot"
    echo "ERROR: WHERE THE HELL AM I?"
    exit 1
fi


########################################################################
# Final message
########################################################################
set +x

echo "===========================================================================================" >> $sMultibuildLogFile
echo "======= PASSED: ALL TESTED" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile

echo ""
echo ""
echo "------------------------------------------------------------------------"
echo "            MULTIBUILD TEST ALL TESTS COMPLETED SUCCESSFULLY"
echo "------------------------------------------------------------------------"

cat $sMultibuildLogFile

