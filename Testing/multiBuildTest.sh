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
echo "Log file:  " $sMultibuildLogFile
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

echo $sTempSource
echo $sTempBuild

########################################################################
# Test: simpel GNU Make build system
########################################################################
if (( $bGCC == 1 )); then
    #unzip $1 || { exit 1; }
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    make -j || { echo "Legacy build failed!"; exit 1; }
    ./tasgrid -test || { echo "Legacy tasgrid failed!"; exit 1; }
    ./tasdream -test || { echo "Legacy tasdream failed!"; exit 1; }
    ./testTSG.py || { echo "Legacy python test failed!"; exit 1; }
    make examples || { echo "Legacy build examples failed!"; exit 1; }
    ./example_sparse_grids -fast || { echo "Legacy SG examples failed!"; exit 1; }
    ./example_dream -fast || { echo "Legacy dream examples failed!"; exit 1; }
    if (( $bGfortran == 1 )) && (( $bMacOS == 0 )); then
        make fortran
        ./example_sparse_grids_fortran -fast || { echo "Legacy fortran examples failed!"; exit 1; }
    else
        echo "Either using OSX or no gfortran, skipping GNU Make fortran test" >> $sMultibuildLogFile
    fi
    if (( $bOctave == 1 )); then
        make matlab
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
# Test: default install, install.sh
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
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: default install, using install.sh" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


#########################################################################
## Test: default install, pure cmake
#########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/ || { exit 1; }
mkdir Build || { exit 1; }
cd Build || { exit 1; }
cmake -D CMAKE_INSTALL_PREFIX="$sTempBuild/Install" $sTempBuild/Tasmanian || { exit 1; }
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
cd $sTestRoot
echo "======= PASSED: default install, using pure cmake" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


#########################################################################
## Test: pure cmake with lots of warning flags, strict options
#########################################################################
if (( $bNVCC == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/ || { exit 1; }
    mkdir Build || { exit 1; }
    cd Build || { exit 1; }
    cmake -D CMAKE_INSTALL_PREFIX="$sTempBuild/Install" \
          -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow" \
          -D Tasmanian_STRICT_OPTIONS=ON \
          -D Tasmanian_ENABLE_OPENMP=ON \
          -D Tasmanian_ENABLE_BLAS=ON \
          -D Tasmanian_ENABLE_CUBLAS=ON \
          -D Tasmanian_ENABLE_CUDA=ON \
          -D Tasmanian_ENABLE_PYTHON=ON \
          -D Tasmanian_SHARED_LIBRARY=ON \
          -D Tasmanian_STATIC_LIBRARY=ON \
          -D Tasmanian_ENABLE_MATLAB=OFF \
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
    cd $sTestRoot
    echo "======= PASSED: pure cmake with lots of warning flags, all options" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


#########################################################################
## Test: pure cmake with lots of warning flags, cuda + no cublas
#########################################################################
if (( $bNVCC == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/ || { exit 1; }
    mkdir Build || { exit 1; }
    cd Build || { exit 1; }
    cmake -D CMAKE_INSTALL_PREFIX="$sTempBuild/Install" \
          -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow" \
          -D Tasmanian_STRICT_OPTIONS=ON \
          -D Tasmanian_ENABLE_OPENMP=ON \
          -D Tasmanian_ENABLE_BLAS=ON \
          -D Tasmanian_ENABLE_CUBLAS=OFF \
          -D Tasmanian_ENABLE_CUDA=ON \
          -D Tasmanian_ENABLE_PYTHON=ON \
          -D Tasmanian_SHARED_LIBRARY=ON \
          -D Tasmanian_STATIC_LIBRARY=ON \
          -D Tasmanian_ENABLE_MATLAB=OFF \
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
    cd $sTestRoot
    echo "======= PASSED: pure cmake with lots of warning flags, cuda + no cublas" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


########################################################################
# Test: -noomp
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder  -make-j -noomp -verbose -nobashrc || { exit 1; }
if [[ ! -z `./TasInstall/bin/tasgrid -v | grep 'multithreading: Enabled'` ]]; then
    echo "Failed to disable OpenMP"
    exit 1;
fi
if (( $bOctave == 1 )); then
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: no OpenMP install" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: -noblas
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -make-j -noblas -verbose -nobashrc || { exit 1; }
if [[ ! -z `./TasInstall/bin/tasgrid -v | grep cpu-blas` ]]; then
    echo "Failed to disable BLAS"
    exit 1;
fi
if (( $bOctave == 1 )); then
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: no BLAS install" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: -nocublas
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -make-j -nocublas -verbose -nobashrc || { exit 1; }
if [[ ! -z `./TasInstall/bin/tasgrid -v | grep gpu-cublas` ]]; then
    echo "Failed to disable Nvidia cuBLAS"
    exit 1;
fi
if (( $bOctave == 1 )); then
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: no cuBLAS install" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: -cuda
########################################################################
if (( $bNVCC == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -make-j -cuda -verbose -nobashrc || { exit 1; }
    if [[ -z `./TasInstall/bin/tasgrid -v | grep gpu-cuda` ]]; then
        echo "Failed to enable Nvidia CUDA"
        exit 1;
    fi
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
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
# Test: -cuda -nocublas
########################################################################
if (( $bNVCC == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -make-j -cuda -nocublas -verbose -nobashrc || { exit 1; }
    if [[ -z `./TasInstall/bin/tasgrid -v | grep gpu-cuda` ]]; then
        echo "Failed to enable Nvidia CUDA without CUBLAS"
        exit 1;
    fi
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: with CUDA and no CUBLAS install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: with CUDA and no CUBLAS install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


########################################################################
# Test: -nopython
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -make-j -nopython -verbose -nobashrc || { exit 1; }
if [ -f ./TasInstall/python/TasmanianSG.py ]]; then
    echo "Failed to disable Python"
    exit 1;
fi
if (( $bOctave == 1 )); then
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: no python install" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile

########################################################################
# Test: -nospam
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -make-j -nospam -verbose -nobashrc || { exit 1; }
if [[ -f ./TasInstall/python/TasmanianSG.py ]]; then
    echo "Failed to disable Spam (python)"
    exit 1;
fi
if (( $bOctave == 1 )); then
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: no spam (python) install" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: -python=/usr/bin/python3
########################################################################
if (( $bPython3 == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -make-j -python=/usr/bin/python3 -verbose -nobashrc || { exit 1; }
    ./TasInstall/examples/example_sparse_grids.py -fast || { echo "Could not run python3 version of examples"; exit 1; }
    if [[ -z `head -n 1 ./TasInstall/examples/example_sparse_grids.py | grep python3` ]]; then
        echo "Failed to set python3 in the hash-bang of the examples"
        exit 1;
    fi
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
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
# Test: -fortran
########################################################################
if (( $bGfortran == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -make-j -fortran -verbose -nobashrc || { exit 1; }
    if [ ! -f ./TasInstall/lib/libtasmanianfortran.so ] && [ ! -f ./TasInstall/lib/libtasmanianfortran.dylib ]; then
        echo "Failed to enable Fortran"
        exit 1;
    fi
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: fortran install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: fortran install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


########################################################################
# Test: -fortran -noshared -nopython
########################################################################
if (( $bGfortran == 1 )); then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -make-j -fortran -noshared -nospam -verbose -nobashrc || { exit 1; }
    if [ ! -f ./TasInstall/lib/libtasmanianfortran.a ]; then
        echo "Failed to enable Fortran"
        exit 1;
    fi
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: fortran static install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: fortran static install" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


########################################################################
# Test: default install, reject -noshared and -nostatic
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -make-j -noshared -nostatic -verbose -nobashrc
if (( $? == 0 )); then
    echo "Failed to reject install with both -noshared and -nostatic"
    exit 1;
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: reject simultaneous -noshared and -nostatic" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: default install, -nostatic and environment setup
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -nostatic -make-j -verbose -nobashrc || { exit 1; }
# MacOSX always finds the dynamic libs as they are hardcoded with absolute path
if (( $bMacOS == 0 )); then
    $sTempBuild/Tasmanian/TasInstall/bin/tasgrid -v
    if [ $? -eq 0 ]; then
        echo "Tasgrid is supposed to give an error here!"
        exit 1;
    else
        echo "Tasgrid is supposed to fail above"
    fi
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()"
        if [ $? -eq 0 ]; then
            echo "Octave is supposed to give an error here!"
            exit 1;
        else
            echo "Octave is supposed to fail above."
        fi
    fi
    if [ -f $sTempBuild/Tasmanian/TasInstall/libtasmaniansparsegrid.a ]; then
        echo "Failed to apply -nostatic"
        exit 1;
    fi
    if [ ! -f $sTempBuild/Tasmanian/TasInstall/lib/libtasmaniansparsegrid.so ]; then
        echo "Tasmanian shared library is missing"
        exit 1;
    fi
else
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { echo "Failed Octave test build with shared lib (MacOSX)!"; exit 1; }
    fi
    $sTempBuild/Tasmanian/TasInstall/bin/tasgrid -v || { echo "Failed Tasgrid build with shared lib (MacOSX)!"; exit 1; }
fi
# save the environment
SAVE_PATH=$PATH
SAVE_LD_LIB_PATH=$LD_LIBRARY_PATH
# source the new environment (any OS above should work fine)
source $sTempBuild/Tasmanian/TasInstall/config/TasmanianENVsetup.sh
tasgrid -v || { echo "Failed Tasgrid test build with shared lib (loaded setup.sh)!"; exit 1; }
if (( $bOctave == 1 )); then
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { echo "Failed Octave test build with shared lib (loaded setup.sh)!"; exit 1; }
fi
# restore the environment
export PATH=$SAVE_PATH
export LD_LIBRARY_PATH=$SAVE_LD_LIB_PATH
# should be back to the case of MacOSX working and Linux failing to find the shared libs
if (( $bMacOS == 0 )); then
    $sTempBuild/Tasmanian/TasInstall/bin/tasgrid -v
    if [ $? -eq 0 ]; then
        echo "Tasgrid is supposed to give an error here!"
        exit 1;
    else
        echo "Tasgrid is supposed to fail above."
    fi
else
    #echo "MacOSX should find the dynamic libraries regardless of LD_LIBRARY_PATH"
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { echo "Failed Octave test build with shared lib (MacOSX v2)!"; exit 1; }
    $sTempBuild/Tasmanian/TasInstall/bin/tasgrid -v || { echo "Failed Octave test build with shared lib (MacOSX v2)!"; exit 1; }
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: -nostatic" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: default install, -noshared, which is overwritten due to python
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -noshared -make-j -verbose -nobashrc || { exit 1; }
if (( $bOctave == 1 )); then
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
fi
if (( $bMacOS == 0 )); then
    if [ ! -f $sTempBuild/Tasmanian/TasInstall/lib/libtasmaniansparsegrid.so ]; then
        echo "Tasmanian shared library is missing even though enabling python should have overwritten -noshared"
        exit 1;
    fi
else
    if [ ! -f $sTempBuild/Tasmanian/TasInstall/lib/libtasmaniansparsegrid.dylib ]; then
        echo "Tasmanian shared library is missing even though enabling python should have overwritten -noshared"
        exit 1;
    fi
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: -noshared" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: -noshared and -nopython
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -noshared -nopython -make-j -verbose -nobashrc || { exit 1; }
if (( $bOctave == 1 )); then
    octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
fi
if [ -f $sTempBuild/Tasmanian/TasInstall/lib/libtasmaniansparsegrid.so ]; then
    echo "Tasmanian shared library is present even though it should be disabled"
    exit 1;
fi
if [ -d $sTempBuild/Tasmanian/TasInstall/python ]; then
    echo "Failed to disable python and shared libraries"
    exit 1;
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: -noshared -nopython" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: default install, no MATLAB
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
./install ./TasInstall -make-j -verbose -nobashrc || { exit 1; }
octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()"
if [ $? -eq 0 ]; then
    echo "Ocatve is not supposed to work here"
    exit 1;
fi
if [ -d ./TasInstall/matlab ]; then
    echo "Failed to disable MATLAB"
    exit 1;
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: no MATLAB" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: -noinstall
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -noinstall -make-j -verbose -nobashrc || { exit 1; }
if [ -d $sTempBuild/Tasmanian/TasInstall/lib ] || [ -d $sTempBuild/Tasmanian/TasInstall/bin ]; then
    echo "Failed to disable the install command."
    exit 1;
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: -noinstall" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: -notest
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
# break the test, so if it executes it will break the install
echo "(" >> ./Testing/testTSG.in.py
sed -i -e 's/\#TasmanianPostInstallTest/exit\ 1/g' ./Testing/test_post_install.in.sh
./install ./TasInstall ./tsgWorkFolder -notest -make-j -verbose -nobashrc || { echo "Failed to disable testing"; exit 1; }
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: -notest" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: -debug
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder $sDashFort -debug -cuda -make-j -verbose -nobashrc || { echo "Failed to make a debug release"; exit 1; }
if [[ -z `cmake -LA -N Build/ | grep CMAKE_BUILD_TYPE:STRING=Debug` ]]; then
    echo "Failed to set a debug build"
    exit 1;
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: -debug" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Test: -wrong
########################################################################
cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
cd $sTempBuild/Tasmanian || { exit 1; }
mkdir -p tsgWorkFolder
./install ./TasInstall ./tsgWorkFolder -cuda -make-j -verbose -nobashrc -wrong
if [ $? -eq 0 ]; then
    echo "Failed to exit on wrong install input"
    exit 1;
fi
cd $sTempBuild
rm -fr Tasmanian/
cd $sTestRoot
echo "======= PASSED: -wrong" >> $sMultibuildLogFile
echo "===========================================================================================" >> $sMultibuildLogFile


########################################################################
# Alternative compiler tests: cuda versions, clang, etc.
########################################################################
if [ -f /usr/bin/clang++-5.0 ]; then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    if (( $bMacOS == 1 )); then
        ./install ./TasInstall ./tsgWorkFolder -cmake="-DCMAKE_CXX_COMPILER=/usr/bin/clang++-5.0" $sDashFort -noomp -cuda -make-j -verbose -nobashrc || { exit 1; }
    else
        ./install ./TasInstall ./tsgWorkFolder -cmake="-DCMAKE_CXX_COMPILER=/usr/bin/clang++-5.0" $sDashFort -cuda -make-j -verbose -nobashrc || { exit 1; }
    fi
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    if (( $bMacOS == 0 )); then
        if [ -z `./TasInstall/bin/tasgrid -v | grep 'OpenMP multithreading: Enabled'` ]; then
            echo "OpenMP is supposed to work with clang++-5.0, but it failed somehow"
            exit 1;
        fi
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: Clang 5.0" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: Clang 5.0" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi

if [ -f /usr/bin/clang++-4.0 ]; then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    if (( $bMacOS == 1 )); then
        ./install ./TasInstall ./tsgWorkFolder -cmake="-DCMAKE_CXX_COMPILER=/usr/bin/clang++-4.0" $sDashFort -noomp -cuda -make-j -verbose -nobashrc || { exit 1; }
    else
        ./install ./TasInstall ./tsgWorkFolder -cmake="-DCMAKE_CXX_COMPILER=/usr/bin/clang++-4.0" $sDashFort -cuda -make-j -verbose -nobashrc || { exit 1; }
    fi
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    if (( $bMacOS == 0 )); then
        if [ -z `./TasInstall/bin/tasgrid -v | grep 'OpenMP multithreading: Enabled'` ]; then
            echo "OpenMP is supposed to work with clang++-4.0, but it failed somehow"
            exit 1;
        fi
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: Clang 4.0" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: Clang 4.0" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi

if [ -f /usr/bin/g++-7 ]; then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -cmake="-DCMAKE_CXX_COMPILER=/usr/bin/g++-7" $sDashFort -cuda -make-j -verbose -nobashrc || { exit 1; }
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: GCC 7" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: GCC 7" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi

if [ -f /usr/bin/g++-6 ]; then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -cmake="-DCMAKE_CXX_COMPILER=/usr/bin/g++-6" $sDashFort -cuda -make-j -verbose -nobashrc || { exit 1; }
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: GCC 6" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: GCC 6" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi

if [ -f /usr/bin/g++-5 ]; then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -cmake="-DCMAKE_CXX_COMPILER=/usr/bin/g++-5" $sDashFort -cuda -make-j -verbose -nobashrc || { exit 1; }
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: GCC 5" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: GCC 5" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi

if [ -d /usr/local/cuda-8.0 ]; then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -cmake="-DCUDA_TOOLKIT_ROOT_DIR:PATH=/usr/local/cuda-8.0/" $sDashFort -cuda -make-j -verbose -nobashrc || { exit 1; }
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: CUDA 8.0" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: CUDA 8.0" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi

if [ -d /usr/local/cuda-9.0 ]; then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -cmake="-DCUDA_TOOLKIT_ROOT_DIR:PATH=/usr/local/cuda-9.0/" $sDashFort -cuda -make-j -verbose -nobashrc || { exit 1; }
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: CUDA 9.0" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: CUDA 9.0" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi

if [ -d /usr/local/cuda-9.1 ]; then
    cp -r $sTempSource $sTempBuild/Tasmanian || { exit 1; }
    cd $sTempBuild/Tasmanian || { exit 1; }
    mkdir -p tsgWorkFolder
    ./install ./TasInstall ./tsgWorkFolder -cmake="-DCUDA_TOOLKIT_ROOT_DIR:PATH=/usr/local/cuda-9.1/" $sDashFort -cuda -make-j -verbose -nobashrc || { exit 1; }
    if (( $bOctave == 1 )); then
        octave --eval "addpath('$sTempBuild/Tasmanian/TasInstall/matlab/'); tsgCoreTests()" || { exit 1; }
    fi
    cd $sTempBuild
    rm -fr Tasmanian/
    cd $sTestRoot
    echo "======= PASSED: CUDA 9.1" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
else
    echo "======= SKIPPED: CUDA 9.1" >> $sMultibuildLogFile
    echo "===========================================================================================" >> $sMultibuildLogFile
fi


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

