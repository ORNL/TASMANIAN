#include Config/AltBuildSystems/Makefile.in

# disable openmp on non-Linux platforms
UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
OPENMPFLAGS = -fopenmp
else
OPENMPFLAGS = -Wno-unknown-pragmas
endif

# Default C++ compiler
CC = g++
CXXFLAGS = -O3 -std=c++11 $(OPENMPFLAGS) -fPIC

IADD = -I./include -I./SparseGrids/ -I./InterfaceTPL/ -I./DREAM/ -I./Addons/ -I./
LADD = -L./

# set the compile command
COMPILE = $(CC) $(CXXFLAGS) $(IADD) $(LADD)
ECOMPILE = $(CC) $(CXXFLAGS) -I./include

# check if Python exists, disable python with TSG_SKIP_PYTHON=1
ifneq ($(shell echo "import numpy; print(numpy.__version__ )" | python3 2>/dev/null),)
HAS_PYTHON := "Yes"
endif

# object files for each target
SparseGridsObj = \
          ./SparseGrids/tsgIndexSets.o ./SparseGrids/tsgCoreOneDimensional.o ./SparseGrids/tsgIndexManipulator.o ./SparseGrids/tsgGridGlobal.o \
          ./SparseGrids/tsgSequenceOptimizer.o ./SparseGrids/tsgLinearSolvers.o ./SparseGrids/tsgGridSequence.o ./SparseGrids/tsgHardCodedTabulatedRules.o \
          ./SparseGrids/tsgHierarchyManipulator.o ./SparseGrids/tsgGridLocalPolynomial.o ./SparseGrids/tsgRuleWavelet.o ./SparseGrids/tsgGridWavelet.o \
          ./SparseGrids/tsgGridFourier.o ./SparseGrids/tsgDConstructGridGlobal.o ./SparseGrids/tsgAcceleratedDataStructures.o \
          ./SparseGrids/TasmanianSparseGridWrapC.o ./SparseGrids/TasmanianSparseGrid.o ./InterfaceTPL/tsgGpuNull.o

DREAMObj = \
          ./DREAM/tsgDreamState.o ./DREAM/tsgDreamLikelyGaussian.o ./DREAM/tsgDreamSampleWrapC.o

AddonsObj = \
          ./Addons/tsgCLoadNeededValues.o ./Addons/tsgCConstructSurrogate.o ./Addons/tsgCLoadUnstructuredPoints.o ./Addons/tsgCExoticQuadrature.o

TasgridObj = \
          ./Tasgrid/tasgrid_main.o ./SparseGrids/gridtestTestFunctions.o ./SparseGrids/gridtestExternalTests.o ./Tasgrid/tasgridWrapper.o

# test objects
GidtestObj = \
          ./SparseGrids/gridtest_main.o ./SparseGrids/gridtestTestFunctions.o ./SparseGrids/gridtestExternalTests.o \
          ./SparseGrids/gridtestUnitTests.o ./SparseGrids/gridtestTestInterfaceC.o

DREAMTestObj = \
          ./DREAM/dreamtest_main.o ./DREAM/tasdreamExternalTests.o

AddonTestObj = \
          ./Addons/testAddons.o

# example objects
SGExampleObj = \
          ./SparseGrids/Examples/example_sparse_grids.o     ./SparseGrids/Examples/example_sparse_grids_01.o \
          ./SparseGrids/Examples/example_sparse_grids_02.o  ./SparseGrids/Examples/example_sparse_grids_03.o \
          ./SparseGrids/Examples/example_sparse_grids_04.o  ./SparseGrids/Examples/example_sparse_grids_05.o \
          ./SparseGrids/Examples/example_sparse_grids_06.o  ./SparseGrids/Examples/example_sparse_grids_07.o \
          ./SparseGrids/Examples/example_sparse_grids_08.o  ./SparseGrids/Examples/example_sparse_grids_09.o \
          ./SparseGrids/Examples/example_sparse_grids_10.o

DREAMExampleObj = \
          ./DREAM/Examples/example_dream.o     ./DREAM/Examples/example_dream_01.o \
          ./DREAM/Examples/example_dream_02.o  ./DREAM/Examples/example_dream_03.o \
          ./DREAM/Examples/example_dream_04.o  ./DREAM/Examples/example_dream_05.o

# all target
.PHONY: all
all: tasgrid
	@echo ""
	@echo "TASMANIAN Simple GNU Make engine:"
	@echo "Complete build for C/C++, CLI and Python interfaces"
	@echo "MATLAB interface uses 'make matlab'"


# help
.PHONY: help
help:
	@echo "TASMANIAN Simple GNU Make tool"
	@echo " make"
	@echo " make test"
	@echo " make examples"
	@echo " make test_examples"
	@echo " make matlab"


# linear dependence pattern: sparse-grid -> dream -> addons -> python -> tasgrid
tasgrid: Tasmanian.py $(TasgridObj)
	@echo " -- building the tasgrid executable"
	$(COMPILE) $(TasgridObj) $(SparseGridsObj) -o tasgrid

Tasmanian.py: libtasmaniancaddons.so
	@echo " -- building the Python module"
	cp ./InterfacePython/TasmanianConfig.in.py TasmanianConfig.py
	sed -i -e 's|@Tasmanian_VERSION_MAJOR@|'7'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_VERSION_MINOR@|'8'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_license@|'BSD\ 3-Clause\ with\ UT-Battelle\ disclaimer'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_git_hash@|'Tasmanian\ git\ hash\ is\ not\ available\ here'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_libsparsegrid_path@|'`pwd`/libtasmaniansparsegrid.so'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_libdream_path@|'`pwd`/libtasmaniandream.so'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_libcaddons_path@|'`pwd`/libtasmaniancaddons.so'|g' ./TasmanianConfig.py
	cp ./InterfacePython/Tasmanian* .

libtasmaniancaddons.so: libtasmaniandream.so $(AddonsObj)
	@echo " -- building the Addons module"
	cp ./Addons/*.hpp ./include/
	$(COMPILE) -shared $(AddonsObj) -o libtasmaniancaddons.so ./libtasmaniandream.so ./libtasmaniansparsegrid.so

libtasmaniandream.so: libtasmaniansparsegrid.so $(DREAMObj)
	@echo " -- building the DREAM module"
	cp ./DREAM/*.hpp ./include/
	ar rcs libtasmaniandream.a $(DREAMObj)
	$(COMPILE) -shared $(DREAMObj) -o libtasmaniandream.so ./libtasmaniansparsegrid.so

libtasmaniansparsegrid.so: $(SparseGridsObj)
	@echo " -- building the Sparse Grid module"
	cp ./SparseGrids/*.h ./include/
	cp ./SparseGrids/*.hpp ./include/
	ar rcs libtasmaniansparsegrid.a $(SparseGridsObj)
	$(COMPILE) -shared $(SparseGridsObj) -o libtasmaniansparsegrid.so

%.o: %.cpp ./include
	$(COMPILE) -c $< -o $@

include:
	mkdir -p ./include
	cp ./Config/AltBuildSystems/TasmanianConfig.hpp ./include/TasmanianConfig.hpp
	cp ./Config/AltBuildSystems/tasgridLogs.hpp ./include/tasgridLogs.hpp
	cp ./Config/Tasmanian.hpp ./include/


# test target, also compiles the tests
.PHONY: test
test: gridtest dreamtest addontester
	./gridtest
	./dreamtest
	./addontester
ifdef HAS_PYTHON
	cp ./SparseGrids/GaussPattersonRule.table .
	cp ./Config/AltBuildSystems/testConfigureData.py InterfacePython/
	PYTHONPATH=$(PYTHONPATH):./InterfacePython python3 ./InterfacePython/testTSG.py && { echo "SUCCESS: Test completed successfully"; }
else
	@echo "WARNING: Cannot find python3 with numpy in the system path, skipping the Python tests."
endif

addontester: $(AddonTestObj)
	$(COMPILE) $(AddonTestObj) -o addontester libtasmaniansparsegrid.a libtasmaniandream.a

dreamtest: $(DREAMTestObj)
	$(COMPILE) $(DREAMTestObj) -o dreamtest libtasmaniansparsegrid.a libtasmaniandream.a

gridtest: $(GidtestObj)
	$(COMPILE) $(GidtestObj) -o gridtest libtasmaniansparsegrid.a libtasmaniandream.a


# examples
.PHONY: examples
examples: example_sparse_grids example_dream
	cp ./InterfacePython/example_sparse_grids.in.py ./InterfacePython/example_sparse_grids.py
	sed -i -e 's|@Tasmanian_string_python_hashbang@|'\/usr\/bin\/env\ python3'|g' ./InterfacePython/example_sparse_grids.py
	sed -i -e 's|@Tasmanian_python_example_import@|'sys.path.append\(\"`pwd`\"\)'|g' ./InterfacePython/example_sparse_grids.py
	cp ./InterfacePython/example_dream.in.py ./InterfacePython/example_dream.py
	sed -i -e 's|@Tasmanian_string_python_hashbang@|'\/usr\/bin\/env\ python3'|g' ./InterfacePython/example_dream.py
	sed -i -e 's|@Tasmanian_python_example_import@|'sys.path.append\(\"`pwd`\"\)'|g' ./InterfacePython/example_dream.py
	@echo "Compiled the examples, run with the commands"
	@echo "./example_sparse_grids"
	@echo "./example_dream"
	@echo "python3 ./InterfacePython/example_sparse_grids.py"
	@echo "python3 ./InterfacePython/example_dream.py"

example_sparse_grids: $(SGExampleObj)
	$(ECOMPILE) $(SGExampleObj) -o example_sparse_grids ./libtasmaniansparsegrid.a

example_dream: $(DREAMExampleObj)
	$(ECOMPILE) $(DREAMExampleObj) -o example_dream ./libtasmaniansparsegrid.a ./libtasmaniandream.a

SparseGrids/Examples/example_sparse_grids%.o: SparseGrids/Examples/example_sparse_grids%.cpp
	$(ECOMPILE) -c $< -o $@

DREAM/Examples/example_dream%.o: DREAM/Examples/example_dream%.cpp
	$(ECOMPILE) -c $< -o $@


# examples testing
.PHONY: test_examples
test_examples: examples
	./example_sparse_grids
	./example_dream
ifdef HAS_PYTHON
	python3 ./InterfacePython/example_sparse_grids.py
	python3 ./InterfacePython/example_dream.py
endif


# Matlab
.PHONY: matlab
matlab:
	cp Config/AltBuildSystems/tsgGetPaths.m InterfaceMATLAB/tsgGetPaths.m
	mkdir -p tsgMatlabWorkFolder
	sed -i -e 's|ENTER\ THE\ PATH\ TO\ MATLAB\ WORK\ FOLDER|'`pwd`/tsgMatlabWorkFolder/'|g' InterfaceMATLAB/tsgGetPaths.m
	sed -i -e 's|ENTER\ THE\ PATH\ TO\ tasgrid\ EXECUTABLE|'`pwd`/tasgrid'|g' InterfaceMATLAB/tsgGetPaths.m
	@echo ""
	@echo "Setting MATLAB work forlder to "`pwd`"/tsgMatlabWorkFolder/"
	@echo "use MATLAB command: addpath('"`pwd`"/InterfaceMATLAB/');"
	@echo ""


# Clean
.PHONY: clean
clean:
	rm -fr libtasmaniansparsegrid.so
	rm -fr libtasmaniansparsegrid.a
	rm -fr libtasmaniandream.so
	rm -fr libtasmaniandream.a
	rm -fr libtasmaniancaddons.so
	rm -fr tasgrid
	rm -fr gridtest
	rm -fr dreamtest
	rm -fr addontester
	rm -fr Tasmanian*.py
	rm -fr *.pyc
	rm -fr GaussPattersonRule.table
	rm -fr __pycache__
	rm -fr ./InterfacePython/__pycache__
	rm -fr ./InterfacePython/*.pyc
	rm -fr ./InterfacePython/example_sparse_grids.py
	rm -fr ./InterfacePython/example_dream.py
	rm -fr ./InterfacePython/testConfigureData.py
	rm -fr InterfaceMATLAB/tsgGetPaths.m
	rm -fr example_sparse_grids
	rm -fr example_dream
	rm -fr include
	rm -fr ./SparseGrids/*.o
	rm -fr ./DREAM/*.o
	rm -fr ./Addons/*.o
	rm -fr ./Tasgrid/*.o
	rm -fr ./InterfaceTPL/*.o
	rm -fr ./SparseGrids/Examples/*.o
	rm -fr ./DREAM/Examples/*.o
	rm -fr ./tsgMatlabWorkFolder
	rm -fr ./testSave
	rm -fr ./refTestFlename.grid
