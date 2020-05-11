include Config/AltBuildSystems/Makefile.in

IADD = -I./include $(CommonIADD)
LADD = -L./ $(CommonLADD)
LIBS = ./libtasmaniandream.a ./libtasmaniansparsegrid.a $(CommonLIBS)
FFLIBS90 = ./libtasmanianfortran90.a ./libtasmaniansparsegrid.a ./libtasmaniandream.a $(CommonLIBS)

TSG_SOURCE = $(wildcard ./SparseGrids/*.cpp) $(wildcard ./SparseGrids/*.hpp) $(wildcard ./SparseGrids/*.h)
TDR_SOURCE = $(wildcard ./DREAM/*.cpp) $(wildcard ./DREAM/*.hpp)

# headers to be ignored in this build .in., .windows., executables, config files, etc.
CMAKE_IN_HEADERS = $(wildcard ./DREAM/*.in.hpp) \
                   $(wildcard ./DREAM/tasdream*.hpp) \
                   $(wildcard ./SparseGrids/*.*.hpp) \
                   $(wildcard ./SparseGrids/*.*.h) \
                   $(wildcard ./SparseGrids/tasgrid*.hpp) \
                   $(wildcard ./Addons/test*.hpp) \

# header names with ./include as opposed to ./SparseGrids
HEADERS = $(patsubst ./DREAM/%,./include/%,$(filter-out $(CMAKE_IN_HEADERS),$(wildcard ./DREAM/*.hpp))) \
          $(patsubst ./SparseGrids/%,./include/%,$(filter-out $(CMAKE_IN_HEADERS),$(wildcard ./SparseGrids/*.h*))) \
          $(patsubst ./Addons/%,./include/%,$(filter-out $(CMAKE_IN_HEADERS),$(wildcard ./Addons/*.hpp))) \
          ./include/TasmanianConfig.hpp ./include/Tasmanian.hpp

ALL_TARGETS = GaussPattersonRule.table \
              Tasmanian.py TasmanianConfig.py \
              example_sparse_grids.py InterfacePython/testConfigureData.py testTSG.py sandbox.py \
              $(wildcard ./DREAM/Examples/example_dream*.cpp) \
              $(wildcard ./SparseGrids/Examples/example_sparse_grids*.cpp) \
              libtasmaniansparsegrid.so libtasmaniandream.so \
              libtasmaniancaddons.so $(HEADERS)

DREAM_EXAMPLES_OBJ = $(patsubst ./DREAM/Examples/%,%,$(patsubst %.cpp,%.o,$(wildcard ./DREAM/Examples/example_dream*.cpp)))
SG_EXAMPLES_OBJ = $(patsubst ./SparseGrids/Examples/%,%,$(patsubst %.cpp,%.o,$(wildcard ./SparseGrids/Examples/example_sparse_grids*.cpp)))

CONFIGURED_HEADERS = ./SparseGrids/TasmanianConfig.hpp ./SparseGrids/tasgridLogs.hpp

# all target
.PHONY: all
all: $(ALL_TARGETS)
	@echo ""
	@echo "TASMANIAN Simple GNU Make engine:"
	@echo "Complete build for C/C++, CLI and Python interfaces"
	@echo "MATLAB interface: use 'make matlab'"
	@echo ""
	@echo "      or alternatively manually edit ./InterfaceMATLAB/tsgGetPaths.m"
	@echo "      and use MATLAB command addpath('$(shell pwd)/InterfaceMATLAB/')"
	@echo ""

# There are seemingly superfluous dependencies, but those allow for the use of make -j command without conflict
# most of the work is done by the make in the subfolders and
# the extra dependencies ensure only one subfoler make is called at a time

# Common headers
SparseGrids/TasmanianConfig.hpp: ./Config/AltBuildSystems/TasmanianConfig.hpp
	cp ./Config/AltBuildSystems/TasmanianConfig.hpp ./SparseGrids/TasmanianConfig.hpp

include/TasmanianConfig.hpp: ./SparseGrids/TasmanianConfig.hpp
	mkdir -p ./include
	cp ./SparseGrids/TasmanianConfig.hpp ./include/

./SparseGrids/tasgridLogs.hpp: ./Config/AltBuildSystems/tasgridLogs.hpp
	cp ./Config/AltBuildSystems/tasgridLogs.hpp ./SparseGrids/tasgridLogs.hpp

include/Tasmanian.hpp: ./Config/Tasmanian.hpp
	mkdir -p ./include
	cp ./Config/Tasmanian.hpp ./include/

# Sparse Grids
libtasmaniansparsegrid.so: $(TSG_SOURCE) $(CONFIGURED_HEADERS)
	cd SparseGrids; make

# DREAM
libtasmaniandream.so: libtasmaniansparsegrid.so $(TDR_SOURCE) $(CONFIGURED_HEADERS)
	cd DREAM; make

# Addons
libtasmaniancaddons.so: libtasmaniandream.so libtasmaniansparsegrid.so
	cd Addons; make

# Headers
# many calls to mkdir, consider reducing
include/TasmanianDREAM.hpp: ./DREAM/TasmanianDREAM.hpp
	mkdir -p ./include
	cp ./DREAM/TasmanianDREAM.hpp ./include/

include/tsg%.hpp: ./DREAM/tsg%.hpp
	mkdir -p ./include
	cp ./$< ./$@

include/TasmanianSparseGrid.h: ./SparseGrids/TasmanianSparseGrid.h
	mkdir -p ./include
	cp ./SparseGrids/TasmanianSparseGrid.h ./include/

include/TasmanianSparseGrid.hpp: ./SparseGrids/TasmanianSparseGrid.hpp
	mkdir -p ./include
	cp ./SparseGrids/TasmanianSparseGrid.hpp ./include/

include/tsg%.hpp: ./SparseGrids/tsg%.hpp
	mkdir -p ./include
	cp ./$< ./$@

include/TasmanianAddons.hpp: ./Addons/TasmanianAddons.hpp
	mkdir -p ./include
	cp ./Addons/TasmanianAddons.hpp ./include/

include/tsg%.hpp: ./Addons/tsg%.hpp
	mkdir -p ./include
	cp ./$< ./$@

# Python, tables, etc.
GaussPattersonRule.table: ./SparseGrids/GaussPattersonRule.table
	cp ./SparseGrids/GaussPattersonRule.table .

Tasmanian.py: ./InterfacePython/Tasmanian.py
	cp ./InterfacePython/Tasmanian* .

TasmanianConfig.py: ./InterfacePython/TasmanianConfig.in.py
	cp ./InterfacePython/TasmanianConfig.in.py TasmanianConfig.py
	sed -i -e 's|@Tasmanian_VERSION_MAJOR@|'7'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_VERSION_MINOR@|'1'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_license@|'BSD\ 3-Clause\ with\ UT-Battelle\ disclaimer'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_git_hash@|'Tasmanian\ git\ hash\ is\ not\ available\ here'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_libsparsegrid_path@|'`pwd`/libtasmaniansparsegrid.so'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_libdream_path@|'`pwd`/libtasmaniandream.so'|g' ./TasmanianConfig.py
	sed -i -e 's|@Tasmanian_libcaddons_path@|'`pwd`/libtasmaniancaddons.so'|g' ./TasmanianConfig.py

example_sparse_grids.py: ./InterfacePython/example_sparse_grids.in.py
	cp ./InterfacePython/example_sparse_grids.in.py example_sparse_grids.py
	sed -i -e 's|@Tasmanian_string_python_hashbang@|'\/usr\/bin\/env\ python'|g' ./example_sparse_grids.py
	sed -i -e 's|@Tasmanian_python_example_import@|'sys.path.append\(\"`pwd`\"\)'|g' ./example_sparse_grids.py
	cp ./InterfacePython/example_dream.in.py example_dream.py
	sed -i -e 's|@Tasmanian_string_python_hashbang@|'\/usr\/bin\/env\ python'|g' ./example_dream.py
	sed -i -e 's|@Tasmanian_python_example_import@|'sys.path.append\(\"`pwd`\"\)'|g' ./example_dream.py
	cp ./InterfacePython/example_sparse_grids_* .
	cp ./InterfacePython/example_dream_* .

InterfacePython/testConfigureData.py: ./Config/AltBuildSystems/testConfigureData.py
	cp ./Config/AltBuildSystems/testConfigureData.py InterfacePython/

testTSG.py: ./InterfacePython/testTSG.py
	cp ./InterfacePython/testTSG.py .

sandbox.py: ./InterfacePython/sandbox.py
	cp ./InterfacePython/sandbox.py .

example_spa%.cpp: ./SparseGrids/Examples/example_spa%.cpp
	cp ./$< ./$@

example_dre%.cpp: ./DREAM/Examples/example_dre%.cpp
	cp ./$< ./$@

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

# Python 3
.PHONY: python3
python3: TasmanianSG.py testTSG.py example_sparse_grids.py
	sed -i -e 's|\#\!\/usr\/bin\/env\ python|\#\!\/usr\/bin\/env\ python3|g' example_sparse_grids.py
	sed -i -e 's|\#\!\/usr\/bin\/env\ python|\#\!\/usr\/bin\/env\ python3|g' example_dream.py
	sed -i -e 's|\#\!\/usr\/bin\/env\ python|\#\!\/usr\/bin\/env\ python3|g' testTSG.py

# Fortran
#fortran: example_sparse_grids_f90 libtasmanianfortran90.a libtasmanianfortran90.so
.PHONY: fortran
fortran: ./include/tasmaniansg.mod forswigtester example_sparse_grids_fortran
	./forswigtester

./include/tasmaniansg.mod: libtasmaniancaddons.so libtasmaniandream.a
	$(CC) $(OPTC) $(IADD) -I./include/ -c InterfaceSwig/generated/tasmanianFORTRAN_wrap.cxx -o InterfaceSwig/tasmanianFORTRAN_wrap.o
	$(FF) $(OPTF) $(IADD) -c InterfaceSwig/generated/tasmanian.f90 -o InterfaceSwig/tasmanian.o
	ar rcs libtasmanianfortran.a InterfaceSwig/tasmanianFORTRAN_wrap.o InterfaceSwig/tasmanian.o
	$(FF) $(OPTLFF) $(LADD) InterfaceSwig/tasmanianFORTRAN_wrap.o InterfaceSwig/tasmanian.o -shared -o libtasmanianfortran.so ./libtasmaniandream.so ./libtasmaniansparsegrid.so $(LIBS) -lstdc++
	mv tasmanian.mod ./include/

forswigtester: ./include/tasmaniansg.mod
	$(FF) $(OPTLFF) $(LADD) InterfaceSwig/FortranTests/*.f90 -o forswigtester -I./include/ libtasmanianfortran.a $(LIBS) -lstdc++

example_sparse_grids_fortran: ./include/tasmaniansg.mod
	$(FF) $(OPTLFF) $(LADD) InterfaceSwig/FortranExamples/*.f90 -o example_sparse_grids_fortran -I./include/ libtasmanianfortran.a $(LIBS) -lstdc++


# Testing and examples
.PHONY: test
test: $(ALL_TARGETS)
	./gridtest
	./dreamtest
	./addontester
	PYTHONPATH=$(PYTHONPATH):./InterfacePython ./testTSG.py && { echo "SUCCESS: Test completed successfully"; }
	./example_sparse_grids.py

.PHONY: examples
examples: $(ALL_TARGETS) example_dream example_sparse_grids
	echo "Done examples"

example_dream: $(ALL_TARGETS) $(DREAM_EXAMPLES_OBJ)
	$(CC) $(OPTL) $(LADD) $(DREAM_EXAMPLES_OBJ) -o example_dream $(LIBS)

example_sparse_grids: $(ALL_TARGETS) $(SG_EXAMPLES_OBJ)
	$(CC) $(OPTL) $(LADD) $(SG_EXAMPLES_OBJ) -o example_sparse_grids $(LIBS)

%.o: %.cpp
	$(CC) $(OPTC) $(IADD) -c $<

# Clean
.PHONY: clean
clean:
	rm -fr libtasmaniansparsegrid.so
	rm -fr libtasmaniansparsegrid.a
	rm -fr libtasmaniandream.so
	rm -fr libtasmaniandream.a
	rm -fr libtasmaniancaddons.so
	rm -fr libtasmanianfortran.so
	rm -fr libtasmanianfortran.a
	rm -fr InterfaceSwig/*.o
	rm -fr forswigtester
	rm -fr f03_test_file
	rm -fr tasgrid
	rm -fr gridtest
	rm -fr dreamtest
	rm -fr addontester
	rm -fr TasmanianSG.py
	rm -fr TasmanianAddons.py
	rm -fr example_sparse_grids.py
	rm -fr GaussPattersonRule.table
	rm -fr *.pyc
	rm -fr ./InterfacePython/*.pyc
	rm -fr __pycache__
	rm -fr ./InterfacePython/__pycache__
	rm -fr example_sparse_grids*
	rm -fr example_dream*
	rm -fr example_sparse_grids_fortran
	rm -fr refTestFlename.grid
	rm -fr fortranio.test
	rm -fr testSave
	rm -fr testTSG.py
	rm -fr sandbox.py
	rm -fr *.py
	rm -fr include
	rm -fr tasmanianFortranTestCustomRule.table
	rm -fr ./SparseGrids/TasmanianConfig.hpp
	rm -fr ./SparseGrids/tasgridLogs.hpp
	rm -fr ./InterfacePython/testConfigureData.py
	cd SparseGrids; make clean
	cd DREAM; make clean
	cd Addons; make clean
	cd InterfaceFortran; make clean
