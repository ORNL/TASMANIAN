include ../Config/Makefile.in

IADD = $(CommonIADD)
LADD = $(CommonLADD)
LIBS = $(CommonLIBS)

# Uncomment below to enable clang + OpenMP
# Note: this needs libiomp to be symlinked as libomp in the current folder
#CC = clang++-3.7
#COMPILE_OPTIONS = -fopenmp=libomp
# End cland + OpenMP

# Uncomment below for intel icc compiler
#CC = icc
#COMPILE_OPTIONS = -fopenmp -diag-disable 11074 -diag-disable 11076
#OPTC = -O3 -mtune=native $(COMPILE_OPTIONS) -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic
# End icc


LHEADERS = tsgIndexSets.hpp tsgCoreOneDimensional.hpp tsgIndexManipulator.hpp tsgGridGlobal.hpp tsgCacheLagrange.hpp tsgSequenceOptimizer.hpp \
           tsgEnumerates.hpp tsgOneDimensionalWrapper.hpp tsgGridSequence.hpp tsgGridCore.hpp tsgLinearSolvers.hpp \
           tsgRuleLocalPolynomial.hpp tsgHardCodedTabulatedRules.hpp tsgGridLocalPolynomial.hpp \
           tsgRuleWavelet.hpp tsgGridWavelet.hpp \
           tasgridTestFunctions.hpp tasgridExternalTests.hpp tasgridWrapper.hpp \
           tsgAcceleratedDataStructures.hpp \
           TasmanianSparseGrid.hpp

LIBOBJ = tsgIndexSets.o tsgCoreOneDimensional.o tsgIndexManipulator.o tsgGridGlobal.o tsgCacheLagrange.o tsgSequenceOptimizer.o tsgOneDimensionalWrapper.o \
         tsgGridCore.o tsgLinearSolvers.o tsgGridSequence.o tsgRuleLocalPolynomial.o tsgHardCodedTabulatedRules.o \
         tsgGridLocalPolynomial.o tsgRuleWavelet.o tsgGridWavelet.o \
         tsgAcceleratedDataStructures.o \
         TasmanianSparseGrid.o

WROBJ = tasgrid_main.o tasgridTestFunctions.o tasgridExternalTests.o tasgridWrapper.o

LIBNAME = libtasmaniansparsegrid.a
SHAREDNAME = libtasmaniansparsegrid.so

EXECNAME = tasgrid


%.o: %.cpp $(LHEADERS)
	$(CC) $(OPTC) $(IADD) -c $< -o $@

all: $(LIBNAME) $(EXECNAME) $(SHAREDNAME)

$(SHAREDNAME): $(LIBOBJ)
	$(CC) $(OPTL) $(LIBOBJ) -shared -o $(SHAREDNAME) $(LIBS)

$(LIBNAME): $(LIBOBJ)
	ar rcs $(LIBNAME) $(LIBOBJ)

$(EXECNAME):  $(LIBNAME) $(WROBJ)
	$(CC) $(OPTL) $(LADD) -L. $(WROBJ) -o $(EXECNAME) $(LIBNAME) $(LIBS)


clean:
	rm -fr *.o
	rm -fr $(LIBNAME)
	rm -fr $(EXECNAME)
	rm -fr $(SHAREDNAME)
