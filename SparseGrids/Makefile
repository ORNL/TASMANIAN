include ../Config/AltBuildSystems/Makefile.in

# use the common libs
IADD = $(CommonIADD)
LADD = $(CommonLADD)
LIBS = $(CommonLIBS)


LHEADERS = tsgIndexSets.hpp tsgCoreOneDimensional.hpp tsgIndexManipulator.hpp tsgGridGlobal.hpp tsgCacheLagrange.hpp tsgSequenceOptimizer.hpp \
           tsgEnumerates.hpp tsgOneDimensionalWrapper.hpp tsgGridSequence.hpp tsgGridCore.hpp tsgLinearSolvers.hpp \
           tsgHierarchyManipulator.hpp \
           tsgRuleLocalPolynomial.hpp tsgHardCodedTabulatedRules.hpp tsgGridLocalPolynomial.hpp tsgGridFourier.hpp \
           tsgRuleWavelet.hpp tsgCudaLoadStructures.hpp tsgGridWavelet.hpp \
           tsgCudaLinearAlgebra.hpp tsgCudaBasisEvaluations.hpp tsgAcceleratedDataStructures.hpp \
           tsgDConstructGridGlobal.hpp \
           tasgridTestFunctions.hpp tasgridExternalTests.hpp tasgridWrapper.hpp tasgridUnitTests.hpp \
           TasmanianSparseGrid.hpp

LIBOBJ = tsgIndexSets.o tsgCoreOneDimensional.o tsgIndexManipulator.o tsgGridGlobal.o tsgSequenceOptimizer.o \
         tsgLinearSolvers.o tsgGridSequence.o tsgHardCodedTabulatedRules.o tsgHierarchyManipulator.o\
         tsgGridLocalPolynomial.o tsgRuleWavelet.o tsgGridWavelet.o tsgGridFourier.o \
         tsgDConstructGridGlobal.o \
         tsgAcceleratedDataStructures.o $(TASMANIAN_CUDA_KERNELS) \
         TasmanianSparseGridWrapC.o TasmanianSparseGrid.o tsgGpuNull.o

WROBJ = tasgrid_main.o tasgridTestFunctions.o tasgridExternalTests.o tasgridWrapper.o

GTONJ = gridtest_main.o tasgridTestFunctions.o tasgridExternalTests.o tasgridUnitTests.o tasgridTestInterfaceC.o

LIBNAME = libtasmaniansparsegrid.a
SHAREDNAME = libtasmaniansparsegrid.so

EXECNAME = tasgrid

TESTNAME = gridtest

%.o: %.cpp $(LHEADERS)
	$(CC) $(OPTC) $(IADD) -c $< -o $@

all: ../$(LIBNAME) ../$(EXECNAME) ../$(TESTNAME) ../$(SHAREDNAME)

tsgGpuNull.o: ../InterfaceTPL/tsgGpuNull.cpp
	$(CC) $(OPTC) $(IADD) -c ../InterfaceTPL/tsgGpuNull.cpp -o tsgGpuNull.o

../$(SHAREDNAME): $(LIBOBJ) ../$(LIBNAME) ../$(EXECNAME) ../$(TESTNAME)
	cd ..; $(CC) $(OPTL) ./SparseGrids/tsg*.o ./SparseGrids/Tasmanian*.o -shared -o $(SHAREDNAME) $(LIBS)

../$(LIBNAME): $(LIBOBJ)
	ar rcs ../$(LIBNAME) $(LIBOBJ)

../$(EXECNAME): ../$(LIBNAME) $(WROBJ)
	$(CC) $(OPTL) $(LADD) -L. $(WROBJ) -o ../$(EXECNAME) ../$(LIBNAME) $(LIBS)

../$(TESTNAME): ../$(LIBNAME) $(GTONJ)
	$(CC) $(OPTL) $(LADD) -L. $(GTONJ) -o ../$(TESTNAME) ../$(LIBNAME) $(LIBS)

clean:
	rm -fr *.o
