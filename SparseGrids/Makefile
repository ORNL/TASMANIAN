include ../Config/AltBuildSystems/Makefile.in

# use the common libs
IADD = $(CommonIADD)
LADD = $(CommonLADD)
LIBS = $(CommonLIBS)


LHEADERS = tsgIndexSets.hpp tsgCoreOneDimensional.hpp tsgIndexManipulator.hpp tsgGridGlobal.hpp tsgCacheLagrange.hpp tsgSequenceOptimizer.hpp \
           tsgEnumerates.hpp tsgOneDimensionalWrapper.hpp tsgGridSequence.hpp tsgGridCore.hpp tsgLinearSolvers.hpp \
           tsgRuleLocalPolynomial.hpp tsgHardCodedTabulatedRules.hpp tsgGridLocalPolynomial.hpp tsgGridFourier.hpp \
           tsgRuleWavelet.hpp tsgCudaMacros.hpp tsgGridWavelet.hpp \
           tsgCudaLinearAlgebra.hpp tsgCudaBasisEvaluations.hpp tsgAcceleratedDataStructures.hpp \
           tasgridTestFunctions.hpp tasgridExternalTests.hpp tasgridWrapper.hpp tasgridUnitTests.hpp \
           TasmanianSparseGrid.hpp

LIBOBJ = tsgIndexSets.o tsgCoreOneDimensional.o tsgIndexManipulator.o tsgGridGlobal.o tsgSequenceOptimizer.o tsgOneDimensionalWrapper.o \
         tsgGridCore.o tsgLinearSolvers.o tsgGridSequence.o tsgHardCodedTabulatedRules.o \
         tsgGridLocalPolynomial.o tsgRuleWavelet.o tsgGridWavelet.o tsgGridFourier.o \
         tsgAcceleratedDataStructures.o $(TASMANIAN_CUDA_KERNELS) \
         TasmanianSparseGrid.o

WROBJ = tasgrid_main.o tasgridTestFunctions.o tasgridExternalTests.o tasgridWrapper.o

GTONJ = gridtest_main.o tasgridTestFunctions.o tasgridExternalTests.o tasgridUnitTests.o tasgridTestInterfaceC.o

LIBNAME = libtasmaniansparsegrid.a
SHAREDNAME = libtasmaniansparsegrid.so

EXECNAME = tasgrid

TESTNAME = gridtest

%.cu.o: %.cu $(LHEADERS)
	$(NVCC) $(NVCC_OPT) -c $< -o $@

%.o: %.cpp $(LHEADERS)
	$(CC) $(OPTC) $(IADD) -c $< -o $@

all: $(LIBNAME) $(EXECNAME) $(TESTNAME) $(SHAREDNAME)

$(SHAREDNAME): $(LIBOBJ)
	$(CC) $(OPTL) $(LIBOBJ) -shared -o $(SHAREDNAME) $(LIBS)

$(LIBNAME): $(LIBOBJ)
	ar rcs $(LIBNAME) $(LIBOBJ)

$(EXECNAME):  $(LIBNAME) $(WROBJ)
	$(CC) $(OPTL) $(LADD) -L. $(WROBJ) -o $(EXECNAME) $(LIBNAME) $(LIBS)

$(TESTNAME):  $(LIBNAME) $(GTONJ)
	$(CC) $(OPTL) $(LADD) -L. $(GTONJ) -o $(TESTNAME) $(LIBNAME) $(LIBS)

clean:
	rm -fr *.o
	rm -fr $(LIBNAME)
	rm -fr $(EXECNAME)
	rm -fr $(TESTNAME)
	rm -fr $(SHAREDNAME)
