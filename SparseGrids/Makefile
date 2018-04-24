include ../Config/AltBuildSystems/Makefile.in

# use the common libs
IADD = $(CommonIADD)
LADD = $(CommonLADD)
LIBS = $(CommonLIBS)


LHEADERS = tsgIndexSets.hpp tsgCoreOneDimensional.hpp tsgIndexManipulator.hpp tsgGridGlobal.hpp tsgCacheLagrange.hpp tsgSequenceOptimizer.hpp \
           tsgEnumerates.hpp tsgOneDimensionalWrapper.hpp tsgGridSequence.hpp tsgGridCore.hpp tsgLinearSolvers.hpp \
           tsgRuleLocalPolynomial.hpp tsgHardCodedTabulatedRules.hpp tsgGridLocalPolynomial.hpp \
           tsgRuleWavelet.hpp tsgCudaMacros.hpp tsgGridWavelet.hpp \
           tsgCudaLinearAlgebra.hpp tsgCudaBasisEvaluations.hpp tsgAcceleratedDataStructures.hpp \
           tasgridTestFunctions.hpp tasgridExternalTests.hpp tasgridWrapper.hpp \
           TasmanianSparseGrid.hpp

LIBOBJ = tsgIndexSets.o tsgCoreOneDimensional.o tsgIndexManipulator.o tsgGridGlobal.o tsgSequenceOptimizer.o tsgOneDimensionalWrapper.o \
         tsgGridCore.o tsgLinearSolvers.o tsgGridSequence.o tsgRuleLocalPolynomial.o tsgHardCodedTabulatedRules.o \
         tsgGridLocalPolynomial.o tsgRuleWavelet.o tsgGridWavelet.o \
         tsgAcceleratedDataStructures.o $(TASMANIAN_CUDA_KERNELS) \
         TasmanianSparseGrid.o

WROBJ = tasgrid_main.o tasgridTestFunctions.o tasgridExternalTests.o tasgridWrapper.o

LIBNAME = libtasmaniansparsegrid.a
SHAREDNAME = libtasmaniansparsegrid.so

EXECNAME = tasgrid

%.cu.o: %.cu $(LHEADERS)
	$(NVCC) $(NVCC_OPT) -c $< -o $@

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
