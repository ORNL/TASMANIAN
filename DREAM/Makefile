include ../Config/AltBuildSystems/Makefile.in

IADD = -I../SparseGrids $(CommonIADD)
LADD = -L../ -L. $(CommonLADD)
LIBS = ../libtasmaniansparsegrid.a $(CommonLIBS)


LHEADERS = TasmanianDREAM.hpp tsgDreamState.hpp tsgDreamSample.hpp tsgDreamLikelihoodCore.hpp \
           tsgDreamLikelyGaussian.hpp tsgDreamCoreRandom.hpp \
           tsgDreamCorePDF.hpp tsgDreamEnumerates.hpp

LIBOBJ = tsgDreamState.o tsgDreamLikelyGaussian.o tsgDreamSampleWrapC.o

WROBJ = dreamtest_main.o tasdreamExternalTests.o

LIBNAME = libtasmaniandream.a
SHAREDNAME = libtasmaniandream.so

EXECNAME = dreamtest


%.o: %.cpp $(LHEADERS)
	$(CC) $(OPTC) $(IADD) -c $< -o $@

all: ../$(LIBNAME) ../$(EXECNAME) ../$(SHAREDNAME)

../$(SHAREDNAME): $(LIBOBJ) ../$(LIBNAME) ../$(EXECNAME)
	cd ..; $(CC) $(OPTL) $(LADD) ./DREAM/tsg*.o -shared -o $(SHAREDNAME) ./libtasmaniansparsegrid.so $(CommonLIBS)

../$(LIBNAME): $(LIBOBJ)
	ar rcs ../$(LIBNAME) $(LIBOBJ)

../$(EXECNAME): ../$(LIBNAME) $(WROBJ)
	$(CC) $(OPTL) $(LADD) $(WROBJ) -o ../$(EXECNAME) ../$(LIBNAME) $(LIBS)


clean:
	rm -fr *.o
