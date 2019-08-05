include ../Config/AltBuildSystems/Makefile.in

IADD = -I../SparseGrids $(CommonIADD)
LADD = -L../ $(CommonLADD)
LIBS = ../libtasmaniansparsegrid.a $(CommonLIBS)


LHEADERS = TasmanianDREAM.hpp tsgDreamState.hpp tsgDreamSample.hpp tsgDreamLikelihoodCore.hpp \
           tsgDreamLikelyGaussian.hpp tsgDreamInternalBlas.hpp tsgDreamCoreRandom.hpp \
           tsgDreamCorePDF.hpp tsgDreamEnumerates.hpp

LIBOBJ = tsgDreamState.o tsgDreamLikelyGaussian.o

WROBJ = dreamtest_main.o tasdreamExternalTests.o

LIBNAME = libtasmaniandream.a
SHAREDNAME = libtasmaniandream.so

EXECNAME = dreamtest


%.o: %.cpp $(LHEADERS)
	$(CC) $(OPTC) $(IADD) -c $< -o $@

all: $(LIBNAME) $(EXECNAME) $(SHAREDNAME)

$(SHAREDNAME): $(LIBOBJ)
	$(CC) $(OPTL) $(LADD) $(LIBOBJ) -shared -o $(SHAREDNAME) ../libtasmaniansparsegrid.so $(CommonLIBS)

$(LIBNAME): $(LIBOBJ)
	ar rcs $(LIBNAME) $(LIBOBJ)

$(EXECNAME):  $(LIBNAME) $(WROBJ)
	$(CC) $(OPTL) $(LADD) -L. $(WROBJ) -o $(EXECNAME) $(LIBNAME) $(LIBS)


clean:
	rm -fr *.o
	rm -fr $(LIBNAME)
	rm -fr $(EXECNAME)
	rm -fr $(SHAREDNAME)
