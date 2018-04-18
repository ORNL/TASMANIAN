include ../Config/AltBuildSystems/Makefile.in

IADD = -I../SparseGrids $(CommonIADD)
LADD = -L../ $(CommonLADD)
LIBS = ../libtasmaniansparsegrid.a $(CommonLIBS)


LHEADERS = tdrCorePDF.hpp tdrEnumerates.hpp tasdreamExternalTests.hpp tasdreamTestPDFs.hpp

LIBOBJ = TasmanianDREAM.o tdrCorePDF.o tasdreamExternalTests.o tasdreamTestPDFs.o

WROBJ = tasdream_main.o tasdreamBenchmark.o

LIBNAME = libtasmaniandream.a
SHAREDNAME = libtasmaniandream.so

EXECNAME = tasdream


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
