include ../Config/Makefile.in

IADD = -I../SparseGrids $(CommonIADD)
LADD = -L../ $(CommonLADD)
LIBS = ../libtasmaniansparsegrid.a $(CommonLIBS)


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


LHEADERS = tdrCorePDF.hpp tdrEnableMPI.hpp tdrEnumerates.hpp tasdreamExternalTests.hpp tasdreamTestPDFs.hpp 

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
