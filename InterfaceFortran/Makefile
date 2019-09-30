include ../Config/AltBuildSystems/Makefile.in

IADD = -I../SparseGrids -I../DREAM -I../Addons $(CommonIADD)
LADD = -L../ $(CommonLADD)
LIBS = ../libtasmaniansparsegrid.a ../libtasmaniandream.a $(CommonLIBS) -lstdc++

CSOURCE = tsgC2Fortran.cpp
FSOURCE = TasmanianSG.f90 tsgC2FortranBridge.f90

LIBNAME = libtasmanianfortran90.a
SHAREDNAME = libtasmanianfortran90.so

TESTNAME = fortester90

all: ../$(LIBNAME) ../$(SHAREDNAME) ../$(TESTNAME)

../$(LIBNAME): tsgC2Fortran.o TasmanianSG.o tsgC2FortranBridge.o
	ar rcs ../$(LIBNAME) tsgC2Fortran.o TasmanianSG.o tsgC2FortranBridge.o

../$(SHAREDNAME): tsgC2Fortran.o TasmanianSG.o tsgC2FortranBridge.o
	$(FF) $(OPTLFF) $(LADD) tsgC2Fortran.o TasmanianSG.o tsgC2FortranBridge.o -shared -o ../$(SHAREDNAME) ../libtasmaniansparsegrid.so $(LIBS)

tsgC2Fortran.o: tsgC2Fortran.cpp
	$(CC) $(OPTC) $(IADD) -c tsgC2Fortran.cpp -o tsgC2Fortran.o

TasmanianSG.o: TasmanianSG.f90
	$(FF) $(OPTF) $(IADD) -c TasmanianSG.f90 -o TasmanianSG.o

tsgC2FortranBridge.o: tsgC2FortranBridge.f90
	$(FF) $(OPTF) $(IADD) -c tsgC2FortranBridge.f90 -o tsgC2FortranBridge.o

../$(TESTNAME): ../$(LIBNAME) tasmaniansg.mod
	$(FF) $(OPTF) $(IADD) -c fortester.f90 -o fortester90.o
	$(FF) $(OPTLFF) $(LADD) fortester90.o -o ../$(TESTNAME) ../$(LIBNAME) $(LIBS)

# clean
.PHONY: clean
clean:
	rm -fr *.o
	rm -fr *.mod
	rm -fr *.so
	rm -fr *.a
