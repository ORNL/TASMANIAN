include ../Config/AltBuildSystems/Makefile.in

IADD = -I../SparseGrids -I../DREAM $(CommonIADD)
LADD = -L../ $(CommonLADD)
LIBS = ./libtasmaniansparsegrid.so ./libtasmaniandream.so $(CommonLIBS)

LHEADERS = $(wildcard ./*.hpp)

LIBOBJ = tsgCLoadNeededValues.o tsgCConstructSurrogate.o tsgCLoadUnstructuredPoints.o tsgCExoticQuadrature.o

MODULENAME = libtasmaniancaddons.so

TESTNAME = addontester

all: $(LIBOBJ) ../$(MODULENAME) ../$(TESTNAME)

%.o: %.cpp $(LHEADERS)
	$(CC) $(OPTC) $(IADD) -c $< -o $@

../$(MODULENAME): $(LIBOBJ)
	cd ..; $(CC) $(OPTL) $(LADD) ./Addons/tsg*.o -shared -o $(MODULENAME) $(LIBS)

../$(TESTNAME): testAddons.o
	$(CC) $(OPTL) $(LADD) testAddons.o -o ../$(TESTNAME) ../libtasmaniandream.a ../libtasmaniansparsegrid.a

clean:
	rm -fr *.o
