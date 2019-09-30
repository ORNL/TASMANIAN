include ../Config/AltBuildSystems/Makefile.in

IADD = -I../SparseGrids -I../DREAM $(CommonIADD)
LADD = -L../ $(CommonLADD)
LIBS = $(shell pwd)/../libtasmaniansparsegrid.so $(shell pwd)/../libtasmaniandream.so $(CommonLIBS)

LHEADERS = $(wildcard ./*.hpp)

LIBOBJ = tsgCLoadNeededValues.o tsgCConstructSurrogate.o

MODULENAME = libtasmaniancaddons.so

TESTNAME = addontester

all: $(LIBOBJ) ../$(MODULENAME) ../$(TESTNAME)

%.o: %.cpp $(LHEADERS)
	$(CC) $(OPTC) $(IADD) -c $< -o $@

../$(MODULENAME): $(LIBOBJ)
	$(CC) $(OPTL) $(LADD) $(LIBOBJ) -shared -o ../$(MODULENAME) $(LIBS)

../$(TESTNAME): testAddons.o
	$(CC) $(OPTL) $(LADD) testAddons.o -o ../$(TESTNAME) ../libtasmaniandream.a ../libtasmaniansparsegrid.a

clean:
	rm -fr *.o
