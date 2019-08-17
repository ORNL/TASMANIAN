include ../Config/AltBuildSystems/Makefile.in

IADD = -I../SparseGrids -I../DREAM $(CommonIADD)
LADD = -L../ $(CommonLADD)
LIBS = $(shell pwd)/../libtasmaniansparsegrid.so $(shell pwd)/../libtasmaniandream.so $(CommonLIBS)

LHEADERS = $(wildcard ./*.hpp)

LIBOBJ = tsgCLoadNeededValues.o tsgCConstructSurrogate.o

MODULENAME = libtasmaniancaddons.so

all: $(LIBOBJ) $(MODULENAME)

%.o: %.cpp $(LHEADERS)
	$(CC) $(OPTC) $(IADD) -c $< -o $@

$(MODULENAME): $(LIBOBJ)
	$(CC) $(OPTL) $(LADD) $(LIBOBJ) -shared -o $(MODULENAME) $(LIBS)

clean:
	rm -fr *.o
	rm -fr $(MODULENAME)
