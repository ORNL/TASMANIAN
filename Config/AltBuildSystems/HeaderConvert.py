#!/usr/bin/env python3

# this script takes TasmanianSparseGridWrapC.cpp and parces out
# all the methods, then puts them in TasmanianSparseGrid.h
# methods must be defined on a single line and start with type + tsg
# e.g., void tsg... or void *tsg

with open("../../SparseGrids/TasmanianSparseGrid.hpp") as tsghpp:
    lTsgHpp = tsghpp.readlines()

with open("../../SparseGrids/TasmanianSparseGridWrapC.cpp") as tsgcpp:
    lTsgCpp = tsgcpp.readlines()

lTsgH =[]

iStage = 0
for sLine in lTsgHpp:
    sLine = sLine.replace('\n', '')
    if (iStage == 0): # copyright statement
        lTsgH.append(sLine + "\n")
        if ("*/" in sLine):
            iStage = 1
    elif (iStage == 1): # ifndef statement
        lTsgH.append("\n")
        lTsgH.append("#ifndef __TASMANIAN_SPARSE_GRID_H\n")
        lTsgH.append("#define __TASMANIAN_SPARSE_GRID_H\n")
        lTsgH.append("\n")
        lTsgH.append("// ------------ C Interface for TasmanianSparseGrid -------------- //\n")
        lTsgH.append("// NOTE: you need to explicitly call the constructor and destructor\n")
        lTsgH.append("//       in all cases, void *grid is a pointer to a C++ class\n")
        lTsgH.append("\n")
        iStage = 2

iStage = 0
for sLine in lTsgCpp:
    if (iStage == 0):
        if ('extern "C"' in sLine):
            iStage = 1
    elif (iStage == 1):
        # if this is the definition of a function
        if (sLine.startswith("//")):
            lTsgH.append(sLine)
        elif (("void* tsg" in sLine) or
            ("void *tsg" in sLine) or
            ("void * tsg" in sLine) or
            ("void tsg" in sLine) or
            ("int* tsg" in sLine) or
            ("int *tsg" in sLine) or
            ("int * tsg" in sLine) or
            ("int tsg" in sLine) or
            ("double* tsg" in sLine) or
            ("double *tsg" in sLine) or
            ("double * tsg" in sLine) or
            ("double tsg" in sLine) or
            ("char* tsg" in sLine) or
            ("char *tsg" in sLine) or
            ("char * tsg" in sLine) or
            ("char tsg" in sLine)):
            sDef = sLine.split("{")[0] + ";"

            lTsgH.append(sDef + "\n")
            #print(sDef)

lTsgH.append("\n")
lTsgH.append("#endif\n")

with open("TasmanianSparseGrid.h", 'w') as outp:
    outp.writelines(lTsgH)
