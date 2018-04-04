#!/usr/bin/env python

# This is a very fragile script, be VERY careful when you run it
# in short, this script takes ../SparseGrids/TasmanianSparseGrid.hpp/cpp files
# parses them and outputs corresponding .h, windows.h, and windows.hpp files
# thus, only the ../SparseGrids/TasmanianSparseGrid.hpp/cpp are ever edited
# while with one go of the script the windows and C headers are generated
# the output is written to this folder, careful not to overwrite things


with open("../SparseGrids/TasmanianSparseGrid.hpp") as tsghpp:
    lTsgHpp = tsghpp.readlines()

with open("../SparseGrids/TasmanianSparseGrid.cpp") as tsgcpp:
    lTsgCpp = tsgcpp.readlines()

lTsgH =[]
lTsgHWin =[]
lTsgHppWin = []

iStage = 0
for sLine in lTsgHpp:
    sLine = sLine.replace('\n', '')
    if (iStage == 0): # copyright statement
        lTsgH.append(sLine + "\n")
        lTsgHWin.append(sLine + "\n")
        lTsgHppWin.append(sLine + "\n")
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

        lTsgHWin.append("\n")
        lTsgHWin.append("#ifndef __TASMANIAN_SPARSE_GRID_H\n")
        lTsgHWin.append("#define __TASMANIAN_SPARSE_GRID_H\n")
        lTsgHWin.append("\n")
        lTsgHWin.append("// ------------ C Interface for TasmanianSparseGrid -------------- //\n")
        lTsgHWin.append("// NOTE: you need to explicitly call the constructor and destructor\n")
        lTsgHWin.append("//       in all cases, void *grid is a pointer to a C++ class\n")
        lTsgHWin.append("\n")
        lTsgHWin.append("#ifdef TSG_DLL\n")
        lTsgHWin.append("#define TSG_API __declspec(dllexport)\n")
        lTsgHWin.append("#else\n")
        lTsgHWin.append("#ifdef TSG_DYNAMIC\n")
        lTsgHWin.append("#define TSG_API __declspec(dllimport)\n")
        lTsgHWin.append("#else\n")
        lTsgHWin.append("#define TSG_API\n")
        lTsgHWin.append("#endif\n")
        lTsgHWin.append("#endif\n")
        lTsgHWin.append("\n")

        lTsgHppWin.append("\n")
        lTsgHppWin.append("#ifndef __TASMANIAN_SPARSE_GRID_HPP\n")
        lTsgHppWin.append("#define __TASMANIAN_SPARSE_GRID_HPP\n")
        lTsgHppWin.append("\n")
        lTsgHppWin.append("#ifdef TSG_DLL\n")
        lTsgHppWin.append("#define TSG_API __declspec(dllexport)\n")
        lTsgHppWin.append("#else\n")
        lTsgHppWin.append("#ifdef TSG_DYNAMIC\n")
        lTsgHppWin.append("#define TSG_API __declspec(dllimport)\n")
        lTsgHppWin.append("#else\n")
        lTsgHppWin.append("#define TSG_API\n")
        lTsgHppWin.append("#endif\n")
        lTsgHppWin.append("#endif\n")

        iStage = 2
    elif (iStage == 2): # skilp until #define
        if ("#define" in sLine):
            iStage = 3
    elif (iStage == 3): # copy headers until class definition
        if ("class TasmanianSparseGrid{" in sLine):
            lTsgHppWin.append("class TSG_API TasmanianSparseGrid{\n") # mod the class definition
            iStage = 4
        else:
            lTsgHppWin.append(sLine + "\n")
    elif (iStage == 4): # copy the rest of the class definition
        lTsgHppWin.append(sLine + "\n")
        if ("};" in sLine):
            iStage = 5
            lTsgHppWin.append("\n")
    else: # nothing left to do with this file
        pass

iStage = 0
for sLine in lTsgCpp:
    if (iStage == 0):
        if ('extern "C"' in sLine):
            iStage = 1
    elif (iStage == 1):
        # if this is the definition of a function
        if (sLine.startswith("//")):
            lTsgH.append(sLine)
            lTsgHWin.append(sLine)
            lTsgHppWin.append(sLine)
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
            lTsgHppWin.append('extern "C" TSG_API ' + sDef + '\n')
            lTsgHWin.append('TSG_API ' + sDef + '\n')
            #print(sDef)

lTsgHppWin.append("\n")
lTsgHppWin.append("}\n")
lTsgHppWin.append("\n")
lTsgHppWin.append("#endif\n")

lTsgH.append("\n")
lTsgH.append("#endif\n")

lTsgHWin.append("\n")
lTsgHWin.append("#endif\n")

with open("TasmanianSparseGrid.windows.hpp", 'w') as outp:
    outp.writelines(lTsgHppWin)

with open("TasmanianSparseGrid.h", 'w') as outp:
    outp.writelines(lTsgH)

with open("TasmanianSparseGrid.windows.h", 'w') as outp:
    outp.writelines(lTsgHWin)
