#!@Tasmanian_string_python_hashbang@

# necessary import for every use of TASMANIAN
#
import TasmanianSG
import numpy as np
import math

#from random import uniform

#import matplotlib.pyplot as plt
#import matplotlib.colors as cols
#from mpl_toolkits.mplot3d import Axes3D

###############################################################################
# This file is included for testing, debugging and development reasons
# This file will reside in the cmake build folder, but not the install folder
# Basically, add testing code to this file and it will be automatically
# included as part of the cmake
# You can also use the file without cmake
###############################################################################

print("Add code to this file to test, debug, or develop features")

import testConfigureData as tdata # needed for Gauss-Patterson table file
ct = TasmanianSG.CustomTabulated()
ct.read(tdata.sGaussPattersonTableFile)
grid = TasmanianSG.TasmanianSparseGrid()
grid.makeGlobalGridCustom(2, 1, 1, "level", ct)

print(grid.getNumDimensions())
print(grid.getNumOutputs())
print(grid.getRule())
print(grid.getCustomRuleDescription())
print(grid.getNumPoints())
print(grid.getPoints())
print(grid.getQuadratureWeights())
