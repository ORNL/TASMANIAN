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

ct = TasmanianSG.makeCustomTabulatedFromFile(tdata.sGaussPattersonTableFile)
print(ct.getDescription())
for i in range(ct.getNumLevels()):
    print(ct.getNumPoints(i))
    print()

sub_ct = TasmanianSG.makeCustomTabulatedSubset(ct, 0, 2, "Test")
print(sub_ct.getDescription())
for i in range(sub_ct.getNumLevels()):
    print(sub_ct.getNumPoints(i))
    print()

