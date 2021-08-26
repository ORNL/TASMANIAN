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

import TasmanianAddons
from ctypes import CFUNCTYPE, c_double
depth = 10
dimension = 1
shift = 1
nref = 51
weight_fn = lambda x : 1.0 if x == 0.0 else math.sin(x) / x

g = TasmanianAddons.createExoticQuadratureGrid(depth, dimension, shift, weight_fn, nref)
pts = g.getPoints()
print(pts)
wts = g.getQuadratureWeights()
print(wts)
