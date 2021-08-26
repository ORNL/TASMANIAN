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

# Test Exotic Quadrature.
import TasmanianAddons as TA

depth = 20
dimension = 1
shift = 1
nref = 101
weight_fn = lambda x: 1.0 if x == 0.0 else math.sin(x) / x
description = "Exotic Quadrature"
is_symmetric = False

grid = TasmanianSG.TasmanianSparseGrid()
TA.loadExoticQuadrature(depth, dimension, shift, weight_fn, nref, description, is_symmetric, grid)
pts = grid.getPoints()
print(pts)
wts = grid.getQuadratureWeights()
print(wts)

