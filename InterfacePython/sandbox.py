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

# print("Add code to this file to test, debug, or develop features")
import Tasmanian
num_levels = 3
num_nodes = np.array([2, 3, 4])
precision = np.array([4, 5, 2])
nodes = [np.array([0.0001, 0.7]), np.array([-0.5, 0.5, 0.7]), np.array([0.1, 0.2, 0.3, 0.4])]
weights = [np.array([1.5, 0.5]), np.array([1.0, 0.9, 0.1]), np.array([2.0, 0.0, 0.0, 0.0])]
description = "Exotic Quadrature"
ct = Tasmanian.makeCustomTabulatedFromData(num_levels, num_nodes, precision, nodes, weights, description)

print(ct.getNumLevels())
print(ct.getNumPoints(1))
print(ct.getIExact(1))
print(ct.getQExact(1))

print('Grabbing weights and nodes...')
w, x = ct.getWeightsNodes(1)
print(w)
print(x)
