#!@Tasmanian_string_python_hashbang@

# necessary import for every use of TASMANIAN
#
import Tasmanian
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

x0 = np.array([1.0, 2.0, 3.0], np.float64)
stepsize = 1/4
adaptive_stepsize = 0.001 * stepsize
func = lambda x : x[0] * x[0] / 2 + x[1] * x[1] + 2.0 * x[2] * x[2]
grad = lambda x : np.array([x[0], 2 * x[1], 4.0 * x[2]], np.float64)

# Constant GD
gds = Tasmanian.Optimization.GradientDescentState(x0, adaptive_stepsize)
status = Tasmanian.Optimization.GradientDescent(grad, stepsize, 100, 0.01, gds)
print("Constant GD")
print(status)
print(gds.getX())
print()

# Adaptive GD
gds = Tasmanian.Optimization.GradientDescentState(x0, adaptive_stepsize)
status = Tasmanian.Optimization.AdaptiveGradientDescent(func, grad, 2.0, 2.0, 100, 0.01, gds)
print("Adaptive GD")
print(status)
print(gds.getX())
print()

# Adaptive PGD
proj = lambda x : np.array([min(max(0.5, xi), 2.0) for xi in x], np.float64)
gds = Tasmanian.Optimization.GradientDescentState(x0, adaptive_stepsize)
status = Tasmanian.Optimization.AdaptiveProjectedGradientDescent(func, grad, proj, 2.0, 2.0, 100, 0.01, gds)
print("Adaptive PGD")
print(status)
print(gds.getX())
print()
