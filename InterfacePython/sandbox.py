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
import sys, math, matplotlib
matplotlib.rcParams['text.usetex'] = True

import numpy as np
import matplotlib.pyplot as plt

import TasmanianSG as TSG
import TasmanianAddons as TA

# Test Exotic Quadrature.
def getExoticIntegral(depth, f, weight_fn, shift, dimension, nref=101, description="Exotic Quadrature", is_symmetric=False):
    grid = TSG.TasmanianSparseGrid()
    TA.loadExoticQuadrature(depth, dimension, shift, weight_fn, nref, description, is_symmetric, grid)
    n_pts = grid.getNumPoints()
    pts = grid.getPoints()
    quad_wts = grid.getQuadratureWeights()
    f_vals = np.apply_along_axis(f, 1, pts)
    return n_pts, np.sum(f_vals * quad_wts)

def getGaussLegendreIntegral(depth, f, weight_fn, dimension):
    grid = TSG.TasmanianSparseGrid()
    grid.makeGlobalGrid(dimension, 1, depth, "qptotal", "gauss-legendre")
    n_pts = grid.getNumPoints()
    pts = grid.getPoints()
    quad_wts = grid.getQuadratureWeights()
    f_vals = np.apply_along_axis(f, 1, pts)
    weight_fn_vals = np.apply_along_axis(weight_fn, 1, pts)
    return n_pts, np.sum(f_vals * weight_fn_vals * quad_wts)

def testSincInstance(f, freq, shift, dimension, true_integral, max_level = 25):
    num_eq_arr = np.array([])
    err_eq_arr = np.array([])
    num_gl_arr = np.array([])
    err_gl_arr = np.array([])
    sinc = np.vectorize(lambda x: 1.0 if x == 0.0 else math.sin(freq * x) / (freq * x))
    weight_fn = lambda xs: np.prod(sinc(xs))
    # Exotic Quadrature
    for i in range(max_level):
        num_eq, val_eq = getExoticIntegral(i + dimension, f, weight_fn, shift, dimension, is_symmetric=True)
        num_eq_arr = np.append(num_eq_arr, num_eq)
        err_eq_arr = np.append(err_eq_arr, math.log10(abs(val_eq - true_integral)+1e-16))
    # Gauss-Legendre
    num_gl = 0
    j = 0
    while num_gl < num_eq:
        num_gl, val_gl = getGaussLegendreIntegral(j+dimension, f, weight_fn, dimension)
        num_gl_arr = np.append(num_gl_arr, num_gl)
        err_gl_arr = np.append(err_gl_arr, math.log10(abs(val_gl - true_integral)+1e-16))
        j += 1
    plt.plot(num_eq_arr, err_eq_arr)
    plt.plot(num_gl_arr, err_gl_arr)
    plt.xlabel(r'Number of Points')
    plt.ylabel(r'$\log_{10}$ Error')
    plt.title(r'$\log_{10}$ Error with $d=' + str(int(dimension)) + r', \rho(x)={\rm sinc}(' + str(int(freq)) +
              r'x), f(\bf{x})=\exp(-\|\bf{x}\|^2)$')
    plt.legend(['Exotic Quadrature', 'Gauss-Legendre'])
    plt.show()

# Test the integral of exp(-x*x) * sinc(x).
f = lambda xs: math.exp(-xs.dot(xs))
# testSincInstance(f, 1.0, 1.0, 1, 1.4321357541271255)
# testSincInstance(f, 10.0, 1.0, 1, 0.32099682841103033)
# testSincInstance(f, 100.0, 1.0, 1, 0.031353648322695503)
# testSincInstance(f, 10.0, 1.0, 2, 0.32099682841103033 ** 2)
testSincInstance(f, 100.0, 1.0, 2, 0.031353648322695503 ** 2)
# testSincInstance(f, 10.0, 1.0, 3, 0.32099682841103033 ** 3)
# testSincInstance(f, 100.0, 1.0, 3, 0.031353648322695503 ** 3)
# testSincInstance(f, 10.0, 1.0, 4, 0.32099682841103033 ** 4, max_level=15)
# testSincInstance(f, 100.0, 1.0, 4, 0.031353648322695503 ** 4, max_level=15)
