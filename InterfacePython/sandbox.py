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
import numpy as np
import matplotlib.pyplot as plt

def getExoticIntegral(depth, f, weight_fn, shift, nref=101, description="Exotic Quadrature", is_symmetric=False):
    grid = TasmanianSG.TasmanianSparseGrid()
    TA.loadExoticQuadrature(depth, 1, shift, weight_fn, nref, description, is_symmetric, grid)
    pts = grid.getPoints()
    wts = grid.getQuadratureWeights()
    return np.sum(f(pts[:, 0]) * wts)

def getGaussLegendreIntegral(depth, f, weight_fn):
    grid = TasmanianSG.makeGlobalGrid(1, 0, depth, "qptotal", "gauss-legendre")
    pts = grid.getPoints()
    wts = grid.getQuadratureWeights()
    return np.sum(f(pts[:, 0]) * weight_fn(pts[:, 0]) * wts)

def testSincInstance(f, freq, shift, true_integral):
    num_eq = np.array([])
    err_eq = np.array([])
    num_gl = np.array([])
    err_gl= np.array([])
    weight_fn = np.vectorize(lambda x: 1.0 if x == 0.0 else math.sin(freq * x) / (freq * x))
    num_problems = 40
    for i in range(num_problems):
        num_eq = np.append(num_eq, i+2)
        val_eq = getExoticIntegral(i+1, f, weight_fn, shift)
        err_eq = np.append(err_eq, math.log10(abs(val_eq - true_integral)+1e-12))
    for j in range(num_problems * 2):
        num_gl = np.append(num_gl, (j+2)/2 if j % 2 == 0 else (j+3)/2)
        val_gl = getGaussLegendreIntegral(j+1, f, weight_fn)
        err_gl = np.append(err_gl, math.log10(abs(val_gl - true_integral)+1e-12))
    plt.plot(num_eq, err_eq)
    plt.plot(num_gl, err_gl)
    plt.xlabel('Number of Points')
    plt.ylabel('Log10 Error')
    plt.title('Log10 Error with ρ(x)=sinc(' + str(int(freq)) + 'x), f(x)=exp(-x^2)')
    plt.legend(['Exotic Quadrature', 'Gauss-Legendre'])
    plt.show()

# Test the integral of exp(-x*x) * sinc(x).
f = np.vectorize(lambda x: math.exp(-x * x))
# testSincInstance(f, 1.0, 1.0, 1.4321357541271255)
# testSincInstance(f, 10.0, 1.0, 0.32099682841103033)
testSincInstance(f, 100.0, 1.0, 0.031353648322695503)

