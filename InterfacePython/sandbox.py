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

pss = Tasmanian.Optimization.ParticleSwarmState(2, 10)
print([pss.isPositionInitialized(),
       pss.isVelocityInitialized(),
       pss.isBestPositionInitialized(),
       pss.isCacheInitialized()])

pss.initializeParticlesInsideBox(np.array([-1.0, -2.0]),
                                 np.array([1.0, 2.0]))

print([pss.isPositionInitialized(),
       pss.isVelocityInitialized(),
       pss.isBestPositionInitialized(),
       pss.isCacheInitialized()])

print(pss.getParticlePositions())

def obj_fn(x_batch, y):
    rows, _ = x_batch.shape
    for i in range(rows):
        y[i] = sum(abs(x_batch[i,:]))

def inside_fn(x):
    return True

Tasmanian.Optimization.ParticleSwarm(obj_fn, 1000, inside_fn, pss, 0.5, 2.0, 2.0)

print(pss.getParticlePositions())

print([pss.isPositionInitialized(),
       pss.isVelocityInitialized(),
       pss.isBestPositionInitialized(),
       pss.isCacheInitialized()])

