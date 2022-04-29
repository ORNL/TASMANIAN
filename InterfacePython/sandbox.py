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

import TasmanianOPT as OPT

iNumDimensions = 2
iNumParticles = 3
PSS = OPT.ParticleSwarmState(iNumDimensions, iNumParticles)

#print(PSS.getNumDimensions())

pp = PSS.getParticlePositions()
print("init_pp = \n", pp, "\n")

lower = np.asarray([-3.0, -1.0])
upper = np.asarray([3.0, 2.0])
PSS.initializeParticlesInsideBox(lower, upper)

#print(PSS.getParticlePositions())
#print(PSS.getParticleVelocities())

pp = PSS.getParticlePositions()
print("box_pp = \n", pp, "\n")

def l1_fn(x_batch, fval_batch):
    for i in range(fval_batch.shape[0]):
        fval_batch[i] = abs(x_batch[i,0]) + abs(x_batch[i,1])

def inside_fn(x):
    return -3.0 <= x[0] and x[0] <= 3.0 and -1.0 <= x[1] and x[1] <= 2.0

fval = np.empty((1,))
for i in range(20):
    OPT.ParticleSwarm(l1_fn, 1, inside_fn, PSS, 0.5, 2.0, 2.0)
    bp = PSS.getBestPosition()
    l1_fn(np.asarray([bp]), fval)
    print("i = ", i, ", f(bpp) = ", fval[0], "\n")
