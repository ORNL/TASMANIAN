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
iNumParticles = 10
PSS = OPT.ParticleSwarmState(iNumDimensions, iNumParticles)

print("iNumDimensions = ", PSS.getNumDimensions())
print("iNumParticles = ", PSS.getNumParticles(), "\n")

pp = PSS.getParticlePositions()
print("init_pp = \n", pp, "\n")

pv = PSS.getParticleVelocities()
print("init_pv = \n", pv, "\n")

bpp = PSS.getBestParticlePositions()
print("init_bpp = \n", bpp, "\n")

bp = PSS.getBestPosition()
print("init_bp = \n", bp, "\n")

a1 = np.arange(iNumParticles * iNumDimensions, dtype=np.float64)
a1.resize((iNumParticles, iNumDimensions))
a2 = np.arange(-iNumParticles * iNumDimensions, dtype=np.float64, step=-1)
a2.resize((iNumParticles, iNumDimensions))
a3 = np.arange(start=(iNumParticles + 1) * iNumDimensions, stop=2 * (iNumParticles + 1) * iNumDimensions, dtype=np.float64)
a3.resize((iNumParticles + 1, iNumDimensions))

PSS.setParticlePositions(a1)
pp = PSS.getParticlePositions()
print("set_pp = \n", pp, "\n")

PSS.setParticleVelocities(a2)
pv = PSS.getParticleVelocities()
print("set_pv = \n", pv, "\n")

PSS.setBestParticlePositions(a3)
bpp = PSS.getBestParticlePositions()
print("set_bpp = \n", bpp, "\n")

bp = PSS.getBestPosition()
print("set_bp = \n", bp, "\n")

PSS.clearBestParticles()
bpp = PSS.getBestParticlePositions()
print("clear_bpp = \n", bpp, "\n")

lower = np.asarray([-3.0, -1.0])
upper = np.asarray([3.0, 2.0])
PSS.initializeParticlesInsideBox(lower, upper)

pp = PSS.getParticlePositions()
print("box_pp = \n", pp, "\n")

pv = PSS.getParticleVelocities()
print("box_pv = \n", pv, "\n")

bp = PSS.getBestPosition()
print("box_bp = \n", bp, "\n")

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
    print("f(bpp) = ", fval[0])
