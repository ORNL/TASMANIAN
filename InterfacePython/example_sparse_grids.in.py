#!@Tasmanian_string_python_hashbang@

##############################################################################################################################################################################
# Copyright (c) 2017, Miroslav Stoyanov
#
# This file is part of
# Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
#    and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
# THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
# COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
# THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
# IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
##############################################################################################################################################################################

# The sys.path.append() command is necessary only if Tasmanian is
# not included in the system PYTHONPATH
# If PYTHONPATH is set (e.g., source TasmanianENVsetup.sh) you can jump
# straight to import Tasmanian
import sys
@Tasmanian_python_example_import@
import Tasmanian

import example_sparse_grids_01
import example_sparse_grids_02
import example_sparse_grids_03
import example_sparse_grids_04
import example_sparse_grids_05
import example_sparse_grids_06
import example_sparse_grids_07
import example_sparse_grids_08
import example_sparse_grids_09
import example_sparse_grids_10

if __name__ == "__main__":
    example_sparse_grids_01.example_01()
    example_sparse_grids_02.example_02()
    example_sparse_grids_03.example_03()
    example_sparse_grids_04.example_04()
    example_sparse_grids_05.example_05()
    example_sparse_grids_06.example_06()
    example_sparse_grids_07.example_07()
    example_sparse_grids_08.example_08()
    example_sparse_grids_09.example_09()
    example_sparse_grids_10.example_10()
    print("")


exit(0)

import numpy as np

# imports specifically needed by the examples
import math
from random import uniform
from datetime import datetime

bFastTests = False
if (len(sys.argv) > 1):
    bFastTests = True

print("Tasmanian version: {0:s}".format(Tasmanian.__version__))
print("Tasmanian license: {0:s}".format(Tasmanian.__license__))

grid = Tasmanian.SparseGrid()

iDim = 2
iLevel = 6

# EXAMPLE 1: integrate: f(x,y) = exp(-x^2) * cos(y) over [-1,1] x [-1,1]
# using classical Smolyak grid with Clenshaw-Curtis points and weights
grid = Tasmanian.SparseGrid()

print("\n-------------------------------------------------------------------------------------------------")
print("Example 1:  integrate f(x,y) = exp(-x^2) * cos(y), using clenshaw-curtis level nodes")

iDim = 2
iLevel = 6

grid.makeGlobalGrid(iDim, 0, iLevel, "level", "clenshaw-curtis")
aPoints = grid.getPoints()
aWeights = grid.getQuadratureWeights()

fIntegral = sum([aWeights[i] * math.exp(-aPoints[i][0]*aPoints[i][0]) * math.cos(aPoints[i][1]) for i in range(aPoints.shape[0])])
fError = math.fabs(fIntegral - 2.513723354063905e+00)

print("      at level: {0:d}".format(iLevel))
print("      the grid has: {0:d}".format(aPoints.shape[0]))
print("      integral: {0:1.16e}".format(fIntegral))
print("         error: {0:1.16e}\n".format(fError))

iLevel = 7

grid.makeGlobalGrid(iDim, 0, iLevel, "level", "clenshaw-curtis")
aPoints = grid.getPoints()
aWeights = grid.getQuadratureWeights()

fIntegral = sum([aWeights[i] * math.exp(-aPoints[i][0]*aPoints[i][0]) * math.cos(aPoints[i][1]) for i in range(aPoints.shape[0])])
fError = math.fabs(fIntegral - 2.513723354063905e+00)

print("      at level: {0:d}".format(iLevel))
print("      the grid has: {0:d}".format(aPoints.shape[0]))
print("      integral: {0:1.16e}".format(fIntegral))
print("         error: {0:1.16e}".format(fError))

#
# EXAMPLE 2: integrate: f(x,y) = exp(-x^2) * cos(y) over (x,y) in [-5,5] x [-2,3]
# using Gauss-Patterson rules chosen to integrate exactly polynomials of
# total degree up to degree specified by prec

print("\n-------------------------------------------------------------------------------------------------")
print("Example 2: integrate f(x,y) = exp(-x^2) * cos(y) over [-5,5] x [-2,3] using  Gauss-Patterson nodes")

iDim = 2
iPrec = 20

grid.makeGlobalGrid(iDim, 0, iPrec, "qptotal", "gauss-patterson")

grid.setDomainTransform(np.array([[-5.0, 5.0], [-2.0, 3.0]]))

aPoints = grid.getPoints()
aWeights = grid.getQuadratureWeights()

fIntegral = sum([aWeights[i] * math.exp(-aPoints[i][0]*aPoints[i][0]) * math.cos(aPoints[i][1]) for i in range(aPoints.shape[0])])
fError = math.fabs(fIntegral - 1.861816427518323e+00)

print("      at precision: {0:d}".format(iPrec))
print("      the grid has: {0:d}".format(aPoints.shape[0]))
print("      integral: {0:1.16e}".format(fIntegral))
print("         error: {0:1.16e}".format(fError))

iPrec = 40

grid.makeGlobalGrid(iDim, 0, iPrec, "qptotal", "gauss-patterson")

grid.setDomainTransform(np.array([[-5.0, 5.0], [-2.0, 3.0]]))

aPoints = grid.getPoints()
aWeights = grid.getQuadratureWeights()

fIntegral = sum([aWeights[i] * math.exp(-aPoints[i][0]*aPoints[i][0]) * math.cos(aPoints[i][1]) for i in range(aPoints.shape[0])])
fError = math.fabs(fIntegral - 1.861816427518323e+00)

print("      at precision: {0:d}".format(iPrec))
print("      the grid has: {0:d}".format(aPoints.shape[0]))
print("      integral: {0:1.16e}".format(fIntegral))
print("         error: {0:1.16e}".format(fError))

#
# EXAMPLE 3: integrate: f(x,y) = exp(-x^2) * cos(y) over (x,y) in [-5,5] x [-2,3]
# using different rules

print("\n-------------------------------------------------------------------------------------------------")
print("Example 3: integrate f(x,y) = exp(-x^2) * cos(y) over [-5,5] x [-2,3] using different rules\n")

iDim = 2
fExact = 1.861816427518323e+00
mRange = np.array([[-5.0, 5.0], [-2.0, 3.0]])

print("           Clenshaw-Curtis      Gauss-Legendre     Gauss-Patterson")
print(" precision    points     error    points     error    points     error")

liPrecisions = [9, 13, 17, 21, 25, 29]
for iPrec in liPrecisions:
    grid1 = Tasmanian.TasmanianSparseGrid()
    grid1.makeGlobalGrid(iDim, 0, iPrec, "qptotal", "clenshaw-curtis")
    grid2 = Tasmanian.TasmanianSparseGrid()
    grid2.makeGlobalGrid(iDim, 0, iPrec, "qptotal", "gauss-legendre")
    grid3 = Tasmanian.TasmanianSparseGrid()
    grid3.makeGlobalGrid(iDim, 0, iPrec, "qptotal", "gauss-patterson")

    grid1.setDomainTransform(mRange)
    grid2.setDomainTransform(mRange)
    grid3.setDomainTransform(mRange)

    aPoints = grid1.getPoints()
    aWeights = grid1.getQuadratureWeights()
    iNumP1 = aPoints.shape[0]
    fIntegral = sum([aWeights[i] * math.exp(-aPoints[i][0]*aPoints[i][0]) * math.cos(aPoints[i][1]) for i in range(iNumP1)])
    fError1 = math.fabs(fIntegral - fExact)

    aPoints = grid2.getPoints()
    aWeights = grid2.getQuadratureWeights()
    iNumP2 = aPoints.shape[0]
    fIntegral = sum([aWeights[i] * math.exp(-aPoints[i][0]*aPoints[i][0]) * math.cos(aPoints[i][1]) for i in range(iNumP2)])
    fError2 = math.fabs(fIntegral - fExact)

    aPoints = grid3.getPoints()
    aWeights = grid3.getQuadratureWeights()
    iNumP3 = aPoints.shape[0]
    fIntegral = sum([aWeights[i] * math.exp(-aPoints[i][0]*aPoints[i][0]) * math.cos(aPoints[i][1]) for i in range(iNumP3)])
    fError3 = math.fabs(fIntegral - fExact)

    print("    {0:2d}       {1:3d}  {2:1.2e}      {3:4d}  {4:1.2e}       {5:3d}  {6:1.2e}".format(iPrec, iNumP1, fError1, iNumP2, fError2, iNumP3, fError3))

#
# EXAMPLE 4: interpolate: f(x,y) = exp(-x^2) * cos(y)
# with a rule that exactly interpolates polynomials of total degree up to
# degree specified by prec

print("\n-------------------------------------------------------------------------------------------------")
print("Example 4: interpolate f(x,y) = exp(-x^2) * cos(y), using clenshaw-curtis iptotal rule\n")

iDim = 2
iOut = 1
iPrec = 10
aX = np.array([0.3, 0.7])
fExact = math.exp(-aX[0]*aX[0]) * math.cos(aX[1])

grid.makeGlobalGrid(iDim, iOut, iPrec, "iptotal", "clenshaw-curtis")

aPoints = grid.getPoints()

aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI][0] = math.exp(- aPoints[iI][0] * aPoints[iI][0]) * math.cos(aPoints[iI][1])

grid.loadNeededPoints(aVals)

aY = grid.evaluate(aX)

print("  using polynomials of total degree up to: {0:1d}".format(iPrec))
print("                 the grid has: {0:1d} points".format(aPoints.shape[0]))
print("         interpolant at (0.3,0.7): {0:1.16e}".format(aY[0]))
print("                    error: {0:1.16e}\n".format(math.fabs(aY[0] - fExact)))

iPrec = 12
grid.makeGlobalGrid(iDim, iOut, iPrec, "iptotal", "clenshaw-curtis")

aPoints = grid.getPoints()

aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI][0] = math.exp(- aPoints[iI][0] * aPoints[iI][0]) * math.cos(aPoints[iI][1])

grid.loadNeededPoints(aVals)

aY = grid.evaluate(aX)

print("  using polynomials of total degree up to: {0:1d}".format(iPrec))
print("                 the grid has: {0:1d} points".format(aPoints.shape[0]))
print("         interpolant at (0.3,0.7): {0:1.16e}".format(aY[0]))
print("                    error: {0:1.16e}".format(math.fabs(aY[0] - fExact)))

if (bFastTests):
    exit(0)

# random points to be used for error checking
aPnts = np.empty([1000, 4])
for iI in range(1000):
    for iJ in range(4):
        aPnts[iI][iJ] = uniform(-1.0, 1.0)

#
# EXAMPLE 5:
# interpolate: f(x1,x2,x3,x4) = exp(-x1^2) * cos(x2) * exp(-x3^2) * cos(x4)
# with Global and Sequence Leja rules

iDim = 4
iOut = 1
iPrec = 15

aTres = np.empty([1000,])
for iI in range(1000):
    aTres[iI] = math.exp(- aPnts[iI][0] ** 2) * math.cos(aPnts[iI][1]) * math.exp(- aPnts[iI][2] ** 2) * math.cos(aPnts[iI][3])

clock_start = datetime.now()
grid.makeGlobalGrid(iDim, iOut, iPrec, "iptotal", "leja")
clock_end = datetime.now()
stage1 = clock_end - clock_start

aPoints = grid.getPoints()

print("\n-------------------------------------------------------------------------------------------------")
print("Example 5: interpolate f(x1,x2,x3,x4) = exp(-x1^2) * cos(x2) * exp(-x3^2) * cos(x4)")
print("       comparign the performance of Global and Sequence grids with leja nodes")
print("       using polynomials of total degree up to: {0:1d}".format(iPrec))
print("                    grids have: {0:1d} points".format(aPoints.shape[0]))
print("       both grids are evaluated at 1000 random points \n")

aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.exp(- aPoints[iI][0] ** 2) * math.cos(aPoints[iI][1]) * math.exp(- aPoints[iI][2] ** 2) * math.cos(aPoints[iI][3])

clock_start = datetime.now()
grid.loadNeededPoints(aVals)
clock_end = datetime.now()
stage2 = clock_end - clock_start

clock_start = datetime.now()
aRes = grid.evaluateBatch(aPnts)
clock_end = datetime.now()
stage3 = clock_end - clock_start

clock_start = datetime.now()
grid.makeSequenceGrid(iDim, iOut, iPrec, "iptotal", "leja")
clock_end = datetime.now()
sstage1 = clock_end - clock_start

clock_start = datetime.now()
grid.loadNeededPoints(aVals)
clock_end = datetime.now()
sstage2 = clock_end - clock_start

clock_start = datetime.now()
aRes = grid.evaluateBatch(aPnts)
clock_end = datetime.now()
sstage3 = clock_end - clock_start

print("          Stage    Global Grid  Sequence Grid")
print("      make grid {0:14d} {1:14d}   microseconds".format(stage1.microseconds, sstage1.microseconds))
print("    load values {0:14d} {1:14d}   microseconds".format(stage2.microseconds, sstage2.microseconds))
print("       evaluate {0:14d} {1:14d}   microseconds".format(stage3.microseconds, sstage3.microseconds))

print("\n NOTE: due to the somewhat random nature of garbage collection in Python,")
print("       the actual timing above may be inacurate")

#
# EXAMPLE 6:
# interpolate: f(x,y) = exp(-x^2) * cos(y)
# using different refinement schemes

aPnts = aPnts[:,:2]
aTres = np.empty([1000,])
for iI in range(1000):
    aTres[iI] = math.exp(- aPnts[iI][0] ** 2) * math.cos(aPnts[iI][1])

iDim = 2
iOut = 1

grid1.makeGlobalGrid(iDim, iOut, 3, "iptotal", "leja")
grid2.makeGlobalGrid(iDim, iOut, 3, "iptotal", "leja")
grid3.makeSequenceGrid(iDim, iOut, 3, "iptotal", "leja")

aPoints = grid1.getPoints()
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.exp(- aPoints[iI][0] ** 2) * math.cos(aPoints[iI][1])

grid1.loadNeededPoints(aVals)
grid2.loadNeededPoints(aVals)
grid3.loadNeededPoints(aVals)

print("\n-------------------------------------------------------------------------------------------------")
print("Example 6: interpolate: f(x,y) = exp(-x^2) * cos(y)")
print("   using leja nodes and different refinement schemes")
print("   the error is estimated as the maximum from 1000 random points\n")

print("          Toral Degree          Curved         Surplus")
print(" iteration    points     error    points     error    points     error")

for iK in range(10):
    grid1.setAnisotropicRefinement("iptotal", 10, 0)
    aPoints = grid1.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(- aPoints[iI][0] ** 2) * math.cos(aPoints[iI][1])
    grid1.loadNeededPoints(aVals)

    aRes = grid1.evaluateBatch(aPnts)
    fError1 = max(np.fabs(aRes[:,0] - aTres))

    grid2.setAnisotropicRefinement("ipcurved", 10, 0)
    aPoints = grid2.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(- aPoints[iI][0] ** 2) * math.cos(aPoints[iI][1])
    grid2.loadNeededPoints(aVals)

    aRes = grid2.evaluateBatch(aPnts)
    fError2 = max(np.fabs(aRes[:,0] - aTres))

    grid3.setSurplusRefinement(1.E-10, 0)
    aPoints = grid3.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(- aPoints[iI][0] ** 2) * math.cos(aPoints[iI][1])
    grid3.loadNeededPoints(aVals)

    aRes = grid3.evaluateBatch(aPnts)
    fError3 = max(np.fabs(aRes[:,0] - aTres))

    print(" {0:9d} {1:9d}  {2:1.2e} {3:9d}  {4:1.2e} {5:9d}  {6:1.2e}".format(iK+1, grid1.getNumPoints(), fError1, grid2.getNumPoints(), fError2, grid3.getNumPoints(), fError3))

#
# EXAMPLE 7:
# interpolate: f(x,y) = exp(-x^2) * cos(y)
# using localp and semilocalp grids

iDim = 2
iOut = 1
iDepth = 7

print("\n-------------------------------------------------------------------------------------------------")
print("Example 7: interpolate: f(x,y) = exp(-x^2) * cos(y)")
print("       using localp and semi-localp rules with depth {0:1d}".format(iDepth))
print("       the error is estimated as the maximum from 1000 random points\n")

grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, 2, "localp")

aPoints = grid.getPoints()
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.exp(- aPoints[iI][0] ** 2) * math.cos(aPoints[iI][1])
grid.loadNeededPoints(aVals)

aRes = grid.evaluateBatch(aPnts)
fError1 = max(np.fabs(aRes[:,0] - aTres))

grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, 2, "semi-localp")

# localp and semi-localp grids have the same points, no need to recompute aVals
grid.loadNeededPoints(aVals)

aRes = grid.evaluateBatch(aPnts)
fError2 = max(np.fabs(aRes[:,0] - aTres))

print("     Number of points: {0:1d}".format(aPoints.shape[0]))
print("     Error for    localp: {0:1.4e}".format(fError1))
print("     Error for   semi-localp: {0:1.4e}".format(fError2))
print(" Note: semi-localp wins this competition because the function is very smooth")

#
# EXAMPLE 8:
# interpolate: f(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y)
# using localp and semilocalp grids

aTres = np.empty([1000,])
for iI in range(1000):
    aTres[iI] = math.cos(0.5 * math.pi * aPnts[iI][0]) * math.cos(0.5 * math.pi * aPnts[iI][1])

iDim = 2
iOut = 1
iDepth = 7

print("\n-------------------------------------------------------------------------------------------------")
print("Example 8: interpolate f(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y)")
print("       using localp and localp-zero rules with depth {0:1d}".format(iDepth))
print("       the error is estimated as the maximum from 1000 random points\n")

grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, 2, "localp")

aPoints = grid.getPoints()
iNumP1 = aPoints.shape[0]
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.cos(0.5 * math.pi * aPoints[iI][0]) * math.cos(0.5 * math.pi * aPoints[iI][1])
grid.loadNeededPoints(aVals)

aRes = grid.evaluateBatch(aPnts)
fError1 = max(np.fabs(aRes[:,0] - aTres))

grid.makeLocalPolynomialGrid(iDim, iOut, iDepth-1, 2, "localp-zero")

aPoints = grid.getPoints()
iNumP2 = aPoints.shape[0]
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.cos(0.5 * math.pi * aPoints[iI][0]) * math.cos(0.5 * math.pi * aPoints[iI][1])
grid.loadNeededPoints(aVals)

aRes = grid.evaluateBatch(aPnts)
fError2 = max(np.fabs(aRes[:,0] - aTres))

print(" For localp    Number of points: {0:1d}   Error: {1:1.16e}".format(iNumP1, fError1))
print(" For localp-zero   Number of points: {0:1d}   Error: {1:1.16e}".format(iNumP2, fError2))
print(" Note: localp-zero wins this competition because the function is zero at the boundary")

#
# EXAMPLE 9:
# interpolate: f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))
# using different refinement schemes

aTres = np.empty([1000,])
for iI in range(1000):
    aTres[iI] = math.exp(-aPnts[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPnts[iI][1]))

iDim = 2
iOut = 1
iDepth = 2
fTol = 1.E-5

grid1.makeLocalPolynomialGrid(iDim, iOut, iDepth, -1, "localp")

aPoints = grid1.getPoints()
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
grid1.loadNeededPoints(aVals)

grid2.makeLocalPolynomialGrid(iDim, iOut, iDepth, -1, "localp")
grid2.loadNeededPoints(aVals)

print("\n-------------------------------------------------------------------------------------------------")
print("Example 9: interpolate f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))")
print("   the error is estimated as the maximum from 1000 random points")
print("   tolerance is set at 1.E-5 and maximal order polynomials are used\n")

print("               Classic         FDS")
print(" precision    points     error    points     error")

for iK in range(7):
    grid1.setSurplusRefinement(fTol, -1, "classic")
    aPoints = grid1.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
    grid1.loadNeededPoints(aVals)

    aRes = grid1.evaluateBatch(aPnts)
    fError1 = max(np.fabs(aRes[:,0] - aTres))

    grid2.setSurplusRefinement(fTol, -1, "fds")
    aPoints = grid2.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
    grid2.loadNeededPoints(aVals)

    aRes = grid2.evaluateBatch(aPnts)
    fError2 = max(np.fabs(aRes[:,0] - aTres))

    print(" {0:9d} {1:9d}  {2:1.2e} {3:9d}  {4:1.2e}".format(iK+1, grid1.getNumPoints(), fError1, grid2.getNumPoints(), fError2))

#
# EXAMPLE 10:
# interpolate: f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))
# using local polynomails and wavelets

grid1.makeLocalPolynomialGrid(iDim, iOut, iDepth=3, iOrder=1, sRule="localp")
aPoints = grid1.getPoints()
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
grid1.loadNeededPoints(aVals)

grid2.makeWaveletGrid(iDim, iOut, iDepth=1, iOrder=1)
aPoints = grid2.getPoints()
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
grid2.loadNeededPoints(aVals)

print("\n-------------------------------------------------------------------------------------------------")
print("Example 10: interpolate f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))")
print("   the error is estimated as the maximum from 1000 random points")
print("   using local polynomials and wavelets\n")

print("           Polynomials        Wavelets")
print(" precision    points     error    points     error")

for iK in range(8):
    grid1.setSurplusRefinement(fTol, -1, "fds")
    aPoints = grid1.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
    grid1.loadNeededPoints(aVals)

    aRes = grid1.evaluateBatch(aPnts)
    fError1 = max(np.fabs(aRes[:,0] - aTres))

    grid2.setSurplusRefinement(fTol, -1, "fds")
    aPoints = grid2.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
    grid2.loadNeededPoints(aVals)

    aRes = grid2.evaluateBatch(aPnts)
    fError2 = max(np.fabs(aRes[:,0] - aTres))

    print(" {0:9d} {1:9d}  {2:1.2e} {3:9d}  {4:1.2e}".format(iK+1, grid1.getNumPoints(), fError1, grid2.getNumPoints(), fError2))


# random points to be used for error checking
aPnts = np.empty([1000, 3])
for iI in range(1000):
    for iJ in range(3):
        aPnts[iI][iJ] = uniform(-1.0, 1.0)


#
# EXAMPLE 11:
#
# interpolate: f(x,y,z) = 1/((1+4x^2)*(1+5y^2)*(1+6z^2))
# using classical and conformal transformation
#
aTres = 1.0 / ((1.0 + 4.0*aPnts[:,0]**2)*(1.0 + 5.0*aPnts[:,1]**2)*(1.0 + 6.0*aPnts[:,2]**2));

iDim = 3;
iPrec = 12;
iOuts = 1;
grid.makeGlobalGrid(iDim, iOuts, iPrec, 'iptotal', 'clenshaw-curtis');
aPoints = grid.getNeededPoints()
aVals = 1.0 / ((1.0 + 4.0*aPoints[:,0]**2)*(1.0 + 5.0*aPoints[:,1]**2)*(1.0 + 6.0*aPoints[:,2]**2));
grid.loadNeededPoints(aVals.reshape([aPoints.shape[0], 1]));
aRes = grid.evaluateBatch(aPnts);
fError1 = max(np.fabs(aRes[:,0] - aTres))
iNump1 = aPoints.shape[0]

grid.makeGlobalGrid(iDim, iOuts, iPrec, 'iptotal', 'clenshaw-curtis');
grid.setConformalTransformASIN(np.array([4, 4, 4]))
aPoints = grid.getNeededPoints()
aVals = 1.0 / ((1.0 + 4.0*aPoints[:,0]**2)*(1.0 + 5.0*aPoints[:,1]**2)*(1.0 + 6.0*aPoints[:,2]**2));
grid.loadNeededPoints(aVals.reshape([aPoints.shape[0], 1]));
aRes = grid.evaluateBatch(aPnts);
fError2 = max(np.fabs(aRes[:,0] - aTres))

grid.makeLocalPolynomialGrid(iDim, iOuts, iPrec-4, 2);
aPoints = grid.getNeededPoints()
aVals = 1.0 / ((1.0 + 4.0*aPoints[:,0]**2)*(1.0 + 5.0*aPoints[:,1]**2)*(1.0 + 6.0*aPoints[:,2]**2));
grid.loadNeededPoints(aVals.reshape([aPoints.shape[0], 1]));
aRes = grid.evaluateBatch(aPnts);
fError3 = max(np.fabs(aRes[:,0] - aTres))
iNump2 = aPoints.shape[0]

grid.makeLocalPolynomialGrid(iDim, iOuts, iPrec-4, 2);
grid.setConformalTransformASIN(np.array([4, 4, 4]))
aPoints = grid.getNeededPoints()
aVals = 1.0 / ((1.0 + 4.0*aPoints[:,0]**2)*(1.0 + 5.0*aPoints[:,1]**2)*(1.0 + 6.0*aPoints[:,2]**2));
grid.loadNeededPoints(aVals.reshape([aPoints.shape[0], 1]));
aRes = grid.evaluateBatch(aPnts);
fError4 = max(np.fabs(aRes[:,0] - aTres))
iNump2 = aPoints.shape[0]


print("\n----------------------------------------------------------------------------")
print("Example 11: interpolate f(x,y,z) = 1/((1+4x^2)*(1+5y^2)*(1+6z^2))")
print("            using conformal transformation")
print("            the error is estimated as the maximum from 1000 random points")
print("");
print("Grid Type    nodes     error regular   error conformal")
print("Global         {0:1d}         {1:1.3e}         {2:1.3e}".format(iNump1, fError1, fError2))
print("Localp        {0:1d}         {1:1.3e}         {2:1.3e}".format(iNump2, fError3, fError4))
print("")
print("Note: conformal maps address specific problems with the region of analyticity of a function")
print("      the map can accelerate or slow down convergence depending on the problem")

print("\n-------------------------------------------------------------------------------------------------\n")
