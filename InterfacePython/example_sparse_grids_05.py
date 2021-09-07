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

import numpy as np
import Tasmanian

def example_05():

    print("\n---------------------------------------------------------------------------------------------------\n")
    print("Example 5: interpolate f(x,y) = exp(-x^2) * cos(y), using leja rule")
    print("           employ adaptive refinement to increase accuracy per samples")

    iNumInputs = 2
    iNumOutputs = 1
    def model(aX):
        return np.ones((1,)) * np.exp(-aX[0] * aX[0]) * np.cos(aX[1])

    iTestGridSize = 33
    dx = np.linspace(-1.0, 1.0, iTestGridSize) # sample on a uniform grid
    aMeshX, aMeshY = np.meshgrid(dx, dx)
    aTestPoints = np.column_stack([aMeshX.reshape((iTestGridSize**2, 1)),
                                   aMeshY.reshape((iTestGridSize**2, 1))])
    aReferenceValues = np.exp(-aTestPoints[:,0]**2) * np.cos(aTestPoints[:,1])
    aReferenceValues = aReferenceValues.reshape((aReferenceValues.shape[0], 1))

    def testGrid(grid, aTestPoints, aReferenceValues):
        aResult = grid.evaluateBatch(aTestPoints)
        for i in range(20):
            aX = aTestPoints[i,:]
        return np.max(np.abs(aResult - aReferenceValues))


    iInitialLevel = 5

    grid_isotropic = Tasmanian.SparseGrid()
    grid_iptotal   = Tasmanian.SparseGrid()
    grid_icurved   = Tasmanian.SparseGrid()
    grid_surplus   = Tasmanian.SparseGrid()
    grid_isotropic.makeGlobalGrid(iNumInputs, iNumOutputs, iInitialLevel, "level", "leja")
    grid_iptotal.copyGrid(grid_isotropic)
    grid_icurved.copyGrid(grid_isotropic)
    grid_surplus.copyGrid(grid_isotropic)

    iNumThreads = 1
    iBudget = 100
    print("{0:>22s}{1:>22s}{2:>22s}{3:>22s}".format("isotropic", "iptotal", "ipcurved", "surplus"))
    print("{0:>8s}{1:>14s}{0:>8s}{1:>14s}{0:>8s}{1:>14s}{0:>8s}{1:>14s}".format("points", "error"))

    bBelowBudget = True
    while(bBelowBudget):
        sInfo = ""
        if (grid_isotropic.getNumLoaded() < iBudget):
            Tasmanian.loadNeededValues(lambda x, tid : model(x), grid_isotropic, iNumThreads)
            sInfo += "{0:>8d}{1:>14s}".format(grid_isotropic.getNumLoaded(),
                     "{0:1.4e}".format(testGrid(grid_isotropic, aTestPoints, aReferenceValues)))

            iLevel = 0
            while(grid_isotropic.getNumNeeded() == 0):
                grid_isotropic.updateGlobalGrid(iLevel, "level")
                iLevel += 1
        else:
            sInfo += "{0:>22s}".format("")

        if (grid_iptotal.getNumLoaded() < iBudget):
            Tasmanian.loadNeededValues(lambda x, tid : model(x), grid_iptotal, iNumThreads)
            sInfo += "{0:>8d}{1:>14s}".format(grid_iptotal.getNumLoaded(),
                     "{0:1.4e}".format(testGrid(grid_iptotal, aTestPoints, aReferenceValues)))

            grid_iptotal.setAnisotropicRefinement("iptotal", 10, 0);
        else:
            sInfo += "{0:>22s}".format("")

        if (grid_icurved.getNumLoaded() < iBudget):
            Tasmanian.loadNeededValues(lambda x, tid : model(x), grid_icurved, iNumThreads)
            sInfo += "{0:>8d}{1:>14s}".format(grid_icurved.getNumLoaded(),
                     "{0:1.4e}".format(testGrid(grid_icurved, aTestPoints, aReferenceValues)))

            grid_icurved.setAnisotropicRefinement("ipcurved", 10, 0);
        else:
            sInfo += "{0:>22s}".format("")

        if (grid_surplus.getNumLoaded() < iBudget):
            Tasmanian.loadNeededValues(lambda x, tid : model(x), grid_surplus, iNumThreads)
            sInfo += "{0:>8d}{1:>14s}".format(grid_surplus.getNumLoaded(),
                     "{0:1.4e}".format(testGrid(grid_surplus, aTestPoints, aReferenceValues)))

            grid_surplus.setSurplusRefinement(1.E-8, 0)
        else:
            sInfo += "{0:>22s}".format("")

        print(sInfo)
        bBelowBudget = (grid_isotropic.getNumLoaded() < iBudget
                       or grid_icurved.getNumLoaded() < iBudget
                       or grid_icurved.getNumLoaded() < iBudget
                       or grid_surplus.getNumLoaded() < iBudget)


if (__name__ == "__main__"):
    example_05()
