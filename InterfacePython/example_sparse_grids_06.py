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

def example_06():

    print("\n---------------------------------------------------------------------------------------------------\n")
    print("Example 6: interpolate f(x,y) = exp(-5 * x^2) * cos(z), using the rleja rule")
    print("           employs adaptive construction\n")

    iNumInputs = 2
    iNumOutputs = 1
    iNumSamplesPerBatch = 1 # the model is set for a single sample per-batch
    def model(aX):
        # note that the model has to return a 2-D numpy.ndarray
        return np.ones((1,1)) * np.exp(-5.0 * aX[0, 0] ** 2.0) * np.cos(aX[0, 1])

    iTestGridSize = 33
    dx = np.linspace(-1.0, 1.0, iTestGridSize) # sample on a uniform grid
    aMeshX, aMeshY = np.meshgrid(dx, dx)
    aTestPoints = np.column_stack([aMeshX.reshape((iTestGridSize**2, 1)),
                                   aMeshY.reshape((iTestGridSize**2, 1))])
    aReferenceValues = np.exp(-5.0 * aTestPoints[:,0]**2) * np.cos(aTestPoints[:,1])
    aReferenceValues = aReferenceValues.reshape((aReferenceValues.shape[0], 1))

    def testGrid(grid, aTestPoints, aReferenceValues):
        aResult = grid.evaluateBatch(aTestPoints)
        for i in range(20):
            aX = aTestPoints[i,:]
        return np.max(np.abs(aResult - aReferenceValues))

    grid = Tasmanian.SparseGrid()
    grid.makeSequenceGrid(iNumInputs, iNumOutputs, 2, "level", "rleja")

    iNumThreads = 1

    print("{0:>8s}{1:>14s}".format("points", "error"))
    for i in range(4):
        iBudget = 50 * (i + 1)

        Tasmanian.constructAnisotropicSurrogate(lambda x, tid : model(x),
                                                iBudget, iNumThreads, iNumSamplesPerBatch, grid,
                                                "iptotal", 0)

        print("{0:>8d}{1:>14s}".format(grid.getNumPoints(),
              "{0:1.4e}".format(testGrid(grid, aTestPoints, aReferenceValues))))


if (__name__ == "__main__"):
    example_06()
