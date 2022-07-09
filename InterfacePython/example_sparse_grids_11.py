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
from random import uniform
import Tasmanian

def example_11():

    print("\n---------------------------------------------------------------------------------------------------\n")
    print("Example 11: Using unstructured data\n")

    iNumInputs = 2 # using two inputs for testing

    # test the error on a uniform dense grid with 10K points
    iTestGridSize = 100
    dx = np.linspace(-1.0, 1.0, iTestGridSize) # sample on a uniform grid
    aMeshX, aMeshY = np.meshgrid(dx, dx)
    aTestPoints = np.column_stack([aMeshX.reshape((iTestGridSize**2, 1)),
                                   aMeshY.reshape((iTestGridSize**2, 1))])

    def get_error(grid, model, aTestPoints):
        aGridResult = grid.evaluateBatch(aTestPoints)
        aModelResult = np.empty((aTestPoints.shape[0], 1), np.float64)
        for i in range(aTestPoints.shape[0]):
            aModelResult[i,:] = model(aTestPoints[i,:])
        return np.max(np.abs(aModelResult[:,0] - aGridResult[:,0]))

    def model(aX):
        return np.ones((1,)) * np.exp( - aX[0]**2 - aX[1]**2 )

    grid = Tasmanian.makeGlobalGrid(iNumInputs, 1, 4, "level", "clenshaw-curtis")

    # generate random data
    iNumData = 2000
    inputs = np.zeros((iNumData, iNumInputs))
    outputs = np.zeros((iNumData, 1))

    for i in range(iNumData):
        inputs[i,0] = uniform(-1.0, 1.0)
        inputs[i,1] = uniform(-1.0, 1.0)
        outputs[i,:] = model(inputs[i,:])

    if (not grid.isAccelerationAvailable("cpu-blas") and
        not grid.isAccelerationAvailable("gpu-cuda") and
        not grid.isAccelerationAvailable("gpu-magma")):
        print("Skipping example 11, BLAS, CUDA, or MAGMA acceleration required.")
        return

    Tasmanian.loadUnstructuredDataL2(inputs, outputs, 1.E-4, grid)

    print("Using construction from unstructured (random) data")
    print("    approximatino error = {0:1.6E}".format(get_error(grid, model, aTestPoints)))
    print("  note: compared to the C++ code, this example does not fix the random seed")
    print("        the error will be slightly different every time you run the example")
    print("")


if (__name__ == "__main__"):
    example_11()
