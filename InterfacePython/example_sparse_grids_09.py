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

def example_09():

    print("\n---------------------------------------------------------------------------------------------------\n")
    print("Example 9: comparison between local polynomial refinement strategies\n")

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

    def sharp_model(aX):
        return np.ones((1,)) * np.exp(-aX[0]) / (1.0 + 100.0 * np.exp(-10.0 * aX[1]))

    grid_classic = Tasmanian.makeLocalPolynomialGrid(iNumInputs, 1, 2,
                                                     iOrder = -1, sRule = "localp")
    grid_fds = Tasmanian.copyGrid(grid_classic)

    print("Using batch refinement:")
    print("               classic                  fds")
    print("  points         error  points        error")

    fTolerance = 1.E-5;
    while((grid_classic.getNumNeeded() > 0) or (grid_fds.getNumNeeded() > 0)):
        Tasmanian.loadNeededValues(lambda x, tid: sharp_model(x), grid_classic, 4);
        Tasmanian.loadNeededValues(lambda x, tid: sharp_model(x), grid_fds, 4);

        # print the results at this stage
        print("{0:>8d}{1:>14.4e}{2:>8d}{3:>14.4e}".format(
            grid_classic.getNumLoaded(), get_error(grid_classic, sharp_model, aTestPoints),
            grid_fds.getNumLoaded(), get_error(grid_fds, sharp_model, aTestPoints)))

        # setting refinement for each grid
        grid_classic.setSurplusRefinement(fTolerance, 0, "classic");
        grid_fds.setSurplusRefinement(fTolerance, 0, "fds");

    # reset the grid and repeat using construction
    grid_classic = Tasmanian.makeLocalPolynomialGrid(iNumInputs, 1, 2,
                                                     iOrder = -1, sRule = "localp")
    grid_fds.copyGrid(grid_classic)

    print("\nUsing construction:")
    print("               classic                  fds")
    print("  points         error  points        error")

    iNumThreads = 1
    for budget in range(50, 800, 100):
        Tasmanian.constructSurplusSurrogate(
            lambda x, tid : sharp_model(x.reshape((2,))).reshape((1,1)),
            budget, iNumThreads, 1, grid_classic, fTolerance, "classic")
        Tasmanian.constructSurplusSurrogate(
            lambda x, tid : sharp_model(x.reshape((2,))).reshape((1,1)),
            budget, iNumThreads, 1, grid_fds, fTolerance, "fds")

        print("{0:>8d}{1:>14.4e}{2:>8d}{3:>14.4e}".format(
            grid_classic.getNumLoaded(), get_error(grid_classic, sharp_model, aTestPoints),
            grid_fds.getNumLoaded(), get_error(grid_fds, sharp_model, aTestPoints)))


if (__name__ == "__main__"):
    example_09()
