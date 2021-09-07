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

def example_08():

    print("\n---------------------------------------------------------------------------------------------------\n")
    print("Example 8: interpolate different functions demonstrating the different")
    print("           local polynomial rules\n")

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

    def smooth_model(aX):
        return np.ones((1,)) * np.exp(-aX[0]**2) * np.cos(aX[1])

    iOrder = 2
    grid_localp     = Tasmanian.makeLocalPolynomialGrid(iNumInputs, 1, 7, iOrder, "localp")
    grid_semilocalp = Tasmanian.makeLocalPolynomialGrid(iNumInputs, 1, 7, iOrder, "semi-localp")

    Tasmanian.loadNeededValues(lambda x, tid : smooth_model(x), grid_localp, 4)
    Tasmanian.loadNeededValues(lambda x, tid : smooth_model(x), grid_semilocalp, 4)

    print("Using smooth model: f(x, y) = exp(-x*x) * cos(y)")
    print(" rule_localp,     points = {0:1d}   error = {1:1.4e}".format(
        grid_localp.getNumPoints(), get_error(grid_localp, smooth_model, aTestPoints)))
    print(" rule_semilocalp, points = {0:1d}   error = {1:1.4e}".format(
        grid_semilocalp.getNumPoints(), get_error(grid_semilocalp, smooth_model, aTestPoints)))
    print(" If the model is smooth, rule_semilocalp has an advantage.\n")

    def zero_model(aX):
        return np.ones((1,)) * np.cos(0.5 * np.pi * aX[0]) * np.cos(0.5 * np.pi * aX[1])

    grid_localp0 = Tasmanian.makeLocalPolynomialGrid(iNumInputs, 1, 6, iOrder, "localp-zero")

    Tasmanian.reloadLoadedPoints(lambda x, tid : zero_model(x), grid_localp, 4)
    Tasmanian.loadNeededValues(lambda x, tid : zero_model(x), grid_localp0, 4)

    print("Using homogeneous model: f(x, y) = cos(pi * x / 2) * cos(pi * y / 2)")
    print(" rule_localp,  points = {0:1d}   error = {1:1.4e}".format(
        grid_localp.getNumPoints(), get_error(grid_localp, zero_model, aTestPoints)))
    print(" rule_localp0, points = {0:1d}   error = {1:1.4e}".format(
        grid_localp0.getNumPoints(), get_error(grid_localp0, zero_model, aTestPoints)))
    print(" The rule_localp0 uses basis tuned for models with zero boundary.")


if (__name__ == "__main__"):
    example_08()
