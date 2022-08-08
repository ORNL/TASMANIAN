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

from Tasmanian import Optimization as Opt, makeGlobalGrid
import numpy as np

def example_02():
    print("\n---------------------------------------------------------------------------------------------------\n")
    print("EXAMPLE 2: use various Gradient Descent algorithms to minimize a convex quadratic")

    # Define the objective function `func` and and its gradient `grad`.
    func = lambda x : 2.0 * (x[0] - 1.0) * (x[0] - 1.0) + (x[1] - 2.0) * (x[1] - 2.0) / 2.0
    grad = lambda x : np.array([4.0 * (x[0] - 1.0), x[1] - 2.0])

    # Create an empty Gradient Descent State.
    x0 = np.array([0.0, 0.0])
    state = Opt.GradientDescentState(x0, 0.0)

    # Run the Gradient Descent (GD) algorithm and check the output. Note that the adaptive stepsize in the state is NOT used.
    iNumIterations = 200;
    fTolerance = 1E-3
    dInfo = Opt.GradientDescent(grad, 1/8.0, iNumIterations, fTolerance, state)

    aGlobal = np.array([1.0, 2.0])
    aSolution = state.getX()

    sResult = ""
    fL2Error = 0.0
    for i in range(2):
        sResult = "{0:1s}{1:13.5f}".format(sResult, aSolution[i])
        fL2Error += (aSolution[i] - aGlobal[i]) ** 2
    fL2Error = np.sqrt(fL2Error)
    print("\nUsing the Gradient Descent algorithm on the EXACT objective function, the computed solution is:")
    print(" iterations: {0:14}".format(dInfo['performed_iterations']))
    print("   computed: {0:1s}".format(sResult))
    print("   L2 error: {0:14e}".format(fL2Error))

    # Create a surrogate model for the quadratic.
    # NOTE: local polynomial grids are not well suited for first-order methods like gradient descent.
    grid = makeGlobalGrid(2, 1, 2, "iptotal", "leja")
    needed_points = grid.getNeededPoints()
    needed_values = np.resize(np.apply_along_axis(func, 1, needed_points), [grid.getNumNeeded(), 1])
    grid.loadNeededValues(needed_values)
    surrogate_func = lambda x : grid.evaluate(x)
    surrogate_grad = lambda x : grid.differentiate(x)

    # Run the GD algorithm on the surrogate and check the output.
    state = Opt.GradientDescentState(x0, 0.0)
    dInfo = Opt.GradientDescent(surrogate_grad, 1/8.0, iNumIterations, fTolerance, state)

    aSolution = state.getX()

    sResult = ""
    fL2Error = 0.0
    for i in range(2):
        sResult = "{0:1s}{1:13.5f}".format(sResult, aSolution[i])
        fL2Error += (aSolution[i] - aGlobal[i]) ** 2
    fL2Error = np.sqrt(fL2Error)
    print("\nUsing the Gradient Descent algorithm on the SURROGATE objective function, the computed solution is:")
    print(" iterations: {0:14}".format(dInfo['performed_iterations']))
    print("   computed: {0:1s}".format(sResult))
    print("   L2 error: {0:14e}".format(fL2Error))

    # Run the adaptive GD algorithm and check the output. Note that the adaptive stepsize in the state IS used here.
    state = Opt.GradientDescentState(x0, 1.0)
    dInfo = Opt.AdaptiveGradientDescent(func, grad, 1.25, 1.25, iNumIterations, fTolerance, state)

    aSolution = state.getX()

    sResult = ""
    fL2Error = 0.0
    for i in range(2):
        sResult = "{0:1s}{1:13.5f}".format(sResult, aSolution[i])
        fL2Error += (aSolution[i] - aGlobal[i]) ** 2
    fL2Error = np.sqrt(fL2Error)
    print("\nUsing the ADAPTIVE Gradient Descent algorithm on the EXACT objective function, the computed solution is:")
    print(" iterations: {0:14}".format(dInfo['performed_iterations']))
    print("   computed: {0:1s}".format(sResult))
    print("   L2 error: {0:14e}".format(fL2Error))

    # Run the projected adaptive GD algorithm and check the output. The domain constraint is implicitly enforced by the
    # projection function `proj`.
    state = Opt.GradientDescentState(x0, 1.0)
    proj = lambda x : np.array([min(0.5, max(0.0, x[0])),
                                min(1.5, max(0.0, x[1]))]) # projection onto the box [0.0, 0.5] x [0.0, 1.5]
    dInfo = Opt.AdaptiveProjectedGradientDescent(func, grad, proj, 1.25, 1.25, iNumIterations, fTolerance, state)

    aGlobal = np.array([0.5, 1.5])
    aSolution = state.getX()

    sResult = ""
    fL2Error = 0.0
    for i in range(2):
        sResult = "{0:1s}{1:13.5f}".format(sResult, aSolution[i])
        fL2Error += (aSolution[i] - aGlobal[i]) ** 2
    fL2Error = np.sqrt(fL2Error)
    print("\nUsing the ADAPTIVE PROJECTED Gradient Descent algorithm on the EXACT objective function, the computed solution is:")
    print(" iterations: {0:14}".format(dInfo['performed_iterations']))
    print("   computed: {0:1s}".format(sResult))
    print("   L2 error: {0:14e}".format(fL2Error))


if __name__ == "__main__":
    example_02()
