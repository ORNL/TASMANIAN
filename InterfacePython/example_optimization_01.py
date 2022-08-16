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

from Tasmanian import Optimization as Opt, DREAM, makeLocalPolynomialGrid
import numpy as np

def example_01():
    print("\n---------------------------------------------------------------------------------------------------\n")
    print("EXAMPLE 1: use the Particle Swarm algorithm to minimize the six-hump camel function")

    # Create an empty Particle Swarm state.
    iNumDimensions = 2
    iNumParticles = 50
    randUnif = DREAM.RandomGenerator(sType="default", iSeed=777)
    state = Opt.ParticleSwarmState(iNumDimensions, iNumParticles)

    # Load the state with uniformly initialized particles in the domain [-3.0, 3.0] x [-2.0, 2.0].
    state.initializeParticlesInsideBox(np.array([-3.0, -2.0]), np.array([3.0, 2.0]), random01=randUnif)

    # Define the batched objective function `func` and a function `inside` that returns True if its input is inside the domain.
    # NOTE: `shc` is the unbatched six-hump camel function.
    shc = lambda x : (4 - 2.1*x[0]*x[0] + x[0]*x[0]*x[0]*x[0]/3)*x[0]*x[0] + x[0]*x[1] + (-4.0 + 4.0*x[1]*x[1])*x[1]*x[1]
    func = lambda x_batch : np.apply_along_axis(shc, 1, x_batch)
    inside = lambda x : bool((-3 <= x[0]) and (x[0] <= 3) and (-2 <= x[1]) and (x[1] <= 2))

    # Run the Particle Swarm (PS) algorithm and check the output.
    iNumIterations = 200;
    Opt.ParticleSwarm(func, inside, 0.5, 2, 2, iNumIterations, state, random01=randUnif)

    aLocal1 = np.array([-0.08984201368301331, +0.7126564032704135])
    aLocal2 = np.array([+0.08984201368301331, -0.7126564032704135])
    aSolution = state.getBestPosition()

    sResult = ""
    fL2Error = 0.0
    for i in range(iNumDimensions):
        sResult = "{0:1s}{1:13.5f}".format(sResult, aSolution[i])
        fL2Error += min((aSolution[i] - aLocal1[i]) ** 2,
                        (aSolution[i] - aLocal2[i]) ** 2)
    fL2Error = np.sqrt(fL2Error)
    print("\nUsing the Particle Swarm algorithm on the EXACT objective function, the computed solution is:")
    print(" computed: {0:1s}".format(sResult))
    print(" L2 error: {0:14e}".format(fL2Error))

    # Create a surrogate model for the six-hump camel function.
    grid = makeLocalPolynomialGrid(2, 1, 10, iOrder=1)
    grid.setDomainTransform(np.array([[-3.0, 3.0], [-2.0, 2.0]])) # set the non-canonical domain
    needed_points = grid.getNeededPoints()
    needed_values = np.resize(func(needed_points), [grid.getNumNeeded(), 1])
    grid.loadNeededValues(needed_values)
    surrogate_func = lambda x_batch : np.squeeze(grid.evaluateBatch(x_batch))

    # Run the PS algorithm on the surrogate and check the output.
    state = Opt.ParticleSwarmState(iNumDimensions, iNumParticles)
    state.initializeParticlesInsideBox(np.array([-3.0, -2.0]), np.array([3.0, 2.0]), random01=randUnif)
    Opt.ParticleSwarm(surrogate_func, inside, 0.5, 2, 2, iNumIterations, state, random01=randUnif)

    aSolution = state.getBestPosition()

    sResult = ""
    fL2Error = 0.0
    for i in range(iNumDimensions):
        sResult = "{0:1s}{1:13.5f}".format(sResult, aSolution[i])
        fL2Error += min((aSolution[i] - aLocal1[i]) ** 2,
                        (aSolution[i] - aLocal2[i]) ** 2)
    fL2Error = np.sqrt(fL2Error)
    print("\nUsing the Particle Swarm algorithm on the SURROGATE objective function, the computed solution is:")
    print(" computed: {0:1s}".format(sResult))
    print(" L2 error: {0:14e}".format(fL2Error))

if __name__ == "__main__":
    example_01()
