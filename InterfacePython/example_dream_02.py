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

from Tasmanian import DREAM
from Tasmanian import SparseGrid
from Tasmanian import makeGlobalGrid
from Tasmanian import makeSequenceGrid
from Tasmanian import loadNeededPoints
import numpy

def example_02():
    print("\n---------------------------------------------------------------------------------------------------\n")
    print("EXAMPLE 2: set the inference problem: identify model parameters x_0 and x_1")
    print("           from data (noise free example)")
    print("           model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(M_PI*t + 0.3*M_PI)")
    print("           t in [0,1], t is discretized with 32 equidistant nodes")
    print("           the likelihood is exp(- 16 * (f(x) - d)^2)")
    print("           using a sparse grid to interpolate the likelihood")
    print("     NOTE: See the comments in the C++ DREAM Example 2\n")

    iNumDimensions = 2;
    iNumChains = 100;
    iNumBurnupIterations = 3000;
    iNumCollectIterations = 100;

    def model(aX):
        # using 32 nodes in the interior of (0.0, 1.0)
        t = numpy.linspace(1.0 / 64.0, 1.0 - 1.0 / 64.0, 32)
        return numpy.sin(numpy.pi * aX[0] * t + aX[1])

    def likelihood(aY, aData):
        # using log-form likelihood
        # numpy.ones((1,)) converts the scalar to a 1D numpy array
        return - 0.5 * len(aData) * numpy.sum((aY - aData)**2) * numpy.ones((1,))

    aData = model([1.0, 0.3 * numpy.pi])

    domain_lower = numpy.array([0.5, -0.1])
    domain_upper = numpy.array([8.0,  1.7])

    grid = makeGlobalGrid(iNumDimensions, 1, 10, 'iptotal', 'clenshaw-curtis')
    grid.setDomainTransform(numpy.column_stack([domain_lower, domain_upper]))

    loadNeededPoints(lambda x, tid : likelihood(model(x), aData), grid, iNumThreads = 2)

    state = DREAM.State(iNumChains, grid)
    state.setState(DREAM.genUniformSamples(domain_lower, domain_upper, state.getNumChains()))

    DREAM.Sample(iNumBurnupIterations, iNumCollectIterations,
                 lambda x : grid.evaluateBatch(x).reshape((x.shape[0],)), # link to the SparseGrid
                 DREAM.Domain(grid),
                 state,
                 DREAM.IndependentUpdate("uniform", 0.5),
                 DREAM.DifferentialUpdate(90),
                 typeForm = DREAM.typeLogform)

    aMean, aVariance = state.getHistoryMeanVariance()

    print("Inferred values (using 10th order polynomial sparse grid):")
    print("    frequency:{0:13.6f}   error:{1:14.6e}".format(aMean[0],
                                                             numpy.abs(aMean[0] - 1.0)))
    print("   correction:{0:13.6f}   error:{1:14.6e}\n".format(aMean[1],
                                                             numpy.abs(aMean[0] - 0.3 * numpy.pi)))

    # repeat the example using the faster Sequence grid and 30th order polynomials
    grid = makeSequenceGrid(iNumDimensions, 1, 30, 'iptotal', 'leja')
    grid.setDomainTransform(numpy.column_stack([domain_lower, domain_upper]))

    loadNeededPoints(lambda x, tid : likelihood(model(x), aData), grid, iNumThreads = 2)

    state = DREAM.State(iNumChains, grid)
    state.setState(DREAM.genUniformSamples(domain_lower, domain_upper, state.getNumChains()))

    DREAM.Sample(iNumBurnupIterations, iNumCollectIterations,
                 lambda x : grid.evaluateBatch(x).reshape((x.shape[0],)), # link to the SparseGrid
                 DREAM.Domain(grid),
                 state,
                 DREAM.IndependentUpdate("uniform", 0.5),
                 DREAM.DifferentialUpdate(90),
                 typeForm = DREAM.typeLogform)

    aMean, aVariance = state.getHistoryMeanVariance()

    print("Inferred values (using 30th order polynomial sparse grid):")
    print("    frequency:{0:13.6f}   error:{1:14.6e}".format(aMean[0],
                                                             numpy.abs(aMean[0] - 1.0)))
    print("   correction:{0:13.6f}   error:{1:14.6e}".format(aMean[1],
                                                             numpy.abs(aMean[0] - 0.3 * numpy.pi)))

if __name__ == "__main__":
    example_02()
