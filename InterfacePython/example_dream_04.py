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
from Tasmanian import makeSequenceGrid
from Tasmanian import loadNeededValues
import numpy

def example_04():
    print("\n---------------------------------------------------------------------------------------------------\n")
    print("EXAMPLE 4: similar to Example 3, but the data is noisy and the multi-modal posterior comes")
    print("           from a symmetry in the model")
    print("           set inference problem: identify x_0, x_1, x_2 and x_3 (four) model parameters")
    print("                                  from noisy data")
    print("           model: f(x) = x_0*exp(x_1*t) + x_2*exp(x_3*t),  data: d = exp(t) + 0.4*exp(3*t)")
    print("           note the symmetry between x_1 and x_3 which results in a bi-modal posterior")
    print("           t in [0,1], t is discretized with 64 equidistant nodes")
    print("           likelihood is exp(-64 * (f(x) - d)^2)")
    print("           using a sparse grid to interpolate the model\n")

    iNumDimensions = 4
    iNumOutputs = 64
    iNumChains = 50
    iNumBurnupIterations  = 1000
    iNumCollectIterations = 1000

    def model(aX):
        # using 64 nodes in the interior of (0.0, 1.0)
        t = numpy.linspace(1.0 / 128.0, 1.0 - 1.0 / 128.0, 64)
        return aX[0] * numpy.exp(aX[1] * t) + aX[2] * numpy.exp(aX[3] * t)

    aData = (model([1.0, 1.0, 0.4, 3.0]) # add noise
        + DREAM.genGaussianSamples([0.0 for i in range(iNumOutputs)],
                                   [1.0 / float(iNumOutputs) for i in range(iNumOutputs)],
                                   1).reshape((iNumOutputs,)))

    grid = makeSequenceGrid(iNumDimensions, iNumOutputs, 15, 'iptotal', 'leja')

    domain_lower = numpy.array([0.2, 0.5, 0.2, 0.5])
    domain_upper = numpy.array([1.2, 4.0, 1.2, 4.0])
    grid.setDomainTransform(numpy.column_stack([domain_lower, domain_upper]))

    loadNeededValues(lambda x, tid : model(x), grid, 1)

    likely = DREAM.LikelihoodGaussIsotropic(1.0 / float(iNumOutputs), aData)

    state = DREAM.State(iNumChains, grid)
    state.setState(DREAM.genUniformSamples(domain_lower, domain_upper, state.getNumChains()))

    DREAM.Sample(iNumBurnupIterations, iNumCollectIterations,
                 DREAM.Posterior(grid, likely, prior = "uniform", typeForm = DREAM.typeLogform),
                 DREAM.Domain(grid),
                 state,
                 DREAM.IndependentUpdate("gaussian", 0.01),
                 DREAM.DifferentialUpdate(100),
                 typeForm = DREAM.typeLogform)

    aSamples = state.getHistory()

    low_scale = numpy.average(numpy.array(
        [aSamples[i,0] if aSamples[i,1] < aSamples[i,3] else aSamples[i,2]
         for i in range(aSamples.shape[0])]
    ))
    low_rate = numpy.average(numpy.array(
        [aSamples[i,1] if aSamples[i,1] < aSamples[i,3] else aSamples[i,3]
         for i in range(aSamples.shape[0])]
    ))
    high_scale = numpy.average(numpy.array(
        [aSamples[i,0] if aSamples[i,1] >= aSamples[i,3] else aSamples[i,2]
         for i in range(aSamples.shape[0])]
    ))
    high_rate = numpy.average(numpy.array(
        [aSamples[i,1] if aSamples[i,1] >= aSamples[i,3] else aSamples[i,3]
         for i in range(aSamples.shape[0])]
    ))

    print("Acceptance rate: {0:1.3f}".format(state.getAcceptanceRate()))
    print("Inferred values:")
    print(" low   rate:{0:13.6f}   error:{1:14.6e}".format(low_rate,
                                                     numpy.abs(low_rate - 1.0)))
    print(" low  scale:{0:13.6f}   error:{1:14.6e}\n".format(low_scale,
                                                     numpy.abs(low_scale - 1.0)))
    print(" high  rate:{0:13.6f}   error:{1:14.6e}".format(high_rate,
                                                     numpy.abs(high_rate - 3.0)))
    print(" high scale:{0:13.6f}   error:{1:14.6e}".format(high_scale,
                                                    numpy.abs(high_scale - 0.4)))

if __name__ == "__main__":
    example_04()
