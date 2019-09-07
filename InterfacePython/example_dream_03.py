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
from Tasmanian import loadNeededPoints
import numpy

def example_03():
    print("\n---------------------------------------------------------------------------------------------------\n")
    print("EXAMPLE 3: set the inference problem: identify x_0 and x_1 model parameters")
    print("           from data (noise free example)")
    print("           model: f(x) = sin(x_0*M_PI*t + x_1),")
    print("           data: d = sin(5*M_PI*t + 0.3*M_PI) + sin(10*M_PI*t + 0.1*M_PI)")
    print("           compared to Example 2, the data is a superposition of two signals")
    print("           and the posterior is multi-modal")
    print("             -- problem setup --")
    print("           t in [0,1], t is discretized with 32 equidistant nodes")
    print("           the likelihood is exp(- 16 * (f(x) - d)^2)")
    print("           using a sparse grid to interpolate the model\n")

    iNumDimensions = 2
    iNumOutputs = 32
    iNumChains = 500
    iNumBurnupIterations = 1000
    iNumCollectIterations = 300

    def model(aX):
        # using 32 nodes in the interior of (0.0, 1.0)
        t = numpy.linspace(1.0 / 64.0, 1.0 - 1.0 / 64.0, 32)
        return numpy.sin(numpy.pi * aX[0] * t + aX[1])

    aData = model([5.0, 0.3 * numpy.pi]) + model([10.0, 0.1 * numpy.pi])

    grid = makeSequenceGrid(iNumDimensions, iNumOutputs, 30, 'iptotal', 'leja')

    domain_lower = numpy.array([ 1.0, -0.1])
    domain_upper = numpy.array([12.0,  1.7])
    grid.setDomainTransform(numpy.column_stack([domain_lower, domain_upper]))

    loadNeededPoints(lambda x, tid : model(x), grid, 2)

    likely = DREAM.LikelihoodGaussIsotropic(1.0 / float(iNumOutputs), aData)

    state = DREAM.State(iNumChains, grid)
    state.setState(DREAM.genUniformSamples(domain_lower, domain_upper, state.getNumChains()))

    DREAM.Sample(iNumBurnupIterations, iNumCollectIterations,
                 DREAM.Posterior(grid, likely, prior = "uniform", typeForm = DREAM.typeLogform),
                 DREAM.Domain(grid),
                 state,
                 DREAM.IndependentUpdate("gaussian", 0.01),
                 DREAM.DifferentialUpdate(90),
                 typeForm = DREAM.typeLogform)

    aSamples = state.getHistory()

    # split the frequencies into low and high for the two modes of the posterior
    # here 6.5 is the mid-point of the domain
    low_frequency = numpy.average(numpy.array(
        [aSamples[i,0] for i in range(aSamples.shape[0]) if aSamples[i,0] < 6.5]
    ))
    low_correction = numpy.average(numpy.array(
        [aSamples[i,1] for i in range(aSamples.shape[0]) if aSamples[i,0] < 6.5]
    ))
    high_frequency = numpy.average(numpy.array(
        [aSamples[i,0] for i in range(aSamples.shape[0]) if aSamples[i,0] >= 6.5]
    ))
    high_correction = numpy.average(numpy.array(
        [aSamples[i,1] for i in range(aSamples.shape[0]) if aSamples[i,0] >= 6.5]
    ))

    print("Inferred values:")
    print(" low   frequency:{0:13.6f}   error:{1:14.6e}".format(low_frequency,
                                                     numpy.abs(low_frequency - 5.0)))
    print(" low  correction:{0:13.6f}   error:{1:14.6e}\n".format(low_correction,
                                                     numpy.abs(low_correction - 0.3 * numpy.pi)))
    print(" high  frequency:{0:13.6f}   error:{1:14.6e}".format(high_frequency,
                                                     numpy.abs(high_frequency - 10.0)))
    print(" high correction:{0:13.6f}   error:{1:14.6e}".format(high_correction,
                                                    numpy.abs(high_correction - 0.1 * numpy.pi)))

if __name__ == "__main__":
    example_03()
