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
import numpy

def example_05():
    print("\n---------------------------------------------------------------------------------------------------\n")
    print("EXAMPLE 5: infer the frequency and magnitude of two signals from noisy data")
    print("           the model has 5 parameters: f(x_1 ... x_5) = sum x_k sin(k * pi * t)")
    print("           data = 2.0 * sin(2 * pi * t) + sin(4 * pi * t) + noise")
    print("           t in [0, 1], t is discretized using 64 equidistant nodes")
    print("           using two different likelihood functions, corresponding to l-2 and l-1 norms")
    print("           looking for the mode of the posterior, i.e., the optimal fit to the data\n")

    iNumDimensions = 5
    iNumChains = 50
    iNumBurnupIterations  = 1000
    iNumCollectIterations = 1000
    fScaleT = 1.0 / 64.0 # scale in T

    def model(aX):
        # the model assumes 2D array aX
        t = numpy.linspace(1.0 / 128.0, 1.0 - 1.0 / 128.0, 64)
        aSinMatrix = numpy.column_stack([numpy.sin(float(k) * numpy.pi * t) for k in range(1, 6)])
        return numpy.matmul(aSinMatrix, aX.transpose()).transpose()

    aNominalWeights = numpy.array([0.0, 2.0, 0.0, 1.0, 0.0])
    aData = model(aNominalWeights).reshape((64,))
    aData = aData + DREAM.genUniformSamples([-fScaleT for i in range(64)],
                                            [fScaleT for i in range(64)], 1).reshape((64,))

    likely = DREAM.LikelihoodGaussIsotropic(fScaleT, aData)

    domain_lower = numpy.array([0.0 for i in range(iNumDimensions)])
    domain_upper = numpy.array([3.0 for i in range(iNumDimensions)])

    state = DREAM.State(iNumChains, iNumDimensions)
    state.setState(DREAM.genUniformSamples(domain_lower, domain_upper, state.getNumChains()))

    DREAM.Sample(iNumBurnupIterations, iNumCollectIterations,
                 DREAM.Posterior(model, likely, prior = "uniform", typeForm = DREAM.typeLogform),
                 DREAM.Domain("hypercube", domain_lower, domain_upper),
                 state,
                 DREAM.IndependentUpdate("uniform", 0.0),
                 DREAM.DifferentialUpdate(100),
                 typeForm = DREAM.typeLogform)

    aSolution = state.getApproximateMode()

    sResult = ""
    sErrors = ""
    for i in range(iNumDimensions):
        sResult = "{0:1s}{1:13.5f}".format(sResult, aSolution[i])
        sErrors = "{0:1s}{1:13.5e}".format(sErrors, numpy.abs(aSolution[i] - aNominalWeights[i]))
    print("Using Gaussian likelihood, the computed solution is:")
    print(" computed: {0:1s}".format(sResult))
    print("    error: {0:1s}".format(sErrors))

    ###############################################################################################
    # repeat the example, but change the likelihood
    ###############################################################################################
    def l1likelihood(aX, aData):
        # computes the likelihood, the 32 comes from the discretization in t
        aY = model(aX)
        aResult = numpy.empty((aY.shape[0],), numpy.float64)
        for i in range(aY.shape[0]):
            aResult[i] = -32.0 * numpy.sum(numpy.abs(aY[i,:] - aData))
        return aResult

    state = DREAM.State(iNumChains, iNumDimensions) # reset the state
    state.setState(DREAM.genUniformSamples(domain_lower, domain_upper, state.getNumChains()))

    DREAM.Sample(iNumBurnupIterations, iNumCollectIterations,
                 DREAM.PosteriorEmbeddedLikelihood(lambda x : l1likelihood(x, aData),
                                                   prior = "uniform",
                                                   typeForm = DREAM.typeLogform),
                 DREAM.Domain("hypercube", domain_lower, domain_upper),
                 state,
                 DREAM.IndependentUpdate("uniform", 0.0),
                 DREAM.DifferentialUpdate(100),
                 typeForm = DREAM.typeLogform)

    aSolution = state.getApproximateMode()

    sResult = ""
    sErrors = ""
    for i in range(iNumDimensions):
        sResult = "{0:1s}{1:13.5f}".format(sResult, aSolution[i])
        sErrors = "{0:1s}{1:13.5e}".format(sErrors, numpy.abs(aSolution[i] - aNominalWeights[i]))
    print("\nUsing l-1 likelihood, the computed solution is:")
    print(" computed: {0:1s}".format(sResult))
    print("    error: {0:1s}".format(sErrors))

if __name__ == "__main__":
    example_05()
