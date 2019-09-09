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

def example_01():
    print("\n---------------------------------------------------------------------------------------------------\n")
    print("EXAMPLE 1: make your own probability distribution")
    print("           sample from the Gaussian distribution: f(x) = exp(-x^2)")
    print("           ignoring scaling constants, using 3000 samples")
    print("    See the comments in example_dream_01.cpp\n")

    iNumDimensions = 1
    iNumChains = 30
    iNumBurnupIterations = 200
    iNumCollectIterations = 1000 # total samples are iNumChains x iNumCollectIterations

    state = DREAM.State(iNumChains, GridOrIntDimensions = iNumDimensions)

    # initialize with uniform samples on [-2, 2]
    state.setState(DREAM.genUniformSamples((-2.0,), (2.0,), state.getNumChains()))

    DREAM.Sample(iNumBurnupIterations, iNumCollectIterations,
                 lambda x : numpy.exp( - x[:,0]**2 ), # Gaussian mean 0.0 variance 0.5
                 DREAM.Domain("unbounded"),
                 state,
                 DREAM.IndependentUpdate("uniform", 0.5),
                 DREAM.DifferentialUpdate(90))

    aMean, aVariance = state.getHistoryMeanVariance()

    print("Using regular form:")
    print("mean    {0:13.6f}   error {1:14.6e}".format(aMean[0], numpy.abs(aMean[0])))
    print("variance{0:13.6f}   error {1:14.6e}".format(aVariance[0], numpy.abs(aVariance[0] -0.5)))

    # reset the state
    state = DREAM.State(iNumChains, GridOrIntDimensions = iNumDimensions)

    # initialize with uniform samples on [-2, 2]
    state.setState(DREAM.genUniformSamples((-2.0,), (2.0,), state.getNumChains()))

    DREAM.Sample(iNumBurnupIterations, iNumCollectIterations,
                 lambda x : - x[:,0]**2, # Gaussian mean 0.0 variance 0.5
                 DREAM.Domain("unbounded"),
                 state,
                 DREAM.IndependentUpdate("uniform", 0.5),
                 DREAM.DifferentialUpdate(90),
                 typeForm = DREAM.typeLogform)

    aMean, aVariance = state.getHistoryMeanVariance()

    print("Using logarithm form:")
    print("mean    {0:13.6f}   error {1:14.6e}".format(aMean[0], numpy.abs(aMean[0])))
    print("variance{0:13.6f}   error {1:14.6e}".format(aVariance[0], numpy.abs(aVariance[0] -0.5)))


if __name__ == "__main__":
    example_01()
