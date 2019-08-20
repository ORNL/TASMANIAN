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

from ctypes import c_char_p, c_int, c_double, c_void_p, POINTER, cdll, CFUNCTYPE
import numpy as np
import sys

from TasmanianConfig import __path_libdream__
from TasmanianConfig import TasmanianInputError as InputError

from TasmanianSG import TasmanianSparseGrid as SparseGrid

from TasmanianDreamLikely import *
from TasmanianDreamState import DreamState as State

type_dream_pdf     = CFUNCTYPE(None, c_int, c_int, POINTER(c_double), POINTER(c_double))
type_dream_domain  = CFUNCTYPE(c_int, c_int, POINTER(c_double))
type_dream_iupdate = CFUNCTYPE(None, c_int, POINTER(c_double))
type_dream_dupdate = CFUNCTYPE(c_double)
type_dream_random  = CFUNCTYPE(c_double)

pLibDTSG = cdll.LoadLibrary(__path_libdream__)

pLibDTSG.tsgGenUniformSamples.argtypes  = [c_int, c_int, POINTER(c_double), POINTER(c_double), c_char_p, c_int, type_dream_random]
pLibDTSG.tsgGenGaussianSamples.argtypes = [c_int, c_int, POINTER(c_double), POINTER(c_double), c_char_p, c_int, type_dream_random]

pLibDTSG.tsgDreamSample.argtypes = [c_int, c_int, c_int,
                                    type_dream_pdf, c_void_p,
                                    c_void_p, POINTER(c_double), POINTER(c_double), type_dream_domain,
                                    c_char_p, c_double, type_dream_iupdate,
                                    c_int, type_dream_dupdate,
                                    c_char_p, c_int, type_dream_random]


def pdfWrapper(callableProbability, iNumSamples, iNumDims, pX, pY):
    aX = np.ctypeslib.as_array(pX, (iNumSamples, iNumDims))
    aY = np.ctypeslib.as_array(pY, (iNumSamples,))
    aResult = callableProbability(aX)
    aY[0:iNumSamples] = aResult[0:iNumSamples]

def domainWrapper(callableDomain, iNumDims, pX):
    aX = np.ctypeslib.as_array(pX, (iNumDims,))
    if callableDomain(aX):
        return 1;
    else:
        return 0

def iupdateWrapper(callableIUpdate, iNumDims, pX):
    aX = np.ctypeslib.as_array(pX, (iNumDims,))
    aResult = callableIUpdate(aX)
    aX[0:iNumDims] = aResult[0:iNumDims]


class Domain(object):
    def __init__(self, sType, *args):
        self.TasmanianDreamDomain = True
        if isinstance(sType, SparseGrid):
            self.pGrid = args[0].pGrid
            self.pLower = None
            self.pUpper = None
            self.pCallable = type_dream_domain(lambda : 1)
        elif sType == "hypercube":
            self.pGrid = None
            self.pLower = np.ctypeslib.as_ctypes(args[0])
            self.pUpper = np.ctypeslib.as_ctypes(args[1])
            self.pCallable = type_dream_domain(lambda : 1)
        elif sType == "unbounded":
            self.pGrid = None
            self.pLower = None
            self.pUpper = None
            self.pCallable = type_dream_domain(lambda n, x : domainWrapper(lambda x: True, n, x))
        elif sType == "custom":
            self.pGrid = None
            self.pLower = None
            self.pUpper = None
            self.pCallable = type_dream_domain(lambda n, x : domainWrapper(args[0], n, x))
        else:
            raise InputError("Domain", "unknown domain type")


class IndependentUpdate(object):
    def __init__(self, sType, callableOrMagnitude = 0.0):
        self.TasmanianDreamIndependentUpdate = True
        if not ((sys.version_info.major == 3 and isinstance(sType, str))
           or (sys.version_info.major == 2 and isinstance(sType, basestring))):
            self.sType = "null"
            self.pCallable = type_dream_iupdate(lambda n, x, : iupdateWrapper(sType, n, x))
        else:
            self.sType = sType
            self.fMagnitude = callableOrMagnitude
            self.pCallable = type_dream_iupdate(lambda : 1)
        if (sys.version_info.major == 3):
            self.sType = bytes(self.sType, encoding='utf8')


class DifferentialUpdate(object):
    def __init__(self, callableOrMagnitude):
        self.TasmanianDreamDifferentialUpdate = True
        if (((sys.version_info.major == 3) and isinstance(callableOrMagnitude, int))
            or ((sys.version_info.major == 2) and isinstance(callableOrMagnitude, (int, long)))):
            self.iPercent = callableOrMagnitude
            self.pCallable = type_dream_dupdate(lambda : 1)
        else:
            self.iPercent = -1
            self.pCallable = type_dream_dupdate(callableOrMagnitude)


class RandomGenerator(object):
    def __init__(self, sType = "default", iSeed = -1, callableRNG = lambda : 1):
        self.TasmanianDreamRandomGenerator = True
        self.sType = sType
        if (sys.version_info.major == 3):
            self.sType = bytes(self.sType, encoding='utf8')
        self.iSeed = iSeed
        self.pCallable = type_dream_random(callableRNG)

typePriorCallable = 0
typePriorUniform  = 1

class Posterior(object):
    def __init__(self, model, likelihood, prior, typeForm = typeRegform):
        self.TasmanianPosterior = True

        self.typeForm = typeForm
        if self.typeForm not in [typeRegform, typeLogform]:
            raise InputError("typeForm", "unknown sampling form, must use typeRegform or typeLogform")

        self.model = model
        if hasattr(self.model, "TasmanianSparseGridObject"):
            self.bUseSparsegrid = True
        elif callable(self.model):
            self.bUseSparsegrid = False
        else:
            raise InputError("model", "model should be either an instance of Tasmanian.SparseGrid() or a callable object.")

        self.likelihood = likelihood
        if hasattr(self.likelihood, "TasmanianLikelihood"):
            self.bUseTasLikely = True
        elif callable(self.likelihood):
            self.bUseTasLikely = False
        else:
            raise InputError("likelihood", "likelihood should be either derived from TasmanianLikelihood() or a callable object.")

        if prior == "uniform":
            self.typePrior = typePriorUniform
        elif callable(prior):
            self.typePrior = typePriorCallable
            self.prior = prior
        else:
            raise InputError("prior", "prior should be either 'uniform' or a callable object.")

    def doModel(self, aX):
        if self.bUseSparsegrid:
            return self.model.evaluateBatch(aX)
        else:
            return self.model(aX)

    def doLikelihood(self, aModelValues):
        if self.bUseTasLikely:
            return self.likelihood.getLikelihood(self.typeForm, aModelValues)
        else:
            return self.likelihood(self.typeForm, aModelValues)

    def doPrior(self, aLikely, aX):
        if self.typePrior == typePriorUniform:
            aResult = aLikely
        else:
            # reg-form multiplies, log-form adds
            if self.typeForm == typeRegform:
                aResult = aLikely * self.prior(aX)
            else:
                aResult = aLikely + self.prior(aX)
        return aResult

    def distribution(self):
        return lambda x : self.doPrior(self.doLikelihood(self.doModel(x)), x).reshape((x.shape[0],))

def PosteriorEmbeddedLikelihood(model, prior, typeForm = typeRegform):
    return Posterior(model, lambda t, x : x, prior, typeForm)

def genUniformSamples(lfLower, lfUpper, iNumSamples, random01 = RandomGenerator(sType = "default")):
    '''
    Wrapper around TasDREAM::genUniformSamples()

    Using the one dimensional numpy.ndarrays lfLower and lfUpper,
    which describe the upper and lower bounds of a hypercube,
    generate iNumSamples from a uniform distribution using
    the random01 engine.
    The lengths of lfLower and lfUpper must match.

    Returns a two dimensional numpy.ndarray of shape = (iNumSamples, len(lfLower)),
    '''
    if len(lfLower) != len(lfUpper):
        raise InputError("lfUpper", "The length of lfLower and lfUpper do not match.")

    aLower = np.array(lfLower)
    aUpper = np.array(lfUpper)
    iNumDims = len(lfLower)

    aResult = np.empty((iNumDims * iNumSamples,), np.float64)

    pLibDTSG.tsgGenUniformSamples(iNumDims, iNumSamples, np.ctypeslib.as_ctypes(aLower), np.ctypeslib.as_ctypes(aUpper),
                                  c_char_p(random01.sType), random01.iSeed, random01.pCallable,
                                  np.ctypeslib.as_ctypes(aResult))
    return aResult.reshape((iNumSamples, iNumDims))


def tsgGenGaussianSamples(lfMean, lfDeviation, iNumSamples, random01 = RandomGenerator(sType = "default")):
    '''
    Wrapper around TasDREAM::tsgGenGaussianSamples()

    Using the one dimensional numpy.ndarrays lfMean and lfVariance,
    which describe the mean and standard deviation of a Gaussian distribution,
    generate iNumSamples from a normal (Gaussian) distribution using
    the random01 engine.
    The lengths of lfMean and lfDeviation must match.

    Returns a two dimensional numpy.ndarray of shape = (iNumSamples, len(lfLower)),
    '''
    if len(lfMean) != len(lfDeviation):
        raise InputError("lfDeviation", "The length of lfMean and lfDeviation do not match.")

    aMean = np.array(lfMean)
    aDeviation = np.array(lfDeviation)
    iNumDims = len(lfMean)

    aResult = np.empty((iNumDims * iNumSamples,), np.float64)

    pLibDTSG.tsgGenGaussianSamples(iNumDims, iNumSamples, np.ctypeslib.as_ctypes(aMean), np.ctypeslib.as_ctypes(aDeviation),
                                   c_char_p(random01.sType), random01.iSeed, random01.pCallable,
                                   np.ctypeslib.as_ctypes(aResult))
    return aResult.reshape((iNumSamples, iNumDims))


def Sample(iNumBurnup, iNumCollect,
           probability_distibution, domain_description, dream_state,
           independent_update = IndependentUpdate("none"),
           differential_update = DifferentialUpdate(100),
           random01 = RandomGenerator(sType = "default"),
           typeForm = typeRegform):
    '''
    Wrapper to TasDREAM::SampleDREAM().
    '''
    if not hasattr(domain_description, "TasmanianDreamDomain"):
        raise InputError("domain_description", "domain_description must be an instance of DREAM.Domain()")
    if not hasattr(dream_state, "TasmanainDreamState"):
        raise InputError("dream_state", "dream_state must be an instance of DREAM.State()")
    if not hasattr(independent_update, "TasmanianDreamIndependentUpdate"):
        raise InputError("independent_update", "domain_description must be an instance of DREAM.IndependentUpdate()")
    if not hasattr(differential_update, "TasmanianDreamDifferentialUpdate"):
        raise InputError("differential_update", "domain_description must be an instance of DREAM.DifferentialUpdate()")
    if not hasattr(random01, "TasmanianDreamRandomGenerator"):
        raise InputError("random01", "domain_description must be an instance of DREAM.RandomGenerator()")
    if typeForm not in [typeRegform, typeLogform]:
        raise InputError("typeForm", "unknown sampling form, must use typeRegform or typeLogform")

    if hasattr(probability_distibution, "TasmanianPosterior"):
        pLibDTSG.tsgDreamSample(typeForm, c_int(iNumBurnup), c_int(iNumCollect),
                                type_dream_pdf(lambda m, n, x, y : pdfWrapper(probability_distibution.distribution(), m, n, x, y)), dream_state.pStatePntr,
                                domain_description.pGrid, domain_description.pLower, domain_description.pUpper, domain_description.pCallable,
                                c_char_p(independent_update.sType), independent_update.fMagnitude, independent_update.pCallable,
                                differential_update.iPercent, differential_update.pCallable,
                                c_char_p(random01.sType),
                                random01.iSeed,
                                random01.pCallable)
    elif callable(probability_distibution):
        pLibDTSG.tsgDreamSample(typeForm, c_int(iNumBurnup), c_int(iNumCollect),
                                type_dream_pdf(lambda m, n, x, y : pdfWrapper(probability_distibution, m, n, x, y)), dream_state.pStatePntr,
                                domain_description.pGrid, domain_description.pLower, domain_description.pUpper, domain_description.pCallable,
                                c_char_p(independent_update.sType), independent_update.fMagnitude, independent_update.pCallable,
                                differential_update.iPercent, differential_update.pCallable,
                                c_char_p(random01.sType),
                                random01.iSeed,
                                random01.pCallable)
    else:
        raise InputError("probability_distibution", "probability_distibution must a callable object that takes a 2D numpy.ndarray and returns a 1D ndarray")



