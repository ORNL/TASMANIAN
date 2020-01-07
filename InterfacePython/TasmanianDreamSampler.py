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
from TasmanianConfig import enableVerboseErrors
from TasmanianConfig import TasmanianInputError as InputError

from TasmanianSG import TasmanianSparseGrid as SparseGrid

from TasmanianDreamLikely import *
from TasmanianDreamState import DreamState as State

type_dream_pdf     = CFUNCTYPE(None, c_int, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_int))
type_dream_domain  = CFUNCTYPE(c_int, c_int, POINTER(c_double))
type_dream_iupdate = CFUNCTYPE(None, c_int, POINTER(c_double), POINTER(c_int))
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
                                    c_char_p, c_int, type_dream_random, POINTER(c_int)]


def pdfWrapper(callableProbability, iNumSamples, iNumDims, pX, pY, pErrInfo):
    '''
    Internal use only.

    Wraps around a user-provided callable method callableProbability
    to facilitate callback from the C/C++ code.
    The wrapper corresponds to a single output pdf/likelihood method.

    iNumSamples is the number of samples to process for this call.
    iNumDims is the number of dimensions of each sample vector
    pX is a POINTER(double) with size iNumSamples * iNumDims
    pY is a POINTER(double) with size iNumSamples
    '''
    pErrInfo[0] = 1
    aX = np.ctypeslib.as_array(pX, (iNumSamples, iNumDims))
    aY = np.ctypeslib.as_array(pY, (iNumSamples,))
    aResult = callableProbability(aX)
    if aResult.shape != (iNumSamples,) and enableVerboseErrors:
        print("ERROR: incorrect output from the probability distribution, should be (iNumSamples,)")
        return
    aY[0:iNumSamples] = aResult[0:iNumSamples]
    pErrInfo[0] = 0

def domainWrapper(callableDomain, iNumDims, pX):
    '''
    Internal use only.

    Wraps around a user-provided callable method callableDomain
    to facilitate callback from the C/C++ code.
    The wrapper corresponds to the test of a single point.

    iNumDims is the dimensions of the point
    pX is a POINTER(double) with size iNumDims
    returns 1 if callableDomain() returns true and 0 otherwise.
    '''
    aX = np.ctypeslib.as_array(pX, (iNumDims,))
    if callableDomain(aX):
        return 1;
    else:
        return 0

def iupdateWrapper(callableIUpdate, iNumDims, pX, pErrInfo):
    '''
    Internal use only.

    Wraps around a user-provided callable method callableIUpdate
    to facilitate callback from the C/C++ code.
    The wrapper corresponds to the independent update of a single point.

    iNumDims is the dimensions of the point
    pX is a POINTER(double) with size iNumDims
    pX will be overwritten with pX + update
    '''
    pErrInfo[0] = 1
    aX = np.ctypeslib.as_array(pX, (iNumDims,))
    aResult = callableIUpdate()
    if aResult.shape != (iNumDims,) and enableVerboseErrors:
        print("ERROR: incorrect output from the independent update, should be (iNumDimensions,)")
        return
    aX[0:iNumDims] = aResult[0:iNumDims]
    pErrInfo[0] = 0


class Domain(object):
    '''
    Domain objects specify the domain of sampling and is used to test
    each candidate sample. The objects roughly correspond to the inside()
    parameter given to TasDREAM::SampleDREAM() but with extra flexibility
    to reduce the cost of C++ to Python callbacks.

    The domain can be specified in one of the following ways:
    1) Domain(grid)
       The grid is an instance of Tasmanian.SparseGrid() and the domain
       will match the domain of the sparse grid. The grid domain is usually
       a hypercube with lower and upper limit, but in some cases, e.g.,
       using gauss-hermite and gauss-laguerre rules, the domain can be
       unbounded in some directions.
    2) Domain("hypercube", aLower, aUpper)
       The domain will be set to a hypercube, similar to TasDREAM::hypercube().
       The aLower and aUpper must be one dimensional numpy.ndarrays with
       size equal to the dimensions of the domain and the corresponding
       values must specify the lower and upper limit in each direction.
    3) Domain("unbounded")
       Samples will be accepted and rejected only based on the probability
       density, i.e., the sampling domain includes all real valued vectors.
    4) Domain("custom", inside)
       The inside parameter must be a callable object that accepts an single
       numpy.ndarray with size equal to the number of dimensions, and must
       return True or False depending on whether the vector is inside
       the domain.
    '''
    def __init__(self, sType, *args):
        '''
        Constructor that translates the input arguments into inputs
        to the internal Tasmanian API.

        See help(Tasmanian.DREAM.Domain)
        '''
        self.TasmanianDreamDomain = True
        if hasattr(sType, "TasmanianSparseGridObject"):
            self.pGrid = sType.pGrid
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
            raise InputError("Domain", "DREAM.Domain() was given an unknown domain type")


class IndependentUpdate(object):
    '''
    Each step of the DREAM procedure uses two updates to generate new proposed
    vectors for each chain. The independent update is a zero-mean random vector
    pulled from a distribution that does not depend on the other DREAM chains.
    Thus is the Python equivalent to the independent_update() parameter
    to TasDREAM::SampleDREAM().

    The independent update can be specified in one of the following ways:
    1) IndependentUpdate(sType, fMagnitude)
       The sType is a string and fMagnitude is a float64 corresponding to the
       dist and magnitude of the TasDREAM::SampleDREAM() overload.
       The sType must be either "gaussian" or "uniform".
       Note that the magnitude can be zero.
    2) IndependentUpdate(iupdate)
       The iupdate must be a callable object that accepts a numpy.ndarray
       with length equal to the sampling dimension and returns same vector
       but perturbed with the update, e.g., lambda x : x + update().
    '''
    def __init__(self, sType, callableOrMagnitude = 0.0):
        '''
        Constructor that translates the input arguments into inputs
        to the internal Tasmanian API.

        See help(Tasmanian.DREAM.IndependentUpdate)
        '''
        self.TasmanianDreamIndependentUpdate = True
        if not ((sys.version_info.major == 3 and isinstance(sType, str))
           or (sys.version_info.major == 2 and isinstance(sType, basestring))):
            if not callable(sType):
                raise InputError("sType", "DREAM.IndependentUpdate() the sType must be either a string or a callable object")
            self.sType = "null"
            self.fMagnitude = 0.0
            self.pCallable = type_dream_iupdate(lambda n, x, err : iupdateWrapper(sType, n, x, err))
        else:
            self.sType = sType
            self.fMagnitude = callableOrMagnitude
            self.pCallable = type_dream_iupdate(lambda : 1)
        if (sys.version_info.major == 3):
            self.sType = bytes(self.sType, encoding='utf8')


class DifferentialUpdate(object):
    '''
    Each step of the DREAM procedure uses two updates to generate new proposed
    vectors for each chain. The differential update is defined as a scale
    multiplied by the difference between two of the DREAM chains.
    This class sets the scale, i.e., the magnitude of the differential update.
    Thus is the Python equivalent to the differential_update() parameter
    to TasDREAM::SampleDREAM().

    The differential update can be specified in one of the following ways:
    1) DifferentialUpdate(iMagnitude)
       The iMagnitude is an integer that defines the percent scale of the
       differential update, i.e., the update is float(iMagnitude) / 100.0
       This is equivalent to using TasDREAM::const_percent<iMagnitude>
       in the C++ call to TasDREAM::SampleDREAM().
    2) DifferentialUpdate(dupdate)
       The dupdate must be a callable object that accepts no inputs
       and returns a float64 indicating the magnitude.

    Note: by symmetry, positive and negative scales with the same magnitude
          are statistically equivalent.
    '''
    def __init__(self, callableOrMagnitude):
        '''
        Constructor that translates the input arguments into inputs
        to the internal Tasmanian API.

        See help(Tasmanian.DREAM.DifferentialUpdate)
        '''
        self.TasmanianDreamDifferentialUpdate = True
        if (((sys.version_info.major == 3) and isinstance(callableOrMagnitude, int))
            or ((sys.version_info.major == 2) and isinstance(callableOrMagnitude, (int, long)))):
            self.iPercent = callableOrMagnitude
            self.pCallable = type_dream_dupdate(lambda : 1)
        else:
            if not callable(callableOrMagnitude):
                raise InputError("callableOrMagnitude", "DREAM.DifferentialUpdate() the callableOrMagnitude must be either an integer or a callable object")
            self.iPercent = -1
            self.pCallable = type_dream_dupdate(callableOrMagnitude)


class RandomGenerator(object):
    '''
    A class that helps specify the random number generator to be used on the
    C++ side of the call to TasDREAM::SampleDREAM().

    The random number engine can be defined in one of several ways:
    1) RandomGenerator("default") or just RandomGenerator()
       The C++ code will use the rand() command with seed initialized from
       the system time.
    2) RandomGenerator("minstd_rand")
       The rand() command in C++ is implementation dependent, i.e., the random
       number may be different between compilers even if the same seed is used.
       The std::minstd_rand generator is defined in the standard and the output
       is implementation independent, which makes it a good option for
       testing and debugging (probably not very good for actual applications).
       The random seed used here will still be taken from the system time.
    3) RandomGenerator("minstd_rand", iSeed)
       In addition to the standard std::minstd_rand engine this also specifies
       the seed.
    4) RandomGenerator(callableRNG = random01)
       The random01 parameter must be a callable object that returns random
       numbers distributed uniformly in (0, 1), e.g.,
       from random import uniform
       RandomGenerator(callableRNG = lambda : uniform(0.0, 1.0))
    '''
    def __init__(self, sType = "default", iSeed = -1, callableRNG = lambda : 1):
        '''
        Constructor that translates the input arguments into inputs
        to the internal Tasmanian API.

        See help(Tasmanian.DREAM.RandomGenerator)
        '''
        if sType not in ["default", "minstd_rand"]:
            raise InputError("sType", "DREAM.RandomGenerator() the sType is invalid")
        if not callable(callableRNG):
            raise InputError("callableRNG", "DREAM.RandomGenerator() the callableRNG must be a callable object")
        self.TasmanianDreamRandomGenerator = True
        self.sType = sType
        if (sys.version_info.major == 3):
            self.sType = bytes(self.sType, encoding='utf8')
        self.iSeed = iSeed
        self.pCallable = type_dream_random(callableRNG)

typePriorCallable = 0 # internal API, the values are synchronized with the C/C++ code
typePriorUniform  = 1

class Posterior(object):
    '''
    A class that helps combine the components of a Bayesian posterior into
    a single callable object. Strictly speaking, this is a factory class that
    generates probability distribution callable objects used in DREAM.Sample().
    The intended use of the class is embedded as in the C++ function
    TasDREAM::posterior().

    The probability distribution takes a two dimensional numpy.ndarray
    object with shape (num_samples, num_dimensions) and returns a one
    dimensional numpy.ndarray with length (num_samples,).
    Each entry of the result corresponds to the (unscaled) probability
    density at the corresponding vector.

    The usage of this class is:
    DREAM.Sample(..., Posterior(model, likelihood, prior), ...)
    or if log-form sampling is to be used:
    DREAM.Sample(..., Posterior(model, likelihood, prior, typeForm), ...)

    Combining the three components gives us:
        result = likelihood(typeForm, model(x)) * prior(x)
    of if using typeLogform:
        result = likelihood(typeForm, model(x)) + prior(x)

    The model, likelihood, and prior, can be callable objects so that

    model: input.shape  = (num_samples, num_dimensions)
           output.shape = (num_samples, num_outputs)

    likelihood: input.shape  = (num_samples, num_outputs)
                output.shape = (num_samples,)

    prior: input.shape  = (num_samples, num_dimensions)
           output.shape = (num_samples,)

    In addition:
    model: can be an instance of Tasmanian.SparseGrid()
    likelihood: can be an instance of DREAM.LikelihoodGaussIsotropic()
                of LikelihoodGaussAnisotropic().
    prior: can be the string "uniform"

    The typeForm must be either typeRegform (default) or typeLogform
    and must be synchronized with the form used in DREAM.Sample().
    The likelihood need only implement the form used in DREAM.Sample().
    The form here indicates whether to add or multiply by the log of
    the values of the prior.
    '''
    def __init__(self, model, likelihood, prior, typeForm = typeRegform):
        '''
        Constructor that translates the input arguments into inputs
        to the internal Tasmanian API.

        See help(Tasmanian.DREAM.Posterior)
        '''
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
        '''
        Internal API, applies the model part of the posterior.
        '''
        if self.bUseSparsegrid:
            return self.model.evaluateBatch(aX)
        else:
            return self.model(aX)

    def doLikelihood(self, aModelValues):
        '''
        Internal API, applies the likelihood part of the posterior.
        '''
        if self.bUseTasLikely:
            return self.likelihood.getLikelihood(self.typeForm, aModelValues)
        else:
            return self.likelihood(self.typeForm, aModelValues)

    def doPrior(self, aLikely, aX):
        '''
        Internal API, computes the prior part of the posterior.
        '''
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
        '''
        Internal API, combines the input objects into a single lambda.
        '''
        return lambda x : self.doPrior(self.doLikelihood(self.doModel(x)), x).reshape((x.shape[0],))

def PosteriorEmbeddedLikelihood(model, prior, typeForm = typeRegform):
    '''
    Shortcut to Tasmanian.DREAM.Posterior() when the model and the likelihood
    are combined into a single object (callable or sparse grid).
    The prior and typeForm are the same.
    '''
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


def genGaussianSamples(lfMean, lfDeviation, iNumSamples, random01 = RandomGenerator(sType = "default")):
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
    The Markov Chain Monte Carlo algorithms generate chains of samples
    with distribution that converges to the target probability density.
    The current state of the chains (i.e., the last vector) is updated by
    first generating a set of proposals and then accepting/rejecting randomly
    but with the constraint that the only accepted vectors belong to the
    specified domain. The proposals are constructed from two components,
    the independent component that is iid random with zero mean and
    the differential component that is takes from the scaled difference
    between the state of two randomly selected chains.

    iNumBurnup: is a non-negative integer that indicates the number of
                iterations to discard, i.e., samples that have no statistical
                significance because the chain state has not converged to
                the target probability distribution.
    iNumCollect: is a non-negative integer that indicates the number of
                iterations to store in the state history.

    probability_distibution: is either an instance of DREAM.Posterior() or
                a callable object that accepts a two dimensional numpy.ndarray
                and then returns a one dimensional numpy.ndarray.
                See help(Tasmanian.DREAM.Posterior)
                Can also use a lambda to capture a Tasmanian sparse grid, i.e.,
                lambda x : grid.evaluateBatch(x).reshape((x.shape[0],))

    domain_description: is an instance of DREAM.Domain()
                See help(Tasmanian.DREAM.Domain)

    dream_state: is an instance of DREAM.State() with number of dimensions
                that match the dimensions used by the probability_distibution.
                See help(Tasmanian.DREAM.State)
                On input the state chains have to be initialized, on output
                the history will be populated with the iNumCollect samples.
                (if the history is not empty, the new samples will be appended)

    independent_update: is an instance of DREAM.IndependentUpdate()
                See help(Tasmanian.DREAM.IndependentUpdate)

    differential_update: is an instance of DREAM.DifferentialUpdate()
                See help(Tasmanian.DREAM.DifferentialUpdate)

    random01: is an instance of DREAM.RandomGenerator()
                See help(Tasmanian.DREAM.RandomGenerator)

    typeForm: is either DREAM.typeRegform or DREAM.typeLogform
               The form modifies the accept/reject criteria based on whether
               the probability_distibution returns the actual distribution
               of interest or the logarithm of that (which is sometimes more
               stable and/or easier to approximate with a sparse grid).

    Note: the typeForm must match the return of the probability_distibution.
    '''
    if not hasattr(domain_description, "TasmanianDreamDomain"):
        raise InputError("domain_description", "domain_description must be an instance of DREAM.Domain()")
    if not hasattr(dream_state, "TasmanainDreamState"):
        raise InputError("dream_state", "dream_state must be an instance of DREAM.State()")
    if not hasattr(independent_update, "TasmanianDreamIndependentUpdate"):
        raise InputError("independent_update", "independent_update must be an instance of DREAM.IndependentUpdate()")
    if not hasattr(differential_update, "TasmanianDreamDifferentialUpdate"):
        raise InputError("differential_update", "differential_update must be an instance of DREAM.DifferentialUpdate()")
    if not hasattr(random01, "TasmanianDreamRandomGenerator"):
        raise InputError("random01", "random01 must be an instance of DREAM.RandomGenerator()")
    if typeForm not in [typeRegform, typeLogform]:
        raise InputError("typeForm", "unknown sampling form, must use typeRegform or typeLogform")

    pErrorCode = (c_int * 1)()

    if hasattr(probability_distibution, "TasmanianPosterior"):
        pLibDTSG.tsgDreamSample(typeForm, c_int(iNumBurnup), c_int(iNumCollect),
                                type_dream_pdf(lambda m, n, x, y, err : pdfWrapper(probability_distibution.distribution(), m, n, x, y, err)),
                                dream_state.pStatePntr,
                                domain_description.pGrid, domain_description.pLower, domain_description.pUpper, domain_description.pCallable,
                                c_char_p(independent_update.sType), independent_update.fMagnitude, independent_update.pCallable,
                                differential_update.iPercent, differential_update.pCallable,
                                c_char_p(random01.sType),
                                random01.iSeed,
                                random01.pCallable,
                                pErrorCode)
    elif callable(probability_distibution):
        pLibDTSG.tsgDreamSample(typeForm, c_int(iNumBurnup), c_int(iNumCollect),
                                type_dream_pdf(lambda m, n, x, y, err : pdfWrapper(probability_distibution, m, n, x, y, err)), dream_state.pStatePntr,
                                domain_description.pGrid, domain_description.pLower, domain_description.pUpper, domain_description.pCallable,
                                c_char_p(independent_update.sType), independent_update.fMagnitude, independent_update.pCallable,
                                differential_update.iPercent, differential_update.pCallable,
                                c_char_p(random01.sType),
                                random01.iSeed,
                                random01.pCallable,
                                pErrorCode)
    else:
        raise InputError("probability_distibution", "probability_distibution must be a callable object that takes a 2D numpy.ndarray and returns a 1D ndarray")

    if pErrorCode[0] != 0:
        raise InputError("Sample", "An error occurred during the call to Tasmanian.")

