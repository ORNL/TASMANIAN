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
        if (((sys.version_info.major == 3) and isinstance(callableOrMagnitude, int))
            or ((sys.version_info.major == 2) and isinstance(callableOrMagnitude, (int, long)))):
            self.iPercent = callableOrMagnitude
            self.pCallable = type_dream_dupdate(lambda : 1)
        else:
            self.iPercent = -1
            self.pCallable = type_dream_dupdate(callableOrMagnitude)


class RandomGenerator(object):
    def __init__(self, sType = "default", iSeed = -1, callableRNG = lambda : 1):
        self.sType = sType
        if (sys.version_info.major == 3):
            self.sType = bytes(self.sType, encoding='utf8')
        self.iSeed = iSeed
        self.pCallable = type_dream_random(callableRNG)


def Sample(iNumBurnup, iNumCollect,
           probability_distibution, domain_description, dream_state,
           independent_update = IndependentUpdate("none"),
           differential_update = DifferentialUpdate(100),
           random01 = RandomGenerator(),
           typeForm = typeRegform):
    '''
    Wrapper to TasDREAM::SampleDREAM().
    '''
    pLibDTSG.tsgDreamSample(typeForm, iNumBurnup, iNumCollect,
                            type_dream_pdf(lambda m, n, x, y : pdfWrapper(probability_distibution, m, n, x, y)), dream_state.pStatePntr,
                            domain_description.pGrid, domain_description.pLower, domain_description.pUpper, domain_description.pCallable,
                            c_char_p(independent_update.sType), independent_update.fMagnitude, independent_update.pCallable,
                            differential_update.iPercent, differential_update.pCallable,
                            c_char_p(random01.sType), random01.iSeed, random01.pCallable)



