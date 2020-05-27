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

from ctypes import c_char_p, c_int, c_double, c_void_p, POINTER, cdll, create_string_buffer
import numpy as np
import sys

from TasmanianConfig import __path_libdream__
from TasmanianConfig import TasmanianInputError as InputError

pLibDTSG = cdll.LoadLibrary(__path_libdream__)

pLibDTSG.tsgGetNumOutputsLikelihood.restype = c_int

pLibDTSG.tsgDeleteLikelihood.argtypes = [c_void_p]
pLibDTSG.tsgGetNumOutputsLikelihood.argtypes = [c_void_p]
pLibDTSG.tsgGetLikelihood.argtypes = [c_void_p, c_int, POINTER(c_double), c_int, POINTER(c_double)]

pLibDTSG.tsgMakeLikelihoodGaussIsotropic.restype = c_void_p
pLibDTSG.tsgMakeLikelihoodGaussIsotropic.argtypes = [c_int, c_double, POINTER(c_double), c_int]

pLibDTSG.tsgMakeLikelihoodGaussAnisotropic.restype = c_void_p
pLibDTSG.tsgMakeLikelihoodGaussAnisotropic.argtypes = [c_int, POINTER(c_double), POINTER(c_double), c_int]

typeRegform = 0
typeLogform = 1

class TasmanianLikelihood(object):
    '''
    Generic wrapper to TasDREAM::TasmanianLikelihood derived class.
    '''
    def __init__(self):
        '''
        The derived classes will have a member pClassPntr, init to nullptr here.
        '''
        self.TasmanianLikelihood = True
        self.pClassPntr = c_void_p(None)

    def __del__(self):
        '''
        Delete the grid on exit
        '''
        if self.pClassPntr:
            pLibDTSG.tsgDeleteLikelihood(self.pClassPntr)

    def getNumOutputs(self):
        '''
        Returns the number of outputs set in the likelihood.
        '''
        return pLibDTSG.tsgGetNumOutputs(self.pClassPntr)

    def getLikelihood(self, typeForm, llfModel):
        '''
        Returns the likelihood for the model outputs.

        typeForm: indicates to use the regular or log-form
                use DREAM.typeRegform or DREAM.typeLogform
        llfModel: is a two dimensional numpy.ndarray with
                .shape[0] > 0
                .shape[1] == getNumOutputs()

        Returns a one dimensional ndarray with size .shape[0]

        '''
        iNumSamples = llfModel.shape[0]
        iNumOutputs = llfModel.shape[1]
        aResult = np.empty((iNumSamples,), np.float64)
        pLibDTSG.tsgGetLikelihood(self.pClassPntr, typeForm,
                                  np.ctypeslib.as_ctypes(llfModel.reshape((iNumSamples*iNumOutputs,))),
                                  iNumSamples, np.ctypeslib.as_ctypes(aResult))
        return aResult


class LikelihoodGaussIsotropic(TasmanianLikelihood):
    '''
    Wrapper around TasDREAM::LikelihoodGaussIsotropic
    '''
    def __init__(self, fVariance, lfData, iNumSamples = 1):
        '''
        Make a new likelihood with the given data.

        fVariance:   float64 is the variance
        lfData:      one dimensional ndarray with the data
        iNumSamples: number of samples used for the data
        '''
        super(LikelihoodGaussIsotropic, self).__init__()
        iNumOutputs = len(lfData)
        self.pClassPntr = pLibDTSG.tsgMakeLikelihoodGaussIsotropic(
            iNumOutputs, fVariance, np.ctypeslib.as_ctypes(lfData), iNumSamples)


class LikelihoodGaussAnisotropic(TasmanianLikelihood):
    '''
    Wrapper around TasDREAM::LikelihoodGaussIsotropic
    '''
    def __init__(self, lfVariance, lfData, iNumSamples = 1):
        '''
        Make a new likelihood with the given data.

        lfVariance:  one dimensional ndarray with the data variance
        lfData:      one dimensional ndarray with the data
        iNumSamples: number of samples used for the data
        '''
        super(LikelihoodGaussAnisotropic, self).__init__()
        if (len(lfVariance) != len(lfData)):
            raise InputError("lfVariance", "Mismatch in size of variance and data vectors.")
        iNumOutputs = len(lfData)
        self.pClassPntr = pLibDTSG.tsgMakeLikelihoodGaussAnisotropic(
            iNumOutputs, np.ctypeslib.as_ctypes(lfVariance), np.ctypeslib.as_ctypes(lfData), iNumSamples)
