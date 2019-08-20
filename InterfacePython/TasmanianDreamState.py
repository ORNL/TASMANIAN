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

from TasmanianSG import TasmanianSparseGrid as SparseGrid

pLibDTSG = cdll.LoadLibrary(__path_libdream__)

pLibDTSG.tsgMakeDreamState.restype = c_void_p
pLibDTSG.tsgMakeDreamState.argtypes = [c_int, c_int]

pLibDTSG.tsgDeleteDreamState.argtypes = [c_void_p]

pLibDTSG.tsgDreamStateGetDims.restype = c_int
pLibDTSG.tsgDreamStateGetDims.argtypes = [c_void_p]

pLibDTSG.tsgDreamStateGetChains.restype = c_int
pLibDTSG.tsgDreamStateGetChains.argtypes = [c_void_p]

pLibDTSG.tsgDreamStateGetNumHistory.restype = c_int
pLibDTSG.tsgDreamStateGetNumHistory.argtypes = [c_void_p]

pLibDTSG.tsgDreamStateSet.argtypes = [c_void_p, POINTER(c_double)]
pLibDTSG.tsgDreamStateGetHistory.argtypes = [c_void_p, POINTER(c_double)]
pLibDTSG.tsgDreamStateGetHistoryPDF.argtypes = [c_void_p, POINTER(c_double)]
pLibDTSG.tsgDreamStateGetMeanVar.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
pLibDTSG.tsgDreamStateGetMode.argtypes = [c_void_p, POINTER(c_double)]

pLibDTSG.tsgDreamStateGetRate.restype = c_double
pLibDTSG.tsgDreamStateGetRate.argtypes = [c_void_p]

class DreamState:
    '''
    Wrapper to class TasDREAM::TasmanianDREAM
    '''
    def __init__(self, iNumChains, GridOrIntDimensions):
        '''
        Constructs a new state with given number of chains and dimensions.
        The dimensions can be given as an integer or inferred from
        a Tasmanian.SparseGrid object.

        iNumChains: positive integer indicating the number of chains.
        GridOrIntDimensions: either a positive integer or a non-empty
            Tasmanian.SparseGrid object that gives the dimensions
        '''
        self.TasmanainDreamState = True
        if (isinstance(GridOrIntDimensions, SparseGrid)):
            self.pStatePntr = pLibDTSG.tsgMakeDreamState(iNumChains, GridOrIntDimensions.getNumDimensions())
        else:
            self.pStatePntr = pLibDTSG.tsgMakeDreamState(iNumChains, GridOrIntDimensions)

    def __del__(self):
        '''
        Deletes an instance of the DREAM state.
        '''
        pLibDTSG.tsgDeleteDreamState(self.pStatePntr)

    def getNumDimensions(self):
        '''
        Return the number of dimensions.
        '''
        return pLibDTSG.tsgDreamStateGetDims(self.pStatePntr)

    def getNumChains(self):
        '''
        Return the number of chains.
        '''
        return pLibDTSG.tsgDreamStateGetChains(self.pStatePntr)

    def getNumHistory(self):
        '''
        Return the number of samples saved in the history.
        '''
        return pLibDTSG.tsgDreamStateGetNumHistory(self.pStatePntr)

    def setState(self, llfNewState):
        '''
        Set a new state for the DREAM chains.

        llfNewState: is a two dimensional numpy.ndarray with
            .shape[0] = .getNumChains()
            .shape[1] = .getNumDimensions()
        '''
        iNumChains = self.getNumChains()
        iNumDims = self.getNumDimensions()
        if (llfNewState.shape[0] != iNumChains):
            raise InputError("llfNewState", "llfNewState.shape[0] should match the number of chains")
        if (llfNewState.shape[1] != iNumDims):
            raise InputError("llfNewState", "llfNewState.shape[1] should match the number of dimensions")
        pLibDTSG.tsgDreamStateSet(self.pStatePntr,
                                  np.ctypeslib.as_ctypes(llfNewState.reshape((iNumChains * iNumDims,))))

    def getHistory(self):
        '''
        Returns the saved history in a two dimensional numpy.ndarray
        '''
        iNumDims = self.getNumDimensions()
        iNumHistory = self.getNumHistory()
        aResult = np.empty((iNumDims * iNumHistory,), np.float64)
        pLibDTSG.tsgDreamStateGetHistory(self.pStatePntr, np.ctypeslib.as_ctypes(aResult))
        return aResult.reshape((iNumHistory, iNumDims))

    def getHistoryPDF(self):
        '''
        Returns the saved pdf in a one dimensional numpy.ndarray
        '''
        iNumHistory = self.getNumHistory()
        aResult = np.empty((NumHistory,), np.float64)
        pLibDTSG.tsgDreamStateGetHistoryPDF(self.pStatePntr, np.ctypeslib.as_ctypes(aResult))
        return aResult

    def getHistoryMeanVariance(self):
        '''
        Returns two one dimensional numpy.ndarrays corresponding
        to the mean and variance of the history.
        '''
        iNumDims = self.getNumDimensions()
        aMean = np.empty((iNumDims,))
        aVar  = np.empty((iNumDims,))
        pLibDTSG.tsgDreamStateGetMeanVar(self.pStatePntr,
                                         np.ctypeslib.as_ctypes(aMean),
                                         np.ctypeslib.as_ctypes(aVar))
        return aMean, aVar

    def getApproximateMode(self):
        '''
        Returns the approximate mode, i.e., the sample with highest probability density.
        '''
        iNumDims = self.getNumDimensions()
        aMode = np.empty((iNumDims,))
        pLibDTSG.tsgDreamStateGetMode(self.pStatePntr, np.ctypeslib.as_ctypes(aMode))
        return aMode

    def getAcceptanceRate(self):
        '''
        Return the acceptance rate accumulated with the history.
        '''
        return pLibDTSG.tsgDreamStateGetRate(self.pStatePntr)
