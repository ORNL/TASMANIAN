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

from TasmanianConfig import TasmanianInputError, __path_libsparsegrid__
pLibTSG = cdll.LoadLibrary(__path_libsparsegrid__)

# The ctypes library requires that we manually specify the return types of functions. In C, this is done by header files, so the
# declarations below serve as a header.

# Return types for TasmanianSparseGrid C class methods.
pLibTSG.tsgConstructTasmanianSparseGrid.restype = c_void_p
pLibTSG.tsgGetVersion.restype = c_char_p
pLibTSG.tsgGetLicense.restype = c_char_p
pLibTSG.tsgGetVersionMajor.restype = c_int
pLibTSG.tsgGetVersionMinor.restype = c_int
pLibTSG.tsgIsOpenMPEnabled.restype = c_int
pLibTSG.tsgIsCudaEnabled.restype = c_int
pLibTSG.tsgIsHipEnabled.restype = c_int
pLibTSG.tsgIsDpcppEnabled.restype = c_int
pLibTSG.tsgGetNumDimensions.restype = c_int
pLibTSG.tsgGetNumOutputs.restype = c_int
pLibTSG.tsgGetNumLoaded.restype = c_int
pLibTSG.tsgGetNumNeeded.restype = c_int
pLibTSG.tsgGetNumPoints.restype = c_int
pLibTSG.tsgRead.restype = c_int
pLibTSG.tsgGetAlpha.restype = c_double
pLibTSG.tsgGetBeta.restype = c_double
pLibTSG.tsgGetOrder.restype = c_int
pLibTSG.tsgGetCustomRuleDescription.restype = c_char_p
pLibTSG.tsgGetLoadedPoints.restype = POINTER(c_double)
pLibTSG.tsgGetNeededPoints.restype = POINTER(c_double)
pLibTSG.tsgGetPoints.restype = POINTER(c_double)
pLibTSG.tsgGetQuadratureWeights.restype = POINTER(c_double)
pLibTSG.tsgGetInterpolationWeights.restype = POINTER(c_double)
pLibTSG.tsgBatchGetInterpolationWeights.restype = POINTER(c_double)
pLibTSG.tsgIsSetDomainTransfrom.restype = c_int
pLibTSG.tsgIsSetConformalTransformASIN.restype = c_int
pLibTSG.tsgEvaluateSparseHierarchicalFunctionsGetNZ.restype = c_int
pLibTSG.tsgIsUsingConstruction.restype = c_int
pLibTSG.tsgGetCandidateConstructionPointsVoidPntr.restype = c_void_p
pLibTSG.tsgGetCandidateConstructionPointsSurplusVoidPntr.restype = c_void_p
pLibTSG.tsgGetCandidateConstructionPointsPythonGetNP.restype = c_int
pLibTSG.tsgGetAccelerationType.restype = c_char_p
pLibTSG.tsgIsAccelerationAvailable.restype = c_int
pLibTSG.tsgGetGPUID.restype = c_int
pLibTSG.tsgGetNumGPUs.restype = c_int
pLibTSG.tsgGetGPUMemory.restype = c_int

# Argument types for TasmanianSparseGrid C class methods.
pLibTSG.tsgDestructTasmanianSparseGrid.argtypes = [c_void_p]
pLibTSG.tsgCopySubGrid.argtypes = [c_void_p, c_void_p, c_int, c_int]
pLibTSG.tsgWrite.argtypes = [c_void_p, c_char_p]
pLibTSG.tsgWriteBinary.argtypes = [c_void_p, c_char_p]
pLibTSG.tsgRead.argtypes = [c_void_p, c_char_p]
pLibTSG.tsgMakeGlobalGrid.argtypes = [c_void_p, c_int, c_int, c_int, c_char_p, c_char_p, POINTER(c_int), c_double, c_double, c_char_p, POINTER(c_int)]
pLibTSG.tsgMakeSequenceGrid.argtypes = [c_void_p, c_int, c_int, c_int, c_char_p, c_char_p, POINTER(c_int), POINTER(c_int)]
pLibTSG.tsgMakeLocalPolynomialGrid.argtypes = [c_void_p, c_int, c_int, c_int, c_int, c_char_p, POINTER(c_int)]
pLibTSG.tsgMakeWaveletGrid.argtypes = [c_void_p, c_int, c_int, c_int, c_int, POINTER(c_int)]
pLibTSG.tsgMakeFourierGrid.argtypes = [c_void_p, c_int, c_int, c_int, c_char_p, POINTER(c_int), POINTER(c_int)]
pLibTSG.tsgMakeGridFromCustomTabulated.argtypes = [c_void_p, c_int, c_int, c_int, c_char_p, c_void_p, POINTER(c_int), POINTER(c_int)]
pLibTSG.tsgUpdateGlobalGrid.argtypes = [c_void_p, c_int, c_char_p, POINTER(c_int), POINTER(c_int)]
pLibTSG.tsgUpdateSequenceGrid.argtypes = [c_void_p, c_int, c_char_p, POINTER(c_int), POINTER(c_int)]
pLibTSG.tsgUpdateFourierGrid.argtypes = [c_void_p, c_int, c_char_p, POINTER(c_int), POINTER(c_int)]
pLibTSG.tsgGetAlpha.argtypes = [c_void_p]
pLibTSG.tsgGetBeta.argtypes = [c_void_p]
pLibTSG.tsgGetOrder.argtypes = [c_void_p]
pLibTSG.tsgGetNumDimensions.argtypes = [c_void_p]
pLibTSG.tsgGetNumOutputs.argtypes = [c_void_p]
pLibTSG.tsgCopyRuleChars.argtypes = [c_void_p, c_int, c_char_p, POINTER(c_int)] # char is not really const
pLibTSG.tsgGetCustomRuleDescription.argtypes = [c_void_p]
pLibTSG.tsgGetNumLoaded.argtypes = [c_void_p]
pLibTSG.tsgGetNumNeeded.argtypes = [c_void_p]
pLibTSG.tsgGetNumPoints.argtypes = [c_void_p]
pLibTSG.tsgGetLoadedPoints.argtypes = [c_void_p]
pLibTSG.tsgGetNeededPoints.argtypes = [c_void_p]
pLibTSG.tsgGetPoints.argtypes = [c_void_p]
pLibTSG.tsgGetLoadedPointsStatic.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgGetNeededPointsStatic.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgGetPointsStatic.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgGetQuadratureWeights.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgGetQuadratureWeightsStatic.argtypes = [c_void_p]
pLibTSG.tsgGetInterpolationWeights.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgGetInterpolationWeightsStatic.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
pLibTSG.tsgLoadNeededValues.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgGetLoadedValuesStatic.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgEvaluate.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
pLibTSG.tsgEvaluateFast.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
pLibTSG.tsgIntegrate.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgEvaluateBatch.argtypes = [c_void_p, POINTER(c_double), c_int, POINTER(c_double)]
pLibTSG.tsgBatchGetInterpolationWeights.argtypes = [c_void_p, POINTER(c_double), c_int]
pLibTSG.tsgBatchGetInterpolationWeightsStatic.argtypes = [c_void_p, POINTER(c_double), c_int, POINTER(c_double)]
pLibTSG.tsgIsGlobal.argtypes = [c_void_p]
pLibTSG.tsgIsSequence.argtypes = [c_void_p]
pLibTSG.tsgIsLocalPolynomial.argtypes = [c_void_p]
pLibTSG.tsgIsWavelet.argtypes = [c_void_p]
pLibTSG.tsgIsFourier.argtypes = [c_void_p]
pLibTSG.tsgSetDomainTransform.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
pLibTSG.tsgIsSetDomainTransfrom.argtypes = [c_void_p]
pLibTSG.tsgClearDomainTransform.argtypes = [c_void_p]
pLibTSG.tsgGetDomainTransform.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
pLibTSG.tsgSetConformalTransformASIN.argtypes = [c_void_p, POINTER(c_int)]
pLibTSG.tsgIsSetConformalTransformASIN.argtypes = [c_void_p]
pLibTSG.tsgClearConformalTransform.argtypes = [c_void_p]
pLibTSG.tsgGetConformalTransformASIN.argtypes = [c_void_p, POINTER(c_int)]
pLibTSG.tsgClearLevelLimits.argtypes = [c_void_p]
pLibTSG.tsgGetLevelLimits.argtypes = [c_void_p, POINTER(c_int)]
pLibTSG.tsgSetAnisotropicRefinement.argtypes = [c_void_p, c_char_p, c_int, c_int, POINTER(c_int)]
pLibTSG.tsgEstimateAnisotropicCoefficientsStatic.argtypes = [c_void_p, c_char_p, c_int, POINTER(c_int)]
pLibTSG.tsgSetGlobalSurplusRefinement.argtypes = [c_void_p, c_double, c_int, POINTER(c_int)]
pLibTSG.tsgSetLocalSurplusRefinement.argtypes = [c_void_p, c_double, c_char_p, c_int, POINTER(c_int), POINTER(c_double)]
pLibTSG.tsgClearRefinement.argtypes = [c_void_p]
pLibTSG.tsgMergeRefinement.argtypes = [c_void_p]
pLibTSG.tsgRemovePointsByHierarchicalCoefficient.argtypes = [c_void_p, c_double, c_int, POINTER(c_double)]
pLibTSG.tsgRemovePointsByHierarchicalCoefficientHardCutoff.argtypes = [c_void_p, c_int, c_int, POINTER(c_double)]
pLibTSG.tsgEvaluateHierarchicalFunctions.argtypes = [c_void_p, POINTER(c_double), c_int, POINTER(c_double)]
pLibTSG.tsgGetHierarchicalSupportStatic.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgSetHierarchicalCoefficients.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgEvaluateSparseHierarchicalFunctionsGetNZ.argtypes = [c_void_p, POINTER(c_double), c_int]
pLibTSG.tsgEvaluateSparseHierarchicalFunctionsStatic.argtypes = [c_void_p, POINTER(c_double), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_double)]
pLibTSG.tsgGetHierarchicalCoefficientsStatic.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgIntegrateHierarchicalFunctionsStatic.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgBeginConstruction.argtypes = [c_void_p]
pLibTSG.tsgIsUsingConstruction.argtypes = [c_void_p]
pLibTSG.tsgGetCandidateConstructionPointsVoidPntr.argtypes = [c_void_p, c_char_p, c_int, POINTER(c_int), POINTER(c_int)]
pLibTSG.tsgGetCandidateConstructionPointsSurplusVoidPntr.argtypes = [c_void_p, c_double, c_char_p, c_int, POINTER(c_int), POINTER(c_double)]
pLibTSG.tsgGetCandidateConstructionPointsPythonGetNP.argtypes = [c_void_p, c_void_p]
pLibTSG.tsgGetCandidateConstructionPointsPythonStatic.argtypes = [c_void_p, POINTER(c_double)]
pLibTSG.tsgGetCandidateConstructionPointsPythonDeleteVect.argtypes = [c_void_p]
pLibTSG.tsgLoadConstructedPoint.argtypes = [c_void_p, POINTER(c_double), c_int, POINTER(c_double)]
pLibTSG.tsgFinishConstruction.argtypes = [c_void_p]
pLibTSG.tsgPrintStats.argtypes = [c_void_p]
pLibTSG.tsgEnableAcceleration.argtypes = [c_void_p, c_char_p]
pLibTSG.tsgEnableAccelerationGPU.argtypes = [c_void_p, c_char_p, c_int]
pLibTSG.tsgGetAccelerationType.argtypes = [c_void_p]
pLibTSG.tsgIsAccelerationAvailable.argtypes = [c_char_p]
pLibTSG.tsgSetGPUID.argtypes = [c_void_p, c_int]
pLibTSG.tsgGetGPUID.argtypes = [c_void_p]
pLibTSG.tsgGetGPUMemory.argtypes = [c_int]
pLibTSG.tsgGetGPUName.argtypes = [c_int, c_int, c_char_p, POINTER(c_int)] # not really const here

# Return types for CustomTabulated C class methods.
pLibTSG.tsgConstructCustomTabulated.restype = c_void_p;
pLibTSG.tsgReadCustomTabulated.restype = c_int
pLibTSG.tsgGetNumLevelsCustomTabulated.restype = c_int;
pLibTSG.tsgGetNumPointsCustomTabulated.restype = c_int;
pLibTSG.tsgGetIExactCustomTabulated.restype = c_int;
pLibTSG.tsgGetQExactCustomTabulated.restype = c_int;
pLibTSG.tsgGetDescriptionCustomTabulated.restype = c_char_p;
pLibTSG.tsgMakeCustomTabulatedFromData.restype = c_void_p;

# Argument types for CustomTabulated C class methods.
pLibTSG.tsgDestructCustomTabulated.argtypes = [c_void_p];
pLibTSG.tsgReadCustomTabulated.argtypes = [c_void_p, c_char_p]
pLibTSG.tsgWriteCustomTabulated.argtypes = [c_void_p, c_char_p]
pLibTSG.tsgGetNumLevelsCustomTabulated.argtypes = [c_void_p];
pLibTSG.tsgGetNumPointsCustomTabulated.argtypes = [c_void_p, c_int];
pLibTSG.tsgGetIExactCustomTabulated.argtypes = [c_void_p, c_int];
pLibTSG.tsgGetQExactCustomTabulated.argtypes = [c_void_p, c_int];
pLibTSG.tsgGetDescriptionCustomTabulated.argtypes = [c_void_p];
pLibTSG.tsgGetWeightsNodesStaticCustomTabulated.argtypes = [c_void_p, c_int, POINTER(c_double), POINTER(c_double)]
pLibTSG.tsgMakeCustomTabulatedFromData.argtype = [c_int, POINTER(c_int), POINTER(c_int), POINTER(c_double), POINTER(c_double), c_char_p]

# Specifications for other C functions.
pLibTSG.tsgPythonGetGlobalPolynomialSpace.restype = POINTER(c_int)
pLibTSG.tsgPythonGetGlobalPolynomialSpace.argtypes = [c_void_p, c_int, POINTER(c_int)]
pLibTSG.tsgDeleteInts.argtypes = [POINTER(c_int)]

bTsgPlotting = True
try:
    import matplotlib.pyplot as tsgPlot
except:
    bTsgPlotting = False
    tsgPlot = []

lsTsgGlobalRules = ["clenshaw-curtis", "clenshaw-curtis-zero", "chebyshev", "chebyshev-odd", "gauss-legendre", "gauss-legendre-odd", "gauss-patterson", "leja", "leja-odd", "rleja", "rleja-odd", "rleja-double2", "rleja-double4", "rleja-shifted", "rleja-shifted-even", "rleja-shifted-double", "max-lebesgue", "max-lebesgue-odd", "min-lebesgue", "min-lebesgue-odd", "min-delta", "min-delta-odd", "gauss-chebyshev1", "gauss-chebyshev1-odd", "gauss-chebyshev2", "gauss-chebyshev2-odd", "fejer2", "gauss-gegenbauer", "gauss-gegenbauer-odd", "gauss-jacobi", "gauss-jacobi-odd", "gauss-laguerre", "gauss-laguerre-odd", "gauss-hermite", "gauss-hermite-odd", "custom-tabulated"]
lsTsgGlobalTypes = ["level", "curved", "iptotal", "ipcurved", "qptotal", "qpcurved", "hyperbolic", "iphyperbolic", "qphyperbolic", "tensor", "iptensor", "qptensor"]
lsTsgCurvedTypes = ["curved", "ipcurved", "qpcurved"]
lsTsgRefineTypes = ["classic", "parents", "direction", "fds", "stable"]
lsTsgSequenceRules = ["leja", "rleja", "rleja-shifted", "max-lebesgue", "min-lebesgue", "min-delta"]
lsTsgLocalRules = ["localp", "semi-localp", "localp-zero", "localp-boundary"]
lsTsgAccelTypes = ["none", "cpu-blas", "gpu-default", "gpu-cublas", "gpu-cuda", "gpu-rocblas", "gpu-hip", "gpu-magma"]

class TasmanianSimpleSparseMatrix:
    def __init__(self):
        self.aPntr = []
        self.aIndx = []
        self.aVals = []
        self.iNumRows = 0
        self.iNumCols = 0

    def getDenseForm(self):
        if ((self.iNumRows == 0) or (self.iNumCols == 0)):
            return np.empty([0,0], np.float64)
        aMat = np.zeros([self.iNumRows, self.iNumCols], np.float64 if not np.iscomplexobj(self.aVals) else np.complex128)
        for iI in range(self.iNumRows):
            iJ = self.aPntr[iI]
            while(iJ < self.aPntr[iI+1]):
                aMat[iI, self.aIndx[iJ]] = self.aVals[iJ]
                iJ += 1
        return aMat


class TasmanianSparseGrid:
    def __init__(self):
        '''
        constructor, creates an empty grid instance
        '''
        self.TasmanianSparseGridObject = True
        self.pGrid = pLibTSG.tsgConstructTasmanianSparseGrid()

    def __del__(self):
        '''
        destructor, calls the C++ destructor and releases all memory
        used by this instance of the class
        Make sure to call this to avoid memory leaks
        '''
        pLibTSG.tsgDestructTasmanianSparseGrid(self.pGrid)

    def stringBufferToString(self, pName, iNumChars):
        if (sys.version_info.major == 3):
            S = [s for s in pName]
            sName = ""
            for iI in range(iNumChars):
                sName += str(S[iI], encoding='utf8')
        else:
            S = [s for s in pName]
            sName = ""
            for iI in range(iNumChars):
                sName += S[iI]
        return sName

    def getVersion(self):
        '''
        returns the hardcoded version string from the library

        '''
        sVersion = pLibTSG.tsgGetVersion()
        if (sys.version_info.major == 3):
            sVersion = str(sVersion, encoding='utf8')
        return sVersion

    def getLicense(self):
        '''
        returns the hardcoded license string from the library

        '''
        sLicense = pLibTSG.tsgGetLicense()
        if (sys.version_info.major == 3):
            sLicense = str(sLicense, encoding='utf8')
        return sLicense

    def getVersionMajor(self):
        '''
        returns the hardcoded version major int
        '''
        return pLibTSG.tsgGetVersionMajor()

    def getVersionMinor(self):
        '''
        returns the hardcoded version minor int
        '''
        return pLibTSG.tsgGetVersionMinor()

    def isOpenMPEnabled(self):
        '''
        returns True if the library has been compiled with OpenMP support
        '''
        return (pLibTSG.tsgIsOpenMPEnabled() != 0)

    def isCudaEnabled(self):
        '''
        returns True if the library has been compiled with CUDA support
        '''
        return (pLibTSG.tsgIsCudaEnabled() != 0)

    def isHipEnabled(self):
        '''
        returns True if the library has been compiled with HIP support
        '''
        return (pLibTSG.tsgIsHipEnabled() != 0)

    def isDpcppEnabled(self):
        '''
        returns True if the library has been compiled with DPC++ support
        '''
        return (pLibTSG.tsgIsDpcppEnabled() != 0)

    def read(self, sFilename):
        '''
        reads the grid from a file
        discards any existing grid held by this class

        sFilename: string indicating a grid file where a grid was
                   already written using write from Python or any other
                   Tasmanian interface

        output: boolean
                True: the read was successful
                False: the read failed,
                       check the CLI output for an error message

        '''
        effective_filename = bytes(sFilename, encoding='utf8') if sys.version_info.major == 3 else sFilename
        if pLibTSG.tsgRead(self.pGrid, c_char_p(effective_filename)) == 0:
            raise TasmanianInputError("sFilename", "ERROR: {0:1s} does not appear to be a valid Tasmanian file.".format(sFilename))

    def write(self, sFilename, bUseBinaryFormat = True):
        '''
        writes the grid to a file

        sFilename: string indicating a grid file where a grid will
                   be written

        bUseBinaryFormat: boolean
                True: write to a binary file
                False: write to an ASCII file

        '''
        if (sys.version_info.major == 3):
            sFilename = bytes(sFilename, encoding='utf8')
        if (bUseBinaryFormat):
            pLibTSG.tsgWriteBinary(self.pGrid, c_char_p(sFilename))
        else:
            pLibTSG.tsgWrite(self.pGrid, c_char_p(sFilename))

    def makeGlobalGrid(self, iDimension, iOutputs, iDepth, sType, sRule, liAnisotropicWeights=[], fAlpha=0.0, fBeta=0.0, sCustomFilename="", liLevelLimits=[]):
        '''
        creates a new sparse grid using a global rule
        discards any existing grid held by this class

        iDimension: int (positive)
                    the number of inputs

        iOutputs: int (non-negative)
                  the number of outputs

        iDepth: int (non-negative)
                controls the density of the grid, i.e.,
                the offset for the tensor selection, the meaning of
                iDepth depends on sType
                Example 1: sType == 'iptotal' will give a grid that
                           interpolates exactly all polynomials of
                           degree up to and including iDepth
                Example 2: sType == 'qptotal' will give a grid that
                           integrates exactly all polynomials of degree
                           up to and including iDepth

        sType: string identifying the tensor selection strategy
              'level'     'curved'     'hyperbolic'     'tensor'
              'iptotal'   'ipcurved'   'iphyperbolic'   'iptensor'
              'qptotal'   'qpcurved'   'qphyperbolic'   'qptensor'

        sRule: string (defines the 1-D rule that induces the grid)

           Interpolation rules

              Note: the quadrature induced by those rules is constructed
                    by integrating the interpolant

              'clenshaw-curtis'    'clenshaw-curtis-zero'      'fejer2'
              'rleja'    'rleja-odd'  'rleja-double2'   'rleja-double4'
              'rleja-shifted'   'rleja-shifted-even'
              'max-lebesgue'    'max-lebesgue-odd'
              'min-lebesgue'    'min-lebesgue-odd'
              'leja'            'leja-odd'
              'min-delta'       'min-delta-odd'

              'chebyshev'       'chebyshev-odd'
                approximation using roots of Chebyshev polynomials
                non-nested case (in contrast to Clenshaw-Curtis nodes)

           Quadrature rules, the weights target exactness with respect
                             to the highest polynomial degree possible

               'gauss-legendre'  'gauss-legendre-odd'
                approximation using roots of polynomials orthogonal in
                measure Uniform

               'gauss-patterson'  (a.k.a. nested Gauss-Legendre)
                Note: the nodes and weights are hard-coded hence there
                is a limit on the highest possible depth

               'gauss-chebyshev1'  'gauss-chebyshev1-odd'
               'gauss-chebyshev2'  'gauss-chebyshev2-odd'
                 approximation using roots of polynomials orthogonal in
                 measures  1/sqrt(1-x^2) and sqrt(1-x^2)  (respectively)

              'gauss-gegenbauer'  'gauss-gegenbauer-odd'
                approximation using roots of polynomials orthogonal in
                measure (1-x^2)^alpha

              'gauss-jacobi'
                approximation using roots of polynomials orthogonal in
                measure (1-x)^alpha * (1+x)^beta

              'gauss-laguerre'
                approximation using roots of polynomials orthogonal in
                measure x^alpha * epx(-x)

              'gauss-hermite'  'gauss-hermite-odd'
                approximation using roots of polynomials orthogonal in
                measure |x|^alpha * epx(-x^2)

        liAnisotropicWeights: list or numpy.ndarray of weights
                              length must be iDimension or 2*iDimension
                              the first iDimension weights
                                                       must be positive
                              see the manual for details

        fAlpha, fBeta: floats
              fAlpha : the alpha parameter for Gegenbauer, Jacobi,
                       Hermite and Laguerre rules
              fBeta  : the beta parameter for Jacobi rules

        sCustomRule: string giving the path to the file with
                     custom-tabulated rule

        '''
        if (iDimension <= 0):
            raise TasmanianInputError("iDimension", "ERROR: dimension should be a positive integer")
        if (iOutputs < 0):
            raise TasmanianInputError("iOutputs", "ERROR: outputs should be a non-negative integer")
        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")
        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")
        if (sRule not in lsTsgGlobalRules):
            raise TasmanianInputError("sRule", "ERROR: invalid global rule, see TasmanianSG.lsTsgGlobalRules for list of accepted global rules")
        pAnisoWeights = None
        if (len(liAnisotropicWeights) > 0):
            if (sType in lsTsgCurvedTypes):
                iNumWeights = 2*iDimension
            else:
                iNumWeights = iDimension
            if (len(liAnisotropicWeights) != iNumWeights):
                raise TasmanianInputError("liAnisotropicWeights", "ERROR: wrong number of liAnisotropicWeights, sType '{0:s}' needs {1:1d} weights but len(liAnisotropicWeights) == {2:1d}".format(sType, iNumWeights, len(liAnisotropicWeights)))
            else:
                aAWeights = np.array([liAnisotropicWeights[i] for i in range(iNumWeights)], np.int32)
                pAnisoWeights = np.ctypeslib.as_ctypes(aAWeights)

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')
            sRule = bytes(sRule, encoding='utf8')
            if (sCustomFilename):
                sCustomFilename = bytes(sCustomFilename, encoding='utf8')

        pCustomRule = None
        if (sCustomFilename):
            pCustomRule = c_char_p(sCustomFilename)

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        pLibTSG.tsgMakeGlobalGrid(self.pGrid, iDimension, iOutputs, iDepth, c_char_p(sType), c_char_p(sRule), pAnisoWeights, c_double(fAlpha), c_double(fBeta), pCustomRule, pLevelLimits)

    def makeSequenceGrid(self, iDimension, iOutputs, iDepth, sType, sRule, liAnisotropicWeights=[], liLevelLimits=[]):
        '''
        creates a new sparse grid using a sequence rule
        discards any existing grid held by this class

        iDimension: int (positive)
              the number of inputs

        iOutputs: int (non-negative)
              the number of outputs

        iDepth: int (non-negative)
                controls the density of the grid, i.e.,
                the offset for the tensor selection, the meaning of
                iDepth depends on sType
                Example 1: sType == 'iptotal' will give a grid that
                           interpolates exactly all polynomials of
                           degree up to and including iDepth
                Example 2: sType == 'qptotal' will give a grid that
                           integrates exactly all polynomials of degree
                           up to and including iDepth

        sType: string identifying the tensor selection strategy
              'level'     'curved'     'hyperbolic'     'tensor'
              'iptotal'   'ipcurved'   'iphyperbolic'   'iptensor'
              'qptotal'   'qpcurved'   'qphyperbolic'   'qptensor'

        sRule: string (defines the 1-D rule that induces the grid)
              'leja'       'rleja'      'rleja-shifted'
              'max-lebesgue'   'min-lebesgue'   'min-delta'

        liAnisotropicWeights: list or numpy.ndarray of weights
                              length must be iDimension or 2*iDimension
                              the first iDimension weights
                                                       must be positive
                              see the manual for details

        '''
        if (iDimension <= 0):
            raise TasmanianInputError("iDimension", "ERROR: dimension should be a positive integer")
        if (iOutputs < 0):
            raise TasmanianInputError("iOutputs", "ERROR: outputs should be a non-negative integer")
        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")
        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")
        if (sRule not in lsTsgSequenceRules):
            raise TasmanianInputError("sRule", "ERROR: invalid sequence rule, see TasmanianSG.lsTsgSequenceRules for list of accepted sequence rules")
        pAnisoWeights = None
        if (len(liAnisotropicWeights) > 0):
            if (sType in lsTsgCurvedTypes):
                iNumWeights = 2*iDimension
            else:
                iNumWeights = iDimension
            if (len(liAnisotropicWeights) != iNumWeights):
                raise TasmanianInputError("liAnisotropicWeights", "ERROR: wrong number of liAnisotropicWeights, sType '{0:s}' needs {1:1d} weights but len(liAnisotropicWeights) == {2:1d}".format(sType, iNumWeights, len(liAnisotropicWeights)))
            else:
                pAnisoWeights = (c_int*iNumWeights)()
                for iI in range(iNumWeights):
                    pAnisoWeights[iI] = liAnisotropicWeights[iI]

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')
            sRule = bytes(sRule, encoding='utf8')

        pLibTSG.tsgMakeSequenceGrid(self.pGrid, iDimension, iOutputs, iDepth, c_char_p(sType), c_char_p(sRule), pAnisoWeights, pLevelLimits)

    def makeLocalPolynomialGrid(self, iDimension, iOutputs, iDepth, iOrder=1, sRule="localp", liLevelLimits=[]):
        '''
        creates a new sparse grid using a local polynomial rule
        discards any existing grid held by this class

        iDimension: int (positive)
              the number of inputs

        iOutputs: int (non-negative)
              the number of outputs

        iDepth: int (non-negative)
                controls the density of the grid, i.e.,
                the number of levels to use

        iOrder: int (must be -1 or bigger)
                -1 indicates largest possible order
                 1 means linear, 2 means quadratic, etc.
                 0 means piece-wise constant, it has different hierarchy
                   then the other orders, most notably the 1D rule
                   triples the number of points per level (as opposed
                   to double for the other cases)

        sRule: string (defines the 1-D rule that induces the grid)
              'localp'  'localp-zero'  'semi-localp'  'localp-boundary'

        '''
        if (iDimension <= 0):
            raise TasmanianInputError("iDimension", "ERROR: dimension should be a positive integer")
        if (iOutputs < 0):
            raise TasmanianInputError("iOutputs", "ERROR: outputs should be a non-negative integer")
        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")
        if (iOrder < -1):
            raise TasmanianInputError("iOrder", "ERROR: order should be a non-negative integer")
        if (sRule not in lsTsgLocalRules):
            raise TasmanianInputError("sRule", "ERROR: invalid local polynomial rule, see TasmanianSG.lsTsgLocalRules for list of accepted sequence rules")

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        if (sys.version_info.major == 3):
            sRule = bytes(sRule, encoding='utf8')

        pLibTSG.tsgMakeLocalPolynomialGrid(self.pGrid, iDimension, iOutputs, iDepth, iOrder, c_char_p(sRule), pLevelLimits)

    def makeWaveletGrid(self, iDimension, iOutputs, iDepth, iOrder=1, liLevelLimits=[]):
        '''
        creates a new sparse grid using a wavelet rule
        discards any existing grid held by this class

        iDimension: int (positive)
              the number of inputs

        iOutputs: int (non-negative)
              the number of outputs

        iDepth: int (non-negative)
                controls the density of the grid, i.e.,
                the number of levels to use

        iOrder: int (must be 1 or 3)
              only wavelets of order 1 and 3 are implemented

        '''
        if (iDimension <= 0):
            raise TasmanianInputError("iDimension", "ERROR: dimension should be a positive integer")
        if (iOutputs < 0):
            raise TasmanianInputError("iOutputs", "ERROR: outputs should be a non-negative integer")
        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")
        if (iOrder not in [1, 3]):
            raise TasmanianInputError("iOrder", "ERROR: order should be either 1 or 3 (only linear and cubic wavelets are available)")

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        pLibTSG.tsgMakeWaveletGrid(self.pGrid, iDimension, iOutputs, iDepth, iOrder, pLevelLimits)

    def makeFourierGrid(self, iDimension, iOutputs, iDepth, sType, liAnisotropicWeights=[], liLevelLimits=[]):
        '''
        creates a new sparse grid using a Fourier rule
        discards any existing grid held by this class

        iDimension: int (positive)
              the number of inputs

        iOutputs: int (non-negative)
              the number of outputs

        iDepth: int (non-negative)
                controls the density of the grid, i.e.,
                the offset for the tensor selection, the meaning of
                iDepth depends on sType

        sType: string identifying the tensor selection strategy
              'level'     'curved'     'hyperbolic'     'tensor'
              'iptotal'   'ipcurved'   'iphyperbolic'   'iptensor'
              'qptotal'   'qpcurved'   'qphyperbolic'   'qptensor'

        liAnisotropicWeights: list or numpy.ndarray of weights
                              length must be iDimension or 2*iDimension
                              the first iDimension weights
                                                       must be positive
                              see the manual for details

        '''
        if (iDimension <= 0):
            raise TasmanianInputError("iDimension", "ERROR: dimension should be a positive integer")
        if (iOutputs < 0):
            raise TasmanianInputError("iOutputs", "ERROR: outputs should be a non-negative integer")
        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")
        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")
        pAnisoWeights = None
        if (len(liAnisotropicWeights) > 0):
            if (sType in lsTsgCurvedTypes):
                iNumWeights = 2*iDimension
            else:
                iNumWeights = iDimension
            if (len(liAnisotropicWeights) != iNumWeights):
                raise TasmanianInputError("liAnisotropicWeights", "ERROR: wrong number of liAnisotropicWeights, sType '{0:s}' needs {1:1d} weights but len(liAnisotropicWeights) == {2:1d}".format(sType, iNumWeights, len(liAnisotropicWeights)))
            else:
                pAnisoWeights = (c_int*iNumWeights)()
                for iI in range(iNumWeights):
                    pAnisoWeights[iI] = liAnisotropicWeights[iI]

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')

        pLibTSG.tsgMakeFourierGrid(self.pGrid, iDimension, iOutputs, iDepth, c_char_p(sType), pAnisoWeights, pLevelLimits)

    def makeGlobalGridCustom(self, iDimension, iOutputs, iDepth, sType, pCustomTabulated, liAnisotropicWeights=[], liLevelLimits=[]):
        '''
        Creates a new sparse grid using a CustomTabulated object. Discards any existing grid held by this class.

        Uses the same inputs as in makeGlobalGrid(), but with an additional input named pCustomTabulated. This new input should be
        a Python CustomTabulated object, which will be used to populate the weights and nodes of the output grid.

        See the manual for more details.
        '''
        if not hasattr(pCustomTabulated, "CustomTabulatedObject"):
            raise TasmanianInputError("pCustomTabulated", "ERROR: pCustomTabulated must be an instance of CustomTabulated")
        self.__testMakeGlobalGrid(iDimension, iOutputs, iDepth, sType, liAnisotropicWeights, liLevelLimits)
        pAnisoWeights = None
        if (len(liAnisotropicWeights) > 0):
            iNumWeights = 2*iDimension if sType in lsTsgCurvedTypes else iDimension
            pAnisoWeights = (c_int*iNumWeights)()
            for iI in range(iNumWeights):
                pAnisoWeights[iI] = liAnisotropicWeights[iI]
        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]
        effective_sType = bytes(sType, encoding='utf8') if (sys.version_info.major == 3) else sType;
        pLibTSG.tsgMakeGridFromCustomTabulated(c_void_p(self.pGrid), c_int(iDimension), c_int(iOutputs), c_int(iDepth),
                                               c_char_p(effective_sType), c_void_p(pCustomTabulated.pCustomTabulated),
                                               pAnisoWeights, pLevelLimits)

    def copyGrid(self, pGrid, iOutputsBegin = 0, iOutputsEnd = -1):
        '''
        accepts an instance of TasmanianSparseGrid class and creates
        a hard copy of the class and all included data
        original class is not modified

        pGrid: instance of TasmanianSparseGrid class
            the source for the copy

        iOutputsBegin: integer indicating the first output to copy

        iOutputsEnd: integer one bigger than the last output to copy
                     if set to -1, all outputs from iOutputsBegin to
                     the end will be copied

        Examples:

        grid.copyGrid(other, 0, -1) # copy all outputs (default)
        grid.copyGrid(other, 0, other.getNumOutputs()) # also copy all
        grid.copyGrid(other, 0, 3) # copy outputs 0, 1, and 2
        grid.copyGrid(other, 1, 4) # copy outputs 1, 2, and 3

        '''
        if (not isinstance(pGrid, TasmanianSparseGrid)):
            raise TasmanianInputError("pGrid", "ERROR: pGrid must be an instance of TasmanianSparseGrid")

        pLibTSG.tsgCopySubGrid(self.pGrid, pGrid.pGrid, iOutputsBegin, iOutputsEnd)

    def updateGlobalGrid(self, iDepth, sType, liAnisotropicWeights=[], liLevelLimits=[]):
        '''
        adds the points defined by depth, type and anisotropy
        to the existing grid

        basically, the same as calling makeGlobalGrid with sRule,
                   fAlpha and fBeta of this grid and the new iDepth,
                   sType and liAnisotropicWeights then adding the
                   resulting points to the current grid

        inputs: see help(TasmanianSG.TasmanianSparseGrid.makeGlobalGrid)

        '''
        if (not self.isGlobal()):
            raise TasmanianInputError("updateGlobalGrid", "ERROR: calling updateGlobalGrid for a grid that is not global")

        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")

        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")

        iDimension = self.getNumDimensions()
        pAnisoWeights = None
        if (len(liAnisotropicWeights) > 0):
            if (sType in lsTsgCurvedTypes):
                iNumWeights = 2*iDimension
            else:
                iNumWeights = iDimension
            if (len(liAnisotropicWeights) != iNumWeights):
                raise TasmanianInputError("liAnisotropicWeights", "ERROR: wrong number of liAnisotropicWeights, sType '{0:s}' needs {1:1d} weights but len(liAnisotropicWeights) == {2:1d}".format(sType, iNumWeights, len(liAnisotropicWeights)))
            else:
                pAnisoWeights = (c_int*iNumWeights)()
                for iI in range(iNumWeights):
                    pAnisoWeights[iI] = liAnisotropicWeights[iI]

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        pLibTSG.tsgUpdateGlobalGrid(self.pGrid, iDepth, sType, pAnisoWeights, pLevelLimits)

    def updateSequenceGrid(self, iDepth, sType, liAnisotropicWeights=[], liLevelLimits=[]):
        '''
        adds the points defined by depth, type and anisotropy
        to the existing grid

        basically, the same as calling makeSequenceGrid() with sRule,
                   of this grid and the new iDepth, sType and
                   liAnisotropicWeights then adding the resulting points
                   to the current grid

        inputs: see help(Tasmanian.SparseGrid.makeGlobalGrid)

        '''
        if (not self.isSequence()):
            raise TasmanianInputError("updateSequenceGrid", "ERROR: calling updateSequenceGrid for a grid that is not a sequence grid")

        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")

        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")

        iDimension = self.getNumDimensions()
        pAnisoWeights = None
        if (len(liAnisotropicWeights) > 0):
            if (sType in lsTsgCurvedTypes):
                iNumWeights = 2*iDimension
            else:
                iNumWeights = iDimension
            if (len(liAnisotropicWeights) != iNumWeights):
                raise TasmanianInputError("liAnisotropicWeights", "ERROR: wrong number of liAnisotropicWeights, sType '{0:s}' needs {1:1d} weights but len(liAnisotropicWeights) == {2:1d}".format(sType, iNumWeights, len(liAnisotropicWeights)))
            else:
                pAnisoWeights = (c_int*iNumWeights)()
                for iI in range(iNumWeights):
                    pAnisoWeights[iI] = liAnisotropicWeights[iI]

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        pLibTSG.tsgUpdateSequenceGrid(self.pGrid, iDepth, sType, pAnisoWeights, pLevelLimits)

    def updateFourierGrid(self, iDepth, sType, liAnisotropicWeights=[], liLevelLimits=[]):
        '''
        adds the points defined by depth, type and anisotropy
        to the existing grid

        basically, the same as calling makeFourierGrid() with sRule,
                   of this grid and the new iDepth, sType and
                   liAnisotropicWeights then adding the resulting points
                   to the current grid

        inputs: see help(Tasmanian.SparseGrid.makeGlobalGrid)

        '''
        if (not self.isFourier()):
            raise TasmanianInputError("updateFourierGrid", "ERROR: calling updateFourierGrid for a grid that is not a Fourier grid")

        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")

        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see Tasmanian.TasmanianSG.lsTsgGlobalTypes for list of accepted types")

        iDimension = self.getNumDimensions()
        pAnisoWeights = None
        if (len(liAnisotropicWeights) > 0):
            if (sType in lsTsgCurvedTypes):
                iNumWeights = 2*iDimension
            else:
                iNumWeights = iDimension
            if (len(liAnisotropicWeights) != iNumWeights):
                raise TasmanianInputError("liAnisotropicWeights", "ERROR: wrong number of liAnisotropicWeights, sType '{0:s}' needs {1:1d} weights but len(liAnisotropicWeights) == {2:1d}".format(sType, iNumWeights, len(liAnisotropicWeights)))
            else:
                pAnisoWeights = (c_int*iNumWeights)()
                for iI in range(iNumWeights):
                    pAnisoWeights[iI] = liAnisotropicWeights[iI]

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        pLibTSG.tsgUpdateFourierGrid(self.pGrid, iDepth, sType, pAnisoWeights, pLevelLimits)

    def getAlpha(self):
        '''
        returns the value of fAlpha in the call to makeGlobalGrid
        if makeGlobalGrid has not been called, returns 0.0

        '''
        return pLibTSG.tsgGetAlpha(self.pGrid)

    def getBeta(self):
        '''
        returns the value of fBeta in the call to makeGlobalGrid
        if makeGlobalGrid has not been called, returns 0.0

        '''
        return pLibTSG.tsgGetBeta(self.pGrid)

    def getOrder(self):
        '''
        returns the value of iOrder in the call to
        makeLocalPolynomialGrid or makeWaveletGrid

        if makeLocalPolynomialGrid and makeWaveletGrid
        have not been called, returns -1

        '''
        return pLibTSG.tsgGetOrder(self.pGrid)

    def getNumDimensions(self):
        '''
        returns the value of iDimension in the make***Grid command
        if no grid has been made, it returns 0

        '''
        return pLibTSG.tsgGetNumDimensions(self.pGrid)

    def getNumOutputs(self):
        '''
        returns the value of iOutputs in the make***Grid command
        if no grid has been made, it returns 0

        '''
        return pLibTSG.tsgGetNumOutputs(self.pGrid)

    def getRule(self):
        '''
        returns the value of sRule in the make***Grid command
        if makeWaveletGrid is used, returns "wavelet"
        if no grid has been made, it returns "unknown"

        '''
        pName = create_string_buffer(128);
        iNumChars = np.array([0], np.int32)
        pLibTSG.tsgCopyRuleChars(self.pGrid, 128, pName, np.ctypeslib.as_ctypes(iNumChars))
        return self.stringBufferToString(pName, iNumChars[0])

    def getCustomRuleDescription(self):
        '''
        returns the description provided in the custom rule file
        if not using a custom grid, returns ""

        '''
        if ("custom-tabulated" in self.getRule()):
            sRule = pLibTSG.tsgGetCustomRuleDescription(self.pGrid)
            if (sys.version_info.major == 3):
                sRule = str(sRule, encoding='utf8')
            return sRule
        else:
            return ""

    def getNumLoaded(self):
        '''
        returns the number of points loaded in the existing interpolant

        '''
        return pLibTSG.tsgGetNumLoaded(self.pGrid)

    def getNumNeeded(self):
        '''
        returns the number of points needed to form the interpolant or
        form the next interpolant following a refinement

        '''
        return pLibTSG.tsgGetNumNeeded(self.pGrid)

    def getNumPoints(self):
        '''
        if points have been loaded, returns the same as getNumLoaded()
        otherwise, returns the same as getNumNeeded()

        '''
        return pLibTSG.tsgGetNumPoints(self.pGrid)

    def getLoadedPoints(self):
        '''
        returns the points loaded in the existing interpolant

        output: a 2-D numpy.ndarray of size getNumLoaded() X iDimension
            reach row corresponds to one point
            if (getNumLoaded() == 0): returns numpy.empty([0,0])

        '''
        iNumDims = self.getNumDimensions()
        iNumPoints = self.getNumLoaded()
        if (iNumPoints == 0):
            return np.empty([0, 0], np.float64)
        aPoints = np.empty([iNumPoints * iNumDims], np.float64)
        pLibTSG.tsgGetLoadedPointsStatic(self.pGrid, np.ctypeslib.as_ctypes(aPoints))
        return aPoints.reshape([iNumPoints, iNumDims])

    def getNeededPoints(self):
        '''
        returns the points needed to form the interpolant or the next
        level of refinement following a set***Refinement() call

        output: 2-D numpy.ndarray of size getNumNeeded() X iDimension
            reach row corresponds to one point
            if (getNumNeeded() == 0): returns numpy.empty([0,0])

        '''
        iNumDims = self.getNumDimensions()
        iNumPoints = self.getNumNeeded()
        if (iNumPoints == 0):
            return np.empty([0, 0], np.float64)
        aPoints = np.empty([iNumPoints * iNumDims], np.float64)
        pLibTSG.tsgGetNeededPointsStatic(self.pGrid, np.ctypeslib.as_ctypes(aPoints))
        return aPoints.reshape([iNumPoints, iNumDims])

    def getPoints(self):
        '''
        if points have been loaded, gives the same as getLoadedPoints()
        otherwise, returns the same as getNeededPoints()

        '''
        iNumDims = self.getNumDimensions()
        iNumPoints = self.getNumPoints()
        if (iNumPoints == 0):
            return np.empty([0, 0], np.float64)
        aPoints = np.empty([iNumPoints * iNumDims], np.float64)
        pLibTSG.tsgGetPointsStatic(self.pGrid, np.ctypeslib.as_ctypes(aPoints))
        return aPoints.reshape([iNumPoints, iNumDims])

    def getQuadratureWeights(self):
        '''
        returns the quadrature weights associated with
        the points in getPoints()

        output: a 1-D numpy.ndarray of length getNumPoints()
                the order of the weights matches
                the order in getPoints()

        '''
        iNumPoints = self.getNumPoints()
        if (iNumPoints == 0):
            return np.empty([0], np.float64)
        aWeights = np.empty([iNumPoints], np.float64)
        pLibTSG.tsgGetQuadratureWeightsStatic(self.pGrid, np.ctypeslib.as_ctypes(aWeights))
        return aWeights

    def getInterpolationWeights(self, lfX):
        '''
        returns the interpolation weights associated with the points
        in getPoints()

        lfX: a 1-D numpy.ndarray with length iDimensions
             the entries indicate the points for evaluating the weights

        output: a 1-D numpy.ndarray of length getNumPoints()
            the order of the weights matches the order in getPoints()

        '''
        iNumX = len(lfX)
        if (iNumX != self.getNumDimensions()):
            raise TasmanianInputError("lfX", "ERROR: len(lfX) should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumX))
        iNumPoints = self.getNumPoints()
        if (iNumPoints == 0):
            return np.empty([0], np.float64)
        aWeights = np.empty([iNumPoints], np.float64)
        pLibTSG.tsgGetInterpolationWeightsStatic(self.pGrid, np.ctypeslib.as_ctypes(lfX), np.ctypeslib.as_ctypes(aWeights))
        return aWeights

    def getInterpolationWeightsBatch(self, llfX):
        '''
        returns the interpolation weights associated with the points
        in getPoints()

        finds multiple weights with a single library call
        uses OpenMP if enabled in libtasmaniansparsegrids.so

        llfX: a 2-D numpy.ndarray with second dimension iDimensions
              each row in the array is a single requested point

        output: a 2-D numpy.ndarray
                with dimensions llfX.shape[0] X getNumPoints()
                each row corresponds to the weight for one row of llfX

        '''
        if (len(llfX.shape) != 2):
            raise TasmanianInputError("llfX", "ERROR: llfX should be a 2-D numpy.ndarray instread it has dimension {0:1d}".format(len(llfX.shape)))
        iNumX = llfX.shape[0]
        if (iNumX == 0):
            return np.empty([0, self.getNumPoints()], np.float64)
        iNumDim = llfX.shape[1]
        if (iNumDim != self.getNumDimensions()):
            raise TasmanianInputError("llfX", "ERROR: llfX.shape[1] should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumDim))
        iNumPoints = self.getNumPoints()
        if (iNumPoints == 0):
            return np.empty([0, 0], np.float64)
        aWeights = np.empty([iNumX, iNumPoints], np.float64)
        pLibTSG.tsgBatchGetInterpolationWeightsStatic(self.pGrid,
            np.ctypeslib.as_ctypes(llfX.reshape([iNumX * iNumDim])), iNumX,
            np.ctypeslib.as_ctypes(aWeights.reshape([iNumX * iNumPoints])))
        return aWeights

    def loadNeededValues(self, llfVals):
        '''
        loads the values of the target function at the needed points
        if there are no needed points, this reset the currently loaded
        values

        llfVals: a 2-D numpy.ndarray
                 with dimensions getNumNeeded() X iOutputs
                 each row corresponds to the values of the outputs at
                 the corresponding needed point. The order and leading
                 dimension must match the points obtained form
                 getNeededPoints()

        '''
        if (len(llfVals.shape) != 2):
            raise TasmanianInputError("llfVals", "ERROR: llfVals should be a 2-D numpy.ndarray, instead it has {0:1d} dimensions".format(len(llfVals.shape)))
        if (self.getNumNeeded() == 0):
            if (llfVals.shape[0] != self.getNumLoaded()):
                raise TasmanianInputError("llfVals", "ERROR: leading dimension of llfVals is {0:1d} but the number of current points is {1:1d}".format(llfVals.shape[0], self.getNumLoaded()))
        elif (llfVals.shape[0] != self.getNumNeeded()):
            raise TasmanianInputError("llfVals", "ERROR: leading dimension of llfVals is {0:1d} but the number of needed points is {1:1d}".format(llfVals.shape[0], self.getNumNeeded()))
        if (llfVals.shape[1] != self.getNumOutputs()):
            raise TasmanianInputError("llfVals", "ERROR: second dimension of llfVals is {0:1d} but the number of outputs is set to {1:1d}".format(llfVals.shape[1], self.getNumOutputs()))
        iNumPoints = llfVals.shape[0]
        iNumDims = llfVals.shape[1]
        pLibTSG.tsgLoadNeededValues(self.pGrid, np.ctypeslib.as_ctypes(llfVals.reshape([iNumPoints * iNumDims])))

    def loadNeededPoints(self, llfVals):
        '''
        Alias of loadNeededValues().
        '''
        self.loadNeededValues(llfVals)

    def getLoadedValues(self):
        '''
        Returns the model values as given to Tasmanian by the loadNeededPoints() method.
        The ordering will match the current internal ordering, e.g., mixing the different
        model values from different refinement iterations.

        Returns a two dimensional numpy.ndarray with size getNumPoints() by getNumOutputs()
        '''
        iNumPoints  = self.getNumPoints()
        iNumOutputs = self.getNumOutputs()
        if (iNumPoints == 0 or iNumOutputs == 0):
            return np.empty([0, 0], np.float64)
        aVals = np.empty((iNumPoints * iNumOutputs,), np.float64)
        pLibTSG.tsgGetLoadedValuesStatic(self.pGrid, np.ctypeslib.as_ctypes(aVals))
        return aVals.reshape((iNumPoints, iNumOutputs))

    def evaluateThreadSafe(self, lfX):
        '''
        evaluates the intepolant at a single points of interest and
        returns the result
        This is the thread safe version, but it does not use
        acceleration of any type

        this should be called after the grid has been created and after
        values have been loaded

        lfX: a 1-D numpy.ndarray with length iDimensions
             the entries indicate the points for evaluating the weights

        output: returns a 1-D numpy.ndarray of length iOutputs
            the values of the interpolant at lfX

        '''
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("evaluateThreadSafe", "ERROR: cannot call evaluate for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (len(lfX.shape) != 1):
            raise TasmanianInputError("lfX", "ERROR: lfX should be 1D numpy array")
        iNumX = lfX.shape[0]
        if (iNumX != self.getNumDimensions()):
            raise TasmanianInputError("lfX", "ERROR: lfX should have lenth {0:1d} instead it has length {1:1d}".format(self.getNumDimensions(),iNumX))
        iNumOutputs = self.getNumOutputs()
        aY = np.empty([iNumOutputs], np.float64)
        pLibTSG.tsgEvaluate(self.pGrid, np.ctypeslib.as_ctypes(lfX), np.ctypeslib.as_ctypes(aY))
        return aY

    def evaluate(self, lfX):
        '''
        evaluates the intepolant at a single points of interest and
        returns the result
        This is the accelerated version using the selected acceleration
        type, but it is potentially not thread safe

        this should be called after the grid has been created and after
        values have been loaded

        lfX: a 1-D numpy.ndarray with length iDimensions
             the entries indicate the points for evaluating the weights

        output: returns a 1-D numpy.ndarray of length iOutputs
            the values of the interpolant at lfX

        '''
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("evaluate", "ERROR: cannot call evaluate for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (len(lfX.shape) != 1):
            raise TasmanianInputError("lfX", "ERROR: lfX should be 1D numpy array")
        iNumX = lfX.shape[0]
        if (iNumX != self.getNumDimensions()):
            raise TasmanianInputError("lfX", "ERROR: lfX should have lenth {0:1d} instead it has length {1:1d}".format(self.getNumDimensions(),iNumX))
        iNumOutputs = self.getNumOutputs()
        aY = np.empty([iNumOutputs], np.float64)
        pLibTSG.tsgEvaluateFast(self.pGrid, np.ctypeslib.as_ctypes(lfX), np.ctypeslib.as_ctypes(aY))
        return aY

    def evaluateBatch(self, llfX):
        '''
        evaluates the intepolant at the points of interest and returns
        the result

        this should be called after the grid has been created and after
        values have been loaded

        llfX: a 2-D numpy.ndarray
              with second dimension equal to iDimensions
              each row in the array is a single requested point

        output: a 2-D numpy.ndarray
                with dimensions llfX.shape[0] X iOutputs
                each row corresponds to the value of the interpolant
                for one row of llfX

        '''
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("evaluateBatch", "ERROR: cannot call evaluateBatch for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (len(llfX.shape) != 2):
            raise TasmanianInputError("llfX", "ERROR: llfX should be a 2-D numpy.ndarray instread it has dimension {0:1d}".format(len(llfX.shape)))
        iNumX = llfX.shape[0]
        if (iNumX == 0):
            return np.empty([0, self.getNumOutputs()], np.float64)
        iNumDim = llfX.shape[1]
        if (iNumDim != self.getNumDimensions()):
            raise TasmanianInputError("llfX", "ERROR: llfX.shape[1] should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumDim))
        iNumOutputs = self.getNumOutputs()
        aY = np.empty([iNumX, iNumOutputs], np.float64)
        # np.ctypeslib.as_ctypes(llfX.reshape([iNumX*iNumDim,])) messes up, the first 4 entries randomly get set to machine eps (10^-310) and 0
        lfX = llfX.reshape([iNumX*iNumDim,])
        pLibTSG.tsgEvaluateBatch(self.pGrid, np.ctypeslib.as_ctypes(lfX), iNumX, np.ctypeslib.as_ctypes(aY.reshape([iNumX*iNumOutputs,])))
        return aY

    def integrate(self):
        '''
        returns the integral of the interpolant

        output: returns a 1-D numpy.ndarray of length iOutputs
            the integral of the interpolant

        '''
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("integrate", "ERROR: cannot call integrate for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        iNumOutputs = self.getNumOutputs()
        aQ = np.empty([iNumOutputs], np.float64)
        pLibTSG.tsgIntegrate(self.pGrid, np.ctypeslib.as_ctypes(aQ))
        return aQ

    def isGlobal(self):
        '''
        returns True if using a global grid

        '''
        return (pLibTSG.tsgIsGlobal(self.pGrid) != 0)

    def isSequence(self):
        '''
        returns True if using a sequence grid

        '''
        return (pLibTSG.tsgIsSequence(self.pGrid) != 0)

    def isLocalPolynomial(self):
        '''
        returns True if using a local polynomial grid

        '''
        return (pLibTSG.tsgIsLocalPolynomial(self.pGrid) != 0)

    def isWavelet(self):
        '''
        returns True if using a local wavelet grid

        '''
        return (pLibTSG.tsgIsWavelet(self.pGrid) != 0)

    def isFourier(self):
        '''
        returns True if using a Fourier grid
        '''
        return (pLibTSG.tsgIsFourier(self.pGrid) != 0)

    def setDomainTransform(self, llfTransform):
        '''
        sets the lower and upper bound for each dimension

        Note: gauss-laguerre and gauss-hermite rules are defined on
              unbounded domain, in which case  this sets the
              shift and scale parameters, consult the manual

        llfTransform: a 2-D numpy.ndarray of size iDimension X 2
                      transform specifies the lower and upper bound
                      of the domain in each direction.

                      For gauss-laguerre and gauss-hermite grids, the
                      transform gives the a and b parameters of the
                      weights
                       exp(-b (x - a))
                       exp(-b (x - a)^2)

        '''
        lShape = llfTransform.shape
        if (len(lShape) != 2):
            raise TasmanianInputError("llfTransform", "ERROR: llfTransform should be a 2-D numpy.ndarray")
        if (lShape[0] != self.getNumDimensions()):
            raise TasmanianInputError("llfTransform", "ERROR: the first dimension of llfTransform is {0:1d} and it should match iDimension: {1:1d}".format(lShape[0], self.getNumDimensions()))
        if (lShape[1] != 2):
            raise TasmanianInputError("llfTransform", "ERROR: the second dimension of llfTransform is {0:1d} and it should be 2".format(lShape[1]))
        iNumDimensions = llfTransform.shape[0]
        # NOTE: copy is done to convert 2-D ndarray to two 1-D arrays
        pA = (c_double*iNumDimensions)()
        pB = (c_double*iNumDimensions)()
        for iI in range(iNumDimensions):
            pA[iI] = llfTransform[iI][0]
            pB[iI] = llfTransform[iI][1]
        pLibTSG.tsgSetDomainTransform(self.pGrid, pA, pB)

    def isSetDomainTransfrom(self):
        '''
        returns True if the grid is defined for non-canonical domain
        returns False if using a canonical domain

        '''
        return (pLibTSG.tsgIsSetDomainTransfrom(self.pGrid) != 0)

    def clearDomainTransform(self):
        '''
        resets the domain to canonical
        loaded values will be kept, however, the values now correspond
        to canonical points and may be invalid for your application

        '''
        pLibTSG.tsgClearDomainTransform(self.pGrid)

    def getDomainTransform(self):
        '''
        returns llfTransform from the call to setDomainTransform()

        if setDomainTransform() has not been called or if the
        transformed has been cleared by clearDomainTransform(),
        then this returns an empty matrix

        '''
        if (not self.isSetDomainTransfrom()):
            return np.empty([0,2], np.float64)

        iNumDimensions = self.getNumDimensions()
        pA = (c_double*iNumDimensions)()
        pB = (c_double*iNumDimensions)()
        pLibTSG.tsgGetDomainTransform(self.pGrid, pA, pB)
        llfTransform = np.empty([iNumDimensions, 2], np.float64)
        for iI in range(iNumDimensions):
            llfTransform[iI][0] = pA[iI]
            llfTransform[iI][1] = pB[iI]
        return llfTransform

    def setConformalTransformASIN(self, liTruncation):
        '''
        sets conformal domain transform based on truncated
        Maclaurin series of arcsin()

        liTruncation: 1-D numpy.ndarray of non-negative integers
                      indicating the truncation order in each direction
                      0 indicates no transform applied to this direction

        '''
        lShape = liTruncation.shape
        if (len(lShape) != 1):
            raise TasmanianInputError("liTruncation", "ERROR: liTruncation should be a 1-D numpy.ndarray")
        if (lShape[0] != self.getNumDimensions()):
            raise TasmanianInputError("liTruncation", "ERROR: the length of liTruncation is {0:1d} and it should match iDimension: {1:1d}".format(lShape[0], self.getNumDimensions()))
        iNumDimensions = lShape[0]
        pTruncation = (c_int*iNumDimensions)()
        for iI in range(iNumDimensions):
            pTruncation[iI] = liTruncation[iI] # this converts Python longs to c_int
        pLibTSG.tsgSetConformalTransformASIN(self.pGrid, pTruncation)

    def isSetConformalTransformASIN(self):
        '''
        returns True if conformal transform is set
        returns False otherwise

        see: setConformalTransformASIN()

        '''
        return (pLibTSG.tsgIsSetConformalTransformASIN(self.pGrid) != 0)

    def clearConformalTransform(self):
        '''
        resets the conformal domain transform
        loaded values will be kept, however, the values now correspond
        to canonical points and may be invalid for your application

        '''
        pLibTSG.tsgClearConformalTransform(self.pGrid)

    def getConformalTransformASIN(self):
        '''
        returns liTruncation from the call to setConformalTransformASIN()

        if setConformalTransformASIN() has not been called or if the
        transformed has been cleared by clearConformalTransform(),
        then this returns an empty matrix

        '''
        if (not self.isSetConformalTransformASIN()):
            return np.empty([0,], np.int)

        iNumDimensions = self.getNumDimensions()
        pTruncation = (c_int*iNumDimensions)()
        pLibTSG.tsgGetConformalTransformASIN(self.pGrid, pTruncation)
        liTruncation = np.empty([iNumDimensions,], np.int)
        for iI in range(iNumDimensions):
            liTruncation[iI] = pTruncation[iI] # convert c_int to python long
        return liTruncation

    def clearLevelLimits(self):
        '''
        clears the limits set by the last make***Grid or refine command
        if no limits are set, this has no effect
        '''
        pLibTSG.tsgClearLevelLimits(self.pGrid)

    def getLevelLimits(self):
        '''
        returns the limits set by the last call to make***Grid or refine
        returns a vector of integers corresponding to the limits for each
        direction, -1 indicates no limit
        '''
        iNumDimensions = self.getNumDimensions()
        pTruncation = (c_int*iNumDimensions)()
        pLibTSG.tsgGetLevelLimits(self.pGrid, pTruncation)
        liLimits = np.empty([iNumDimensions,], np.int)
        for iI in range(iNumDimensions):
            liLimits[iI] = pTruncation[iI] # convert c_int to python long
        return liLimits

    def setAnisotropicRefinement(self, sType, iMinGrowth, iOutput, liLevelLimits = []):
        '''
        estimates anisotropic coefficients from the current set of
        loaded points and updates the grid with the best points
        according to the estimate

        sType: string identifying the estimate to use (see the Manual)
               recommended: 'iptotal'   'ipcurved'

        iMinGrowth: int (positive)
                minimum number of new points to include in the new grid

        iOutput: int (indicates the output to use)
             selects which output to use for refinement
             sequence grids accept -1 to indicate all outputs

        '''
        if (self.getNumOutputs() == 0):
            raise TasmanianInputError("setAnisotropicRefinement", "ERROR: cannot set refinement for grid with iOutput = 0")
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("setAnisotropicRefinement", "ERROR: cannot call setAnisotropicRefinement for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (iMinGrowth <= 0):
            raise TasmanianInputError("iMinGrowth", "ERROR: the number of growth should be positive integer")
        if (iOutput == -1):
            if (not self.isSequence()):
                raise TasmanianInputError("iOutput", "ERROR: iOutput = -1 can be used only for sequence grids")
        if (iOutput < -1):
            raise TasmanianInputError("iOutput", "ERROR: iOutput should be -1 or a non-negative integer")
        if (iOutput >= self.getNumOutputs()):
            raise TasmanianInputError("iOutput", "ERROR: iOutput cannot exceed the index of the last output {0:1d}".format(self.getNumOutputs() - 1))
        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            iDimension = self.getNumDimensions()
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')
        pLibTSG.tsgSetAnisotropicRefinement(self.pGrid, c_char_p(sType), iMinGrowth, iOutput, pLevelLimits)

    def estimateAnisotropicCoefficients(self, sType, iOutput):
        '''
        returns the estimate of the anisotropic coefficients from the
        current set of loaded points
        see the manual

        sType: string identifying the estimate to use (see the Manual)
               recommended: 'iptotal'   'ipcurved'


        iOutput: int (indicates the output to use)
             selects which output to use for refinement
             sequence grids accept -1 to indicate all outputs

        outputs: 1-D numpy.ndarray
                 of length getNumDimensions() or 2*getNumDimensions()
                 the first set of getNumDimensions() entries correspond
                                                 to the xi coefficients
                 the second set of getNumDimensions() entries correspond
                                                 to the eta coefficients

        '''
        if (self.getNumOutputs() == 0):
            raise TasmanianInputError("estimateAnisotropicCoefficients", "ERROR: cannot set refinement for grid with iOutput = 0")
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("estimateAnisotropicCoefficients", "ERROR: cannot call estimateAnisotropicCoefficients for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (iOutput == -1):
            if (not self.isSequence()):
                raise TasmanianInputError("iOutput", "ERROR: iOutput = -1 can be used only for sequence grids")
        if (iOutput < -1):
            raise TasmanianInputError("iOutput", "ERROR: iOutput should be -1 or a non-negative integer")
        if (iOutput >= self.getNumOutputs()):
            raise TasmanianInputError("iOutput", "ERROR: iOutput cannot exceed the index of the last output {0:1d}".format(self.getNumOutputs() - 1))
        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")

        iNumCoeffs = self.getNumDimensions()
        if ("curved" in sType):
            iNumCoeffs = iNumCoeffs * 2
        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')

        aCoeff = np.empty([iNumCoeffs], np.int32)
        pLibTSG.tsgEstimateAnisotropicCoefficientsStatic(self.pGrid, c_char_p(sType), iOutput, np.ctypeslib.as_ctypes(aCoeff))

        return aCoeff

    def setSurplusRefinement(self, fTolerance, iOutput, sCriteria = "", liLevelLimits = [], llfScaleCorrection = []):
        '''
        using hierarchical surplusses as an error indicator, the surplus
        refinement adds points to the grid to improve accuracy

        when using sequence grids: this algorithm corresponds to the
                                   greedy Knapsack problem

        when using local polynomial or wavelet grids, this call
                                 corresponds to local spatial refinement

        fTolerance: float (non-negative)
                    the relative error tolerance, i.e.,
                    we refine only for points associated with surplus
                    that exceeds the tolerance

        iOutput: int (indicates the output to use)
                 selects which output to use for refinement
                 sequence and local polynomial grids accept -1 to
                 indicate all outputs

        sCriteria: hierarchical and direction refinement strategy
                   'classic'  'parents'   'direction'   'fds'   'stable'
                  applicable only for Local Polynomial and Wavelet grids

        llfScaleCorrection: 2-D numpy.ndarray of non-negative numbers
                            Instead of comparing the normalized surpluses to
                            the tolerance, the scaled surplus will be used.
                            The correction allows to manually guide the
                            refinement process.
                            The surplus of the j-th output of the i-th point
                            will be scaled by llfScaleCorrection[i][j].
                            llfScaleCorrection.shape[0] must be equal to
                            getNumLoaded()
                            If empty, the scale is assumed 1.0
                            llfScaleCorrection.shape[1] must be equal to the
                            number of outputs used in the process, which is
                            equal to getNumOutputs() for iOutput == -1,
                            or 1 if iOutput > -1.
        '''
        if (self.isGlobal()):
            if (self.getRule() not in lsTsgSequenceRules):
                raise TasmanianInputError("setSurplusRefinement", "ERROR: setSurplusRefinement cannot be used with global grids with non-sequence rule")
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("setSurplusRefinement", "ERROR: cannot call setSurplusRefinement for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (fTolerance < 0.0):
            raise TasmanianInputError("fTolerance", "ERROR: fTolerance must be non-negative")
        iActiveOutputs = self.getNumOutputs()
        pScaleCorrection = None
        if (len(llfScaleCorrection) > 0):
            if (iOutput > -1):
                iActiveOutputs = 1
            if (len(llfScaleCorrection.shape) != 2):
                raise TasmanianInputError("llfScaleCorrection", "ERROR: llfScaleCorrection must be a 2-D numpy.ndarray, instead it has {0:1d} dimensions".format(len(llfScaleCorrection.shape)))
            if (llfScaleCorrection.shape[0] != self.getNumLoaded()):
                raise TasmanianInputError("llfScaleCorrection", "ERROR: leading dimension of llfScaleCorrection is {0:1d} but the number of current points is {1:1d}".format(llfScaleCorrection.shape[0], self.getNumLoaded()))
            if (llfScaleCorrection.shape[1] != iActiveOutputs):
                raise TasmanianInputError("llfScaleCorrection", "ERROR: second dimension of llfScaleCorrection is {0:1d} but the refinement is set to use {1:1d}".format(llfScaleCorrection.shape[0], iActiveOutputs))
            pScaleCorrection = np.ctypeslib.as_ctypes(llfScaleCorrection.reshape([self.getNumLoaded() * iActiveOutputs]))

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            iDimension = self.getNumDimensions()
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")
            pLevelLimits = (c_int*iDimension)()
            for iI in range(iDimension):
                pLevelLimits[iI] = liLevelLimits[iI]

        if (len(sCriteria) == 0):
            if (not self.isSequence() and not self.isGlobal()):
                raise TasmanianInputError("sCriteria", "ERROR: sCriteria must be specified")
            pLibTSG.tsgSetGlobalSurplusRefinement(self.pGrid, c_double(fTolerance), iOutput, pLevelLimits)
        else:
            if (self.isSequence()):
                raise TasmanianInputError("sCriteria", "ERROR: sCriteria cannot be used for sequence grids")
            if (not sCriteria in lsTsgRefineTypes):
                raise TasmanianInputError("sCriteria", "ERROR: invalid criteria, see TasmanianSG.lsTsgRefineTypes for the list of accepted types")
            if (sys.version_info.major == 3):
                sCriteria = bytes(sCriteria, encoding='utf8')
            pLibTSG.tsgSetLocalSurplusRefinement(self.pGrid, c_double(fTolerance), c_char_p(sCriteria), iOutput, pLevelLimits, pScaleCorrection)

    def clearRefinement(self):
        '''
        clear the last call to set***Refinement,
        only works if called before the points are loaded, i.e.,
        before loadNeededPoints()

        if getNumNeeded() == 0, this call will have no effect

        '''
        pLibTSG.tsgClearRefinement(self.pGrid)

    def mergeRefinement(self):
        '''
        combines the loaded and needed points into a single grid
        it also invalidates any currently loaded values, i.e., the
        grid cannot be used for internal integration or interpolation
        until loadNeededPoints() or setHierarchicalCoefficients()
        is called (even if those have been called before)

        if getNumNeeded() == 0, this call will have no effect

        '''
        pLibTSG.tsgMergeRefinement(self.pGrid)

    def beginConstruction(self):
        '''
        start dynamic construction procedure
        '''
        pLibTSG.tsgBeginConstruction(self.pGrid)

    def isUsingConstruction(self):
        '''
        check if using dynamic construction
        '''
        return (pLibTSG.tsgIsUsingConstruction(self.pGrid) != 0)

    def getCandidateConstructionPoints(self, sType, liAnisotropicWeightsOrOutput, liLevelLimits = []):
        '''
        returns the sorted points for the construction
        '''
        if (not self.isUsingConstruction()):
            raise TasmanianInputError("getCandidateConstructionPoints", "ERROR: calling getCandidateConstructionPoints() before beginConstruction()")
        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")
        iNumDims = self.getNumDimensions()
        pAnisoWeights = None
        iOutput = -1

        if (((sys.version_info.major == 3) and isinstance(liAnisotropicWeightsOrOutput, int))
            or ((sys.version_info.major == 2) and isinstance(liAnisotropicWeightsOrOutput, (int, long)))):
            iOutput = liAnisotropicWeightsOrOutput
        elif (isinstance(liAnisotropicWeightsOrOutput, (list, np.ndarray))):
            if (len(liAnisotropicWeightsOrOutput) > 0):
                if (sType in lsTsgCurvedTypes):
                    iNumWeights = 2*iNumDims
                else:
                    iNumWeights = iNumDims
                if (len(liAnisotropicWeightsOrOutput) != iNumWeights):
                    raise TasmanianInputError("liAnisotropicWeightsOrOutput", "ERROR: wrong number of liAnisotropicWeights, sType '{0:s}' needs {1:1d} weights but len(liAnisotropicWeights) == {2:1d}".format(sType, iNumWeights, len(liAnisotropicWeightsOrOutput)))
                else:
                    aAWeights = np.array([liAnisotropicWeightsOrOutput[i] for i in range(iNumWeights)], np.int32)
                    pAnisoWeights = np.ctypeslib.as_ctypes(aAWeights)
        else:
            raise TasmanianInputError("liAnisotropicWeightsOrOutput", "ERROR: liAnisotropicWeightsOrOutput should be either an integer or numpy.ndarray")

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iNumDims):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to the grid dimension")
            pLevelLimits = (c_int*iNumDims)()
            for iI in range(iNumDims):
                pLevelLimits[iI] = liLevelLimits[iI]

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')

        pVector = pLibTSG.tsgGetCandidateConstructionPointsVoidPntr(self.pGrid, c_char_p(sType), iOutput, pAnisoWeights, pLevelLimits)

        iNumPoints = pLibTSG.tsgGetCandidateConstructionPointsPythonGetNP(self.pGrid, pVector)
        if (iNumPoints == 0):
            pLibTSG.tsgGetCandidateConstructionPointsPythonDeleteVect(pVector)
            return np.empty([0, 0], np.float64)
        aPoints = np.empty([iNumPoints * iNumDims], np.float64)

        pLibTSG.tsgGetCandidateConstructionPointsPythonStatic(pVector, np.ctypeslib.as_ctypes(aPoints))
        pLibTSG.tsgGetCandidateConstructionPointsPythonDeleteVect(pVector)

        return aPoints.reshape([iNumPoints, iNumDims])

    def getCandidateConstructionPointsSurplus(self, fTolerance, sRefinementType, iOutput = -1, liLevelLimits = [], aScaleCorrection = []):
        '''
        returns the sorted points for the construction
        '''
        if (not self.isUsingConstruction()):
            raise TasmanianInputError("getCandidateConstructionPointsSurplus", "ERROR: calling getCandidateConstructionPointsSurplus() before beginConstruction()")
        iNumDims = self.getNumDimensions()

        if (not sRefinementType in lsTsgRefineTypes):
            raise TasmanianInputError("sRefinementType", "ERROR: calling getCandidateConstructionPointsSurplus() called with incorrect type, see TasmanianSG.lsTsgRefineTypes")
        if (sys.version_info.major == 3):
            sRefinementType = bytes(sRefinementType, encoding='utf8')

        pLevelLimits = None
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iNumDims):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to the grid dimension")
            pLevelLimits = (c_int*iNumDims)()
            for iI in range(iNumDims):
                pLevelLimits[iI] = liLevelLimits[iI]

        pScale = None
        if (len(aScaleCorrection) > 0):
            pScale = np.ctypeslib.as_ctypes(aScaleCorrection.reshape([np.prod(aScaleCorrection.shape),]))

        pVector = pLibTSG.tsgGetCandidateConstructionPointsSurplusVoidPntr(self.pGrid, c_double(fTolerance), c_char_p(sRefinementType), iOutput, pLevelLimits, pScale)

        iNumPoints = pLibTSG.tsgGetCandidateConstructionPointsPythonGetNP(self.pGrid, pVector)
        if (iNumPoints == 0):
            pLibTSG.tsgGetCandidateConstructionPointsPythonDeleteVect(pVector)
            return np.empty([0, 0], np.float64)
        aPoints = np.empty([iNumPoints * iNumDims], np.float64)

        pLibTSG.tsgGetCandidateConstructionPointsPythonStatic(pVector, np.ctypeslib.as_ctypes(aPoints))
        pLibTSG.tsgGetCandidateConstructionPointsPythonDeleteVect(pVector)

        return aPoints.reshape([iNumPoints, iNumDims])


    def loadConstructedPoint(self, lfX, lfY):
        '''
        load the currently computed point or points

        lfX: should be 1D (single point case) or 2D (multi-point case)
             numpy.ndarray
             1D case: the lenght should be iDimensions indicating the
             computed point
             2D case: lfX.shape[1] shold be iDimensions and lfX.shape[0]
             should be the number of computed points
        lfY: should be a 1D or 2D numpy.ndarray matching lfY
             1D case: the length should be iOutputs indicated the
             corresponding model value
             2D case: lfY.shape[1] should be iOutputs and lfY.shape[0]
             must match lfX.shape[0]

        if lfX or lfY are not numpy.ndarrays, an attempt would be made
        to cast them using numpy.array() funciton (e.g., lfX and lfY
        could be 1D lists or tuples).
        '''
        if (not self.isUsingConstruction()):
            raise TasmanianInputError("loadConstructedPoint", "ERROR: calling loadConstructedPoint() before beginConstruction()")
        iNumDims = self.getNumDimensions()
        iNumOuts = self.getNumOutputs()
        if (not isinstance(lfX, np.ndarray)):
            lfX = np.array(lfX)
        if (not isinstance(lfY, np.ndarray)):
            lfY = np.array(lfY)
        if (len(lfX.shape) == 1):
            if (lfX.shape[0] != iNumDims):
                raise TasmanianInputError("lfX", "ERROR: lfX should be 1D numpy.ndarray with length equal to the grid dimension")
            if (len(lfY.shape) != 1 or lfY.shape[0] != iNumOuts):
                raise TasmanianInputError("lfY", "ERROR: lfY should be 1D numpy.ndarray with length equal to the model outputs")
            pLibTSG.tsgLoadConstructedPoint(self.pGrid, np.ctypeslib.as_ctypes(lfX), 1, np.ctypeslib.as_ctypes(lfY))
        elif (len(lfX.shape) == 2):
            iNumX = lfX.shape[0]
            if (iNumX == 0):
                return
            if (len(lfY.shape) != 2):
                raise TasmanianInputError("lfY", "ERROR: if lfX is 2D array, then lfY should be 2D array as well")
            if (lfX.shape[1] != iNumDims):
                raise TasmanianInputError("lfX", "ERROR: lfX should be 2D numpy.ndarray with shape[1] equal to the grid dimension")
            if (lfY.shape[1] != iNumOuts):
                raise TasmanianInputError("lfY", "ERROR: lfY should be 2D numpy.ndarray with shape[1] equal to the model outputs")
            if (lfY.shape[0] != iNumX):
                raise TasmanianInputError("lfY", "ERROR: lfY should provide the same number of entries as lfX, i.e., shape[0] must match")
            pLibTSG.tsgLoadConstructedPoint(self.pGrid, np.ctypeslib.as_ctypes(lfX.reshape(iNumX * iNumDims)), iNumX, np.ctypeslib.as_ctypes(lfY.reshape(iNumX * iNumOuts)))
        else:
            raise TasmanianInputError("lfX", "ERROR: lfX should be 1D or 2D numpy.ndarray")


    def finishConstruction(self):
        '''
        end the dynamic construction procedure
        '''
        pLibTSG.tsgFinishConstruction(self.pGrid)

    def removePointsByHierarchicalCoefficient(self, fTolerance, iOutput = -1, aScaleCorrection = [], iNumKeep = -1):
        '''
        removes any points in the grid with relative surplus that
        exceeds the tolerance or keeps the set number of points
        with largest surplus

        fTolerance: float (positive)
                    the relative surplus tolerance, i.e.,
                    we keep only for points associated with surplus
                    that exceeds the tolerance
                    if iNumKeep is positive, then fTolerance is ignored

        iOutput: int (indicates the output to use)
                 selects which output to consider
                 accept -1 to indicate all outputs

        lfScaleCorrection: numpy array of doubles, either 1D or 2D
                           if iOutputs = -1 and getNumOutputs() > 1,
                           then using 2D array with shape
                           getNumLoaded() X getNumOutputs() with
                           one weight per hierarchical coefficient
                           if iOutputs > -1, then using 1D array with
                           one weight per point

        iNumKeep: int (positive or equal to -1)
                  indicates the number of points to keep
                  if set to -1 then fTolerance is used as a cutoff
                  if positive then the given number of points will be kept
        '''
        if (not self.isLocalPolynomial()):
            raise TasmanianInputError("removePointsByHierarchicalCoefficient", "ERROR: calling removePointsByHierarchicalCoefficient for a grid that isn't local polynomial")
        if (iNumKeep == -1 and fTolerance <= 0.0):
            raise TasmanianInputError("fTolerance", "ERROR: fTolerance must be a positive integer")
        if (iOutput < -1):
            raise TasmanianInputError("iOutput", "ERROR: iOutput should be -1 or a non-negative integer")
        if (iOutput >= self.getNumOutputs()):
            raise TasmanianInputError("iOutput", "ERROR: iOutput cannot exceed the index of the last output {0:1d}".format(self.getNumOutputs() - 1))
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("removePointsByHierarchicalCoefficient", "ERROR: calling removePointsByHierarchicalCoefficient when no points are loades")
        if (iNumKeep == 0 or iNumKeep < -1 or iNumKeep > self.getNumLoaded()):
            raise TasmanianInputError("iNumKeep", "ERROR: iNumKeep should be either -1 or positive without exceeding the number of loaded points.")
        if (len(aScaleCorrection) > 0):
            lShape = aScaleCorrection.shape
            if (len(lShape) != 2):
                raise TasmanianInputError("aScaleCorrection", "ERROR: aScaleCorrection should be a 2D array")
            if (lShape[0] != self.getNumLoaded()):
                raise TasmanianInputError("aScaleCorrection", "ERROR: aScaleCorrection.shape[0] should match getNumLoaded()")
            if (iOutput == -1 and lShape[1] != self.getNumOutputs()):
                raise TasmanianInputError("aScaleCorrection", "ERROR: aScaleCorrection.shape[1] should match getNumOutputs()")
            if (iOutput != -1 and lShape[1] != 1):
                raise TasmanianInputError("aScaleCorrection", "ERROR: aScaleCorrection.shape[1] should be 1")

        if (len(aScaleCorrection) == 0):
            if (iNumKeep == -1):
                pLibTSG.tsgRemovePointsByHierarchicalCoefficient(self.pGrid, fTolerance, iOutput, None)
            else:
                pLibTSG.tsgRemovePointsByHierarchicalCoefficientHardCutoff(self.pGrid, iNumKeep, iOutput, None)
        else:
            lShape = aScaleCorrection.shape
            iNumWeights = lShape[0]
            if (iOutput == -1):
                iNumWeights *= lShape[1]
            if (iNumKeep == -1):
                pLibTSG.tsgRemovePointsByHierarchicalCoefficient(self.pGrid, fTolerance, iOutput,           np.ctypeslib.as_ctypes(aScaleCorrection.reshape([iNumWeights,])))
            else:
                pLibTSG.tsgRemovePointsByHierarchicalCoefficientHardCutoff(self.pGrid, iNumKeep, iOutput,           np.ctypeslib.as_ctypes(aScaleCorrection.reshape([iNumWeights,])))

    def getHierarchicalCoefficients(self):
        '''
        For global grids, this just returns the values loaded using the
        call to loadNeededPoints().
        In all other cases, this returns the list of hierarchical
        coefficients, i.e., surpluses.

        returns a 2-D numpy array getNumPoints() by getNumOutputs()
        '''
        iNumOuts = self.getNumOutputs()
        if (iNumOuts == 0):
            return np.empty([0,0], np.float64)
        iNumPoints = self.getNumLoaded()
        if (iNumPoints == 0):
            return np.empty([0,iNumOuts], np.float64)

        if (not self.isFourier()):
            aSurp = np.empty([iNumOuts * iNumPoints], np.float64)
            pLibTSG.tsgGetHierarchicalCoefficientsStatic(self.pGrid, np.ctypeslib.as_ctypes(aSurp))
        else:
            aSurp = np.empty([2 * iNumOuts * iNumPoints], np.float64)
            pLibTSG.tsgGetHierarchicalCoefficientsStatic(self.pGrid, np.ctypeslib.as_ctypes(aSurp))
            aSurp = aSurp[:(iNumOuts * iNumPoints)] + 1j * aSurp[(iNumOuts * iNumPoints):]

        return aSurp.reshape([iNumPoints, iNumOuts])

    def evaluateHierarchicalFunctions(self, llfX):
        '''
        evaluates the hierarchical functions at a set of points in the
        domain and return a 2-D numpy.ndarray with the result

        llfX: a 2-D numpy.ndarray with llfX.shape[1] == iDimensions
              the entries indicate the points for evaluating the weights

        output: returns a 2-D numpy.ndarray of
                shape == [llfX.shape[0], getNumPoints()]
                the values of the basis functions at the points
        '''
        if (len(llfX.shape) != 2):
            raise TasmanianInputError("llfX", "ERROR: calling evaluateHierarchicalFunctions llfX should be a 2-D numpy array")
        if (llfX.shape[1] != self.getNumDimensions()):
            raise TasmanianInputError("llfX", "ERROR: calling evaluateHierarchicalFunctions llfX.shape[1] is not equal to getNumDimensions()")
        iNumX = llfX.shape[0]
        # see evaluateBatch()
        lfX = llfX.reshape([llfX.shape[0] * llfX.shape[1]])

        if not self.isFourier():
            aResult = np.empty([iNumX * self.getNumPoints()], np.float64)
            pLibTSG.tsgEvaluateHierarchicalFunctions(self.pGrid, np.ctypeslib.as_ctypes(lfX), iNumX, np.ctypeslib.as_ctypes(aResult))
        else:
            aResult = np.empty([2 * iNumX * self.getNumPoints()], np.float64)
            pLibTSG.tsgEvaluateHierarchicalFunctions(self.pGrid, np.ctypeslib.as_ctypes(lfX), iNumX, np.ctypeslib.as_ctypes(aResult))
            aResult = aResult[0::2] + 1j * aResult[1::2]

        return aResult.reshape([iNumX, self.getNumPoints()])

    def getHierarchicalSupport(self):
        '''
        returns the support of the hierarchical basis in a 2-D numpy.ndarray

        - only local-polynomial and wavelet grids have restricted support
        - the support of all basis is restricted to the domain even
          if the support includes additional ares

        output: shape == [getNumPoints(), getNumDimensions()]
                the support entries correspond to the output of getPoints()
                if x represents a point in the domain and
                if abs( x[j] - getPoints()[i, j] ) > getSupport()[i, j]
                then the i-th entry in evaluateHierarchicalFunctions(x)
                corresponding to x is guaranteed to be zero
        '''
        if (self.getNumPoints() == 0):
            return np.empty([0,0], np.float64)
        aResult = np.empty((self.getNumPoints() * self.getNumDimensions()), np.float64)
        pLibTSG.tsgGetHierarchicalSupportStatic(self.pGrid, np.ctypeslib.as_ctypes(aResult))
        return aResult.reshape((self.getNumPoints(), self.getNumDimensions()))

    def evaluateSparseHierarchicalFunctions(self, llfX):
        '''
        evaluates the hierarchical functions at a set of points in the
        domain. The distinction between this function and
        evaluateHierarchicalFunctions() lies in the type of the returned
        result, namely a sparse vs a dense matrix.

        The motivation for this function is that Local Polynomial and
        Wavelet grids usually result in sparse matrices

        llfX: a 2-D numpy.ndarray with llfX.shape[1] == iDimensions
              the entries indicate the points for evaluating

        output: returns a TasmanianSimpleSparseMatrix class
                which is a simple class with three fields:
                aPntr, aIndx, and aVals which are numpy.ndarray of types
                int32, int32, and float64
                iNumRows and iNumCols are meta fields and have values
                iNumRows = llfX.shape[0]
                iNumCols = self.getNumPoints()
                The sparse matrix is compressed along the llfX.shape[0]
                dimension, i.e., using column compressed format
        '''
        if (len(llfX.shape) != 2):
            raise TasmanianInputError("llfX", "ERROR: calling evaluateSparseHierarchicalFunctions(), llfX should be a 2-D numpy array")
        if (llfX.shape[1] != self.getNumDimensions()):
            raise TasmanianInputError("llfX", "ERROR: calling evaluateSparseHierarchicalFunctions(), llfX.shape[1] is not equal to getNumDimensions()")
        iNumX = llfX.shape[0]
        pMat = TasmanianSimpleSparseMatrix()
        iNumNZ = pLibTSG.tsgEvaluateSparseHierarchicalFunctionsGetNZ(self.pGrid, np.ctypeslib.as_ctypes(llfX.reshape([llfX.shape[0] * llfX.shape[1]])), iNumX)
        pMat.aPntr = np.empty([iNumX+1,], np.int32)
        pMat.aIndx = np.empty([iNumNZ,], np.int32)
        pMat.aVals = np.empty([iNumNZ if not self.isFourier() else 2 * iNumNZ,], np.float64)
        pMat.iNumRows = iNumX
        pMat.iNumCols = self.getNumPoints()
        # see evaluateBatch()
        lfX = llfX.reshape([llfX.shape[0] * llfX.shape[1]])
        pLibTSG.tsgEvaluateSparseHierarchicalFunctionsStatic(self.pGrid, np.ctypeslib.as_ctypes(lfX), iNumX,
                                                             np.ctypeslib.as_ctypes(pMat.aPntr), np.ctypeslib.as_ctypes(pMat.aIndx), np.ctypeslib.as_ctypes(pMat.aVals))

        return pMat

    def setHierarchicalCoefficients(self, llfCoefficients):
        '''
        Local polynomial, Wavelet, and Sequence grids construct
        approximation using hierarchical coefficients based on the
        loaded values. This function does the opposite, the hierarchical
        coefficients are loaded directly and the values are computed
        based on the coefficients. The coefficients can be computed,
        e.g., by solving least-squares or compressed sensing problem
                   min || A c - f ||
        where A is a matrix returned by evaluateHierarchicalFunctions()
        or evaluateSparseHierarchicalFunctions() for a set of points
        llfX; f are the values of the target function at the llfX
        points; and c is the vector with corresponding hierarchical
        coefficients.

        If there is a pending refinement, i.e., getNumLoaded() != 0 and
        getNumNeeded() != 0, then the refinement is discarded (since it
        was computed based on the old and now obsolete values)

        llfCoefficients: a 2-D numpy.ndarray
                         with dimensions getNumPoints() X iOutputs
                         each row corresponds to the values of the
                         coefficients at the corresponding point.
                         The order and leading dimension must match the
                         points obtained form getPoints(), the same
                         order as the second dimension of
                         evaluateHierarchicalFunctions()
        '''
        if (len(llfCoefficients.shape) != 2):
            raise TasmanianInputError("llfCoefficients", "ERROR: llfCoefficients should be a 2-D numpy.ndarray, instead it has {0:1d} dimensions".format(len(llfCoefficients.shape)))
        if (llfCoefficients.shape[0] != self.getNumPoints()):
            raise TasmanianInputError("llfCoefficients", "ERROR: leading dimension of llfCoefficients is {0:1d} but the number of current points is {1:1d}".format(llfCoefficients.shape[0], self.getNumNeeded()))
        if (llfCoefficients.shape[1] != self.getNumOutputs()):
            raise TasmanianInputError("llfCoefficients", "ERROR: second dimension of llfCoefficients is {0:1d} but the number of outputs is set to {1:1d}".format(llfCoefficients.shape[1], self.getNumOutputs()))
        if (self.isFourier() and (not np.iscomplexobj(llfCoefficients))):
            raise TasmanianInputError("llfCoefficients", "ERROR: using Fourier grid but llfCoefficients is not complex")

        iNumPoints = llfCoefficients.shape[0]
        iNumDims = llfCoefficients.shape[1]

        if self.isFourier():
            llfCoefficientsTmp = np.vstack((np.real(llfCoefficients), np.imag(llfCoefficients)))
            llfCoefficients = llfCoefficientsTmp.reshape([2 * iNumPoints * iNumDims,])
        else:
            llfCoefficients = llfCoefficients.reshape([iNumPoints * iNumDims,])

        pLibTSG.tsgSetHierarchicalCoefficients(self.pGrid, np.ctypeslib.as_ctypes(llfCoefficients))

    def integrateHierarchicalFunctions(self):
        '''
        Computes the integrals of the hierarchical basis functions,
        i.e., the same functions computed by evaluateHierarchicalFunctions().

        returns a one dimensional np.ndarray with size self.getNumPoints()
        '''
        iNumPoints = self.getNumPoints()
        if (iNumPoints == 0):
            return np.empty([0,], np.float64)

        aIntegrals = np.zeros([iNumPoints,])
        pLibTSG.tsgIntegrateHierarchicalFunctionsStatic(self.pGrid, np.ctypeslib.as_ctypes(aIntegrals))

        return aIntegrals

    def getGlobalPolynomialSpace(self, bInterpolation):
        '''
        returns a matrix corresponding to the polynomial space that is
        integrated or interpolated exactly by the current grid

        bInterpolation: boolean
                indicates whether to give the space associated
                with integration or interpolation

        output: is a 2-D numpy.ndarray of integers
            output.shape[0] indicates the cardinality of the space
            output.shape[1] is equal to iDimension
            each row corresponds to a multi-index associated with a
            polynomial in a hierarchical tensor basis, for example,
            monomials

            see the manual for details

        '''
        iInterp = 0
        if (bInterpolation):
            iInterp = 1
        pNumIndexes = (c_int*1)()
        pIndexes = pLibTSG.tsgPythonGetGlobalPolynomialSpace(self.pGrid, iInterp, pNumIndexes)
        iNumDimensions = self.getNumDimensions()
        lliPolynomials = np.empty([pNumIndexes[0], iNumDimensions], np.int)
        for iI in range(pNumIndexes[0]):
            for iJ in range(iNumDimensions):
                lliPolynomials[iI][iJ] = pIndexes[iI*iNumDimensions + iJ]
        pLibTSG.tsgDeleteInts(pIndexes)
        return lliPolynomials

    def enableAcceleration(self, sAccelerationType, iGPUID = None):
        '''
        Enables the use of accelerated backend libraries and extensions,
        such as BLAS and CUDA.
        Each acceleration type requires corresponding CMake compile
        options, otherwise the backend will fallback to the closest
        available options.

        sAccelerationType: string

          'none'
              core fallback mode, relies on sequential implementation
              if compiled with Tasmanian_ENABLE_OPENMP this will use
              simple "omp parallel for" to take advantage of multiple cpu cores

          'cpu-blas'
              uses BLAS level 2 and 3 functions for acceleration of batch
              evaluations
              requires Tasmanian_ENABLE_BLAS=ON
              this is the default mode, if available

          'gpu-default'
              uses CUDA kernels, cuBlas, cuSparse and MAGMA libraries for
              accelerated matrix operations, e.g., cublasDgemm
              refer to TasGrid::TypeAcceleration for more details

          'gpu_cublas'
              uses the Nvidia cuBlas and cuSparse libraries

          'gpu-cuda'
              uses custom CUDA kernels in addition to the accelerated
              linear algebra libraries

           'gpu-magma'
              uses the custom CUDA kernels and the MAGMA library in place
              of the default Nvidia libraries

        iGPU: integer
              indicates the GPU device to use, if set to None then device
              zero will be used first or the device set with setGPUID()
        '''
        if (sAccelerationType not in lsTsgAccelTypes):
            raise TasmanianInputError("sAccelerationType", "ERROR: invalid acceleration type")
        if (sys.version_info.major == 3):
            sAccelerationType = bytes(sAccelerationType, encoding='utf8')
        if (iGPUID is None):
            pLibTSG.tsgEnableAcceleration(self.pGrid, c_char_p(sAccelerationType))
        else:
            if ((iGPUID < 0) or (iGPUID >= self.getNumGPUs())):
                raise TasmanianInputError("iGPUID", "ERROR: invalid GPU ID number")
            pLibTSG.tsgEnableAcceleration(self.pGrid, c_char_p(sAccelerationType), iGPUID)

    def getAccelerationType(self):
        '''
        returns the type of acceleration set by enableAcceleration
        '''
        sAccType = pLibTSG.tsgGetAccelerationType(self.pGrid)
        if (sys.version_info.major == 3):
            sAccType = str(sAccType, encoding='utf8')
        return sAccType

    def isAccelerationAvailable(self, sAccelerationType):
        '''
           returns True if the library has been compiled with support
           for sAccelerationType.
           Even if this returns False, you can use the type for
           enableAcceleration, but the library will default to the next
           available type (see the Manual)
        '''
        if (sAccelerationType not in lsTsgAccelTypes):
            raise TasmanianInputError("sAccelerationType", "ERROR: invalid acceleration type")
        if (sys.version_info.major == 3):
            sAccelerationType = bytes(sAccelerationType, encoding='utf8')
        return (pLibTSG.tsgIsAccelerationAvailable(sAccelerationType) != 0)

    def setGPUID(self, iGPUID):
        '''
        when using cuda on a machine with multiple GPUs, this helps set
        the GPU for this grid
        NOTE: each instance of the sparse grids class holds a separate
              instance of iGPUID and different grids can be assigned to
              different GPUs (on multi-gpu system)
        iGPUID can be changed at any time, however, this will cause
        some of the internal cache to be invalidated and it may lead
        to extraneous data movement

        calling read or make***Grid will reset the selected GPU

        defaults to 0

        this doesn't do anything unless enableAcceleration is called
        using a "gpu-" acceleration type
        '''
        if ((iGPUID < 0) or (iGPUID >= self.getNumGPUs())):
            raise TasmanianInputError("iGPUID", "ERROR: invalid GPU ID number")
        pLibTSG.tsgSetGPUID(self.pGrid, iGPUID)

    def getGPUID(self):
        '''
        returns the GPU ID set using setGPUID
        '''
        return pLibTSG.tsgGetGPUID(self.pGrid)

    def getNumGPUs(self):
        '''
        returns the number of available GPUs according to cuda

        this is one of several functions designed to allow basic
        management of multi-gpu setup with only Tasmanian module
        '''
        return pLibTSG.tsgGetNumGPUs()

    def getGPUMemory(self, iGPUID):
        '''
        returns the total memory (in MegaBytes, 1024**2 bytes) of the
        corresponding GPU

        this is one of several functions designed to allow basic
        management of multi-gpu setup with only Tasmanian module
        '''
        if ((iGPUID < 0) or (iGPUID >= self.getNumGPUs())):
            raise TasmanianInputError("iGPUID", "ERROR: invalid GPU ID number")
        return pLibTSG.tsgGetGPUMemory(iGPUID)

    def getGPUName(self, iGPUID):
        '''
        return the cuda name ID of the corresponding GPU

        this is one of several functions designed to allow basic
        management of multi-gpu setup with only Tasmanian module
        '''
        if ((iGPUID < 0) or (iGPUID >= self.getNumGPUs())):
            raise TasmanianInputError("iGPUID", "ERROR: invalid GPU ID number")
        pName = create_string_buffer(256)
        iNumChars = np.array([0], np.int32)
        pLibTSG.tsgGetGPUName(iGPUID, 256, pName, np.ctypeslib.as_ctypes(iNumChars))
        return self.stringBufferToString(pName, iNumChars[0])

    def printStats(self):
        '''
        calls the library printStats() function, which displays basic
        information about this instance of the grid
        '''
        pLibTSG.tsgPrintStats(self.pGrid)

    def plotPoints2D(self, pAxisObject=tsgPlot, sStyle="bo", iMarkerSize=3):
        '''
        plots the points in a 2D plot using matplotlib.pyplot
        applicable only for grids with iDimensions == 2

        pAxisObject: axis object from the matplotlib.pyplot package

        sStyle: string
                the matplotlib.pyplot style, e.g.,
                'ko' will make black cirlces, 'rx' will use red crosses

        iMarkerSize: positive integer
                     the marker size for plotting the points
        '''
        if (not bTsgPlotting):
            raise TasmanianInputError("plotPoints2D", "ERROR: could not load matplotlib.pyplot")

        if (self.getNumDimensions() != 2):
            raise TasmanianInputError("plotPoints2D", "ERROR: cannot plot a grid with other than 2 dimensions")

        aPoints = self.getPoints()

        fXmin = min(aPoints[:,0])
        fXmax = max(aPoints[:,0])
        fYmin = min(aPoints[:,1])
        fYmax = max(aPoints[:,1])
        if (fXmin == fXmax):
            fXmin = fXmin - 0.1
            fXmax = fXmax + 0.1
        if (fYmin == fYmax):
            fYmin = fYmin - 0.1
            fYmax = fYmax + 0.1

        pAxisObject.plot(aPoints[:,0], aPoints[:,1], sStyle, markersize=iMarkerSize)
        pAxisObject.axis([fXmin - 0.1 * np.fabs(fXmin), fXmax + 0.1 * np.fabs(fYmax), fYmin - 0.1 * np.fabs(fYmin), fYmax + 0.1 * np.fabs(fYmax)])

    def plotResponse2D(self, iOutput=0, iNumDim0=100, iNumDim1=100, pAxisObject=tsgPlot, sCmap="jet"):
        '''
        plots the response in a 2D plot using matplotlib.pyplot
        applicable only for grids with iDimensions == 2

        iOutput is the output to use for plotting

        iNumDim0, iNumDim1: positive integers
               the points for the plot are selected on a dense grid with
               number of points iNumDim0 and iNumDim1 in dimensions
               0 and 1 respectively

        pAxisObject: axis object from the matplotlib.pyplot package

        sCmap: string indicating the map to use, e.g., "jet" or "heat"
        '''
        if (not bTsgPlotting):
            raise TasmanianInputError("plotResponse2D", "ERROR: could not load matplotlib.pyplot")
        if (iOutput < 0):
            raise TasmanianInputError("iOutput", "ERROR: iOutput should be a non-negative integer")
        if (iOutput >= self.getNumOutputs()):
            raise TasmanianInputError("iOutput", "ERROR: iOutput cannot exceed the index of the last output {0:1d}".format(self.getNumOutputs() - 1))
        if (self.getNumDimensions() != 2):
            raise TasmanianInputError("plotResponse2D", "ERROR: cannot plot a grid with other than 2 dimensions")
        if (iNumDim0 < 1):
            raise TasmanianInputError("iNumDim0", "ERROR: the number of points should be at least 1")
        if (iNumDim1 < 1):
            raise TasmanianInputError("iNumDim1", "ERROR: the number of points should be at least 1")

        aPoints = self.getPoints()

        fXmin = min(aPoints[:,0])
        fXmax = max(aPoints[:,0])
        fYmin = min(aPoints[:,1])
        fYmax = max(aPoints[:,1])
        if (fXmin == fXmax):
            fXmin = fXmin - 0.1
            fXmax = fXmax + 0.1
        if (fYmin == fYmax):
            fYmin = fYmin - 0.1
            fYmax = fYmax + 0.1

        x = np.linspace(fXmin, fXmax, iNumDim0)
        y = np.linspace(fYmax, fYmin, iNumDim1) # flip the order of y to match the top-to-bottom pixel indexing

        XX, YY = np.meshgrid(x, y)
        ZZ = self.evaluateBatch(np.vstack((XX.reshape((iNumDim0*iNumDim1,)), YY.reshape((iNumDim0*iNumDim1,)))).T)
        ZZ = ZZ[:,iOutput].reshape((iNumDim0,iNumDim1))

        pAxisObject.imshow(ZZ, cmap=sCmap, extent=[fXmin, fXmax, fYmin, fYmax])

    def __testMakeGlobalGrid(self, iDimension, iOutputs, iDepth, sType, liAnisotropicWeights, liLevelLimits):
        '''
        Tests the correctness of makeGlobalGrid() and its variants.
        '''
        if (iDimension <= 0):
            raise TasmanianInputError("iDimension", "ERROR: dimension should be a positive integer")
        if (iOutputs < 0):
            raise TasmanianInputError("iOutputs", "ERROR: outputs should be a non-negative integer")
        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")
        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")
        if (len(liAnisotropicWeights) > 0):
            iNumWeights = 2*iDimension if sType in lsTsgCurvedTypes else iDimension
            if (len(liAnisotropicWeights) != iNumWeights):
                raise TasmanianInputError("liAnisotropicWeights", "ERROR: wrong number of liAnisotropicWeights, sType '{0:s}' needs {1:1d} "
                                          "weights but len(liAnisotropicWeights) == {2:1d}".format(sType, iNumWeights, len(liAnisotropicWeights)))
        if (len(liLevelLimits) > 0):
            if (len(liLevelLimits) != iDimension):
                raise TasmanianInputError("liLevelLimits", "ERROR: invalid number of level limits, must be equal to iDimension")

def makeGlobalGrid(iDimension, iOutputs, iDepth, sType, sRule, liAnisotropicWeights=[], fAlpha=0.0, fBeta=0.0, sCustomFilename="", liLevelLimits=[]):
    '''
    Factory method equivalent to TasmanianSparseGrid.makeGlobalGrid().
    '''
    grid = TasmanianSparseGrid()
    grid.makeGlobalGrid(iDimension, iOutputs, iDepth, sType, sRule, liAnisotropicWeights, fAlpha, fBeta, sCustomFilename, liLevelLimits)
    return grid

def makeSequenceGrid(iDimension, iOutputs, iDepth, sType, sRule, liAnisotropicWeights=[], liLevelLimits=[]):
    '''
    Factory method equivalent to TasmanianSparseGrid.makeSequenceGrid().
    '''
    grid = TasmanianSparseGrid()
    grid.makeSequenceGrid(iDimension, iOutputs, iDepth, sType, sRule, liAnisotropicWeights, liLevelLimits)
    return grid

def makeFourierGrid(iDimension, iOutputs, iDepth, sType, liAnisotropicWeights=[], liLevelLimits=[]):
    '''
    Factory method equivalent to TasmanianSparseGrid.makeFourierGrid().
    '''
    grid = TasmanianSparseGrid()
    grid.makeFourierGrid(iDimension, iOutputs, iDepth, sType, liAnisotropicWeights, liLevelLimits)
    return grid

def makeLocalPolynomialGrid(iDimension, iOutputs, iDepth, iOrder=1, sRule="localp", liLevelLimits=[]):
    '''
    Factory method equivalent to TasmanianSparseGrid.makeLocalPolynomialGrid().
    '''
    grid = TasmanianSparseGrid()
    grid.makeLocalPolynomialGrid(iDimension, iOutputs, iDepth, iOrder, sRule, liLevelLimits)
    return grid

def makeWaveletGrid(iDimension, iOutputs, iDepth, iOrder=1, liLevelLimits=[]):
    '''
    Factory method equivalent to TasmanianSparseGrid.makeWaveletGrid().
    '''
    grid = TasmanianSparseGrid()
    grid.makeWaveletGrid(iDimension, iOutputs, iDepth, iOrder, liLevelLimits)
    return grid

def makeGlobalGridCustom(iDimension, iOutputs, iDepth, sType, pCustomTabulated, liAnisotropicWeights=[], liLevelLimits=[]):
    '''
    Factory method equivalent to TasmanianSparseGrid.makeGlobalGridCustom().
    '''
    grid = TasmanianSparseGrid()
    grid.makeGlobalGridCustom(iDimension, iOutputs, iDepth, sType, pCustomTabulated, liAnisotropicWeights, liLevelLimits)
    return grid

def copyGrid(source, iOutputsBegin = 0, iOutputsEnd = -1):
    '''
    Factory method equivalent to TasmanianSparseGrid.copyGrid().
    '''
    grid = TasmanianSparseGrid()
    grid.copyGrid(source, iOutputsBegin, iOutputsEnd)
    return grid

def check_np_arr(np_data_name, np_data, expected_length, expected_dimension):
    '''
    Utility function that checks if a NumPy array conforms to given input lengths and dimensions.
    '''
    if len(np_data.shape) != expected_dimension:
        raise TasmanianInputError(np_data_name, "ERROR: the dimension of "  + np_data_name +
                                  " does not match the expected dimension of " + str(expected_dimension))
    if len(np_data) != expected_length:
        raise TasmanianInputError(np_data_name, "ERROR: the length of "  + np_data_name +
                                  " does not match the expected length of " + str(expected_length))

class CustomTabulated:
    def __init__(self):
        '''
        Constructor that creates an empty CustomTabulated instance.
        '''
        self.CustomTabulatedObject = True
        self.pCustomTabulated = pLibTSG.tsgConstructCustomTabulated()

    def __del__(self):
        '''
        Destructor that calls the C++ destructor and releases all memory used by this instance of the class.
        Make sure to call this to avoid memory leaks.
        '''
        pLibTSG.tsgDestructCustomTabulated(self.pCustomTabulated)

    def read(self, sFilename):
        '''
        Reads the CustomTabulated object from a file and discards any existing grid held by this class.

        sFilename: string indicating a CustomTabulated that was already written using write from Python or any other
                   Tasmanian interfaces.

        output: boolean
            True: the read was successful
            False: the read failed, check the CLI output for an error message
        '''
        effective_filename = bytes(sFilename, encoding='utf8') if sys.version_info.major == 3 else sFilename
        if pLibTSG.tsgReadCustomTabulated(self.pCustomTabulated, c_char_p(effective_filename)) == 0:
            raise TasmanianInputError("sFilename", "ERROR: {0:1s} does not appear to be a valid Tasmanian file.".format(sFilename))

    def write(self, sFilename):
        '''
        Writes the CustomTabulated object to a file in ASCII format.

        sFilename: string indicating a location where the CustomTabulated instance will be written to.
        '''
        effective_filename = bytes(sFilename, encoding='utf8') if sys.version_info.major == 3 else sFilename
        pLibTSG.tsgWriteCustomTabulated(self.pCustomTabulated, c_char_p(effective_filename))

    def getNumLevels(self):
        return pLibTSG.tsgGetNumLevelsCustomTabulated(self.pCustomTabulated)

    def getNumPoints(self, level):
        return pLibTSG.tsgGetNumPointsCustomTabulated(self.pCustomTabulated, c_int(level))

    def getIExact(self, level):
        return pLibTSG.tsgGetIExactCustomTabulated(self.pCustomTabulated, c_int(level))

    def getQExact(self, level):
        return pLibTSG.tsgGetQExactCustomTabulated(self.pCustomTabulated, c_int(level))

    def getDescription(self):
        return pLibTSG.tsgGetDescriptionCustomTabulated(self.pCustomTabulated).decode('UTF-8')

    def getWeightsNodes(self, level):
        '''
        Outputs the weights and nodes (in that order) of the CustomTabulated instance for a given level.

        level: int indicating the level where the weights/nodes will be pulled from.

        output: int (weights), int (nodes)
        '''
        iNumPoints = self.getNumPoints(level)
        if (iNumPoints == 0):
            return np.empty([0], np.float64), np.empty([0], np.float64)
        pWeights, pNodes = (iNumPoints * c_double)(), (iNumPoints * c_double)()
        pLibTSG.tsgGetWeightsNodesStaticCustomTabulated(self.pCustomTabulated, c_int(level), pWeights, pNodes)
        return np.array(pWeights), np.array(pNodes)

def makeCustomTabulatedFromFile(filename):
    '''
    Wrapper that calls the CustomTabulated constructor which reads its data from a file.

    filename: string indicating the location of the file.

    output: CustomTabulated
    '''
    ct = CustomTabulated()
    ct.read(filename)
    return ct

def makeCustomTabulatedFromData(num_levels, num_nodes, precision, nodes, weights, description):
    '''
    Wrapper that calls the CustomTabulated constructor which takes ownership of its input data.

    num_levels: int indicating the number of levels of the instance.
    num_nodes: 1D NumPy array of length num_levels whose i-th entry is the number of nodes at level i.
    precision: 1D NumPy array of length num_levels whose i-th entry is the (quadrature/integration) precision at level i.
    nodes: Python list of length num_levels whose i-th entry is a 1D NumPy array containing the nodes at level i.
    weights: Python list of length num_levels whose i-th entry is a 1D NumPy array containing the quadrature weights at level i.
    description: string that briefly describes the instance.

    output: CustomTabulated
    '''
    # nodes and weights are expected to be Python lists.
    if len(nodes) != num_levels:
        raise TasmanianInputError("nodes", "ERROR: the length of nodes does not match the expected length " + str(num_levels))
    if len(weights) != num_levels:
        raise TasmanianInputError("weights", "ERROR: the length of weights does not match the expected length " + str(num_levels))
    check_np_arr("num_nodes", num_nodes, num_levels, 1);
    check_np_arr("precision", precision, num_levels, 1);
    for i in range(num_levels):
        check_np_arr("nodes["+str(i)+"]", nodes[i], num_nodes[i], 1);
        check_np_arr("weights["+str(i)+"]", weights[i], num_nodes[i], 1);
    ct = CustomTabulated()
    effective_description = bytes(description, encoding='utf8') if sys.version_info.major == 3 else description
    # create the C arrays for num_nodes, precision, nodes, and weights by copying.
    pNumNodes, pPrecision, pNodes, pWeights = None, None, None, None
    if (num_levels > 0):
        pNumNodes = (c_int * num_levels)()
        pPrecision = (c_int * num_levels)()
        for iI in range(num_levels):
            pNumNodes[iI] = num_nodes[iI]
            pPrecision[iI] = precision[iI]
        total_num_nodes = np.sum(num_nodes)
        pNodes = (c_double * total_num_nodes)()
        node_idx = 0
        for lfJ in nodes:
            for fK in lfJ:
                pNodes[node_idx] = fK
                node_idx += 1
        pWeights = (c_double * total_num_nodes)()
        weight_idx = 0
        for lfJ in weights:
            for fK in lfJ:
                pWeights[weight_idx] = fK
                weight_idx += 1
    ct.pCustomTabulated = pLibTSG.tsgMakeCustomTabulatedFromData(c_int(num_levels), pNumNodes, pPrecision, pNodes, pWeights,
                                                                 c_char_p(effective_description))
    return ct
