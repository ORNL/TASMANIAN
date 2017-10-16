#!/usr/bin/env python

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

from ctypes import c_char_p, c_int, c_double, c_void_p, POINTER, cdll
import numpy as np
import sys

bTsgPlotting = True
try:
    import matplotlib.pyplot as tsgPlot
except:
    bTsgPlotting = False

__version__ = "5.1"
__license__ = "BSD 3-Clause with UT-Battelle disclamer"
__author__ = "Miroslav Stoyanov"


lsTsgGlobalRules = ["clenshaw-curtis", "clenshaw-curtis-zero", "chebyshev", "chebyshev-odd", "gauss-legendre", "gauss-legendre-odd", "gauss-patterson", "leja", "leja-odd", "rleja", "rleja-odd", "rleja-double2", "rleja-double4", "rleja-shifted", "rleja-shifted-even", "rleja-shifted-double", "max-lebesgue", "max-lebesgue-odd", "min-lebesgue", "min-lebesgue-odd", "min-delta", "min-delta-odd", "gauss-chebyshev1", "gauss-chebyshev1-odd", "gauss-chebyshev2", "gauss-chebyshev2-odd", "fejer2", "gauss-gegenbauer", "gauss-gegenbauer-odd", "gauss-jacobi", "gauss-jacobi-odd", "gauss-laguerre", "gauss-laguerre-odd", "gauss-hermite", "gauss-hermite-odd", "custom-tabulated"]
lsTsgGlobalTypes = ["level", "curved", "iptotal", "ipcurved", "qptotal", "qpcurved", "hyperbolic", "iphyperbolic", "qphyperbolic", "tensor", "iptensor", "qptensor"]
lsTsgCurvedTypes = ["curved", "ipcurved", "qpcurved"]
lsTsgSequenceRules = ["leja", "rleja", "rleja-shifted", "max-lebesgue", "min-lebesgue", "min-delta"]
lsTsgLocalRules = ["localp", "semi-localp", "localp-zero"]
lsTsgAccelTypes = ["none", "cpu-blas", "gpu-fullmem"]

class TasmanianInputError(Exception):
    '''Exception raised for incorret input to Tasmanian

    Attributes:
    sVariable -- string containing the variable name with incorrect value
    sMessage -- message regarding the error

    '''
    def __init__(self, sVar, sMess):
        self.sVariable = sVar
        self.sMessage = sMess        

    def printInfo(self):
        '''
        prints information for the incorect function or variable

        '''
        print("Incorrect input for: {0:s}".format(tsgError.sVariable))
        print(tsgError.sMessage)


class TasmanianSparseGrid:
    def __init__(self, tasmanian_library=0):
        '''
        constructor, creates an empty grid instance

        tasmanian_library: indicates the libtasmaniansparsegrids.so file

        if tasmanian_library is an int: look for the library in the
                                        default locaiton
                                        "libtasmaniansparsegrid.so"
                                        or
                                        the lirary install prefix given
                                        to cmake

        if tasmanian_library is a string: use it as path, for example:
                TasmanianSparseGrid(
                           tasmanian_library =
                           "opt/Tasmanian/lib/libtasmaniansparsegrid.so"
                          )

        othwerwise: tasmanian_library must be an instance of ctypes.cdll
                    this is useful when creating lots of instances of
                    TasmanianSparseGrid in order to avoid having each
                    instance load a sepatare copy of the common library

        '''
        if (isinstance(tasmanian_library, int)):
            self.pLibTSG = cdll.LoadLibrary("libtasmaniansparsegrid.dll")
        elif (isinstance(tasmanian_library, basestring)):
            self.pLibTSG = cdll.LoadLibrary(tasmanian_library)
        else:
            self.pLibTSG = tasmanian_library

        # ctypes requires that we manually specify the return types of functions
        # in C this is done by header files, so this serves as a header
        self.pLibTSG.tsgConstructTasmanianSparseGrid.restype = c_void_p
        self.pLibTSG.tsgGetVersion.restype = c_char_p
        self.pLibTSG.tsgGetLicense.restype = c_char_p
        self.pLibTSG.tsgGetVersionMajor.restype = c_int
        self.pLibTSG.tsgGetVersionMinor.restype = c_int
        self.pLibTSG.tsgIsCudaEnabled.restype = c_int
        self.pLibTSG.tsgIsBLASEnabled.restype = c_int
        self.pLibTSG.tsgIsOpenMPEnabled.restype = c_int
        self.pLibTSG.tsgGetNumDimensions.restype = c_int
        self.pLibTSG.tsgGetNumOutputs.restype = c_int
        self.pLibTSG.tsgGetNumLoaded.restype = c_int
        self.pLibTSG.tsgGetNumNeeded.restype = c_int
        self.pLibTSG.tsgGetNumPoints.restype = c_int
        self.pLibTSG.tsgRead.restype = c_int
        self.pLibTSG.tsgGetAlpha.restype = c_double
        self.pLibTSG.tsgGetBeta.restype = c_double
        self.pLibTSG.tsgGetOrder.restype = c_int
        self.pLibTSG.tsgGetRule.restype = c_char_p
        self.pLibTSG.tsgGetCustomRuleDescription.restype = c_char_p
        self.pLibTSG.tsgGetLoadedPoints.restype = POINTER(c_double)
        self.pLibTSG.tsgGetNeededPoints.restype = POINTER(c_double)
        self.pLibTSG.tsgGetPoints.restype = POINTER(c_double)
        self.pLibTSG.tsgGetQuadratureWeights.restype = POINTER(c_double)
        self.pLibTSG.tsgGetInterpolationWeights.restype = POINTER(c_double)
        self.pLibTSG.tsgBatchGetInterpolationWeights.restype = POINTER(c_double)
        self.pLibTSG.tsgIsSetDomainTransfrom.restype = c_int
        self.pLibTSG.tsgIsSetConformalTransformASIN.restype = c_int
        self.pLibTSG.tsgEstimateAnisotropicCoefficients.restype = POINTER(c_int)
        self.pLibTSG.tsgEvalHierarchicalFunctions.restype = POINTER(c_double)
        self.pLibTSG.tsgBatchEvalHierarchicalFunctions.restype = POINTER(c_double)
        self.pLibTSG.tsgGetSurpluses.restype = POINTER(c_double)
        self.pLibTSG.tsgGetGlobalPolynomialSpace.restype = POINTER(c_int)
        self.pLibTSG.tsgGetAccelerationType.restype = POINTER(c_int)
        self.pLibTSG.tsgGetGPUID.restype = POINTER(c_int)
        self.pLibTSG.tsgGetNumGPUs.restype = c_int
        self.pLibTSG.tsgGetGPUmemory.restype = c_int
        self.pLibTSG.tsgGetGPUname.restype = c_char_p

        self.pLibTSG.tsgDestructTasmanianSparseGrid.argtypes = [c_void_p]
        self.pLibTSG.tsgCopyGrid.argtypes = [c_void_p, c_void_p]
        self.pLibTSG.tsgErrorLogCerr.argtypes = [c_void_p]
        self.pLibTSG.tsgDisableErrorLog.argtypes = [c_void_p]
        self.pLibTSG.tsgWrite.argtypes = [c_void_p, c_char_p]
        self.pLibTSG.tsgRead.argtypes = [c_void_p, c_char_p]
        self.pLibTSG.tsgMakeGlobalGrid.argtypes = [c_void_p, c_int, c_int, c_int, c_char_p, c_char_p, POINTER(c_int), c_double, c_double, c_char_p]
        self.pLibTSG.tsgMakeSequenceGrid.argtypes = [c_void_p, c_int, c_int, c_int, c_char_p, c_char_p, POINTER(c_int)]
        self.pLibTSG.tsgMakeLocalPolynomialGrid.argtypes = [c_void_p, c_int, c_int, c_int, c_int, c_char_p]
        self.pLibTSG.tsgMakeWaveletGrid.argtypes = [c_void_p, c_int, c_int, c_int, c_int]
        self.pLibTSG.tsgUpdateGlobalGrid.argtypes = [c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.pLibTSG.tsgUpdateSequenceGrid.argtypes = [c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.pLibTSG.tsgGetAlpha.argtypes = [c_void_p]
        self.pLibTSG.tsgGetBeta.argtypes = [c_void_p]
        self.pLibTSG.tsgGetOrder.argtypes = [c_void_p]
        self.pLibTSG.tsgGetNumDimensions.argtypes = [c_void_p]
        self.pLibTSG.tsgGetNumOutputs.argtypes = [c_void_p]
        self.pLibTSG.tsgGetRule.argtypes = [c_void_p]
        self.pLibTSG.tsgGetCustomRuleDescription.argtypes = [c_void_p]
        self.pLibTSG.tsgGetNumLoaded.argtypes = [c_void_p]
        self.pLibTSG.tsgGetNumNeeded.argtypes = [c_void_p]
        self.pLibTSG.tsgGetNumPoints.argtypes = [c_void_p]
        self.pLibTSG.tsgGetLoadedPoints.argtypes = [c_void_p]
        self.pLibTSG.tsgGetNeededPoints.argtypes = [c_void_p]
        self.pLibTSG.tsgGetPoints.argtypes = [c_void_p]
        self.pLibTSG.tsgGetQuadratureWeights.argtypes = [c_void_p]
        self.pLibTSG.tsgGetInterpolationWeights.argtypes = [c_void_p, POINTER(c_double)]
        self.pLibTSG.tsgLoadNeededPoints.argtypes = [c_void_p, POINTER(c_double)]
        self.pLibTSG.tsgEvaluate.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
        self.pLibTSG.tsgEvaluateFast.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
        self.pLibTSG.tsgIntegrate.argtypes = [c_void_p, POINTER(c_double)]
        self.pLibTSG.tsgEvaluateBatch.argtypes = [c_void_p, POINTER(c_double), c_int, POINTER(c_double)]
        self.pLibTSG.tsgBatchGetInterpolationWeights.argtypes = [c_void_p, POINTER(c_double), c_int]
        self.pLibTSG.tsgIsGlobal.argtypes = [c_void_p]
        self.pLibTSG.tsgIsSequence.argtypes = [c_void_p]
        self.pLibTSG.tsgIsLocalPolynomial.argtypes = [c_void_p]
        self.pLibTSG.tsgIsWavelet.argtypes = [c_void_p]
        self.pLibTSG.tsgSetDomainTransform.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
        self.pLibTSG.tsgIsSetDomainTransfrom.argtypes = [c_void_p]
        self.pLibTSG.tsgClearDomainTransform.argtypes = [c_void_p]
        self.pLibTSG.tsgGetDomainTransform.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
        self.pLibTSG.tsgSetConformalTransformASIN.argtypes = [c_void_p, POINTER(c_int)]
        self.pLibTSG.tsgIsSetConformalTransformASIN.argtypes = [c_void_p]
        self.pLibTSG.tsgClearConformalTransform.argtypes = [c_void_p]
        self.pLibTSG.tsgGetConformalTransformASIN.argtypes = [c_void_p, POINTER(c_int)]
        self.pLibTSG.tsgSetAnisotropicRefinement.argtypes = [c_void_p, c_char_p, c_int, c_int]
        self.pLibTSG.tsgEstimateAnisotropicCoefficients.argtypes = [c_void_p, c_char_p, c_int, POINTER(c_int)]
        self.pLibTSG.tsgSetGlobalSurplusRefinement.argtypes = [c_void_p, c_double, c_int]
        self.pLibTSG.tsgSetLocalSurplusRefinement.argtypes = [c_void_p, c_double, c_char_p, c_int]
        self.pLibTSG.tsgClearRefinement.argtypes = [c_void_p]
        self.pLibTSG.tsgRemovePointsBySurplus.argtypes = [c_void_p, c_double, c_int]
        self.pLibTSG.tsgEvalHierarchicalFunctions.argtypes = [c_void_p, POINTER(c_double)]
        self.pLibTSG.tsgBatchEvalHierarchicalFunctions.argtypes = [c_void_p, POINTER(c_double), c_int]
        self.pLibTSG.tsgSetHierarchicalCoefficients.argtypes = [c_void_p, POINTER(c_double)]
        self.pLibTSG.tsgGetSurpluses.argtypes = [c_void_p]
        self.pLibTSG.tsgGetGlobalPolynomialSpace.argtypes = [c_void_p, c_int, POINTER(c_int)]
        self.pLibTSG.tsgPrintStats.argtypes = [c_void_p]
        self.pLibTSG.tsgEnableAcceleration.argtypes = [c_void_p, c_char_p]
        self.pLibTSG.tsgGetAccelerationType.argtypes = [c_void_p]
        self.pLibTSG.tsgSetGPUID.argtypes = [c_void_p, c_int]
        self.pLibTSG.tsgGetGPUID.argtypes = [c_void_p]
        self.pLibTSG.tsgGetGPUmemory.argtypes = [c_int]
        self.pLibTSG.tsgGetGPUname.argtypes = [c_int]

        self.pLibTSG.tsgDeleteDoubles.argtypes = [POINTER(c_double)]
        self.pLibTSG.tsgDeleteInts.argtypes = [POINTER(c_int)]
        self.pLibTSG.tsgDeleteChars.argtypes = [c_char_p]

        self.pGrid = self.pLibTSG.tsgConstructTasmanianSparseGrid()

    def __del__(self):
        '''
        destructor, calls the C++ destructor and releases all memory
        used by this instance of the class
        Make sure to call this to avoid memory leaks

        '''
        self.pLibTSG.tsgDestructTasmanianSparseGrid(self.pGrid)

    def getVersion(self):
        '''
        returns the hardcoded version string from the library

        '''
        sVersion = self.pLibTSG.tsgGetVersion()
        if (sys.version_info.major == 3):
            sVersion = str(sVersion, encoding='utf8')
        return sVersion

    def getLicense(self):
        '''
        returns the hardcoded license string from the library

        '''
        sLicense = self.pLibTSG.tsgGetLicense()
        if (sys.version_info.major == 3):
            sLicense = str(sLicense, encoding='utf8')
        return sLicense
        
    def getVersionMajor(self):
        '''
        returns the hardcoded version major int
        '''
        return self.pLibTSG.tsgGetVersionMajor()
        
    def getVersionMinor(self):
        '''
        returns the hardcoded version minor int
        '''
        return self.pLibTSG.tsgGetVersionMinor()
        
    def isCudaEnabled(self):
        '''
        returns True if the library has been compiled with CUDA support
        '''
        return (self.pLibTSG.tsgIsCudaEnabled() == 0)
        
    def isBLASEnabled(self):
        '''
        returns True if the library has been compiled with BLAS support
        '''
        return (self.pLibTSG.tsgIsBLASEnabled() == 0)
        
    def isOpenMPEnabled(self):
        '''
        returns True if the library has been compiled with OpenMP support
        '''
        return (self.pLibTSG.tsgIsOpenMPEnabled() == 0)
    
    def setErrorLogCerr(self):
        '''
        sets the error log from the C++ library to the cerr stream
        '''
        self.pLibTSG.tsgErrorLogCerr(self.pGrid)
        
    def disableLog(self):
        '''
        disables the errors from the C++ library
        '''
        self.pLibTSG.tsgDisableErrorLog(self.pGrid)

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
        if (sys.version_info.major == 3):
            sFilename = bytes(sFilename, encoding='utf8')
        return (self.pLibTSG.tsgRead(self.pGrid, c_char_p(sFilename)) == 0)

    def write(self, sFilename):
        '''
        writes the grid to a file

        sFilename: string indicating a grid file where a grid will
                   be written

        '''
        if (sys.version_info.major == 3):
            sFilename = bytes(sFilename, encoding='utf8')
        self.pLibTSG.tsgWrite(self.pGrid, c_char_p(sFilename))

    def makeGlobalGrid(self, iDimension, iOutputs, iDepth, sType, sRule, liAnisotropicWeights=[], fAlpha=0.0, fBeta=0.0, sCustomFilename=""):
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
                measure (1-x)^alpha*(1+x)^beta

              'gauss-laguerre'
                approximation using roots of polynomials orthogonal in
                measure x^alpha*epx(-x)

              'gauss-hermite'  'gauss-hermite-odd'
                approximation using roots of polynomials orthogonal in
                measure |x|^alpha*epx(-x^2)

        liAnisotropicWeights: list or numpy.array of weights
                              length must be iDimension or 2*iDimension
                              the first iDimension wegiths
                                                       must be positive
                              see the manual for details

        fAlpha, fBeta: floats
              fAlpha : the alpha parameter for Gegenbauer, Jacobi,
                       Hermite and Laguerre rules
              fBeta  : the beta parameter for Jacobi rules

        sCustomRule: string giving the parth to the file with
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
                pAnisoWeights = (c_int*iNumWeights)()
                for iI in range(iNumWeights):
                    pAnisoWeights[iI] = liAnisotropicWeights[iI]
        pCustomRule = None
        if (sCustomFilename):
            pCustomRule = c_char_p(sCustomFilename)

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')
            sRule = bytes(sRule, encoding='utf8')

        self.pLibTSG.tsgMakeGlobalGrid(self.pGrid, iDimension, iOutputs, iDepth, c_char_p(sType), c_char_p(sRule), pAnisoWeights, c_double(fAlpha), c_double(fBeta), pCustomRule)

    def makeSequenceGrid(self, iDimension, iOutputs, iDepth, sType, sRule, liAnisotropicWeights=[]):
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

        liAnisotropicWeights: list or numpy.array of weights
                              length must be iDimension or 2*iDimension
                              the first iDimension wegiths
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

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')
            sRule = bytes(sRule, encoding='utf8')

        self.pLibTSG.tsgMakeSequenceGrid(self.pGrid, iDimension, iOutputs, iDepth, c_char_p(sType), c_char_p(sRule), pAnisoWeights)

    def makeLocalPolynomialGrid(self, iDimension, iOutputs, iDepth, iOrder=1, sRule="localp"):
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
                   tripples the number of points per level (as opposed
                   to double for the other cases)

        sRule: string (defines the 1-D rule that induces the grid)
              'localp'   'localp-zero'   'semi-localp'

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

        if (sys.version_info.major == 3):
            sRule = bytes(sRule, encoding='utf8')

        self.pLibTSG.tsgMakeLocalPolynomialGrid(self.pGrid, iDimension, iOutputs, iDepth, iOrder, c_char_p(sRule))

    def makeWaveletGrid(self, iDimension, iOutputs, iDepth, iOrder=1):
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

        self.pLibTSG.tsgMakeWaveletGrid(self.pGrid, iDimension, iOutputs, iDepth, iOrder)

    def copyGrid(self, pGrid):
            '''
            accepts an instance of TasmanianSparseGrid class and creates
            a hard copy of the class and all included data
            original class is not modified

            pGrid: instance of TasmanianSparseGrid class
                   the source for the copy

            '''
            if (not isinstance(pGrid, TasmanianSparseGrid)):
                    raise TasmanianInputError("pGrid", "ERROR: pGrid must be an instance of TasmanianSparseGrid")
            self.pLibTSG.tsgCopyGrid(self.pGrid, pGrid.pGrid)

    def updateGlobalGrid(self, iDepth, sType, liAnisotropicWeights=[]):
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

        self.pLibTSG.tsgUpdateGlobalGrid(self.pGrid, iDepth, sType, pAnisoWeights)

    def updateSequenceGrid(self, iDepth, sType, liAnisotropicWeights=[]):
        '''
        adds the points defined by depth, type and anisotropy
        to the existing grid

        basically, the same as calling makeGlobalGrid with sRule,
                   of this grid and the new iDepth, sType and
                   liAnisotropicWeights then adding the resulting points
                   to the current grid

        inputs: see help(TasmanianSG.TasmanianSparseGrid.makeGlobalGrid)

        '''
        if (not self.isSequence()):
            raise TasmanianInputError("updateSequenceGrid", "ERROR: calling updateSequenceGrid for a grid that is not a sequence grid")
            
        if (iDepth < 0):
            raise TasmanianInputError("iDepth", "ERROR: depth should be a non-negative integer")

        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")

        iDimension = self.getNumDimensions()
        pAnisoWeights = c_void_p(0)
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

        self.pLibTSG.tsgUpdateSequenceGrid(self.pGrid, iDepth, sType, pAnisoWeights)

    def getAlpha(self):
        '''
        returns the value of fAlpha in the call to makeGlobalGrid
        if makeGlobalGrid has not been called, returns 0.0

        '''
        return self.pLibTSG.tsgGetAlpha(self.pGrid)

    def getBeta(self):
        '''
        returns the value of fBeta in the call to makeGlobalGrid
        if makeGlobalGrid has not been called, returns 0.0

        '''
        return self.pLibTSG.tsgGetBeta(self.pGrid)

    def getOrder(self):
        '''
        returns the value of iOrder in the call to
        makeLocalPolynomialGrid or makeWaveletGrid

        if makeLocalPolynomialGrid and makeWaveletGrid
        have not been called, returns -1

        '''
        return self.pLibTSG.tsgGetOrder(self.pGrid)

    def getNumDimensions(self):
        '''
        returns the value of iDimension in the make***Grid command
        if no grid has been made, it returns 0

        '''
        return self.pLibTSG.tsgGetNumDimensions(self.pGrid)

    def getNumOutputs(self):
        '''
        returns the value of iOutputs in the make***Grid command
        if no grid has been made, it returns 0

        '''
        return self.pLibTSG.tsgGetNumOutputs(self.pGrid)

    def getRule(self):
        '''
        returns the value of sRule in the make***Grid command
        if makeWaveletGrid is used, returns "wavelet"
        if no grid has been made, it returns "unknown"

        '''
        sRule = self.pLibTSG.tsgGetRule(self.pGrid)
        if (sys.version_info.major == 3):
            sRule = str(sRule, encoding='utf8')
        return sRule

    def getCustomRuleDescription(self):
        '''
        returns the description provided in the custom rule file
        if not using a custom grid, returns ""

        '''
        if ("custom-tabulated" in self.getRule()):
            sRule = self.pLibTSG.tsgGetCustomRuleDescription(self.pGrid)
            if (sys.version_info.major == 3):
                sRule = str(sRule, encoding='utf8')
            return sRule
        else:
            return ""

    def getNumLoaded(self):
        '''
        returns the number of points loaded in the existing interpolant

        '''
        return self.pLibTSG.tsgGetNumLoaded(self.pGrid)

    def getNumNeeded(self):
        '''
        returns the number of points needed to form the interpolant or
        form the next interpolant following a refinement

        '''
        return self.pLibTSG.tsgGetNumNeeded(self.pGrid)

    def getNumPoints(self):
        '''
        if points have been loaded, returns the same as getNumLoaded()
        otherwise, returns the same as getNumNeeded()

        '''
        return self.pLibTSG.tsgGetNumPoints(self.pGrid)

    def getLoadedPoints(self):
        '''
        returns the points loaded in the existing interpolant

        output: a 2-D numpy.array of size getNumLoaded() X iDimension
            reach row correspoinds to one point
            if (getNumLoaded() == 0): returns numpy.empty([0,0])

        '''
        pPoints = self.pLibTSG.tsgGetLoadedPoints(self.pGrid)
        iDimensions = self.getNumDimensions()
        iNumPoints = self.getNumLoaded()
        if (iNumPoints == 0):
            return np.empty([0, 0], np.float64)
        llfPoints = np.empty([iNumPoints, iDimensions], np.float64)
        for iI in range(iNumPoints):
            for iJ in range(iDimensions):
                llfPoints[iI][iJ] = pPoints[iI*iDimensions + iJ]
        self.pLibTSG.tsgDeleteDoubles(pPoints)
        return llfPoints

    def getNeededPoints(self):
        '''
        returns the points needed to form the interpolant or the next
        level of refinement following a set***Refinement() call

        output: 2-D numpy.array of size getNumNeeded() X iDimension
            reach row correspoinds to one point
            if (getNumNeeded() == 0): returns numpy.empty([0,0])

        '''
        pPoints = self.pLibTSG.tsgGetNeededPoints(self.pGrid)
        iDimensions = self.getNumDimensions()
        iNumPoints = self.getNumNeeded()
        llfPoints = np.empty([iNumPoints, iDimensions], np.float64)
        for iI in range(iNumPoints):
            for iJ in range(iDimensions):
                llfPoints[iI][iJ] = pPoints[iI*iDimensions + iJ]
        self.pLibTSG.tsgDeleteDoubles(pPoints)
        return llfPoints

    def getPoints(self):
        '''
        if points have been loaded, gives the same as getLoadedPoints()
        otherwise, returns the same as getNeededPoints()

        '''
        pPoints = self.pLibTSG.tsgGetPoints(self.pGrid)
        iDimensions = self.getNumDimensions()
        iNumPoints = self.getNumPoints()
        llfPoints = np.empty([iNumPoints, iDimensions], np.float64)
        for iI in range(iNumPoints):
            for iJ in range(iDimensions):
                llfPoints[iI][iJ] = pPoints[iI*iDimensions + iJ]
        self.pLibTSG.tsgDeleteDoubles(pPoints)
        return llfPoints

    def getQuadratureWeights(self):
        '''
        returns the quadrature weights associated with
        the points in getPoints()

        output: a 1-D numpy.array of length getNumPoints()
                the order of the weights matches
                the order in getPoints()

        '''
        pWeights = self.pLibTSG.tsgGetQuadratureWeights(self.pGrid)
        iNumPoints = self.getNumPoints()
        lfWeights = np.empty([iNumPoints,], np.float64)
        for iI in range(iNumPoints):
            lfWeights[iI] = pWeights[iI]
        self.pLibTSG.tsgDeleteDoubles(pWeights)
        return lfWeights

    def getInterpolationWeights(self, lfX):
        '''
        returns the interpolation weights associated with the points
        in getPoints()

        lfX: a 1-D numpy.array with length iDimensions
             the entries indicate the points for evaluating the weights

        output: a 1-D numpy.array of length getNumPoints()
            the order of the weights matches the order in getPoints()

        '''
        iNumX = len(lfX)
        if (iNumX != self.getNumDimensions()):
                raise TasmanianInputError("lfX", "ERROR: len(lfX) should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumX))
        pX = (c_double*iNumX)()
        for iI in range(iNumX):
            pX[iI] = lfX[iI]
        pWeights = self.pLibTSG.tsgGetInterpolationWeights(self.pGrid, pX)
        iNumPoints = self.getNumPoints()
        lfWeights = np.empty([iNumPoints], np.float64)
        for iI in range(iNumPoints):
            lfWeights[iI] = pWeights[iI]
        self.pLibTSG.tsgDeleteDoubles(pWeights)
        return lfWeights

    def getInterpolationWeightsBatch(self, llfX):
        '''
        returns the interpolation weights associated with the points
        in getPoints()

        finds multiple weights with a single library call
        uses OpenMP if enabled in libtasmaniansparsegrids.so

        llfX: a 2-D numpy.array with second dimension iDimensions
              each row in the array is a single requested point

        output: a 2-D numpy.array
                with dimensions llfX.shape[0] X getNumPoints()
                each row corresponds to the weight for one row of llfX

        '''
        if (len(llfX.shape) != 2):
            raise TasmanianInputError("llfX", "ERROR: llfX should be a 2-D numpy.array instread it has dimension {0:1d}".format(len(llfX.shape)))
        iNumX = llfX.shape[0]
        if (iNumX == 0):
            return numpy.empty([0, self.getNumPoints()], np.float64)
        iNumDim = llfX.shape[1]
        if (iNumDim != self.getNumDimensions()):
            raise TasmanianInputError("llfX", "ERROR: llfX.shape[1] should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumDim))
        pX = (c_double*(iNumX*iNumDim))()
        for iI in range(iNumX):
            for iJ in range(iNumDim):
                pX[iI*iNumDim + iJ] = llfX[iI][iJ]
        iNumPoints = self.getNumPoints()
        pWeights = self.pLibTSG.tsgBatchGetInterpolationWeights(self.pGrid, pX, iNumX)
        llfWeights = np.empty([iNumX, iNumPoints], np.float64)
        for iI in range(iNumX):
            for iJ in range(iNumPoints):
                llfWeights[iI][iJ] = pWeights[iI*iNumPoints + iJ]
        self.pLibTSG.tsgDeleteDoubles(pWeights)
        return llfWeights

    def loadNeededPoints(self, llfVals):
        '''
        loads the values of the target function at the needed points
        if there are no needed points, this reset the currently loaded
        values

        llfVals: a 2-D numpy.array
                 with dimensions getNumNeeded() X iOutputs
                 each row corresponds to the values of the outputs at
                 the corresponding needed point. The order and leading
                 dimension must match the points obtained form
                 getNeededPoints()

        '''
        if (len(llfVals.shape) != 2):
            raise TasmanianInputError("llfVals", "ERROR: llfVals should be a 2-D numpy.array, instead it has {0:1d} dimensions".format(len(llfVals.shape)))
        if (llfVals.shape[0] != self.getNumNeeded()):
            if (self.getNumNeeded() == 0):
                if (llfVals.shape[0] != self.getNumLoaded()):
                    raise TasmanianInputError("llfVals", "ERROR: leading dimension of llfVals is {0:1d} but the number of current points is {1:1d}".format(llfVals.shape[0], self.getNumNeeded()))
            else:
                raise TasmanianInputError("llfVals", "ERROR: leading dimension of llfVals is {0:1d} but the number of needed points is {1:1d}".format(llfVals.shape[0], self.getNumNeeded()))
        if (llfVals.shape[1] != self.getNumOutputs()):
            raise TasmanianInputError("llfVals", "ERROR: second dimension of llfVals is {0:1d} but the number of outputs is set to {1:1d}".format(llfVals.shape[1], self.getNumOutputs()))
        iNumPoints = llfVals.shape[0]
        iNumDims = llfVals.shape[1]
        pVals = (c_double*(iNumPoints*iNumDims))()
        for iI in range(iNumPoints):
            for iJ in range(iNumDims):
                pVals[iI*iNumDims + iJ] = llfVals[iI][iJ]
        self.pLibTSG.tsgLoadNeededPoints(self.pGrid, pVals)

    def evaluateThreadSafe(self, lfX):
    #def evaluate(self, lfX):
        '''
        evaluates the intepolant at a single points of interest and
        returns the result
        This is the thread safe version, but it does not use 
        acceleration of any type

        this should be called after the grid has been created and after
        values have been loaded

        lfX: a 1-D numpy.array with length iDimensions
             the entries indicate the points for evaluating the weights

        output: returns a 1-D numpy.array of length iOutputs
            the values of the interpolant at lfX

        '''
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("evaluate", "ERROR: cannot call evaluate for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        iNumX = len(lfX)
        if (iNumX != self.getNumDimensions()):
            raise TasmanianInputError("lfX", "ERROR: lfX should have lenth {0:1d} instead it has length {1:1d}".format(self.getNumDimensions(),iNumX))
        pX = (c_double*iNumX)()
        for iI in range(iNumX):
            pX[iI] = lfX[iI]
        iNumOutputs = self.getNumOutputs()
        pY = (c_double*iNumOutputs)()
        self.pLibTSG.tsgEvaluate(self.pGrid, pX, pY)
        lfY = np.empty([iNumOutputs,], np.float64)
        for iI in range(iNumOutputs):
            lfY[iI] = pY[iI]
        return lfY
    
    def evaluate(self, lfX):
        '''
        evaluates the intepolant at a single points of interest and
        returns the result
        This is the accelerated version using the selected acceleration
        type, but it is potentially not thread safe

        this should be called after the grid has been created and after
        values have been loaded

        lfX: a 1-D numpy.array with length iDimensions
             the entries indicate the points for evaluating the weights

        output: returns a 1-D numpy.array of length iOutputs
            the values of the interpolant at lfX

        '''
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("evaluate", "ERROR: cannot call evaluate for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        iNumX = len(lfX)
        if (iNumX != self.getNumDimensions()):
            raise TasmanianInputError("lfX", "ERROR: lfX should have lenth {0:1d} instead it has length {1:1d}".format(self.getNumDimensions(),iNumX))
        pX = (c_double*iNumX)()
        for iI in range(iNumX):
            pX[iI] = lfX[iI]
        iNumOutputs = self.getNumOutputs()
        pY = (c_double*iNumOutputs)()
        self.pLibTSG.tsgEvaluateFast(self.pGrid, pX, pY)
        lfY = np.empty([iNumOutputs,], np.float64)
        for iI in range(iNumOutputs):
            lfY[iI] = pY[iI]
        return lfY

    def evaluateBatch(self, llfX):
        '''
        evaluates the intepolant at the points of interest and returns
        the result

        this should be called after the grid has been created and after
        values have been loaded

        llfX: a 2-D numpy.array
              with second dimension equal to iDimensions
              each row in the array is a single requested point

        output: a 2-D numpy.array
                with dimensions llfX.shape[0] X iOutputs
                each row corresponds to the value of the interpolant
                for one row of llfX

        '''
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("evaluateBatch", "ERROR: cannot call evaluateBatch for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (len(llfX.shape) != 2):
            raise TasmanianInputError("llfX", "ERROR: llfX should be a 2-D numpy.array instread it has dimension {0:1d}".format(len(llfX.shape)))
        iNumX = llfX.shape[0]
        if (iNumX == 0):
            return numpy.empty([0, self.getNumPoints()], np.float64)
        iNumDim = llfX.shape[1]
        if (iNumDim != self.getNumDimensions()):
            raise TasmanianInputError("llfX", "ERROR: llfX.shape[1] should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumDim))
        iNumDim = llfX.shape[1]
        pX = (c_double*(iNumX*iNumDim))()
        for iI in range(iNumX):
            for iJ in range(iNumDim):
                pX[iI*iNumDim + iJ] = llfX[iI][iJ]
        iNumOutputs = self.getNumOutputs()
        pY = (c_double*(iNumOutputs*iNumX))()
        self.pLibTSG.tsgEvaluateBatch(self.pGrid, pX, iNumX, pY)
        llfY = np.empty([iNumX, iNumOutputs], np.float64)
        for iI in range(iNumX):
            for iJ in range(iNumOutputs):
                llfY[iI][iJ] = pY[iI*iNumOutputs + iJ]
        return llfY

    def integrate(self):
        '''
        returns the integral of the interpolant

        output: returns a 1-D numpy.array of length iOutputs
            the integral of the interpolant

        '''
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("integrate", "ERROR: cannot call integrate for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        iNumOutputs = self.getNumOutputs()
        pQ = (c_double*iNumOutputs)()
        self.pLibTSG.tsgIntegrate(self.pGrid, pQ)
        lfQ = np.empty([iNumOutputs,], np.float64)
        for iI in range(iNumOutputs):
            lfQ[iI] = pQ[iI]
        return lfQ

    def isGlobal(self):
        '''
        returns True if using a global grid

        '''
        return (self.pLibTSG.tsgIsGlobal(self.pGrid) == 0)

    def isSequence(self):
        '''
        returns True if using a sequence grid

        '''
        return (self.pLibTSG.tsgIsSequence(self.pGrid) == 0)

    def isLocalPolynomial(self):
        '''
        returns True if using a local polynomial grid

        '''
        return (self.pLibTSG.tsgIsLocalPolynomial(self.pGrid) == 0)

    def isWavelet(self):
        '''
        returns True if using a local wavelet grid

        '''
        return (self.pLibTSG.tsgIsWavelet(self.pGrid) == 0)

    def setDomainTransform(self, llfTransform):
        '''
        sets the lower and upper bound for each dimension

        Note: gauss-laguerre and gauss-hermite rules are defined on
              unbounded domain, in which case  this sets the
              shift and scale parameters

        llfTransform: a 2-D numpy.array of size iDimension X 2
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
            raise TasmanianInputError("llfTransform", "ERROR: llfTransform should be a 2-D numpy.array")
        if (lShape[0] != self.getNumDimensions()):
            raise TasmanianInputError("llfTransform", "ERROR: the first dimension of llfTransform is {0:1d} and it should match iDimension: {1:1d}".format(lShape[0], self.getNumDimensions()))
        if (lShape[1] != 2):
            raise TasmanianInputError("llfTransform", "ERROR: the second dimension of llfTransform is {0:1d} and it should be 2".format(lShape[1]))
        iNumDimensions = llfTransform.shape[0]
        pA = (c_double*iNumDimensions)()
        pB = (c_double*iNumDimensions)()
        for iI in range(iNumDimensions):
            pA[iI] = llfTransform[iI][0]
            pB[iI] = llfTransform[iI][1]
        self.pLibTSG.tsgSetDomainTransform(self.pGrid, pA, pB)

    def isSetDomainTransfrom(self):
        '''
        returns True if the grid is defined for non-canonical domain
        returns False if using a canonical domain

        '''
        return (self.pLibTSG.tsgIsSetDomainTransfrom(self.pGrid) == 0)

    def clearDomainTransform(self):
        '''
        resets the domain to canonical
        loaded values will be kept, however, the values now correspond
        to canonical points and may be invalid for your applicaiton

        '''
        self.pLibTSG.tsgClearDomainTransform(self.pGrid)

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
        self.pLibTSG.tsgGetDomainTransform(self.pGrid, pA, pB)
        llfTransform = np.empty([iNumDimensions, 2], np.float64)
        for iI in range(iNumDimensions):
            llfTransform[iI][0] = pA[iI]
            llfTransform[iI][1] = pB[iI]
        return llfTransform
        
    def setConformalTransformASIN(self, liTruncation):
        '''
        sets conformal domain transform based on truncated 
        Maclaurin series of arcsin()
        
        liTruncation: 1-D numpy.array of non-negative integers
                      indicating the truncation order in each direction
                      0 indicates no transform applied to this direction
        
        '''
        lShape = liTruncation.shape
        if (len(lShape) != 1):
            raise TasmanianInputError("liTruncation", "ERROR: liTruncation should be a 1-D numpy.array")
        if (lShape[0] != self.getNumDimensions()):
            raise TasmanianInputError("liTruncation", "ERROR: the length of liTruncation is {0:1d} and it should match iDimension: {1:1d}".format(lShape[0], self.getNumDimensions()))
        iNumDimensions = lShape[0]
        pTruncation = (c_int*iNumDimensions)()
        for iI in range(iNumDimensions):
            pTruncation[iI] = liTruncation[iI]
        self.pLibTSG.tsgSetConformalTransformASIN(self.pGrid, pTruncation)
            
    def isSetConformalTransformASIN(self):
        '''
        returns True if conformal transform is set
        returns False otherwise
        
        see: setConformalTransformASIN()

        '''
        return (self.pLibTSG.tsgIsSetConformalTransformASIN(self.pGrid) == 0)
            
    def clearConformalTransform(self):
        '''
        resets the conformal domain transform
        loaded values will be kept, however, the values now correspond
        to canonical points and may be invalid for your applicaiton

        '''
        self.pLibTSG.tsgClearConformalTransform(self.pGrid)
            
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
        self.pLibTSG.tsgGetConformalTransformASIN(self.pGrid, pTruncation)
        liTruncation = np.empty([iNumDimensions,], np.int)
        for iI in range(iNumDimensions):
            liTruncation[iI] = pTruncation[iI]
        return liTruncation

    def setAnisotropicRefinement(self, sType, iMinGrowth, iOutput):
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
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("setAnisotropicRefinement", "ERROR: cannot call setAnisotropicRefinement for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (self.getNumOutputs() == 0):
            raise TasmanianInputError("setAnisotropicRefinement", "ERROR: cannot set refinement for grid with iOutput = 0")
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("setAnisotropicRefinement", "ERROR: cannot set refinement before values are loaded")
        if (not iMinGrowth > 0):
            raise TasmanianInputError("iMinGrowth", "ERROR: the number of growth should be positive integer")
        if (iOutput == -1):
            if (not self.isSequence()):
                raise TasmanianInputError("iOutput", "ERROR: iOutput = -1 can be used only for sequence grids")
        if (iOutput < -1):
            raise TasmanianInputError("iOutput", "ERROR: iOutput should be -1 or a non-negative integer")
        if (iOutput >= self.getNumOutputs()):
            raise TasmanianInputError("iOutput", "ERROR: iOutput cannot exceed the index of the last output {0:1d}".format(getNumOutputs() - 1))
        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")

        if (sys.version_info.major == 3):
            sType = bytes(sType, encoding='utf8')
        self.pLibTSG.tsgSetAnisotropicRefinement(self.pGrid, c_char_p(sType), iMinGrowth, iOutput)

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

        outputs: 1-D numpy.array
                 of length getNumDimensions() or 2*getNumDimensions()
                 the first set of getNumDimensions() entries correspond
                                                 to the xi coefficients
                 the second set of getNumDimensions() entries correspond
                                                 to the eta coefficients

        '''
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("estimateAnisotropicCoefficients", "ERROR: cannot call estimateAnisotropicCoefficients for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (getNumOutputs() == 0):
            raise TasmanianInputError("setAnisotropicRefinement", "ERROR: cannot set refinement for grid with iOutput = 0")
        if (getNumLoaded() == 0):
            raise TasmanianInputError("setAnisotropicRefinement", "ERROR: cannot set refinement before values are loaded")
        if (iOutput == -1):
            if (not self.isSequence()):
                raise TasmanianInputError("iOutput", "ERROR: iOutput = -1 can be used only for sequence grids")
        if (iOutput < -1):
            raise TasmanianInputError("iOutput", "ERROR: iOutput should be -1 or a non-negative integer")
        if (iOutput >= getNumOutputs()):
            raise TasmanianInputError("iOutput", "ERROR: iOutput cannot exceed the index of the last output {0:1d}".format(getNumOutputs() - 1))
        if (sType not in lsTsgGlobalTypes):
            raise TasmanianInputError("sType", "ERROR: invalid type, see TasmanianSG.lsTsgGlobalTypes for list of accepted types")

        pNumCoefficients = (c_int*1)()
        pCoefficients = self.pLibTSG.tsgEstimateAnisotropicCoefficients(self.pGrid, c_char_p(sType), iOutput, pNumCoefficients)
        liCoefficients = np.empty([pNumCoefficients[0], ], np.int)
        for iI in range(pNumCoefficients[0]):
            liCoefficients[iI] = pCoefficients[iI]
        self.pLibTSG.tsgDeleteInts(pCoefficients)
        return liCoefficients

    def setSurplusRefinement(self, fTolerance, iOutput, sCriteria=""):
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

        sCriteria: hierarhical and direction refinement strategy
                   'classic'  'parents'   'direction'   'fds'
                  applicable only for Local Polynomial and Wavelet grids

        '''
        if (self.isGlobal()):
            raise TasmanianInputError("setSurplusRefinement", "ERROR: setSurplusRefinement cannot be used with global grids")
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("setSurplusRefinement", "ERROR: cannot call setSurplusRefinement for a grid before any points are loaded, i.e., call loadNeededPoints first!")
        if (fTolerance < 0.0):
            raise TasmanianInputError("fTolerance", "ERROR: fTolerance must be non-negative")
        if (len(sCriteria) == 0):
            if (not self.isSequence()):
                raise TasmanianInputError("sCriteria", "ERROR: sCriteria must be specified")
            self.pLibTSG.tsgSetGlobalSurplusRefinement(self.pGrid, c_double(fTolerance), iOutput)
        else:
            if (self.isSequence()):
                raise TasmanianInputError("sCriteria", "ERROR: sCriteria cannot be used for sequence grids")
            if (sys.version_info.major == 3):
                sCriteria = bytes(sCriteria, encoding='utf8')
            self.pLibTSG.tsgSetLocalSurplusRefinement(self.pGrid, c_double(fTolerance), c_char_p(sCriteria), iOutput)

    def clearRefinement(self):
        '''
        clear the last call to set***Refinement,
        only works if called before the points are loaded, i.e.,
        before loadNeededPoints()

        if getNumNeeded() == 0, this call will have no effect

        '''
        self.pLibTSG.tsgClearRefinement(self.pGrid)

    def removePointsBySurplus(self, fTolerance, iOutput):
        '''
        removes any points in the grid with relative surplus that
        exceeds the tolerance
        applies only to local polynomial grids

        fTolerance: float (positive)
                    the relative surplus tolerance, i.e.,
                    we keep only for points associated with surplus
                    that exceeds the tolerance

        iOutput: int (indicates the output to use)
                 selects which output to consider
                 accept -1 to indicate all outputs

        '''
        if (not self.isLocalPolynomial()):
            raise TasmanianInputError("removePointsBySurplus", "ERROR: calling removePointsBySurplus for a grid that isn't local polynomial")
        if (fTolerance <= 0.0):
            raise TasmanianInputError("fTolerance", "ERROR: fTolerance must be a positive integer")
        if (iOutput < -1):
            raise TasmanianInputError("iOutput", "ERROR: iOutput should be -1 or a non-negative integer")
        if (iOutput >= self.getNumOutputs()):
            raise TasmanianInputError("iOutput", "ERROR: iOutput cannot exceed the index of the last output {0:1d}".format(self.getNumOutputs() - 1))
        if (self.getNumLoaded() == 0):
            raise TasmanianInputError("removePointsBySurplus", "ERROR: calling removePointsBySurplus when no points are loades")
        self.pLibTSG.tsgRemovePointsBySurplus(self.pGrid, fTolerance, iOutput)

    def evalHierarchicalFunctions(self, lfX):
        '''
        EXPERIMENTAL CAPABILITY

        evaluates the hierarchical functions at a single points of
        interest and returns the result

        this should be called after the grid has been created and after
        values have been loaded (even if you load dummy values)

        lfX: a 1-D numpy.array with length iDimensions
             the entries indicate the points for evaluating the weights

        output: returns a 1-D numpy.array of length iOutputs
                the values of the interpolant at lfX
        '''
        iNumX = len(lfX)
        pX = (c_double*iNumX)()
        for iI in range(iNumX):
            pX[iI] = lfX[iI]
        pVals = self.pLibTSG.tsgEvalHierarchicalFunctions(self.pGrid, pX)
        iNumPoints = self.getNumPoints()
        lfVals = np.empty([iNumPoints], np.float64)
        for iI in range(iNumPoints):
            lfVals[iI] = pVals[iI]
        self.pLibTSG.tsgDeleteDoubles(pVals)
        return lfVals

    def evalBatchHierarchicalFunctions(self, llfX):
        '''
        EXPERIMENTAL CAPABILITY

        evaluates the hierarchical functions at a multiple points of
        interest and returns the result. This function uses OpenMP if
        enabled in the tasmanian library

        this should be called after the grid has been created and after
        values have been loaded (even if you load dummy values)

        llfX: a 2-D numpy.array with second dimension being iDimensions
              the entries indicate the points for evaluating the weights

        output: returns a 2-D numpy.array of size
                llfX.shape[0] x iOutputs
                the values of the interpolant at llfX points
        '''
        if (len(llfX.shape) != 2):
            raise TasmanianInputError("llfX", "ERROR: llfX should be a 2-D numpy.array instread it has dimension {0:1d}".format(len(llfX.shape)))
        iNumX = llfX.shape[0]
        if (iNumX == 0):
            return numpy.empty([0, self.getNumPoints()], np.float64)
        iNumDim = llfX.shape[1]
        if (iNumDim != self.getNumDimensions()):
            raise TasmanianInputError("llfX", "ERROR: llfX.shape[1] should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumDim))
        pX = (c_double*(iNumX*iNumDim))()
        for iI in range(iNumX):
            for iJ in range(iNumDim):
                pX[iI*iNumDim + iJ] = llfX[iI][iJ]
        iNumPoints = self.getNumPoints()
        pVals = self.pLibTSG.tsgBatchEvalHierarchicalFunctions(self.pGrid, pX, iNumX)
        llfVals = np.empty([iNumX, iNumPoints], np.float64)
        for iI in range(iNumX):
            for iJ in range(iNumPoints):
                llfVals[iI][iJ] = pVals[iI*iNumPoints + iJ]
        self.pLibTSG.tsgDeleteDoubles(pVals)
        return llfVals

    def setHierarchicalCoefficients(self, llfCoefficients):
        '''
        EXPERIMENTAL CAPABILITY

        loads the hierarchical surplus coefficients associated with the
        target function

        this must be calles after dummy values have been loaded (see
        loadNeededPoints())

        the coefficients will overwrite the existing values

        llfCoefficients: a 2-D numpy.array
                         with dimensions getNumPoints() X iOutputs
                         each row corresponds to the values of the
                         outputs at the corresponding needed point.
                         The order and leading dimension must match the
                         points obtained form getPoints()
        '''
        if (len(llfCoefficients.shape) != 2):
            raise TasmanianInputError("llfCoefficients", "ERROR: llfCoefficients should be a 2-D numpy.array, instead it has {0:1d} dimensions".format(len(llfCoefficients.shape)))
        if (llfCoefficients.shape[0] != self.getNumNeeded()):
            if (self.getNumNeeded() == 0):
                if (llfCoefficients.shape[0] != self.getNumLoaded()):
                    raise TasmanianInputError("llfCoefficients", "ERROR: leading dimension of llfCoefficients is {0:1d} but the number of current points is {1:1d}".format(llfCoefficients.shape[0], self.getNumNeeded()))
            else:
                raise TasmanianInputError("llfCoefficients", "ERROR: leading dimension of llfCoefficients is {0:1d} but the number of needed points is {1:1d}".format(llfCoefficients.shape[0], self.getNumNeeded()))
        if (llfCoefficients.shape[1] != self.getNumOutputs()):
            raise TasmanianInputError("llfCoefficients", "ERROR: second dimension of llfCoefficients is {0:1d} but the number of outputs is set to {1:1d}".format(llfCoefficients.shape[1], self.getNumOutputs()))
        iNumPoints = llfCoefficients.shape[0]
        iNumDims = llfCoefficients.shape[1]
        pCoefficients = (c_double*(iNumPoints*iNumDims))()
        for iI in range(iNumPoints):
            for iJ in range(iNumDims):
                pCoefficients[iI*iNumDims + iJ] = llfCoefficients[iI][iJ]
        self.pLibTSG.tsgSetHierarchicalCoefficients(self.pGrid, pCoefficients)

    def getSurpuses(self):
        '''
        EXPERIMENTAL CAPABILITY
        '''
        iNumOut = self.getNumOutputs()
        iNumPoints = self.getNumPoints()
        pSurp = self.pLibTSG.tsgGetSurpluses(self.pGrid)
        llfS = np.empty([iNumPoints, iNumOut], np.float64)
        for iI in range(iNumPoints):
            for iJ in range(iNumOut):
                llfS[iI][iJ] = pSurp[iI*iNumOut + iJ]
        return llfS

    def getGlobalPolynomialSpace(self, bInterpolation):
        '''
        returns a matrix corresponding to the polynomial space that is
        integrated or interpolated exactly by the current grid

        bInterpolation: boolean
                indicates whether to give the space associated
                with integration or interpolation

        output: is a 2-D numpy.array of integers
            output.shape[0] indicates the cardinality of the space
            output.shape[1] is equal to iDimension
            each row corresponds to a multi-index associated with a
            polynomial in a hierarchical tensor basis, for example,
            monomials

            see the manual for details

        '''
        iInterp = 1
        if (bInterpolation):
            iInterp = 0
        pNumIndexes = (c_int*1)()
        pIndexes = self.pLibTSG.tsgGetGlobalPolynomialSpace(self.pGrid, iInterp, pNumIndexes)
        iNumDimensions = self.getNumDimensions()
        lliPolynomials = np.empty([pNumIndexes[0], iNumDimensions], np.int)
        for iI in range(pNumIndexes[0]):
            for iJ in range(iNumDimensions):
                lliPolynomials[iI][iJ] = pIndexes[iI*iNumDimensions + iJ]
        self.pLibTSG.tsgDeleteInts(pIndexes)
        return lliPolynomials

    def enableAcceleration(self, sAccelerationType):
        '''
        enables the use of acceleration, for example GPU
        must specify acceleration type and perhaps additional 
        parameters, if needed by the acceleration type
        currently, acceleration pertains primarily to batch evaluations
        of global interpolants
        to use acceleration, you must compile using cmake and enable
        a corresponding flag

        sAccelerationType: string

          'none' 
              fallback mode, relies on sequential implementation
              if compiled with Tasmanian_ENABLE_OPENMP this will use
              simple "omp parallel for" to take advantage of multicore

          'cpu-blas'
              uses BLAS dgemm function for acceleration of batch
              evaluations
              requires Tasmanian_ENABLE_BLAS switched ON
              if enabled, this is the default mode

          'gpu-fullmem'
              uses CUDA cuBlas library for accelerated matrix 
              multiplication, i.e., cublasDgemm
              requires Tasmanian_ENABLE_CUDA
              this mode assumes that all data structures will fit in
              GPU memory, mainly the matrix of loaded values, 
              the matrix of interpolation weights, and the result matrix
              total storage (in doubles) =
                  getNumOutputs() * getNumPoints()
                + getNumPoints() * num_x
                + getNumOutputs() * num_x
              where num_x is the size of the batch
              may have to split the batch into smaller units to work
              around memory constraints

        '''
        if (sAccelerationType not in lsTsgAccelTypes):
            raise TasmanianInputError("sAccelerationType", "ERROR: invalid acceleration type")
        if (sys.version_info.major == 3):
            sAccelerationType = bytes(sAccelerationType, encoding='utf8')
        self.pLibTSG.tsgEnableAcceleration(self.pGrid, c_char_p(sAccelerationType))
        
    def getAccelerationType(self):
        '''
        returns the type of acceleration set by enableAcceleration
        '''
        return lsTsgAccelTypes[self.pLibTSG.tsgGetAccelerationType(self.pGrid)]

    def setGPUID(self, iGPUID):
        '''
        when using cuda on a machine with multiple GPUs, this helps set
        the GPU for this grid
        NOTE: each instance of the sparse grids class holds a separate
              instance of iGPUID and different grids can be assigned to
              different GPUs (on multigpu system)
        iGPUID can be changed at any time, however, this will cause
        some of the internal cache to be invalidated and it may lead
        to extraneous data movement
        
        calling read or make***Grid will reset the selected GPU
        
        defaults to 0
        
        this doesn't do anything unless enableAcceleration is called
        using a "gpu-" acceleration type
        '''
        self.pLibTSG.tsgSetGPUID(self.pGrid, iGPUID)
        
    def getGPUID(self):
        '''
        returns the GPU ID set using setGPUID
        '''
        return self.pLibTSG.tsgGetGPUID(self.pGrid)

    def getNumGPUs(self):
        '''
        returns the number of available GPUs according to cuda
        
        this is one of several functions designed to allow basic 
        management of multi-gpu setup with only Tasmanian module
        '''
        return self.pLibTSG.tsgGetNumGPUs()

    def getGPUmemory(self, iGPUID):
        '''
        returns the total memory (in MegaBytes, 1024**2 bytes) of the
        corresponding GPU
        
        this is one of several functions designed to allow basic 
        management of multi-gpu setup with only Tasmanian module
        '''
        if ((iGPUID < 0) or (iGPUID >= self.getNumGPUs())):
            raise TasmanianInputError("iGPUID", "ERROR: invalid GPU ID number")
        return self.pLibTSG.tsgGetGPUmemory(iGPUID)

    def getGPUname(self, iGPUID):
        '''
        return the cuda name ID of the corresponding GPU
        
        this is one of several functions designed to allow basic 
        management of multi-gpu setup with only Tasmanian module
        '''
        if ((iGPUID < 0) or (iGPUID >= self.getNumGPUs())):
            raise TasmanianInputError("iGPUID", "ERROR: invalid GPU ID number")
        # the pointer returned here has to be deleted
        pName = self.pLibTSG.tsgGetGPUname(iGPUID)
        if (sys.version_info.major == 3):
            sName = str(pName, encoding='utf8')
        else:
            sName = pName # create a copy to release the pointer
        return sName

    def printStats(self):
        '''
        calls the library printStats() function, which displays basic
        information about this instance of the grid
        '''
        self.pLibTSG.tsgPrintStats(self.pGrid)

    def plotPoints2D(self, sStyle="bo", iNumFigure=1, iMarkerSize=3, bShow=True):
        '''
        plots the points in a 2D plot using matplotlib.pyplot
        applicable only for grids with iDimensions == 2

        returns an instance of matplotlib.pyplot

        '''
        if (not bTsgPlotting):
            raise TasmanianInputError("plotPoints2D", "ERROR: could not load matplotlib.pyplot")
        if (self.getNumDimensions() != 2):
            raise TasmanianInputError("plotPoints2D", "ERROR: cannot plot a grid with other than 2 dimensions")

        tsgPlot.figure(iNumFigure)
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

        tsgPlot.plot(aPoints[:,0], aPoints[:,1], sStyle, markersize=iMarkerSize)
        tsgPlot.axis([fXmin - 0.1 * np.fabs(fXmin), fXmax + 0.1 * np.fabs(fYmax), fYmin - 0.1 * np.fabs(fYmin), fYmax + 0.1 * np.fabs(fYmax)])
        if (bShow):
            tsgPlot.show()
        return tsgPlot

    def plotResponse2D(self, iOutput=0, iNumDim0=100, iNumDim1=100, sStyle="bo", iNumFigure=1, iMarkerSize=3, bShow=True):
        '''
        plots the response in a 2D plot using matplotlib.pyplot
        applicable only for grids with iDimensions == 2

        iOutput is the output to use for plotting
        iNumDim0, iNumDim1: positive integers
               the points for the plot are selected on a dense grid with
               number of points iNumDim0 and iNumDim1 in dimensions
               0 and 1

        returns an instance of matplotlib.pyplot

        '''
        if (iOutput < 0):
            raise TasmanianInputError("iOutput", "ERROR: iOutput should be a non-negative integer")
        if (iOutput >= self.getNumOutputs()):
            raise TasmanianInputError("iOutput", "ERROR: iOutput cannot exceed the index of the last output {0:1d}".format(self.getNumOutputs() - 1))
        if (not bTsgPlotting):
            raise TasmanianInputError("plotPoints2D", "ERROR: could not load matplotlib.pyplot")
        if (self.getNumDimensions() != 2):
            raise TasmanianInputError("plotPoints2D", "ERROR: cannot plot a grid with other than 2 dimensions")

        tsgPlot.figure(iNumFigure)
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
        y = np.linspace(fYmin, fYmax, iNumDim1)

        XX, YY = np.meshgrid(x, y)
        ZZ = self.evaluateBatch(np.vstack((XX.reshape((iNumDim0*iNumDim1,)), YY.reshape((iNumDim0*iNumDim1,)))).T)
        ZZ = ZZ[:,iOutput].reshape((iNumDim0,iNumDim1))

        tsgPlot.imshow(ZZ, cmap="jet", extent=[fXmin, fXmax, fYmin, fYmax])
        if (bShow):
            tsgPlot.show()
        return tsgPlot
