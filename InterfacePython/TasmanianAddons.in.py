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

from ctypes import c_int, c_double, c_void_p, POINTER, CFUNCTYPE, cdll
import numpy as np
import sys

pLibCTSG = cdll.LoadLibrary("@Tasmanian_libcaddons_path@")

type_lpnmodel = CFUNCTYPE(None, c_int, POINTER(c_double), c_int, POINTER(c_double), c_int)

pLibCTSG.tsgLoadNeededPoints.argtypes = [c_int, type_lpnmodel, c_void_p, c_int]


def tsgLnpModelWrapper(oUserModel, iSizeX, pX, iSizeY, pY, iThreadID):
    '''
    DO NOT CALL DIRECTLY
    This is callback from C++, see TasGrid::loadNeededPoints()

    Creates an interface between a user callable object and
    the Tasmanian callback from C++.
    The callback passes in two raw-arrays with given sizes
    (that correspond to grid inputs and outputs)
    and the id of the running thread.
    The raw-array is wrapped in a numpy array structure and given
    to the user, the result is written to the output array.

    oUserModel: is callable object, e.g., function or lambda

    iSizeX: number of entries in x, equal to getNumDimensions()
    pX:     raw-array corresponding to the model inputs
    iSizeY: number of entries in y, equal to getNumOutputs()
    pY:     raw-array corresponding to the model outputs

    iThreadID: the id of the running thread
    '''
    aX = np.ctypeslib.as_array(pX, (iSizeX,))
    aY = np.ctypeslib.as_array(pY, (iSizeY,))
    aResult = oUserModel(aX, iThreadID)
    aY[0:iSizeY] = aResult[0:iSizeY]


def loadNeededPoints(callableModel, grid, iNumThreads = 1):
    '''
    Wrapper to TasGrid::loadNeededPoints(), non-overwrite version.

    If the grid has needed points, the callableModel will be called
    for each grid point (i.e., model input) and the resulting values
    will be loaded in the grid.

    callableModel: is callable object, e.g., function or lambda
        The object must accept two inputs and give one output:
        aY = callableModel(aX, iThreadID)

        aX: is a one dimensional numpy.ndarray with size equal
            to the number of model inputs

        iThreadID: is the ID of the thread executing the model,
            always between 0 and iNumThreads -1
            Two simultaneous calls to callableModel() will always
            have different ids.

        Return: aY must be a one dimensional numpy.ndarray with
            size equal to the number of model outputs

        Note: if iNumThreads > 1, then callableModel() must be thread-safe.

    grid: must be an instance of Tasmanian.SparseGrid()
        model values will be loaded in the grid

    iNumThreads: integer, if greater than 1 the model will be called
        in parallel from multiple threads.
        See TasGrid::loadNeededPoints().

    '''
    iOverwrite = 0 # do not overwrite
    pLibCTSG.tsgLoadNeededPoints(iOverwrite,
                                 type_lpnmodel(lambda nx, x, ny, y, tid : tsgLnpModelWrapper(callableModel, nx, x, ny, y, tid)),
                                 grid.pGrid, iNumThreads)

def reloadLoadedPoints(callableModel, grid, iNumThreads = 1):
    '''
    Wrapper to TasGrid::loadNeededPoints(), overwrite version.

    Clears any pending refinement (i.e., needed points) and overwrites the model
    values associated with the existing loaded points.

    The inputs are identical to Tasmanian.loadNeededPoints().

    '''
    iOverwrite = 1 # do overwrite
    pLibCTSG.tsgLoadNeededPoints(iOverwrite,
                                 type_lpnmodel(lambda nx, x, ny, y, tid : tsgLnpModelWrapper(callableModel, nx, x, ny, y, tid)),
                                 grid.pGrid, iNumThreads)
