# Copyright (c) 2022, Miroslav Stoyanov & Weiwei Kong
#
# This file is part of
# Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
# conditions are met:
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
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND
# IMPLIED. THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
# THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL
# ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE. THE USER ASSUMES
# RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING
# FROM OR ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.

from ctypes import c_char_p, c_int, c_double, c_void_p, POINTER, CDLL, cast, CFUNCTYPE, Structure, RTLD_GLOBAL
from numpy.ctypeslib import as_ctypes
import numpy as np
import sys

from TasmanianConfig import __path_libdream__
from TasmanianConfig import TasmanianInputError as InputError

type_optim_obj_fn_single  = CFUNCTYPE(c_double, c_int, POINTER(c_double), POINTER(c_int))
type_optim_grad_fn_single = CFUNCTYPE(None, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_int))
type_optim_proj_fn_single = CFUNCTYPE(None, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_int))

# OptimizationStatus class (mirrors the C++ struct).
class optimization_status(Structure):
    pass
optimization_status._fields_ = [("performed_iterations", c_int), ("residual", c_double)]
def clean_status(input_status):
    return dict((field, getattr(input_status, field)) for field, _ in input_status._fields_)

pLibDTSG = CDLL(__path_libdream__, mode = RTLD_GLOBAL)

pLibDTSG.tsgGradientDescentState_Construct.restype = c_void_p
pLibDTSG.tsgGradientDescentState_Construct.argtypes = [c_int, POINTER(c_double), c_double]

pLibDTSG.tsgGradientDescentState_Destruct.argtypes = [c_void_p]

pLibDTSG.tsgGradientDescentState_GetNumDimensions.restype = c_int
pLibDTSG.tsgGradientDescentState_GetNumDimensions.argtypes = [c_void_p]
pLibDTSG.tsgGradientDescentState_GetAdaptiveStepsize.restype = c_double
pLibDTSG.tsgGradientDescentState_GetAdaptiveStepsize.argtypes = [c_void_p]
pLibDTSG.tsgGradientDescentState_GetX.argtypes = [c_void_p, POINTER(c_double)]

pLibDTSG.tsgGradientDescentState_SetAdaptiveStepsize.argtypes = [c_void_p, c_double]
pLibDTSG.tsgGradientDescentState_SetX.argtypes = [c_void_p, POINTER(c_double)]

pLibDTSG.tsgGradientDescent_AdaptProj.argtypes = [type_optim_obj_fn_single, type_optim_grad_fn_single, type_optim_proj_fn_single,
                                                  c_double, c_double, c_int, c_double, c_void_p, POINTER(c_int)]
pLibDTSG.tsgGradientDescent_AdaptProj.restype = optimization_status
pLibDTSG.tsgGradientDescent_Adapt.argtypes = [type_optim_obj_fn_single, type_optim_grad_fn_single, c_double, c_double, c_int,
                                              c_double, c_void_p, POINTER(c_int)]
pLibDTSG.tsgGradientDescent_Adapt.restype = optimization_status
pLibDTSG.tsgGradientDescent_Const.argtypes = [type_optim_grad_fn_single, c_double, c_int, c_double, c_void_p, POINTER(c_int)]
pLibDTSG.tsgGradientDescent_Const.restype = optimization_status

class GradientDescentState:
    '''
    Wrapper to class TasOptimization::GradientDescentState
    '''
    def __init__(self, aX0, fInitialStepsize):
        '''
        Constructs a new state with a given number of dimensions and particles.

        aX0              : a one-dimensional NumPy array representing the starting point.
        fInitialStepsize : initial adaptive stepsize.
        '''
        self.TasmanianGradientDescentState = True
        if (len(aX0.shape) != 1):
            raise InputError("aX0", "aX0 should be one-dimensional")
        self.pStatePntr = c_void_p(pLibDTSG.tsgGradientDescentState_Construct(aX0.shape[0], as_ctypes(aX0), fInitialStepsize))

    def __del__(self):
        '''
        Deletes an instance of the particle swarm state.
        '''
        pLibDTSG.tsgGradientDescentState_Destruct(self.pStatePntr)

    def getNumDimensions(self):
        '''
        Return the number of dimensions.
        '''
        return pLibDTSG.tsgGradientDescentState_GetNumDimensions(self.pStatePntr)

    def getAdaptiveStepsize(self):
        '''
        Return the adaptive stepsize.
        '''
        return pLibDTSG.tsgGradientDescentState_GetAdaptiveStepsize(self.pStatePntr)

    def getX(self):
        '''
        Return the current candidate point.
        '''
        iNumDims = self.getNumDimensions()
        aResult = np.zeros((iNumDims,), np.float64)
        pLibDTSG.tsgGradientDescentState_GetX(self.pStatePntr, as_ctypes(aResult))
        return aResult

    def setAdaptiveStepsize(self, fNewStepsize):
        '''
        Set new adaptive stepsize

        fNewStepsize : a double representing the new adaptive stepsize.
        '''
        pLibDTSG.tsgGradientDescentState_SetAdaptiveStepsize(self.pStatePntr, fNewStepsize)

    def setX(self, aXNew):
        '''
        Set new particle positions from a NumPy array.

        aXNew : a one-dimensional NumPy array representing the new candidate point.
        '''
        iNumDims = self.getNumDimensions()
        if (len(aXNew.shape) != 1):
            raise InputError("aXNew", "aXNew should be one-dimensional")
        if (aXNew.shape[0] != iNumDims):
            raise InputError("aXNew", "aXNew.shape[0] should match the number of dimensions")
        pLibDTSG.tsgGradientDescentState_SetX(self.pStatePntr, as_ctypes(aXNew))

# Helper methods for converting Python lambdas to C lambdas.
def convert_py_obj_fn_single(pObjectiveFunction, sCallerName):
    def cpp_obj_fn(num_dim, x_single_ptr, err_arr):
        err_arr[0] = 1
        aX = np.ctypeslib.as_array(x_single_ptr, (num_dim, ))
        aResult = pObjectiveFunction(aX)
        if not isinstance(aResult, float):
            print("ERROR: incorrect output from the objective function given to " + sCallerName +
                  "(), should be a 'float' but received '" + type(aResult).__name__ + "'")
            return 0
        err_arr[0] = 0
        return aResult
    return cpp_obj_fn

def convert_py_grad_fn_single(pGradientFunction, sCallerName):
    def cpp_grad_fn(num_dim, x_single_ptr, grad_ptr, err_arr):
        err_arr[0] = 1
        aX = np.ctypeslib.as_array(x_single_ptr, (num_dim,))
        aResult = pGradientFunction(aX)
        if aResult.shape != (num_dim, ):
            print("ERROR: incorrect output from the gradient function given to " + sCallerName +
                  "(), should be a NumPy array with shape (iNumDims,)")
            return
        aGrad = np.ctypeslib.as_array(grad_ptr, (num_dim,))
        aGrad[0:num_dim] = aResult[0:num_dim]
        err_arr[0] = 0
    return cpp_grad_fn

def convert_py_proj_fn_single(pProjectionFunction, sCallerName):
    def cpp_proj_fn(num_dim, x_single_ptr, proj_ptr, err_arr):
        err_arr[0] = 1
        aX = np.ctypeslib.as_array(x_single_ptr, (num_dim, ))
        aResult = pProjectionFunction(aX)
        if aResult.shape != (num_dim, ):
            print("ERROR: incorrect output from the projection function given to " + sCallerName +
                  "(), should be a NumPy array with shape (iNumDims,)")
            return
        aProj =  np.ctypeslib.as_array(proj_ptr, (num_dim,))
        aProj[0:num_dim] = aResult[0:num_dim]
        err_arr[0] = 0
    return cpp_proj_fn

# Main algorithms.
def AdaptiveProjectedGradientDescent(pObjectiveFunction, pGradientFunction, pProjectionFunction, fIncreaseCoeff, fDecreaseCoeff,
                                     iMaxIterations, fTolerance, oGradientDescentState):
    '''
    Wrapper around TasOptimization::GradientDescent().

    Runs at most iMaxIterations of the adaptive stepsize projected gradient descent algorithm on an input state oGradientDescentState
    to minimize the function pObjectiveFunction with gradient pGradientFunction and projection function pProjectionFunction. The
    parameters fIncreaseCoeff and fDecreaseCoeff control the evolution of the algorithm. The parameter fTolerance is a first-order
    stationarity tolerance used for early termination (see the C++ documentation for more details).

    Returns a dictionary containing the results of the run.

    pObjectiveFunction    : a Python lambda representing the objective function; it should take in a one-dimensional NumPy array
                            (x_single) and produce a double representing be the result of evaluating the objective function on
                            x_single; it is expected that x_single.shape[0] = .getNumDimensions().
    pGradientFunction     : a Python lambda representing the gradient of the objective function; it should take in a
                            one-dimensional NumPy array (x_single) and produce a one-dimensional NumPy representing the gradient
                            (grad) at x_single; it is expected that x_single.shape[0] = grad.shape[0] = .getNumDimensions().
    pProjectionFunction   : a Python lambda representing the domain projection function; it should take in a one-dimensional NumPy
                            array (x_single) and produce a one-dimensional NumPy representing the projection (proj) of x_single
                            onto the domain; it is expected that x_single.shape[0] = proj.shape[0] = .getNumDimensions().
    fIncreaseCoeff        : a positive double controlling how fast the stepsize is increased; it is expected to be greater than 1.
    fDecreaseCoeff        : a positive double controlling how fast the stepsize is decreased; it is expected to be greater than 1.
    iMaxIterations        : a positive integer representing the maximum number iterations the algorithm is run.
    fTolerance            : a positive double representing the stationarity tolerance for early tolerance.
    oGradientDescentState : an instance of the Python GradientDescentState class; it will contain the results of applying the algorithm.
    '''
    cpp_obj_fn = convert_py_obj_fn_single(pObjectiveFunction, "GradientDescent")
    cpp_grad_fn = convert_py_grad_fn_single(pGradientFunction, "GradientDescent")
    cpp_proj_fn = convert_py_proj_fn_single(pProjectionFunction, "GradientDescent")
    pErrorCode = (c_int * 1)()
    oStatus = pLibDTSG.tsgGradientDescent_AdaptProj(type_optim_obj_fn_single(cpp_obj_fn), type_optim_grad_fn_single(cpp_grad_fn),
                                                    type_optim_proj_fn_single(cpp_proj_fn), c_double(fIncreaseCoeff), c_double(fDecreaseCoeff),
                                                    c_int(iMaxIterations), c_double(fTolerance), oGradientDescentState.pStatePntr, pErrorCode)
    if pErrorCode[0] != 0:
        raise InputError("GradientDescent", "An error occurred during the call to Tasmanian.")
    return clean_status(oStatus)

def AdaptiveGradientDescent(pObjectiveFunction, pGradientFunction, fIncreaseCoeff, fDecreaseCoeff, iMaxIterations, fTolerance,
                            oGradientDescentState):
    '''
    Wrapper around TasOptimization::GradientDescent().

    Runs at most iMaxIterations of the adaptive stepsize (unconstrained) gradient descent algorithm on an input state oGradientDescentState
    to minimize the function pObjectiveFunction with gradient pGradientFunction and projection function pProjectionFunction. The
    parameters fIncreaseCoeff and fDecreaseCoeff control the evolution of the algorithm. The parameter fTolerance is a first-order
    stationarity tolerance used for early termination (see the C++ documentation for more details).

    Returns a dictionary containing the results of the run.

    pObjectiveFunction    : a Python lambda representing the objective function; it should take in a one-dimensional NumPy array
                            (x_single) and produce a double representing be the result of evaluating the objective function on
                            x_single; it is expected that x_single.shape[0] = .getNumDimensions().
    pGradientFunction     : a Python lambda representing the gradient of the objective function; it should take in a
                            one-dimensional NumPy array (x_single) and produce a one-dimensional NumPy representing the gradient
                            (grad) at x_single; it is expected that x_single.shape[0] = grad.shape[0] = .getNumDimensions().
    fIncreaseCoeff        : a positive double controlling how fast the stepsize is increased; it is expected to be greater than 1.
    fDecreaseCoeff        : a positive double controlling how fast the stepsize is decreased; it is expected to be greater than 1.
    iMaxIterations        : a positive integer representing the maximum number iterations the algorithm is run.
    fTolerance            : a positive double representing the stationarity tolerance for early tolerance.
    oGradientDescentState : an instance of the Python GradientDescentState class; it will contain the results of applying the algorithm.
    '''
    cpp_obj_fn = convert_py_obj_fn_single(pObjectiveFunction, "GradientDescent")
    cpp_grad_fn = convert_py_grad_fn_single(pGradientFunction, "GradientDescent")
    pErrorCode = (c_int * 1)()
    oStatus = pLibDTSG.tsgGradientDescent_Adapt(type_optim_obj_fn_single(cpp_obj_fn), type_optim_grad_fn_single(cpp_grad_fn),
                                                c_double(fIncreaseCoeff), c_double(fDecreaseCoeff), c_int(iMaxIterations),
                                                c_double(fTolerance), oGradientDescentState.pStatePntr, pErrorCode)
    if pErrorCode[0] != 0:
        raise InputError("GradientDescent", "An error occurred during the call to Tasmanian.")
    return clean_status(oStatus)

def GradientDescent(pGradientFunction, fStepsize, iMaxIterations, fTolerance, oGradientDescentState):
    '''
    Wrapper around TasOptimization::GradientDescent().

    Runs at most iMaxIterations of the constant stepsize (unconstrained) gradient descent algorithm on an input state oGradientDescentState
    to minimize the function pObjectiveFunction with gradient pGradientFunction and projection function pProjectionFunction.
    The parameter fTolerance is a first-order stationarity tolerance used for early termination (see the C++ documentation for
    more details).

    Returns a dictionary containing the results of the run.

    pGradientFunction     : a Python lambda representing the gradient of the objective function; it should take in a
                            one-dimensional NumPy array (x_single) and produce a one-dimensional NumPy representing the gradient
                            (grad) at x_single; it is expected that x_single.shape[0] = grad.shape[0] = .getNumDimensions().
    fStepsize             : a positive double representing the stepsize.
    iMaxIterations        : a positive integer representing the maximum number iterations the algorithm is run.
    fTolerance            : a positive double representing the stationarity tolerance for early tolerance.
    oGradientDescentState : an instance of the Python GradientDescentState class; it will contain the results of applying the algorithm.
    '''
    cpp_grad_fn = convert_py_grad_fn_single(pGradientFunction, "GradientDescent")
    pErrorCode = (c_int * 1)()
    oStatus = pLibDTSG.tsgGradientDescent_Const(type_optim_grad_fn_single(cpp_grad_fn), c_double(fStepsize), c_int(iMaxIterations),
                                                c_double(fTolerance), oGradientDescentState.pStatePntr, pErrorCode)
    if pErrorCode[0] != 0:
        raise InputError("GradientDescent", "An error occurred during the call to Tasmanian.")
    return clean_status(oStatus)

