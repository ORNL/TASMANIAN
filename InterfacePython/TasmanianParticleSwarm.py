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

from ctypes import c_char_p, c_int, c_bool, c_double, c_void_p, POINTER, cdll, cast, CFUNCTYPE
from numpy.ctypeslib import as_ctypes, ndpointer
import numpy as np
import sys

from TasmanianConfig import __path_libdream__
from TasmanianConfig import TasmanianInputError as InputError

import TasmanianDREAM as DREAM

c_double_p        = POINTER(c_double)
np_double_arr     = ndpointer(c_double, flags="C_CONTIGUOUS")
type_optim_obj_fn = CFUNCTYPE(None, c_int, c_int, c_double_p, c_double_p)
type_optim_dom_fn = CFUNCTYPE(c_bool, c_int, c_double_p)
type_dream_random = CFUNCTYPE(c_double)

pLibDTSG = cdll.LoadLibrary(__path_libdream__)

pLibDTSG.tsgParticleSwarmState_Construct.restype = c_void_p
pLibDTSG.tsgParticleSwarmState_Construct.argtypes = [c_int, c_int]

pLibDTSG.tsgParticleSwarmState_Destruct.argtypes = [c_void_p]

pLibDTSG.tsgParticleSwarmState_GetNumDimensions.restype = c_int
pLibDTSG.tsgParticleSwarmState_GetNumDimensions.argtype = [c_void_p]
pLibDTSG.tsgParticleSwarmState_GetNumParticles.restype = c_int
pLibDTSG.tsgParticleSwarmState_GetNumParticles.argtype = [c_void_p]
pLibDTSG.tsgParticleSwarmState_GetParticlePositions.argtypes = [c_void_p, np_double_arr]
pLibDTSG.tsgParticleSwarmState_GetParticleVelocities.argtypes = [c_void_p, np_double_arr]
pLibDTSG.tsgParticleSwarmState_GetBestParticlePositions.argtypes = [c_void_p, np_double_arr]
pLibDTSG.tsgParticleSwarmState_GetBestPosition.argtypes = [c_void_p, np_double_arr]

pLibDTSG.tsgParticleSwarmState_SetParticlePositions.argtypes = [c_void_p, np_double_arr]
pLibDTSG.tsgParticleSwarmState_SetParticleVelocities.argtypes = [c_void_p, np_double_arr]
pLibDTSG.tsgParticleSwarmState_SetBestParticlePositions.argtypes = [c_void_p, np_double_arr]

pLibDTSG.tsgParticleSwarmState_ClearBestParticles.argtypes = [c_void_p]
pLibDTSG.tsgParticleSwarmState_ClearCache.argtypes = [c_void_p]
pLibDTSG.tsgParticleSwarmState_InitializeParticlesInsideBox.argtypes = [c_void_p, np_double_arr, np_double_arr, c_char_p, c_int, type_dream_random]
pLibDTSG.tsgParticleSwarm.argtypes = [type_optim_obj_fn, c_int, type_optim_dom_fn, c_void_p, c_double, c_double, c_double, c_char_p, c_int, type_dream_random]

class ParticleSwarmState:
    '''
    Wrapper to class TasOptimization::ParticleSwarmState
    '''
    def __init__(self, iNumDimensions, iNumParticles):
        '''
        Constructs a new state with a given number of dimensions and particles.

        iNumDimensions : positive integer indicating the number of dimensions.
        iNumParticles  : positive integer indicating the number of particles.
        '''
        self.TasmanianParticleSwarmState = True
        self.pStatePntr = c_void_p(pLibDTSG.tsgParticleSwarmState_Construct(iNumDimensions, iNumParticles))

    def __del__(self):
        '''
        Deletes an instance of the particle swarm state.
        '''
        pLibDTSG.tsgParticleSwarmState_Destruct(self.pStatePntr)

    def getNumDimensions(self):
        '''
        Return the number of dimensions.
        '''
        return pLibDTSG.tsgParticleSwarmState_GetNumDimensions(self.pStatePntr)

    def getNumParticles(self):
        '''
        Return the number of particles.
        '''
        return pLibDTSG.tsgParticleSwarmState_GetNumParticles(self.pStatePntr)

    def getParticlePositions(self):
        '''
        Return the particle positions as a NumPy array.
        '''
        iNumDims = self.getNumDimensions()
        iNumPart = self.getNumParticles()
        aResult = np.zeros((iNumDims * iNumPart,), np.float64)
        pLibDTSG.tsgParticleSwarmState_GetParticlePositions(self.pStatePntr, aResult)
        return aResult.reshape((iNumPart, iNumDims))

    def getParticleVelocities(self):
        '''
        Return the particle velocities as a NumPy array.
        '''
        iNumDims = self.getNumDimensions()
        iNumPart = self.getNumParticles()
        aResult = np.zeros((iNumDims * iNumPart,), np.float64)
        pLibDTSG.tsgParticleSwarmState_GetParticleVelocities(self.pStatePntr, aResult)
        return aResult.reshape((iNumPart, iNumDims))

    def getBestParticlePositions(self):
        '''
        Return the best particle positions as a NumPy array.
        '''
        iNumDims = self.getNumDimensions()
        iNumPart = self.getNumParticles()
        aResult = np.zeros((iNumDims * (iNumPart + 1),), np.float64)
        pLibDTSG.tsgParticleSwarmState_GetBestParticlePositions(self.pStatePntr, aResult)
        return aResult.reshape((iNumPart + 1, iNumDims))

    def getBestPosition(self):
        '''
        Return the best particle position as a NumPy array.
        '''
        iNumDims = self.getNumDimensions()
        aResult = np.zeros((iNumDims,), np.float64)
        pLibDTSG.tsgParticleSwarmState_GetBestPosition(self.pStatePntr, aResult)
        return aResult

    def setParticlePositions(self, llfNewPPosns):
        '''
        Set new particle positions from a NumPy array.

        llfNewPPosns : a two-dimensional numpy.ndarray with
            .shape[0] = .getNumParticles()
            .shape[1] = .getNumDimensions()
        '''
        iNumPart = self.getNumParticles()
        iNumDims = self.getNumDimensions()
        if (llfNewPPosns.shape[0] != iNumPart):
            raise InputError("llfNewPPosns", "llfNewPPosns.shape[0] should match the number of particles")
        if (llfNewPPosns.shape[1] != iNumDims):
            raise InputError("llfNewPPosns", "llfNewPPosns.shape[1] should match the number of dimensions")
        pLibDTSG.tsgParticleSwarmState_SetParticlePositions(self.pStatePntr, llfNewPPosns)

    def setParticleVelocities(self, llfNewPVelcs):
        '''
        Set new particle velocities from a NumPy array.

        llfNewPVelcs : a two-dimensional numpy.ndarray with
            .shape[0] = .getNumParticles()
            .shape[1] = .getNumDimensions()
        '''
        iNumPart = self.getNumParticles()
        iNumDims = self.getNumDimensions()
        if (llfNewPVelcs.shape[0] != iNumPart):
            raise InputError("llfNewPVelcs", "llfNewPVelcs.shape[0] should match the number of particles")
        if (llfNewPVelcs.shape[1] != iNumDims):
            raise InputError("llfNewPVelcs", "llfNewPVelcs.shape[1] should match the number of dimensions")
        pLibDTSG.tsgParticleSwarmState_SetParticleVelocities(self.pStatePntr, llfNewPVelcs)

    def setBestParticlePositions(self, llfNewBPPosns):
        '''
        Set new best particle positions from a NumPy array.

        llfNewBPPosns : a two-dimensional numpy.ndarray with
            .shape[0] = .getNumParticles() + 1
            .shape[1] = .getNumDimensions()
        '''
        iNumPart = self.getNumParticles()
        iNumDims = self.getNumDimensions()
        if (llfNewBPPosns.shape[0] != iNumPart + 1):
            raise InputError("llfNewBPPosns", "llfNewBPPosns.shape[0] should match the number of particles + 1")
        if (llfNewBPPosns.shape[1] != iNumDims):
            raise InputError("llfNewBPPosns", "llfNewBPPosns.shape[1] should match the number of dimensions")
        llfNewBPPosns.resize(((iNumPart + 1) * iNumDims,))
        pLibDTSG.tsgParticleSwarmState_SetBestParticlePositions(self.pStatePntr, llfNewBPPosns)

    def clearBestParticles(self):
        '''
        Clears the vector of best particles.
        '''
        pLibDTSG.tsgParticleSwarmState_ClearBestParticles(self.pStatePntr)

    def clearCache(self):
        '''
        Clears the internal cache of the managed ParticleSwarm C++ object.
        '''
        pLibDTSG.tsgParticleSwarmState_ClearCache(self.pStatePntr)

    def initializeParticlesInsideBox(self, lfBoxLower, lfBoxUpper, random01 = DREAM.RandomGenerator(sType = "default")):
        '''
        Initialize the particle swarm state with randomized particles (determined by gen_random01) inside a box bounded by the
        parameters box_lower and box_upper. The i-th entry in these parameters respectively represent the lower and upper
        bounds on the particle position values for the i-th dimension. The parameter random01 controls the distribution of
        the particles.

        lfBoxLower : a one-dimensional NumPy array with lfBoxLower.shape[0] = .getNumDimensions()
        lfBoxUpper : a one-dimensional NumPy array with lfBoxLower.shape[0] = .getNumDimensions()
        random01   : a DREAM.RandomGenerator instance that produces floats in [0,1]
        '''
        iNumDims = self.getNumDimensions()
        if (lfBoxLower.shape[0] != iNumDims):
            raise InputError("lfBoxLower", "lfBoxLower.shape[0] should match the number of dimensions")
        if (lfBoxUpper.shape[0] != iNumDims):
            raise InputError("lfBoxUpper", "lfBoxUpper.shape[0] should match the number of dimensions")
        pLibDTSG.tsgParticleSwarmState_InitializeParticlesInsideBox(self.pStatePntr, lfBoxLower, lfBoxUpper, c_char_p(random01.sType),
                                                                    c_int(random01.iSeed), type_dream_random(random01.pCallable))

def ParticleSwarm(pObjectiveFunction, iNumIterations, pInside, oParticleSwarmState, fInertiaWeight, fCognitiveCoeff,
                  fSocialCoeff, random01 = DREAM.RandomGenerator(sType = "default")):
    '''
    Wrapper around TasOptimization::ParticleSwarm().

    Runs iNumIterations of the particle swarm algorithm on an input state cParticleSwarmState to minimize the function
    pObjectiveFunction over a domain specified by pInside. The parameters fInertiaWeight, fCognitiveCoeff, and fSocialCoeff
    control the evolution of the algorithm. The parameter random01 controls the distribution of the particles.

    pObjectiveFunction  : a Python lambda representing the objective function; it should take in two 2D NumPy arrays (x_batch, fval_batch)
                          and produce no outputs; when it is called, it should write the result of applying the objective function
                          to every row of x_batch to fval_batch; it is expected that x_batch.shape[1] = .getNumDimensions().
    iNumIterations      : a positive integer representing the number iterations the algorithm is run.
    pInside             : a Python lambda representing the function domain; it should take in one 2D NumPy array (x) and produce an
                          integer (isInside); when called, it should return True if x is in the domain and False otherwise;
                          it is expected that x.shape[0] = .getNumDimensions().
    oParticleSwarmState : an instance of the Python ParticleSwarmState class; it will contain the results of applying the algorithm.
    fInertiaWeight      : a double that controls the speed of the particles; should be in (0,1).
    fCognitiveCoeff     : a double that controls how much a particle favors its own trajectory; usually in [1,3].
    fSocialCoeff        : a double that controls how much a particle favors the swarm's trajectories; usually in [1,3].
    random01            : a DREAM.RandomGenerator instance that produces floats in [0,1]

    NOTE: np.frombuffer() is used below in order to avoid copying the data owned by the C pointers into the generated NumPy arrays.
    '''
    def cpp_obj_fn(num_dim, num_batch, x_batch_ptr, fval_ptr):
        py_x_batch_ptr = cast(x_batch_ptr, POINTER(c_double * num_batch * num_dim))
        np_x_batch = np.frombuffer(py_x_batch_ptr.contents)
        np_x_batch.resize((num_batch, num_dim))
        py_fval_ptr = cast(fval_ptr, POINTER(c_double * num_batch))
        np_fval = np.frombuffer(py_fval_ptr.contents)
        pObjectiveFunction(np_x_batch, np_fval)

    def cpp_dom_fn(num_dim, x_ptr):
        py_x_ptr = cast(x_ptr, POINTER(c_double * num_dim))
        np_x = np.frombuffer(py_x_ptr.contents)
        return pInside(np_x)

    pLibDTSG.tsgParticleSwarm(type_optim_obj_fn(cpp_obj_fn), c_int(iNumIterations), type_optim_dom_fn(cpp_dom_fn),
                              oParticleSwarmState.pStatePntr, c_double(fInertiaWeight), c_double(fCognitiveCoeff),
                              c_double(fSocialCoeff), c_char_p(random01.sType), c_int(random01.iSeed),
                              type_dream_random(random01.pCallable))
