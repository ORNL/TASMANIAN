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

import unittest
import Tasmanian
DREAM = Tasmanian.DREAM
Opt = Tasmanian.Optimization
import os
import numpy as np

import testCommon

ttc = testCommon.TestTasCommon()

class TestTasClass(unittest.TestCase):
    '''
    Tests for the Tasmanian Optimization python bindings module
    '''
    def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")

    def testNothing(self):
        pass

    def checkParticleSwarmState(self):
        iNumDimensions = 2
        iNumParticles = 5
        # Initialization tests.
        state = Opt.ParticleSwarmState(iNumDimensions, iNumParticles)
        self.assertEqual(iNumDimensions, state.getNumDimensions())
        self.assertEqual(iNumParticles, state.getNumParticles())
        self.assertFalse(state.isPositionInitialized())
        self.assertFalse(state.isVelocityInitialized())
        self.assertFalse(state.isBestPositionInitialized())
        self.assertFalse(state.isCacheInitialized())
        # Loading/Unloading tests.
        pp = np.arange(iNumParticles * iNumDimensions, dtype=np.float64).reshape(iNumParticles, iNumDimensions)
        state.setParticlePositions(pp)
        self.assertTrue(state.isPositionInitialized())
        np.testing.assert_array_equal(pp, state.getParticlePositions())
        pv =  np.arange(iNumParticles * iNumDimensions, dtype=np.float64).reshape(iNumParticles, iNumDimensions)
        state.setParticleVelocities(pv)
        np.testing.assert_array_equal(pv, state.getParticleVelocities())
        bpp =  np.arange((iNumParticles+1) * iNumDimensions, dtype=np.float64).reshape(iNumParticles+1, iNumDimensions)
        state.setBestParticlePositions(bpp)
        self.assertTrue(state.isBestPositionInitialized())
        np.testing.assert_array_equal(bpp, state.getBestParticlePositions())
        state.clearBestParticles()
        self.assertFalse(state.isBestPositionInitialized())
        np.testing.assert_array_equal(np.zeros([iNumParticles + 1, iNumDimensions]), state.getBestParticlePositions())
        # Generation tests.
        state = Opt.ParticleSwarmState(iNumDimensions, iNumParticles)
        state.initializeParticlesInsideBox(np.array([-1.0, 1.0]), np.array([2.0, 3.0]),
                                           random01=DREAM.RandomGenerator(sType="default", iSeed=777))
        self.assertTrue(state.isPositionInitialized())
        self.assertTrue(state.isVelocityInitialized())
        pp = state.getParticlePositions()
        self.assertTrue(np.all((-1.0 <= pp[:,0]) & (pp[:,0] <= 2.0)))
        self.assertTrue(np.all((1.0 <= pp[:,1]) & (pp[:,1] <= 3.0)))

    def checkParticleSwarm(self):
        iNumDimensions = 2
        iNumParticles = 50
        state = Opt.ParticleSwarmState(iNumDimensions, iNumParticles)
        state.initializeParticlesInsideBox(np.array([-3.0, -2.0]), np.array([3.0, 2.0]),
                                           random01=DREAM.RandomGenerator(sType="default", iSeed=777))
        # Six-hump-camel function.
        shc = lambda x : (4 - 2.1*x[0]*x[0] + x[0]*x[0]*x[0]*x[0]/3)*x[0]*x[0] + x[0]*x[1] + (-4.0 + 4.0*x[1]*x[1])*x[1]*x[1]
        f = lambda x_batch : np.apply_along_axis(shc, 1, x_batch)
        inside = lambda x : bool((-3 <= x[0]) and (x[0] <= 3) and (-2 <= x[1]) and (x[1] <= 2))
        # Main call + tests.
        Opt.ParticleSwarm(f, inside, 0.5, 2, 2, 1, state, random01=DREAM.RandomGenerator(sType="default", iSeed=777))
        self.assertTrue(state.isCacheInitialized())
        state.clearCache()
        self.assertFalse(state.isCacheInitialized())
        iNumIterations = 200
        Opt.ParticleSwarm(f, inside, 0.5, 2, 2, iNumIterations, state)
        self.assertTrue(np.allclose(np.array([-0.08984201368301331, +0.7126564032704135]), state.getBestPosition()) or
                        np.allclose(np.array([+0.08984201368301331, -0.7126564032704135]), state.getBestPosition()) )

    def checkGradientDescentState(self):
        x0 = np.array([0.0, 1.0])
        stepsize = 2.0;
        # Initialization tests.
        state = Opt.GradientDescentState(x0, stepsize)
        self.assertEqual(len(x0), state.getNumDimensions())
        np.testing.assert_array_equal(x0, state.getX())
        self.assertEqual(stepsize, state.getAdaptiveStepsize())
        # Loading/Unloading tests.
        x1 = np.array([3.0, 4.0])
        state.setX(x1)
        np.testing.assert_array_equal(x1, state.getX())
        state.setAdaptiveStepsize(5.0)
        self.assertEqual(5.0, state.getAdaptiveStepsize())

    def checkGradientDescent(self):
        func = lambda x : 2.0 * x[0] * x[0] + x[1] * x[1] / 2.0
        grad = lambda x : np.array([4.0 * x[0], x[1]])
        lambda0 = 3.0
        xBar = np.array([0.0, 0.0])
        xBarConstr = np.array([-0.5, 0.5])
        tol = 1E-6
        iter_limit = 100
        # Main call + tests.
        state = Opt.GradientDescentState(np.array([-1.0, 1.0]), lambda0)
        result = Opt.GradientDescent(grad, 1/4.0, iter_limit, tol, state)
        self.assertLessEqual(result['performed_iterations'], iter_limit)
        self.assertTrue(np.allclose(xBar, state.getX(), atol=tol * 10))
        self.assertEqual(lambda0, state.getAdaptiveStepsize())
        state = Opt.GradientDescentState(np.array([-1.0, 1.0]), lambda0)
        result = Opt.AdaptiveGradientDescent(func, grad, 1.25, 1.25, iter_limit, tol, state)
        self.assertLessEqual(result['performed_iterations'], iter_limit)
        self.assertTrue(np.allclose(xBar, state.getX(), atol=tol * 10))
        proj = lambda x : np.array([min(max(-3.0, x[0]), -0.5), min(max(0.5, x[1]), 3.0)])
        state = Opt.GradientDescentState(np.array([-2.0, 2.0]), lambda0)
        result = Opt.AdaptiveProjectedGradientDescent(func, grad, proj, 1.25, 1.25, iter_limit, tol, state)
        self.assertLessEqual(result['performed_iterations'], iter_limit)
        self.assertTrue(np.allclose(xBarConstr, state.getX(), atol=tol))

    def performOptTests(self):
        self.checkParticleSwarmState()
        self.checkParticleSwarm()
        self.checkGradientDescentState()
        self.checkGradientDescent()

if __name__ == "__main__":
    tester = TestTasClass()
    tester.performOptTests()
