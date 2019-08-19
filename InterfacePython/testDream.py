import unittest
import Tasmanian
from Tasmanian import DREAM
import os
import numpy as np

import testCommon

ttc = testCommon.TestTasCommon()

class TestTasClass(unittest.TestCase):
    '''
    Miscelaneous tests that don't quite fit in other categories.
    '''
    def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")

    def testNothing(self):
        pass

    def checkGaussian3D(self):
        iNumDimensions = 3
        iNumSamples = 1000
        iNumChains = 20
        iNumIterations = int(iNumSamples / iNumChains) + 2
        iNumBurnup = iNumIterations

        state = DREAM.State(iNumChains, iNumDimensions)
        state.setState(DREAM.tsgGenGaussianSamples([2.0, 2.0, 2.0], [3.0, 3.0, 3.0], iNumChains, DREAM.RandomGenerator("minstd_rand", 42)))

        DREAM.Sample(iNumBurnup, iNumSamples,
                     lambda x : np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ), # Gaussian mean 2.0 variance 9.0
                     DREAM.Domain("unbounded"),
                     state,
                     DREAM.IndependentUpdate("gaussian", 3.0),
                     DREAM.DifferentialUpdate(0),
                     DREAM.RandomGenerator("minstd_rand", 42))

        aMean, aVar = state.getHistoryMeanVariance()
        self.assertTrue(np.max(np.abs(aMean - 2.0)) < 0.2, "Failed the mean test for Gaussian 3D")
        self.assertTrue(np.max(np.abs(aVar - 9.0)) < 1.0, "Failed the variance test for Gaussian 3D")

        state = DREAM.State(iNumChains, iNumDimensions)
        state.setState(DREAM.tsgGenGaussianSamples([1.5, 2.0, 2.5], [0.5, 1.0, 2.0], iNumChains, DREAM.RandomGenerator("minstd_rand", 42)))

        likely = DREAM.LikelihoodGaussAnisotropic(np.array([0.25, 1.0, 4.0]), np.array([1.5, 2.0, 2.5]), 1)

        DREAM.Sample(iNumBurnup, iNumSamples,
                     lambda x : likely.getLikelihood(DREAM.typeLogform, x), # Gaussian mean 2.0 variance 9.0
                     DREAM.Domain("unbounded"),
                     state,
                     DREAM.IndependentUpdate("gaussian", 3.0),
                     DREAM.DifferentialUpdate(0),
                     DREAM.RandomGenerator("minstd_rand", 72),
                     DREAM.typeLogform)

        aMean, aVar = state.getHistoryMeanVariance()
        self.assertTrue(np.max(np.abs(aMean - np.array([1.5, 2.0, 2.5]))) < 0.4, "Failed the mean test for anisotropic Gaussian 3D")
        self.assertTrue(np.max(np.abs(aVar - np.array([0.25, 1.0, 4.0]))) < 1.0, "Failed the variance test for anisotropic Gaussian 3D")


    def checkAnalytic(self):
        self.checkGaussian3D()

    def performDreamTests(self):
        self.checkAnalytic()

if (__name__ == "__main__"):
    tester = TestTasClass()
    tester.performDreamTests()
