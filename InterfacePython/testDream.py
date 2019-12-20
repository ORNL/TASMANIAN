import unittest
import Tasmanian
DREAM = Tasmanian.DREAM
import os
import numpy as np

import testCommon

ttc = testCommon.TestTasCommon()

class TestTasClass(unittest.TestCase):
    '''
    Tests for the Tasmanian DREAM python bindings module
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
        state.setState(DREAM.genGaussianSamples([2.0, 2.0, 2.0], [3.0, 3.0, 3.0], iNumChains, DREAM.RandomGenerator("minstd_rand", 42)))

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
        state.setState(DREAM.genGaussianSamples([1.5, 2.0, 2.5], [0.5, 1.0, 2.0], iNumChains, DREAM.RandomGenerator("minstd_rand", 42)))

        likely = DREAM.LikelihoodGaussAnisotropic(np.array([0.25, 1.0, 4.0]), np.array([1.5, 2.0, 2.5]), 1)

        DREAM.Sample(iNumBurnup, iNumSamples,
                     lambda x : likely.getLikelihood(DREAM.typeLogform, x), # Gaussian mean 2.0 variance 9.0
                     DREAM.Domain("unbounded"),
                     state,
                     DREAM.IndependentUpdate("gaussian", 3.0),
                     DREAM.DifferentialUpdate(97),
                     DREAM.RandomGenerator("minstd_rand", 72),
                     DREAM.typeLogform)

        aMean, aVar = state.getHistoryMeanVariance()
        self.assertTrue(np.max(np.abs(aMean - np.array([1.5, 2.0, 2.5]))) < 0.4, "Failed the mean test for anisotropic Gaussian 3D")
        self.assertTrue(np.max(np.abs(aVar - np.array([0.25, 1.0, 4.0]))) < 1.0, "Failed the variance test for anisotropic Gaussian 3D")

    def checkGrid(self):
        iNumDimensions = 2
        iNumSamples = 1000
        iNumChains = 20
        iNumIterations = int(iNumSamples / iNumChains) + 2
        iNumBurnup = iNumIterations

        grid = Tasmanian.SparseGrid()
        grid.makeSequenceGrid(2, 1, 24, "level", "rleja")
        Tasmanian.loadNeededPoints(lambda x, tid : np.ones((1,)) * np.exp( -0.5 * np.sum((x - 0.3)**2) / 0.01 ), grid, 1)

        state = DREAM.State(iNumChains, grid)
        state.setState(DREAM.genGaussianSamples([0.3, 0.3], [0.1, 0.1], iNumChains, DREAM.RandomGenerator("minstd_rand", 55)))

        DREAM.Sample(iNumBurnup, iNumSamples,
                     DREAM.PosteriorEmbeddedLikelihood(grid, "uniform"),
                     DREAM.Domain("hypercube", np.array((0.0, 0.0)), np.array((1.0, 1.0))),
                     state,
                     DREAM.IndependentUpdate("uniform", 0.1),
                     DREAM.DifferentialUpdate(98),
                     DREAM.RandomGenerator("minstd_rand", 33))

        aMean, aVar = state.getHistoryMeanVariance()
        self.assertTrue(np.max(np.abs(aMean - np.array([0.3, 0.3]))) < 0.05, "Failed the mean test for grid likelihood")
        self.assertTrue(np.max(np.abs(aVar - np.array([0.01, 0.01]))) < 0.005, "Failed the variance test for grid likelihood")

    def checkModel(self):
        iNumDimensions = 2
        iNumSamples = 1000
        iNumChains = 20
        iNumIterations = int(iNumSamples / iNumChains) + 2
        iNumBurnup = iNumIterations

        likely = DREAM.LikelihoodGaussIsotropic(0.01, np.array([1.0, 2.0]))

        state = DREAM.State(iNumChains, iNumDimensions)
        state.setState(DREAM.genUniformSamples([-1.0, -1.0], [1.0, 1.0], iNumChains, DREAM.RandomGenerator("minstd_rand", 67)))

        DREAM.Sample(iNumBurnup, iNumSamples,
                     DREAM.Posterior(lambda x : np.matmul(np.cos(x), np.array([[1.0, 0.0], [0.0, 2.0]])),
                                     likely, "uniform", DREAM.typeLogform),
                     DREAM.Domain("hypercube", np.array((0.0, 0.0)), np.array((1.0, 1.0))),
                     state,
                     DREAM.IndependentUpdate("uniform", 0.1),
                     DREAM.DifferentialUpdate(98),
                     DREAM.RandomGenerator("minstd_rand", 88),
                     DREAM.typeLogform)

        aMode = state.getApproximateMode()
        self.assertTrue(np.max(np.abs(aMode)) < 0.01, "Failed the mode test for custom model and Tasmanian likelihood")

    def performDreamTests(self):
        self.checkGaussian3D() # analytic model/likelihood
        self.checkGrid()
        self.checkModel() # custom model

if __name__ == "__main__":
    tester = TestTasClass()
    tester.performDreamTests()
