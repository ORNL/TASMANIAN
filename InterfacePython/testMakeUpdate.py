import unittest
import Tasmanian
import TasmanianSG # needed to test all rules
import numpy as np
import os, sys
from ctypes import cdll

from random import shuffle

import testConfigureData as tdata # needed to keep track of version tests and ways to load the library
import testCommon

ttc = testCommon.TestTasCommon()

class TestTasClass(unittest.TestCase):
    '''
    Test make/update, check points, check library paths
    '''
    def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")

    def testNothing(self):
        pass

    def checkPathsVersions(self):
        '''
        Check different ways to specify the library path and check the
        version number.
        '''
        grid = Tasmanian.TasmanianSparseGrid()

        sVersion = "{0:1d}.{1:1d}".format(grid.getVersionMajor(), grid.getVersionMinor())
        self.assertEqual(sVersion, Tasmanian.__version__, "version mismatch")
        sLicense = grid.getLicense()
        self.assertEqual(sLicense, Tasmanian.__license__, "license mismatch")

        iVM = int(sVersion.split('.')[0])
        iVm = int(sVersion.split('.')[1])
        self.assertEqual(iVM, grid.getVersionMajor(), "version major mismatch")
        self.assertEqual(iVm, grid.getVersionMinor(), "version minor mismatch")

    def checkMakeAgainstKnown(self):
        '''
        Make grids and check against pen-and-paper answers, check weights too.
        '''
        grid = Tasmanian.makeGlobalGrid(1, 0, 4, 'level', 'gauss-hermite', [], 2.0)
        aW = grid.getQuadratureWeights()
        self.assertTrue(np.abs(sum(aW) - 0.5 * np.sqrt(np.pi)) < 1.E-14, "Gauss-Hermite Alpha")

        grid.makeGlobalGrid(2, 0, 2, 'level', 'leja', [2, 1])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0, 0.0]])
        aP = grid.getPoints()
        np.testing.assert_equal(aA, aP, 'Anisotropy Global not equal', True)

        grid.makeGlobalGrid(2, 0, 4, 'ipcurved', 'leja', [20, 10, 0, 7])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [0.0, np.sqrt(1.0/3.0)], [1.0, 0.0], [1.0, 1.0], [-1.0, 0.0]])
        aP = grid.getPoints()
        np.testing.assert_almost_equal(aA, aP, 12, 'Anisotropy Global not equal', True)

        grid.makeSequenceGrid(2, 1, 2, 'level', 'leja', [2, 1])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0, 0.0]])
        aP = grid.getPoints()
        np.testing.assert_equal(aA, aP, 'Anisotropy Sequence not equal', True)

        grid = Tasmanian.makeSequenceGrid(2, 1, 4, 'ipcurved', 'leja', [20, 10, 0, 7])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [0.0, np.sqrt(1.0/3.0)], [1.0, 0.0], [1.0, 1.0], [-1.0, 0.0]])
        aP = grid.getPoints()
        np.testing.assert_almost_equal(aA, aP, 12, 'Anisotropy Sequence not equal', True)

        grid = Tasmanian.makeFourierGrid(2, 1, 2, 'level', [1, 2])
        aA = np.array([[0.0, 0.0], [0.0, 1.0/3.0], [0.0, 2.0/3.0], [1.0/3.0, 0.0], [2.0/3.0, 0.0], [1.0/9.0, 0.0], [2.0/9.0, 0.0], [4.0/9.0, 0.0], [5.0/9.0, 0.0], [7.0/9.0, 0.0], [8.0/9.0, 0.0]])
        aP = grid.getPoints()
        np.testing.assert_almost_equal(aA, aP, 12, 'Anisotropy Fourier not equal', True)

        # this is a very important test, checks the curved rule and covers the non-lower-set index selection code
        grid.makeGlobalGrid(2, 1, 1, 'ipcurved', 'rleja', [10, 10, -21, -21])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, 0.5], [0.0, 0.25], [0.0, 0.75], [0.0, 0.125],
                       [1.0, 0.0], [1.0, 1.0], [1.0, 0.5], [1.0, 0.25], [1.0, 0.75], [1.0, 0.125],
                       [0.5, 0.0], [0.5, 1.0], [0.5, 0.5], [0.5, 0.25], [0.5, 0.75], [0.5, 0.125],
                       [0.25, 0.0], [0.25, 1.0], [0.25, 0.5], [0.25, 0.25], [0.25, 0.75],
                       [0.75, 0.0], [0.75, 1.0], [0.75, 0.5], [0.75, 0.25],
                       [0.125, 0.0], [0.125, 1.0], [0.125, 0.5]])
        aA = np.cos(np.pi * aA)
        aP = grid.getPoints()
        np.testing.assert_almost_equal(aA, aP, 12, 'Anisotropy heavily curved not equal', True) # 12 is the number of dec places

        # test weights
        aX = np.array([[0.0, -0.3], [-0.44, 0.7], [0.82, -0.01]])
        grid.makeGlobalGrid(2, 0, 4, 'level', 'chebyshev')
        aW = grid.getInterpolationWeightsBatch(aX)
        aX1 = aX[0,:].reshape([2,])
        aX2 = aX[1,:].reshape([2,])
        aX3 = aX[2,:].reshape([2,])
        aW1 = grid.getInterpolationWeights(aX1)
        aW2 = grid.getInterpolationWeights(aX2)
        aW3 = grid.getInterpolationWeights(aX3)
        aWW = np.row_stack([aW1, aW2, aW3])
        np.testing.assert_equal(aW, aWW, "Batch weights", True)

        grid.makeGlobalGrid(3, 0, 1, 'level', 'chebyshev')
        aTrans = np.array([[-1.0, -4.0], [2.0, 5.0], [-3.0, 5]])
        grid.setDomainTransform(aTrans)
        aA = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        aB = np.array([[-2.5, 3.5, 1.0], [-2.5, 3.5, -3.0], [-2.5, 3.5, 5.0], [-2.5,  2.0, 1.0], [-2.5, 5.0, 1.0], [-1.0, 3.5, 1.0], [-4.0, 3.5, 1.0]])
        np.testing.assert_equal(aB, grid.getPoints(), "Original equal", True)
        grid.clearDomainTransform()
        np.testing.assert_equal(aA, grid.getPoints(), "Tansformed equal", True)

        grid.makeGlobalGrid(3, 0, 1, 'level', 'fejer2')
        aA = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -0.707106781186548], [0.0, 0.0, 0.707106781186548], [0.0, -0.707106781186548, 0.0], [0.0, 0.707106781186548, 0.0], [-0.707106781186548, 0.0, 0.0], [0.707106781186548, 0.0, 0.0]])
        np.testing.assert_almost_equal(aA, grid.getPoints(), 14, "Original equal", True)
        grid.setConformalTransformASIN(np.array([4, 6, 3]))
        aA = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -0.60890205], [0.0, 0.0, 0.60890205], [0.0, -0.57892643, 0.0], [0.0, 0.57892643, 0.0], [-0.59587172, 0.0, 0.0], [0.59587172, 0.0, 0.0]])
        np.testing.assert_almost_equal(aA, grid.getPoints(), 7, "Transformed equal", True)
        grid.clearConformalTransform()
        aA = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -0.707106781186548], [0.0, 0.0, 0.707106781186548], [0.0, -0.707106781186548, 0.0], [0.0, 0.707106781186548, 0.0], [-0.707106781186548, 0.0, 0.0], [0.707106781186548, 0.0, 0.0]])
        np.testing.assert_almost_equal(aA, grid.getPoints(), 14, "Original equal", True)

        # number of points
        grid = Tasmanian.makeLocalPolynomialGrid(2, 1, 2, 2, "localp")
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN, 13, "num needed")
        self.assertEqual(iNL, 0, "num loaded")
        self.assertEqual(iNP, iNN, "num points")
        ttc.loadExpN2(grid)
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN, 0, "num needed")
        self.assertEqual(iNL, 13, "num loaded")
        self.assertEqual(iNP, iNL, "num points")
        grid.setSurplusRefinement(0.10, 0, "classic")
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN, 8, "num needed")
        self.assertEqual(iNL, 13, "num loaded")
        self.assertEqual(iNP, iNL, "num points")
        grid.setSurplusRefinement(0.05, 0, "classic")
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN, 16, "num needed")
        self.assertEqual(iNL, 13, "num loaded")
        self.assertEqual(iNP, iNL, "num points")
        ttc.loadExpN2(grid)
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN,  0, "num needed")
        self.assertEqual(iNL, 29, "num loaded")
        self.assertEqual(iNP, iNL, "num points")

    def checkUpdate(self):
        '''
        Check the update and make working together.
        '''
        grid = Tasmanian.SparseGrid()
        for iI in range(2):
            if (iI == 0):
                sMake = "grid.makeGlobalGrid(2, 1, 2, 'level', 'leja', [2, 1])"
                sUpda1 = "grid.updateGlobalGrid(2, 'level', [2, 1])"
                sUpda2 = "grid.updateGlobalGrid(3, 'level', [2, 1])"
                sUpda3 = "grid.updateGlobalGrid(4, 'ipcurved', [20, 10, 0, 7])"
            else:
                sMake = "grid.makeSequenceGrid(2, 1, 2, 'level', 'leja', [2, 1])"
                sUpda1 = "grid.updateSequenceGrid(2, 'level', [2, 1])"
                sUpda2 = "grid.updateSequenceGrid(3, 'level', [2, 1])"
                sUpda3 = "grid.updateSequenceGrid(4, 'ipcurved', [20, 10, 0, 7])"
            exec(sMake)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 4, "num needed")
            self.assertEqual(iNL, 0, "num loaded")
            self.assertEqual(iNP, iNN, "num points")
            ttc.loadExpN2(grid)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 0, "num needed")
            self.assertEqual(iNL, 4, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            exec(sUpda1)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 0, "num needed")
            self.assertEqual(iNL, 4, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            exec(sUpda2)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 2, "num needed")
            self.assertEqual(iNL, 4, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0, 0.0]])
            aP = grid.getPoints()
            np.testing.assert_equal(aA, aP, "Update Global not equal", True)
            aP = grid.getLoadedPoints()
            np.testing.assert_equal(aA, aP, "Update Global not equal", True)
            aA = np.array([[0.0, np.sqrt(1.0/3.0)], [1.0, 1.0]])
            aP = grid.getNeededPoints()
            np.testing.assert_equal(aA, aP, "Update Global not equal", True)
            ttc.loadExpN2(grid)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 0, "num needed")
            self.assertEqual(iNL, 6, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            exec(sUpda3)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 1, "num needed")
            self.assertEqual(iNL, 6, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            aA = np.array([[-1.0, 0.0]])
            aP = grid.getNeededPoints()
            np.testing.assert_equal(aA, aP, "Update Global not equal", True)

        grid.makeFourierGrid(2, 1, 2, "level")
        ttc.loadExpN2(grid)
        grid.updateFourierGrid(3, "level")
        self.assertEqual(grid.getNumNeeded(), 60, "failed at updateFourierGrid()")

    def checkDefaults(self):
        '''
        Test that default values are set correctly for alpha/beta/order,
        empty returns, empty copy.
        '''
        grid = Tasmanian.SparseGrid()
        for sRule in TasmanianSG.lsTsgLocalRules:
            grid.makeLocalPolynomialGrid(3, 1, 3, 0, sRule)
            self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
            self.assertEqual(grid.getBeta(), 0.0, "failed beta")
            self.assertEqual(grid.getOrder(), 0, "failed order")
            grid.makeLocalPolynomialGrid(3, 1, 3, 1, sRule)
            self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
            self.assertEqual(grid.getBeta(), 0.0, "failed beta")
            self.assertEqual(grid.getOrder(), 1, "failed order")
            grid.makeLocalPolynomialGrid(3, 1, 3, 2, sRule)
            self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
            self.assertEqual(grid.getBeta(), 0.0, "failed beta")
            self.assertEqual(grid.getOrder(), 2, "failed order")

        grid = Tasmanian.makeWaveletGrid(3, 1, 2, 1)
        self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), 1, "failed order")
        grid.makeWaveletGrid(3, 1, 2, 3)
        self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), 3, "failed order")

        try:
            grid.copyGrid([])
            self.assertTrue(False, "failed to raise exception on copy grid")
        except Tasmanian.InputError as TsgError:
            TsgError.bShowOnExit = False

        # default alpha/beta and order
        grid.makeGlobalGrid(2, 0, 2, 'level', 'leja', [2, 1])
        self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), -1, "failed order")
        grid.makeGlobalGrid(1, 0, 4, 'level', 'gauss-hermite', [], 2.0)
        self.assertEqual(grid.getAlpha(), 2.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), -1, "failed order")
        grid.makeGlobalGrid(1, 0, 4, 'level', 'gauss-jacobi', [], 3.0, 2.0)
        self.assertEqual(grid.getAlpha(), 3.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 2.0, "failed beta")
        self.assertEqual(grid.getOrder(), -1, "failed order")

        grid.makeSequenceGrid(2, 1, 1, 'level', 'leja', [2, 1])
        self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), -1, "failed order")

        # test empty returns
        dummy_ans = np.empty([0, 0], np.float64)
        grid_dummy = Tasmanian.SparseGrid()
        aX = np.empty([0,], np.float64)
        np.testing.assert_equal(dummy_ans, grid_dummy.getPoints(), "Empty read", True)
        np.testing.assert_equal(dummy_ans, grid_dummy.getNeededPoints(), "Empty read", True)
        np.testing.assert_equal(dummy_ans, grid_dummy.getLoadedPoints(), "Empty read", True)
        dummy_ans = np.empty([0,], np.float64)
        np.testing.assert_equal(dummy_ans, grid_dummy.getQuadratureWeights(), "Empty read", True)
        np.testing.assert_equal(dummy_ans, grid_dummy.getInterpolationWeights(aX), "Empty read", True)
        dummy_ans = np.empty([0, 0], np.float64)
        aX = np.empty([0, 0], np.float64)
        np.testing.assert_equal(dummy_ans, grid_dummy.getInterpolationWeightsBatch(aX), "Empty read", True)
        aX = np.empty([1, 0], np.float64)
        np.testing.assert_equal(dummy_ans, grid_dummy.getInterpolationWeightsBatch(aX), "Empty read", True)
        grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev')
        ttc.loadExpN2(grid)
        dummy_ans = np.empty([0, 1], np.float64)
        aX = np.empty([0, 2], np.float64)
        np.testing.assert_equal(dummy_ans, grid.evaluateBatch(aX), "Empty batch eval", True)

    def checkLevelLimits(self):
        '''
        Check setting level limits.
        '''
        grid = Tasmanian.SparseGrid()

        # Level Limits
        grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis', liLevelLimits = [1, 4])
        aPoints = grid.getPoints()
        ttc.checkPoints(aPoints[:,0], lMustHave = [-1.0, 1.0], lMustNotHave = [1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0)])
        ttc.checkPoints(aPoints[:,1], lMustHave = [-1.0, 1.0, 1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0)], lMustNotHave = [])
        ttc.loadExpN2(grid)
        grid.setAnisotropicRefinement('iptotal', 20, 0, [1, 4])
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        ttc.checkPoints(aPoints[:,0], lMustHave = [], lMustNotHave = [1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0)])

        grid.makeSequenceGrid(3, 1, 3, 'level', 'leja', liLevelLimits = [3, 2, 1])
        aPoints = grid.getPoints()
        ttc.checkPoints(aPoints[:,0], lMustHave = [0.0, -1.0, 1.0, 1.0/np.sqrt(3.0)], lMustNotHave = [])
        ttc.checkPoints(aPoints[:,1], lMustHave = [0.0, -1.0, 1.0], lMustNotHave = [1.0/np.sqrt(3.0)])
        ttc.checkPoints(aPoints[:,2], lMustHave = [0.0, 1.0], lMustNotHave = [-1.0, 1.0/np.sqrt(3.0)])
        ttc.loadExpN2(grid)
        grid.setAnisotropicRefinement('iptotal', 5, 0)
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        ttc.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [1.0/np.sqrt(3.0)])
        ttc.checkPoints(aPoints[:,2], lMustHave = [], lMustNotHave = [-1.0, 1.0/np.sqrt(3.0)])
        grid.clearRefinement()
        grid.setAnisotropicRefinement('iptotal', 10, 0, [3, 2, 2])
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        ttc.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [1.0/np.sqrt(3.0)])
        ttc.checkPoints(aPoints[:,2], lMustHave = [1.0], lMustNotHave = [1.0/np.sqrt(3.0)])
        grid.clearRefinement()
        grid.setSurplusRefinement(1.E-8, 0, '', [3, 2, 1])
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        ttc.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [1.0/np.sqrt(3.0)])
        ttc.checkPoints(aPoints[:,2], lMustHave = [], lMustNotHave = [-1.0, 1.0/np.sqrt(3.0)])

        # check that nodes from level 2 (+-0.5) and 3 (+-0.25, +-0.75) appear only in the proper dimension
        grid.makeLocalPolynomialGrid(3, 1, 3, 1, 'localp', [1, 2, 3])
        aPoints = grid.getPoints()
        ttc.checkPoints(aPoints[:,0], lMustHave = [0.0, -1.0, 1.0], lMustNotHave = [0.5, -0.5, -0.75, -0.25, 0.25, 0.75])
        ttc.checkPoints(aPoints[:,1], lMustHave = [0.0, -1.0, 1.0, 0.5, -0.5], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])
        ttc.checkPoints(aPoints[:,2], lMustHave = [0.0, -1.0, 1.0, 0.5, -0.5, -0.75, -0.25, 0.25, 0.75], lMustNotHave = [])
        ttc.loadExpN2(grid)
        grid.setSurplusRefinement(1.E-8, 0, 'classic')
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        ttc.checkPoints(aPoints[:,0], lMustHave = [], lMustNotHave = [0.5, -0.5, -0.75, -0.25, 0.25, 0.75])
        ttc.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])
        ttc.checkPoints(aPoints[:,2], lMustHave = [], lMustNotHave = [])
        grid.clearRefinement()
        grid.setSurplusRefinement(1.E-8, 0, 'classic', [2, 2, 3])
        aPoints = grid.getNeededPoints()
        ttc.checkPoints(aPoints[:,0], lMustHave = [0.5, -0.5], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])
        ttc.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])

        grid.makeWaveletGrid(2, 1, 3, 1, [0, 2])
        aPoints = grid.getPoints()
        ttc.checkPoints(aPoints[:,0], lMustHave = [0.0, -1.0, 1.0], lMustNotHave = [-0.5, 0.5])
        ttc.checkPoints(aPoints[:,1], lMustHave = [0.0, -1.0, 1.0, -0.5, 0.5, -0.75, -0.25, 0.25, 0.75], lMustNotHave = [-0.125, 0.125])
        ttc.loadExpN2(grid)
        grid.setSurplusRefinement(1.E-8, 0, 'classic') # grid is effectively full tensor, cannot refine within the level limits
        self.assertTrue((grid.getNumNeeded() == 0), 'added points even though it should not have');
        grid.setSurplusRefinement(1.E-8, 0, 'classic', [1, 2])
        aPoints = grid.getNeededPoints()
        ttc.checkPoints(aPoints[:,0], lMustHave = [0.5, -0.5], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])
        ttc.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [-0.125, 0.125])

        # level limits I/O
        grid.makeLocalPolynomialGrid(3, 1, 3, 1, 'localp', [1, 2, 3])
        aLimits = grid.getLevelLimits()
        np.testing.assert_almost_equal(aLimits, np.array([1, 2, 3]), 14, "Could not read level limits", True)

        grid.clearLevelLimits()
        aLimits = grid.getLevelLimits()
        np.testing.assert_almost_equal(aLimits, np.array([-1, -1, -1]), 14, "Could not read level limits", True)

    def performMakeUpdateTests(self):
        self.checkPathsVersions()
        self.checkMakeAgainstKnown()
        self.checkUpdate()
        self.checkDefaults()
        self.checkLevelLimits()
