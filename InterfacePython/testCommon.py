import unittest
import TasmanianSG
import sys, os
import numpy as np

class TestTasCommon(unittest.TestCase):
    '''
    The class contains functions and untilities that will be used in many
    places throughout the testing suite. The class needs to run a "dummy"
    test from within the class constructor, otherwise an object cannot
    be initialized from outside of this module.
    '''
    def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")

    def testNothing(self):
        pass

    def compareGrids(self, gridA, gridB, bTestRuleNames = True):
        '''
        Compaes two grids by checking points, weights, and several evaluate.
        The evaluates are done on a canonical [-1, 1] interval.

        The test passes if the grids are mathematically idential.

        bTestRuleNames should be True in most cases, i.e., check if the names
        of the rule in each grid matches. The exception is when comparing the
        GaussPatterson grids build with custom file and sRule = "gauss-patterson"
        '''
        # test basic number of points data
        self.assertEqual(gridA.getNumDimensions(), gridB.getNumDimensions(), "error in getNumDimensions()")
        self.assertEqual(gridA.getNumOutputs(), gridB.getNumOutputs(), "error in getNumOutputs()")
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), "error in getNumPoints()")
        self.assertEqual(gridA.getNumLoaded(), gridB.getNumLoaded(), "error in getNumLoaded()")
        self.assertEqual(gridA.getNumNeeded(), gridB.getNumNeeded(), "error in getNumNeeded()")
        self.assertEqual(gridA.isUsingConstruction(), gridB.isUsingConstruction(), "error in isUsingConstruction()")

        if (gridA.getNumPoints() == 0): # emptry grid, nothing else to check
            return

        # load test points (canonical domain only), make sure to avoid grid points, e.g., 1/2, 1/4, etc.
        if (gridA.getNumDimensions() == 1):
            mX1 = np.array([1.0/3.0])
            mX2 = np.array([-1.0/3.0])
            mX3 = np.array([-1.0/5.0])
            mX4 = np.array([1.0/7.0])
            mX5 = np.array([-1.0/7.0])
        elif (gridA.getNumDimensions() == 2):
            mX1 = np.array([1.0/3.0, 1.0/6.0])
            mX2 = np.array([-1.0/3.0, 1.0/6.0])
            mX3 = np.array([-1.0/5.0, -1.0/7.0])
            mX4 = np.array([1.0/7.0, 1.0/5.0])
            mX5 = np.array([-1.0/7.0, -1.0/13.0])
        elif (gridA.getNumDimensions() == 3):
            mX1 = np.array([1.0/3.0, 1.0/6.0, -1.0/5.0])
            mX2 = np.array([-1.0/3.0, 1.0/6.0, 1.0/6.0])
            mX3 = np.array([-1.0/5.0, -1.0/7.0, 1.0/3.0])
            mX4 = np.array([1.0/7.0, 1.0/5.0, 2.0/3.0])
            mX5 = np.array([-1.0/7.0, -1.0/13.0, -2.0/3.0])

        aBatchPoints = np.row_stack([mX1, mX2, mX3, mX4, mX5])

        pA = gridA.getPoints()
        pB = gridB.getPoints()
        np.testing.assert_equal(pA, pB, "Points not equal", True)

        pA = gridA.getLoadedPoints()
        pB = gridB.getLoadedPoints()
        np.testing.assert_equal(pA, pB, "Loaded not equal", True)

        pA = gridA.getNeededPoints()
        pB = gridB.getNeededPoints()
        np.testing.assert_equal(pA, pB, "Needed not equal", True)

        # test rule data
        self.assertEqual(gridA.getAlpha(), gridB.getAlpha(), "error in getAlpha()")
        self.assertEqual(gridA.getBeta(), gridB.getBeta(), "error in getBeta()")
        self.assertEqual(gridA.getOrder(), gridB.getOrder(), "error in getOrder()")

        if (bTestRuleNames):
            self.assertEqual(gridA.getRule(), gridB.getRule(), "error in getRule()")
            self.assertEqual(gridA.getCustomRuleDescription(), gridB.getCustomRuleDescription(), "error in getCustomRuleDescription()")

        self.assertEqual(gridA.isGlobal(), gridB.isGlobal(), "error in isGlobal()")
        self.assertEqual(gridA.isSequence(), gridB.isSequence(), "error in isSequence()")
        self.assertEqual(gridA.isLocalPolynomial(), gridB.isLocalPolynomial(), "error in isLocalPolynomial()")
        self.assertEqual(gridA.isWavelet(), gridB.isWavelet(), "error in isWavelet()")
        self.assertEqual(gridA.isFourier(), gridB.isFourier(), "error in isFourier()")

        # weights
        pA = gridA.getQuadratureWeights()
        pB = gridB.getQuadratureWeights()
        np.testing.assert_equal(pA, pB, "Quadrature not equal", True)

        pA = gridA.getInterpolationWeights(mX1)
        pB = gridB.getInterpolationWeights(mX1)
        np.testing.assert_equal(pA, pB, "Interpolation test 1 not equal", True)

        pA = gridA.getInterpolationWeights(mX2)
        pB = gridB.getInterpolationWeights(mX2)
        np.testing.assert_equal(pA, pB, "Interpolation test 2 not equal", True)

        pA = gridA.getInterpolationWeights(mX3)
        pB = gridB.getInterpolationWeights(mX3)
        np.testing.assert_equal(pA, pB, "Interpolation test 3 not equal", True)

        # evaluate (values have been loaded)
        if (gridA.getNumLoaded() > 0):
            pA = gridA.evaluate(mX4)
            pB = gridB.evaluate(mX4)
            np.testing.assert_almost_equal(pA, pB, 12, "Interpolation test 4 not equal", True)

            pA = gridA.evaluate(mX5)
            pB = gridB.evaluate(mX5)
            np.testing.assert_almost_equal(pA, pB, 12, "Interpolation test 5 not equal", True)

            pA = gridA.integrate()
            pB = gridB.integrate()
            np.testing.assert_almost_equal(pA, pB, 12, "Integration test not equal", True)

            pA = gridA.evaluateBatch(aBatchPoints)
            pB = gridB.evaluateBatch(aBatchPoints)
            np.testing.assert_almost_equal(pA, pB, 12, "Interpolation test 6 (batch) not equal", True)

            pA = gridA.getHierarchicalCoefficients()
            pB = gridB.getHierarchicalCoefficients()
            np.testing.assert_almost_equal(pA, pB, 12, "getHierarchicalCoefficients() not equal", True)

        # domain transforms
        self.assertEqual(gridA.isSetDomainTransfrom(), gridB.isSetDomainTransfrom(), "error in isSetDomainTransfrom()")

        pA = gridA.getDomainTransform()
        pB = gridB.getDomainTransform()
        np.testing.assert_equal(pA, pB, "Domain test no equal", True)

        self.assertEqual(gridA.isSetConformalTransformASIN(), gridB.isSetConformalTransformASIN(), "error in isSetConformalTransformASIN()")

        pA = gridA.getLevelLimits()
        pB = gridB.getLevelLimits()
        np.testing.assert_equal(pA, pB, "Level limit test no equal", True)

        pA = gridA.getConformalTransformASIN()
        pB = gridB.getConformalTransformASIN()
        np.testing.assert_equal(pA, pB, "Conformal transform ASIN not equal", True)

    def compareCustomTabulated(self, ctA, ctB):
        '''
        Compaes two CustomTabulated instances by checking nodes, weights, and other metadata.

        The test passes if the two instances are mathematically identical.
        '''
        # Test metadata.
        self.assertEqual(ctA.getDescription(), ctB.getDescription(), "error in getDescription()")
        self.assertEqual(ctA.getNumLevels(), ctB.getNumLevels(), "error in getNumLevels()")
        for level in range(ctA.getNumLevels()):
            self.assertEqual(ctA.getNumPoints(level), ctB.getNumPoints(level), "error in getNumPoints() at level " + str(level))
            self.assertEqual(ctA.getIExact(level), ctB.getIExact(level), "error in getIExact() at level " + str(level))
            self.assertEqual(ctA.getQExact(level), ctB.getQExact(level), "error in getQExact() at level " + str(level))
            wA, nA = ctA.getWeightsNodes(level)
            wB, nB = ctB.getWeightsNodes(level)
            np.testing.assert_equal(wA, wB, "Weights from getWeightsNodes() are not equal at level " + str(level), True)
            np.testing.assert_equal(nA, nB, "Nodes from getWeightsNodes() are not equal at level " + str(level), True)

    def loadExpN2(self, grid):
        '''
        If there are needed points, load the points with exp(-sum(x)**2) where
        x is each of the needed points.

        The function is needed to test whether read/write correctly works on the
        loaded values, so some values have to be easily loaded even if they
        are not meaningful as in the convergence/correctness tests.
        '''
        if (grid.getNumNeeded() == 0): return
        mPoints = grid.getNeededPoints()
        iOut = grid.getNumOutputs()
        iDim = grid.getNumDimensions()
        grid.loadNeededPoints(np.column_stack([np.exp(-np.sum(mPoints**2, axis=1)) for i in range(iOut)]))

    def checkPoints(self, aPoints, lMustHave, lMustNotHave):
        '''
        aPoints is 1D array of points
        lMustHave and lMustNotHave are the points to test
        Example: checkPoints(aPoints[:,0], [0.0, 1.0, -1.0], [0.5, -0.5])
        '''
        for x in lMustHave:
            self.assertTrue((np.sum(np.where(np.abs(aPoints - x) < 0.001, 1, 0)) > 0), 'did not properly limit level, did not find {0:1.5e}'.format(x))
        for x in lMustNotHave:
            self.assertTrue((np.sum(np.where(np.abs(aPoints - x) < 0.001, 1, 0)) == 0), 'did not properly limit level, did find {0:1.5e}'.format(x))
