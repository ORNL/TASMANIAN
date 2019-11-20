import unittest
import TasmanianSG
import numpy as np

import testCommon # needed to compare grids after merge refinement

ttc = testCommon.TestTasCommon()

class TestTasClass(unittest.TestCase):
    '''
    Test methods for directly setting hierarchical matrices and coeffs,
    useful to build grids from unstructured/random data.
    '''
    def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")

    def testNothing(self):
        pass

    def checkAgainstKnown(self):
        '''
        Check hierarchical functions against known (pen-and-paper)
        solutions.
        '''
        grid = TasmanianSG.TasmanianSparseGrid()

        # evalHierarchicalBasis (all done in batch), dense mode
        grid.makeLocalPolynomialGrid(2, 2, 1, 1, 'localp')
        aT = np.empty([0,2], np.float64)
        np.testing.assert_equal(aT, grid.getHierarchicalCoefficients(), "coeff mismatch", True)
        grid.makeLocalPolynomialGrid(2, 0, 1, 1, 'localp')
        aT = np.empty([0,0], np.float64)

        np.testing.assert_equal(aT, grid.getHierarchicalCoefficients(), "coeff mismatch", True)
        grid.makeLocalPolynomialGrid(2, 2, 1, 1, 'localp')
        aPoints = grid.getPoints()
        aV = np.column_stack([np.exp(2.0 * aPoints[:,0] + aPoints[:,1]), np.sin(3.0 * aPoints[:,0] + aPoints[:,1])])
        grid.loadNeededPoints(aV)
        aS = grid.getHierarchicalCoefficients()
        aT = np.array([[1.0, 0.0], [np.exp(-1.0)-1.0, np.sin(-1.0)], [np.exp(1.0)-1.0, np.sin(1.0)], [np.exp(-2.0)-1.0, np.sin(-3.0)], [np.exp(2.0)-1.0, np.sin(3.0)]])
        np.testing.assert_almost_equal(aS, aT, 14, "Surplusses equal", True)

        grid.makeGlobalGrid(3, 1, 4, 'level', 'fejer2')
        aPoints = grid.getPoints()
        aVan = grid.evaluateHierarchicalFunctions(aPoints)
        aIdentity = np.zeros([grid.getNumPoints(), grid.getNumPoints()], np.float64)
        for iI in range(grid.getNumPoints()):
            aIdentity[iI,iI] = 1.0
        np.testing.assert_almost_equal(aVan, aIdentity, 14, "evaluateHierarchicalFunctions", True)

        grid.makeGlobalGrid(2, 1, 1, 'level', 'clenshaw-curtis')
        aPoints = np.array([[0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666]])
        f0 = (1.0 - aPoints[:,0]**2) + (1.0 - aPoints[:,1]**2) - 1.0
        f1 = 0.5 * aPoints[:,1] * (aPoints[:,1] - 1.0)
        f2 = 0.5 * aPoints[:,1] * (aPoints[:,1] + 1.0)
        f3 = 0.5 * aPoints[:,0] * (aPoints[:,0] - 1.0)
        f4 = 0.5 * aPoints[:,0] * (aPoints[:,0] + 1.0)
        aResult = np.column_stack([f0, f1, f2, f3, f4])
        aVan = grid.evaluateHierarchicalFunctions(aPoints)
        np.testing.assert_almost_equal(aVan, aResult, 14, "evaluateHierarchicalFunctions", True)

        grid.makeSequenceGrid(2, 1, 2, 'level', 'leja')
        aPoints = np.array([[0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666]])
        f0 = np.ones([aPoints.shape[0],])
        f1 = aPoints[:,1]
        f2 = 0.5 * aPoints[:,1] * (aPoints[:,1] - 1.0)
        f3 = aPoints[:,0]
        f4 = aPoints[:,0] * aPoints[:,1]
        f5 = 0.5 * aPoints[:,0] * (aPoints[:,0] - 1.0)
        aResult = np.column_stack([f0, f1, f2, f3, f4, f5])
        aVan = grid.evaluateHierarchicalFunctions(aPoints)
        np.testing.assert_almost_equal(aVan, aResult, 14, "evaluateHierarchicalFunctions", True)

        grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp')
        aPoints = np.array([[0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666]])

        f0 = np.ones([aPoints.shape[0],])
        f1 = -aPoints[:,1]
        f2 = np.copy(aPoints[:,1])
        f3 = -aPoints[:,0]
        f4 = np.copy(aPoints[:,0])
        for iI in range(aPoints.shape[0]):
            if (aPoints[iI,1] > 0.0): f1[iI] = 0.0
            if (aPoints[iI,1] < 0.0): f2[iI] = 0.0
            if (aPoints[iI,0] > 0.0): f3[iI] = 0.0
            if (aPoints[iI,0] < 0.0): f4[iI] = 0.0
        aResult = np.column_stack([f0, f1, f2, f3, f4])
        aVan = grid.evaluateHierarchicalFunctions(aPoints)
        np.testing.assert_almost_equal(aVan, aResult, 14, "evaluateHierarchicalFunctions", True)

        grid.makeFourierGrid(2, 1, 1, "level")
        aPoints = np.array([[0.25, 0.5]]);
        aResult = np.array([[1.0 + 0.0 * 1j, -1.0 + 0.0 * 1j, -1.0 + 0.0 * 1j, 0.0 - 1.0 * 1j, 0.0 + 1.0 * 1j]])
        aVan = grid.evaluateHierarchicalFunctions(aPoints)
        np.testing.assert_almost_equal(aVan, aResult, 14, "evaluateHierarchicalFunctions", True)

        # sparse hierarchical functions
        pSparse = TasmanianSG.TasmanianSimpleSparseMatrix()
        np.testing.assert_almost_equal(np.empty([0,0], np.float64), pSparse.getDenseForm(), 14, "TasmanianSimpleSparseMatrix.getDense()", True)

        grid.makeLocalPolynomialGrid(2, 1, 4, 1, 'localp')
        aPoints = np.array([[0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666]])
        aDense = grid.evaluateHierarchicalFunctions(aPoints)
        pSparse = grid.evaluateSparseHierarchicalFunctions(aPoints)
        np.testing.assert_almost_equal(aDense, pSparse.getDenseForm(), 14, "evaluateSparseHierarchicalFunctions", True)

        grid.makeWaveletGrid(2, 1, 4, 1)
        aPoints = np.array([[0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666]])
        aDense = grid.evaluateHierarchicalFunctions(aPoints)
        pSparse = grid.evaluateSparseHierarchicalFunctions(aPoints)
        np.testing.assert_almost_equal(aDense, pSparse.getDenseForm(), 14, "evaluateSparseHierarchicalFunctions", True)

    def checkSetCoeffsMergeRefine(self):
        '''
        Set the coefficients and use mergeRefinement()
        '''
        aPoints = np.array([[0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666]])
        # set hierarchical coeffs/merge refinement
        gridA = TasmanianSG.TasmanianSparseGrid()
        gridB = TasmanianSG.TasmanianSparseGrid()
        gridA.makeGlobalGrid(2, 1, 4, 'level', 'chebyshev')
        gridB.makeGlobalGrid(2, 1, 4, 'level', 'chebyshev')
        ttc.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)

        gridA = TasmanianSG.TasmanianSparseGrid()
        gridB = TasmanianSG.TasmanianSparseGrid()
        gridA.makeGlobalGrid(2, 1, 4, 'level', 'fejer2')
        gridB.makeGlobalGrid(2, 1, 4, 'level', 'fejer2')
        ttc.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)
        gridA.setAnisotropicRefinement('iptotal', 10, 0)
        gridB.setAnisotropicRefinement('iptotal', 10, 0)
        ttc.loadExpN2(gridA)
        gridB.mergeRefinement()
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), 'different number of points after merge refinement')
        aRes = gridB.evaluateBatch(aPoints)
        np.testing.assert_almost_equal(aRes, np.zeros(aRes.shape), 14, "not zero after merged refinement", True)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)

        gridA.makeSequenceGrid(2, 1, 5, 'ipcurved', 'min-delta')
        gridB.makeSequenceGrid(2, 1, 5, 'ipcurved', 'min-delta')
        ttc.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)
        gridA.setAnisotropicRefinement('iptotal', 30, 0)
        gridB.setAnisotropicRefinement('iptotal', 30, 0)
        ttc.loadExpN2(gridA)
        gridB.mergeRefinement()
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), 'different number of points after merge refinement')
        aRes = gridB.evaluateBatch(aPoints)
        np.testing.assert_almost_equal(aRes, np.zeros(aRes.shape), 14, "not zero after merged refinement", True)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)

        gridA.makeLocalPolynomialGrid(2, 1, 4, 2, 'semi-localp')
        gridB.makeLocalPolynomialGrid(2, 1, 4, 2, 'semi-localp')
        ttc.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)
        gridA.setSurplusRefinement(1.E-4, 0, 'classic')
        gridB.setSurplusRefinement(1.E-4, 0, 'classic')
        ttc.loadExpN2(gridA)
        gridB.mergeRefinement()
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), 'different number of points after merge refinement')
        aRes = gridB.evaluateBatch(aPoints)
        np.testing.assert_almost_equal(aRes, np.zeros(aRes.shape), 14, "not zero after merged refinement", True)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)

        gridA.makeWaveletGrid(2, 1, 2, 1)
        gridB.makeWaveletGrid(2, 1, 2, 1)
        ttc.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)
        gridA.setSurplusRefinement(1.E-3, 0, 'classic')
        gridB.setSurplusRefinement(1.E-3, 0, 'classic')
        ttc.loadExpN2(gridA)
        gridB.mergeRefinement()
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), 'different number of points after merge refinement')
        aRes = gridB.evaluateBatch(aPoints)
        np.testing.assert_almost_equal(aRes, np.zeros(aRes.shape), 14, "not zero after merged refinement", True)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)

        gridA.makeFourierGrid(2, 1, 3, 'level')
        gridB.makeFourierGrid(2, 1, 3, 'level')
        ttc.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        ttc.compareGrids(gridA, gridB)

    def checkIntegrals(self):
        grid = TasmanianSG.TasmanianSparseGrid()

        lTests = ['grid.makeGlobalGrid(2, 1, 4, "level", "clenshaw-curtis")',
                  'grid.makeSequenceGrid(2, 1, 4, "level", "rleja")',
                  'grid.makeFourierGrid(2, 1, 3, "level")',
                  'grid.makeLocalPolynomialGrid(2, 1, 4)',
                  'grid.makeWaveletGrid(2, 1, 2)']

        for t in lTests:
            exec(t)
            grid.setDomainTransform(np.array([[-3.0, 5.0], [-4.0, 2.0]]))
            ttc.loadExpN2(grid)
            aIntegrals = grid.integrateHierarchicalFunctions()
            aSurps = grid.getHierarchicalCoefficients()
            fResult = np.real(np.sum(aIntegrals * aSurps.reshape((aSurps.shape[0],))))
            fExpected = grid.integrate()
            np.testing.assert_almost_equal(fResult, fExpected[0], 12, 'integrals different by too much')

        # check support
        grid.makeLocalPolynomialGrid(1, 1, 3)
        aSup = grid.getHierarchicalSupport()
        aResult = (1.0, 1.0, 1.0, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25)
        np.testing.assert_almost_equal(aSup, np.array(aResult).reshape((9, 1)), 12, 'support different by too much')

    def checkValues(self):
        lTests = ['grid.makeGlobalGrid(2, 2, 3, "level", "clenshaw-curtis")',
                  'grid.makeSequenceGrid(2, 2, 3, "level", "rleja")',
                  'grid.makeFourierGrid(2, 2, 2, "level")',
                  'grid.makeLocalPolynomialGrid(2, 2, 2)',
                  'grid.makeWaveletGrid(2, 2, 2)']

        for t in lTests:
            grid = TasmanianSG.TasmanianSparseGrid()
            aResult = grid.getLoadedValues()
            self.assertTrue(len(aResult) == 0, "failed to return a empty array")
            exec(t)
            iNP = grid.getNumNeeded()
            iOuts = 2
            aReferece = np.array([float(i) for i in range(iNP * iOuts)]).reshape(iNP, iOuts)
            grid.loadNeededPoints(aReferece)
            aResult = grid.getLoadedValues()
            np.testing.assert_almost_equal(aResult, aReferece, 14, 'values in getLoadedValues() differ by too much')

    def performUnstructuredDataTests(self):
        self.checkAgainstKnown()
        self.checkSetCoeffsMergeRefine()
        self.checkIntegrals()
        self.checkValues()
