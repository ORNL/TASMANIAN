import unittest
import TasmanianSG
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

    def checkPolynomialSpace(self):
        '''
        Construct grids with specific parameters and check the polynomial
        space against known (pen-and-paper) exact spaces.
        '''
        grid = TasmanianSG.TasmanianSparseGrid()

        aA = np.array([[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2], [3, 0], [4, 0]])
        grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis')
        aP = grid.getGlobalPolynomialSpace(True)
        np.testing.assert_equal(aP, aA, "poly space mismatch", True)

        aA = np.array([[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [1, 5], [2, 0], [2, 1], [2, 2], [2, 3], [3, 0], [3, 1], [3, 2], [3, 3], [4, 0], [4, 1], [5, 0], [5, 1]])
        grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis')
        aP = grid.getGlobalPolynomialSpace(False)
        np.testing.assert_equal(aP, aA, "poly space mismatch", True)

        aA = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [2, 0]])
        grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja')
        aP = grid.getGlobalPolynomialSpace(True)
        np.testing.assert_equal(aP, aA, "poly space mismatch", True)

        grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja')
        aP = grid.getGlobalPolynomialSpace(False)
        np.testing.assert_equal(aP, aA, "poly space mismatch", True)

        aA = np.array([[0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [1, 3], [2, 0], [2, 1], [3, 0], [3, 1]])
        grid.makeGlobalGrid(2, 1, 1, 'level', 'leja-odd')
        aP = grid.getGlobalPolynomialSpace(False)
        np.testing.assert_equal(aP, aA, "poly space mismatch", True)

        grid.makeSequenceGrid(2, 1, 2, 'level', 'leja')
        aP = grid.getGlobalPolynomialSpace(False)
        np.testing.assert_equal(aP, aA, "poly space mismatch", True)

    def checkPlotting(self):
        '''
        If matplotlib is available and there is an active display, then
        make plot commands but do not show the final plots.
        If matplotlib is missing, check the exceptions.
        '''
        grid = TasmanianSG.TasmanianSparseGrid()

        if (TasmanianSG.bTsgPlotting):
            if ((os.name == "posix") and (os.environ.get('DISPLAY') is None)):
                print("NOTE: there is no display, cannot run plotting tests.")
                return
            grid.makeGlobalGrid(2, 1, 20, 'level', 'leja')
            ttc.loadExpN2(grid)
            grid.plotResponse2D()
            grid.plotPoints2D()

            grid.makeGlobalGrid(2, 1, 0, 'level', 'leja')
            ttc.loadExpN2(grid)
            grid.plotResponse2D()
            grid.plotPoints2D()

            grid.makeGlobalGrid(3, 1, 3, 'level', 'leja')
            try:
                grid.plotPoints2D()
                self.assertTrue(False, "failed to flag plot exception when when using a grid with other than 2D")
            except TasmanianSG.TasmanianInputError as TSGError:
                TSGError.bShowOnExit = False
                self.assertEqual(TSGError.sVariable, "plotPoints2D", "error raising exception for plotPoints2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotResponse2D()
                self.assertTrue(False, "failed to flag plot exception when when using a grid with other than 2D")
            except TasmanianSG.TasmanianInputError as TSGError:
                TSGError.bShowOnExit = False
                self.assertEqual(TSGError.sVariable, "plotResponse2D", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))

            grid.makeGlobalGrid(2, 2, 2, 'level', 'leja')
            try:
                grid.plotResponse2D(iOutput = -1)
                self.assertTrue(False, "failed to flag plot exception when when using wrong output")
            except TasmanianSG.TasmanianInputError as TSGError:
                TSGError.bShowOnExit = False
                self.assertEqual(TSGError.sVariable, "iOutput", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotResponse2D(iOutput = 3)
                self.assertTrue(False, "failed to flag plot exception when when using wrong output")
            except TasmanianSG.TasmanianInputError as TSGError:
                TSGError.bShowOnExit = False
                self.assertEqual(TSGError.sVariable, "iOutput", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotResponse2D(iNumDim0 = -1)
                self.assertTrue(False, "failed to flag plot exception when when using wrong output")
            except TasmanianSG.TasmanianInputError as TSGError:
                TSGError.bShowOnExit = False
                self.assertEqual(TSGError.sVariable, "iNumDim0", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotResponse2D(iNumDim1 = -1)
                self.assertTrue(False, "failed to flag plot exception when when using wrong output")
            except TasmanianSG.TasmanianInputError as TSGError:
                TSGError.bShowOnExit = False
                self.assertEqual(TSGError.sVariable, "iNumDim1", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
        else:
            grid.makeGlobalGrid(2, 1, 3, 'level', 'leja')
            ttc.loadExpN2(grid)
            try:
                grid.plotResponse2D()
                self.assertTrue(False, "failed to flag plotexception  when missing matplotlib")
            except TasmanianSG.TasmanianInputError as TSGError:
                TSGError.bShowOnExit = False
                self.assertEqual(TSGError.sVariable, "plotResponse2D", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotPoints2D()
                self.assertTrue(False, "failed to flag plot exception when missing matplotlib")
            except TasmanianSG.TasmanianInputError as TSGError:
                TSGError.bShowOnExit = False
                self.assertEqual(TSGError.sVariable, "plotPoints2D", "error raising exception for plotPoints2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))

    def checkDerivatives(self):
        '''
        Check the consistency of derivatives in terms of dimensions and ordering.
        '''
        func1 = lambda x : 2.0 * x[0] * x[0] + x[1] * x[1] / 2.0
        grid1 = TasmanianSG.makeGlobalGrid(2, 1, 4, "iptotal", "gauss-legendre")
        points1 = grid1.getNeededPoints()
        values1 = np.resize(np.apply_along_axis(func1, 1, points1), [grid1.getNumNeeded(), 1])
        grid1.loadNeededValues(values1)
        grad1 = grid1.differentiate(np.array([3.0, 4.0]))
        self.assertEqual(grad1.shape, (2,))
        self.assertTrue(np.allclose(grad1, np.array([12.0, 4.0])))

        func2 = lambda x : [2.0 * x[0] * x[0] + x[1] * x[1] / 2.0 + x[2] * x[2], x[0] * x[1] * x[2]]
        grid2 = TasmanianSG.makeGlobalGrid(3, 2, 4, "iptotal", "gauss-legendre")
        points2 = grid2.getNeededPoints()
        values2 = np.resize(np.apply_along_axis(func2, 1, points2), [grid2.getNumNeeded(), 2])
        grid2.loadNeededValues(values2)
        grad2 = grid2.differentiate(np.array([1.0, 2.0, 3.0]))
        self.assertEqual(grad2.shape, (2, 3))
        self.assertTrue(np.allclose(grad2, np.array([[4.0, 2.0, 6.0], [6.0, 3.0, 2.0]])))

    def performMiscTests(self):
        self.checkPolynomialSpace()
        self.checkPlotting()
        self.checkDerivatives()
