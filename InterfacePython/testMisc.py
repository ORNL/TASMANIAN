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

    def performMiscTests(self):
        self.checkPolynomialSpace()
        self.checkPlotting()
