import unittest
import Tasmanian
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

    def checkLoadNeeded(self):
        '''
        Check the load needed construction procedure.
        '''
        gridA = Tasmanian.SparseGrid()
        gridB = Tasmanian.SparseGrid()

        for i in range(4):
            if i < 2:
                gridA.makeGlobalGrid(2, 1, 3, 'level', 'fejer2')
                gridB.copyGrid(gridA)
                Tasmanian.loadNeededPoints(lambda x, i : np.ones((1,)) * np.exp(-np.sum(x**2)), gridA, i)
            else:
                gridA.makeLocalPolynomialGrid(2, 1, 2, 2)
                gridB.copyGrid(gridA)
                gridA.loadNeededPoints(np.ones((gridA.getNumPoints(), 1)))
                Tasmanian.reloadLoadedPoints(lambda x, i : np.ones((1,)) * np.exp(-np.sum(x**2)), gridA, i - 2)

            ttc.loadExpN2(gridB)
            ttc.compareGrids(gridA, gridB)

    def performAddonTests(self):
        self.checkLoadNeeded()
