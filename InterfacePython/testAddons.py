import unittest
import Tasmanian
import os
import numpy as np

import testCommon

ttc = testCommon.TestTasCommon()

def batch_exp(x):
    aResult = np.zeros((x.shape[0], 1))
    for i in range(x.shape[0]):
        aResult[i, 0] = np.exp(-np.sum(x[i,:]))
    return aResult

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

    def checkConstruction(self):
        gridA = Tasmanian.SparseGrid()
        gridB = Tasmanian.SparseGrid()

        for i in range(12):
            if i < 2:
                gridA.makeGlobalGrid(2, 1, 3, 'level', 'fejer2')
                gridB.copyGrid(gridA)
                Tasmanian.constructAnisotropicSurrogate(lambda x, i : np.ones((1,1)) * np.exp(-np.sum(x**2)),
                                                        gridA.getNumPoints(), i, 1, gridA, "hyperbolic", 0)
            elif i < 4:
                gridA.makeGlobalGrid(2, 1, 3, 'level', 'min-lebesgue')
                gridB.copyGrid(gridA)
                Tasmanian.constructAnisotropicSurrogate(lambda x, y, i : np.ones((1,1)) * np.exp(-np.sum(x**2)),
                                                        gridA.getNumPoints(), i - 2, 1, gridA, "ipcurved", (6, 5, 1, 1),
                                                        bUseInitialGuess = True, sCheckpointFilename = "cpf_addons")
                os.remove("cpf_addons")
            elif i < 6:
                gridA.makeGlobalGrid(2, 1, 3, 'level', 'min-delta', liLevelLimits = (5, 5))
                gridB.copyGrid(gridA)
                gridA.clearLevelLimits()
                Tasmanian.constructAnisotropicSurrogate(lambda x, i : np.ones((1,1)) * np.exp(-np.sum(x**2)),
                                                        gridA.getNumPoints(), i - 4, 1, gridA, "iplevel", (1, 2),
                                                        liLevelLimits = (5, 5), bUseInitialGuess = False)
            elif i < 8:
                gridA.makeGlobalGrid(2, 1, 3, 'level', 'min-delta', liLevelLimits = (10, 10))
                gridB.copyGrid(gridA)
                gridA.clearLevelLimits()
                Tasmanian.constructAnisotropicSurrogate(lambda x, y, i : np.ones((1,1)) * np.exp(-np.sum(x**2)),
                                                        gridA.getNumPoints(), i - 6, 1, gridA, "iplevel", 0,
                                                        liLevelLimits = (10, 10), bUseInitialGuess = True)
            elif i < 10:
                gridA.makeLocalPolynomialGrid(2, 1, 2, 2, liLevelLimits = (7, 7))
                gridB.copyGrid(gridA)
                gridA.clearLevelLimits()
                Tasmanian.constructSurplusSurrogate(lambda x, y, i : np.ones((1,1)) * np.exp(-np.sum(x**2)),
                                                    gridA.getNumPoints(), i - 8, 1, gridA, 1.E-5, "stable", -1,
                                                    liLevelLimits = (7, 7), bUseInitialGuess = True)
            elif i < 12:
                gridA.makeLocalPolynomialGrid(2, 1, 2, 3)
                gridB.copyGrid(gridA)
                Tasmanian.constructSurplusSurrogate(lambda x, i : np.ones((1,1)) * np.exp(-np.sum(x**2)),
                                                    gridA.getNumPoints(), i - 10, 1, gridA, 1.E-5, "fds", 0,
                                                    bUseInitialGuess = False, sCheckpointFilename = "cpf_addons")
                os.remove("cpf_addons")

            ttc.loadExpN2(gridB)
            gridA.finishConstruction()
            ttc.compareGrids(gridA, gridB)

    def checkBatchConstruct(self):
        gridA = Tasmanian.SparseGrid()
        gridB = Tasmanian.SparseGrid()

        gridA.makeGlobalGrid(2, 1, 6, 'level', 'fejer2')
        gridB.copyGrid(gridA)
        Tasmanian.constructAnisotropicSurrogate(lambda x, i : batch_exp(x),
                                                gridA.getNumPoints(), 2, 1, gridA, "iptotal", 0)
        Tasmanian.constructAnisotropicSurrogate(lambda x, i : batch_exp(x),
                                                gridB.getNumPoints(), 2, 10, gridB, "iptotal", 0)
        ttc.compareGrids(gridA, gridB)

        gridA.makeLocalPolynomialGrid(2, 1, 3)
        gridB.copyGrid(gridA)
        Tasmanian.constructSurplusSurrogate(lambda x, i : batch_exp(x),
                                            gridA.getNumPoints(), 2, 1, gridA, 1.E-11, "classic")
        Tasmanian.constructSurplusSurrogate(lambda x, i : batch_exp(x),
                                            gridB.getNumPoints(), 2, 10, gridB, 1.E-11, "classic")
        ttc.compareGrids(gridA, gridB)

    def checkUnstructuredL2(self):
        grid = Tasmanian.makeLocalPolynomialGrid(2, 1, 5)
        if not grid.isAccelerationAvailable("cpu-blas"): return

        x = np.linspace(-0.9, 0.9, 40)
        xx, yy = np.meshgrid(x, x)

        points = np.column_stack([xx.reshape((xx.size,)), yy.reshape((yy.size,))])
        model = np.exp(points[:,0] - points[:,1]).reshape((points.shape[0], 1))

        Tasmanian.loadUnstructuredDataL2(points, model, 1.E-5, grid)

        x = np.linspace(-0.8, 0.8, 20)
        xx, yy = np.meshgrid(x, x)
        points = np.column_stack([xx.reshape((xx.size,)), yy.reshape((yy.size,))])
        model = np.exp(points[:,0] - points[:,1])

        surrogate = grid.evaluateBatch(points)

        self.assertLess(np.max(np.abs(model - surrogate[:,0])), 5.E-3, "Constructed from unstructured too inaccurate.")

    def performAddonTests(self):
        self.checkLoadNeeded()
        self.checkConstruction()
        self.checkBatchConstruct()
        self.checkUnstructuredL2()
