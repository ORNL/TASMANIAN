import unittest
import TasmanianSG
import numpy as np

from random import shuffle

import testCommon

ttc = testCommon.TestTasCommon()

class TestTasClass(unittest.TestCase):
    '''
    Test the refinement capabilities:
    * set different refinements
    * estimate anisotropic coefficients
    * read/write refinements
    '''
    def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")

    def testNothing(self):
        pass

    def checkSetClear(self):
        '''
        Set refinement and clear refinement
        '''
        grid = TasmanianSG.TasmanianSparseGrid()
        grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'semi-localp')
        ttc.loadExpN2(grid)
        self.assertEqual(grid.getNumNeeded(), 0, "num needed")
        grid.setSurplusRefinement(0.0001, 0, 'classic')
        self.assertTrue((grid.getNumNeeded() > 0), "num needed")
        grid.clearRefinement()
        self.assertEqual(grid.getNumNeeded(), 0, "num needed")

    def checkAnisoCoeff(self):
        '''
        Check anisotropic coefficients
        '''
        grid = TasmanianSG.TasmanianSparseGrid()
        grid.makeGlobalGrid(2, 1, 9, 'level', 'rleja')
        aP = grid.getPoints()
        aV = np.exp(aP[:,0] + aP[:,1]**2)
        grid.loadNeededPoints(aV.reshape([aP.shape[0], 1]))
        aC = grid.estimateAnisotropicCoefficients('iptotal', 0)
        self.assertEqual(len(aC.shape), 1, 'dimensions of the estimated anisotropic weight')
        self.assertEqual(aC.shape[0], 2, 'dimensions of the estimated anisotropic weight')
        self.assertLess(np.abs(float(aC[0]) / float(aC[1]) - 2.0), 0.2, 'wrong anisotropic weights estimated')
        aC = grid.estimateAnisotropicCoefficients('ipcurved', 0)
        self.assertEqual(len(aC.shape), 1, 'dimensions of the estimated anisotropic weight')
        self.assertEqual(aC.shape[0], 4, 'dimensions of the estimated anisotropic weight')
        self.assertLess(np.abs(float(aC[0]) / float(aC[1]) - 2.0), 0.2, 'wrong anisotropic weights estimated, alpha curved')
        self.assertLess(aC[2], 0.0, 'wrong anisotropic weights estimated, beta 1')
        self.assertLess(aC[3], 0.0, 'wrong anisotropic weights estimated, beta 2')

    def checkLocalpSurplus(self):
        '''
        Check surplus refinement for local polynomial grids
        '''
        grid = TasmanianSG.TasmanianSparseGrid()
        grid.makeLocalPolynomialGrid(2, 1, 4, 1, 'semi-localp')
        ttc.loadExpN2(grid)
        aPoints = grid.getPoints()
        aScale = np.array([[1.0 if aPoints[i,0] > 0.0 else 0.0] for i in range(aPoints.shape[0])])
        grid.setSurplusRefinement(1.E-9, 0, 'classic', [], aScale)
        aNeeded = grid.getNeededPoints()
        for iI in range(aNeeded.shape[0]):
            self.assertLess(0.0, aNeeded[iI, 0], 'wrong set of needed points after rescaling')


    def checkFileIO(self):
        '''
        Read/Write regular refinement.
        '''
        grid = TasmanianSG.TasmanianSparseGrid()

        sRefinementIOGrids = ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis')",
                              "grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja')"]
        for sTest in sRefinementIOGrids:
            exec(sTest)
            ttc.loadExpN2(grid)
            grid.setAnisotropicRefinement('iptotal', 20, 0)
            self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
            gridB = TasmanianSG.TasmanianSparseGrid()
            gridB.copyGrid(grid)
            ttc.compareGrids(grid, gridB)

            grid.write("refTestFlename.grid", bUseBinaryFormat = True)
            gridB.makeLocalPolynomialGrid(1, 1, 0, 1)
            gridB.read("refTestFlename.grid")
            ttc.compareGrids(grid, gridB)

            grid.write("refTestFlename.grid", bUseBinaryFormat = False)
            gridB.makeLocalPolynomialGrid(1, 1, 0, 1)
            gridB.read("refTestFlename.grid")
            ttc.compareGrids(grid, gridB)

    def checkConstruction(self):
        '''
        Test read/write when using construction.
        '''
        llTest = ["gridA.makeGlobalGrid(3, 2, 2, 'level', 'clenshaw-curtis'); gridB.makeGlobalGrid(3, 2, 2, 'level', 'clenshaw-curtis')",
                  "gridA.makeSequenceGrid(3, 2, 4, 'level', 'leja'); gridB.makeSequenceGrid(3, 2, 4, 'level', 'leja')",
                  "gridA.makeLocalPolynomialGrid(3, 2, 2); gridB.makeLocalPolynomialGrid(3, 2, 2)",
                  "gridA.makeWaveletGrid(3, 2, 2); gridB.makeWaveletGrid(3, 2, 2)",
                  "gridA.makeFourierGrid(3, 2, 2, 'level'); gridB.makeFourierGrid(3, 2, 2, 'level')",]

        for sMakeGrids in llTest:
            for sFormat in [False, True]: # test binary and ascii format
                gridA = TasmanianSG.TasmanianSparseGrid()
                gridB = TasmanianSG.TasmanianSparseGrid()
                gridC = TasmanianSG.TasmanianSparseGrid()

                exec(sMakeGrids)

                gridA.beginConstruction()
                gridB.beginConstruction()
                #gridA.printStats()

                gridB.write("testSave", bUseBinaryFormat = sFormat)
                gridB.makeSequenceGrid(1, 1, 0, "level", "rleja") # clean the grid
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)
                gridC.copyGrid(gridA)
                ttc.compareGrids(gridA, gridC)

                for t in range(5): # use 5 iterations
                    if (gridA.isLocalPolynomial() or gridA.isWavelet()):
                        aPointsA = gridA.getCandidateConstructionPointsSurplus(1.E-4, "fds")
                        aPointsB = gridB.getCandidateConstructionPointsSurplus(1.E-4, "fds")
                        aPointsC = gridC.getCandidateConstructionPointsSurplus(1.E-4, "fds")
                    else:
                        aPointsA = gridA.getCandidateConstructionPoints("level", 0)
                        aPointsB = gridB.getCandidateConstructionPoints("level", 0)
                        aPointsC = gridC.getCandidateConstructionPoints("level", 0)
                    np.testing.assert_almost_equal(aPointsA, aPointsB, decimal=11)
                    np.testing.assert_almost_equal(aPointsA, aPointsC, decimal=11)

                    iNumPoints = int(aPointsA.shape[0] / 2)
                    if (iNumPoints > 32): iNumPoints = 32

                    # use the first samples (up to 32) and shuffle the order
                    # add one of the samples further in the list
                    liSamples = list(range(iNumPoints + 1))
                    shuffle(liSamples)
                    for iI in range(len(liSamples)):
                        if (liSamples[iI] == iNumPoints):
                            liSamples[iI] = iNumPoints + 1
                    #liSamples = map(lambda i: i if i < iNumPoints else iNumPoints + 1, liSamples)

                    for iI in liSamples: # compute and load the samples
                        aPoint = aPointsA[iI, :]
                        aValue = np.array([np.exp(aPoint[0] + aPoint[1]), 1.0 / ((aPoint[0] - 1.3) * (aPoint[1] - 1.6) * (aPoint[2] - 2.0))])

                        gridA.loadConstructedPoint(aPoint, aValue)
                        gridB.loadConstructedPoint(aPoint, aValue)
                        gridC.loadConstructedPoint(aPoint, aValue)

                    # using straight construction or read/write should produce the same result
                    ttc.compareGrids(gridA, gridC)
                    gridB.write("testSave", bUseBinaryFormat = sFormat)
                    gridB.makeSequenceGrid(1, 1, 0, "level", "rleja")
                    gridB.read("testSave")
                    ttc.compareGrids(gridA, gridB)
                    gridC.copyGrid(gridA)
                    ttc.compareGrids(gridA, gridC)

                gridA.finishConstruction()
                gridB.finishConstruction()

                gridB.write("testSave", bUseBinaryFormat = sFormat)
                gridB.makeSequenceGrid(1, 1, 0, "level", "rleja")
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)
                gridC.copyGrid(gridA)
                ttc.compareGrids(gridA, gridC)

        # check multi-point load
        gridA = TasmanianSG.TasmanianSparseGrid()
        gridA.makeLocalPolynomialGrid(3, 2, 4);
        ttc.loadExpN2(gridA)

        gridB = TasmanianSG.TasmanianSparseGrid()
        gridB.makeLocalPolynomialGrid(3, 2, 0)

        gridB.beginConstruction()
        aX = gridA.getPoints()
        aY = gridA.evaluateBatch(aX)
        gridB.loadConstructedPoint(aX, aY)
        gridB.finishConstruction()
        ttc.compareGrids(gridA, gridB)

        # check some mem-leaks and crashes (correctness is elsewhere)
        gridA = TasmanianSG.TasmanianSparseGrid()
        gridA.makeLocalPolynomialGrid(2, 5, 0)
        gridA.beginConstruction()
        gridA.loadConstructedPoint(np.empty([0, 2]), np.empty([0, 5])) # empty input, check for crash

        gridA.makeLocalPolynomialGrid(2, 1, 1)
        gridA.loadNeededPoints(np.ones([5, 1]))
        gridA.beginConstruction()
        aPoints = gridA.getCandidateConstructionPointsSurplus(1.E-4, "classic") # should generate empty output
        np.testing.assert_almost_equal(aPoints, np.empty([0, 0]), 14, "failed to generate empty list of construction points", True)

        gridA.makeLocalPolynomialGrid(2, 1, 0)
        gridA.loadNeededPoints(np.ones([1, 1]))
        gridA.beginConstruction()
        aPoints = gridA.getCandidateConstructionPointsSurplus(1.E-4, "classic", 0, [], np.array([[1.E-6]])) # should generate empty output
        np.testing.assert_almost_equal(aPoints, np.empty([0, 0]), 14, "failed to generate empty list of construction points", True)

        gridA.makeGlobalGrid(2, 1, 1, "tensor", "clenshaw-curtis")
        gridA.loadNeededPoints(np.ones([9, 1]))
        gridA.beginConstruction()
        aPoints = gridA.getCandidateConstructionPoints("ipcurved", [5, 5, 2, 2], [1, 1]) # should generate empty output
        np.testing.assert_almost_equal(aPoints, np.empty([0, 0]), 14, "failed to generate empty list of construction points", True)

    def checkRemovePoints(self):
        '''
        tests removePointsByHierarchicalCoefficient()
        '''
        grid = TasmanianSG.makeLocalPolynomialGrid(2, 1, 1)
        aPoints = grid.getNeededPoints()
        grid.loadNeededValues(np.exp(-aPoints[:,0]**2 -0.5*aPoints[:,1]**2).reshape((grid.getNumNeeded(), 1)))

        reduced = TasmanianSG.TasmanianSparseGrid()

        reduced.copyGrid(grid)
        reduced.removePointsByHierarchicalCoefficient(0.6)
        self.assertEqual(reduced.getNumPoints(), 3, "failed to remove points with threshold 0.6")
        np.testing.assert_almost_equal(reduced.getLoadedPoints(), np.array([[0.0, 0.0], [-1.0, 0.0], [1.0, 0.0]]), 14, "failed reduce 1", True)

        reduced.copyGrid(grid)
        reduced.removePointsByHierarchicalCoefficient(0.7)
        self.assertEqual(reduced.getNumPoints(), 1, "failed to remove points with threshold 0.7")
        np.testing.assert_almost_equal(reduced.getLoadedPoints(), np.array([[0.0, 0.0]]), 14, "failed reduce 2", True)

        reduced.copyGrid(grid)
        reduced.removePointsByHierarchicalCoefficient(0.0, iNumKeep = 3)
        self.assertEqual(reduced.getNumPoints(), 3, "failed to remove points down to 3")
        np.testing.assert_almost_equal(reduced.getLoadedPoints(), np.array([[0.0, 0.0], [-1.0, 0.0], [1.0, 0.0]]), 14, "failed reduce 3", True)

        reduced.copyGrid(grid)
        reduced.removePointsByHierarchicalCoefficient(0.0, iNumKeep = 1)
        self.assertEqual(reduced.getNumPoints(), 1, "failed to remove points down to 1")
        np.testing.assert_almost_equal(reduced.getLoadedPoints(), np.array([[0.0, 0.0]]), 14, "failed reduce 4", True)

        reduced.copyGrid(grid)
        reduced.removePointsByHierarchicalCoefficient(0.0, aScaleCorrection = np.array([[1.0], [1.0], [1.0], [0.1], [0.1]]), iNumKeep = 3)
        self.assertEqual(reduced.getNumPoints(), 3, "failed to remove corrected points")
        np.testing.assert_almost_equal(reduced.getLoadedPoints(), np.array([[0.0, 0.0], [0.0, -1.0], [0.0, 1.0]]), 14, "failed reduce 5", True)

    def performRefinementTest(self):
        self.checkSetClear()
        self.checkAnisoCoeff()
        self.checkLocalpSurplus()
        self.checkFileIO()
        self.checkConstruction()
        self.checkRemovePoints()
