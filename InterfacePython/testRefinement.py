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

            grid.write("refTestFlename.grid")
            gridB.makeLocalPolynomialGrid(1, 1, 0, 1)
            gridB.read("refTestFlename.grid")
            ttc.compareGrids(grid, gridB)

    def checkConstruction(self):
        '''
        Test read/write when using construction.
        '''
        for sFormat in [False, True]: # test binary and ascii format
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            gridA.makeGlobalGrid(3, 2, 2, "level", "clenshaw-curtis")
            gridB.makeGlobalGrid(3, 2, 2, "level", "clenshaw-curtis")

            gridA.beginConstruction()
            gridB.beginConstruction()

            gridB.write("testSave", bUseBinaryFormat = sFormat)
            gridB.makeSequenceGrid(1, 1, 0, "level", "rleja") # clean the grid
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            for t in range(5): # use 5 iterations
                aPointsA = gridA.getCandidateConstructionPoints("level", 0)
                aPointsB = gridB.getCandidateConstructionPoints("level", 0)
                np.testing.assert_almost_equal(aPointsA, aPointsB, decimal=11)

                iNumPoints = int(aPointsA.shape[0] / 2)
                if (iNumPoints > 32): iNumPoints = 32

                # use the first samples (up to 32) and shuffle the order
                # add one of the samples further in the list
                liSamples = list(range(iNumPoints + 1))
                shuffle(liSamples)
                liSamples = map(lambda i: i if i < iNumPoints else iNumPoints + 1, liSamples)

                for iI in liSamples: # compute and load the samples
                    aPoint = aPointsA[iI, :]
                    aValue = np.array([np.exp(aPoint[0] + aPoint[1]), 1.0 / ((aPoint[0] - 1.3) * (aPoint[1] - 1.6) * (aPoint[2] - 2.0))])

                    gridA.loadConstructedPoint(aPoint, aValue)
                    gridB.loadConstructedPoint(aPoint, aValue)

                # using straight construction or read/write should produce the same result
                gridB.write("testSave", bUseBinaryFormat = sFormat)
                gridB.makeSequenceGrid(1, 1, 0, "level", "rleja")
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)

            gridA.finishConstruction()
            gridB.finishConstruction()

            gridB.write("testSave", bUseBinaryFormat = sFormat)
            gridB.makeSequenceGrid(1, 1, 0, "level", "rleja")
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

    def performRefinementTest(self):
        self.checkSetClear()
        self.checkAnisoCoeff()
        self.checkFileIO()
        self.checkConstruction()
