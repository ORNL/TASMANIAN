import unittest
import TasmanianSG
import sys, os
import numpy as np

from random import shuffle

import testConfigureData as tdata # needed for Gauss-Patterson table file
import testCommon
ttc = testCommon.TestTasCommon()

class TestTasClass(unittest.TestCase):
    '''
    Test the read/write functionality in different formats, also the
    reading of various meta-data (i.e., version or gpu info)
    '''
    def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")

    def testNothing(self):
        pass

    def checkMetaIO(self):
        '''
        Test the version I/O, available acceleration, gpu info, print stats
        '''
        grid = TasmanianSG.TasmanianSparseGrid()
        print("Tasmanian Sparse Grids version: {0:1s}".format(grid.getVersion()))
        print("                 version major: {0:1d}".format(grid.getVersionMajor()))
        print("                 version minor: {0:1d}".format(grid.getVersionMinor()))
        print("                       License: {0:1s}".format(grid.getLicense()))
        if (grid.isOpenMPEnabled()):
            sStatus = "Enabled"
        else:
            sStatus = "Disabled"
        print("                        OpenMP: {0:1s}".format(sStatus))
        sAvailableAcceleration = ""
        for s in ["cpu-blas", "gpu-cublas", "gpu-cuda", "gpu-magma"]:
            if (grid.isAccelerationAvailable(s)):
                sAvailableAcceleration += " " + s
        print("        Available acceleration:{0:1s}".format(sAvailableAcceleration))

        print("                Available GPUs:")
        if (grid.getNumGPUs() > 0):
            for iGPU in range(grid.getNumGPUs()):
                sName = ""
                if ("@CMAKE_CXX_COMPILER_ID@" == "MSVC"):
                    sName = "Unavailable under Windows"
                else:
                    sName = grid.getGPUName(iGPU)
                sMem = grid.getGPUMemory(iGPU)
                print("     {0:2d}: {1:20s} with{2:6d}MB RAM".format(iGPU, sName, sMem))
        else:
            print("            none")

        # covers the printStats() in python and C++
        grid.printStats() # empty grid
        grid.makeGlobalGrid(2, 0, 1, "level", "gauss-gegenbauer", [], 3.0)
        grid.printStats()
        grid.makeGlobalGrid(2, 0, 1, "level", "gauss-jacobi", [], 3.0, 3.0)
        grid.printStats()
        grid.makeGlobalGrid(2, 1, 1, "level", "custom-tabulated", [], 0.0, 0.0, tdata.sGaussPattersonTableFile)
        grid.printStats()
        grid.makeSequenceGrid(2, 1, 3, "level", "rleja")
        ttc.loadExpN2(grid)
        grid.printStats()
        grid.makeLocalPolynomialGrid(1, 1, 3)
        grid.setDomainTransform(np.array([[2.0, 3.0]]))
        grid.printStats()
        grid.makeWaveletGrid(2, 1, 1)
        grid.printStats()
        grid.makeFourierGrid(3, 1, 1, "level")
        grid.enableAcceleration("gpu-cuda")
        grid.printStats()

    def checkReadWriteGlobal(self):
        '''
        Test reading and writing of Global grids.
        '''
        # iDimension, iOutputs, iDepth, sType, sRule, fAlpha, fBeta, useTransform, loadFunciton, limitLevels
        lGrids = [[3, 2, 2, "level", "leja", 0.0, 0.0, False, False, False],
                  [2, 1, 4, "level", "clenshaw-curtis", 0.0, 0.0, False, False, False],
                  [3, 1, 3, "level", "rleja", 0.0, 0.0, True, False, False],
                  [3, 1, 2, "iptotal", "chebyshev", 0.0, 0.0, False, True, False],
                  [3, 1, 3, "level", "leja", 0.0, 0.0, True, True, False],
                  [3, 1, 3, "level", "leja", 0.0, 0.0, True, True, True],
                  [2, 1, 5, "qptotal", "gauss-hermite", 1.0, 3.0, False, False, False],
                  [3, 1, 2, "level", "gauss-laguerre", 3.0, 0.0, False, False, True],]
        aTransform = np.array([[0.0,1.0],[0.0,1.0],[-2.0,-1.0]])

        for lT in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            if (lT[7]):
                gridA.makeGlobalGrid(lT[0], lT[1], lT[2], lT[3], lT[4], [], lT[5], lT[6], [], [3, 2, 1])
            else:
                gridA.makeGlobalGrid(lT[0], lT[1], lT[2], lT[3], lT[4], [], lT[5], lT[6])
            #gridA.printStats()
            if (lT[7]):
                gridA.setDomainTransform(aTransform)
            if (lT[8]):
                ttc.loadExpN2(gridA)

            gridA.write("testSave")
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridB.makeGlobalGrid(1, 0, 1, "level", "rleja")
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.copyGrid(gridA)
            ttc.compareGrids(gridA, gridB)

        # test an error message from wrong read
        print("Attempting a bogus read to see if error would be properly registered")
        self.assertFalse(gridB.read("testSaveBlah"), "Failed to flag a fake read")
        print("GOOD: error was registered")

        # custom rule test
        gridA = TasmanianSG.TasmanianSparseGrid()
        gridB = TasmanianSG.TasmanianSparseGrid()
        gridA.makeGlobalGrid(2, 0, 4, 'level', 'custom-tabulated', [], 0.0, 0.0, tdata.sGaussPattersonTableFile)
        gridB = TasmanianSG.TasmanianSparseGrid()
        gridB.makeGlobalGrid(2, 0, 4, 'level', 'gauss-patterson')
        ttc.compareGrids(gridA, gridB, bTestRuleNames = False)
        gridA.write("testSave")
        gridB.makeGlobalGrid(2, 0, 4, 'level', 'clenshaw-curtis')
        gridB.read("testSave")
        ttc.compareGrids(gridA, gridB)
        gridA.write("testSave", bUseBinaryFormat = True)
        gridB.makeGlobalGrid(3, 0, 4, 'level', 'leja')
        gridB.read("testSave")
        ttc.compareGrids(gridA, gridB)

    def checkReadWriteSequence(self):
        '''
        Test reading and writing of Sequence grids.
        '''
        # iDimension, iOutputs, iDepth, sType, sRule, useTransform, loadFunciton, limitLevels
        lGrids = [[3, 2, 2, "level", "leja", False, False, False],
                  [2, 1, 4, "level", "max-lebesgue", False, False, False],
                  [3, 1, 3, "level", "rleja", True, False, False],
                  [3, 1, 3, "level", "rleja", True, False, True],
                  [3, 1, 2, "iptotal", "min-delta", False, True, False],
                  [3, 1, 3, "level", "leja", True, True, False],
                  [3, 1, 3, "level", "leja", True, True, True],]

        for lT in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            if (lT[7]):
                gridA.makeSequenceGrid(lT[0], lT[1], lT[2], lT[3], lT[4], [], [2, 3, 1])
            else:
                gridA.makeSequenceGrid(lT[0], lT[1], lT[2], lT[3], lT[4])
            if (lT[5]):
                gridA.setDomainTransform(np.array([[0.0,1.0],[0.0,1.0],[-2.0,-1.0]]))
            if (lT[6]):
                ttc.loadExpN2(gridA)

            gridA.write("testSave")
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridB.makeGlobalGrid(1, 0, 1, "level", "rleja")
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.copyGrid(gridA)
            ttc.compareGrids(gridA, gridB)

    def checkReadWriteLocalp(self):
        '''
        Test reading and writing of Localp grids.
        '''
        # iDimension, iOutputs, iDepth, iorder, sRule, useTransform, loadFunciton, limitLevels
        lGrids = [[3, 2, 2, 0, "localp", False, False, False],
                  [3, 0, 2, 0, "localp-boundary", False, False, False],
                  [2, 1, 4, 1, "semi-localp", False, False, False],
                  [3, 1, 3, 2, "localp", True, False, False],
                  [3, 1, 2, 3, "localp-zero", False, True, False],
                  [3, 1, 2, 3, "localp-zero", False, True, True],
                  [3, 1, 3, 4, "semi-localp", True, True, False],
                  [3, 1, 3, -1, "semi-localp", True, True, False],
                  [3, 1, 3, -1, "semi-localp", True, True, True],]

        for lT in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            if (lT[7]):
                gridA.makeLocalPolynomialGrid(lT[0], lT[1], lT[2], lT[3], lT[4], [3, 1, 2])
            else:
                gridA.makeLocalPolynomialGrid(lT[0], lT[1], lT[2], lT[3], lT[4])
            if (lT[5]):
                gridA.setDomainTransform(np.array([[0.0,1.0],[0.0,1.0],[-2.0,-1.0]]))
            if (lT[6]):
                ttc.loadExpN2(gridA)

            gridA.write("testSave")
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridB.makeGlobalGrid(1, 0, 1, "level", "rleja")
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.copyGrid(gridA)
            ttc.compareGrids(gridA, gridB)

    def checkReadWriteWavelet(self):
        '''
        Test reading and writing of Wavelet grids.
        '''
        # iDimension, iOutputs, iDepth, iOrder, useTransform, loadFunciton
        lGrids = [[3, 2, 2, 1, False, False, False],
                  [3, 0, 2, 1, False, False, False],
                  [2, 1, 4, 1, False, False, False],
                  [3, 1, 1, 3, True, False, False],
                  [3, 1, 1, 3, True, False, True],
                  [3, 1, 2, 1, False, True, False],
                  [3, 1, 2, 3, True, True, True],
                  [3, 1, 2, 3, True, True, False],]

        for lT in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            if (lT[6]):
                gridA.makeWaveletGrid(lT[0], lT[1], lT[2], lT[3], [1, 1, 2])
            else:
                gridA.makeWaveletGrid(lT[0], lT[1], lT[2], lT[3])
            if (lT[4]):
                gridA.setDomainTransform(np.array([[0.0,1.0],[0.0,1.0],[-2.0,-1.0]]))
            if (lT[5]):
                ttc.loadExpN2(gridA)

            gridA.write("testSave")
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridB.makeGlobalGrid(1, 0, 1, "level", "rleja")
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.copyGrid(gridA)
            ttc.compareGrids(gridA, gridB)

    def checkReadWriteFourier(self):
        '''
        Test reading and writing of Fourier grids.
        '''
        # iDimension, iOutputs, iDepth, useTransform, loadFunction, useLevelLimits
        lGrids = [[3, 2, 2, False, False, False],
                  [2, 1, 4, False, False, False],
                  [3, 1, 1, True, False, False],
                  [3, 1, 1, True, False, True],
                  [3, 1, 2, False, True, False],
                  [3, 1, 2, True, True, True],
                  [3, 1, 2, True, True, False],]

        for lT in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            if (lT[5]):
                gridA.makeFourierGrid(lT[0], lT[1], lT[2], 'level', [1, 1, 2])
            else:
                gridA.makeFourierGrid(lT[0], lT[1], lT[2], 'level')
            if (lT[3]):
                gridA.setDomainTransform(np.array([[0.0,1.0],[0.0,1.0],[-2.0,-1.0]]))
            if (lT[4]):
                ttc.loadExpN2(gridA)

            gridA.write("testSave")
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridB.makeGlobalGrid(1, 0, 1, "level", "rleja")
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.copyGrid(gridA)
            ttc.compareGrids(gridA, gridB)

    def checkReadWriteMisc(self):
        '''
        Test reading and writing of domain transforms and testing all rules.
        '''
        lGrids = ['gridA.makeGlobalGrid(3, 2, 4, "level", "clenshaw-curtis"); gridA.setDomainTransform(aTransform); gridA.setConformalTransformASIN(np.array([3,4,5]))',
                  'gridA.makeGlobalGrid(3, 2, 4, "level", "gauss-legendre"); gridA.setConformalTransformASIN(np.array([3,5,1]))',
                  'gridA.makeSequenceGrid(2, 2, 5, "level", "leja"); gridA.setConformalTransformASIN(np.array([0,4]))',
                  'gridA.makeLocalPolynomialGrid(3, 1, 4, 2, "localp"); gridA.setDomainTransform(aTransform); gridA.setConformalTransformASIN(np.array([5,3,0]))',
                  'gridA.getNumPoints()']
        aTransform = np.array([[0.0,1.0],[0.0,1.0],[-2.0,-1.0]])

        for sGrid in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            exec(sGrid)
            gridA.write("testSave")
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.read("testSave")
            ttc.compareGrids(gridA, gridB)

            gridB.makeSequenceGrid(1, 1, 0, "level", "leja");
            gridB.makeLocalPolynomialGrid(1, 1, 0)
            gridB.copyGrid(gridA)
            ttc.compareGrids(gridA, gridB)

        # Make a grid with every possible rule (catches false-positive and memory crashes)
        for sType in TasmanianSG.lsTsgGlobalTypes:
            for sRule in TasmanianSG.lsTsgGlobalRules:
                if ("custom-tabulated" in sRule):
                    gridA.makeGlobalGrid(2, 0, 2, sType, sRule, sCustomFilename = tdata.sGaussPattersonTableFile)
                else:
                    gridA.makeGlobalGrid(2, 0, 2, sType, sRule)
                gridA.write("testSave", bUseBinaryFormat = False)
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)
                gridB.makeGlobalGrid(1, 0, 0, "level", "clenshaw-curtis")
                gridA.write("testSave", bUseBinaryFormat = True)
                gridB.read("testSave")

        for sType in TasmanianSG.lsTsgGlobalTypes:
            for sRule in TasmanianSG.lsTsgSequenceRules:
                gridA.makeSequenceGrid(2, 1, 3, sType, sRule)
                gridA.write("testSave")
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)
                gridB.makeGlobalGrid(1, 0, 0, "level", "clenshaw-curtis")
                gridA.write("testSave", bUseBinaryFormat = True)
                gridB.read("testSave")

    def checkReadWriteConstruction(self):
        '''
        Test read/write when using non- construction.
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

    def performIOTest(self):
        self.checkMetaIO()

        self.checkReadWriteGlobal()
        self.checkReadWriteSequence()
        self.checkReadWriteLocalp()
        self.checkReadWriteWavelet()
        self.checkReadWriteFourier()

        self.checkReadWriteConstruction()

        self.checkReadWriteMisc()
