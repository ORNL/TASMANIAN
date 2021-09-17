import unittest
import TasmanianSG
import sys, os
import numpy as np

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
        sGPUbackend = "none"
        if (grid.isCudaEnabled()): sGPUbackend = "CUDA"
        if (grid.isHipEnabled()): sGPUbackend = "ROCm/HIP"
        if (grid.isDpcppEnabled()): sGPUbackend = "oneAPI/DPC++"
        print("                   GPU backend: {0:1s}".format(sGPUbackend))
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

            gridA.write("testSave", bUseBinaryFormat = False)
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
        try:
            gridB.read("Test_If_Bogus_Filename_Produces_an_Error")
        except TasmanianSG.TasmanianInputError as TSGError:
            TSGError.bShowOnExit = False
            self.assertEqual(TSGError.sVariable, "sFilename", "Reading a bogus file properly failed, but the error information is wrong.")

        # custom rule test
        gridA = TasmanianSG.TasmanianSparseGrid()
        gridB = TasmanianSG.TasmanianSparseGrid()
        gridA.makeGlobalGrid(2, 0, 4, 'level', 'custom-tabulated', [], 0.0, 0.0, tdata.sGaussPattersonTableFile)
        gridB = TasmanianSG.TasmanianSparseGrid()
        gridB.makeGlobalGrid(2, 0, 4, 'level', 'gauss-patterson')
        ttc.compareGrids(gridA, gridB, bTestRuleNames = False)
        gridA.write("testSave", bUseBinaryFormat = False)
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

            gridA.write("testSave", bUseBinaryFormat = False)
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

            gridA.write("testSave", bUseBinaryFormat = False)
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

            gridA.write("testSave", bUseBinaryFormat = False)
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

            gridA.write("testSave", bUseBinaryFormat = False)
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

    def checkCopySubgrid(self):
        grid = TasmanianSG.TasmanianSparseGrid() # test grid
        gridTotal = TasmanianSG.TasmanianSparseGrid() # test grid
        gridRef1 = TasmanianSG.TasmanianSparseGrid() # test grid
        gridRef2 = TasmanianSG.TasmanianSparseGrid() # test grid
        gridRef3 = TasmanianSG.TasmanianSparseGrid() # test grid
        gridRef4 = TasmanianSG.TasmanianSparseGrid() # test grid
        # outputs:   source       0         0, 1, 2     2, 3, 4       5
        lGrids = ['gridTotal', 'gridRef1', 'gridRef2', 'gridRef3', 'gridRef4']
        lOutputs = [6, 1, 3, 3, 1];
        lMake = ['.makeGlobalGrid(2, iOutputs, 4, "level", "clenshaw-curtis")',
                 '.makeSequenceGrid(2, iOutputs, 4, "level", "rleja")',
                 '.makeLocalPolynomialGrid(2, iOutputs, 4, 2)',
                 '.makeWaveletGrid(2, iOutputs, 3)',
                 '.makeFourierGrid(2, iOutputs, 4, "level")']

        for sMake in lMake:
            for iI in range(len(lOutputs)):
                iOutputs = lOutputs[iI]
                exec(lGrids[iI] + sMake)

            gridTotalValues = []
            gridRef1Values  = []
            gridRef2Values  = []
            gridRef3Values  = []
            gridRef4Values  = []

            aPoints = gridTotal.getPoints()
            for iI in range(aPoints.shape[0]):
                lModel = [i * np.exp(aPoints[iI][0] + aPoints[iI][1]) for i in range(1, 7)]
                gridTotalValues.append(lModel)
                gridRef1Values.append(lModel[0:1])
                gridRef2Values.append(lModel[0:3])
                gridRef3Values.append(lModel[2:5])
                gridRef4Values.append(lModel[5:6])

            for sGrid in lGrids:
                exec(sGrid + ".loadNeededPoints(np.row_stack(" + sGrid + "Values))")

            grid.copyGrid(gridTotal)
            ttc.compareGrids(grid, gridTotal)
            grid.copyGrid(gridTotal, 0, 1)
            ttc.compareGrids(grid, gridRef1)
            grid.copyGrid(gridTotal, 0, 3)
            ttc.compareGrids(grid, gridRef2)
            grid.copyGrid(gridTotal, 2, 5)
            ttc.compareGrids(grid, gridRef3)
            grid.copyGrid(gridTotal, 5, 6)
            ttc.compareGrids(grid, gridRef4)

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
            gridA.write("testSave", bUseBinaryFormat = False)
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
                gridA.write("testSave", bUseBinaryFormat = False)
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)
                gridB.makeGlobalGrid(1, 0, 0, "level", "clenshaw-curtis")
                gridA.write("testSave", bUseBinaryFormat = True)
                gridB.read("testSave")
                ttc.compareGrids(gridA, gridB)

    def checkReadWriteCustomTabulated(self):
        '''
        Test reading and writing of the CustomTabulated class.
        '''
        description = "testCT"
        def create_nodes(j):
            return np.linspace(-1.0, 1.0, num=j)
        def create_weights(j):
            weights = np.linspace(0.0, 1.0, num=j)
            weights = 2.0 * weights / weights.sum()
            return weights
        # Read and write from explicitly given data.
        for i in range(4):
            num_levels = i
            num_nodes = np.array([3*(j+1) for j in range(i)])
            precision = np.array([2*(j+1)-1 for j in range(i)])
            nodes = [create_nodes(j) for j in num_nodes]
            weights = [create_weights(j) for j in num_nodes]
            ctA = TasmanianSG.makeCustomTabulatedFromData(num_levels, num_nodes, precision, nodes, weights, description)            # .
            ctB = TasmanianSG.CustomTabulated()
            ctA.write("testSave")
            ctB.read("testSave")
            for i in range(num_levels):
                read_weights, read_nodes = ctB.getWeightsNodes(i)
                np.testing.assert_almost_equal(read_weights, weights[i])
                np.testing.assert_almost_equal(read_nodes, nodes[i])
        # Read and write from a file.
        ctA = TasmanianSG.makeCustomTabulatedFromFile(tdata.sGaussPattersonTableFile)
        ctB = TasmanianSG.CustomTabulated()
        ctA.write("testSave")
        ctB.read("testSave")
        ttc.compareCustomTabulated(ctA, ctB)
        for i in range(ctB.getNumLevels()):
            read_weights, read_nodes = ctB.getWeightsNodes(i)
            grid = TasmanianSG.makeGlobalGrid(1, 0, i, "level", "gauss-patterson")
            np.testing.assert_almost_equal(read_weights, grid.getQuadratureWeights())
            np.testing.assert_almost_equal(read_nodes, grid.getPoints().flatten())
        # Test an error message from wrong read.
        try:
            ctB.read("Test_If_Bogus_Filename_Produces_an_Error")
        except TasmanianSG.TasmanianInputError as TSGError:
            TSGError.bShowOnExit = False
            self.assertEqual(TSGError.sVariable, "sFilename", "Reading a bogus file properly failed, but the error information is wrong.")

    def checkGlobalGridCustom(self):
        '''
        Test makeGlobalGridCustom(), which creates a grid from a CustomTabulated instance.
        '''
        gridA = TasmanianSG.makeGlobalGrid(1, 1, 3, "level", "custom-tabulated", [], 0.0, 0.0, tdata.sGaussPattersonTableFile)
        ct = TasmanianSG.makeCustomTabulatedFromFile(tdata.sGaussPattersonTableFile)
        gridB = TasmanianSG.makeGlobalGridCustom(1, 1, 3, "level", ct)
        ttc.compareGrids(gridA, gridB)
        gridA = TasmanianSG.makeGlobalGrid(2, 1, 3, "level", "custom-tabulated", [], 0.0, 0.0, tdata.sGaussPattersonTableFile)
        gridB = TasmanianSG.makeGlobalGridCustom(2, 1, 3, "level", ct)
        ttc.compareGrids(gridA, gridB)

    def performIOTest(self):
        self.checkMetaIO()

        self.checkReadWriteGlobal()
        self.checkReadWriteSequence()
        self.checkReadWriteLocalp()
        self.checkReadWriteWavelet()
        self.checkReadWriteFourier()

        self.checkCopySubgrid()
        self.checkReadWriteMisc()

        self.checkReadWriteCustomTabulated()
        self.checkGlobalGridCustom()
