#!@Tasmanian_string_python_hashbang@

import unittest
import TasmanianSG
import math
import sys, os
import numpy as np
from random import uniform
from ctypes import cdll

grid = TasmanianSG.TasmanianSparseGrid()

# python-coverage run testTSG.py
# python-coverage html
# python-coverage report

class TestTasmanian(unittest.TestCase):
    def compareGrids(self, gridA, gridB, bTestRuleNames = True):
        self.assertEqual(gridA.getNumDimensions(), gridB.getNumDimensions(), "error in getNumDimensions()")
        self.assertEqual(gridA.getNumOutputs(), gridB.getNumOutputs(), "error in getNumOutputs()")
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), "error in getNumPoints()")
        self.assertEqual(gridA.getNumLoaded(), gridB.getNumLoaded(), "error in getNumLoaded()")
        self.assertEqual(gridA.getNumNeeded(), gridB.getNumNeeded(), "error in getNumNeeded()")

        mX1 = np.array([1.0/3.0, 1.0/6.0])
        mX2 = np.array([-1.0/3.0, 1.0/6.0])
        mX3 = np.array([-1.0/5.0, -1.0/7.0])
        mX4 = np.array([1.0/7.0, 1.0/5.0])
        mX5 = np.array([-1.0/7.0, -1.0/13.0])

        if (gridA.getNumDimensions() == 3):
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

        self.assertEqual(gridA.getAlpha(), gridB.getAlpha(), "error in getAlpha()")
        self.assertEqual(gridA.getBeta(), gridB.getBeta(), "error in getBeta()")
        self.assertEqual(gridA.getOrder(), gridB.getOrder(), "error in getOrder()")

        if (bTestRuleNames):
            self.assertEqual(gridA.getRule(), gridB.getRule(), "error in getRule()")
            self.assertEqual(gridA.getCustomRuleDescription(), gridB.getCustomRuleDescription(), "error in getCustomRuleDescription()")

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

        if (gridA.getNumLoaded() > 0):
            pA = gridA.evaluate(mX4)
            pB = gridB.evaluate(mX4)
            np.testing.assert_almost_equal(pA, pB, 15, "Interpolation test 4 not equal", True)

            pA = gridA.evaluate(mX5)
            pB = gridB.evaluate(mX5)
            np.testing.assert_almost_equal(pA, pB, 15, "Interpolation test 5 not equal", True)

            pA = gridA.integrate()
            pB = gridB.integrate()
            np.testing.assert_almost_equal(pA, pB, 15, "Integration test not equal", True)

            pA = gridA.evaluateBatch(aBatchPoints)
            pB = gridB.evaluateBatch(aBatchPoints)
            np.testing.assert_almost_equal(pA, pB, 15, "Interpolation test 6 (batch) not equal", True)

            pA = gridA.getHierarchicalCoefficients()
            pB = gridB.getHierarchicalCoefficients()
            np.testing.assert_almost_equal(pA, pB, 15, "getHierarchicalCoefficients() not equal", True)

        self.assertEqual(gridA.isGlobal(), gridB.isGlobal(), "error in isGlobal()")
        self.assertEqual(gridA.isSequence(), gridB.isSequence(), "error in isSequence()")
        self.assertEqual(gridA.isLocalPolynomial(), gridB.isLocalPolynomial(), "error in isLocalPolynomial()")
        self.assertEqual(gridA.isWavelet(), gridB.isWavelet(), "error in isWavelet()")

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

    def loadExpN2(self, grid):
        if (grid.getNumNeeded() == 0): return
        mPoints = grid.getNeededPoints()
        iOut = grid.getNumOutputs()
        iDim = grid.getNumDimensions()
        grid.loadNeededPoints(np.column_stack([np.exp(-np.sum(mPoints**2, axis=1)) for i in range(iOut)]))

    def testBasicIO(self):
        print("\nTesting core I/O test")

        # test library meta I/O
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
        for s in ["cpu-blas", "gpu-cublas", "gpu-cuda", "gpu-default"]:
            if (grid.isAccelerationAvailable(s)):
                sAvailableAcceleration += " " + s
        print("        Available acceleration:{0:1s}".format(sAvailableAcceleration))

        print("                Available GPUs:")
        if (grid.getNumGPUs() > 0):
            for iGPU in range(grid.getNumGPUs()):
                sName = grid.getGPUName(iGPU)
                sMem = grid.getGPUMemory(iGPU)
                print("     {0:2d}: {1:20s} with{2:6d}MB RAM".format(iGPU, sName, sMem))
        else:
            print("            none")

        #grid.printStats() # covers the printStats(), but silence for a release as it prints an empty grid

        # test I/O for Global Grids
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
                self.loadExpN2(gridA)

            gridA.write("testSave")
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridB.makeGlobalGrid(1, 0, 1, "level", "rleja")
            gridB.copyGrid(gridA)
            self.compareGrids(gridA, gridB)

        # test an error message from wrong read
        #print("Attempting a bogus read to see if error would be properly registered")
        gridB.disableLog()
        self.assertFalse(gridB.read("testSaveBlah"), "Failed to flag a fake read")
        gridB.setErrorLogCerr()

        # custom rule test
        grid.makeGlobalGrid(2, 0, 4, 'level', 'custom-tabulated', [], 0.0, 0.0, "GaussPattersonRule.table")
        grid1 = TasmanianSG.TasmanianSparseGrid()
        grid1.makeGlobalGrid(2, 0, 4, 'level', 'gauss-patterson')
        self.compareGrids(grid, grid1, bTestRuleNames = False)
        grid.write("testSave")
        grid1.makeGlobalGrid(2, 0, 4, 'level', 'clenshaw-curtis')
        grid1.read("testSave")
        self.compareGrids(grid, grid1)
        grid.write("testSave", bUseBinaryFormat = True)
        grid1.makeGlobalGrid(3, 0, 4, 'level', 'leja')
        grid1.read("testSave")
        self.compareGrids(grid, grid1)

        # test I/O for Sequence Grids
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
                gridA.setDomainTransform(aTransform)
            if (lT[6]):
                self.loadExpN2(gridA)

            gridA.write("testSave")
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridB.makeGlobalGrid(1, 0, 1, "level", "rleja")
            gridB.copyGrid(gridA)
            self.compareGrids(gridA, gridB)

        # test I/O for Local Polynomial Grids
        # iDimension, iOutputs, iDepth, iorder, sRule, useTransform, loadFunciton, limitLevels
        lGrids = [[3, 2, 2, 0, "localp", False, False, False],
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
                gridA.setDomainTransform(aTransform)
            if (lT[6]):
                self.loadExpN2(gridA)

            gridA.write("testSave")
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridB.makeGlobalGrid(1, 0, 1, "level", "rleja")
            gridB.copyGrid(gridA)
            self.compareGrids(gridA, gridB)

        # test I/O for Local Wavelet Grids
        # iDimension, iOutputs, iDepth, iOrder, useTransform, loadFunciton
        lGrids = [[3, 2, 2, 1, False, False, False],
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
                gridA.setDomainTransform(aTransform)
            if (lT[5]):
                self.loadExpN2(gridA)

            gridA.write("testSave")
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridB.makeGlobalGrid(1, 0, 1, "level", "rleja")
            gridB.copyGrid(gridA)
            self.compareGrids(gridA, gridB)

        lGrids = ['gridA.makeGlobalGrid(3, 2, 4, "level", "clenshaw-curtis"); gridA.setDomainTransform(aTransform); gridA.setConformalTransformASIN(np.array([3,4,5]))',
                  'gridA.makeGlobalGrid(3, 2, 4, "level", "gauss-legendre"); gridA.setConformalTransformASIN(np.array([3,5,1]))',
                  'gridA.makeSequenceGrid(2, 2, 5, "level", "leja"); gridA.setConformalTransformASIN(np.array([0,4]))',
                  'gridA.makeLocalPolynomialGrid(3, 1, 4, 2, "localp"); gridA.setDomainTransform(aTransform); gridA.setConformalTransformASIN(np.array([5,3,0]))',]

        for sGrid in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            exec(sGrid)
            gridA.write("testSave")
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridA.write("testSave", bUseBinaryFormat = True)
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)

            gridB.makeSequenceGrid(1, 1, 0, "level", "leja");
            gridB.copyGrid(gridA)
            self.compareGrids(gridA, gridB)

        # Make a grid with every possible rule (catches false-positive and memory crashes)
        for sType in TasmanianSG.lsTsgGlobalTypes:
            for sRule in TasmanianSG.lsTsgGlobalRules:
                if ("custom-tabulated" in sRule):
                    gridA.makeGlobalGrid(2, 0, 2, sType, sRule, sCustomFilename = "GaussPattersonRule.table")
                else:
                    gridA.makeGlobalGrid(2, 0, 2, sType, sRule)
                gridA.write("testSave", bUseBinaryFormat = False)
                gridB.read("testSave")
                self.compareGrids(gridA, gridB)
                gridB.makeGlobalGrid(1, 0, 0, "level", "clenshaw-curtis")
                gridA.write("testSave", bUseBinaryFormat = True)
                gridB.read("testSave")

        for sType in TasmanianSG.lsTsgGlobalTypes:
            for sRule in TasmanianSG.lsTsgSequenceRules:
                gridA.makeSequenceGrid(2, 1, 3, sType, sRule)
                gridA.write("testSave")
                gridB.read("testSave")
                self.compareGrids(gridA, gridB)
                gridB.makeGlobalGrid(1, 0, 0, "level", "clenshaw-curtis")
                gridA.write("testSave", bUseBinaryFormat = True)
                gridB.read("testSave")

    def testAcceleratedEvaluate(self):
        print("\nTesting accelerated evaluate consistency")

        iGPUID = @Tasmanian_TESTS_GPU_ID@

        # acceleration meta-data
        grid = TasmanianSG.TasmanianSparseGrid()
        grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'semi-localp')
        for accel in TasmanianSG.lsTsgAccelTypes:
            grid.enableAcceleration(accel)
            sA = grid.getAccelerationType()
            bTest = ((accel not in sA) or (sA not in accel))
            self.assertFalse(bTest, "set/get Acceleration")

        if (@Tasmanian_cmake_synctest_enable@):
            bHasBlas = ("@Tasmanian_ENABLE_BLAS@" == "ON")
            bHasCuBlas = ("@Tasmanian_ENABLE_CUBLAS@" == "ON")
            bHasCuda = ("@Tasmanian_ENABLE_CUDA@" == "ON")
            self.assertTrue((grid.isAccelerationAvailable("cpu-blas") == bHasBlas), "failed to match blas")
            self.assertTrue((grid.isAccelerationAvailable("gpu-cublas") == bHasCuBlas), "failed to match cublas")
            self.assertTrue((grid.isAccelerationAvailable("gpu-cuda") == bHasCuda), "failed to match cuda")
            self.assertTrue((grid.isAccelerationAvailable("gpu-default") == (bHasCuBlas or bHasCuda)), "failed to match cuda")

        self.assertTrue((grid.getGPUID() == 0), "did not default to gpu 0")

        if (grid.getNumGPUs() > 1):
            grid.setGPUID(1)
            self.assertTrue((grid.getGPUID() == 1), "did not set to gpu 1")
            sName = grid.getGPUName(1) # mostly checks for memory leaks and crashes

        if (grid.getNumGPUs() > 0):
            grid.setGPUID(0)
            self.assertTrue((grid.getGPUID() == 0), "did not set to gpu 0")
            sName = grid.getGPUName(0) # mostly checks for memory leaks and crashes

        # consistency with evaluate, not a timing test
        grid = TasmanianSG.TasmanianSparseGrid()

        aTestPointsCanonical = np.array([[ uniform(-1.0, 1.0) for j in range(2) ] for i in range(100) ])
        aTestPointsTransformed = np.array([[ uniform(3.0, 5.0) for j in range(2) ] for i in range(100) ])
        aDomainTransform = np.array([[3.0, 5.0],[3.0, 5.0]])

        iFastEvalSubtest = 6

        lTests = [ 'grid.makeGlobalGrid(2, 2, 4, "level", "clenshaw-curtis")',
                   'grid.makeGlobalGrid(2, 2, 4, "level", "chebyshev")',
                   'grid.makeSequenceGrid(2, 2, 4, "level", "leja")',
                   'grid.makeLocalPolynomialGrid(2, 3, 4, 1, "localp")',
                   'grid.makeLocalPolynomialGrid(2, 2, 4, 2, "localp")',
                   'grid.makeLocalPolynomialGrid(2, 1, 4, 3, "localp")',
                   'grid.makeLocalPolynomialGrid(2, 1, 4, 4, "semi-localp")',
                   'grid.makeWaveletGrid(2, 1, 3, 1)',
                   'grid.makeWaveletGrid(2, 1, 3, 3)' ]
        lTests = ['grid.makeLocalPolynomialGrid(2, 3, 4, 1, "localp")']

        iNumGPUs = grid.getNumGPUs()
        lsAccelTypes = ["none", "cpu-blas", "gpu-cuda", "gpu-cublas"]

        for sTest in lTests:
            for iI in range(2):
                iC = 0
                iGPU = 0
                if (iGPUID > -1):
                    iGPU = iGPUID
                while (iC < len(lsAccelTypes)):
                    sAcc = lsAccelTypes[iC]
                    exec(sTest)
                    if (iI == 0):
                        aTestPoints = aTestPointsCanonical
                    else:
                        aTestPoints = aTestPointsTransformed
                        grid.setDomainTransform(aDomainTransform)

                    grid.enableAcceleration(sAcc)
                    if ((sAcc == "gpu-cuda") or (sAcc == "gpu-cublas")):
                        if (iGPU < grid.getNumGPUs()): # without cuda or cublas, NumGPUs is 0 and cannot set GPU
                            grid.setGPUID(iGPU)

                    self.loadExpN2(grid)

                    aRegular = np.array([grid.evaluateThreadSafe(aTestPoints[i,:]) for i in range(aTestPoints.shape[0]) ])
                    aBatched = grid.evaluateBatch(aTestPoints)
                    np.testing.assert_almost_equal(aRegular, aBatched, 14, "Batch evaluation test not equal: {0:1s}, acceleration: {1:1s}, gpu: {2:1d}".format(sTest, sAcc, iGPU), True)

                    aFast = np.array([ grid.evaluate(aTestPoints[i,:]) for i in range(iFastEvalSubtest) ])
                    np.testing.assert_almost_equal(aRegular[0:iFastEvalSubtest,:], aFast, 14, "Batch evaluation test not equal: {0:1s}, acceleration: {1:1s}, gpu: {2:1d}".format(sTest, sAcc, iGPU), True)

                    if ((sAcc == "gpu-cuda") or (sAcc == "gpu-cublas")):
                        if (iGPUID == -1):
                            iGPU += 1
                            if (iGPU >= iNumGPUs):
                                iC += 1
                                iGPU = 0
                        else:
                            iC += 1
                    else:
                        iC += 1

    def testBasicException(self):
        print("\nTesting error handling")
        grid = TasmanianSG.TasmanianSparseGrid()

        # notError tests here are needed to ensure that the multi-statement commands fail for the right function
        llTests = [["grid.makeGlobalGrid(-1, 1,  4, 'level', 'clenshaw-curtis')", "iDimension"],
                   ["grid.makeGlobalGrid(2, -1,  4, 'level', 'clenshaw-curtis')", "iOutputs"],
                   ["grid.makeGlobalGrid(2,  1, -4, 'level', 'clenshaw-curtis')", "iDepth"],
                   ["grid.makeGlobalGrid(2,  1,  4, 'wrong', 'clenshaw-curtis')", "sType"],
                   ["grid.makeGlobalGrid(2,  1,  4, 'level', 'clenshaw-wrong')", "sRule"],
                   ["grid.makeGlobalGrid(2,  1,  4, 'level', 'clenshaw-curtis', [1,2,3])", "liAnisotropicWeights"],
                   ["grid.makeGlobalGrid(2,  1,  4, 'level', 'clenshaw-curtis', [1,2], liLevelLimits = [1, 2, 3])", "liLevelLimits"],
                   ["grid.makeGlobalGrid(2,  1,  4, 'level', 'clenshaw-curtis', [1,2], liLevelLimits = [1, 2])", "notError"],
                   ["grid.makeSequenceGrid(-1, 1,  4, 'level', 'leja')", "iDimension"],
                   ["grid.makeSequenceGrid(2, -1,  4, 'level', 'leja')", "iOutputs"],
                   ["grid.makeSequenceGrid(2,  1, -4, 'level', 'leja')", "iDepth"],
                   ["grid.makeSequenceGrid(2,  1,  4, 'wrong', 'leja')", "sType"],
                   ["grid.makeSequenceGrid(2,  1,  4, 'level', 'weja')", "sRule"],
                   ["grid.makeSequenceGrid(2,  1,  4, 'level', 'leja', [1, 2, 3])", "liAnisotropicWeights"],
                   ["grid.makeSequenceGrid(2,  1,  4, 'level', 'leja', [1, 2], [1, 2, 3])", "liLevelLimits"],
                   ["grid.makeSequenceGrid(2,  1,  4, 'level', 'leja', liLevelLimits = [1, 2])", "notError"],
                   ["grid.makeLocalPolynomialGrid(-1, 1,  4,  2, 'localp')", "iDimension"],
                   ["grid.makeLocalPolynomialGrid(2, -1,  4,  2, 'localp')", "iOutputs"],
                   ["grid.makeLocalPolynomialGrid(2,  1, -4,  2, 'localp')", "iDepth"],
                   ["grid.makeLocalPolynomialGrid(2,  1,  4, -2, 'localp')", "iOrder"],
                   ["grid.makeLocalPolynomialGrid(2,  1,  4,  2, 'lowrong')", "sRule"],
                   ["grid.makeLocalPolynomialGrid(2,  1,  4,  2, 'localp', [1, 2, 3])", "liLevelLimits"],
                   ["grid.makeLocalPolynomialGrid(2,  1,  4,  2, 'localp', [1, 2])", "notError"],
                   ["grid.makeWaveletGrid(-1, 1,  4,  1)", "iDimension"],
                   ["grid.makeWaveletGrid(2, -1,  4,  1)", "iOutputs"],
                   ["grid.makeWaveletGrid(2,  1, -4,  3)", "iDepth"],
                   ["grid.makeWaveletGrid(2,  1,  4,  2)", "iOrder"],
                   ["grid.makeWaveletGrid(2,  1,  4,  1, [1, 2, 3])", "liLevelLimits"],
                   ["grid.makeWaveletGrid(2,  1,  4,  1, [2, 1])", "notError"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja')", "notError"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev')", "notError"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2]))", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateGlobalGrid(1,'iptotal')", "updateGlobalGrid"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev'); grid.updateGlobalGrid(-1,'iptotal')", "iDepth"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev'); grid.updateGlobalGrid(4,'wrong')", "sType"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev'); grid.updateGlobalGrid(4,'iptotal',[1,2,3])", "liAnisotropicWeights"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev'); grid.updateGlobalGrid(4,'iptotal',liLevelLimits = [1,2,3])", "liLevelLimits"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev'); grid.updateGlobalGrid(4,'iptotal',liLevelLimits = [1,2])", "notError"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(4,'iptotal')", "updateSequenceGrid"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(-1,'iptotal')", "iDepth"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(4,'wrong')", "sType"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(4,'iptotal',[1,2,3])", "liAnisotropicWeights"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(4,'iptotal',liLevelLimits = [1,2,3])", "liLevelLimits"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(4,'iptotal',liLevelLimits = [2,3])", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); aW = grid.getInterpolationWeights(np.array([1,2,3]))", "lfX"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); aW = grid.getInterpolationWeightsBatch(np.array([1,2,3]))", "llfX"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); aW = grid.getInterpolationWeightsBatch(np.array([[1,2,3],[1,2,3]]))", "llfX"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([1,]))", "llfVals"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,3]))", "llfVals"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([5,2]))", "llfVals"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2])); grid.loadNeededPoints(np.ones([5,2]))", "llfVals"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.evaluate(np.zeros([1,2]))", "evaluate"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.evaluateThreadSafe(np.zeros(np.array([1,2])))", "evaluateThreadSafe"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2])); grid.evaluate(np.zeros([1,3]))", "lfX"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2])); grid.evaluate(np.zeros([3]))", "lfX"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2])); grid.evaluateThreadSafe(np.zeros([1,3]))", "lfX"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2])); grid.evaluateThreadSafe(np.zeros([3]))", "lfX"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.evaluateBatch(np.zeros([1,2]))", "evaluateBatch"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2])); grid.evaluateBatch(np.zeros([1,3]))", "llfX"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2])); grid.evaluateBatch(np.zeros([2,]))", "llfX"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.integrate()", "integrate"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'gauss-legendre'); grid.setDomainTransform(np.zeros([2,]))", "llfTransform"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'gauss-legendre'); grid.setDomainTransform(np.zeros([2,1]))", "llfTransform"],
                   ["grid.makeGlobalGrid(3, 1, 2, 'level', 'gauss-legendre'); grid.setDomainTransform(np.zeros([2,2]))", "llfTransform"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'gauss-legendre')", "notError"],
                   ["grid.makeGlobalGrid(3, 1, 2, 'level', 'gauss-legendre'); grid.setDomainTransform(np.zeros([3,2]))", "notError"],
                   ["grid.makeGlobalGrid(3, 1, 2, 'level', 'clenshaw-curtis');", "notError"],
                   ["grid.makeGlobalGrid(3, 1, 2, 'level', 'clenshaw-curtis'); grid.setConformalTransformASIN(np.array([0,2,4]))", "notError"],
                   ["grid.makeGlobalGrid(3, 1, 2, 'level', 'clenshaw-curtis'); grid.setConformalTransformASIN(np.array([0,2]))", "liTruncation"],
                   ["grid.makeGlobalGrid(3, 1, 2, 'level', 'clenshaw-curtis'); grid.setConformalTransformASIN(np.array([[0,2,3],[1,2,3]]))", "liTruncation"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'rleja'); grid.setAnisotropicRefinement('iptotal', 10, 1);", "setAnisotropicRefinement"],
                   ["grid.makeGlobalGrid(2, 0, 2, 'level', 'clenshaw-curtis'); grid.setAnisotropicRefinement('iptotal', 10, 1);", "setAnisotropicRefinement"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'fejer2'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', -2, 1);", "iMinGrowth"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -1);", "iOutput"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, 0);", "notError"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, 0, [2, 3]);", "notError"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, 0, [2, 3, 3]);", "liLevelLimits"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -1);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -2);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, 5);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.setAnisotropicRefinement('wrong', 10, 0);", "sType"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -1, [3, 4]);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -1, [3, 4, 5]);", "liLevelLimits"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'rleja'); grid.estimateAnisotropicCoefficients('iptotal', 1);", "estimateAnisotropicCoefficients"],
                   ["grid.makeGlobalGrid(2, 0, 2, 'level', 'clenshaw-curtis'); grid.estimateAnisotropicCoefficients('iptotal', 1);", "estimateAnisotropicCoefficients"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); self.loadExpN2(grid); grid.estimateAnisotropicCoefficients('iptotal', -1);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.estimateAnisotropicCoefficients('iptotal', -1);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.estimateAnisotropicCoefficients('ipcurved', -1);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.estimateAnisotropicCoefficients('iptotal', -2);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.estimateAnisotropicCoefficients('iptotal', 5);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); self.loadExpN2(grid); grid.estimateAnisotropicCoefficients('wrong', 0);", "sType"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); grid.setSurplusRefinement(1.E-4, 0, 'classic');", "setSurplusRefinement"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.setSurplusRefinement(1.E-4, 0, 'classic');", "setSurplusRefinement"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); self.loadExpN2(grid); grid.setSurplusRefinement(-1.E-4, 0, 'classic');", "fTolerance"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); self.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic');", "sCriteria"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); self.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); self.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, '', [2, 3, 4]);", "liLevelLimits"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); self.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, '', [2, 3]);", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); self.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0);", "sCriteria"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); self.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic');", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); self.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic', [2, 3, 4]);", "liLevelLimits"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); self.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic', [2, 3]);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.removePointsByHierarchicalCoefficient(1.E-4, 0);", "removePointsByHierarchicalCoefficient"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.removePointsByHierarchicalCoefficient(-1.E-4, 0);", "fTolerance"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.removePointsByHierarchicalCoefficient(1.E-4, -2);", "iOutput"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.removePointsByHierarchicalCoefficient(1.E-4, 3);", "iOutput"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.removePointsByHierarchicalCoefficient(1.E-4, 0);", "removePointsByHierarchicalCoefficient"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); self.loadExpN2(grid); grid.removePointsByHierarchicalCoefficient(1.E-4, 0);", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([3,]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([10,]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([11,2]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([13,2]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([13,3]));", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, 0,  np.ones([13,3]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, 0,  np.ones([11,]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, 0,  np.ones([13,]));", "notError"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); grid.evaluateHierarchicalFunctions(np.array([[1.0, 1.0], [0.5, 0.3]]));", "notError"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); grid.evaluateHierarchicalFunctions(np.array([[1.0,], [0.5,]]));", "llfX"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); grid.evaluateHierarchicalFunctions(np.array([1.0, 1.0]));", "llfX"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.evaluateSparseHierarchicalFunctions(np.array([[1.0, 1.0], [0.5, 0.3]]));", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.evaluateSparseHierarchicalFunctions(np.array([[1.0,], [0.5,]]));", "llfX"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.evaluateSparseHierarchicalFunctions(np.array([1.0, 1.0]));", "llfX"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.setHierarchicalCoefficients(np.array([1.0, 1.0]));", "llfCoefficients"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.setHierarchicalCoefficients(np.array([[1.0, 1.0],]));", "llfCoefficients"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.setHierarchicalCoefficients(np.array([[1.0, 1.0],[1.0, 1.0],[1.0, 1.0],[1.0, 1.0],[1.0, 1.0],]));", "llfCoefficients"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.setHierarchicalCoefficients(np.array([[1.0,],[1.0,],[1.0,],[1.0,],[1.0,]]));", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.enableAcceleration('gpu-wrong');", "sAccelerationType"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.enableAcceleration('gpu-default');", "notError"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.isAccelerationAvailable('cpu-wrong');", "sAccelerationType"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.isAccelerationAvailable('cpu-blas');", "notError"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.getGPUMemory(-1);", "iGPUID"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.getGPUMemory(1000000);", "iGPUID"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.getGPUMemory(grid1.getNumGPUs());", "iGPUID"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.getGPUName(-1);", "iGPUID"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.getGPUName(1000000);", "iGPUID"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.getGPUName(grid1.getNumGPUs());", "iGPUID"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.setGPUID(-1);", "iGPUID"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.setGPUID(1000000);", "iGPUID"],
                   ["grid1 = TasmanianSG.TasmanianSparseGrid(); grid1.setGPUID(grid1.getNumGPUs());", "iGPUID"],]

        for lTest in llTests:
            try:
                exec(lTest[0])
                self.assertEqual(lTest[1], "notError", "failed to raise exception for invalid '{0:1s}' using test\n '{1:1s}'".format(lTest[1],lTest[0]))
            except TasmanianSG.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, lTest[1], "error raising exception for '{0:1s}' using test\n '{1:1s}'\n Error.sVariable = '{2:1s}'".format(lTest[1],lTest[0],TSGError.sVariable))

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

    def testFullCoverageA(self):
        print("\nTesting core make/update grid")

        if (@Tasmanian_cmake_synctest_enable@):
            grid = TasmanianSG.TasmanianSparseGrid("@Tasmanian_libsparsegrid_path@")
            sVersion = grid.getVersion()
            self.assertEqual(sVersion, TasmanianSG.__version__, "version mismatch")

            pLibTSG = cdll.LoadLibrary("@Tasmanian_libsparsegrid_path@")
            grid = TasmanianSG.TasmanianSparseGrid(pLibTSG)
            sVersion = grid.getVersion()
            self.assertEqual(sVersion, TasmanianSG.__version__, "version mismatch")

        grid = TasmanianSG.TasmanianSparseGrid()

        sVersion = grid.getVersion()
        self.assertEqual(sVersion, TasmanianSG.__version__, "version mismatch")
        sLicense = grid.getLicense()
        self.assertEqual(sLicense, TasmanianSG.__license__, "license mismatch")

        iVM = int(sVersion.split('.')[0])
        iVm = int(sVersion.split('.')[1])
        self.assertEqual(iVM, grid.getVersionMajor(), "version major mismatch")
        self.assertEqual(iVm, grid.getVersionMinor(), "version minor mismatch")

        grid.makeGlobalGrid(1, 0, 4, 'level', 'gauss-hermite', [], 2.0)
        aW = grid.getQuadratureWeights()
        self.assertTrue(np.abs(sum(aW) - 0.5 * np.sqrt(np.pi)) < 1.E-14, "Gauss-Hermite Alpha")

        grid.makeGlobalGrid(2, 0, 2, 'level', 'leja', [2, 1])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0, 0.0]])
        aP = grid.getPoints()
        np.testing.assert_equal(aA, aP, 'Anisotropy Global not equal', True)

        grid.makeGlobalGrid(2, 0, 4, 'ipcurved', 'leja', [20, 10, 0, 7])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [0.0, math.sqrt(1.0/3.0)], [1.0, 0.0], [1.0, 1.0], [-1.0, 0.0]])
        aP = grid.getPoints()
        np.testing.assert_equal(aA, aP, 'Anisotropy Global not equal', True)

        grid.makeSequenceGrid(2, 1, 2, 'level', 'leja', [2, 1])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0, 0.0]])
        aP = grid.getPoints()
        np.testing.assert_equal(aA, aP, 'Anisotropy Sequence not equal', True)

        grid.makeSequenceGrid(2, 1, 4, 'ipcurved', 'leja', [20, 10, 0, 7])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [0.0, math.sqrt(1.0/3.0)], [1.0, 0.0], [1.0, 1.0], [-1.0, 0.0]])
        aP = grid.getPoints()
        np.testing.assert_equal(aA, aP, 'Anisotropy Sequence not equal', True)

        # this is a very important test, checks the curved rule and covers the non-lower-set index selection code
        grid.makeGlobalGrid(2, 1, 1, 'ipcurved', 'rleja', [10, 10, -21, -21])
        aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, 0.5], [0.0, 0.25], [0.0, 0.75], [0.0, 0.125],
                       [1.0, 0.0], [1.0, 1.0], [1.0, 0.5], [1.0, 0.25], [1.0, 0.75], [1.0, 0.125],
                       [0.5, 0.0], [0.5, 1.0], [0.5, 0.5], [0.5, 0.25], [0.5, 0.75], [0.5, 0.125],
                       [0.25, 0.0], [0.25, 1.0], [0.25, 0.5], [0.25, 0.25], [0.25, 0.75],
                       [0.75, 0.0], [0.75, 1.0], [0.75, 0.5], [0.75, 0.25],
                       [0.125, 0.0], [0.125, 1.0], [0.125, 0.5]])
        aA = np.cos(math.pi * aA)
        aP = grid.getPoints()
        np.testing.assert_almost_equal(aA, aP, 12, 'Anisotropy heavily curved not equal', True) # 12 is the number of dec places

        for sRule in TasmanianSG.lsTsgLocalRules:
            grid.makeLocalPolynomialGrid(3, 1, 3, 0, sRule)
            self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
            self.assertEqual(grid.getBeta(), 0.0, "failed beta")
            self.assertEqual(grid.getOrder(), 0, "failed order")
            grid.makeLocalPolynomialGrid(3, 1, 3, 1, sRule)
            self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
            self.assertEqual(grid.getBeta(), 0.0, "failed beta")
            self.assertEqual(grid.getOrder(), 1, "failed order")
            grid.makeLocalPolynomialGrid(3, 1, 3, 2, sRule)
            self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
            self.assertEqual(grid.getBeta(), 0.0, "failed beta")
            self.assertEqual(grid.getOrder(), 2, "failed order")

        grid.makeWaveletGrid(3, 1, 2, 1)
        self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), 1, "failed order")
        grid.makeWaveletGrid(3, 1, 2, 3)
        self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), 3, "failed order")

        try:
            grid.copyGrid([])
            self.assertTrue(False, "failed to raise exception on copy grid")
        except TasmanianSG.TasmanianInputError as TsgError:
            with open(os.devnull, 'w') as devnul:
                sys.stdout = devnul
                TsgError.printInfo() # Miro: silence this for a release
                sys.stdout = sys.__stdout__

        for iI in range(2):
            if (iI == 0):
                sMake = "grid.makeGlobalGrid(2, 1, 2, 'level', 'leja', [2, 1])"
                sUpda1 = "grid.updateGlobalGrid(2, 'level', [2, 1])"
                sUpda2 = "grid.updateGlobalGrid(3, 'level', [2, 1])"
                sUpda3 = "grid.updateGlobalGrid(4, 'ipcurved', [20, 10, 0, 7])"
            else:
                sMake = "grid.makeSequenceGrid(2, 1, 2, 'level', 'leja', [2, 1])"
                sUpda1 = "grid.updateSequenceGrid(2, 'level', [2, 1])"
                sUpda2 = "grid.updateSequenceGrid(3, 'level', [2, 1])"
                sUpda3 = "grid.updateSequenceGrid(4, 'ipcurved', [20, 10, 0, 7])"
            exec(sMake)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 4, "num needed")
            self.assertEqual(iNL, 0, "num loaded")
            self.assertEqual(iNP, iNN, "num points")
            self.loadExpN2(grid)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 0, "num needed")
            self.assertEqual(iNL, 4, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            exec(sUpda1)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 0, "num needed")
            self.assertEqual(iNL, 4, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            exec(sUpda2)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 2, "num needed")
            self.assertEqual(iNL, 4, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            aA = np.array([[0.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0, 0.0]])
            aP = grid.getPoints()
            np.testing.assert_equal(aA, aP, "Update Global not equal", True)
            aP = grid.getLoadedPoints()
            np.testing.assert_equal(aA, aP, "Update Global not equal", True)
            aA = np.array([[0.0, math.sqrt(1.0/3.0)], [1.0, 1.0]])
            aP = grid.getNeededPoints()
            np.testing.assert_equal(aA, aP, "Update Global not equal", True)
            self.loadExpN2(grid)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 0, "num needed")
            self.assertEqual(iNL, 6, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            exec(sUpda3)
            iNN = grid.getNumNeeded()
            iNL = grid.getNumLoaded()
            iNP = grid.getNumPoints()
            self.assertEqual(iNN, 1, "num needed")
            self.assertEqual(iNL, 6, "num loaded")
            self.assertEqual(iNP, iNL, "num points")
            aA = np.array([[-1.0, 0.0]])
            aP = grid.getNeededPoints()
            np.testing.assert_equal(aA, aP, "Update Global not equal", True)

        grid.makeGlobalGrid(2, 0, 2, 'level', 'leja', [2, 1])
        self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), -1, "failed order")
        grid.makeGlobalGrid(1, 0, 4, 'level', 'gauss-hermite', [], 2.0)
        self.assertEqual(grid.getAlpha(), 2.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), -1, "failed order")
        grid.makeGlobalGrid(1, 0, 4, 'level', 'gauss-jacobi', [], 3.0, 2.0)
        self.assertEqual(grid.getAlpha(), 3.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 2.0, "failed beta")
        self.assertEqual(grid.getOrder(), -1, "failed order")

        grid.makeSequenceGrid(2, 1, 1, 'level', 'leja', [2, 1])
        self.assertEqual(grid.getAlpha(), 0.0, "failed alpha")
        self.assertEqual(grid.getBeta(), 0.0, "failed beta")
        self.assertEqual(grid.getOrder(), -1, "failed order")

        grid.makeLocalPolynomialGrid(2, 1, 2, 2, "localp")
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN, 13, "num needed")
        self.assertEqual(iNL, 0, "num loaded")
        self.assertEqual(iNP, iNN, "num points")
        self.loadExpN2(grid)
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN, 0, "num needed")
        self.assertEqual(iNL, 13, "num loaded")
        self.assertEqual(iNP, iNL, "num points")
        grid.setSurplusRefinement(0.10, 0, "classic")
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN, 8, "num needed")
        self.assertEqual(iNL, 13, "num loaded")
        self.assertEqual(iNP, iNL, "num points")
        grid.setSurplusRefinement(0.05, 0, "classic")
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN, 16, "num needed")
        self.assertEqual(iNL, 13, "num loaded")
        self.assertEqual(iNP, iNL, "num points")
        self.loadExpN2(grid)
        iNN = grid.getNumNeeded()
        iNL = grid.getNumLoaded()
        iNP = grid.getNumPoints()
        self.assertEqual(iNN,  0, "num needed")
        self.assertEqual(iNL, 29, "num loaded")
        self.assertEqual(iNP, iNL, "num points")

        # test empty returns
        dummy_ans = np.empty([0, 0], np.float64)
        grid_dummy = TasmanianSG.TasmanianSparseGrid()
        aX = np.empty([0,], np.float64)
        np.testing.assert_equal(dummy_ans, grid_dummy.getPoints(), "Empty read", True)
        np.testing.assert_equal(dummy_ans, grid_dummy.getNeededPoints(), "Empty read", True)
        np.testing.assert_equal(dummy_ans, grid_dummy.getLoadedPoints(), "Empty read", True)
        dummy_ans = np.empty([0,], np.float64)
        np.testing.assert_equal(dummy_ans, grid_dummy.getQuadratureWeights(), "Empty read", True)
        np.testing.assert_equal(dummy_ans, grid_dummy.getInterpolationWeights(aX), "Empty read", True)
        dummy_ans = np.empty([0, 0], np.float64)
        aX = np.empty([0, 0], np.float64)
        np.testing.assert_equal(dummy_ans, grid_dummy.getInterpolationWeightsBatch(aX), "Empty read", True)
        aX = np.empty([1, 0], np.float64)
        np.testing.assert_equal(dummy_ans, grid_dummy.getInterpolationWeightsBatch(aX), "Empty read", True)
        grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev')
        self.loadExpN2(grid)
        dummy_ans = np.empty([0, 1], np.float64)
        aX = np.empty([0, 2], np.float64)
        np.testing.assert_equal(dummy_ans, grid.evaluateBatch(aX), "Empty batch eval", True)

        # test weights
        aX = np.array([[0.0, -0.3], [-0.44, 0.7], [0.82, -0.01]])
        grid.makeGlobalGrid(2, 0, 4, 'level', 'chebyshev')
        aW = grid.getInterpolationWeightsBatch(aX)
        aX1 = aX[0,:].reshape([2,])
        aX2 = aX[1,:].reshape([2,])
        aX3 = aX[2,:].reshape([2,])
        aW1 = grid.getInterpolationWeights(aX1)
        aW2 = grid.getInterpolationWeights(aX2)
        aW3 = grid.getInterpolationWeights(aX3)
        aWW = np.row_stack([aW1, aW2, aW3])
        np.testing.assert_equal(aW, aWW, "Batch weights", True)

        grid.makeGlobalGrid(3, 0, 1, 'level', 'chebyshev')
        aTrans = np.array([[-1.0, -4.0], [2.0, 5.0], [-3.0, 5]])
        grid.setDomainTransform(aTrans)
        aA = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        aB = np.array([[-2.5, 3.5, 1.0], [-2.5, 3.5, -3.0], [-2.5, 3.5, 5.0], [-2.5,  2.0, 1.0], [-2.5, 5.0, 1.0], [-1.0, 3.5, 1.0], [-4.0, 3.5, 1.0]])
        np.testing.assert_equal(aB, grid.getPoints(), "Original equal", True)
        grid.clearDomainTransform()
        np.testing.assert_equal(aA, grid.getPoints(), "Tansformed equal", True)

        grid.makeGlobalGrid(3, 0, 1, 'level', 'fejer2')
        aA = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -0.707106781186548], [0.0, 0.0, 0.707106781186548], [0.0, -0.707106781186548, 0.0], [0.0, 0.707106781186548, 0.0], [-0.707106781186548, 0.0, 0.0], [0.707106781186548, 0.0, 0.0]])
        np.testing.assert_almost_equal(aA, grid.getPoints(), 14, "Original equal", True)
        grid.setConformalTransformASIN(np.array([4, 6, 3]))
        aA = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -0.60890205], [0.0, 0.0, 0.60890205], [0.0, -0.57892643, 0.0], [0.0, 0.57892643, 0.0], [-0.59587172, 0.0, 0.0], [0.59587172, 0.0, 0.0]])
        np.testing.assert_almost_equal(aA, grid.getPoints(), 7, "Transformed equal", True)
        grid.clearConformalTransform()
        aA = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -0.707106781186548], [0.0, 0.0, 0.707106781186548], [0.0, -0.707106781186548, 0.0], [0.0, 0.707106781186548, 0.0], [-0.707106781186548, 0.0, 0.0], [0.707106781186548, 0.0, 0.0]])
        np.testing.assert_almost_equal(aA, grid.getPoints(), 14, "Original equal", True)

        # Level Limits
        grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis', liLevelLimits = [1, 4])
        aPoints = grid.getPoints()
        self.checkPoints(aPoints[:,0], lMustHave = [-1.0, 1.0], lMustNotHave = [1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0)])
        self.checkPoints(aPoints[:,1], lMustHave = [-1.0, 1.0, 1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0)], lMustNotHave = [])
        self.loadExpN2(grid)
        grid.setAnisotropicRefinement('iptotal', 20, 0, [1, 4])
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        self.checkPoints(aPoints[:,0], lMustHave = [], lMustNotHave = [1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0)])

        grid.makeSequenceGrid(3, 1, 3, 'level', 'leja', liLevelLimits = [3, 2, 1])
        aPoints = grid.getPoints()
        self.checkPoints(aPoints[:,0], lMustHave = [0.0, -1.0, 1.0, 1.0/np.sqrt(3.0)], lMustNotHave = [])
        self.checkPoints(aPoints[:,1], lMustHave = [0.0, -1.0, 1.0], lMustNotHave = [1.0/np.sqrt(3.0)])
        self.checkPoints(aPoints[:,2], lMustHave = [0.0, 1.0], lMustNotHave = [-1.0, 1.0/np.sqrt(3.0)])
        self.loadExpN2(grid)
        grid.setAnisotropicRefinement('iptotal', 30, 0)
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        self.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [1.0/np.sqrt(3.0)])
        self.checkPoints(aPoints[:,2], lMustHave = [], lMustNotHave = [-1.0, 1.0/np.sqrt(3.0)])
        grid.clearRefinement()
        grid.setAnisotropicRefinement('iptotal', 30, 0, [3, 2, 2])
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        self.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [1.0/np.sqrt(3.0)])
        self.checkPoints(aPoints[:,2], lMustHave = [-1.0], lMustNotHave = [1.0/np.sqrt(3.0)])
        grid.clearRefinement()
        grid.setSurplusRefinement(1.E-8, 0, '', [3, 2, 1])
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        self.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [1.0/np.sqrt(3.0)])
        self.checkPoints(aPoints[:,2], lMustHave = [], lMustNotHave = [-1.0, 1.0/np.sqrt(3.0)])

        # check that nodes from level 2 (+-0.5) and 3 (+-0.25, +-0.75) appear only in the proper dimension
        grid.makeLocalPolynomialGrid(3, 1, 3, 1, 'localp', [1, 2, 3])
        aPoints = grid.getPoints()
        self.checkPoints(aPoints[:,0], lMustHave = [0.0, -1.0, 1.0], lMustNotHave = [0.5, -0.5, -0.75, -0.25, 0.25, 0.75])
        self.checkPoints(aPoints[:,1], lMustHave = [0.0, -1.0, 1.0, 0.5, -0.5], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])
        self.checkPoints(aPoints[:,2], lMustHave = [0.0, -1.0, 1.0, 0.5, -0.5, -0.75, -0.25, 0.25, 0.75], lMustNotHave = [])
        self.loadExpN2(grid)
        grid.setSurplusRefinement(1.E-8, 0, 'classic')
        self.assertTrue(grid.getNumNeeded() > 0, 'did not refine')
        aPoints = grid.getNeededPoints()
        self.checkPoints(aPoints[:,0], lMustHave = [], lMustNotHave = [0.5, -0.5, -0.75, -0.25, 0.25, 0.75])
        self.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])
        self.checkPoints(aPoints[:,2], lMustHave = [], lMustNotHave = [])
        grid.clearRefinement()
        grid.setSurplusRefinement(1.E-8, 0, 'classic', [2, 2, 3])
        aPoints = grid.getNeededPoints()
        self.checkPoints(aPoints[:,0], lMustHave = [0.5, -0.5], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])
        self.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])

        grid.makeWaveletGrid(2, 1, 3, 1, [0, 2])
        aPoints = grid.getPoints()
        self.checkPoints(aPoints[:,0], lMustHave = [0.0, -1.0, 1.0], lMustNotHave = [-0.5, 0.5])
        self.checkPoints(aPoints[:,1], lMustHave = [0.0, -1.0, 1.0, -0.5, 0.5, -0.75, -0.25, 0.25, 0.75], lMustNotHave = [-0.125, 0.125])
        self.loadExpN2(grid)
        grid.setSurplusRefinement(1.E-8, 0, 'classic') # grid is effectively full tensor, cannot refine within the level limits
        self.assertTrue((grid.getNumNeeded() == 0), 'added points even though it should not have');
        grid.setSurplusRefinement(1.E-8, 0, 'classic', [1, 2])
        aPoints = grid.getNeededPoints()
        self.checkPoints(aPoints[:,0], lMustHave = [0.5, -0.5], lMustNotHave = [-0.75, -0.25, 0.25, 0.75])
        self.checkPoints(aPoints[:,1], lMustHave = [], lMustNotHave = [-0.125, 0.125])

        # level limits I/O
        grid.makeLocalPolynomialGrid(3, 1, 3, 1, 'localp', [1, 2, 3])
        aLimits = grid.getLevelLimits()
        np.testing.assert_almost_equal(aLimits, np.array([1, 2, 3]), 14, "Could not read level limits", True)

        grid.clearLevelLimits()
        aLimits = grid.getLevelLimits()
        np.testing.assert_almost_equal(aLimits, np.array([-1, -1, -1]), 14, "Could not read level limits", True)

    def testFullCoverageB(self):
        print("\nTesting core refine grid")
        grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'semi-localp')
        self.loadExpN2(grid)
        self.assertEqual(grid.getNumNeeded(), 0, "num needed")
        grid.setSurplusRefinement(0.0001, 0, 'classic')
        self.assertTrue((grid.getNumNeeded() > 0), "num needed")
        grid.clearRefinement()
        self.assertEqual(grid.getNumNeeded(), 0, "num needed")

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

    def testFullCoverageC(self):
        print("\nTesting core learning from random samples")

        # evalHierarchicalBasis (all done in batch)
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
        aT = np.array([[1.0, 0.0], [math.exp(-1.0)-1.0, math.sin(-1.0)], [math.exp(1.0)-1.0, math.sin(1.0)], [math.exp(-2.0)-1.0, math.sin(-3.0)], [math.exp(2.0)-1.0, math.sin(3.0)]])
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

        # sparse hierarchical functions
        pSparse = TasmanianSG.TasmanianSimpleSparseMatrix()
        np.testing.assert_almost_equal(np.empty([0,0], np.float64), pSparse.getDenseForm(), 14, "TasmanianSimpleSparseMatrix.getDense()", True)

        grid.makeGlobalGrid(2, 1, 4, 'iptotal', 'fejer2')
        aPoints = np.array([[0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666]])
        aDense = grid.evaluateHierarchicalFunctions(aPoints)
        pSparse = grid.evaluateSparseHierarchicalFunctions(aPoints)
        np.testing.assert_almost_equal(aDense, pSparse.getDenseForm(), 14, "evaluateSparseHierarchicalFunctions", True)

        grid.makeSequenceGrid(2, 1, 4, 'iptotal', 'min-delta')
        aPoints = np.array([[0.33, 0.25], [-0.27, 0.39], [0.97, -0.76], [-0.44, 0.21], [-0.813, 0.03], [-0.666, 0.666]])
        aDense = grid.evaluateHierarchicalFunctions(aPoints)
        pSparse = grid.evaluateSparseHierarchicalFunctions(aPoints)
        np.testing.assert_almost_equal(aDense, pSparse.getDenseForm(), 14, "evaluateSparseHierarchicalFunctions", True)

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

        gridA = TasmanianSG.TasmanianSparseGrid()
        gridB = TasmanianSG.TasmanianSparseGrid()
        gridA.makeGlobalGrid(2, 1, 4, 'level', 'chebyshev')
        gridB.makeGlobalGrid(2, 1, 4, 'level', 'chebyshev')
        self.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        self.compareGrids(gridA, gridB)

        gridA = TasmanianSG.TasmanianSparseGrid()
        gridB = TasmanianSG.TasmanianSparseGrid()
        gridA.makeGlobalGrid(2, 1, 4, 'level', 'fejer2')
        gridB.makeGlobalGrid(2, 1, 4, 'level', 'fejer2')
        self.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        self.compareGrids(gridA, gridB)
        gridA.setAnisotropicRefinement('iptotal', 10, 0)
        gridB.setAnisotropicRefinement('iptotal', 10, 0)
        self.loadExpN2(gridA)
        gridB.mergeRefinement()
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), 'different number of points after merge refinement')
        aRes = gridB.evaluateBatch(aPoints)
        np.testing.assert_almost_equal(aRes, np.zeros(aRes.shape), 14, "not zero after merged refinement", True)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        self.compareGrids(gridA, gridB)

        gridA.makeSequenceGrid(2, 1, 5, 'ipcurved', 'min-delta')
        gridB.makeSequenceGrid(2, 1, 5, 'ipcurved', 'min-delta')
        self.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        self.compareGrids(gridA, gridB)
        gridA.setAnisotropicRefinement('iptotal', 30, 0)
        gridB.setAnisotropicRefinement('iptotal', 30, 0)
        self.loadExpN2(gridA)
        gridB.mergeRefinement()
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), 'different number of points after merge refinement')
        aRes = gridB.evaluateBatch(aPoints)
        np.testing.assert_almost_equal(aRes, np.zeros(aRes.shape), 14, "not zero after merged refinement", True)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        self.compareGrids(gridA, gridB)

        gridA.makeLocalPolynomialGrid(2, 1, 4, 2, 'semi-localp')
        gridB.makeLocalPolynomialGrid(2, 1, 4, 2, 'semi-localp')
        self.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        self.compareGrids(gridA, gridB)
        gridA.setSurplusRefinement(1.E-4, 0, 'classic')
        gridB.setSurplusRefinement(1.E-4, 0, 'classic')
        self.loadExpN2(gridA)
        gridB.mergeRefinement()
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), 'different number of points after merge refinement')
        aRes = gridB.evaluateBatch(aPoints)
        np.testing.assert_almost_equal(aRes, np.zeros(aRes.shape), 14, "not zero after merged refinement", True)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        self.compareGrids(gridA, gridB)

        gridA.makeWaveletGrid(2, 1, 2, 1)
        gridB.makeWaveletGrid(2, 1, 2, 1)
        self.loadExpN2(gridA)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        self.compareGrids(gridA, gridB)
        gridA.setSurplusRefinement(1.E-3, 0, 'classic')
        gridB.setSurplusRefinement(1.E-3, 0, 'classic')
        self.loadExpN2(gridA)
        gridB.mergeRefinement()
        self.assertEqual(gridA.getNumPoints(), gridB.getNumPoints(), 'different number of points after merge refinement')
        aRes = gridB.evaluateBatch(aPoints)
        np.testing.assert_almost_equal(aRes, np.zeros(aRes.shape), 14, "not zero after merged refinement", True)
        gridB.setHierarchicalCoefficients(gridA.getHierarchicalCoefficients())
        self.compareGrids(gridA, gridB)

    def testFullCoverageZ(self):
        print("\nTesting plotting and other misc")
        if (TasmanianSG.bTsgPlotting):
            if ((os.name == "posix") and (os.environ.get('DISPLAY') is None)):
                print("NOTE: there is no display, cannot run plotting tests.")
                return
            grid.makeGlobalGrid(2, 1, 20, 'level', 'leja')
            self.loadExpN2(grid)
            grid.plotResponse2D()
            grid.plotPoints2D()

            grid.makeGlobalGrid(2, 1, 0, 'level', 'leja')
            self.loadExpN2(grid)
            grid.plotResponse2D()
            grid.plotPoints2D()

            grid.makeGlobalGrid(3, 1, 3, 'level', 'leja')
            try:
                grid.plotPoints2D()
                self.assertTrue(False, "failed to flag plot exception when when using a grid with other than 2D")
            except TasmanianSG.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, "plotPoints2D", "error raising exception for plotPoints2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotResponse2D()
                self.assertTrue(False, "failed to flag plot exception when when using a grid with other than 2D")
            except TasmanianSG.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, "plotResponse2D", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))

            grid.makeGlobalGrid(2, 2, 2, 'level', 'leja')
            try:
                grid.plotResponse2D(iOutput = -1)
                self.assertTrue(False, "failed to flag plot exception when when using wrong output")
            except TasmanianSG.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, "iOutput", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotResponse2D(iOutput = 3)
                self.assertTrue(False, "failed to flag plot exception when when using wrong output")
            except TasmanianSG.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, "iOutput", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotResponse2D(iNumDim0 = -1)
                self.assertTrue(False, "failed to flag plot exception when when using wrong output")
            except TasmanianSG.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, "iNumDim0", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotResponse2D(iNumDim1 = -1)
                self.assertTrue(False, "failed to flag plot exception when when using wrong output")
            except TasmanianSG.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, "iNumDim1", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
        else:
            grid.makeGlobalGrid(2, 1, 3, 'level', 'leja')
            self.loadExpN2(grid)
            try:
                grid.plotResponse2D()
                self.assertTrue(False, "failed to flag plotexception  when missing matplotlib")
            except TasmanianSG.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, "plotResponse2D", "error raising exception for plotResponse2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))
            try:
                grid.plotPoints2D()
                self.assertTrue(False, "failed to flag plot exception when missing matplotlib")
            except TasmanianSG.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, "plotPoints2D", "error raising exception for plotPoints2D() using test\n Error.sVariable = '{0:1s}'".format(TSGError.sVariable))

if __name__ == '__main__':
    unittest.main()


