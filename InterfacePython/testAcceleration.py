import unittest
import TasmanianSG
import sys, os
import numpy as np

from random import uniform

import testConfigureData as tdata # needed to keep track of configured acceleration methods
import testCommon

ttc = testCommon.TestTasCommon()

class TestTasClass(unittest.TestCase):
    '''
    Test the acceleration options:
    * consistency between CMake options and options reported by C++/Python
    * ability to switch between GPUs, read properties set GPUID, etc.
    * consistency between accelerated and reference evaluations
    '''
    def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")

    def testNothing(self):
        pass

    def checkMeta(self):
        '''
        Check if the CMake options are propagated to the C++ library and the
        correct acceleration options are passed between C++ and Python.
        '''
        grid = TasmanianSG.TasmanianSparseGrid()

        if (tdata.bEnableSyncTests):
            self.assertTrue((grid.isAccelerationAvailable("cpu-blas") == tdata.bHasBlas), "failed to match blas")
            self.assertTrue((grid.isAccelerationAvailable("gpu-cublas") == tdata.bHasCuBlas), "failed to match cublas")
            self.assertTrue((grid.isAccelerationAvailable("gpu-cuda") == tdata.bHasCuda), "failed to match cuda")
            self.assertTrue((grid.isAccelerationAvailable("gpu-default") == (tdata.bHasCuBlas or tdata.bHasCuda)), "failed to match cuda")

            # acceleration meta-data
            lsAvailableAcc = []
            if (tdata.bHasBlas): lsAvailableAcc.append("cpu-blas")
            if (tdata.bHasCuBlas): lsAvailableAcc.append("gpu-cublas")
            if (tdata.bHasCuda): lsAvailableAcc.append("gpu-cuda")
            grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'semi-localp')
            for accel in lsAvailableAcc:
                grid.enableAcceleration(accel)
                sA = grid.getAccelerationType()
                bTest = ((accel not in sA) or (sA not in accel))
                self.assertFalse(bTest, "set/get Acceleration")

            lsAvailableAcc = []
            if (tdata.bHasCuBlas): lsAvailableAcc.append(("gpu-rocblas", "gpu-cublas"))
            if (tdata.bHasCuda): lsAvailableAcc.append(("gpu-hip", "gpu-cuda"))
            for accel in lsAvailableAcc:
                grid.enableAcceleration(accel[0])
                sA = grid.getAccelerationType()
                bTest = ((accel[1] not in sA) or (sA not in accel[1]))
                self.assertFalse(bTest, "set/get Acceleration with HIP - CUDA aliases")


    def checkMultiGPU(self):
        '''
        Check setting and resetting the GPU ID and reading the device names.
        '''
        grid = TasmanianSG.TasmanianSparseGrid()

        self.assertTrue((grid.getGPUID() == 0), "did not default to gpu 0")

        if (grid.getNumGPUs() > 1):
            grid.setGPUID(1)
            self.assertTrue((grid.getGPUID() == 1), "did not set to gpu 1")
            if (not tdata.bUsingMSVC):
                sName = grid.getGPUName(1) # mostly checks for memory leaks and crashes

        if (grid.getNumGPUs() > 0):
            grid.setGPUID(0)
            self.assertTrue((grid.getGPUID() == 0), "did not set to gpu 0")
            if (not tdata.bUsingMSVC):
                sName = grid.getGPUName(0) # mostly checks for memory leaks and crashes

    def checkEvaluateConsistency(self):
        '''
        Check for consistency between accelerated and reference evaluations.
        In short, set a grid, load points, evaluate in different ways.
        Test all visible GPUs and combinations cuBlas/MAGMA, etc.
        '''
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
                   'grid.makeWaveletGrid(2, 1, 3, 3)',
                   'grid.makeFourierGrid(2, 1, 3, "level")' ]
        lTests = ['grid.makeLocalPolynomialGrid(2, 3, 4, 1, "localp")']

        iNumGPUs = grid.getNumGPUs()
        lsAccelTypes = ["none", "cpu-blas", "gpu-cuda", "gpu-cublas", "gpu-magma"]

        for sTest in lTests:
            for iI in range(2):
                iC = 0
                iGPU = 0
                if (tdata.iGPUID > -1):
                    iGPU = tdata.iGPUID

                while (iC < len(lsAccelTypes)):
                    sAcc = lsAccelTypes[iC]
                    exec(sTest)
                    if (iI == 0):
                        aTestPoints = aTestPointsCanonical
                    else:
                        aTestPoints = aTestPointsTransformed
                        grid.setDomainTransform(aDomainTransform)

                    grid.enableAcceleration(sAcc)
                    if ((sAcc == "gpu-cuda") or (sAcc == "gpu-cublas") or (sAcc == "gpu-magma")):
                        if (iGPU < grid.getNumGPUs()): # without cuda or cublas, NumGPUs is 0 and cannot set GPU
                            grid.setGPUID(iGPU)

                    ttc.loadExpN2(grid)

                    aRegular = np.array([grid.evaluateThreadSafe(aTestPoints[i,:]) for i in range(aTestPoints.shape[0]) ])
                    aBatched = grid.evaluateBatch(aTestPoints)
                    np.testing.assert_almost_equal(aRegular, aBatched, 14, "Batch evaluation test not equal: {0:1s}, acceleration: {1:1s}, gpu: {2:1d}".format(sTest, sAcc, iGPU), True)

                    aFast = np.array([ grid.evaluate(aTestPoints[i,:]) for i in range(iFastEvalSubtest) ])
                    np.testing.assert_almost_equal(aRegular[0:iFastEvalSubtest,:], aFast, 14, "Batch evaluation test not equal: {0:1s}, acceleration: {1:1s}, gpu: {2:1d}".format(sTest, sAcc, iGPU), True)

                    if ((sAcc == "gpu-cuda") or (sAcc == "gpu-cublas") or (sAcc == "gpu-magma")):
                        if (tdata.iGPUID == -1):
                            iGPU += 1
                            if (iGPU >= iNumGPUs):
                                iC += 1
                                iGPU = 0
                        else:
                            iC += 1
                    else:
                        iC += 1

    def performAccelerationTest(self):
        self.checkMeta()
        self.checkMultiGPU()
        self.checkEvaluateConsistency()
