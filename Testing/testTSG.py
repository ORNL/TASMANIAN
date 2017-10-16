#!/usr/bin/env python

import unittest
import TasmanianSG
import math
import numpy as np
from random import uniform

grid = TasmanianSG.TasmanianSparseGrid()

# in principle, should always test the wavelets, but it takes long
_TestWavelets = True

class TestTasmanian(unittest.TestCase):
        
    def tsgReadMatrix(self, sFilename):
        with open(sFilename, 'r') as infile:
            lsLines = infile.readlines()
        
        sLine = lsLines.pop(0)
        liDims = [int(s) for s in sLine.split() if s.isdigit()]
        
        if (len(liDims) != 2):
            return []
            
        M = np.empty(liDims)
        iI = 0
        for sLine in lsLines:
            lfRow = [float(s) for s in sLine.split() if s.isdigit()]
            iJ = 0
            for fV in lfRow:
                M[iI,iJ] = fV
                iJ += 1
                
        return M
        
    def compareGrids(self, gridA, gridB):
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
            #print "{0:2.20e}   {1:2.20e}".format(pA[0], pB[0])
            #np.testing.assert_equal(pA, pB, "Interpolation test 4 not equal", True)
            np.testing.assert_almost_equal(pA, pB, 15, "Interpolation test 4 not equal", True)
            
            pA = gridA.evaluate(mX5)
            pB = gridB.evaluate(mX5)
            np.testing.assert_almost_equal(pA, pB, 15, "Interpolation test 5 not equal", True)
            
            pA = gridA.integrate()
            pB = gridB.integrate()
            np.testing.assert_almost_equal(pA, pB, 15, "Integration test not equal", True)
        
        self.assertEqual(gridA.isGlobal(), gridB.isGlobal(), "error in isGlobal()")
        self.assertEqual(gridA.isSequence(), gridB.isSequence(), "error in isSequence()")
        self.assertEqual(gridA.isLocalPolynomial(), gridB.isLocalPolynomial(), "error in isLocalPolynomial()")
        self.assertEqual(gridA.isWavelet(), gridB.isWavelet(), "error in isWavelet()")
        
        self.assertEqual(gridA.isSetDomainTransfrom(), gridB.isSetDomainTransfrom(), "error in isSetDomainTransfrom()")
        
        pA = gridA.getDomainTransform()
        pB = gridB.getDomainTransform()
        np.testing.assert_equal(pA, pB, "Domain test no equal", True)
        
        self.assertEqual(gridA.isSetConformalTransformASIN(), gridB.isSetConformalTransformASIN(), "error in isSetConformalTransformASIN()")

        pA = gridA.getConformalTransformASIN()
        pB = gridB.getConformalTransformASIN()
        np.testing.assert_equal(pA, pB, "Conformal transform ASIN not equal", True)
        
    def loadExpN2(self, grid):
        if (grid.getNumNeeded() == 0): return
        mPoints = grid.getNeededPoints()
        iOut = grid.getNumOutputs()
        iDim = grid.getNumDimensions()
        lVals = [[math.exp(-sum([mPoints[i,j]**2 for j in range(iDim)])) for k in range(iOut)] for i in range(mPoints.shape[0])]
        grid.loadNeededPoints(np.array(lVals))
    
    def testAABasicLibMeta(self): # the AA is so that this test goes first
        print("\nLibrary core check, read the library meta")
        grid = TasmanianSG.TasmanianSparseGrid()
        print("Tasmanian Sparse Grids version: {0:1s}".format(grid.getVersion()))
        print("                 version major: {0:1d}".format(grid.getVersionMajor()))
        print("                 version minor: {0:1d}".format(grid.getVersionMinor()))
        print("                       License: {0:1s}".format(grid.getLicense()))
        if (grid.isCudaEnabled()):
            sStatus = "Enabled"
        else:
            sStatus = "Disabled"
        print("                   Nvidia CUDA: {0:1s}".format(sStatus))
        if (grid.isBLASEnabled()):
            sStatus = "Enabled"
        else:
            sStatus = "Disabled"
        print("                          BLAS: {0:1s}".format(sStatus))
        if (grid.isOpenMPEnabled()):
            sStatus = "Enabled"
        else:
            sStatus = "Disabled"
        print("                        OpenMP: {0:1s}".format(sStatus))
    
    def testIO(self):
        print("\nRunning core I/O test")
        # test I/O for Global Grids
        # iDimension, iOutputs, iDepth, sType, sRule, fAlpha, fBeta, useTransform, loadFunciton
        lGrids = [[3, 2, 2, "level", "leja", 0.0, 0.0, False, False],
                  [2, 1, 4, "level", "clenshaw-curtis", 0.0, 0.0, False, False],
                  [3, 1, 3, "level", "rleja", 0.0, 0.0, True, False],
                  [3, 1, 2, "iptotal", "chebyshev", 0.0, 0.0, False, True],
                  [3, 1, 3, "level", "leja", 0.0, 0.0, True, True],
                  [2, 1, 5, "qptotal", "gauss-hermite", 1.0, 3.0, False, False],
                  [2, 1, 3, "level", "gauss-laguerre", 3.0, 0.0, False, False],]
        aTransform = np.array([[0.0,1.0],[0.0,1.0],[-2.0,-1.0]])
        
        for lT in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            gridA.makeGlobalGrid(lT[0], lT[1], lT[2], lT[3], lT[4], [], lT[5], lT[6])
            if (lT[7]):
                gridA.setDomainTransform(aTransform)
            if (lT[8]):
                self.loadExpN2(gridA)
            
            gridA.write("testSave")
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)
            
        # test an error message from wrong read
        #print("Attempting a bogus read to see if error would be properly registered")
        gridB.disableLog()
        self.assertFalse(gridB.read("testSaveBlah"), "Failed to flag a fake read")
        gridB.setErrorLogCerr()
        
        # test I/O for Sequence Grids
        # iDimension, iOutputs, iDepth, sType, sRule, useTransform, loadFunciton
        lGrids = [[3, 2, 2, "level", "leja", False, False],
                  [2, 1, 4, "level", "max-lebesgue", False, False],
                  [3, 1, 3, "level", "rleja", True, False],
                  [3, 1, 2, "iptotal", "min-delta", False, True],
                  [3, 1, 3, "level", "leja", True, True],]
                  
        for lT in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            gridA.makeSequenceGrid(lT[0], lT[1], lT[2], lT[3], lT[4])
            if (lT[5]):
                gridA.setDomainTransform(aTransform)
            if (lT[6]):
                self.loadExpN2(gridA)
            
            gridA.write("testSave")
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)
            
        # test I/O for Local Polynomial Grids
        # iDimension, iOutputs, iDepth, sType, sRule, useTransform, loadFunciton
        lGrids = [[3, 2, 2, 0, "localp", False, False],
                  [2, 1, 4, 1, "semi-localp", False, False],
                  [3, 1, 3, 2, "localp", True, False],
                  [3, 1, 2, 3, "localp-zero", False, True],
                  [3, 1, 3, 4, "semi-localp", True, True],
                  [3, 1, 3, -1, "semi-localp", True, True],]
                  
        for lT in lGrids:
            gridA = TasmanianSG.TasmanianSparseGrid()
            gridB = TasmanianSG.TasmanianSparseGrid()

            gridA.makeLocalPolynomialGrid(lT[0], lT[1], lT[2], lT[3], lT[4])
            if (lT[5]):
                gridA.setDomainTransform(aTransform)
            if (lT[6]):
                self.loadExpN2(gridA)
            
            gridA.write("testSave")
            gridB.read("testSave")
            self.compareGrids(gridA, gridB)
            
        # test I/O for Local Wavelet Grids
        # iDimension, iOutputs, iDepth, sType, sRule, useTransform, loadFunciton
        if (_TestWavelets):
            lGrids = [[3, 2, 2, 1, False, False],
                      [2, 1, 4, 1, False, False],
                      [3, 1, 3, 3, True, False],
                      [3, 1, 2, 1, False, True],
                      [3, 1, 3, 3, True, True],]
            
            for lT in lGrids:
                gridA = TasmanianSG.TasmanianSparseGrid()
                gridB = TasmanianSG.TasmanianSparseGrid()
            
                gridA.makeWaveletGrid(lT[0], lT[1], lT[2], lT[3])
                if (lT[4]):
                    gridA.setDomainTransform(aTransform)
                if (lT[5]):
                    self.loadExpN2(gridA)
                
                gridA.write("testSave")
                gridB.read("testSave")
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
        
                
    def testAcceleratedEvaluate(self):
        print("\nAccelerated evaluate consistency test")
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
        
        iNumGPUs = grid.getNumGPUs()
        lsAccelTypes = ["none", "cpu-blas", "gpu-fullmem"]
        
        for sTest in lTests:
            for iI in range(2):
                iC = 0
                iGPU = 0
                while (iC < len(lsAccelTypes)):
                    sAcc = lsAccelTypes[iC]
                    exec(sTest)
                    if (iI == 0):
                        aTestPoints = aTestPointsCanonical
                    else:
                        aTestPoints = aTestPointsTransformed
                        grid.setDomainTransform(aDomainTransform)
                    
                    grid.enableAcceleration(sAcc)
                    if (sAcc == "gpu-fullmem"):
                        grid.setGPUID(iGPU)
                    
                    self.loadExpN2(grid)
                    
                    aRegular = np.array([ grid.evaluateThreadSafe(aTestPoints[i,:]) for i in range(aTestPoints.shape[0]) ])
                    aBatched = grid.evaluateBatch(aTestPoints)
                    np.testing.assert_almost_equal(aRegular, aBatched, 14, "Batch evaluation test not equal: {0:1s}, acceleration: {1:1s}, gpu: {2:1d}".format(sTest, sAcc, iGPU), True)
                    
                    aFast = np.array([ grid.evaluate(aTestPoints[i,:]) for i in range(iFastEvalSubtest) ])
                    np.testing.assert_almost_equal(aRegular[0:iFastEvalSubtest,:], aFast, 14, "Batch evaluation test not equal: {0:1s}, acceleration: {1:1s}, gpu: {2:1d}".format(sTest, sAcc, iGPU), True)
                    
                    if (sAcc == "gpu-fullmem"):
                        iGPU += 1
                        if (iGPU >= iNumGPUs):
                            iC += 1
                            iGPU = 0
                    else:
                        iC += 1

    def testBasicException(self):
        print("\nTesting error handling")
        grid = TasmanianSG.TasmanianSparseGrid()
        #print("\nAttempting bogus grid construction, should see many errors")
        
        llTests = [["grid.makeGlobalGrid(-1, 1,  4, 'level', 'clenshaw-curtis')", "iDimension"],
                   ["grid.makeGlobalGrid(2, -1,  4, 'level', 'clenshaw-curtis')", "iOutputs"],
                   ["grid.makeGlobalGrid(2,  1, -4, 'level', 'clenshaw-curtis')", "iDepth"],
                   ["grid.makeGlobalGrid(2,  1,  4, 'wrong', 'clenshaw-curtis')", "sType"],
                   ["grid.makeGlobalGrid(2,  1,  4, 'level', 'clenshaw-wrong')", "sRule"],
                   ["grid.makeGlobalGrid(2,  1,  4, 'level', 'clenshaw-curtis', [1,2,3])", "liAnisotropicWeights"],
                   ["grid.makeSequenceGrid(-1, 1,  4, 'level', 'leja')", "iDimension"],
                   ["grid.makeSequenceGrid(2, -1,  4, 'level', 'leja')", "iOutputs"],
                   ["grid.makeSequenceGrid(2,  1, -4, 'level', 'leja')", "iDepth"],
                   ["grid.makeSequenceGrid(2,  1,  4, 'wrong', 'leja')", "sType"],
                   ["grid.makeSequenceGrid(2,  1,  4, 'level', 'weja')", "sRule"],
                   ["grid.makeSequenceGrid(2,  1,  4, 'level', 'leja', [1,2,3])", "liAnisotropicWeights"],
                   ["grid.makeLocalPolynomialGrid(-1, 1,  4,  2, 'localp')", "iDimension"],
                   ["grid.makeLocalPolynomialGrid(2, -1,  4,  2, 'localp')", "iOutputs"],
                   ["grid.makeLocalPolynomialGrid(2,  1, -4,  2, 'localp')", "iDepth"],
                   ["grid.makeLocalPolynomialGrid(2,  1,  4, -2, 'localp')", "iOrder"],
                   ["grid.makeLocalPolynomialGrid(2,  1,  4,  2, 'lowrong')", "sRule"],
                   ["grid.makeWaveletGrid(-1, 1,  4,  1)", "iDimension"],
                   ["grid.makeWaveletGrid(2, -1,  4,  1)", "iOutputs"],
                   ["grid.makeWaveletGrid(2,  1, -4,  3)", "iDepth"],
                   ["grid.makeWaveletGrid(2,  1,  4,  2)", "iOrder"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja')", "notError"], # notError tests here are needed to ensure that the multi-statement commands fail for the right function
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev')", "notError"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2]))", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateGlobalGrid(1,'iptotal')", "updateGlobalGrid"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev'); grid.updateGlobalGrid(-1,'iptotal')", "iDepth"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev'); grid.updateGlobalGrid(4,'wrong')", "sType"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev'); grid.updateGlobalGrid(4,'iptotal',[1,2,3])", "liAnisotropicWeights"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(4,'iptotal')", "updateSequenceGrid"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(-1,'iptotal')", "iDepth"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(4,'wrong')", "sType"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateSequenceGrid(4,'iptotal',[1,2,3])", "liAnisotropicWeights"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); aW = grid.getInterpolationWeights(np.array([1,2,3]))", "lfX"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); aW = grid.getInterpolationWeightsBatch(np.array([1,2,3]))", "llfX"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); aW = grid.getInterpolationWeightsBatch(np.array([[1,2,3],[1,2,3]]))", "llfX"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([1,]))", "llfVals"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,3]))", "llfVals"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([5,2]))", "llfVals"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2])); grid.loadNeededPoints(np.ones([5,2]))", "llfVals"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.evaluate(np.zeros([1,2]))", "evaluate"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededPoints(np.zeros([6,2])); grid.evaluate(np.zeros([1,3]))", "lfX"],
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
                   ["grid.makeGlobalGrid(3, 1, 2, 'level', 'clenshaw-curtis'); grid.setConformalTransformASIN(np.array([[0,2,3],[1,2,3]]))", "liTruncation"],]

        for lTest in llTests:
            try:
                exec(lTest[0])
                self.assertEqual(lTest[1], "notError", "failed to raise exception for invalid '{0:1s}' using test\n '{1:1s}'".format(lTest[1],lTest[0]))
            except TasmanianSG.TasmanianInputError as TSGError: 
                self.assertEqual(TSGError.sVariable, lTest[1], "error raising exception for '{0:1s}' using test\n '{1:1s}'\n Error.sVariable = '{2:1s}'".format(lTest[1],lTest[0],TSGError.sVariable))
                 
        #print("End of the many errors (more errors could follow, but not without warning)\n")

if __name__ == '__main__':
    unittest.main()

