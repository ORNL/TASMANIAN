import unittest
import Tasmanian
DREAM = Tasmanian.DREAM
Opt = Tasmanian.Optimization
import numpy as np

import testCommon

ttc = testCommon.TestTasCommon()

# tests should pass when errors are raised, thus suppress the verbose error messages
# to avoid confusing a verbose message with a real error (i.e., failed test)
Tasmanian.useVerboseErrors(False)

class TestTasClass(unittest.TestCase):
    '''
    Test the python exceptions, create a bunch of incorrect faults and
    check whether the correct exception is raised.
    '''
    def __init__(self):
        unittest.TestCase.__init__(self, "testNothing")

    def testNothing(self):
        pass

    def getSparseGridTests(self):
        # The list llTest contains list-of-lists (could make into tuples)
        # Each sub-list (tuple) has two entries, one is the command that should raise the exception,
        # the second is the "sError" variable in the exception class.
        # notError tests here are needed to ensure that the multi-statement commands fail for the correct function.
        return    [["grid.makeGlobalGrid(-1, 1,  4, 'level', 'clenshaw-curtis')", "iDimension"],
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
                   ["grid.makeFourierGrid(-1, 1,  4, 'level')", "iDimension"],
                   ["grid.makeFourierGrid(2, -1,  4, 'level')", "iOutputs"],
                   ["grid.makeFourierGrid(2,  1, -4, 'level')", "iDepth"],
                   ["grid.makeFourierGrid(2,  1,  4, 'wrong')", "sType"],
                   ["grid.makeFourierGrid(2,  1,  4, 'level', [1, 2, 3])", "liAnisotropicWeights"],
                   ["grid.makeFourierGrid(2,  1,  4, 'level', [1, 2], [1, 2, 3])", "liLevelLimits"],
                   ["grid.makeFourierGrid(2,  1,  4, 'level', liLevelLimits = [1, 2])", "notError"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja')", "notError"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'chebyshev')", "notError"],
                   ["grid.makeSequenceGrid(2, 2, 2, 'level', 'rleja'); grid.loadNeededValues(np.zeros([6,2]))", "notError"],
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
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'rleja'); grid.updateFourierGrid(3,'level')", "updateFourierGrid"],
                   ["grid.makeFourierGrid(2, 1, 2, 'level'); grid.updateFourierGrid(-3,'level')", "iDepth"],
                   ["grid.makeFourierGrid(2, 1, 2, 'level'); grid.updateFourierGrid(3,'wrong')", "sType"],
                   ["grid.makeFourierGrid(2, 1, 2, 'level'); grid.updateFourierGrid(3,'iptotal',[1])", "liAnisotropicWeights"],
                   ["grid.makeFourierGrid(2, 1, 2, 'level'); grid.updateFourierGrid(3,'iptotal',liLevelLimits = [4])", "liLevelLimits"],
                   ["grid.makeFourierGrid(2, 1, 2, 'level'); grid.updateFourierGrid(3,'iptotal', [4, 3], [4, 3])", "notError"],
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
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'fejer2'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', -2, 1);", "iMinGrowth"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -1);", "iOutput"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, 0);", "notError"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, 0, [2, 3]);", "notError"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, 0, [2, 3, 3]);", "liLevelLimits"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -1);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -2);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, 5);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('wrong', 10, 0);", "sType"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -1, [3, 4]);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.setAnisotropicRefinement('iptotal', 10, -1, [3, 4, 5]);", "liLevelLimits"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'rleja'); grid.estimateAnisotropicCoefficients('iptotal', 1);", "estimateAnisotropicCoefficients"],
                   ["grid.makeGlobalGrid(2, 0, 2, 'level', 'clenshaw-curtis'); grid.estimateAnisotropicCoefficients('iptotal', 1);", "estimateAnisotropicCoefficients"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); ttc.loadExpN2(grid); grid.estimateAnisotropicCoefficients('iptotal', -1);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.estimateAnisotropicCoefficients('iptotal', -1);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.estimateAnisotropicCoefficients('ipcurved', -1);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.estimateAnisotropicCoefficients('iptotal', -2);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.estimateAnisotropicCoefficients('iptotal', 5);", "iOutput"],
                   ["grid.makeSequenceGrid(2, 1, 3, 'iptotal', 'leja'); ttc.loadExpN2(grid); grid.estimateAnisotropicCoefficients('wrong', 0);", "sType"],
                   ["grid.makeGlobalGrid(2, 1, 2, 'level', 'clenshaw-curtis'); grid.setSurplusRefinement(1.E-4, 0, 'classic');", "setSurplusRefinement"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.setSurplusRefinement(1.E-4, 0, 'classic');", "setSurplusRefinement"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); ttc.loadExpN2(grid); grid.setSurplusRefinement(-1.E-4, 0, 'classic');", "fTolerance"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic');", "sCriteria"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, '', [2, 3, 4]);", "liLevelLimits"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, '', [2, 3]);", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0);", "sCriteria"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'class');", "sCriteria"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic');", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic', [2, 3, 4]);", "liLevelLimits"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic', [], np.ones([3]));", "llfScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic', [], np.ones([grid.getNumPoints() - 1, 1]));", "llfScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic', [], np.ones([grid.getNumPoints(), 2]));", "llfScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 2, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, -1, 'classic', [], np.ones([grid.getNumPoints(), 2]));", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 2, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, -1, 'classic', [], np.ones([grid.getNumPoints(), 3]));", "llfScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.setSurplusRefinement(1.E-4, 0, 'classic', [2, 3]);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.removePointsByHierarchicalCoefficient(1.E-4, 0);", "removePointsByHierarchicalCoefficient"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.removePointsByHierarchicalCoefficient(-1.E-4, 0);", "fTolerance"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.removePointsByHierarchicalCoefficient(1.E-4, -2);", "iOutput"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.removePointsByHierarchicalCoefficient(1.E-4, 3);", "iOutput"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); grid.removePointsByHierarchicalCoefficient(1.E-4, 0);", "removePointsByHierarchicalCoefficient"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 2, 1, 'localp'); ttc.loadExpN2(grid); grid.removePointsByHierarchicalCoefficient(1.E-4, 0);", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([3,]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([10,]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([11,2]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([13,2]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, -1, np.ones([13,3]));", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, 0,  np.ones([13,3]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, 0,  np.ones([11,]));", "aScaleCorrection"],
                   ["grid.makeLocalPolynomialGrid(2, 3, 2, 1, 'localp'); grid.loadNeededPoints(np.ones([grid.getNumNeeded(), grid.getNumOutputs()])); grid.removePointsByHierarchicalCoefficient(1.E-4, 0,  np.ones([13,1]));", "notError"],
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
                   ["grid.makeFourierGrid(2, 1, 1, 'level'); grid.setHierarchicalCoefficients(np.array([[1.0,],[1.0,],[1.0,],[1.0,],[1.0,]]));", "llfCoefficients"],
                   ["grid.makeGlobalGrid(2, 1, 3, 'level', 'rleja'); grid.getCandidateConstructionPoints('level', -1);", "getCandidateConstructionPoints"],
                   ["grid.makeGlobalGrid(2, 1, 3, 'level', 'rleja'); grid.beginConstruction(); grid.getCandidateConstructionPoints('lev', -1);", "sType"],
                   ["grid.makeGlobalGrid(2, 1, 3, 'level', 'rleja'); grid.beginConstruction(); grid.getCandidateConstructionPoints('level', [0]);", "liAnisotropicWeightsOrOutput"],
                   ["grid.makeGlobalGrid(2, 1, 3, 'level', 'rleja'); grid.beginConstruction(); grid.getCandidateConstructionPoints('level', 'string');", "liAnisotropicWeightsOrOutput"],
                   ["grid.makeGlobalGrid(2, 1, 3, 'level', 'rleja'); grid.beginConstruction(); grid.getCandidateConstructionPoints('level', -1, [2]);", "liLevelLimits"],
                   ["grid.makeGlobalGrid(2, 1, 3, 'level', 'rleja'); grid.beginConstruction(); grid.getCandidateConstructionPoints('level', -1, [2, 1]);", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.getCandidateConstructionPointsSurplus(1.E-5, 'classic');", "getCandidateConstructionPointsSurplus"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.getCandidateConstructionPointsSurplus(1.E-5, 'classic', 0, [2]);", "liLevelLimits"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.getCandidateConstructionPointsSurplus(1.E-5, 'classic', -1, [2, 3]);", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.loadConstructedPoint([0.0, 0.0], [1.0]);", "notError"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.loadConstructedPoint([0.0, 0.0], [1.0]);", "loadConstructedPoint"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.loadConstructedPoint([0.0], [1.0]);", "lfX"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.loadConstructedPoint([0.0, 0.0], [1.0, 2.0]);", "lfY"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.loadConstructedPoint([[0.0, 0.0], [1.0, 0.0]], [1.0, 2.0]);", "lfY"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.loadConstructedPoint([[0.0,], [1.0,]], [[1.0,], [2.0,]]);", "lfX"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.loadConstructedPoint([[0.0, 0.0], [1.0, 0.0]], [[1.0, 2.0], [1.0, 2.0]]);", "lfY"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.loadConstructedPoint([[0.0, 0.0], [1.0, 0.0]], [[1.0,], [2.0,], [3.0,]]);", "lfY"],
                   ["grid.makeLocalPolynomialGrid(2, 1, 1, 1, 'localp'); grid.beginConstruction(); grid.loadConstructedPoint([[[0.0, 0.0], [1.0, 0.0]],], [[1.0,], [2.0,]]);", "lfX"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.enableAcceleration('gpu-wrong');", "sAccelerationType"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.enableAcceleration('gpu-default');", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.enableAcceleration('gpu-default', 0 if grid.isAccelerationAvailable('gpu-cuda') else None);", "notError"],
                   ["grid.makeSequenceGrid(2, 1, 2, 'level', 'leja'); grid.enableAcceleration('gpu-default', -11);", "iGPUID"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.isAccelerationAvailable('cpu-wrong');", "sAccelerationType"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.isAccelerationAvailable('cpu-blas');", "notError"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.getGPUMemory(-1);", "iGPUID"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.getGPUMemory(1000000);", "iGPUID"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.getGPUMemory(grid1.getNumGPUs());", "iGPUID"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.getGPUName(-1);", "iGPUID"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.getGPUName(1000000);", "iGPUID"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.getGPUName(grid1.getNumGPUs());", "iGPUID"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.setGPUID(-1);", "iGPUID"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.setGPUID(1000000);", "iGPUID"],
                   ["grid1 = Tasmanian.SparseGrid(); grid1.setGPUID(grid1.getNumGPUs());", "iGPUID"],
                   ["grid.makeLocalPolynomialGrid(1, 1, 1, 1, 'localp'); Tasmanian.loadNeededPoints(lambda x, tid : x, grid, 1);", "notError"],
                   ["grid.makeLocalPolynomialGrid(1, 1, 1, 1, 'localp'); Tasmanian.loadNeededPoints(lambda x, tid : np.ones((2,)) * x, grid, 1);", "loadNeededValues"],
                   ]

    def getCustomTabulatedTests(self):
        # same format as in getSparseGridTests().
        return    [["Tasmanian.makeCustomTabulatedFromData(1, np.array([2]), np.array([1]), np.array([[0.0]]), np.array([[2.0]]), 'testCT')", "nodes[0]"], 
                   ["Tasmanian.makeCustomTabulatedFromData(1, np.array([2]), np.array([1]), np.array([[0.0, 1.0]]), np.array([[2.0]]), 'testCT')", "weights[0]"], 
                   ["Tasmanian.makeCustomTabulatedFromData(1, np.array([2, 3]), np.array([1]), np.array([[0.0, 1.0]]), np.array([[1.0, 1.0]]), 'testCT')", "num_nodes"], 
                   ["Tasmanian.makeCustomTabulatedFromData(1, np.array([2]), np.array([1, 2]), np.array([[0.0, 1.0]]), np.array([[1.0, 1.0]]), 'testCT')", "precision"], 
                  ]

    def getDreamTests(self):
        # see getSparseGridTests() for comments about the format
        return [["DREAM.Domain('incorrect')", "Domain"],
                ["DREAM.IndependentUpdate(1)", "sType"],
                ["DREAM.DifferentialUpdate(0.0)", "callableOrMagnitude"],
                ["DREAM.RandomGenerator('incorrect')", "sType"],
                ["DREAM.RandomGenerator(callableRNG = 0.0)", "callableRNG"],
                ["DREAM.Posterior(lambda x: x, lambda y: y, lambda x: 1, typeForm = 'incorrect')", "typeForm"],
                ["DREAM.Posterior('incorrect', lambda y: y, lambda x: 1)", "model"],
                ["DREAM.Posterior(lambda x: x, 'incorrect', lambda x: 1)", "likelihood"],
                ["DREAM.Posterior(lambda x: x, lambda y: y, 'prior')", "prior"],
                ["DREAM.genUniformSamples([0, 0], [1], 20)", "lfUpper"],
                ["DREAM.genGaussianSamples([0, 0], [1], 20)", "lfDeviation"],
                ["DREAM.State(3, 2).setState(np.array([[0.1, 0.1],[0.1, 0.2],[0.1, 0.3],]))", "notError"],
                ["DREAM.State(3, 2).setState(np.array([[0.1, 0.1],[0.1, 0.2],]))", "llfNewState"],
                ["DREAM.State(3, 2).setState(np.array([[0.1,],[0.1,],[0.1,],]))", "llfNewState"],
                ["""DREAM.Sample(5, 5, lambda x : np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ),
                    DREAM.Domain('unbounded'),
                    state,
                    DREAM.IndependentUpdate('gaussian', 3.0),
                    DREAM.DifferentialUpdate(0))""", "notError"],
                ["""DREAM.Sample(5, 5, lambda x : np.ones((2, 1)) * np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ),
                    DREAM.Domain('unbounded'),
                    state,
                    DREAM.IndependentUpdate('gaussian', 3.0),
                    DREAM.DifferentialUpdate(0))""", "Sample"],
                ["""DREAM.Sample(5, 5, lambda x : np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ),
                    DREAM.Domain('unbounded'),
                    state,
                    DREAM.IndependentUpdate(lambda : np.ones((3, 1))),
                    DREAM.DifferentialUpdate(0))""", "Sample"],
                ["""DREAM.Sample(5, 5, DREAM.Domain('unbounded'),
                    DREAM.Domain('unbounded'),
                    state,
                    DREAM.IndependentUpdate('gaussian', 3.0),
                    DREAM.DifferentialUpdate(0))""", "probability_distibution"],
                ["""DREAM.Sample(5, 5, lambda x : np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ),
                    (-1.0, 1.0),
                    state,
                    DREAM.IndependentUpdate('gaussian', 3.0),
                    DREAM.DifferentialUpdate(0))""", "domain_description"],
                ["""DREAM.Sample(5, 5, lambda x : np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ),
                    DREAM.Domain('unbounded'),
                    grid,
                    DREAM.IndependentUpdate('gaussian', 3.0),
                    DREAM.DifferentialUpdate(0))""", "dream_state"],
                ["""DREAM.Sample(5, 5, lambda x : np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ),
                    DREAM.Domain('unbounded'),
                    state,
                    DREAM.DifferentialUpdate(0),
                    DREAM.DifferentialUpdate(0))""", "independent_update"],
                ["""DREAM.Sample(5, 5, lambda x : np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ),
                    DREAM.Domain('unbounded'),
                    state,
                    DREAM.IndependentUpdate('gaussian', 3.0),
                    DREAM.IndependentUpdate('gaussian', 3.0))""", "differential_update"],
                ["""DREAM.Sample(5, 5, lambda x : np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ),
                    DREAM.Domain('unbounded'),
                    state,
                    DREAM.IndependentUpdate('gaussian', 3.0),
                    DREAM.DifferentialUpdate(0), lambda : 0.5)""", "random01"],
                ["""DREAM.Sample(5, 5, lambda x : np.exp( -0.5 * np.sum( (x - 2.0)**2, 1 ) / 9.0 ),
                    DREAM.Domain('unbounded'),
                    state,
                    DREAM.IndependentUpdate('gaussian', 3.0),
                    DREAM.DifferentialUpdate(0), typeForm = 'alpha')""", "typeForm"],
                ]

    def getOptimizationTests(self):
        # see getSparseGridTests() for comments about the format
        return [["pss.setParticlePositions(np.array([[1, 2, 3]]))", "llfNewPPosns"],
                ["pss.setParticlePositions(np.array([[1, 2], [2, 1]]))", "llfNewPPosns"],
                ["pss.setParticleVelocities(np.array([[1, 2, 3]]))", "llfNewPVelcs"],
                ["pss.setParticleVelocities(np.array([[1, 2], [2, 1]]))", "llfNewPVelcs"],
                ["pss.setBestParticlePositions(np.array([[1, 2, 3]]))", "llfNewBPPosns"],
                ["pss.initializeParticlesInsideBox(np.array([-1.0, -2.0, -3.0]), np.array([2.0, 1.0]))", "lfBoxLower"],
                ["pss.initializeParticlesInsideBox(np.array([-1.0, -2.0]), np.array([2.0, 1.0, 3.0]))", "lfBoxUpper"],
                ["pss.initializeParticlesInsideBox(np.array([-1.0, -2.0]), np.array([1.0, 3.0]));" +
                 "f = lambda x_batch : np.append(np.apply_along_axis(np.sum, 1, x_batch), np.array([0.5]), axis=0);" +
                 "inside = lambda x : True;" +
                 "Opt.ParticleSwarm(f, inside, 0.5, 2, 2, 1, pss)", "ParticleSwarm"],
                ["pss.initializeParticlesInsideBox(np.array([-1.0, -2.0]), np.array([1.0, 3.0]));" +
                 "f = lambda x_batch : np.apply_along_axis(np.sum, 1, x_batch);" +
                 "inside = lambda x : 0.01;" +
                 "Opt.ParticleSwarm(f, inside, 0.5, 2, 2, 1, pss)", "ParticleSwarm"],
                ]

    def testListedExceptions(self, llTests):
        grid = Tasmanian.SparseGrid()
        state = DREAM.State(10, 2)
        state.setState(DREAM.genGaussianSamples([-1.0, -1.0], [1.0, 1.0], 10, DREAM.RandomGenerator("minstd_rand", 42)))
        pss = Opt.ParticleSwarmState(2, 3)

        for lTest in llTests:
            try:
                exec(lTest[0])
                self.assertEqual(lTest[1], "notError", "failed to raise exception for invalid '{0:1s}' using test\n '{1:1s}'".format(lTest[1],lTest[0]))
            except Tasmanian.TasmanianInputError as TSGError:
                self.assertEqual(TSGError.sVariable, lTest[1], "error raising exception for '{0:1s}' using test\n '{1:1s}'\n Error.sVariable = '{2:1s}'".format(lTest[1],lTest[0],TSGError.sVariable))

    def performExceptionsTest(self):
        self.testListedExceptions(self.getSparseGridTests())
        self.testListedExceptions(self.getCustomTabulatedTests())
        self.testListedExceptions(self.getDreamTests())
        self.testListedExceptions(self.getOptimizationTests())
