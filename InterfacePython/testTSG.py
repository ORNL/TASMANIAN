#!/usr/bin/env python3

# The hash-bang is used only but the GNU Make build system
# The CMake build system uses the PYTHON_EXECUTABLE and ignores the hash-bang.

import unittest
import TasmanianSG
import math
import sys, os
import numpy as np
from random import uniform
from ctypes import cdll

import testCommon
import testBasicIO
import testAcceleration
import testExceptions
import testMakeUpdate
import testRefinement
import testUnstructuredData
import testMisc
import testAddons
import testDream

grid = TasmanianSG.TasmanianSparseGrid()

ttc = testCommon.TestTasCommon()

# python-coverage run testTSG.py
# python-coverage html
# python-coverage report

class TestTasmanian(unittest.TestCase):
    def testBasicIO(self):
        print("\nTesting core I/O test")
        tester = testBasicIO.TestTasClass()
        tester.performIOTest()

    def testAcceleratedEvaluate(self):
        print("\nTesting accelerated evaluate consistency")
        tester = testAcceleration.TestTasClass()
        tester.performAccelerationTest()

    def testBasicException(self):
        print("\nTesting error handling")
        tester = testExceptions.TestTasClass()
        tester.performExceptionsTest()

    def testAMakeUpdate(self):
        print("\nTesting core make/update grid")
        tester = testMakeUpdate.TestTasClass()
        tester.performMakeUpdateTests()

    def testBRefinement(self):
        print("\nTesting core refine grid")
        tester = testRefinement.TestTasClass()
        tester.performRefinementTest()

    def testCUnsructuredData(self):
        print("\nTesting core learning from random samples")
        tester = testUnstructuredData.TestTasClass()
        tester.performUnstructuredDataTests()

    def testZMisc(self):
        print("\nTesting plotting and other misc")
        tester = testMisc.TestTasClass()
        tester.performMiscTests()

    def testAddons(self):
        print("\nTesting addon wrappers")
        tester = testAddons.TestTasClass()
        tester.performAddonTests()

    def testDream(self):
        print("\nTesting DREAM wrappers")
        tester = testDream.TestTasClass()
        tester.performDreamTests()

if __name__ == '__main__':
    unittest.main()


