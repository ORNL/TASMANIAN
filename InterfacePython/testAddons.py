import unittest
import Tasmanian
import os
import numpy as np

import testCommon

ttc = testCommon.TestTasCommon()

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
        pass

    def performAddonTests(self):
        self.checkLoadNeeded()
