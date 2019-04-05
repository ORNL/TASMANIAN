########################################################################
# simple pythons cript to check if numpy module exists
########################################################################

import numpy as np
import ctypes as ct
import sys

a = np.array([0.0, 1.0, 2.0])
b = (ct.c_int * 3)()
