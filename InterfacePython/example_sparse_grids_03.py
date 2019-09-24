##############################################################################################################################################################################
# Copyright (c) 2017, Miroslav Stoyanov
#
# This file is part of
# Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
#    and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
# THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
# COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
# THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
# IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
##############################################################################################################################################################################

import numpy as np
import Tasmanian

def example_03():

    print("\n---------------------------------------------------------------------------------------------------\n")
    print("Example 3: integrate exp(-x1^2 - x2^2) * cos(x3) * cos(x4)")
    print("           for x1, x2 in [-5,5]; x3, x4 in [-2,3]")
    print("           using different rules and total degree polynomial space\n")

    def make_grid(iPrecision, sRule):
        grid = Tasmanian.makeGlobalGrid(4, 0, iPrecision, "qptotal", sRule)
        grid.setDomainTransform(np.array([[-5.0, 5.0], [-5.0, 5.0], [-2.0, 3.0], [-2.0, 3.0]]))
        return grid

    def print_error(grid):
        fExactIntegral = 1.861816427518323e+00 * 1.861816427518323e+00
        aPoints = grid.getPoints()
        aWeights = grid.getQuadratureWeights()

        fApproximateIntegral = np.sum(aWeights * np.exp(-aPoints[:,0]**2 -aPoints[:,1]**2)
                                      * np.cos(aPoints[:,2]) * np.cos(aPoints[:,3]))
        fError = np.abs(fApproximateIntegral - fExactIntegral)
        return "{0:>10d}{1:>10.2e}".format(grid.getNumPoints(), fError)

    print("               Clenshaw-Curtis      Gauss-Legendre    Gauss-Patterson")
    print(" precision    points     error    points     error    points    error")

    for prec in range(5, 41, 5):
        print("{0:>10d}{1:1s}{2:1s}{3:1s}".format(
            prec,
            print_error(make_grid(prec, "clenshaw-curtis")),
            print_error(make_grid(prec, "gauss-legendre-odd")),
            print_error(make_grid(prec, "gauss-patterson"))))

    print("\nAt 311K points the Gauss-Legendre error is O(1.E-1),")
    print("                   Clenshaw-Curtis error is O(1.E-7) at 320K points.")
    print("At 70K points the Gauss-Patterson error is O(1.E-4),")
    print("                  Clenshaw-Curtis needs 158K points to achieve the same.")


if (__name__ == "__main__"):
    example_03()
