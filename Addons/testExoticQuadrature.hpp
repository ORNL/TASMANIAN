/*
 * Copyright (c) 2021, Miroslav Stoyanov & William Kong
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN:
 * TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND
 * DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED. THERE ARE NO EXPRESS OR
 * IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR
 * THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT,
 * TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH
 * THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN
 * INJURY OR DAMAGE. THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES,
 * PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED
 * BY, RESULTING FROM OR ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR
 * DISPOSAL OF THE SOFTWARE.
 */

// TODO: Add full unit tests and documentation.

#include "TasmanianAddons.hpp"
#include "tasgridCLICommon.hpp"
#include <cstring>

// Tests for getRoots().
inline bool testRootSizes() {
    bool passed = true;
    std::vector<std::vector<double>> roots;
    int nref = 5;
    std::vector<double> ref_weights(nref), ref_points(nref);
    TasGrid::OneDimensionalNodes::getGaussLegendre(nref, ref_weights, ref_points);
    // n = 0, 1, 5
    std::vector<int> n_vec = {1, 2, 5};
    for (auto n : n_vec) {
        roots = TasGrid::getRoots(0, ref_weights, ref_points);
        if (roots.size() != 0) {
            std::cout << "ERROR: Test failed in testRootSizes() for n = "
                      << n << std::endl;
            passed = false;
        }
    }
    return passed;
}

// Tests for getExoticQuadrature().
inline bool testBasicAttributes() {
    bool passed = true;
    const int n = 3;
    const int nref = 11;
    const char *descr = "TEST_DESCRIPTION";
    auto sinc = [](double x) -> double {
      return (x == 0.0 ? 1.0 : sin(x) / (x));
    };
    // Different behaviors expected for shift == 0.0 vs shift != 0.0.
    std::vector<double> shifts = {0.0, 1.0};
    for (auto shift : shifts) {
        TasGrid::CustomTabulated ct =
                TasGrid::getExoticQuadrature(n, shift, sinc, nref, descr);
        if (strcmp(ct.getDescription(), descr) != 0) {
            std::cout << "ERROR: Test failed in testCtAttributes() for shift = "
                      << shift << " on getDescription()"<< std::endl;
            passed = false;
        }
        if (ct.getNumLevels() != n) {
            std::cout << "ERROR: Test failed in testCtAttributes() for shift = "
                      << shift << " on getNumLevels()" << std::endl;
            passed = false;
        }
        for (int i=0; i<n; i++) {
            if (shift == 0.0 && ct.getNumPoints(i) != i+1) {
                std::cout << "ERROR: Test failed in testCtAttributes() for shift = "
                          << shift << ", i = " << i
                          << " on getNumPoints()" << std::endl;
                passed = false;
            }
            if (shift != 0.0 && ct.getNumPoints(i) != 2*(i+1)) {
                std::cout << "ERROR: Test failed in testCtAttributes() for shift = "
                          << shift << ", i = " << i
                          << " on getNumPoints()" << std::endl;
                passed = false;
            }
            if (ct.getQExact(i) != 2*(i+1)-1) {
                std::cout << "ERROR: Test failed in testCtAttributes() for shift = "
                          << shift << ", i = " << i
                          << " on getQExact()" << std::endl;
                passed = false;
            }
        }
    }
    return passed;
}

// Integrates the function f(x) * sinc(freq * x) over [-1, 1]  using Exotic
// quadrature at level n and a specified shift.
inline double integrateFnTimesSinc(std::function<double(double)> f,
                                   int n,
                                   double freq,
                                   double shift) {
  auto sinc = [freq](double x)->double {
    return (x == 0.0 ? 1.0 : sin(freq * x) / (freq * x));
  };
  TasGrid::CustomTabulated ct = TasGrid::getExoticQuadrature(n, shift, sinc);
  TasGrid::TasmanianSparseGrid sg;
  sg.makeGlobalGrid(1, 1, n, TasGrid::type_qptotal, std::move(ct));
  std::vector<double> quad_points = sg.getPoints();
  std::vector<double> quad_weights = sg.getQuadratureWeights();
  double integral = 0.0;
  for (size_t i=0; i<quad_weights.size(); i++) {
      integral += f(quad_points[i]) * quad_weights[i];
  }
  return integral;
}

// Test the accuracy of 1D problem instances.
inline bool wrapSincTest1D(std::function<double(double)> f,
                           int level,
                           double freq,
                           double shift,
                           double exact_integral) {
    bool passed = true;
    double precision = 1e-12;
    double approx_integral = integrateFnTimesSinc(f, level, 1.0, 0.0);
    if (std::abs(approx_integral - exact_integral) > precision) {
        std::cout << "ERROR: " << std::setprecision(16)
                  << "Computed integral value " << approx_integral
                  << " does not match exact integral " << exact_integral
                  << std::endl;
        passed = false;
    }
    return passed;
}
inline bool testExpMx2_sinc1_shift0() {
    auto f = [](double x)->double {return std::exp(-x*x);};
    bool passed = wrapSincTest1D(f, 20, 1.0, 0.0, 1.4321357541271255);
    if (not passed) {
        std::cout << "ERROR: Test failed in testExpMx_sinc1_shift0().";
    }
    return passed;
}
inline bool testExpMx2_sinc1_shift1() {
    auto f = [](double x)->double {return std::exp(-x*x);};
    bool passed = wrapSincTest1D(f, 20, 1.0, 1.0, 1.4321357541271255);
    if (not passed) {
        std::cout << "ERROR: Test failed in testExpMx_sinc1_shift0().";
    }
    return passed;
}

// inline void debugSincT(float T, double exact) {
    // int max_n = 15;
    // int nref = 1001;
    // double shift = 1.0;
    // auto f = [](double x)->double{return exp(-x * x);};
    // auto sinc = [T](double x)->double{return(x == 0.0 ? 1.0 : sin(T * x) / (T * x));};

    // std::vector<std::vector<double>> points_cache(max_n), weights_cache(max_n);
    // TasGrid::getExoticGaussLegendreCache(max_n, shift, sinc, nref, weights_cache, points_cache);

    // std::cout << std::fixed << std::setprecision(6) << "Exotic     \t GL" << std::endl;
    // for (int n=1; n<max_n; n++) {
    //     // Exotic Quadrature
    //     double approx_exotic = 0.0;
    //     for (size_t i=0; i<weights_cache[n-1].size(); i++) {
    //         approx_exotic += f(points_cache[n-1][i]) * weights_cache[n-1][i];
    //     }
    //     std::cout << log10(std::abs(approx_exotic - exact));
    //     std::cout << "\t";
    //     // Gauss-Legendre
    //     std::vector<double> gl_weights(n), gl_points(n);
    //     TasGrid::OneDimensionalNodes::getGaussLegendre(n, gl_weights, gl_points);
    //     double approx_GL = 0.0;
    //     for (size_t i=0; i<gl_weights.size(); i++) {
    //         approx_GL +=
    //                 f(gl_points[i]) * sinc(gl_points[i]) * gl_weights[i];
    //     }
    //     std::cout << log10(std::abs(approx_GL - exact));
    //     std::cout << std::endl;
    // }
    // std::cout << endl;
// }

