/*
 * Copyright (c) 2021, Miroslav Stoyanov & William Kong
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
 * conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
 *    and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND
 * IMPLIED. THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE
 * OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL
 * ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE. THE USER ASSUMES
 * RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING
 * FROM OR ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */


#include "TasmanianAddons.hpp"
#include "tasgridCLICommon.hpp"
#include <cstring>

// Test the output of getRoots().
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
            std::cout << "ERROR: Test failed in testRootSizes() for n = " << n << std::endl;
            passed = false;
        }
    }
    return passed;
}

// Test the attributes of the output of getExoticQuadrature().
inline bool testBasicAttributes() {
    bool passed = true;
    const int n = 3;
    const int nref = 11;
    const char *descr = "TEST_DESCRIPTION";
    auto sinc = [](double x) -> double {return (x == 0.0 ? 1.0 : sin(x) / (x));};
    // Different behaviors expected for shift == 0.0 vs shift != 0.0.
    std::vector<double> shifts = {0.0, 1.0};
    for (auto shift : shifts) {
        TasGrid::CustomTabulated ct = TasGrid::getExoticQuadrature(n, shift, sinc, nref, descr);
        if (strcmp(ct.getDescription(), descr) != 0) {
            std::cout << "ERROR: Test failed in testCtAttributes() for shift = " << shift << " on getDescription()"<< std::endl;
            passed = false;
        }
        if (ct.getNumLevels() != n) {
            std::cout << "ERROR: Test failed in testCtAttributes() for shift = " << shift << " on getNumLevels()" << std::endl;
            passed = false;
        }
        for (int i=0; i<n; i++) {
            if (ct.getQExact(i) != 2*(i+1)-1) {
                std::cout << "ERROR: Test failed in testCtAttributes() for shift = " << shift << ", i = " << i
                          << " on getQExact()" << std::endl;
                passed = false;
            }
        }
    }
    return passed;
}


// Given a dimension, integrates the function f(x[1],...,x[dimension]) * sinc(freq * x[1]) * ... * sinc(freq * x[dimension]) over
// [-1, 1]^dimension  using Exotic quadrature at a given depth and a specified shift. Then, compares this integral against a given exact
// integral value.
inline bool wrapSincTest(std::function<double(std::vector<double>)> f, int dimension, int depth, double freq, double shift,
                         double exact_integral, double precision = 1e-10) {

    // Initialize.
    bool passed = true;
    auto sinc = [freq](double x)->double {return (x == 0.0 ? 1.0 : sin(freq * x) / (freq * x));};
    int level = depth % 2 == 1 ? (depth + 1) / 2 : depth / 2 + 1;
    TasGrid::CustomTabulated ct = TasGrid::getExoticQuadrature(level, shift, sinc);
    TasGrid::TasmanianSparseGrid sg;
    sg.makeGlobalGrid(dimension, 1, depth, TasGrid::type_qptotal, std::move(ct));

    // Compute the integral and compare to the reference value.
    std::vector<double> quad_points = sg.getPoints();
    std::vector<double> quad_weights = sg.getQuadratureWeights();
    double approx_integral = 0.0;
    assert(quad_weights.size() * dimension == quad_points.size());
    for (size_t i=0; i<quad_weights.size(); i++) {
        std::vector<double> x(dimension);
        for (int d=0; d<dimension; d++) {
            x[d] = quad_points[dimension * i + d];
        }
        approx_integral += f(x) * quad_weights[i];
    }
    if (std::abs(approx_integral - exact_integral) > precision) {
        std::cout << "ERROR: " << std::setprecision(16) << "Computed integral value " << approx_integral
                  << " does not match exact integral " << exact_integral << " for test problem " << std::endl
                  << "∫_[-1,1]^n f(x[1],...,x[n]) * sinc(freq*x[n]) * ... * sinc(freq*x[n]) dx[1] ... dx[n]" << std::endl
                  << "with inputs: " << std::endl << std::endl << std::left
                  << std::setw(10) << "dimension" << " = " << dimension << std::endl
                  << std::setw(10) << "depth"     << " = " << depth     << std::endl
                  << std::setw(10) << "freq"      << " = " << freq      << std::endl
                  << std::setw(10) << "shift"     << " = " << shift     << std::endl
                  << std::setw(10) << "precision" << " = " << precision << std::endl << std::endl;
        passed = false;
    }
    return passed;

}

// Tests for the accuracy of getExoticQuadrature() when weight_fn(x) is sinc(freq * x) and f(x) is exp(-x'*x).
inline bool testAccuracy() {
    bool passed = true;
    int depth, dimension;
    // 1D problem instances.
    depth = 20;
    dimension = 1;
    auto f1 = [](std::vector<double> x)->double{return std::exp(-x[0]*x[0]);};
    passed = passed && wrapSincTest(f1, dimension, depth, 1.0, 0.0, 1.4321357541271255);
    passed = passed && wrapSincTest(f1, dimension, depth, 1.0, 1.0, 1.4321357541271255);
    passed = passed && wrapSincTest(f1, dimension, depth, 10.0, 1.0, 0.32099682841103033);
    passed = passed && wrapSincTest(f1, dimension, depth, 100.0, 1.0, 0.031353648322695503);
    // 2D problem instances.
    depth = 40;
    dimension = 2;
    auto f2 = [](std::vector<double> x)->double{return std::exp(-x[0]*x[0]-x[1]*x[1]);};
    passed = passed && wrapSincTest(f2, dimension, depth, 1.0, 0.0, 2.051012818249270);
    passed = passed && wrapSincTest(f2, dimension, depth, 1.0, 1.0, 2.051012818249270);
    passed = passed && wrapSincTest(f2, dimension, depth, 10.0, 1.0, 0.1030389638499404);
    passed = passed && wrapSincTest(f2, dimension, depth, 100.0, 1.0, 0.0009830512631432665);
    return passed;
}

// Main wrapper
inline bool testExoticQuadrature() {
    bool passed = true;
    passed = passed && testRootSizes();
    passed = passed && testBasicAttributes();
    passed = passed && testAccuracy();
    return passed;
}

