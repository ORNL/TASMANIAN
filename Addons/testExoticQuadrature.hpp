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

inline void debugGetRoots() {
    // Test 1 (Print Roots)
    int n = 10;
    int nref = 101;
    std::vector<double> ref_weights(nref), ref_points(nref);
    TasGrid::OneDimensionalNodes::getGaussLegendre(nref, ref_weights, ref_points);
    std::vector<double> ref_integral_weights(ref_weights.size());
    std::transform(ref_weights.begin(),
                   ref_weights.end(),
                   ref_integral_weights.begin(),
                   [](double x)->double{
                       return (x == 0.0 ? 2.0 : 1.0 + sin(x) / x);});
    std::vector<std::vector<double>> roots =
            TasGrid::getRoots(n, ref_integral_weights, ref_points);
    for (int j=0; j<roots.size(); j++) {
        std::cout << "n = " << j + 1 << std::endl;
        for (int k=0; k<roots[j].size(); k++) {
            std::cout << roots[j][k] << std::endl;
        }
        std::cout << std::endl;
    }
}

inline void debugGetExoticQuadrature() {
    int n = 7;
    int nref = 201;
    double shift = 1.0;
    auto sinc = [](double x)->double{return(x == 0.0 ? 1.0 : sin(x) / x);};
    TasGrid::CustomTabulated ct = TasGrid::getExoticQuadrature(n, shift, sinc, nref, "tmp");
    std::vector<double> points;
    std::vector<double> weights;
    for (int j=0; j<n; j++) {
        std::cout << "POINTS, n = " << j + 1 << std::endl;
        ct.getWeightsNodes(j, weights, points);
        for (size_t k=0; k<points.size(); k++) {
            std::cout << points[k] << std::endl;
        }
        std::cout << "WEIGHTS, n = " << j + 1 << std::endl;
        for (size_t k=0; k<weights.size(); k++) {
            std::cout << weights[k] << std::endl;
        }
        std::cout << std::endl;
    }
}

inline void debugSincT(float T, double exact) {

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
}