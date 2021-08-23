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

// TODO: Add better documentation.

#ifndef __TASMANIAN_ADDONS_EXOTICQUADRATURE_HPP
#define __TASMANIAN_ADDONS_EXOTICQUADRATURE_HPP

//! \internal
//! \brief Contains algorithms for generating exotic quadrature rules.
//! \endinternal

// Generate the n roots of the n-th degree orthogonal polynomial.
// ref_integral_weights and ref_points are used in the computation of the
// integrals that form the elements of the Jacobi matrix.

// Get Exotic Gauss-Legendre points and weights. Specific to the case where the
// integrand F(x) is of the form:
//     F(x) := f(x) * [weight_fn(x) - shift] + shift * f(x)

#include "tsgAddonsCommon.hpp"

namespace TasGrid {

//! \brief Generate all roots up to the n roots of the n-th degree orthogonal polynomial.
std::vector<std::vector<double>> getRoots(
      const int n,
      const std::vector<double> &ref_integral_weights,
      const std::vector<double> &ref_points);

//! The integrand is assumed to be of the form F(x) = f(x) * weight_fn(x) and
//! the function x->weight_fn(x)+shift is nonnegative.
inline TasGrid::CustomTabulated getExoticQuadrature(const int n,
                                                    const double shift,
                                                    std::function<double(double)> weight_fn,
                                                    const int nref,
                                                    const char* description) {

    // Create the set of reference weights and points for the first term. These
    // reference weights are with respect to the measure induced by the
    // function [weight_fn(x) + shift].
    std::vector<double> ref_weights(nref), ref_points(nref);
    TasGrid::OneDimensionalNodes::getGaussLegendre(nref, ref_weights, ref_points);
    for (size_t k=0; k<nref; k++) {
        ref_weights[k] *= (weight_fn(ref_points[k]) + shift);
    }
    std::vector<std::vector<double>> points_cache = getRoots(n, ref_weights, ref_points);

    // Create the set of weights for the first term.
    std::vector<std::vector<double>> weights_cache(n);
    for (int i=0; i<n; i++) {
        weights_cache[i] = std::vector<double>(points_cache[i].size(), 0.0);
        for (size_t j=0; j<points_cache[i].size(); j++) {
            // Integrate the function l_j(x) according to the measure induced
            // by the function [weight_fn(x) + shift].
            for (int k=0; k<nref; k++) {
                weights_cache[i][j] +=
                        lagrange_eval(j, points_cache[i], ref_points[k]) *
                        ref_weights[k];
            }
        }
    }

    // Create and append the set of points and weights for the second term if
    // the shift is nonzero.
    if (shift != 0.0) {
        for (int i=0; i<n; i++) {
            std::vector<double>
                    correction_points(points_cache[i].size()),
                    correction_weights(points_cache[i].size());
            TasGrid::OneDimensionalNodes::getGaussLegendre(points_cache[i].size(),
                                                           correction_weights,
                                                           correction_points);
            for (auto &w : correction_weights) w *= -shift;
            if (correction_points.size() % 2 == 1) {
                // Zero out for stability.
                correction_points[(correction_points.size() - 1) / 2] = 0.0;
            }
            points_cache[i].insert(points_cache[i].end(),
                                   correction_points.begin(),
                                   correction_points.end());
            weights_cache[i].insert(weights_cache[i].end(),
                                    correction_weights.begin(),
                                    correction_weights.end());
        }
    }

    // Output to a CustomTabulated instance.
    std::vector<int> num_nodes(n), precision(n);
    for (int i=0; i<n; i++) {
        num_nodes[i] = points_cache[i].size();
        precision[i] = num_nodes[i] - 1;
    }
    return TasGrid::CustomTabulated(
        n, std::move(num_nodes), std::move(precision), std::move(points_cache),
        std::move(weights_cache), description);
}

} // namespace(TasGrid)

void getExoticGaussLegendreCache(
        const int n,
        const double shift,
        std::function<double(double)> weight_fn,
        const int nref,
        std::vector<std::vector<double>> &weights_cache,
        std::vector<std::vector<double>> &points_cache);
}
#endif
