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

#include "tsgCoreOneDimensional.hpp"
#include "tsgLinearSolvers.hpp"
#include <cassert>

namespace TasGrid {

// Evaluates a polynomial with roots given by \b roots at the point \b x.
inline double poly_eval(const std::vector<double> &roots, double x) {
    double eval_value = 1.0;
    for (auto r : roots) eval_value *= (x - r);
    return eval_value;
}
// Generate roots.
inline std::vector<std::vector<double>> getRoots(
        const int n,
        const std::vector<double> &ref_weights,
        const std::vector<double> &ref_points) {

    // Compute the roots incrementally.
    assert(ref_points.size() == ref_weights.size());
    std::vector<double>
            poly_m1_vals(ref_points.size(), 0.0),
            poly_vals(ref_points.size(), 1.0);
    std::vector<double> diag, offdiag;
    std::vector<std::vector<double>> roots(n);
    for (int i=0; i<n; i++) {
        // Form the tridiagonal vectors, alpha and beta, of the Jacobi matrix.
        double diag_numr = 0.0;
        double diag_denm = 0.0;
        double sqr_offdiag_numr = 0.0;
        double sqr_offdiag_denm = 0.0;
        for (int j=0; j<ref_points.size(); j++) {
            diag_numr += ref_points[j] * (poly_vals[j] * poly_vals[j]) * ref_weights[j];
            diag_denm += (poly_vals[j] * poly_vals[j]) * ref_weights[j];
            sqr_offdiag_numr += (poly_vals[j] * poly_vals[j]) * ref_weights[j];
            sqr_offdiag_denm +=  (poly_m1_vals[j] * poly_m1_vals[j])  * ref_weights[j];
        }
        diag.push_back(diag_numr / diag_denm);
        if (i >= 1) {
            offdiag.push_back(std::sqrt(sqr_offdiag_numr / sqr_offdiag_denm));
        }
        // Compute the roots. The inputs need to be copied because the
        // tridiagonal eigensolver generates the eigenvalues in place and
        // destroys some of its inputs.
        roots[i] = diag;
        std::vector<double> offdiag_dummy = offdiag;
        TasmanianTridiagonalSolver::getSymmetricEigenvalues(i+1, roots[i], offdiag_dummy);
        if (roots[i].size() % 2 == 1) {
            // Zero out the center for stability.
            roots[i][(roots[i].size() - 1) / 2] = 0.0;
        }
        // Update the values of the polynomials at ref_points.
        std::swap(poly_m1_vals, poly_vals);
        std::transform(ref_points.begin(), ref_points.end(), poly_vals.begin(),
                       [&roots, i](double x){return poly_eval(roots[i], x);});
    }
    return roots;
}

//! \brief Evaluates a Lagrange basis polynomial at a point \b x.
inline double lagrange_eval(int idx, const std::vector<double> &roots, double x) {
    double eval_value = 1.0;
    for (int i=0; i<idx; i++) {
        eval_value *= (x - roots[i]) / (roots[idx] - roots[i]);
    }
    for (size_t i=idx+1; i<roots.size(); i++) {
        eval_value *= (x - roots[i]) / (roots[idx] - roots[i]);
    }
    return eval_value;
}
//! \brief Generate the full cache of Exotic (weighted and possibly shifted)
//! Gauss-Legendre weights and points for an (input) max level \b n, using
//! parameters \b shift, \b weight_fn, and \b nref. By default, the reference
//! points and reference weights are generated by the Gauss-Legendre rule.
inline void getExoticGaussLegendreCache(
        const int n,
        const double shift,
        std::function<double(double)> weight_fn,
        const int nref,
        std::vector<std::vector<double>> &weights_cache,
        std::vector<std::vector<double>> &points_cache) {

    // Create the set of reference weights and points for the first term. These
    // reference weights are with respect to the measure induced by
    // [weight_fn() + shift].
    std::vector<double> ref_weights(nref), ref_points(nref);
    TasGrid::OneDimensionalNodes::getGaussLegendre(nref, ref_weights, ref_points);
    for (size_t k=0; k<nref; k++) {
        ref_weights[k] *= (weight_fn(ref_points[k]) + shift);
    }
    points_cache = getRoots(n, ref_weights, ref_points);

    // Create the set of weights for the first term.
    weights_cache = std::vector<std::vector<double>>(n);
    for (int i=0; i<n; i++) {
        weights_cache[i] = std::vector<double>(points_cache[i].size(), 0.0);
        for (size_t j=0; j<points_cache[i].size(); j++) {
            // Integrate the function l_j(x) relative to the measure given by
            // ref_weights().
            for (int k=0; k<nref; k++) {
                weights_cache[i][j] +=
                        lagrange_eval(j, points_cache[i], ref_points[k]) *
                        ref_weights[k];
            }
        }
    }

    // Create and append the set of points and weights for the second term.
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

} // namespace(TasGrid)
#endif
