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

#ifndef __TASMANIAN_ADDONS_EXOTICQUADRATURE_HPP
#define __TASMANIAN_ADDONS_EXOTICQUADRATURE_HPP

#include "tsgAddonsCommon.hpp"

namespace TasGrid {

// Evaluates a polynomial with roots given by \b roots at the point \b x.
inline double poly_eval(const std::vector<double> &roots, double x) {
    double eval_value = 1.0;
    for (auto r : roots) eval_value *= (x - r);
    return eval_value;
}

//! \brief Generate the roots of all orthogonal polynomials up to degree n, where orthogonality is with respect to the
// measure given by ref_weights. Specializes depending on if the underlying weight function is symmetric or not.
template <bool is_symmetric>
inline std::vector<std::vector<double>> getRoots(const int n, const std::vector<double> &ref_weights,
                                                 const std::vector<double> &ref_points) {

    // Compute the roots incrementally.
    assert(ref_points.size() == ref_weights.size());
    assert(ref_points.size() > 0);
    std::vector<double> poly_m1_vals(ref_points.size(), 0.0), poly_vals(ref_points.size(), 1.0);
    std::vector<double> diag, offdiag;
    std::vector<std::vector<double>> roots(n);
    for (int i=0; i<n; i++) {
        // Form the diagonal and offdiagonal vectors of the Jacobi matrix.
        if (is_symmetric) {
            diag.push_back(0.0);
        } else {
            double diag_numr = 0.0;
            double diag_denm = 0.0;
            for (size_t j=0; j<ref_points.size(); j++) {
            diag_numr += ref_points[j] * (poly_vals[j] * poly_vals[j]) * ref_weights[j];
            diag_denm += (poly_vals[j] * poly_vals[j]) * ref_weights[j];
            }
            diag.push_back(diag_numr / diag_denm);
        }
        if (i >= 1) {
            double sqr_offdiag_numr = 0.0;
            double sqr_offdiag_denm = 0.0;
            for (size_t j=0; j<ref_points.size(); j++) {
                sqr_offdiag_numr += (poly_vals[j] * poly_vals[j]) * ref_weights[j];
                sqr_offdiag_denm += (poly_m1_vals[j] * poly_m1_vals[j])  * ref_weights[j];
            }
            offdiag.push_back(std::sqrt(sqr_offdiag_numr / sqr_offdiag_denm));
        }
        roots[i] = TasmanianTridiagonalSolver::getSymmetricEigenvalues(i+1, diag, offdiag);
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
inline double lagrange_eval(size_t idx, const std::vector<double> &roots, double x) {
    double eval_value = 1.0;
    for (size_t i=0; i<idx; i++) {
        eval_value *= (x - roots[i]) / (roots[idx] - roots[i]);
    }
    for (size_t i=idx+1; i<roots.size(); i++) {
        eval_value *= (x - roots[i]) / (roots[idx] - roots[i]);
    }
    return eval_value;
}

//! \brief Generate the full cache of Exotic (weighted and possibly shifted) quadrature weights and points for an (input) max level
//! \b n, using parameters \b shift, \b weight_fn_vals, \b shifted_weights, \b ref_points, and \b nref.
//!
//! The integral  I := ∫ F(x) * weight_fn(x) dx, for a given level, is approximated as:
//! I ≈ sum{f(roots[i]) * shifted_weights[i]} - shift * sum{f(correction_points[i]) * correction_weights[i]}.
//!
//! It is assumed that shifted_weights[i] is nonnegative and it is computed as shifted_weights[i] = ref_weights[i] *
//! (weight_fn_vals[i] + shift) for every i, where ref_weights are the quadrature weights corresponding to ref_points.
inline TasGrid::CustomTabulated getExoticQuadrature(const int n, const double shift, const std::vector<double> &weight_fn_vals,
                                                    const std::vector<double> &shifted_weights, const std::vector<double> &ref_points,
                                                    const int nref, const char* description, const bool is_symmetric = false) {

    // Create the set of points (roots) for the first term.
    assert(shifted_weights.size() == ref_points.size());
    assert(weight_fn_vals.size() == ref_points.size());
    std::vector<std::vector<double>> points_cache;
    if (is_symmetric) {
        points_cache = getRoots<true>(n, shifted_weights, ref_points);
    } else {
        points_cache = getRoots<false>(n, shifted_weights, ref_points);
    }

    // Create the set of weights for the first term.
    std::vector<std::vector<double>> weights_cache(n);
    for (int i=0; i<n; i++) {
        weights_cache[i] = std::vector<double>(points_cache[i].size(), 0.0);
        for (size_t j=0; j<points_cache[i].size(); j++) {
            // Integrate the function l_j(x) according to the measure induced by the function [weight_fn(x) + shift].
            for (int k=0; k<nref; k++) {
                weights_cache[i][j] += lagrange_eval(j, points_cache[i], ref_points[k]) * shifted_weights[k];
            }
        }
    }

    // Create and append the set of correction points and weights for the second term only if the shift is nonzero.
    if (shift != 0.0) {
        for (int i=0; i<n; i++) {
            const int init_size = points_cache[i].size();
            std::vector<double> correction_points, correction_weights;
            TasGrid::OneDimensionalNodes::getGaussLegendre(init_size, correction_weights, correction_points);
            for (auto &w : correction_weights) w *= -shift;
            if (correction_points.size() % 2 == 1) {
                // Zero out for stability.
                correction_points[(correction_points.size() - 1) / 2] = 0.0;
            }
            // Account for duplicates up to a certain tolerance.
            for (size_t j=0; j<correction_points.size(); j++) {
                int nonunique_idx = -1;
                for (int k=0; k<init_size; k++) {
                    if (std::abs(correction_points[j] - points_cache[i][k]) <= Maths::num_tol) {
                        nonunique_idx = k;
                        break;
                    }
                }
                if (nonunique_idx == -1) {
                    points_cache[i].push_back(correction_points[j]);
                    weights_cache[i].push_back(correction_weights[j]);
                } else {
                    weights_cache[i][nonunique_idx] += correction_weights[j];
                }
            }
        }
    }

    // Output to a CustomTabulated instance.
    std::vector<int> num_nodes(n), precision(n);
    for (int i=0; i<n; i++) {
        num_nodes[i] = points_cache[i].size();
        precision[i] = 2 * (i + 1) - 1;
    }
    return TasGrid::CustomTabulated(n, std::move(num_nodes), std::move(precision), std::move(points_cache),
                                    std::move(weights_cache), description);
}

// Basic function overloads.

// If only the weight function is provided, the reference points and weights are generated by the GL rule and the weight_fn values
// are generated by applying the weight function to these reference points.
inline TasGrid::CustomTabulated getExoticQuadrature(const int n, const double shift, std::function<double(double)> weight_fn,
                                                    const int nref, const char* description, const bool is_symmetric = false) {
    std::vector<double> shifted_weights, ref_points, weight_fn_vals(nref);
    TasGrid::OneDimensionalNodes::getGaussLegendre(nref, shifted_weights, ref_points);
    std::transform(ref_points.begin(), ref_points.end(), weight_fn_vals.begin(), weight_fn);
    if (not std::all_of(weight_fn_vals.begin(), weight_fn_vals.end(), [shift](double x){return (x+shift)>=0.0;})) {
        throw std::invalid_argument("shift needs to be chosen so that weight_fn(x)+shift is nonnegative!");
    }
    for (int k=0; k<nref; k++) {
        shifted_weights[k] *= (weight_fn_vals[k] + shift);
    }
    return getExoticQuadrature(n, shift, weight_fn_vals, shifted_weights, ref_points, nref, description, is_symmetric);
}

// If only the weight function is provided, the reference points are generated by the GL rule.

inline TasGrid::CustomTabulated getExoticQuadrature(const int n, const double shift, std::function<double(double)> weight_fn,
                                                    const bool is_symmetric = true) {
    return getExoticQuadrature(n, shift, weight_fn, 50*n+1, "Exotic quadrature", is_symmetric);
}

} // namespace(TasGrid)
#endif