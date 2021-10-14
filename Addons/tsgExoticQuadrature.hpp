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

// Specifies which algorithm is used in computing the Gauss quadrature points and weights in getGaussNodesAndWeights() (developer use only).
static int gauss_quadrature_version = 1;

/*!
 * \internal
 * \ingroup TasmanianAddonsExoticQuad
 * \brief Evaluates a polynomial with roots given by \b roots at the point \b x.
 * \endinternal
 */
inline double poly_eval(const std::vector<double> &roots, double x) {
    double eval_value = 1.0;
    for (auto r : roots) eval_value *= (x - r);
    return eval_value;
}

/*!
 * \internal
 * \ingroup TasmanianAddonsExoticQuad
 * \brief Evaluates a Lagrange basis polynomial at a point \b x.
 * \endinternal
 */
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

/*!
 * \internal
 * \ingroup TasmanianAddonsExoticQuad
 * \brief Generate the points and weights of a Gaussian quadrature (up to degree n) whose weight function is evaluated on \b
 * ref_points and whose corresponding values are stored in \b ref_weights. Specializes depending on if the underlying weight
 * function is symmetric or not (see \b is_symmetric).
 * \endinternal
 */
template <bool is_symmetric>
inline void getGaussNodesAndWeights(const int n, const std::vector<double> &ref_points, const std::vector<double> &ref_weights,
                                    std::vector<std::vector<double>> &points_cache, std::vector<std::vector<double>> &weights_cache) {
    #ifndef Tasmanian_ENABLE_BLAS
    // If BLAS is not enabled, we must use an implementation that does not require it.
    if (gauss_quadrature_version == 1) {
        gauss_quadrature_version = 2;
    }
    #endif
    // Compute the roots incrementally.
    assert(ref_points.size() == ref_weights.size());
    assert(ref_points.size() > 0);
    double mu0 = 0.0;
    for (auto &w : ref_weights) mu0 += w;
    std::vector<double> poly_m1_vals(ref_points.size(), 0.0), poly_vals(ref_points.size(), 1.0);
    std::vector<double> diag, offdiag;
    points_cache.resize(n);
    weights_cache.resize(n);
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
        // Compute the Gauss points and weights.
        if (gauss_quadrature_version == 1) {
            // Use LAPACK to obtain the Gauss points, and use the reference points and weights to obtain the Gauss weights.
            points_cache[i] = TasmanianTridiagonalSolver::getSymmetricEigenvalues(i+1, diag, offdiag);
            if (points_cache[i].size() % 2 == 1 && is_symmetric) {
                points_cache[i][(points_cache[i].size() - 1) / 2] = 0.0;
            }
            weights_cache[i] = std::vector<double>(points_cache[i].size(), 0.0);
            for (size_t j=0; j<points_cache[i].size(); j++) {
                // Integrate the function l_j(x) according to the measure induced by the function [weight_fn(x) + shift].
                for (size_t k=0; k<ref_points.size(); k++) {
                    weights_cache[i][j] += lagrange_eval(j, points_cache[i], ref_points[k]) * ref_weights[k];
                }
            }
        } else if (gauss_quadrature_version == 2) {
            // Use TasmanianTridiagonalSolver::decompose().
            std::vector<double> dummy_diag = diag;
            std::vector<double> dummy_offdiag = offdiag;
            TasmanianTridiagonalSolver::decompose(dummy_diag, dummy_offdiag, mu0, points_cache[i], weights_cache[i]);
        } else {
            throw std::invalid_argument("ERROR: gauss_quadrature_version must be a valid number!");
        }
        // Update the values of the polynomials at ref_points.
        std::swap(poly_m1_vals, poly_vals);
        std::transform(ref_points.begin(), ref_points.end(), poly_vals.begin(),
                       [&points_cache, i](double x){return poly_eval(points_cache[i], x);});
    }
}

/*!
 * \internal
 * \ingroup TasmanianAddonsExoticQuad
 * \brief Generate the full cache of Exotic (weighted and possibly shifted) quadrature weights and points for an (input) max level
 * \b n, using parameters \b shift, \b weight_fn_vals, \b shifted_weights, and \b ref_points.
 *
 * It is assumed that shifted_weights[i] is nonnegative and is already computed as shifted_weights[i] = ref_weights[i] *
 * (weight_fn_vals[i] + shift) for every i.
 * \endinternal
 */
inline TasGrid::CustomTabulated getShiftedExoticQuadrature(const int n, const double shift, const std::vector<double> &shifted_weights,
                                                           const std::vector<double> &ref_points, const char* description,
                                                           const bool is_symmetric = false) {
    // Create the set of points (roots) and weights for the first term.
    assert(shifted_weights.size() == ref_points.size());
    std::vector<std::vector<double>> points_cache, weights_cache;
    if (is_symmetric) {
        getGaussNodesAndWeights<true>(n, ref_points, shifted_weights, points_cache, weights_cache);
    } else {
        getGaussNodesAndWeights<false>(n, ref_points, shifted_weights, points_cache, weights_cache);
    }

    // Create and append the set of correction points and weights for the second term only if the shift is nonzero.
    if (shift != 0.0) {
        for (int i=0; i<n; i++) {
            const size_t init_size = points_cache[i].size();
            std::vector<double> correction_points, correction_weights;
            TasGrid::OneDimensionalNodes::getGaussLegendre(static_cast<int>(init_size), correction_weights, correction_points);
            for (auto &w : correction_weights) w *= -shift;
            if (correction_points.size() % 2 == 1 && is_symmetric) {
                // Zero out for stability.
                correction_points[(correction_points.size() - 1) / 2] = 0.0;
            }
            // Account for duplicates up to a certain tolerance.
            for (size_t j=0; j<correction_points.size(); j++) {
                size_t nonunique_idx = -1;
                for (size_t k=0; k<init_size; k++) {
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
        num_nodes[i] = static_cast<int>(points_cache[i].size());
        precision[i] = 2 * (i + 1) - 1;
    }
    return TasGrid::CustomTabulated(n, std::move(num_nodes), std::move(precision), std::move(points_cache),
                                    std::move(weights_cache), description);
}

/*!
 * \internal
 * \ingroup TasmanianAddonsExoticQuad
 * \brief Multiplies the reference quadrature weights by the shifted values of the weights function.
 *
 * \param vals are the values of the exotic weight function at the reference quadrature points
 * \param shift is the shift magnitude that ensures all \b vals[i] + \b shift are positive
 * \param ref_weights are the reference quadrature weights, will be overwritten with the result ref_weights[i] * (vals[i] + shift)
 *
 * \throws std::invalid_argument if \b vals[i] + \b shift is negative for some index
 *
 * \endinternal
 */
inline void shiftReferenceWeights(std::vector<double> const &vals, double shift, std::vector<double> &ref_weights){
    assert(ref_weights.size() == vals.size());
    if (not std::all_of(vals.begin(), vals.end(), [shift](double x){return (x+shift)>=0.0;}))
        throw std::invalid_argument("ERROR: shift needs to be chosen so that weight_fn(x)+shift is nonnegative!");
    for (size_t k=0; k<ref_weights.size(); k++)
        ref_weights[k] *= (vals[k] + shift);
}

/*
 * \ingroup TasmanianAddonsExoticQuad
 * \brief Constructs a CustomTabulated object containing the exotic quadrature nodes and quadrature weights for a given level \b n,
 * scalar \b shift, a reference TasmanianSparseGrid \b grid containing loaded values of the weight function, a string \b description,
 * and an optional bool \b is_symmetric that should be set to true if the weight function is symmetric.
 *
 * It is assumed that grid is a one dimensional interpolant/surrogate to the weight function. Only the nodes in grid.getLoadedPoints(),
 * the quadrature weights in grid.getQuadratureWeights(), and the loaded weight function values in grid.getLoadedValues() are used
 * in the computation of the exotic quadrature. Specifically, these vectors are used to orthogonalize the polynomial basis needed
 * to compute the optimal exotic quadrature.
 */
inline TasGrid::CustomTabulated getExoticQuadrature(const int n, const double shift, const TasGrid::TasmanianSparseGrid &grid,
                                                    const char* description, const bool is_symmetric = false) {
    if (grid.getNumLoaded() == 0) {
        throw std::invalid_argument("ERROR: grid needs to contain a nonzero number of loaded points returned by grid.getLoadedPoints()!");
    }
    if (grid.getNumDimensions() != 1) {
        throw std::invalid_argument("ERROR: grid needs to be one dimensional!");
    }
    if (grid.getNumOutputs() != 1) {
        throw std::invalid_argument("ERROR: grid needs to have scalar outputs!");
    }
    std::vector<double> ref_points = grid.getLoadedPoints();
    std::vector<double> shifted_weights = grid.getQuadratureWeights();
    std::vector<double> weight_fn_vals = std::vector<double>(grid.getLoadedValues(), grid.getLoadedValues() + grid.getNumLoaded());
    shiftReferenceWeights(weight_fn_vals, shift, shifted_weights);
    return getShiftedExoticQuadrature(n, shift, shifted_weights, ref_points, description, is_symmetric);
}

/*
 * \ingroup TasmanianAddonsExoticQuad
 * \brief Constructs a CustomTabulated object containing the exotic quadrature nodes and quadrature weights for a given level \b n,
 * scalar \b shift, a 1D function \b weight_fn() that returns the value of the weight function at its input, a positive integer
 * \b nref that specifies how many points are used in orthogonalizing the polynomial basis, a string \b description, and an optional
 * bool \b is_symmetric that should be set to true if the weight function is symmetric.
 *
 * The reference points and weights used for orthogonalizing the polynomial basis are generated by the Gauss-Legendre rule
 * with \b nref points and the weight function values are generated by applying \b weight_fn() to these reference points.
 */
inline TasGrid::CustomTabulated getExoticQuadrature(const int n, const double shift, std::function<double(double)> weight_fn,
                                                    const int nref, const char* description, const bool is_symmetric = false) {
    std::vector<double> shifted_weights, ref_points, weight_fn_vals(nref);
    TasGrid::OneDimensionalNodes::getGaussLegendre(nref, shifted_weights, ref_points);
    std::transform(ref_points.begin(), ref_points.end(), weight_fn_vals.begin(), weight_fn);
    shiftReferenceWeights(weight_fn_vals, shift, shifted_weights);
    return getShiftedExoticQuadrature(n, shift, shifted_weights, ref_points, description, is_symmetric);
}

} // namespace(TasGrid)
#endif
