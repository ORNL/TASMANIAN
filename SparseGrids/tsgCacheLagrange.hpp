/*
 * Copyright (c) 2017, Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
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
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
 * THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
 * THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
 * IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TSG_CACHE_LAGRANGE_HPP
#define __TSG_CACHE_LAGRANGE_HPP

/*!
 * \internal
 * \file tsgCacheLagrange.hpp
 * \brief Cache data structure for evaluate with Global grids.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianAcceleration
 *
 * \endinternal
 */

#include "tsgGridCore.hpp"

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianAcceleration
 * \brief Cache that holds the values of 1D Lagrange polynomials.
 *
 * Global grids have the most complex evaluation procedure,
 * Lagrange polynomials have to be evaluated for each direction and \em each \em level,
 * before the weighted sum can be accumulated to compute the values of the basis functions.
 * The CacheLagrange class computes the values of the Lagrange polynomials
 * and caches them for easy access.
 * \endinternal
 */
template <typename T>
class CacheLagrange{
public:
    /*!
     * \brief Constructor that takes into account a single canonical point \b x.
     *
     * The cache is constructed for each dimension and each level up to \b max_levels,
     * the values of the Lagrange polynomials are computed in two passes resulting in O(n) operations.
     * - \b num_dimensions is the number of dimensions to consider
     * - \b max_levels indicates how many levels to consider in each direction,
     *   heavily anisotropic grids require only a few levels for the "less important" directions
     * - \b rule is the wrapper of the Global grid that contains information about number of points per level
     *   and the actual nodes with the pre-computed Lagrange coefficients
     * - \b holds the coordinates of the canonical point to cache
     */
    CacheLagrange(int num_dimensions, const std::vector<int> &max_levels, const OneDimensionalWrapper &rule, const double x[]) :
            cache(std::vector<std::vector<T>>(num_dimensions, std::vector<T>())), offsets(rule.getPointsCount()){
        for(int dim=0; dim<num_dimensions; dim++){
            cache[dim].resize(offsets[max_levels[dim] + 1]);
            for(int level=0; level <= max_levels[dim]; level++)
                cacheLevel(level, x[dim], rule, &(cache[dim][offsets[level]]));
        }
    }
    //! \brief Destructor, clear all used data.
    ~CacheLagrange() = default;

    //! \brief Computes the values of all Lagrange polynomials for the given level at the given x
    static void cacheLevel(int level, double x, const OneDimensionalWrapper &rule, T *cc){
        const double *nodes = rule.getNodes(level);
        const double *coeff = rule.getCoefficients(level);
        int num_points = rule.getNumPoints(level);

        cc[0] = 1.0;
        T c = 1.0;
        for(int j=0; j<num_points-1; j++){
            c *= (x - nodes[j]);
            cc[j+1] = c;
        }
        c = (rule.getRule() == rule_clenshawcurtis0) ? (x * x - 1.0) : 1.0;
        cc[num_points-1] *= c * coeff[num_points-1];
        for(int j=num_points-2; j>=0; j--){
            c *= (x - nodes[j+1]);
            cc[j] *= c * coeff[j];
        }
    }

    //! \brief Return the Lagrange cache for given \b dimension, \b level and offset local to the level
    T getLagrange(int dimension, int level, int local) const{
        return cache[dimension][offsets[level] + local];
    }

protected:
    /*!
     * \brief Stores the values of the Lagrange polynomias for each dimension and each level
     *
     * cache[i] corresponds to the i-th dimension.
     * The i-th vector has the values of the Lagrange polynomials going level by level,
     * so the value of the Lagrange polynomial for dimension i at level l with local offset j
     * will be cache[i][ offsets[l] + j ]
     */
    std::vector<std::vector<T>> cache;
    //! \brief Stores the offsets for each level within each dimension
    std::vector<int> offsets;
};


/*!
 * \internal
 * \ingroup TasmanianAcceleration
 * \brief Cache that holds the derivatives of 1D Lagrange polynomials. Uses the same interface as CacheLagrange.
 * \endinternal
 */

template <typename T>
class CacheLagrangeDerivative {
public:
    /*!
     * \brief Constructor that takes into account a single canonical point \b x.
     *
     * The cache is constructed for each dimension and each level up to \b max_levels,
     * the derivatives of the Lagrange polynomials are computed in two passes resulting in O(n) operations.
     * - \b num_dimensions is the number of dimensions to consider
     * - \b max_levels indicates how many levels to consider in each direction,
     *   heavily anisotropic grids require only a few levels for the "less important" directions
     * - \b rule is the wrapper of the Global grid that contains information about number of points per level
     *   and the actual nodes with the pre-computed Lagrange coefficients
     * - \b holds the coordinates of the canonical point to cache
     */
    CacheLagrangeDerivative(int num_dimensions, const std::vector<int> &max_levels, const OneDimensionalWrapper &rule, const double x[]) :
            cache(std::vector<std::vector<T>>(num_dimensions, std::vector<T>())), offsets(rule.getPointsCount()) {
        for(int dim=0; dim<num_dimensions; dim++){
            cache[dim].resize(offsets[max_levels[dim] + 1]);
            for(int level=0; level <= max_levels[dim]; level++)
                cacheDerivativeLevel(level, x[dim], rule, &(cache[dim][offsets[level]]));
        }
    }
    //! \brief Destructor, clear all used data.
    ~CacheLagrangeDerivative() = default;

    //! \brief Computes the derivatives of all Lagrange polynomials for the given level at the given x
    static void cacheDerivativeLevel(int level, double x, const OneDimensionalWrapper &rule, T *cc){
        /*
         * The usual j-th Lagrange derivative is of the form Lj(x) = cj * fj(x) * gj(x) where
         *
         *   fj(x) = (x - nodes[0])...(x - nodes[j-1]),
         *   gj(x) = (x - nodes[j+1])...(x - nodes[num_points-1]),
         *
         * and cj is the j-th Lagrange coefficient. The product rule gives Lj'(x) = cj * [fj'(x) * gj(x) + fj(x) * gj'(x)].
         * Hence, we do a forward pass for [fj(x), fj'(x)] and a backward pass for [gj(x), gj'(x)].
         */
        const double *nodes = rule.getNodes(level);
        const double *coeff = rule.getCoefficients(level);
        int num_points = rule.getNumPoints(level);
        // cc first stores fj'(x), aux_f stores fj(x), and aux_g stores gj(x).
        std::vector<T> aux_f(num_points), aux_g(num_points);
        cc[0] = (rule.getRule() == rule_clenshawcurtis0) ? 2.0 * x : 0.0;
        aux_f[0] = (rule.getRule() == rule_clenshawcurtis0) ? x * x - 1.0 : 1.0;
        aux_g[num_points-1] = 1.0;
        for(int i=1; i<num_points; i++) {
            aux_f[i] = aux_f[i-1] * (x - nodes[i-1]);
            aux_g[num_points-1-i] = aux_g[num_points-i] * (x - nodes[num_points-i]);
            cc[i] = aux_f[i-1] + (x - nodes[i-1]) * cc[i-1];
        }
        cc[num_points-1] *= coeff[num_points-1];
        T diff_gj = 0.0;
        for(int i=num_points-2; i>=0; i--) {
            diff_gj = aux_g[i+1] + diff_gj * (x - nodes[i+1]);
            cc[i] = coeff[i] * (cc[i] * aux_g[i] + aux_f[i] * diff_gj);
        }
    }

    //! \brief Return the Lagrange derivative cache for given \b dimension, \b level and offset local to the level
    T getLagrangeDerivative(int dimension, int level, int local) const {
        return cache[dimension][offsets[level] + local];
    }

protected:
    //! \brief See CacheLagrange::cache
    std::vector<std::vector<T>> cache;
    //! \brief See CacheLagrange::offsets
    std::vector<int> offsets;
};

}

#endif
