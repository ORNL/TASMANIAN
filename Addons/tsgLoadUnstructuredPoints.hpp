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

#ifndef __TASMANIAN_ADDONS_LOADUNSTRUCTURED_HPP
#define __TASMANIAN_ADDONS_LOADUNSTRUCTURED_HPP

/*!
 * \internal
 * \file tsgLoadUnstructuredPoints.hpp
 * \brief Templates for using unstructured data.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianAddonsCommon
 *
 * Templates that infer the surrogate coefficients from unstructured data.
 * \endinternal
 */

#include "tsgLoadNeededValues.hpp"

/*!
 * \ingroup TasmanianAddons
 * \addtogroup TasmanianAddonsLoadUnstructured Load from Unstructured Point Set
 *
 * Templates that infer sparse grid coefficients from a set of unstructured points and model values.
 * These work very similar to TasGrid::loadNeededPoints(), but the data is not assumed to align
 * to the sparse grid. While removing some of the restrictions, the unstructured approach
 * always requires more data to achieve the same level of accuracy compared to the carefully chosen
 * points, and the approximation is valid only in the convex hull of the data, i.e.,
 * extrapolating outside of the data-cloud is not mathematically stable.
 * Furthermore, the inference usually relies on solving some linear or non-linear problem
 * which may not be stable and may require additional regularization, especially if the data
 * is not sufficient to provide values to all points, e.g., all data falls outside of the support
 * of some locally supported basis functions.
 */

namespace TasGrid{


template<typename scalar_type>
inline void loadUnstructuredDataL2(double const data_points[], int num_data, double const model_values[],
                                   double tolerance, TasmanianSparseGrid &grid){
    #ifndef Tasmanian_ENABLE_BLAS
    throw std::runtime_error("The loadUnstructuredDataL2() method requires Tasmanian_ENABLE_BLAS=ON in CMake!");
    #endif
    if (grid.getNumNeeded() != 0)
        grid.mergeRefinement();
    int num_equations = (tolerance > 0.0) ? num_data + grid.getNumPoints() : num_data;

    Data2D<scalar_type> basis_matrix(grid.getNumPoints(), num_equations, 0.0);
    Data2D<scalar_type> coefficients(grid.getNumOutputs(), num_equations, 0.0);

    grid.evaluateHierarchicalFunctions(data_points, num_data, reinterpret_cast<double*>(basis_matrix.getStrip(0)));
    if (tolerance > 0.0){
        double correction = std::sqrt(tolerance);
        for(int i=0; i<grid.getNumPoints(); i++)
            basis_matrix.getStrip(i + num_data)[i] = correction;
    }

    auto icoeff = coefficients.begin();
    for(size_t i=0; i<Utils::size_mult(num_data, grid.getNumOutputs()); i++)
        *icoeff++ = model_values[i];

    AccelerationContext acc;
    TasmanianDenseSolver::solvesLeastSquares(&acc, num_equations, grid.getNumPoints(),
                                             basis_matrix.data(), grid.getNumOutputs(), coefficients.data());

    if (std::is_same<scalar_type, std::complex<double>>::value){
        std::vector<double> real_coeffs(Utils::size_mult(2 * grid.getNumOutputs(), grid.getNumPoints()));
        icoeff = coefficients.begin();
        for(size_t i=0; i<Utils::size_mult(grid.getNumOutputs(), grid.getNumPoints()); i++)
            real_coeffs[i] = std::real(*icoeff++);
        icoeff = coefficients.begin();
        for(size_t i=Utils::size_mult(grid.getNumOutputs(), grid.getNumPoints()); i<Utils::size_mult(2 * grid.getNumOutputs(), grid.getNumPoints()); i++)
            real_coeffs[i] = std::imag(*icoeff++);
        grid.setHierarchicalCoefficients(real_coeffs.data());
    }else{
        grid.setHierarchicalCoefficients(reinterpret_cast<double*>(coefficients.data()));
    }
}

inline void loadUnstructuredDataL2(double const data_points[], int num_data, double const model_values[],
                                   double tolerance, TasmanianSparseGrid &grid){
    if (grid.isFourier()){
        loadUnstructuredDataL2<std::complex<double>>(data_points, num_data, model_values, tolerance, grid);
    }else{
        loadUnstructuredDataL2<double>(data_points, num_data, model_values, tolerance, grid);
    }
}

inline void loadUnstructuredDataL2(std::vector<double> const &data_points, std::vector<double> const &model_values,
                                   double tolerance, TasmanianSparseGrid &grid){
    int num_data = static_cast<int>(data_points.size() / grid.getNumDimensions());
    if (model_values.size() < Utils::size_mult(num_data, grid.getNumOutputs()))
        throw std::runtime_error("In loadUnstructuredDataL2(), provided more points than data.");
    loadUnstructuredDataL2(data_points.data(), num_data, model_values.data(), tolerance, grid);
}

}

#endif
