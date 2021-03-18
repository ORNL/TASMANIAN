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
 * In the interior of the data-cloud, the accuracy of the surrogate is very sensitive to
 * the distribution of the points and good accuracy can be achieved only in areas where
 * the points are sufficiently dense.
 *
 * The inference usually relies on solving some linear or non-linear problem
 * which may not be stable and may require additional regularization, especially if the data
 * is not sufficient to provide values for all hierarchical coefficients, e.g., all data falls outside of the support
 * of some locally supported basis functions.
 */

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianAddonsLoadUnstructured
 * \brief Returns true if the grid has a GPU algorithm for computing the hierarchical basis.
 *
 * \endinternal
 */
inline bool hasGPUBasis(TasmanianSparseGrid const &grid){
    return AccelerationMeta::isAvailable(accel_gpu_cuda)
        and not (grid.isLocalPolynomial() and ((grid.getOrder() < 0) or (grid.getOrder() > 2)))
        and not (grid.isWavelet() and grid.getOrder() == 3);
}

/*!
 * \internal
 * \ingroup TasmanianAddonsLoadUnstructured
 * \brief Generates the coefficients and solves the least-squares problem.
 *
 * The template assumes that the arrays are allocated on the GPU device
 * and the \b model_values are set to the correct size that would include padding
 * needed for the regularization (note that the padded may not be initialized).
 * The basis matrix will be constructed on the GPU and must fit in memory.
 * \endinternal
 */
template<typename scalar_type>
void generateCoefficientsGPU(double const data_points[], int num_data, scalar_type model_values[],
                             double tolerance, TasmanianSparseGrid &grid){
    AccelerationContext const *acceleration = grid.getAccelerationContext();
    int num_outputs = grid.getNumOutputs();
    int num_points = grid.getNumPoints();
    int num_equations = (tolerance > 0.0) ? num_data + num_points : num_data;

    GpuVector<scalar_type> basis_matrix(acceleration, num_equations, num_points);
    grid.evaluateHierarchicalFunctionsGPU(data_points, num_data, reinterpret_cast<double*>(basis_matrix.data()));

    if (tolerance > 0.0){
        double correction = std::sqrt(tolerance);
        constexpr long long esize = (std::is_same<scalar_type, double>::value) ? 1 : 2;
        long long num_total = static_cast<long long>(num_points) * static_cast<long long>(num_points) * esize;
        TasGpu::fillDataGPU(acceleration, 0.0, num_total, 1,
                            reinterpret_cast<double*>(basis_matrix.data() + Utils::size_mult(num_data, num_points)));

        long long stride = static_cast<long long>(num_points + 1) * esize;
        TasGpu::fillDataGPU(acceleration, correction, num_points, stride,
                            reinterpret_cast<double*>(basis_matrix.data() + Utils::size_mult(num_data, num_points)));

        TasGpu::fillDataGPU(acceleration, 0.0, esize * num_outputs * num_points, 1,
                            reinterpret_cast<double*>(model_values + Utils::size_mult(num_data, num_outputs)));
    }

    TasmanianDenseSolver::solvesLeastSquaresGPU(acceleration, num_equations, num_points,
                                                basis_matrix.data(), num_outputs, model_values);
}

/*!
 * \internal
 * \ingroup TasmanianAddonsLoadUnstructured
 * \brief Generates the coefficients and solves the least-squares problem.
 *
 * The template assumes that the arrays are allocated on the CPU.
 * The basis matrix is build on the CPU, unless MAGMA acceleration is enabled;
 * in the magma case, the matrix will be build on the GPU by using several
 * rows at a time as to avoid overfilling the GPU memory (out-of-core).
 * \endinternal
 */
template<typename scalar_type>
Data2D<scalar_type> generateCoefficients(double const data_points[], int num_data, double const model_values[], double tolerance, TasmanianSparseGrid &grid){
    int num_dimensions = grid.getNumDimensions();
    int num_outputs = grid.getNumOutputs();
    int num_points = grid.getNumPoints();
    int num_equations = (tolerance > 0.0) ? num_data + num_points : num_data;

    AccelerationContext const *acceleration = grid.getAccelerationContext();

    Data2D<scalar_type> basis_matrix(num_points, num_equations, 0.0);
    Data2D<scalar_type> coefficients(num_outputs, num_equations, 0.0);

    if (acceleration->mode == accel_gpu_magma and hasGPUBasis(grid)){
        int mem_usage = 268435456 / ((std::is_same<double, scalar_type>::value) ? 1 : 2); // 2GB of double numbers
        if (grid.getGPUMemory(acceleration->device) < 2048) mem_usage /= 2; // use minimum of 1GB
        int num_batch = std::min(num_data, mem_usage / num_points);
        if (num_batch < 128){
            // cannot use the GPU with less than 128 entries, might as well use the CPU
            grid.evaluateHierarchicalFunctions(data_points, num_data, reinterpret_cast<double*>(basis_matrix.data()));
        }else{
            acceleration->setDevice();
            GpuVector<double> gpu_points(acceleration, num_dimensions, num_batch);
            GpuVector<scalar_type> gpu_matrix(acceleration, num_points, num_batch);
            for(int i = 0; i < num_data; i += num_batch){
                int num_this_batch = std::min(num_batch, num_data - i);
                gpu_points.load(acceleration, Utils::size_mult(num_this_batch, num_dimensions),
                                data_points + Utils::size_mult(i, num_dimensions));
                grid.evaluateHierarchicalFunctionsGPU(gpu_points.data(), num_this_batch, reinterpret_cast<double*>(gpu_matrix.data()));
                gpu_matrix.unload(acceleration, Utils::size_mult(num_this_batch, num_points), basis_matrix.getStrip(i));
            }
        }
    }else{
        grid.evaluateHierarchicalFunctions(data_points, num_data, reinterpret_cast<double*>(basis_matrix.data()));
    }

    if (tolerance > 0.0){
        double correction = std::sqrt(tolerance);
        for(int i=0; i<grid.getNumPoints(); i++)
            basis_matrix.getStrip(i + num_data)[i] = correction;
    }

    auto icoeff = coefficients.begin();
    for(size_t i=0; i<Utils::size_mult(num_data, grid.getNumOutputs()); i++)
        *icoeff++ = model_values[i];

    TasmanianDenseSolver::solvesLeastSquares(acceleration, num_equations, grid.getNumPoints(),
                                             basis_matrix.data(), grid.getNumOutputs(), coefficients.data());

    return coefficients;
}

/*!
 * \internal
 * \ingroup TasmanianAddonsLoadUnstructured
 * \brief Template implementation that handles the case of Fourier grids vs. all other types.
 *
 * Fourier grids require complex arithmetic while the other grid works with real types.
 * This template operates on the double or std::complex<double> specified by \b scalar_type.
 * \endinternal
 */
template<typename scalar_type>
inline void loadUnstructuredDataL2tmpl(double const data_points[], int num_data, double const model_values[],
                                       double tolerance, TasmanianSparseGrid &grid){

    if (grid.empty()) throw std::runtime_error("Cannot use loadUnstructuredDataL2() with an empty grid.");
    if (grid.getNumNeeded() != 0)
        grid.mergeRefinement();

    AccelerationContext const *acceleration = grid.getAccelerationContext();
    if (acceleration->mode == accel_none)
        throw std::runtime_error("The loadUnstructuredDataL2() method cannot be used with acceleration mode accel_none.");

    int num_dimensions = grid.getNumDimensions();
    int num_outputs = grid.getNumOutputs();
    int num_points = grid.getNumPoints();
    int num_equations = (tolerance > 0.0) ? num_data + num_points : num_data;

    Data2D<scalar_type> coefficients =
        [&]()->Data2D<scalar_type>{
            if (acceleration->mode == accel_gpu_cuda and hasGPUBasis(grid)){
                acceleration->setDevice();
                GpuVector<double> gpu_points(acceleration, num_dimensions, num_data, data_points);
                GpuVector<scalar_type> gpu_values(acceleration, num_outputs, num_equations);
                TasGpu::load_n(acceleration, model_values, Utils::size_mult(num_outputs, num_data), gpu_values.data());
                generateCoefficientsGPU<scalar_type>(gpu_points.data(), num_data, gpu_values.data(), tolerance, grid);
                return Data2D<scalar_type>(num_outputs, num_equations, gpu_values.unload(acceleration));
            }else{
                return generateCoefficients<scalar_type>(data_points, num_data, model_values, tolerance, grid);
            }
        }();

    // do the real-complex split (used in Fourier grids)
    if (std::is_same<scalar_type, std::complex<double>>::value){
        std::vector<double> real_coeffs(Utils::size_mult(2 * grid.getNumOutputs(), grid.getNumPoints()));
        auto icoeff = coefficients.begin();
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

/*!
 * \ingroup TasmanianAddonsLoadUnstructured
 * \brief Construct a sparse grid surrogate using a least-squares fit.
 *
 * Constructs the coefficients of the grid using a least-squares fit \f$ \min_c \| G_c(x) - Y \|_2 \f$,
 * where \b c indicates the hierarchical coefficients of the grid, \b x indicates the data-points
 * and \b Y is the associated model outputs. The problem is set as a linear system of equations
 * so that \f$ G_c(x) = H c \f$ where \b H is the matrix of hierarchical coefficients.
 * Depending on the distribution of the points and/or the type of grid, the system may not be well posed,
 * hence the simple regularization option to add \b tolerance multiplied by the norm of \b c.
 *
 * \param data_points array of data points similar to TasmanianSparseGrid::loadConstructedPoints()
 *          with size grid.getNumDimensions() by \b num_data
 * \param num_data the number of data points
 * \param model_values array of model values corresponding to the \b data_points, the size is grid.getNumOutputs() by \b num_data
 * \param tolerance if positive, will add regularization to the least-squares problem, must be less than the accuracy
 *                  of the surrogate so that it will not affect the approximation error;
 *                  if set to zero or less, no regularization will be used but depending on the type of grid and
 *                  the distribution of the points, the problem may not be well posed
 * \param grid is a non-empty grid, if a refinement has been set, it will be merged with the existing points before inference
 *
 * \throws std::runtime_error if the grid is empty or \b accel_none is used
 * \throws std::runtime_error on various cuda errors, e.g., not enough GPU memory,
 *                            or if the matrix \b G_c does not have a full column rank, e.g., regularization tolerance
 *                            is \b 0 and the data is not sufficient to provide coefficients for all basis functions
 *
 * The problem requires the QR factorization of \b G_c and the solution against multiple right hand sides
 * corresponding to the model outputs. This is computationally very expensive and cannot be done without an
 * appropriate acceleration mode:
 * - \b accel_cpu_blas will perform all operations on the CPU using OpenMP and the BLAS and LAPACK implementations,
 *      OpenMP is optional but recommended
 * - \b accel_gpu_cublas will form \b G_c on the CPU, but then the matrix will be moved to the GPU device and
 *      will perform the solver stage using cuBlas and cuSolverDn
 * - \b accel_gpu_cuda will form \b G_c on the GPU and perform the solve using cuBlas and cuSolverDn,
 *      if the there is no GPU algorithm for the corresponding TasmanianSparseGrid::evaluateHierarchicalFunctionsGPU()
 *      then the CPU will be used identical to the cublas mode
 * - \b accel_gpu_magma will form \b G_c on the GPU device, but will do so in batches as to limit the GPU memory
 *      used at a given time, then the linear algebra will be performed using the MAGMA out-of-core methods
 *
 * The cuda and cublas modes, require that the data and the associated matrices fit in the memory of the GPU device,
 * MAGMA does not have that restriction.
 */
inline void loadUnstructuredDataL2(double const data_points[], int num_data, double const model_values[],
                                   double tolerance, TasmanianSparseGrid &grid){
    if (grid.isFourier()){
        loadUnstructuredDataL2tmpl<std::complex<double>>(data_points, num_data, model_values, tolerance, grid);
    }else{
        loadUnstructuredDataL2tmpl<double>(data_points, num_data, model_values, tolerance, grid);
    }
}

/*!
 * \ingroup TasmanianAddonsLoadUnstructured
 * \brief Overload that used vectors and infers the number of data points from the size of the vectors.
 *
 * See TasGrid::loadUnstructuredDataL2() for details, in addition
 *
 * \throws std::runtime_error if the sizes of the \b data_points and \b model_values do not match
 */
inline void loadUnstructuredDataL2(std::vector<double> const &data_points, std::vector<double> const &model_values,
                                   double tolerance, TasmanianSparseGrid &grid){
    if (grid.empty()) throw std::runtime_error("Cannot use loadUnstructuredDataL2() with an empty grid.");
    int num_data = static_cast<int>(data_points.size() / grid.getNumDimensions());
    if (model_values.size() < Utils::size_mult(num_data, grid.getNumOutputs()))
        throw std::runtime_error("In loadUnstructuredDataL2(), provided more points than data.");
    loadUnstructuredDataL2(data_points.data(), num_data, model_values.data(), tolerance, grid);
}

}

#endif
