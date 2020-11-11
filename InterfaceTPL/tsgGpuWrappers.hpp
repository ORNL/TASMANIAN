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

#ifndef __TASMANIAN_GPU_WRAPPERS_HPP
#define __TASMANIAN_GPU_WRAPPERS_HPP

/*!
 * \internal
 * \file tsgGpuWrappers.hpp
 * \brief Wrappers to GPU functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * The header contains definitions of various operations that can be performed
 * on the GPU devices with the corresponding GPU backend.
 * \endinternal
 */

#include "tsgAcceleratedDataStructures.hpp"

namespace TasGrid{
namespace TasGpu{

/*!
 * \brief Least squares solver with data sitting on the gpu device.
 *
 * Solves the least squares problem \f$ min_x \| A x - B \|^2_2 \f$ where \b A has \b n rows and \b m columns
 * (n >= m) and \b B has \b n rows and \b nrhs columns. Both matrices are stored in row-major format
 * on the GPU device associated with the given \b acceleration.
 * The matrix \b B will be overwritten with the solution.
 */
template<typename scalar_type>
void solveLSmultiGPU(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]);

/*!
 * \brief Identical to TasGpu::solveLSmultiGPU() but the arrays are on the CPU and the MAGMA out-of-core implementation is used.
 */
template<typename scalar_type>
void solveLSmultiOOC(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]);

//! \brief Identical to TasGpu::solveLSmultiGPU() but the data starts with the CPU and gets uploaded to the GPU first.
template<typename scalar_type>
void solveLSmulti(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]){
    GpuVector<scalar_type> gpuA(acceleration, m, n, A);
    GpuVector<scalar_type> gpuB(acceleration, nrhs, n, B);
    solveLSmultiGPU(acceleration, n, m, gpuA.data(), nrhs, gpuB.data());
    gpuB.unload(acceleration, B);
}

//! \brief Factorize \f$ A = P L U \f$, arrays are on the GPU.
void factorizePLU(AccelerationContext const *acceleration, int n, double A[], int_gpu_lapack ipiv[]);
//! \brief Solve A x = b using a PLU factorization.
void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int_gpu_lapack const ipiv[], double b[]);
//! \brief Solve A x = b using a PLU factorization, B is in row-major format.
void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int_gpu_lapack const ipiv[], int nrhs, double B[]);

/*!
 * \brief Wrapper to GPU BLAS that multiplies dense matrices (e.g., cuBlas, MAGMA).
 *
 * Computes \f$ C = \alpha A B + \beta C \f$ where \b A is M by K, \b B is K by N, and \b C is M by N
 * all stored in column major format on the GPU device associated with the \b acceleration.
 * The signature is near identical to BLAS sgemm() or dgemm() but there are no transpose variants
 * and the leading dimensions are inferred as if the matrices have no padding.
 */
template<typename scalar_type>
void denseMultiply(AccelerationContext const *acceleration, int M, int N, int K,
                   typename GpuVector<scalar_type>::value_type alpha, GpuVector<scalar_type> const &A,
                   GpuVector<scalar_type> const &B, typename GpuVector<scalar_type>::value_type beta, scalar_type C[]);

//! \brief Identical to TasGpu::denseMultiply() but both \b B and \b C are array in CPU memory.
template<typename scalar_type>
void denseMultiplyMixed(AccelerationContext const *acceleration, int M, int N, int K, typename GpuVector<scalar_type>::value_type alpha,
                        GpuVector<scalar_type> const &A, scalar_type const B[],
                        typename GpuVector<scalar_type>::value_type beta, scalar_type C[]){
    GpuVector<scalar_type> gpuB(acceleration, K, N, B), gpuC(acceleration, M, N);
    denseMultiply(acceleration, M, N, K, alpha, A, gpuB, beta, gpuC.data());
    gpuC.unload(acceleration, C);
}

/*!
 * \brief Wrapper to GPU methods that multiplies a sparse and a dense matrix.
 *
 * Computes \f$ C = \alpha A B \f$ where \b A is M by K, \b B is K by N, and \b C is M by N,
 * matrices are in column major format with \b B being sparse in column compressed format.
 */
template<typename scalar_type>
void sparseMultiply(AccelerationContext const *acceleration, int M, int N, int K, typename GpuVector<scalar_type>::value_type alpha,
                    const GpuVector<scalar_type> &A, const GpuVector<int> &pntr, const GpuVector<int> &indx,
                    const GpuVector<scalar_type> &vals, scalar_type C[]);

//! \brief Identical to TasGpu::sparseMultiply() but the sparse matrix and the result \b C are in CPU memory.
template<typename T>
void sparseMultiplyMixed(AccelerationContext const *acceleration, int M, int N, int K, typename GpuVector<T>::value_type alpha, const GpuVector<T> &A,
                         const std::vector<int> &pntr, const std::vector<int> &indx, const std::vector<T> &vals, T C[]){
    GpuVector<int> gpu_pntr(acceleration, pntr), gpu_indx(acceleration, indx);
    GpuVector<T> gpu_vals(acceleration, vals), gpu_c(acceleration, M, N);
    sparseMultiply(acceleration, M, N, K, alpha, A, gpu_pntr, gpu_indx, gpu_vals, gpu_c.data());
    gpu_c.unload(acceleration, C);
}

}
}

#endif
