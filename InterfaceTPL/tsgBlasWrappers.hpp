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

#ifndef __TASMANIAN_BLAS_WRAPPERS_HPP
#define __TASMANIAN_BLAS_WRAPPERS_HPP

/*!
 * \file tsgBlasWrappers.hpp
 * \brief Wrappers to BLAS functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * The header contains a inline wrappers that give C++ style of
 * interface to BLAS operations.
 */

/*!
 * \internal
 * \ingroup Tasmanian
 * \addtogroup TasmanianTPLWrappers Wrappers around TPL functionality
 *
 * Tasmanian uses multiple third party libraries (TPL) to gain advanced functionality,
 * such as optimized liner algebra on the CPU and GPU devices.
 * The libraries often come with C or Fortran style of API and the included C++ wrappers
 * help interface with the C++ internals of Tasmanian.
 * The wrappers are put in a set of private headers and should not be included
 * as part of the public API.
 *
 * \endinternal
 */

#ifndef __TASMANIAN_DOXYGEN_SKIP
// Skip the definitions from Doxygen, this serves as a mock-up header for the BLAS API.
extern "C" double dnrm2_(const int *N, const double *x, const int *incx);
extern "C" void dgemv_(const char *transa, const int *M, const int *N, const double *alpha, const double *A, const int *lda, const double *x, const int *incx, const double *beta, const double *y, const int *incy);
extern "C" void dgemm_(const char* transa, const char* transb, const int *m, const int *n, const int *k, const double *alpha, const double *A, const int *lda, const double *B, const int *ldb, const double *beta, const double *C, const int *ldc);
#endif

/*!
 * \brief Wrappers for BLAS methods.
 * \ingroup TasmanianTPLWrappers
 */
namespace TasBLAS{
    inline double norm2(int N, double const x[], int incx){
        return dnrm2_(&N, x, &incx);
    }
    //! \brief BLAS dgemv
    inline void gemv(char trans, int M, int N, double alpha, double const A[], int lda, double const x[], int incx,
                     double beta, double y[], int incy){
        dgemv_(&trans, &M, &N, &alpha, A, &lda, x, &incx, &beta, y, &incy);
    }
    //! \brief BLAS gemm
    inline void gemm(char transa, char transb, int M, int N, int K, double alpha, double const A[], int lda, double const B[], int ldb,
                     double beta, double C[], int ldc){
        dgemm_(&transa, &transb, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
    }

    //! \brief Returns the square of the norm of the vector.
    template<typename T>
    inline auto norm2_2(int N, T const x[]){
        T nrm = norm2(N, x, 1);
        return nrm * nrm;
    }
    /*!
     * \brief Combination of BLAS gemm and gemv
     *
     * Computes \f$ C = \alpha A B + \beta C \f$ where A is M by K, B is K by N, and C is M by N.
     * The method uses both gemm() and gemv() to handle the cases when either dimension is one.
     */
    template<typename T>
    inline void denseMultiply(int M, int N, int K, T alpha, const T A[], const T B[], T beta, T C[]){
        if (M > 1){
            if (N > 1){ // matrix mode
                gemm('N', 'N', M, N, K, alpha, A, M, B, K, beta, C, M);
            }else{ // matrix vector, A * v = C
                gemv('N', M, K, alpha, A, M, B, 1, beta, C, 1);
            }
        }else{ // matrix vector B^T * v = C
            gemv('T', K, N, alpha, B, K, A, 1, beta, C, 1);
        }
    }
}

#endif
