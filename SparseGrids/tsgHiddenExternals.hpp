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

#ifndef __TASMANIAN_SPARSE_GRID_HIDDEN_INTERNALS_HPP
#define __TASMANIAN_SPARSE_GRID_HIDDEN_INTERNALS_HPP

#include "TasmanianConfig.hpp"

namespace TasGrid{

#ifdef Tasmanian_ENABLE_BLAS
#ifndef __TASMANIAN_DOXYGEN_SKIP
// Skip the definitions from Doxygen, this serves as a mock-up header for the BLAS API.
extern "C" void dgemv_(const char *transa, const int *M, const int *N, const double *alpha, const double *A, const int *lda, const double *x, const int *incx, const double *beta, const double *y, const int *incy);
extern "C" void dgemm_(const char* transa, const char* transb, const int *m, const int *n, const int *k, const double *alpha, const double *A, const int *lda, const double *B, const int *ldb, const double *beta, const double *C, const int *ldc);
#endif

//! \internal
//! \brief Wrappers for BLAS methods.
//! \ingroup TasmanianAcceleration
namespace TasBLAS{
    //! \internal
    //! \brief Wrapper for BLAS dense matrix-matrix and matrix-vector.

    //! Common API computing \f$ C = \alpha A B + \beta C \f$ where A is M by K, B is K by N, and C is M by N.
    //! Transposes are not considered (not needed by Tasmanian).
    //! The method switches between \b dgemm_ and \b dgemv_ depending on the appropriate dimensions.
    inline void denseMultiply(int M, int N, int K, double alpha, const double A[], const double B[], double beta, double C[]){
        if (M > 1){
            if (N > 1){ // matrix mode
                char charN = 'N';
                dgemm_(&charN, &charN, &M, &N, &K, &alpha, A, &M, B, &K, &beta, C, &M);
            }else{ // matrix vector, A * v = C
                char charN = 'N'; int blas_one = 1;
                dgemv_(&charN, &M, &K, &alpha, A, &M, B, &blas_one, &beta, C, &blas_one);
            }
        }else{ // matrix vector B^T * v = C
            char charT = 'T'; int blas_one = 1;
            dgemv_(&charT, &K, &N, &alpha, B, &K, A, &blas_one, &beta, C, &blas_one);
        }
    }
}
#endif


// interfaces to external libraries, for convenience purposes mostly
// for now include only BLAS
#ifdef Tasmanian_ENABLE_BLAS
extern "C" void daxpy_(const int *N, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern "C" double ddot_(const int *N, const double *x, const int *incx, const double *y, const int *incy);
extern "C" double dnrm2_(const int *N, const double *x, const int *incx);
//extern "C" void dgemv_(const char *transa, const int *M, const int *N, const double *alpha, const double *A, const int *lda, const double *x, const int *incx, const double *beta, const double *y, const int *incy);
extern "C" void dtrsv_(const char *uplo, const char *trans, const char *diag, const int *N, const double *A, const int *lda, const double *x, const int *incx);
//extern "C" void dgemm_(const char* transa, const char* transb, const int *m, const int *n, const int *k, const double *alpha, const double *A, const int *lda, const double *B, const int *ldb, const double *beta, const double *C, const int *ldc);
#endif // Tasmanian_ENABLE_BLAS

namespace TasBLAS{
#ifdef Tasmanian_ENABLE_BLAS
    // can't find BLAS implementation that is generally reliable
    // if dgemv_ is called within OpenMP region with large matrices, it returns wrong result!
    // Level 1
//     inline double ddot(int N, const double x[], const double y[]){
//         int ione = 1;
//         return ddot_(&N, x, &ione, y, &ione);
//     }
//     inline double ddot(int N, const double x[]){
//         int ione = 1;
//         double nrm = dnrm2_(&N, x, &ione);
//         return nrm * nrm;
//     }
    // Level 2
    inline void dgemv(int M, int N, const double A[], const double x[], double y[], double alpha = 1.0, double beta = 0.0){ // y = A*x, A is M by N
        char charN = 'N'; int blas_one = 1;
        dgemv_(&charN, &M, &N, &alpha, A, &M, x, &blas_one, &beta, y, &blas_one);
    }
    // Level 3
    inline void dgemm(int M, int N, int K, double alpha, const double A[], const double B[], double beta, double C[]){
        char charN = 'N';
        dgemm_(&charN, &charN, &M, &N, &K, &alpha, A, &M, B, &K, &beta, C, &M);
    }
#else
    // non optimal BLAS subroutines, in case there is no BLAS available
//     inline double ddot(int N, const double x[], const double y[]){
//         double sum = 0.0;
//         for(int i=0; i<N; i++) sum += x[i] * y[i];
//         return sum;
//     }
//     inline double ddot(int N, const double x[]){
//         double sum = 0.0;
//         for(int i=0; i<N; i++) sum += x[i] * x[i];
//         return sum;
//     }
    inline void dgemv(int M, int N, const double A[], const double x[], double y[]){ // y = A*x, A is M by N
        for(int i=0; i<M; i++){
            y[i] = A[i] * x[0];
            for(int j=1; j<N; j++) y[i] += A[j*M + i] * x[j];
        }
    }
#endif // Tasmanian_ENABLE_BLAS

    inline void setzero(int N, double *A){ std::fill(A, A + N, 0.0); } // fill A with N zeros
}

}

#endif
