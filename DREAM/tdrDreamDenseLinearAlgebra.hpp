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

#ifndef __TASMANIAN_DREAM_DENSE_LINEAR_ALGEBRA_HPP
#define __TASMANIAN_DREAM_DENSE_LINEAR_ALGEBRA_HPP

#ifdef Tasmanian_ENABLE_BLAS
extern "C" void dtrsv_(const char *uplo, const char *trans, const char *diag, const int *N, const double *A, const int *lda, const double *x, const int *incx);
extern "C" void dtrsm_(const char *side, const char *uplo, const char* transa, const char* diag, const int *m, const int *n, const double *alpha, const double *A, const int *lda, const double *B, const int *ldb);
#endif

#include "TasmanianConfig.hpp"

namespace TasDREAM{

namespace TasBLAS{
#ifdef Tasmanian_ENABLE_BLAS
    inline void dtrsv_UTN(int N, const double A[], double B[]){ // solve A = U^{-T} B (A comes from cholUTU)
        char charU = 'U', charT = 'T', charN = 'N'; int ione = 1;
        dtrsv_(&charU, &charT, &charN, &N, A, &ione, B, &ione);
    }
    inline void dtrsv_UNN(int N, const double A[], double B[]){ // solve A = U^{-1} B (A comes from cholUTU)
        char charU = 'U', charN = 'N'; int ione = 1;
        dtrsv_(&charU, &charN, &charN, &N, A, &ione, B, &ione);
    }
    inline void dtrsm_LUTN(int M, int N, const double A[], const double B[]){ // solve A = U^{-T} B (A comes from cholUTU)
        char charL = 'L', charU = 'U', charT = 'T', charN = 'N'; double alpha = 1.0; int ione = 1;
        dtrsm_(&charL, &charU, &charT, &charN, &M, &N, &alpha, A, &ione, B, &ione);
    }
    inline void dtrsm_LUNN(int M, int N, const double A[], const double B[]){ // solve A = U^{-1} B (A comes from cholUTU)
        char charL = 'L', charU = 'U', charN = 'N'; double alpha = 1.0; int ione = 1;
        dtrsm_(&charL, &charU, &charN, &charN, &M, &N, &alpha, A, &ione, B, &ione);
    }
#else
    // non optimal BLAS subroutines, in case there is no BLAS available
    inline void dtrsv_UTN(int N, const double A[], double B[]){ // solve A = U^{-T} B (A comes from cholUTU)
        B[0] /= A[0];
        for(int i=1; i<N; i++){
            for(int j=0; j<i; j++) B[i] -= A[i*N + j] * B[i];
            B[i] /= A[i*N + i];
        }
    }
    inline void dtrsv_UNN(int N, const double A[], double B[]){ // solve A = U^{-1} B (A comes from cholUTU)
        for(int i=N-1; i>=0; i--){
            B[i] /= A[(i+1)*N-1];
            for(int j=i-1; j>=0; j--) B[j] -= A[i*N + j] * B[i];
        }
    }
    inline void dtrsm_LUTN(int M, int N, const double A[], double B[]){ // solve A = U^{-T} B (A comes from cholUTU)
        //char charL = 'L', charL = 'U', charL = 'T', charN = 'N'; double alpha = 1.0; int ione = 1;
        //dtrsm_(&charL, &charU, &charT, &charN, &M, &N, &alpha, A, &ione, B, &ione);
        for(int b=0; b<N; b++){ // for each column of B
            B[b*M] /= A[0];
            for(int i=1; i<M; i++){
                for(int j=0; j<i; j++) B[b*M + i] -= A[i*M + j] * B[b*M + i];
                B[b*M + i] /= A[i*M + i];
            }
        }
    }
    inline void dtrsm_LUNN(int M, int N, const double A[], double B[]){ // solve A = U^{-1} B (A comes from cholUTU)
        //char charL = 'L', charL = 'U', charN = 'N'; double alpha = 1.0; int ione = 1;
        //dtrsm_(&charL, &charU, &charN, &charN, &M, &N, &alpha, A, &ione, B, &ione);
        for(int b=0; b<N; b++){ // for each column of B
            for(int i=M-1; i>=0; i--){
                B[b*M+i] /= A[(i+1)*M-1];
                for(int j=i-1; j>=0; j--) B[b*M + j] -= A[i*M + j] * B[b*M + i];
            }
        }
    }
#endif // Tasmanian_ENABLE_BLAS

// used to compute anisotropic weights (at most 2 * dims X 2 * dims, not worth adding lapack dependence for such small matrix)
inline void cholUTU(int N, double A[], int lda = -1){
    if (lda == -1) lda = N;
    for(int i=0; i<N; i++){
        for(int k=0; k<i; k++) A[i*lda + i] -= A[i*lda +k]*A[i*lda +k];
        A[i*lda + i] = sqrt(A[i*lda + i]);
        for(int j=i+1; j<N; j++){
            for(int k=0; k<i; k++) A[j*lda + i] -= A[i*lda +k]*A[j*lda +k];
            A[j*lda + i] /= A[i*lda + i];
        }
    }
}

}

}

#endif