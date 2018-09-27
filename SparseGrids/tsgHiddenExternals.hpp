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

namespace TasGrid{

// interfaces to external libraries, for convenience purposes mostly
// for now include only BLAS
#ifdef Tasmanian_ENABLE_BLAS
extern "C" void daxpy_(const int *N, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern "C" double ddot_(const int *N, const double *x, const int *incx, const double *y, const int *incy);
extern "C" double dnrm2_(const int *N, const double *x, const int *incx);
extern "C" void dgemv_(const char *transa, const int *M, const int *N, const double *alpha, const double *A, const int *lda, const double *x, const int *incx, const double *beta, const double *y, const int *incy);
extern "C" void dtrsv_(const char *uplo, const char *trans, const char *diag, const int *N, const double *A, const int *lda, const double *x, const int *incx);
extern "C" void dgemm_(const char* transa, const char* transb, const int *m, const int *n, const int *k, const double *alpha, const double *A, const int *lda, const double *B, const int *ldb, const double *beta, const double *C, const int *ldc);
extern "C" void dtrsm_(const char *side, const char *uplo, const char* transa, const char* diag, const int *m, const int *n, const double *alpha, const double *A, const int *lda, const double *B, const int *ldb);
#endif // Tasmanian_ENABLE_BLAS


class TasBLAS{
public:
#ifdef Tasmanian_ENABLE_BLAS
    // can't find BLAS implementation that is generally reliable
    // if dgemv_ is called within OpenMP region with large matrices, it returns wrong result!
    // Level 1
    inline static void daxpy(int N, double alpha, const double x[], double y[]){
        int ione = 1;
        daxpy_(&N, &alpha, x, &ione, y, &ione);
    }
    inline static double ddot(int N, const double x[], const double y[]){
        int ione = 1;
        return ddot_(&N, x, &ione, y, &ione);
    }
    inline static double ddot(int N, const double x[]){
        int ione = 1;
        double nrm = dnrm2_(&N, x, &ione);
        return nrm * nrm;
    }
    // Level 2
    inline static void dgemv(int M, int N, const double A[], const double x[], double y[], double alpha = 1.0, double beta = 0.0){ // y = A*x, A is M by N
        char charN = 'N'; int blas_one = 1;
        dgemv_(&charN, &M, &N, &alpha, A, &M, x, &blas_one, &beta, y, &blas_one);
    }
    inline static void dgemtv(int M, int N, const double A[], const double x[], double y[]){ // y = A^T*x, A is M by N
        char charT = 'T'; int blas_one = 1; double alpha = 1.0, beta = 0.0;
        dgemv_(&charT, &M, &N, &alpha, A, &M, x, &blas_one, &beta, y, &blas_one);
    }
    inline static void dtrsv_UTN(int N, const double A[], double B[]){ // solve A = U^{-T} B (A comes from cholUTU)
        char charU = 'U', charT = 'T', charN = 'N'; int ione = 1;
        dtrsv_(&charU, &charT, &charN, &N, A, &ione, B, &ione);
    }
    inline static void dtrsv_UNN(int N, const double A[], double B[]){ // solve A = U^{-1} B (A comes from cholUTU)
        char charU = 'U', charN = 'N'; int ione = 1;
        dtrsv_(&charU, &charN, &charN, &N, A, &ione, B, &ione);
    }
    // Level 3
    inline static void dgemm(int M, int N, int K, double alpha, const double A[], const double B[], double beta, double C[]){
        char charN = 'N';
        dgemm_(&charN, &charN, &M, &N, &K, &alpha, A, &M, B, &K, &beta, C, &M);
    }
    inline static void dtrsm_LUTN(int M, int N, const double A[], const double B[]){ // solve A = U^{-T} B (A comes from cholUTU)
        char charL = 'L', charU = 'U', charT = 'T', charN = 'N'; double alpha = 1.0; int ione = 1;
        dtrsm_(&charL, &charU, &charT, &charN, &M, &N, &alpha, A, &ione, B, &ione);
    }
    inline static void dtrsm_LUNN(int M, int N, const double A[], const double B[]){ // solve A = U^{-1} B (A comes from cholUTU)
        char charL = 'L', charU = 'U', charN = 'N'; double alpha = 1.0; int ione = 1;
        dtrsm_(&charL, &charU, &charN, &charN, &M, &N, &alpha, A, &ione, B, &ione);
    }
#else
    // non optimal BLAS subroutines, in case there is no BLAS available
    inline static void daxpy(int N, double alpha, const double x[], double y[]){
        for(int i=0; i<N; i++) y[i] += alpha * x[i];
    }
    inline static double ddot(int N, const double x[], const double y[]){
        double sum = 0.0;
        for(int i=0; i<N; i++) sum += x[i] * y[i];
        return sum;
    }
    inline static double ddot(int N, const double x[]){
        double sum = 0.0;
        for(int i=0; i<N; i++) sum += x[i] * x[i];
        return sum;
    }
    inline static void dgemm(int M, int N, int K, double alpha, const double A[], const double B[], double beta, double C[]){
        for(int j=0; j<N; j++){
            for(int i=0; i<M; i++){
                double sum = 0.0;
                for(int k=0; k<K; k++){
                    sum += A[k*M + i] * B[j*N + k];
                }
                C[j*M + i] = beta*C[j*M + i] + alpha * sum;
            }
        }
    }
    inline static void dgemv(int M, int N, const double A[], const double x[], double y[]){ // y = A*x, A is M by N
        for(int i=0; i<M; i++){
            y[i] = A[i] * x[0];
            for(int j=1; j<N; j++) y[i] += A[j*M + i] * x[j];
        }
    }
    inline static void dgemtv(int M, int N, const double A[], const double x[], double y[]){ // y = A^T*x, A^T is M by N
        for(int i=0; i<M; i++){
            y[i] = A[i*M] * x[0];
            for(int j=1; j<N; j++) y[i] += A[i*M + j] * x[j];
        }
    }
    inline static void dtrsv_UTN(int N, const double A[], double B[]){ // solve A = U^{-T} B (A comes from cholUTU)
        B[0] /= A[0];
        for(int i=1; i<N; i++){
            for(int j=0; j<i; j++) B[i] -= A[i*N + j] * B[i];
            B[i] /= A[i*N + i];
        }
    }
    inline static void dtrsv_UNN(int N, const double A[], double B[]){ // solve A = U^{-1} B (A comes from cholUTU)
        for(int i=N-1; i>=0; i--){
            B[i] /= A[(i+1)*N-1];
            for(int j=i-1; j>=0; j--) B[j] -= A[i*N + j] * B[i];
        }
    }
    inline static void dtrsm_LUTN(int M, int N, const double A[], double B[]){ // solve A = U^{-T} B (A comes from cholUTU)
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
    inline static void dtrsm_LUNN(int M, int N, const double A[], double B[]){ // solve A = U^{-1} B (A comes from cholUTU)
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
    inline static void cholUTU(int N, double A[], int lda = -1){
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


    inline static void setzero(int N, double *A){ std::fill(A, A + N, 0.0); } // fill A with N zeros
    inline static void setzero(int N, int *A){ std::fill(A, A + N, 0); } // fill A with N zeros

};

}

#endif
