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

#include "tsgEnumerates.hpp"

/*!
 * \internal
 * \file tsgBlasWrappers.hpp
 * \brief Wrappers to BLAS functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * The header contains a inline wrappers that give C++ style of
 * interface to BLAS operations.
 * \endinternal
 */

#ifndef __TASMANIAN_DOXYGEN_SKIP
extern "C"{
// Skip the definitions from Doxygen, this serves as a mock-up header for the BLAS API.
// BLAS level 1
double dnrm2_(const int *N, const double *x, const int *incx);
void dswap_(const int *N, double *x, const int *incx, double *y, const int *incy);
void dscal_(const int *N, const double *alpha, const double *x, const int *incx);
// BLAS level 2
void dgemv_(const char *transa, const int *M, const int *N, const double *alpha, const double *A, const int *lda,
            const double *x, const int *incx, const double *beta, const double *y, const int *incy);
void dtrsv_(const char *uplo, const char *trans, const char *diag, const int *N, const double *A, const int *lda,
            double *x, const int *incx);
// BLAS level 3
void dgemm_(const char* transa, const char* transb, const int *m, const int *n, const int *k, const double *alpha,
            const double *A, const int *lda, const double *B, const int *ldb, const double *beta, const double *C, const int *ldc);
void dtrsm_(const char *side, const char *uplo, const char *trans, const char *diag, const int *M, const int *N,
            const double *alpha, const double *A, const int *lda, double *B, const int *ldb);
void ztrsm_(const char *side, const char *uplo, const char *trans, const char *diag, const int *M, const int *N,
            const std::complex<double> *alpha, const std::complex<double> *A, const int *lda, std::complex<double> *B, const int *ldb);
// LAPACK solvers
// General PLU factorize/solve
void dgetrf_(const int *M, const int *N, double *A, const int *lda, int *ipiv, int *info);
void dgetrs_(const char *trans, const int *N, const int *nrhs, const double *A, const int *lda, const int *ipiv, double *B, const int *ldb, int *info);
// General least-squares solve
void dgels_(const char *trans, const int *M, const int *N, const int *nrhs, double *A, const int *lda,
            double *B, const int *ldb, double *work, int *lwork, int *info);
void zgels_(const char *trans, const int *M, const int *N, const int *nrhs, std::complex<double> *A, const int *lda,
            std::complex<double> *B, const int *ldb, std::complex<double> *work, int *lwork, int *info);
// Symmetric tridiagonal eigenvalue compute
void dstebz_(const char *range, const char *order, const int *N, const double *vl, const double *vu, const int *il, const int *iu, const double *abstol,
             const double D[], const double E[], int *M, int *nsplit, double W[], int iblock[], int isplit[], double work[], int iwork[], int *info);
void dsteqr_(const char *compz, const int *N, double D[], double E[], double Z[], const int *ldz, double work[], int *info);
void dsterf_(const int *N, double D[], double E[], int *info);
// General LQ-factorize and multiply by Q
#ifdef Tasmanian_BLAS_HAS_ZGELQ
void dgelq_(const int *M, const int *N, double *A, const int *lda, double *T, int const *Tsize, double *work, int const *lwork, int *info);
void dgemlq_(const char *side, const char *trans, const int *M, const int *N, const int *K, double const *A, int const *lda,
             double const *T, int const *Tsize, double C[], int const *ldc, double *work, int const *lwork, int *info);
void zgelq_(const int *M, const int *N, std::complex<double> *A, const int *lda, std::complex<double> *T, int const *Tsize,
            std::complex<double> *work, int const *lwork, int *info);
void zgemlq_(const char *side, const char *trans, const int *M, const int *N, const int *K, std::complex<double> const *A, int const *lda,
             std::complex<double> const *T, int const *Tsize, std::complex<double> C[], int const *ldc, std::complex<double> *work, int const *lwork, int *info);
#endif
}
#endif

/*!
 * \ingroup TasmanianTPLWrappers
 * \brief Wrappers for BLAS and LAPACK methods (hidden internal namespace).
 */
namespace TasBLAS{
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief BLAS dnrm2
     */
    inline double norm2(int N, double const x[], int incx){
        return dnrm2_(&N, x, &incx);
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief BLAS dswap
     */
    inline void vswap(int N, double x[], int incx, double y[], int incy){
        dswap_(&N, x, &incx, y, &incy);
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief BLAS dscal
     */
    inline void scal(int N, double alpha, double x[], int incx){
        dscal_(&N, &alpha, x, &incx);
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief BLAS dgemv
     */
    inline void gemv(char trans, int M, int N, double alpha, double const A[], int lda, double const x[], int incx,
                     double beta, double y[], int incy){
        dgemv_(&trans, &M, &N, &alpha, A, &lda, x, &incx, &beta, y, &incy);
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief BLAS dtrsv
     */
    inline void trsv(char uplo, char trans, char diag, int N, double const A[], int lda, double x[], int incx){
        dtrsv_(&uplo, &trans, &diag, &N, A, &lda, x, &incx);
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief BLAS gemm
     */
    inline void gemm(char transa, char transb, int M, int N, int K, double alpha, double const A[], int lda, double const B[], int ldb,
                     double beta, double C[], int ldc){
        dgemm_(&transa, &transb, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief BLAS dtrsm
     */
    inline void trsm(char side, char uplo, char trans, char diag, int M, int N, double alpha, double const A[], int lda, double B[], int ldb){
        dtrsm_(&side, &uplo, &trans, &diag, &M, &N, &alpha, A, &lda, B, &ldb);
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief BLAS ztrsm
     */
    inline void trsm(char side, char uplo, char trans, char diag, int M, int N, std::complex<double> alpha,
                     std::complex<double> const A[], int lda, std::complex<double> B[], int ldb){
        ztrsm_(&side, &uplo, &trans, &diag, &M, &N, &alpha, A, &lda, B, &ldb);
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK dgetrf
     */
    inline void getrf(int M, int N, double A[], int lda, int ipiv[]){
        int info = 0;
        dgetrf_(&M, &N, A, &lda, ipiv, &info);
        if (info != 0) throw std::runtime_error(std::string("Lapack dgetrf_ exited with code: ") + std::to_string(info));
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK dgetrs
     */
    inline void getrs(char trans, int N, int nrhs, double const A[], int lda, int const ipiv[], double B[], int ldb){
        int info = 0;
        dgetrs_(&trans, &N, &nrhs, A, &lda, ipiv, B, &ldb, &info);
        if (info != 0) throw std::runtime_error(std::string("Lapack dgetrs_ exited with code: ") + std::to_string(info));
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK dgels
     */
    inline void gels(char trans, int M, int N, int nrhs, double A[], int lda, double B[], int ldb, double work[], int lwork){
        int info = 0;
        dgels_(&trans, &M, &N, &nrhs, A, &lda, B, &ldb, work, &lwork, &info);
        if (info != 0){
            if (lwork > 0)
                throw std::runtime_error(std::string("Lapack dgels_ solve-stage exited with code: ") + std::to_string(info));
            else
                throw std::runtime_error(std::string("Lapack dgels_ infer-worksize-stage exited with code: ") + std::to_string(info));
        }
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK zgels
     */
    inline void gels(char trans, int M, int N, int nrhs, std::complex<double> A[], int lda, std::complex<double> B[], int ldb, std::complex<double> work[], int lwork){
        int info = 0;
        zgels_(&trans, &M, &N, &nrhs, A, &lda, B, &ldb, work, &lwork, &info);
        if (info != 0){
            if (lwork > 0)
                throw std::runtime_error(std::string("Lapack zgels_ solve-stage exited with code: ") + std::to_string(info));
            else
                throw std::runtime_error(std::string("Lapack zgels_ infer-worksize-stage exited with code: ") + std::to_string(info));
        }
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK dstebz
     */
    inline void stebz(char range, char order, int N, double vl, double vu, int il, int iu, double abstol, double D[], double E[],
                      int& M, int& nsplit, double W[], int iblock[], int isplit[], double work[], int iwork[]) {
        int info = 0;
        dstebz_(&range, &order, &N, &vl, &vu, &il, &iu, &abstol, D, E, &M, &nsplit, W, iblock, isplit, work, iwork, &info);
        if (info != 0) {
            if (info <= 3) {
                throw std::runtime_error(
                    std::string(
                        "Lapack dstebz_ failed to converge for some eigenvalues and exited with code: ") +
                    std::to_string(info));
            } else if (info == 4) {
                throw std::runtime_error(
                    std::string("Lapack dstebz_ used a Gershgorin interval that was too small and exited with code: ") +
                    std::to_string(info));
            } else if (info > 4) {
                throw std::runtime_error(
                    std::string("Lapack dstebz_ failed and exited with code: ") +
                    std::to_string(info));
            } else {
                throw std::runtime_error(
                    std::string(
                        "Lapack dstebz_ had an illegal value at argument number: ") +
                    std::to_string(-info));
            }
        }
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK dsteqr
     */
    inline void steqr(char compz, int N, double D[], double E[], double Z[], int ldz, double work[]) {
        int info = 0;
        dsteqr_(&compz, &N, D, E, Z, &ldz, work, &info);
        if (info != 0) {
            if (info > 0) {
                throw std::runtime_error(
                    std::string("Lapack dsteqr_ failed to converge for some eigenvalues and exited with code: ") +
                    std::to_string(info));
            } else {
                throw std::runtime_error(
                    std::string("Lapack dsteqr_ had an illegal value at argument number: ") +
                    std::to_string(-info));
            }
        }
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK dsterf
     */
    inline void sterf(int N, double D[], double E[]) {
        int info = 0;
        dsterf_(&N, D, E, &info);
        if (info != 0) {
            if (info > 0) {
                throw std::runtime_error(
                    std::string("Lapack dsteqr_ failed to converge for some eigenvalues and exited with code: ") +
                    std::to_string(info));
            } else {
                throw std::runtime_error(
                    std::string("Lapack dsteqr_ had an illegal value at argument number: ") +
                    std::to_string(-info));
            }
        }
    }
#ifdef Tasmanian_BLAS_HAS_ZGELQ
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK dgeql
     */
    inline void geql(int M, int N, double A[], int lda, double T[], int Tsize, double work[], int lwork){
        int info = 0;
        dgelq_(&M, &N, A, &lda, T, &Tsize, work, &lwork, &info);
        if (info != 0){
            if (lwork > 0)
                throw std::runtime_error(std::string("Lapack dgeql_ factorize-stage exited with code: ") + std::to_string(info));
            else
                throw std::runtime_error(std::string("Lapack dgeql_ infer-worksize-stage exited with code: ") + std::to_string(info));
        }
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK dgemlq
     */
    inline void gemlq(char side, char trans, int M, int N, int K, double const A[], int lda, double const T[], int Tsize,
                       double C[], int ldc, double work[], int lwork){
        int info = 0;
        dgemlq_(&side, &trans, &M, &N, &K, A, &lda, T, &Tsize, C, &ldc, work, &lwork, &info);
        if (info != 0){
            if (lwork > 0)
                throw std::runtime_error(std::string("Lapack dgemlq_ compute-stage exited with code: ") + std::to_string(info));
            else
                throw std::runtime_error(std::string("Lapack dgemlq_ infer-worksize-stage exited with code: ") + std::to_string(info));
        }
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK zgeql
     */
    inline void geql(int M, int N, std::complex<double> A[], int lda, std::complex<double> T[], int Tsize, std::complex<double> work[], int lwork){
        int info = 0;
        zgelq_(&M, &N, A, &lda, T, &Tsize, work, &lwork, &info);
        if (info != 0){
            if (lwork > 0)
                throw std::runtime_error(std::string("Lapack zgeql_ factorize-stage exited with code: ") + std::to_string(info));
            else
                throw std::runtime_error(std::string("Lapack zgeql_ infer-worksize-stage exited with code: ") + std::to_string(info));
        }
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief LAPACK zgemlq
     */
    inline void gemlq(char side, char trans, int M, int N, int K, std::complex<double> const A[], int lda, std::complex<double> const T[], int Tsize,
                       std::complex<double> C[], int ldc, std::complex<double> work[], int lwork){
        int info = 0;
        zgemlq_(&side, &trans, &M, &N, &K, A, &lda, T, &Tsize, C, &ldc, work, &lwork, &info);
        if (info != 0){
            if (lwork > 0)
                throw std::runtime_error(std::string("Lapack zgemlq_ compute-stage exited with code: ") + std::to_string(info));
            else
                throw std::runtime_error(std::string("Lapack zgemlq_ infer-worksize-stage exited with code: ") + std::to_string(info));
        }
    }
    #endif

    // higher-level methods building on top of one or more BLAS/LAPACK Methods

    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief Returns the square of the norm of the vector.
     */
    template<typename T>
    inline T norm2_2(int N, T const x[]){
        T nrm = norm2(N, x, 1);
        return nrm * nrm;
    }
    /*!
     * \ingroup TasmanianTPLWrappers
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
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief Conjugates a matrix, no op in the real case.
     */
    inline void conj_matrix(int, int, double[]){}
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief Conjugates the matrix, used in the case when 'T' operation is needed by only 'C' is available in the LAPACK standard.
     */
    inline void conj_matrix(int N, int M, std::complex<double> A[]){
        for(size_t i=0; i<static_cast<size_t>(N) * static_cast<size_t>(M); i++) A[i] = std::conj(A[i]);
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief Returns the transpose symbol, 'T' in the real case.
     */
    constexpr inline char get_trans(double){ return 'T'; }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief Returns the conjugate-transpose symbol, 'C' in the complex case.
     */
    constexpr inline char get_trans(std::complex<double>){ return 'C'; }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief Solves the over-determined least squares problem with single right-hand-side.
     *
     * Note that trans must be a capital letter N or T.
     */
    template<typename scalar_type>
    inline void solveLS(char trans, int N, int M, scalar_type A[], scalar_type b[], int nrhs = 1){
        std::vector<scalar_type> work(1);
        int n = (trans == 'N') ? N : M;
        int m = (trans == 'N') ? M : N;
        char effective_trans = (trans == 'N') ? trans : get_trans(static_cast<scalar_type>(0.0));
        conj_matrix(N, M, A); // does nothing in the real case, computes the conjugate in the complex one
        TasBLAS::gels(effective_trans, n, m, nrhs, A, n, b, N, work.data(), -1);
        work.resize(static_cast<size_t>(std::real(work[0])));
        TasBLAS::gels(effective_trans, n, m, nrhs, A, n, b, N, work.data(), static_cast<int>(work.size()));
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief Compute the LQ factorization of the matrix \b A.
     *
     * The assumption here is that the matrix is in column major format, otherwise this computes the QR factorization.
     * In fact, the Tasmanian uses this method to compute QR of a row-major matrix.
     */
    template<typename scalar_type>
    inline void factorizeLQ(int rows, int cols, scalar_type A[], std::vector<scalar_type> &T){
        T.resize(5);
        std::vector<scalar_type> work(1);
        geql(rows, cols, A, rows, T.data(), -1, work.data(), -1);
        T.resize(static_cast<size_t>(std::real(T[0])));
        work.resize(static_cast<size_t>(std::real(work[0])));
        geql(rows, cols, A, rows, T.data(), static_cast<int>(T.size()), work.data(), static_cast<int>(work.size()));
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief Multiplies C by the Q factor computed with factorizeLQ.
     *
     * Computes \f$ C = C Q^T \f$ where Q comes from the call to factorizeLQ.
     * The matrix C has dimensions M by N, A has dimensions K by N.
     */
    template<typename scalar_type>
    inline void multiplyQ(int M, int N, int K, scalar_type const A[], std::vector<scalar_type> const &T, scalar_type C[]){
        std::vector<scalar_type> work(1);
        gemlq('R', get_trans(static_cast<scalar_type>(0.0)), M, N, K, A, K, T.data(), static_cast<int>(T.size()), C, M, work.data(), -1);
        work.resize(static_cast<int>(std::real(work[0])));
        gemlq('R', get_trans(static_cast<scalar_type>(0.0)), M, N, K, A, K, T.data(), static_cast<int>(T.size()), C, M, work.data(), static_cast<int>(work.size()));
    }
    /*!
     * \ingroup TasmanianTPLWrappers
     * \brief Solves the least-squares assuming row-major format, see TasmanianDenseSolver::solvesLeastSquares()
     */
    template<typename scalar_type>
    void solveLSmulti(int n, int m, scalar_type A[], int nrhs, scalar_type B[]){
        if (nrhs == 1){
            TasBLAS::solveLS('T', n, m, A, B);
        }else{
            #ifdef Tasmanian_BLAS_HAS_ZGELQ
            std::vector<scalar_type> T;
            TasBLAS::factorizeLQ(m, n, A, T);
            TasBLAS::multiplyQ(nrhs, n, m, A, T, B);
            TasBLAS::trsm('R', 'L', 'N', 'N', nrhs, m, 1.0, A, m, B, nrhs);
            #else
            auto Bcols = TasGrid::Utils::transpose(nrhs, n, B);
            TasBLAS::solveLS('T', n, m, A, Bcols.data(), nrhs);
            TasGrid::Utils::transpose(n, nrhs, Bcols.data(), B);
            #endif
        }
    }
}

#endif
