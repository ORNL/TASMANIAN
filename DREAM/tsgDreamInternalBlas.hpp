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

#ifndef __TASMANIAN_DREAM_INTERNAL_BLAS_HPP
#define __TASMANIAN_DREAM_INTERNAL_BLAS_HPP

/*!
 * \internal
 * \file tsgDreamInternalBlas.hpp
 * \brief C++ inline functions that make calls to BLAS.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAM
 *
 * Wrappers that provide C++ API to BLAS-Fortran calls.
 */

/*!
 * \internal
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMBLAS C++ to BLAS wrappers
 *
 * Defined several inline functions that provde C++ wrappers for BLAS,
 * specifically translation is provided for the vector vs raw pointer
 * and the pass-by-value vs pass-by-reference API.
 * \endinternal
 */

#ifndef __TASMANIAN_DOXYGEN_SKIP
// avoiding including BLAS headers, use the standard to define the functions, only the BLAS library is needed without include directories
#ifdef Tasmanian_ENABLE_BLAS
//extern "C" void dtrsv_(const char *uplo, const char *trans, const char *diag, const int *N, const double *A, const int *lda, const double *x, const int *incx);
//extern "C" void dtrsm_(const char *side, const char *uplo, const char* transa, const char* diag, const int *m, const int *n, const double *alpha, const double *A, const int *lda, const double *B, const int *ldb);
extern "C" double dnrm2_(const int *N, const double *x, const int *incx);
extern "C" void dgemv_(const char *transa, const int *M, const int *N, const double *alpha, const double *A, const int *lda, const double *x, const int *incx, const double *beta, const double *y, const int *incy);
#endif
#endif

#include "TasmanianConfig.hpp"

namespace TasDREAM{

namespace TasBLAS{
#ifdef Tasmanian_ENABLE_BLAS
//! \internal
//! \brief Wrapper to BLAS 2-norm
//! \ingroup DREAMBLAS
inline double dnrm2squared(int N, const double x[]){
    int ione = 1;
    double nrm = dnrm2_(&N, x, &ione);
    return nrm * nrm;
}

//! \internal
//! \brief Wrapper to BLAS matrix-vector product, \b y = \b alpha * \b A-transpose * \b x + \b beta * \b y
//! \ingroup DREAMBLAS
inline void dgemtv(int M, int N, const double A[], const double x[], double y[], double alpha = 1.0, double beta = 0.0){ // y = A*x, A is M by N
    char charT = 'T'; int blas_one = 1;
    dgemv_(&charT, &M, &N, &alpha, A, &M, x, &blas_one, &beta, y, &blas_one);
}
#endif

}

}

#endif
