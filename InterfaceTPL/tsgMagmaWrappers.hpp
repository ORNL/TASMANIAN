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

#ifndef __TASMANIAN_MAGMA_WRAPPERS_HPP
#define __TASMANIAN_MAGMA_WRAPPERS_HPP

#include "tsgGpuWrappers.hpp"

#ifdef Tasmanian_ENABLE_MAGMA
#include "magma_v2.h"
#include "magmasparse.h"
#else
#error "Cannot use tsgMagmaWrappers.hpp without Tasmanian_ENABLE_MAGMA"
#endif

/*!
 * \file tsgMagmaWrappers.cpp
 * \brief Wrappers to CUDA functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * Realizations of the GPU algorithms using the CUDA backend.
 */

namespace TasGrid{
namespace TasGpu{

/*
 * Common methods
 */

//! \brief Wrapper around magma_sgemv().
inline void gemv(magma_queue_t mqueue, char transa, int M, int N, float alpha, float const A[], int lda,
                 float const x[], int incx, float beta, float y[], int incy){
    magma_sgemv(magma_trans_const(transa), M, N, alpha, A, lda, x, incx, beta, y, incy, mqueue);
}
//! \brief Wrapper around magma_dgemv().
inline void gemv(magma_queue_t mqueue, char transa, int M, int N, double alpha, double const A[], int lda,
                 double const x[], int incx, double beta, double y[], int incy){
    magma_dgemv(magma_trans_const(transa), M, N, alpha, A, lda, x, incx, beta, y, incy, mqueue);
}

//! \brief Wrapper around magma_sgemm().
inline void gemm(magma_queue_t mqueue, char transa, char transb, int M, int N, int K, float alpha, float const A[], int lda,
                 float const B[], int ldb, float beta, float C[], int ldc){
    magma_sgemm(magma_trans_const(transa), magma_trans_const(transb), M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, mqueue);
}
//! \brief Wrapper around magma_dgemm().
inline void gemm(magma_queue_t mqueue, char transa, char transb, int M, int N, int K, double alpha, double const A[], int lda,
                 double const B[], int ldb, double beta, double C[], int ldc){
    magma_dgemm(magma_trans_const(transa), magma_trans_const(transb), M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, mqueue);
}

}
}

#endif
