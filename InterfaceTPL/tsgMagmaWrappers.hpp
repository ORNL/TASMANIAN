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

//! \brief Calls magma_init() and sets the magma device.
inline void initMagma(AccelerationContext const *acceleration){
    if (not acceleration->engine->called_magma_init){
        magma_init();
        acceleration->engine->called_magma_init = std::unique_ptr<int>(new int(1));
    }
    magma_setdevice(acceleration->device);
}

//! \brief Wrapper around magma_dgeqrf_ooc()
inline void geqrf_ooc(int m, int n, double A[], int lda, double tau[]){
    int info = 0;
    double workspace_size;
    magma_dgeqrf_ooc(m, n, A, lda, tau, &workspace_size, -1, &info);
    if (info != 0)
        throw std::runtime_error("magma_dgeqrf_ooc() returned non-zero status: " + std::to_string(info) + " at size-query stage");
    double *workspace = nullptr;
    magma_dmalloc_pinned(&workspace, static_cast<size_t>(workspace_size));
    magma_dgeqrf_ooc(m, n, A, lda, tau, workspace, static_cast<int>(workspace_size), &info);
    if (info != 0)
        throw std::runtime_error("magma_dgeqrf_ooc() returned non-zero status: " + std::to_string(info) + " at execute stage");
    magma_free_pinned(workspace);
}
//! \brief Wrapper around magma_dgeqrf_ooc()
inline void geqrf_ooc(int m, int n, std::complex<double> A[], int lda, std::complex<double> tau[]){
    int info = 0;
    std::complex<double> workspace_size;
    magma_zgeqrf_ooc(m, n, reinterpret_cast<magmaDoubleComplex*>(A), lda, reinterpret_cast<magmaDoubleComplex*>(tau),
                     reinterpret_cast<magmaDoubleComplex*>(&workspace_size), -1, &info);
    if (info != 0)
        throw std::runtime_error("magma_zgeqrf_ooc() returned non-zero status: " + std::to_string(info) + " at size-query stage");
    magmaDoubleComplex *workspace = nullptr;
    magma_zmalloc_pinned(&workspace, static_cast<size_t>(std::real(workspace_size)));
    magma_zgeqrf_ooc(m, n, reinterpret_cast<magmaDoubleComplex*>(A), lda, reinterpret_cast<magmaDoubleComplex*>(tau),
                     workspace, static_cast<int>(std::real(workspace_size)), &info);
    if (info != 0)
        throw std::runtime_error("magma_zgeqrf_ooc() returned non-zero status: " + std::to_string(info) + " at execute stage");
    magma_free_pinned(workspace);
}

//! \brief Wrapper around magma_dormqr()
inline void gemqr_ooc(magma_side_t side, magma_trans_t trans, int m, int n, int k, double A[], int lda, double tau[], double C[], int ldc){
    int info = 0;
    std::vector<double> workspace(1);
    magma_dormqr(side, trans, m, n, k, A, lda, tau, C, ldc, workspace.data(), -1, &info);
    if (info != 0)
        throw std::runtime_error("magma_dormqr() returned non-zero status: " + std::to_string(info) + " at size-query stage");
    int wsize = static_cast<int>(workspace[0]);
    workspace.resize(static_cast<size_t>(wsize));
    magma_dormqr(side, trans, m, n, k, A, lda, tau, C, ldc, workspace.data(), wsize, &info);
    if (info != 0)
        throw std::runtime_error("magma_dormqr() returned non-zero status: " + std::to_string(info) + " at execute stage");
}
//! \brief Wrapper around magma_zunmqr()
inline void gemqr_ooc(magma_side_t side, magma_trans_t trans, int m, int n, int k, std::complex<double> A[], int lda, std::complex<double> tau[],
                      std::complex<double> C[], int ldc){
    int info = 0;
    std::vector<std::complex<double>> workspace(1);
    magma_zunmqr(side, trans, m, n, k, reinterpret_cast<magmaDoubleComplex*>(A), lda, reinterpret_cast<magmaDoubleComplex*>(tau),
                 reinterpret_cast<magmaDoubleComplex*>(C), ldc, reinterpret_cast<magmaDoubleComplex*>(workspace.data()), -1, &info);
    if (info != 0)
        throw std::runtime_error("magma_zunmqr() returned non-zero status: " + std::to_string(info) + " at size-query stage");
    int wsize = static_cast<int>(std::real(workspace[0]));
    workspace.resize(static_cast<size_t>(wsize));
    magma_zunmqr(side, trans, m, n, k, reinterpret_cast<magmaDoubleComplex*>(A), lda, reinterpret_cast<magmaDoubleComplex*>(tau),
                 reinterpret_cast<magmaDoubleComplex*>(C), ldc, reinterpret_cast<magmaDoubleComplex*>(workspace.data()), wsize, &info);
    if (info != 0)
        throw std::runtime_error("magma_zunmqr() returned non-zero status: " + std::to_string(info) + " at execute stage");
}

//! \brief Wrapper around magma_dtrsm_m().
inline void trsm_ooc(magma_side_t side, magma_uplo_t uplo,
                 magma_trans_t trans, magma_diag_t diag, int m, int n,
                 double alpha, double const A[], int lda, double B[], int ldb){
    magma_dtrsm_m(1, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
}
//! \brief Wrapper around magma_ztrsm_m().
inline void trsm_ooc(magma_side_t side, magma_uplo_t uplo,
                 magma_trans_t trans, magma_diag_t diag, int m, int n,
                 std::complex<double> alpha, std::complex<double> const A[], int lda, std::complex<double> B[], int ldb){
    magma_ztrsm_m(1, side, uplo, trans, diag, m, n, *reinterpret_cast<magmaDoubleComplex*>(&alpha),
                  reinterpret_cast<magmaDoubleComplex const*>(A), lda, reinterpret_cast<magmaDoubleComplex*>(B), ldb);
}

//! \brief Performsthe entire transpose operation.
template<typename scalar_type>
void solveLSmultiOOC(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]){
    initMagma(acceleration);
    std::vector<scalar_type> tau(m);
    auto AT = Utils::transpose(m, n, A);
    auto BT = Utils::transpose(nrhs, n, B);
    geqrf_ooc(n, m, AT.data(), n, tau.data());
    gemqr_ooc(MagmaLeft, (std::is_same<scalar_type, double>::value) ? MagmaTrans : MagmaConjTrans, n, nrhs, m, AT.data(), n, tau.data(), BT.data(), n);
    trsm_ooc(MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit, m, nrhs, 1.0, AT.data(), n, BT.data(), n);
    Utils::transpose(n, nrhs, BT.data(), B);
}

}
}

#endif
