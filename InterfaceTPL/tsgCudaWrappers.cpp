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

#include "tsgGpuWrappers.hpp"

#ifdef Tasmanian_ENABLE_CUDA
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cusolverDn.h>
#else
#error "Cannot use tsgCudaWrappers.cpp without Tasmanian_ENABLE_CUDA"
#endif

#ifdef Tasmanian_ENABLE_MAGMA
#include "magma_v2.h"
#include "magmasparse.h"
#endif

/*!
 * \file tsgCudaWrappers.cpp
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
//! \brief Common error checking, prints a message based on the status and includes the \b info.
inline void cucheck(cudaError_t status, std::string info){
    if (status != cudaSuccess){
        std::string message = "ERROR: cuda failed at ";
        message += info + " with error: ";
        message += cudaGetErrorString(status);
        throw std::runtime_error(message);
    }
}

//! \brief Common error checking, prints a message based on the status and includes the \b info.
inline void cucheck(cublasStatus_t status, std::string info){
    if (status != CUBLAS_STATUS_SUCCESS){
        std::string message = "ERROR: cuBlas failed with code: ";
        if (status == CUBLAS_STATUS_NOT_INITIALIZED){
            message += "CUBLAS_STATUS_NOT_INITIALIZED";
        }else if (status == CUBLAS_STATUS_ALLOC_FAILED){
            message += "CUBLAS_STATUS_ALLOC_FAILED";
        }else if (status == CUBLAS_STATUS_INVALID_VALUE){
            message += "CUBLAS_STATUS_INVALID_VALUE";
        }else if (status == CUBLAS_STATUS_ARCH_MISMATCH){
            message += "CUBLAS_STATUS_ARCH_MISMATCH";
        }else if (status == CUBLAS_STATUS_MAPPING_ERROR){
            message += "CUBLAS_STATUS_MAPPING_ERROR";
        }else if (status == CUBLAS_STATUS_EXECUTION_FAILED){
            message += "CUBLAS_STATUS_EXECUTION_FAILED";
        }else if (status == CUBLAS_STATUS_INTERNAL_ERROR){
            message += "CUBLAS_STATUS_INTERNAL_ERROR";
        }else if (status == CUBLAS_STATUS_NOT_SUPPORTED){
            message += "CUBLAS_STATUS_NOT_SUPPORTED";
        }else if (status == CUBLAS_STATUS_LICENSE_ERROR){
            message += "CUBLAS_STATUS_LICENSE_ERROR";
        }else{
            message += "UNKNOWN";
        }
        message += " at ";
        message += info;
        throw std::runtime_error(message);
    }
}
//! \brief Common error checking, prints a message based on the status and includes the \b info.
inline void cucheck(cusolverStatus_t status, std::string info){
    if (status != CUSOLVER_STATUS_SUCCESS){
        std::string message = "ERROR: cuSolver failed with error code " + std::to_string(status) + " at " + info;
        throw std::runtime_error(message);
    }
}

/*
 * cuBLAS section
 */
//! \brief Make cuBlas handle.
inline cublasHandle_t getCuBlasHandle(GpuEngine *engine){
    if (engine->cublasHandle == nullptr){
        cucheck( cublasCreate(reinterpret_cast<cublasHandle_t*>(&engine->cublasHandle)), "create handle" );
        engine->own_cublas_handle = true;
    }
    return reinterpret_cast<cublasHandle_t>(engine->cublasHandle);
}

//! \brief Wrapper around dgeam().
void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
          int m, int n, double alpha, double const A[], int lda,
          double beta, double const B[], int ldb, double C[], int ldc){
    cucheck(cublasDgeam(handle, transa, transb, m, n, &alpha, A, lda, &beta, B, ldb, C, ldc), "cublasDgeam()");
}
//! \brief Wrapper around zgeam().
void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
          int m, int n, std::complex<double> alpha, std::complex<double> const A[], int lda,
          std::complex<double> beta, std::complex<double> const B[], int ldb, std::complex<double> C[], int ldc){
    cucheck(cublasZgeam(handle, transa, transb, m, n, reinterpret_cast<cuDoubleComplex*>(&alpha),
                        reinterpret_cast<cuDoubleComplex const*>(A), lda, reinterpret_cast<cuDoubleComplex*>(&beta),
                        reinterpret_cast<cuDoubleComplex const*>(B), ldb,
                        reinterpret_cast<cuDoubleComplex*>(C), ldc), "cublasZgeam()");
}

//! \brief Wrapper around dtrsv().
inline void trsv(cublasHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag,
                 int n, const double A[], int lda, double x[], int incx){
    cucheck(cublasDtrsv(handle, uplo, trans, diag, n, A, lda, x, incx), "cublasDtrsv()");
}
//! \brief Wrapper around ztrsv().
inline void trsv(cublasHandle_t handle, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag,
                 int n, const std::complex<double> A[], int lda, std::complex<double> x[], int incx){
    cucheck(cublasZtrsv(handle, uplo, trans, diag, n, reinterpret_cast<cuDoubleComplex const*>(A), lda,
                        reinterpret_cast<cuDoubleComplex*>(x), incx), "cublasZtrsv()");
}

//! \brief Wrapper around dtrsm().
inline void trsm(cublasHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
                 cublasOperation_t trans, cublasDiagType_t diag, int m, int n,
                 double alpha, double const A[], int lda, double B[], int ldb){
    cucheck(cublasDtrsm(handle, side, uplo, trans, diag, m, n, &alpha, A, lda, B, ldb), "cublasDtrsm()");
}
//! \brief Wrapper around ztrsm().
inline void trsm(cublasHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
                 cublasOperation_t trans, cublasDiagType_t diag, int m, int n,
                 std::complex<double> alpha, std::complex<double> const A[], int lda, std::complex<double> B[], int ldb){
    cucheck(cublasZtrsm(handle, side, uplo, trans, diag, m, n, reinterpret_cast<cuDoubleComplex*>(&alpha),
                        reinterpret_cast<cuDoubleComplex const*>(A), lda, reinterpret_cast<cuDoubleComplex*>(B), ldb), "cublasZtrsm()");
}

/*
 * cuSolver section
 */
//! \brief Make cuSolver handle.
inline cusolverDnHandle_t getCuSolverDnHandle(GpuEngine *engine){
    if (engine->cusolverDnHandle == nullptr){
        cucheck( cusolverDnCreate(reinterpret_cast<cusolverDnHandle_t*>(&engine->cusolverDnHandle)), "create handle" );
        engine->own_cusolverdn_handle = true;
    }
    return reinterpret_cast<cusolverDnHandle_t>(engine->cusolverDnHandle);
}

//! \brief Wrapper around dgeqrf() - get size.
inline size_t size_geqrf(cusolverDnHandle_t handle, int m, int n, double A[], int lda){
    int result = 0;
    cucheck(cusolverDnDgeqrf_bufferSize(handle, m, n, A, lda, &result), "cusolverDnDgeqrf_bufferSize()");
    return static_cast<size_t>(result);
}
//! \brief Wrapper around zgeqrf() - get size.
inline size_t size_geqrf(cusolverDnHandle_t handle, int m, int n, std::complex<double> A[], int lda){
    int result = 0;
    cucheck(cusolverDnZgeqrf_bufferSize(handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, &result), "cusolverDnZgeqrf_bufferSize()");
    return static_cast<size_t>(result);
}

//! \brief Wrapper around dgeqrf().
inline void geqrf(cusolverDnHandle_t handle, int m, int n, double A[], int lda, double tau[]){
    GpuVector<double> workspace( size_geqrf(handle, m, n, A, lda) );
    GpuVector<int> info(std::vector<int>(1, 0));
    cucheck(cusolverDnDgeqrf(handle, m, n, A, lda, tau, workspace.data(), static_cast<int>(workspace.size()), info.data()), "cusolverDnDgeqrf()");
    if (info.unload()[0] != 0)
        throw std::runtime_error("cusolverDnDgeqrf() returned non-zero status: " + std::to_string(info.unload()[0]));
}
//! \brief Wrapper around zgeqrf().
inline void geqrf(cusolverDnHandle_t handle, int m, int n, std::complex<double> A[], int lda, std::complex<double> tau[]){
    GpuVector<std::complex<double>> workspace( size_geqrf(handle, m, n, A, lda) );
    GpuVector<int> info(std::vector<int>(1, 0));
    cucheck(cusolverDnZgeqrf(handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda,
                             reinterpret_cast<cuDoubleComplex*>(tau), reinterpret_cast<cuDoubleComplex*>(workspace.data()),
                             static_cast<int>(workspace.size()), info.data()), "cusolverDnZgeqrf()");
    if (info.unload()[0] != 0)
        throw std::runtime_error("cusolverDnZgeqrf() returned non-zero status: " + std::to_string(info.unload()[0]));
}

//! \brief Wrapper around dormqr() - get size.
inline size_t size_gemqr(cusolverDnHandle_t handle, cublasSideMode_t side, cublasOperation_t trans,
                         int m, int n, int k, double const A[], int lda,
                         double const tau[], double const C[], int ldc){
    int result = 0;
    cucheck(cusolverDnDormqr_bufferSize(handle, side, trans, m, n, k, A, lda, tau, C, ldc, &result), "cusolverDnDormqr_bufferSize()");
    return static_cast<size_t>(result);
}
//! \brief Wrapper around zunmqr() - get size.
inline size_t size_gemqr(cusolverDnHandle_t handle, cublasSideMode_t side, cublasOperation_t trans,
                         int m, int n, int k, std::complex<double> const A[], int lda,
                         std::complex<double> const tau[], std::complex<double> const C[], int ldc){
    int result = 0;
    cucheck(cusolverDnZunmqr_bufferSize(handle, side, trans, m, n, k, reinterpret_cast<cuDoubleComplex const*>(A), lda,
                                        reinterpret_cast<cuDoubleComplex const*>(tau), reinterpret_cast<cuDoubleComplex const*>(C), ldc, &result),
            "cusolverDnZunmqr_bufferSize()");
    return static_cast<size_t>(result);
}

//! \brief Wrapper around dormqr().
inline void gemqr(cusolverDnHandle_t handle, cublasSideMode_t side, cublasOperation_t trans,
                  int m, int n, int k, double const A[], int lda, double const tau[], double C[], int ldc){
    GpuVector<int> info(std::vector<int>(1, 0));
    GpuVector<double> workspace( size_gemqr(handle, side, trans, m, n, k, A, lda, tau, C, ldc) );
    cucheck(cusolverDnDormqr(handle, side, trans, m, n, k, A, lda, tau, C, ldc,
                             workspace.data(), static_cast<int>(workspace.size()), info.data()), "cusolverDnDormqr()");
    if (info.unload()[0] != 0)
        throw std::runtime_error("cusolverDnDormqr() returned non-zero status: " + std::to_string(info.unload()[0]));
}
//! \brief Wrapper around zunmqr().
inline void gemqr(cusolverDnHandle_t handle, cublasSideMode_t side, cublasOperation_t trans,
                  int m, int n, int k, std::complex<double> const A[], int lda, std::complex<double> const tau[],
                  std::complex<double> C[], int ldc){
    GpuVector<int> info(std::vector<int>(1, 0));
    GpuVector<std::complex<double>> workspace( size_gemqr(handle, side, trans, m, n, k, A, lda, tau, C, ldc) );
    cucheck(cusolverDnZunmqr(handle, side, trans, m, n, k, reinterpret_cast<cuDoubleComplex const*>(A), lda,
                             reinterpret_cast<cuDoubleComplex const*>(tau), reinterpret_cast<cuDoubleComplex*>(C), ldc,
                             reinterpret_cast<cuDoubleComplex*>(workspace.data()), static_cast<int>(workspace.size()), info.data()),
            "cusolverDnZunmqr()");
    if (info.unload()[0] != 0)
        throw std::runtime_error("cusolverDnZunmqr() returned non-zero status: " + std::to_string(info.unload()[0]));
}

/*
 * Algorithm section
 */
//! \brief Solve a least-squares problem using row-major matrix format.
template<typename scalar_type>
void solveLSmultiGPU(GpuEngine *engine, int n, int m, scalar_type A[], int nrhs, scalar_type B[]){
    cublasHandle_t cublash = getCuBlasHandle(engine);
    cusolverDnHandle_t cusolverdnh = getCuSolverDnHandle(engine);

    GpuVector<scalar_type> AT(n, m);
    geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, n, m, 1.0, A, m, 0.0, A, m, AT.data(), n);

    GpuVector<scalar_type> T(m);
    geqrf(cusolverdnh, n, m, AT.data(), n, T.data());

    cublasOperation_t trans = (std::is_same<scalar_type, double>::value) ? CUBLAS_OP_T : CUBLAS_OP_C;

    if (nrhs == 1){
        gemqr(cusolverdnh, CUBLAS_SIDE_LEFT, trans, n, 1, m, AT.data(), n, T.data(), B, n);
        trsv(cublash, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, m, AT.data(), n, B, 1);
    }else{
        GpuVector<scalar_type> BT(n, nrhs);
        geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, n, nrhs, 1.0, B, nrhs, 0.0, B, nrhs, BT.data(), n);

        gemqr(cusolverdnh, CUBLAS_SIDE_LEFT, trans, n, nrhs, m, AT.data(), n, T.data(), BT.data(), n);
        trsm(cublash, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, m, nrhs, 1.0, AT.data(), n, BT.data(), n);

        geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, nrhs, n, 1.0, BT.data(), n, 0.0, BT.data(), n, B, nrhs);
    }
}

template void solveLSmultiGPU<double>(GpuEngine*, int, int, double[], int, double[]);
template void solveLSmultiGPU<std::complex<double>>(GpuEngine*, int, int, std::complex<double>[], int, std::complex<double>[]);

}
}

#endif
