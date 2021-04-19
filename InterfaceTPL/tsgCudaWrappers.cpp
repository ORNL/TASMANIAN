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

#ifndef __TASMANIAN_CUDA_WRAPPERS_CPP
#define __TASMANIAN_CUDA_WRAPPERS_CPP

#include "tsgGpuWrappers.hpp"
#include "tsgCudaWrappers.hpp"

/*!
 * \file tsgCudaWrappers.cpp
 * \brief Wrappers to CUDA functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * Realizations of the GPU algorithms using the CUDA backend.
 */

namespace TasGrid{
/*
 * Meta methods
 */
template<typename T> void GpuVector<T>::resize(AccelerationContext const*, size_t count){
    if (count != num_entries){ // if the current array is not big enough
        clear(); // resets dynamic_mode
        num_entries = count;
        TasGpu::cucheck( cudaMalloc(((void**) &gpu_data), num_entries * sizeof(T)), "cudaMalloc()");
    }
}
template<typename T> void GpuVector<T>::clear(){
    num_entries = 0;
    if (gpu_data != nullptr) // if I own the data and the data is not null
        TasGpu::cucheck( cudaFree(gpu_data), "cudaFree()");
    gpu_data = nullptr;
}
template<typename T> void GpuVector<T>::load(AccelerationContext const *acc, size_t count, const T* cpu_data){
    resize(acc, count);
    TasGpu::cucheck( cudaMemcpy(gpu_data, cpu_data, num_entries * sizeof(T), cudaMemcpyHostToDevice), "cudaMemcpy() to device");
}
template<typename T> void GpuVector<T>::unload(AccelerationContext const*, size_t num, T* cpu_data) const{
    TasGpu::cucheck( cudaMemcpy(cpu_data, gpu_data, num * sizeof(T), cudaMemcpyDeviceToHost), "cudaMemcpy() from device");
}

template void GpuVector<double>::resize(AccelerationContext const*, size_t);
template void GpuVector<double>::clear();
template void GpuVector<double>::load(AccelerationContext const*, size_t, const double*);
template void GpuVector<double>::unload(AccelerationContext const*, size_t, double*) const;

template void GpuVector<std::complex<double>>::resize(AccelerationContext const*, size_t);
template void GpuVector<std::complex<double>>::clear();
template void GpuVector<std::complex<double>>::load(AccelerationContext const*, size_t, const std::complex<double>*);
template void GpuVector<std::complex<double>>::unload(AccelerationContext const*, size_t, std::complex<double>*) const;

template void GpuVector<float>::resize(AccelerationContext const*, size_t);
template void GpuVector<float>::clear();
template void GpuVector<float>::load(AccelerationContext const*, size_t, const float*);
template void GpuVector<float>::unload(AccelerationContext const*, size_t, float*) const;

template void GpuVector<int>::resize(AccelerationContext const*, size_t);
template void GpuVector<int>::clear();
template void GpuVector<int>::load(AccelerationContext const*, size_t, const int*);
template void GpuVector<int>::unload(AccelerationContext const*, size_t, int*) const;

template<> void deleteHandle<AccHandle::Cublas>(int *p){ cublasDestroy(reinterpret_cast<cublasHandle_t>(p)); }
template<> void deleteHandle<AccHandle::Cusparse>(int *p){ cusparseDestroy(reinterpret_cast<cusparseHandle_t>(p)); }
template<> void deleteHandle<AccHandle::Cusolver>(int *p){ cusolverDnDestroy(reinterpret_cast<cusolverDnHandle_t>(p)); }

void GpuEngine::setCuBlasHandle(void *handle){
    cublas_handle = std::unique_ptr<int, HandleDeleter<AccHandle::Cublas>>
        (reinterpret_cast<int*>(handle), HandleDeleter<AccHandle::Cublas>(false));
}
void GpuEngine::setCuSparseHandle(void *handle){
    cusparse_handle = std::unique_ptr<int, HandleDeleter<AccHandle::Cusparse>>
        (reinterpret_cast<int*>(handle), HandleDeleter<AccHandle::Cusparse>(false));
}
void GpuEngine::setCuSolverDnHandle(void *handle){
    cusolver_handle = std::unique_ptr<int, HandleDeleter<AccHandle::Cusolver>>
        (reinterpret_cast<int*>(handle), HandleDeleter<AccHandle::Cusolver>(false));
}

int AccelerationMeta::getNumGpuDevices(){
    int gpu_count = 0;
    cudaGetDeviceCount(&gpu_count);
    return gpu_count;
}
void AccelerationMeta::setDefaultGpuDevice(int deviceID){
    cudaSetDevice(deviceID);
}
unsigned long long AccelerationMeta::getTotalGPUMemory(int deviceID){
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, deviceID);
    return prop.totalGlobalMem;
}
std::string AccelerationMeta::getGpuDeviceName(int deviceID){
    if ((deviceID < 0) || (deviceID >= getNumGpuDevices())) return std::string();

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, deviceID);

    return std::string(prop.name);
}
template<typename T> void AccelerationMeta::recvGpuArray(AccelerationContext const*, size_t num_entries, const T *gpu_data, std::vector<T> &cpu_data){
    cpu_data.resize(num_entries);
    TasGpu::cucheck( cudaMemcpy(cpu_data.data(), gpu_data, num_entries * sizeof(T), cudaMemcpyDeviceToHost), "cudaRecv(type, type)");
}
template<typename T> void AccelerationMeta::delGpuArray(AccelerationContext const*, T *x){
    TasGpu::cucheck( cudaFree(x), "cudaFree() in delCudaArray()");
}

void* AccelerationMeta::createCublasHandle(){
    cublasHandle_t handle;
    cublasCreate(&handle);
    return (void*) handle;
}
void AccelerationMeta::deleteCublasHandle(void *handle){
    cublasDestroy(reinterpret_cast<cublasHandle_t>(handle));
}

template void AccelerationMeta::recvGpuArray<double>(AccelerationContext const*, size_t num_entries, const double*, std::vector<double>&);
template void AccelerationMeta::recvGpuArray<float>(AccelerationContext const*, size_t num_entries, const float*, std::vector<float>&);
template void AccelerationMeta::recvGpuArray<int>(AccelerationContext const*, size_t num_entries, const int*, std::vector<int>&);

template void AccelerationMeta::delGpuArray<double>(AccelerationContext const*, double*);
template void AccelerationMeta::delGpuArray<float>(AccelerationContext const*, float*);
template void AccelerationMeta::delGpuArray<int>(AccelerationContext const*, int*);

namespace TasGpu{
/*
 * cuBLAS section
 */
//! \brief Converts character to cublas operation.
constexpr cublasOperation_t cublas_trans(char trans){
    return (trans == 'N') ? CUBLAS_OP_N : ((trans == 'T') ? CUBLAS_OP_T : CUBLAS_OP_C);
}

//! \brief Wrapper around sgeam().
void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
          int m, int n, float alpha, float const A[], int lda,
          float beta, float const B[], int ldb, float C[], int ldc){
    cucheck(cublasSgeam(handle, transa, transb, m, n, &alpha, A, lda, &beta, B, ldb, C, ldc), "cublasSgeam()");
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

//! \brief Wrapper around sgemv().
inline void gemv(cublasHandle_t handle, cublasOperation_t transa, int M, int N,
                 float alpha, float const A[], int lda, float const x[], int incx, float beta, float y[], int incy){
    cucheck( cublasSgemv(handle, transa, M, N, &alpha, A, lda, x, incx, &beta, y, incy), "cublasSgemv()");
}
//! \brief Wrapper around dgemv().
inline void gemv(cublasHandle_t handle, cublasOperation_t transa, int M, int N,
                 double alpha, double const A[], int lda, double const x[], int incx, double beta, double y[], int incy){
    cucheck( cublasDgemv(handle, transa, M, N, &alpha, A, lda, x, incx, &beta, y, incy), "cublasDgemv()");
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

//! \brief Wrapper around sgemm().
inline void gemm(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, int M, int N, int K,
                 float alpha, float const A[], int lda, float const B[], int ldb, float beta, float C[], int ldc){
    cucheck( cublasSgemm(handle, transa, transb, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc), "cublasSgemm()");
}
//! \brief Wrapper around dgemm().
inline void gemm(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, int M, int N, int K,
                 double alpha, double const A[], int lda, double const B[], int ldb, double beta, double C[], int ldc){
    cucheck( cublasDgemm(handle, transa, transb, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc), "cublasDgemm()");
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
 * cuSparse section
 */
//! \brief Returns the buffer size needed by gemvi.
template<typename T>
size_t size_sparse_gemvi(cusparseHandle_t handle, cusparseOperation_t transa, int M, int N, int nnz){
    static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value, "size_sparse_gemvi() works only with float and double");
    int buff_size = 0;
    if (std::is_same<T, float>::value){
        cucheck( cusparseSgemvi_bufferSize(handle, transa, M, N, nnz, &buff_size), "cusparseSgemvi_bufferSize()");
    }else{
        cucheck( cusparseDgemvi_bufferSize(handle, transa, M, N, nnz, &buff_size), "cusparseDgemvi_bufferSize()");
    }
    return static_cast<size_t>(buff_size);
}
//! \brief Wrapper around sgemvi().
inline void sparse_gemvi(cusparseHandle_t handle, cusparseOperation_t  transa,
               int M, int N, float alpha, float const A[], int lda,
               int nnz, const float x[], int const indx[], float beta, float y[]){
    GpuVector<float> buff(nullptr, size_sparse_gemvi<float>(handle, transa, M, N, nnz) );
    cucheck( cusparseSgemvi(handle, transa, M, N, &alpha, A, lda, nnz, x, indx, &beta, y, CUSPARSE_INDEX_BASE_ZERO, buff.data()), "cusparseSgemvi()");
}
//! \brief Wrapper around dgemvi().
inline void sparse_gemvi(cusparseHandle_t handle, cusparseOperation_t  transa,
               int M, int N, double alpha, double const A[], int lda,
               int nnz, const double x[], int const indx[], double beta, double y[]){
    GpuVector<double> buff(nullptr, size_sparse_gemvi<double>(handle, transa, M, N, nnz) );
    cucheck( cusparseDgemvi(handle, transa, M, N, &alpha, A, lda, nnz, x, indx, &beta, y, CUSPARSE_INDEX_BASE_ZERO, buff.data()), "cusparseDgemvi()");
}

#if (CUDART_VERSION < 11000)
//! \brief Wrapper around sgemv().
inline void sparse_gemv(cusparseHandle_t handle, cusparseOperation_t transa, int M, int N, int nnz,
                        float alpha, cusparseMatDesc &matdesc, float const vals[], int const pntr[], int const indx[],
                        float const x[], float beta, float y[]){
    cucheck( cusparseScsrmv(handle, transa, M, N, nnz, &alpha, matdesc, vals, pntr, indx, x, &beta, y), "cusparseScsrmv()");
}
//! \brief Wrapper around dgemv().
inline void sparse_gemv(cusparseHandle_t handle, cusparseOperation_t transa, int M, int N, int nnz,
                        double alpha, cusparseMatDesc &matdesc, double const vals[], int const pntr[], int const indx[],
                        double const x[], double beta, double y[]){
    cucheck( cusparseDcsrmv(handle, transa, M, N, nnz, &alpha, matdesc, vals, pntr, indx, x, &beta, y), "cusparseDcsrmv()");
}

//! \brief Wrapper around sgemm().
inline void sparse_gemm(cusparseHandle_t handle, cusparseOperation_t transa, cusparseOperation_t transb,
                        int M, int N, int K, int nnz, float alpha,
                        cusparseMatDesc &matdesc, float const vals[], int const pntr[], int const indx[],
                        float const B[], int ldb, float beta, float C[], int ldc){
    cucheck( cusparseScsrmm2(handle, transa, transb, M, N, K, nnz, &alpha, matdesc, vals, pntr, indx, B, ldb, &beta, C, ldc), "cusparseScsrmm2()");
}
//! \brief Wrapper around dgemm().
inline void sparse_gemm(cusparseHandle_t handle, cusparseOperation_t transa, cusparseOperation_t transb,
                        int M, int N, int K, int nnz, double alpha,
                        cusparseMatDesc &matdesc, double const vals[], int const pntr[], int const indx[],
                        double const B[], int ldb, double beta, double C[], int ldc){
    cucheck( cusparseDcsrmm2(handle, transa, transb, M, N, K, nnz, &alpha, matdesc, vals, pntr, indx, B, ldb, &beta, C, ldc), "cusparseDcsrmm2()");
}
#else
//! \brief Wrapper around sgemv().
template<typename scalar_type>
inline void sparse_gemv(cusparseHandle_t handle, int M, int N, int nnz, typename GpuVector<scalar_type>::value_type alpha,
                        scalar_type const vals[], int const pntr[], int const indx[], scalar_type const x[],
                        typename GpuVector<scalar_type>::value_type beta, scalar_type y[]){
    auto matdesc = makeSparseMatDesc(M, N, nnz, pntr, indx, vals);
    auto xdesc = makeSparseDenseVecDesc(N, x);
    auto ydesc = makeSparseDenseVecDesc(M, y);
    cudaDataType_t cuda_type = (std::is_same<scalar_type, float>::value) ? CUDA_R_32F : CUDA_R_64F;

    size_t bsize = 0;
    cucheck( cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matdesc, xdesc, &beta, ydesc, cuda_type, CUSPARSE_MV_ALG_DEFAULT, &bsize),
        "cusparseSpMV_bufferSize()"
    );

    GpuVector<scalar_type> buffer(nullptr, bsize / sizeof(scalar_type));

    cucheck( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matdesc, xdesc, &beta, ydesc, cuda_type, CUSPARSE_MV_ALG_DEFAULT, buffer.data()),
        "cusparseSpMV()"
    );
}
//! \brief Wrapper around dgemm().
template<typename scalar_type>
inline void sparse_gemm(cusparseHandle_t handle, int M, int N, int K, int nnz, typename GpuVector<scalar_type>::value_type alpha,
                        scalar_type const vals[], int const pntr[], int const indx[],
                        scalar_type const B[], int ldb, typename GpuVector<scalar_type>::value_type beta, scalar_type C[], int ldc){
    auto matdesc = makeSparseMatDesc(M, K, nnz, pntr, indx, vals);
    auto bdesc = makeSparseDenseMatDesc(N, K, ldb, B);
    auto cdesc = makeSparseDenseMatDesc(M, N, ldc, C);
    cudaDataType_t cuda_type = (std::is_same<scalar_type, float>::value) ? CUDA_R_32F : CUDA_R_64F;

    size_t bsize = 0;
    cucheck( cusparseSpMM_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_TRANSPOSE,
                                     &alpha, matdesc, bdesc, &beta, cdesc, cuda_type, CUSPARSE_SPMM_ALG_DEFAULT, &bsize), "cusparseSpMM_bufferSize()");

    GpuVector<scalar_type> buffer(nullptr, bsize / sizeof(scalar_type));

    cucheck( cusparseSpMM(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_TRANSPOSE, &alpha, matdesc, bdesc, &beta, cdesc,
                          cuda_type, CUSPARSE_SPMM_ALG_DEFAULT, buffer.data()),
        "cusparseSpMM()"
    );
}
#endif

/*
 * cuSolver section
 */
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
    GpuVector<double> workspace(nullptr, size_geqrf(handle, m, n, A, lda) );
    GpuVector<int> info(nullptr, std::vector<int>(1, 0));
    cucheck(cusolverDnDgeqrf(handle, m, n, A, lda, tau, workspace.data(), static_cast<int>(workspace.size()), info.data()), "cusolverDnDgeqrf()");
    if (info.unload(nullptr)[0] != 0)
        throw std::runtime_error("cusolverDnDgeqrf() returned non-zero status: " + std::to_string(info.unload(nullptr)[0]));
}
//! \brief Wrapper around zgeqrf().
inline void geqrf(cusolverDnHandle_t handle, int m, int n, std::complex<double> A[], int lda, std::complex<double> tau[]){
    GpuVector<std::complex<double>> workspace( nullptr, size_geqrf(handle, m, n, A, lda) );
    GpuVector<int> info(nullptr, std::vector<int>(1, 0));
    cucheck(cusolverDnZgeqrf(handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda,
                             reinterpret_cast<cuDoubleComplex*>(tau), reinterpret_cast<cuDoubleComplex*>(workspace.data()),
                             static_cast<int>(workspace.size()), info.data()), "cusolverDnZgeqrf()");
    if (info.unload(nullptr)[0] != 0)
        throw std::runtime_error("cusolverDnZgeqrf() returned non-zero status: " + std::to_string(info.unload(nullptr)[0]));
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
    GpuVector<int> info(nullptr, std::vector<int>(1, 0));
    GpuVector<double> workspace( nullptr, size_gemqr(handle, side, trans, m, n, k, A, lda, tau, C, ldc) );
    cucheck(cusolverDnDormqr(handle, side, trans, m, n, k, A, lda, tau, C, ldc,
                             workspace.data(), static_cast<int>(workspace.size()), info.data()), "cusolverDnDormqr()");
    if (info.unload(nullptr)[0] != 0)
        throw std::runtime_error("cusolverDnDormqr() returned non-zero status: " + std::to_string(info.unload(nullptr)[0]));
}
//! \brief Wrapper around zunmqr().
inline void gemqr(cusolverDnHandle_t handle, cublasSideMode_t side, cublasOperation_t trans,
                  int m, int n, int k, std::complex<double> const A[], int lda, std::complex<double> const tau[],
                  std::complex<double> C[], int ldc){
    GpuVector<int> info(nullptr, std::vector<int>(1, 0));
    GpuVector<std::complex<double>> workspace( nullptr, size_gemqr(handle, side, trans, m, n, k, A, lda, tau, C, ldc) );
    cucheck(cusolverDnZunmqr(handle, side, trans, m, n, k, reinterpret_cast<cuDoubleComplex const*>(A), lda,
                             reinterpret_cast<cuDoubleComplex const*>(tau), reinterpret_cast<cuDoubleComplex*>(C), ldc,
                             reinterpret_cast<cuDoubleComplex*>(workspace.data()), static_cast<int>(workspace.size()), info.data()),
            "cusolverDnZunmqr()");
    if (info.unload(nullptr)[0] != 0)
        throw std::runtime_error("cusolverDnZunmqr() returned non-zero status: " + std::to_string(info.unload(nullptr)[0]));
}

//! \brief Wrapper around cusolverDnDgetrs().
void getrs(cusolverDnHandle_t handle, cublasOperation_t trans, int n, int nrhs, double const A[], int lda, int const ipiv[], double B[], int ldb){
    GpuVector<int> info(nullptr, std::vector<int>(1, 0));
    cucheck(cusolverDnDgetrs(handle, trans, n, nrhs, A, lda, ipiv, B, ldb, info.data()), "cusolverDnDgetrs()");
    if (info.unload(nullptr)[0] != 0)
        throw std::runtime_error("cusolverDnDgetrs() returned non-zero status: " + std::to_string(info.unload(nullptr)[0]));
}

//! \brief Wrapper around cusolverDnDgetrf().
void factorizePLU(AccelerationContext const *acceleration, int n, double A[], int ipiv[]){
    cusolverDnHandle_t cusolverdnh = getCuSolverDnHandle(acceleration);
    int size = 0;
    cucheck(cusolverDnDgetrf_bufferSize(cusolverdnh, n, n, A, n, &size), "cusolverDnDgetrf_bufferSize()");
    GpuVector<double> workspace(nullptr, static_cast<size_t>(size));
    GpuVector<int> info(nullptr, std::vector<int>(1, 0));
    cucheck(cusolverDnDgetrf(cusolverdnh, n, n, A, n, workspace.data(), ipiv, info.data()), "cusolverDnDgetrf()");
    if (info.unload(nullptr)[0] != 0)
        throw std::runtime_error("cusolverDnDgetrf() returned non-zero status: " + std::to_string(info.unload(nullptr)[0]));
}

void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int const ipiv[], double b[]){
    cusolverDnHandle_t cusolverdnh = getCuSolverDnHandle(acceleration);
    getrs(cusolverdnh, (trans == 'T') ? CUBLAS_OP_T: CUBLAS_OP_N, n, 1, A, n, ipiv, b, n);
}
void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int const ipiv[], int nrhs, double B[]){
    cublasHandle_t cublash = getCuBlasHandle(acceleration);
    cusolverDnHandle_t cusolverdnh = getCuSolverDnHandle(acceleration);
    GpuVector<double> BT(nullptr, n, nrhs);
    geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, n, nrhs, 1.0, B, nrhs, 0.0, B, nrhs, BT.data(), n);
    getrs(cusolverdnh, (trans == 'T') ? CUBLAS_OP_T: CUBLAS_OP_N, n, nrhs, A, n, ipiv, BT.data(), n);
    geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, nrhs, n, 1.0, BT.data(), n, 0.0, BT.data(), n, B, nrhs);
}

/*
 * Algorithm section
 */
template<typename scalar_type>
void solveLSmultiGPU(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]){
    cublasHandle_t cublash = getCuBlasHandle(acceleration);
    cusolverDnHandle_t cusolverdnh = getCuSolverDnHandle(acceleration);

    GpuVector<scalar_type> AT(nullptr, n, m);
    geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, n, m, 1.0, A, m, 0.0, A, m, AT.data(), n);

    GpuVector<scalar_type> T(nullptr, m);
    geqrf(cusolverdnh, n, m, AT.data(), n, T.data());

    cublasOperation_t trans = (std::is_same<scalar_type, double>::value) ? CUBLAS_OP_T : CUBLAS_OP_C;

    if (nrhs == 1){
        gemqr(cusolverdnh, CUBLAS_SIDE_LEFT, trans, n, 1, m, AT.data(), n, T.data(), B, n);
        trsv(cublash, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, m, AT.data(), n, B, 1);
    }else{
        GpuVector<scalar_type> BT(nullptr, n, nrhs);
        geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, n, nrhs, 1.0, B, nrhs, 0.0, B, nrhs, BT.data(), n);

        gemqr(cusolverdnh, CUBLAS_SIDE_LEFT, trans, n, nrhs, m, AT.data(), n, T.data(), BT.data(), n);
        trsm(cublash, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, m, nrhs, 1.0, AT.data(), n, BT.data(), n);

        geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, nrhs, n, 1.0, BT.data(), n, 0.0, BT.data(), n, B, nrhs);
    }
}

template void solveLSmultiGPU<double>(AccelerationContext const*, int, int, double[], int, double[]);
template void solveLSmultiGPU<std::complex<double>>(AccelerationContext const*, int, int, std::complex<double>[], int, std::complex<double>[]);

#ifndef Tasmanian_ENABLE_MAGMA
template<typename scalar_type>
void solveLSmultiOOC(AccelerationContext const*, int, int, scalar_type[], int, scalar_type[]){}
#endif

template void solveLSmultiOOC<double>(AccelerationContext const*, int, int, double[], int, double[]);
template void solveLSmultiOOC<std::complex<double>>(AccelerationContext const*, int, int, std::complex<double>[], int, std::complex<double>[]);


template<typename scalar_type>
void denseMultiply(AccelerationContext const *acceleration, int M, int N, int K, typename GpuVector<scalar_type>::value_type alpha, GpuVector<scalar_type> const &A,
                   GpuVector<scalar_type> const &B, typename GpuVector<scalar_type>::value_type beta, scalar_type C[]){
    cublasHandle_t cublash = getCuBlasHandle(acceleration);
    if (M > 1){
        if (N > 1){ // matrix-matrix mode
            gemm(cublash, CUBLAS_OP_N, CUBLAS_OP_N, M, N, K, alpha, A.data(), M, B.data(), K, beta, C, M);
        }else{ // matrix vector, A * v = C
            gemv(cublash, CUBLAS_OP_N, M, K, alpha, A.data(), M, B.data(), 1, beta, C, 1);
        }
    }else{ // matrix vector B^T * v = C
        gemv(cublash, CUBLAS_OP_T, K, N, alpha, B.data(), K, A.data(), 1, beta, C, 1);
    }
}

template void denseMultiply<float>(AccelerationContext const*, int, int, int, float,
                                   GpuVector<float> const&, GpuVector<float> const&, float, float[]);
template void denseMultiply<double>(AccelerationContext const*, int, int, int, double,
                                    GpuVector<double> const&, GpuVector<double> const&, double, double[]);

#if (CUDART_VERSION < 11000)
template<typename scalar_type>
void sparseMultiply(AccelerationContext const *acceleration, int M, int N, int K, typename GpuVector<scalar_type>::value_type alpha,
                    GpuVector<scalar_type> const &A, GpuVector<int> const &pntr, GpuVector<int> const &indx,
                    GpuVector<scalar_type> const &vals, scalar_type C[]){

    cusparseHandle_t cusparseh = getCuSparseHandle(acceleration);

    if (N > 1){ // dense matrix has many columns
        cusparseMatDesc matdesc;
        if (M > 1){ // dense matrix has many rows, use matrix-matrix algorithm
            GpuVector<scalar_type> tempC(nullptr, M, N);
            sparse_gemm(cusparseh, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_TRANSPOSE, N, M, K, (int) indx.size(),
                        alpha, matdesc, vals.data(), pntr.data(), indx.data(), A.data(), M, 0.0, tempC.data(), N);

            cublasHandle_t cublash = getCuBlasHandle(acceleration);
            geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, M, N, 1.0, tempC.data(), N, 0.0, tempC.data(), N, C, M);
        }else{ // dense matrix has only one row, use sparse matrix times dense vector
            sparse_gemv(cusparseh, CUSPARSE_OPERATION_NON_TRANSPOSE, N, K, (int) indx.size(),
                        alpha, matdesc, vals.data(), pntr.data(), indx.data(), A.data(), 0.0, C);
        }
    }else{ // sparse matrix has only one column, use dense matrix times sparse vector
        // quote from Nvidia CUDA cusparse manual at https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-gemvi
        // "This function requires no extra storage for the general matrices when operation CUSPARSE_OPERATION_NON_TRANSPOSE is selected."
        // Yet, buffer is required when num_nz exceeds 32 even with CUSPARSE_OPERATION_NON_TRANSPOSE
        sparse_gemvi(cusparseh, CUSPARSE_OPERATION_NON_TRANSPOSE, M, K, alpha, A.data(), M, (int) indx.size(), vals.data(), indx.data(), 0.0, C);
    }
}
#else
template<typename scalar_type>
void sparseMultiply(AccelerationContext const *acceleration, int M, int N, int K, typename GpuVector<scalar_type>::value_type alpha,
                    GpuVector<scalar_type> const &A, GpuVector<int> const &pntr, GpuVector<int> const &indx,
                    GpuVector<scalar_type> const &vals, scalar_type C[]){

    cusparseHandle_t cusparseh = getCuSparseHandle(acceleration);

    if (N > 1){ // dense matrix has many columns
        if (M > 1){ // dense matrix has many rows, use matrix-matrix algorithm
            GpuVector<scalar_type> tempC(nullptr, M, N);
            sparse_gemm(cusparseh, N, M, K, (int) indx.size(), alpha, vals.data(), pntr.data(), indx.data(), A.data(), M, 0.0, tempC.data(), N);

            cublasHandle_t cublash = getCuBlasHandle(acceleration);
            geam(cublash, CUBLAS_OP_T, CUBLAS_OP_T, M, N, 1.0, tempC.data(), N, 0.0, tempC.data(), N, C, M);
        }else{ // dense matrix has only one row, use sparse matrix times dense vector
            sparse_gemv(cusparseh, N, K, (int) indx.size(), alpha, vals.data(), pntr.data(), indx.data(), A.data(), 0.0, C);
        }
    }else{
        sparse_gemvi(cusparseh, CUSPARSE_OPERATION_NON_TRANSPOSE, M, K, alpha, A.data(), M, (int) indx.size(), vals.data(), indx.data(), 0.0, C);
    }
}
#endif

template void sparseMultiply<float>(AccelerationContext const*, int, int, int, float, GpuVector<float> const &A,
                                    GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<float> const &vals, float C[]);
template void sparseMultiply<double>(AccelerationContext const*, int, int, int, double, GpuVector<double> const &A,
                                     GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<double> const &vals, double C[]);

template<typename T> void load_n(AccelerationContext const*, T const *cpu_data, size_t num_entries, T *gpu_data){
    TasGpu::cucheck( cudaMemcpy(gpu_data, cpu_data, num_entries * sizeof(T), cudaMemcpyHostToDevice), "cudaMemcpy() load_n to device");
}

template void load_n<int>(AccelerationContext const*, int const*, size_t, int*);
template void load_n<float>(AccelerationContext const*, float const*, size_t, float*);
template void load_n<double>(AccelerationContext const*, double const*, size_t, double*);
template void load_n<std::complex<double>>(AccelerationContext const*, std::complex<double> const*, size_t, std::complex<double>*);

}
}

#endif
