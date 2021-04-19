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
#include "tsgHipWrappers.hpp"

/*!
 * \file tsgHipWrappers.cpp
 * \brief Wrappers to HIP functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * Realizations of the GPU algorithms using the HIP backend.
 */

namespace TasGrid{
/*
 * Meta methods
 */
template<typename T> void GpuVector<T>::resize(AccelerationContext const*, size_t count){
    if (count != num_entries){ // if the current array is not big enough
        clear(); // resets dynamic_mode
        num_entries = count;
        TasGpu::hipcheck( hipMalloc((void**) &gpu_data, num_entries * sizeof(T)), "hipMalloc()");
    }
}
template<typename T> void GpuVector<T>::clear(){
    num_entries = 0;
    if (gpu_data != nullptr) // if I own the data and the data is not null
        TasGpu::hipcheck( hipFree(gpu_data), "hipFree()");
    gpu_data = nullptr;
}
template<typename T> void GpuVector<T>::load(AccelerationContext const*, size_t count, const T* cpu_data){
    resize(nullptr, count);
    TasGpu::hipcheck( hipMemcpy(gpu_data, cpu_data, num_entries * sizeof(T), hipMemcpyHostToDevice), "hipMemcpy() to device");
}
template<typename T> void GpuVector<T>::unload(AccelerationContext const*, size_t num, T* cpu_data) const{
    TasGpu::hipcheck( hipMemcpy(cpu_data, gpu_data, num * sizeof(T), hipMemcpyDeviceToHost), "hipMemcpy() from device");
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

template<> void deleteHandle<AccHandle::Rocblas>(int *p){ rocblas_destroy_handle(reinterpret_cast<rocblas_handle>(p)); }
template<> void deleteHandle<AccHandle::Rocsparse>(int *p){ rocsparse_destroy_handle(reinterpret_cast<rocsparse_handle>(p)); }

void GpuEngine::setRocBlasHandle(void *handle){
    rblas_handle = std::unique_ptr<int, HandleDeleter<AccHandle::Rocblas>>
        (reinterpret_cast<int*>(handle), HandleDeleter<AccHandle::Rocblas>(false));
}
void GpuEngine::setRocSparseHandle(void *handle){
    rsparse_handle = std::unique_ptr<int, HandleDeleter<AccHandle::Rocsparse>>
        (reinterpret_cast<int*>(handle), HandleDeleter<AccHandle::Rocsparse>(false));
}

int AccelerationMeta::getNumGpuDevices(){
    int gpu_count = 0;
    hipGetDeviceCount(&gpu_count);
    return gpu_count;
}
void AccelerationMeta::setDefaultGpuDevice(int deviceID){
    hipSetDevice(deviceID);
}
unsigned long long AccelerationMeta::getTotalGPUMemory(int deviceID){ // int deviceID
    hipDeviceProp_t prop;
    hipGetDeviceProperties(&prop, deviceID);
    return prop.totalGlobalMem;
}
std::string AccelerationMeta::getGpuDeviceName(int deviceID){ // int deviceID
    if ((deviceID < 0) || (deviceID >= getNumGpuDevices())) return std::string();

    hipDeviceProp_t prop;
    hipGetDeviceProperties(&prop, deviceID);

    return std::string(prop.name);
}
template<typename T> void AccelerationMeta::recvGpuArray(AccelerationContext const*, size_t num_entries, const T *gpu_data, std::vector<T> &cpu_data){
    cpu_data.resize(num_entries);
    TasGpu::hipcheck( hipMemcpy(cpu_data.data(), gpu_data, num_entries * sizeof(T), hipMemcpyDeviceToHost), "hip receive");
}
template<typename T> void AccelerationMeta::delGpuArray(AccelerationContext const*, T *x){
    TasGpu::hipcheck( hipFree(x), "hipFree()");
}

template void AccelerationMeta::recvGpuArray<double>(AccelerationContext const*, size_t num_entries, const double*, std::vector<double>&);
template void AccelerationMeta::recvGpuArray<float>(AccelerationContext const*, size_t num_entries, const float*, std::vector<float>&);
template void AccelerationMeta::recvGpuArray<int>(AccelerationContext const*, size_t num_entries, const int*, std::vector<int>&);

template void AccelerationMeta::delGpuArray<double>(AccelerationContext const*, double*);
template void AccelerationMeta::delGpuArray<float>(AccelerationContext const*, float*);
template void AccelerationMeta::delGpuArray<int>(AccelerationContext const*, int*);

namespace TasGpu{
/*
 * rocBLAS section
 */
//! \brief Converts character to cublas operation.
constexpr rocblas_operation cublas_trans(char trans){
    return (trans == 'N') ? rocblas_operation_none : ((trans == 'T') ? rocblas_operation_transpose : rocblas_operation_conjugate_transpose);
}

//! \brief Wrapper around sgeam().
void geam(rocblas_handle handle, rocblas_operation transa, rocblas_operation transb,
          int m, int n, float alpha, float const A[], int lda,
          float beta, float const B[], int ldb, float C[], int ldc){
    hipcheck(rocblas_sgeam(handle, transa, transb, m, n, &alpha, A, lda, &beta, B, ldb, C, ldc), "cublasSgeam()");
}
//! \brief Wrapper around dgeam().
void geam(rocblas_handle handle, rocblas_operation transa, rocblas_operation transb,
          int m, int n, double alpha, double const A[], int lda,
          double beta, double const B[], int ldb, double C[], int ldc){
    hipcheck(rocblas_dgeam(handle, transa, transb, m, n, &alpha, A, lda, &beta, B, ldb, C, ldc), "cublasDgeam()");
}
//! \brief Wrapper around sgemv().
inline void gemv(rocblas_handle handle, rocblas_operation transa, int M, int N,
                 float alpha, float const A[], int lda, float const x[], int incx, float beta, float y[], int incy){
    hipcheck( rocblas_sgemv(handle, transa, M, N, &alpha, A, lda, x, incx, &beta, y, incy), "rocblas_sgemv()");
}
//! \brief Wrapper around dgemv().
inline void gemv(rocblas_handle handle, rocblas_operation transa, int M, int N,
                 double alpha, double const A[], int lda, double const x[], int incx, double beta, double y[], int incy){
    hipcheck( rocblas_dgemv(handle, transa, M, N, &alpha, A, lda, x, incx, &beta, y, incy), "rocblas_sgemv()");
}

//! \brief Wrapper around sgemm().
inline void gemm(rocblas_handle handle, rocblas_operation transa, rocblas_operation transb, int M, int N, int K,
                 float alpha, float const A[], int lda, float const B[], int ldb, float beta, float C[], int ldc){
    hipcheck( rocblas_sgemm(handle, transa, transb, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc), "rocblas_sgemm()");
}
//! \brief Wrapper around dgemm().
inline void gemm(rocblas_handle handle, rocblas_operation transa, rocblas_operation transb, int M, int N, int K,
                 double alpha, double const A[], int lda, double const B[], int ldb, double beta, double C[], int ldc){
    hipcheck( rocblas_dgemm(handle, transa, transb, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc), "rocblas_dgemm()");
}
//! \brief Wrapper around dtrsm().
inline void trsm(rocblas_handle handle, rocblas_side side, rocblas_fill uplo,
                 rocblas_operation trans, rocblas_diagonal diag, int m, int n,
                 double alpha, double const A[], int lda, double B[], int ldb){
    hipcheck(rocblas_dtrsm(handle, side, uplo, trans, diag, m, n, &alpha, A, lda, B, ldb), "rocblas_dtrsm()");
}
//! \brief Wrapper around ztrsm().
inline void trsm(rocblas_handle handle, rocblas_side side, rocblas_fill uplo,
                 rocblas_operation trans, rocblas_diagonal diag, int m, int n,
                 std::complex<double> alpha, std::complex<double> const A[], int lda, std::complex<double> B[], int ldb){
    hipcheck(rocblas_ztrsm(handle, side, uplo, trans, diag, m, n, reinterpret_cast<rocblas_double_complex*>(&alpha),
                           reinterpret_cast<rocblas_double_complex const*>(A), lda,
                           reinterpret_cast<rocblas_double_complex*>(B), ldb), "rocblas_ztrsm()");
}

/*
 * rocSparse section
 */
inline void sparse_gemv(rocsparse_handle handle, rocsparse_operation trans, int M, int N, int nnz,
                        float alpha, float const vals[], int const pntr[], int const indx[],
                        float const x[], float beta, float y[]){
    rocsparseMatDesc desc;
    rocsparseMatInfo info;
    hipcheck( rocsparse_scsrmv_analysis(handle, trans, M, N, nnz, desc, vals, pntr, indx, info), "sgemv-info");
    hipcheck( rocsparse_scsrmv(handle, trans, M, N, nnz, &alpha, desc, vals, pntr, indx, info, x, &beta, y), "sgemv");
}
inline void sparse_gemv(rocsparse_handle handle, rocsparse_operation trans, int M, int N, int nnz,
                        double alpha, double const vals[], int const pntr[], int const indx[],
                        double const x[], double beta, double y[]){
    rocsparseMatDesc desc;
    rocsparseMatInfo info;
    hipcheck( rocsparse_dcsrmv_analysis(handle, trans, M, N, nnz, desc, vals, pntr, indx, info), "dgemv-info");
    hipcheck( rocsparse_dcsrmv(handle, trans, M, N, nnz, &alpha, desc, vals, pntr, indx, info, x, &beta, y), "dgemv");
}
//! \brief Wrapper around dgemm().
inline void sparse_gemm(rocsparse_handle handle, rocsparse_operation transa, rocsparse_operation transb,
                        int M, int N, int K, int nnz, float alpha,
                        float const vals[], int const pntr[], int const indx[],
                        float const B[], int ldb, float beta, float C[], int ldc){
    rocsparseMatDesc desc;
    hipcheck( rocsparse_scsrmm(handle, transa, transb, M, N, K, nnz, &alpha, desc, vals, pntr, indx, B, ldb, &beta, C, ldc), "dgemm()");
}
//! \brief Wrapper around rocsparse_dcsrmm().
inline void sparse_gemm(rocsparse_handle handle, rocsparse_operation transa, rocsparse_operation transb,
                        int M, int N, int K, int nnz, double alpha,
                        double const vals[], int const pntr[], int const indx[],
                        double const B[], int ldb, double beta, double C[], int ldc){
    rocsparseMatDesc desc;
    hipcheck( rocsparse_dcsrmm(handle, transa, transb, M, N, K, nnz, &alpha, desc, vals, pntr, indx, B, ldb, &beta, C, ldc), "dgemm()");
}


/*
 * rocSolver section
 */
//! \brief Wrapper around rocsolver_dgetrs().
void getrs(rocblas_handle handle, rocblas_operation trans, int n, int nrhs, double const A[], int lda, int const ipiv[], double B[], int ldb){
    hipcheck(rocblas_set_pointer_mode(handle, rocblas_pointer_mode_device), "rocblas_set_pointer_mode()");
    hipcheck(rocsolver_dgetrs(handle, trans, n, nrhs, const_cast<double*>(A), lda, ipiv, B, ldb), "rocsolver_dgetrs()");
    rocblas_set_pointer_mode(handle, rocblas_pointer_mode_host);
}

//! \brief Wrapper around rocsolver_dgetrf().
void factorizePLU(AccelerationContext const *acceleration, int n, double A[], int ipiv[]){
    rocblas_handle rochandle = getRocBlasHandle(acceleration);
    GpuVector<int> info(acceleration, std::vector<int>(4, 0));
    hipcheck(rocblas_set_pointer_mode(rochandle, rocblas_pointer_mode_device), "rocblas_set_pointer_mode()");
    hipcheck(rocsolver_dgetrf(rochandle, n, n, A, n, ipiv, info.data()), "rocsolver_dgetrf()");
    rocblas_set_pointer_mode(rochandle, rocblas_pointer_mode_host);
    if (info.unload(nullptr)[0] != 0)
        throw std::runtime_error("rocsolver_dgetrf() returned non-zero status: " + std::to_string(info.unload(nullptr)[0]));
}

void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int const ipiv[], double b[]){
    rocblas_handle rochandle = getRocBlasHandle(acceleration);
    getrs(rochandle, (trans == 'T') ? rocblas_operation_transpose: rocblas_operation_none, n, 1, A, n, ipiv, b, n);
}
void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int const ipiv[], int nrhs, double B[]){
    rocblas_handle rochandle = getRocBlasHandle(acceleration);
    GpuVector<double> BT(nullptr, n, nrhs);
    geam(rochandle, rocblas_operation_transpose, rocblas_operation_transpose, n, nrhs, 1.0, B, nrhs, 0.0, B, nrhs, BT.data(), n);
    getrs(rochandle, (trans == 'T') ? rocblas_operation_transpose: rocblas_operation_none, n, nrhs, A, n, ipiv, BT.data(), n);
    geam(rochandle, rocblas_operation_transpose, rocblas_operation_transpose, nrhs, n, 1.0, BT.data(), n, 0.0, BT.data(), n, B, nrhs);
}

//! \brief Wrapper around rocsolver_dgelqf().
inline void gelqf(rocblas_handle handle, int m, int n, double A[], double tau[]){
    hipcheck(rocsolver_dgelqf(handle, m, n, A, m, tau), "rocsolver_dgelqf()");
}
//! \brief Wrapper around rocsolver_zgelqf().
inline void gelqf(rocblas_handle handle, int m, int n, std::complex<double> A[], std::complex<double> tau[]){
    hipcheck(rocsolver_zgelqf(handle, m, n, reinterpret_cast<rocblas_double_complex*>(A), m, reinterpret_cast<rocblas_double_complex*>(tau)), "rocsolver_zgelqf()");
}

//! \brief Wrapper around rocsolver_dormlq(), does Q^T times C.
inline void gemlq(rocblas_handle handle, int m, int n, int k, double A[], double tau[], double C[]){
    hipcheck(rocsolver_dormlq(handle, rocblas_side_right, rocblas_operation_transpose, m, n, k, A, k, tau, C, m), "rocsolver_dormlq()");
}
//! \brief Wrapper around rocsolver_dunmlq(), does Q^T times C.
inline void gemlq(rocblas_handle handle, int m, int n, int k, std::complex<double> A[], std::complex<double> tau[], std::complex<double> C[]){
    hipcheck(rocsolver_zunmlq(handle, rocblas_side_right, rocblas_operation_conjugate_transpose, m, n, k,
                              reinterpret_cast<rocblas_double_complex*>(A), k,
                              reinterpret_cast<rocblas_double_complex*>(tau), reinterpret_cast<rocblas_double_complex*>(C), m),
             "rocsolver_zunmlq()");
}

/*
 * Algorithm section
 */
template<typename scalar_type>
void solveLSmultiGPU(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]){
    rocblas_handle rochandle = getRocBlasHandle(acceleration);
    GpuVector<scalar_type> tau(acceleration, std::min(n, m));
    gelqf(rochandle, m, n, A, tau.data());
    gemlq(rochandle, nrhs, n, m, A, tau.data(), B);
    trsm(rochandle, rocblas_side_right, rocblas_fill_lower, rocblas_operation_none, rocblas_diagonal_non_unit, nrhs, m, 1.0, A, m, B, nrhs);
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
    rocblas_handle cublash = getRocBlasHandle(acceleration);
    if (M > 1){
        if (N > 1){ // matrix-matrix mode
            gemm(cublash, rocblas_operation_none, rocblas_operation_none, M, N, K, alpha, A.data(), M, B.data(), K, beta, C, M);
        }else{ // matrix vector, A * v = C
            gemv(cublash, rocblas_operation_none, M, K, alpha, A.data(), M, B.data(), 1, beta, C, 1);
        }
    }else{ // matrix vector B^T * v = C
        gemv(cublash, rocblas_operation_transpose, K, N, alpha, B.data(), K, A.data(), 1, beta, C, 1);
    }
}

template void denseMultiply<float>(AccelerationContext const*, int, int, int, float,
                                   GpuVector<float> const&, GpuVector<float> const&, float, float[]);
template void denseMultiply<double>(AccelerationContext const*, int, int, int, double,
                                    GpuVector<double> const&, GpuVector<double> const&, double, double[]);

template<typename scalar_type>
void sparseMultiply(AccelerationContext const *acceleration, int M, int N, int K, typename GpuVector<scalar_type>::value_type alpha,
                    GpuVector<scalar_type> const &A, GpuVector<int> const &pntr, GpuVector<int> const &indx,
                    GpuVector<scalar_type> const &vals, scalar_type C[]){

    rocsparse_handle rocsparseh = getRocSparseHandle(acceleration);

    if (N > 1){
        if (M > 1){
            GpuVector<scalar_type> tempC(nullptr, M, N);
            sparse_gemm(rocsparseh, rocsparse_operation_none, rocsparse_operation_transpose, N, M, K, static_cast<int>(indx.size()),
                        alpha, vals.data(), pntr.data(), indx.data(), A.data(), M, 0.0, tempC.data(), N);

            rocblas_handle rocblash = getRocBlasHandle(acceleration);
            geam(rocblash, rocblas_operation_transpose, rocblas_operation_transpose, M, N, 1.0, tempC.data(), N, 0.0, tempC.data(), N, C, M);
        }else{
            sparse_gemv(rocsparseh, rocsparse_operation_none, N, K, static_cast<int>(indx.size()),
                        alpha, vals.data(), pntr.data(), indx.data(), A.data(), 0.0, C);
        }
    }else{
        GpuVector<scalar_type> tempC(nullptr, M, N);
        int nnz = static_cast<int>(indx.size());
        GpuVector<int> temp_pntr(nullptr, std::vector<int>{0, nnz});
        sparse_gemm(rocsparseh, rocsparse_operation_none, rocsparse_operation_transpose, N, M, K, nnz,
                    alpha, vals.data(), temp_pntr.data(), indx.data(), A.data(), M, 0.0, tempC.data(), N);

        rocblas_handle rocblash = getRocBlasHandle(acceleration);
        geam(rocblash, rocblas_operation_transpose, rocblas_operation_transpose, M, N, 1.0, tempC.data(), N, 0.0, tempC.data(), N, C, M);
    }
}

template void sparseMultiply<float>(AccelerationContext const*, int, int, int, float, GpuVector<float> const &A,
                                    GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<float> const &vals, float C[]);
template void sparseMultiply<double>(AccelerationContext const*, int, int, int, double, GpuVector<double> const &A,
                                     GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<double> const &vals, double C[]);

template<typename T> void load_n(AccelerationContext const*, T const *cpu_data, size_t num_entries, T *gpu_data){
    TasGpu::hipcheck( hipMemcpy(gpu_data, cpu_data, num_entries * sizeof(T), hipMemcpyHostToDevice), "hipMemcpy() load_n to device");
}

template void load_n<int>(AccelerationContext const*, int const*, size_t, int*);
template void load_n<float>(AccelerationContext const*, float const*, size_t, float*);
template void load_n<double>(AccelerationContext const*, double const*, size_t, double*);
template void load_n<std::complex<double>>(AccelerationContext const*, std::complex<double> const*, size_t, std::complex<double>*);

}
}

#endif
