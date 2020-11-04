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

#ifndef __TASMANIAN_DPCPP_WRAPPERS_CPP
#define __TASMANIAN_DPCPP_WRAPPERS_CPP

#include "tsgDpcppWrappers.hpp"

/*!
 * \file tsgDpcppWrappers.cpp
 * \brief Wrappers to DPC++ functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * Realizations of the GPU algorithms using the DPC++/MKL backend.
 */

namespace TasGrid{

template<typename T> void GpuVector<T>::resize(AccelerationContext const *acc, size_t count){
//     if (count != num_entries){ // if the current array is not big enough
//         clear(); // resets dynamic_mode
//         num_entries = count;
//         TasGpu::hipcheck( hipMalloc((void**) &gpu_data, num_entries * sizeof(T)), "hipMalloc()");
//     }
}
template<typename T> void GpuVector<T>::clear(){
//     num_entries = 0;
//     if (gpu_data != nullptr) // if I own the data and the data is not null
//         TasGpu::hipcheck( hipFree(gpu_data), "hipFree()");
//     gpu_data = nullptr;
}
template<typename T> void GpuVector<T>::load(AccelerationContext const *acc, size_t count, const T* cpu_data){
//     resize(count);
//     TasGpu::hipcheck( hipMemcpy(gpu_data, cpu_data, num_entries * sizeof(T), hipMemcpyHostToDevice), "hipMemcpy() to device");
}
template<typename T> void GpuVector<T>::unload(AccelerationContext const *acc, size_t num, T* cpu_data) const{
//    TasGpu::hipcheck( hipMemcpy(cpu_data, gpu_data, num * sizeof(T), hipMemcpyDeviceToHost), "hipMemcpy() from device");
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

GpuEngine::~GpuEngine(){
}

int AccelerationMeta::getNumGpuDevices(){
    return 0;
}
void AccelerationMeta::setDefaultCudaDevice(int deviceID){
}
unsigned long long AccelerationMeta::getTotalGPUMemory(int deviceID){ // int deviceID
    return 0;
}
std::string AccelerationMeta::getCudaDeviceName(int deviceID){ // int deviceID
    return "";
}
template<typename T> void AccelerationMeta::recvCudaArray(size_t num_entries, const T *gpu_data, std::vector<T> &cpu_data){

}
template<typename T> void AccelerationMeta::delCudaArray(T *x){

}

template void AccelerationMeta::recvCudaArray<double>(size_t num_entries, const double*, std::vector<double>&);
template void AccelerationMeta::recvCudaArray<float>(size_t num_entries, const float*, std::vector<float>&);
template void AccelerationMeta::recvCudaArray<int>(size_t num_entries, const int*, std::vector<int>&);

template void AccelerationMeta::delCudaArray<double>(double*);
template void AccelerationMeta::delCudaArray<float>(float*);
template void AccelerationMeta::delCudaArray<int>(int*);

namespace TasGpu{

//! \brief Wrapper around sgeam().
void geam(//rocblas_handle handle, rocblas_operation transa, rocblas_operation transb,
          int m, int n, float alpha, float const A[], int lda,
          float beta, float const B[], int ldb, float C[], int ldc){

}
//! \brief Wrapper around dgeam().
void geam(//rocblas_handle handle, rocblas_operation transa, rocblas_operation transb,
          int m, int n, double alpha, double const A[], int lda,
          double beta, double const B[], int ldb, double C[], int ldc){

}
//! \brief Wrapper around sgemv().
inline void gemv(//rocblas_handle handle, rocblas_operation transa,
                 int M, int N,
                 float alpha, float const A[], int lda, float const x[], int incx, float beta, float y[], int incy){

}
//! \brief Wrapper around dgemv().
inline void gemv(//rocblas_handle handle, rocblas_operation transa,
                 int M, int N,
                 double alpha, double const A[], int lda, double const x[], int incx, double beta, double y[], int incy){

}

//! \brief Wrapper around sgemm().
inline void gemm(//rocblas_handle handle, rocblas_operation transa, rocblas_operation transb,
                 int M, int N, int K,
                 float alpha, float const A[], int lda, float const B[], int ldb, float beta, float C[], int ldc){

}
//! \brief Wrapper around dgemm().
inline void gemm(//rocblas_handle handle, rocblas_operation transa, rocblas_operation transb,
                 int M, int N, int K,
                 double alpha, double const A[], int lda, double const B[], int ldb, double beta, double C[], int ldc){

}
//! \brief Wrapper around dtrsm().
inline void trsm(//rocblas_handle handle, rocblas_side side, rocblas_fill uplo,
                 //rocblas_operation trans, rocblas_diagonal diag,
                 int m, int n,
                 double alpha, double const A[], int lda, double B[], int ldb){

}
//! \brief Wrapper around ztrsm().
inline void trsm(//rocblas_handle handle, rocblas_side side, rocblas_fill uplo,
                 //rocblas_operation trans, rocblas_diagonal diag,
                 int m, int n,
                 std::complex<double> alpha, std::complex<double> const A[], int lda, std::complex<double> B[], int ldb){

}

inline void sparse_gemv(//rocsparse_handle handle, rocsparse_operation trans, int M, int N, int nnz,
                        float alpha, float const vals[], int const pntr[], int const indx[],
                        float const x[], float beta, float y[]){

}
inline void sparse_gemv(//rocsparse_handle handle, rocsparse_operation trans, int M, int N, int nnz,
                        double alpha, double const vals[], int const pntr[], int const indx[],
                        double const x[], double beta, double y[]){

}

inline void sparse_gemm(//rocsparse_handle handle, rocsparse_operation transa, rocsparse_operation transb,
                        int M, int N, int K, int nnz, float alpha,
                        float const vals[], int const pntr[], int const indx[],
                        float const B[], int ldb, float beta, float C[], int ldc){

}

inline void sparse_gemm(//rocsparse_handle handle, rocsparse_operation transa, rocsparse_operation transb,
                        int M, int N, int K, int nnz, double alpha,
                        double const vals[], int const pntr[], int const indx[],
                        double const B[], int ldb, double beta, double C[], int ldc){

}


//! \brief Wrapper around rocsolver_dgetrs().
void getrs(//rocblas_handle handle, rocblas_operation trans,
           int n, int nrhs, double const A[], int lda, int const ipiv[], double B[], int ldb){

}

//! \brief Wrapper around rocsolver_dgetrf().
void factorizePLU(AccelerationContext const *acceleration, int n, double A[], int ipiv[]){

}

void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int const ipiv[], double b[]){

}
void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int const ipiv[], int nrhs, double B[]){

}

//! \brief Wrapper around rocsolver_dgelqf().
inline void gelqf(//rocblas_handle handle,
                  int m, int n, double A[], double tau[]){

}
//! \brief Wrapper around rocsolver_zgelqf().
inline void gelqf(//rocblas_handle handle,
                  int m, int n, std::complex<double> A[], std::complex<double> tau[]){

}

//! \brief Wrapper around rocsolver_dormlq(), does Q^T times C.
inline void gemlq(//rocblas_handle handle,
                  int m, int n, int k, double A[], double tau[], double C[]){

}
//! \brief Wrapper around rocsolver_dunmlq(), does Q^T times C.
inline void gemlq(//rocblas_handle handle,
                  int m, int n, int k, std::complex<double> A[], std::complex<double> tau[], std::complex<double> C[]){

}

/*
 * Algorithm section
 */
template<typename scalar_type>
void solveLSmultiGPU(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]){

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
    if (M > 1){
        if (N > 1){ // matrix-matrix mode
            //gemm(cublash, rocblas_operation_none, rocblas_operation_none, M, N, K, alpha, A.data(), M, B.data(), K, beta, C, M);
        }else{ // matrix vector, A * v = C
            //gemv(cublash, rocblas_operation_none, M, K, alpha, A.data(), M, B.data(), 1, beta, C, 1);
        }
    }else{ // matrix vector B^T * v = C
        //gemv(cublash, rocblas_operation_transpose, K, N, alpha, B.data(), K, A.data(), 1, beta, C, 1);
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

}

template void sparseMultiply<float>(AccelerationContext const*, int, int, int, float, GpuVector<float> const &A,
                                    GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<float> const &vals, float C[]);
template void sparseMultiply<double>(AccelerationContext const*, int, int, int, double, GpuVector<double> const &A,
                                     GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<double> const &vals, double C[]);

template<typename T> void load_n(T const *cpu_data, size_t num_entries, T *gpu_data){

}

template void load_n<int>(int const*, size_t, int*);
template void load_n<float>(float const*, size_t, float*);
template void load_n<double>(double const*, size_t, double*);
template void load_n<std::complex<double>>(std::complex<double> const*, size_t, std::complex<double>*);

}
}

#endif
