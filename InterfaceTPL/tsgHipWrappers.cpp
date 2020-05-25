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
template<typename T> void GpuVector<T>::resize(size_t count){
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
template<typename T> void GpuVector<T>::load(size_t count, const T* cpu_data){
    resize(count);
    TasGpu::hipcheck( hipMemcpy(gpu_data, cpu_data, num_entries * sizeof(T), hipMemcpyHostToDevice), "hipMemcpy() to device");
}
template<typename T> void GpuVector<T>::unload(size_t num, T* cpu_data) const{
    TasGpu::hipcheck( hipMemcpy(cpu_data, gpu_data, num * sizeof(T), hipMemcpyDeviceToHost), "hipMemcpy() from device");
}

template void GpuVector<double>::resize(size_t);
template void GpuVector<double>::clear();
template void GpuVector<double>::load(size_t, const double*);
template void GpuVector<double>::unload(size_t, double*) const;

template void GpuVector<std::complex<double>>::resize(size_t);
template void GpuVector<std::complex<double>>::clear();
template void GpuVector<std::complex<double>>::load(size_t, const std::complex<double>*);
template void GpuVector<std::complex<double>>::unload(size_t, std::complex<double>*) const;

template void GpuVector<float>::resize(size_t);
template void GpuVector<float>::clear();
template void GpuVector<float>::load(size_t, const float*);
template void GpuVector<float>::unload(size_t, float*) const;

template void GpuVector<int>::resize(size_t);
template void GpuVector<int>::clear();
template void GpuVector<int>::load(size_t, const int*);
template void GpuVector<int>::unload(size_t, int*) const;

GpuEngine::~GpuEngine(){
    if (own_rocblas_handle && rocblasHandle != nullptr){
        rocblas_destroy_handle(reinterpret_cast<rocblas_handle>(rocblasHandle));
        rocblasHandle = nullptr;
        own_rocblas_handle = false;
    }
}

int AccelerationMeta::getNumGpuDevices(){
    int gpu_count = 0;
    hipGetDeviceCount(&gpu_count);
    return gpu_count;
}
void AccelerationMeta::setDefaultCudaDevice(int deviceID){
    hipSetDevice(deviceID);
}
unsigned long long AccelerationMeta::getTotalGPUMemory(int deviceID){ // int deviceID
    hipDeviceProp_t prop;
    hipGetDeviceProperties(&prop, deviceID);
    return prop.totalGlobalMem;
}
std::string AccelerationMeta::getCudaDeviceName(int deviceID){ // int deviceID
    if ((deviceID < 0) || (deviceID >= getNumGpuDevices())) return std::string();

    hipDeviceProp_t prop;
    hipGetDeviceProperties(&prop, deviceID);

    return std::string(prop.name);
}
template<typename T> void AccelerationMeta::recvCudaArray(size_t num_entries, const T *gpu_data, std::vector<T> &cpu_data){
    cpu_data.resize(num_entries);
    TasGpu::hipcheck( hipMemcpy(cpu_data.data(), gpu_data, num_entries * sizeof(T), hipMemcpyDeviceToHost), "hip receive");
}
template<typename T> void AccelerationMeta::delCudaArray(T *x){
    TasGpu::hipcheck( hipFree(x), "hipFree()");
}

template void AccelerationMeta::recvCudaArray<double>(size_t num_entries, const double*, std::vector<double>&);
template void AccelerationMeta::recvCudaArray<float>(size_t num_entries, const float*, std::vector<float>&);
template void AccelerationMeta::recvCudaArray<int>(size_t num_entries, const int*, std::vector<int>&);

template void AccelerationMeta::delCudaArray<double>(double*);
template void AccelerationMeta::delCudaArray<float>(float*);
template void AccelerationMeta::delCudaArray<int>(int*);

namespace TasGpu{
/*
 * rocBLAS section
 */
//! \brief Converts character to cublas operation.
constexpr rocblas_operation cublas_trans(char trans){
    return (trans == 'N') ? rocblas_operation_none : ((trans == 'T') ? rocblas_operation_transpose : rocblas_operation_conjugate_transpose);
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

/*
 * rocSparse section
 */

/*
 * rocSolver section
 */

//! \brief Wrapper around cusolverDnDgetrf().
void factorizePLU(AccelerationContext const*, int, double[], int[]){}

void solvePLU(AccelerationContext const*, char, int, double const[], int const[], double []){}
void solvePLU(AccelerationContext const*, char, int, double const[], int const[], int, double[]){}

/*
 * Algorithm section
 */
template<typename scalar_type>
void solveLSmultiGPU(AccelerationContext const*, int, int, scalar_type[], int, scalar_type[]){}

template void solveLSmultiGPU<double>(AccelerationContext const*, int, int, double[], int, double[]);
template void solveLSmultiGPU<std::complex<double>>(AccelerationContext const*, int, int, std::complex<double>[], int, std::complex<double>[]);

template<typename scalar_type>
void solveLSmultiOOC(AccelerationContext const*, int, int, scalar_type[], int, scalar_type[]){}

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
void sparseMultiply(AccelerationContext const*, int, int, int, typename GpuVector<scalar_type>::value_type,
                    GpuVector<scalar_type> const&, GpuVector<int> const&, GpuVector<int> const&,
                    GpuVector<scalar_type> const&, scalar_type[]){
}

template void sparseMultiply<float>(AccelerationContext const*, int, int, int, float, GpuVector<float> const &A,
                                    GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<float> const &vals, float C[]);
template void sparseMultiply<double>(AccelerationContext const*, int, int, int, double, GpuVector<double> const &A,
                                     GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<double> const &vals, double C[]);

template<typename T> void load_n(T const *cpu_data, size_t num_entries, T *gpu_data){
    TasGpu::hipcheck( hipMemcpy(gpu_data, cpu_data, num_entries * sizeof(T), hipMemcpyHostToDevice), "hipMemcpy() load_n to device");
}

template void load_n<int>(int const*, size_t, int*);
template void load_n<float>(float const*, size_t, float*);
template void load_n<double>(double const*, size_t, double*);
template void load_n<std::complex<double>>(std::complex<double> const*, size_t, std::complex<double>*);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Placeholders for the kernels
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T> void dtrans2can(bool, int, int, int, double const[], double const[], T const[], T[]){}
template<typename T> void devalpwpoly(int, TypeOneDRule, int, int, int, const T[], const T[], const T[], T[]){}

template<typename T>
void devalpwpoly_sparse(int, TypeOneDRule, int, int, int, const T[], GpuVector<T> const&, GpuVector<T> const&,
                                GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int>&, GpuVector<int>&, GpuVector<T>&){}

template<typename T>
void devalseq(int, int, std::vector<int> const&, T const[], GpuVector<int> const&, GpuVector<int> const&, GpuVector<T> const&, GpuVector<T> const&, T[]){}
template<typename T>
void devalfor(int, int, std::vector<int> const&, const T[], GpuVector<int> const&, GpuVector<int> const&, T[], typename GpuVector<T>::value_type[]){}

template<typename T>
void devalglo(bool, bool, int, int, int, int, T const[], GpuVector<T> const&, GpuVector<T> const&, GpuVector<T> const&,
              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, T[]){}

void fillDataGPU(double, long long, long long, double[]){}

template void dtrans2can<float>(bool, int, int, int, double const*, double const*, float const*, float*);
template void dtrans2can<double>(bool, int, int, int, double const*, double const*, double const*, double*);

template void devalpwpoly<float>(int, TypeOneDRule, int, int, int, float const*, float const*, float const*, float*);
template void devalpwpoly<double>(int, TypeOneDRule, int, int, int, double const*, double const*, double const*, double*);

template void devalpwpoly_sparse<float>(int, TypeOneDRule, int, int, int, float const*, GpuVector<float> const&, GpuVector<float> const&,
                                        GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                        GpuVector<int>&, GpuVector<int>&, GpuVector<float>&);
template void devalpwpoly_sparse<double>(int, TypeOneDRule, int, int, int, double const*, GpuVector<double> const&, GpuVector<double> const&,
                                         GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                         GpuVector<int>&, GpuVector<int>&, GpuVector<double>&);

template void devalfor<float>(int, int, std::vector<int> const&, float const*, GpuVector<int> const&, GpuVector<int> const&, float*, float*);
template void devalfor<double>(int, int, std::vector<int> const&, double const*, GpuVector<int> const&, GpuVector<int> const&, double*, double*);

template void devalseq<float>(int, int, std::vector<int> const&, float const*, GpuVector<int> const&,
                              GpuVector<int> const&, GpuVector<float> const&, GpuVector<float> const&, float*);
template void devalseq<double>(int, int, std::vector<int> const&, double const*, GpuVector<int> const&,
                               GpuVector<int> const&, GpuVector<double> const&, GpuVector<double> const&, double*);

template void devalglo<float>(bool, bool, int, int, int, int,
                              float const*, GpuVector<float> const&, GpuVector<float> const&, GpuVector<float> const&,
                              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, float*);
template void devalglo<double>(bool, bool, int, int, int, int,
                               double const*, GpuVector<double> const&, GpuVector<double> const&, GpuVector<double> const&,
                               GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                               GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                               GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, double*);

}
}

#endif
