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
    if (count != num_entries){ // if the current array is not big enoug
        clear(); // resets dynamic_mode
        num_entries = count;
        sycl::queue *q = getSyclQueue(acc);
        sycl_queue = acc->engine->internal_queue;
        gpu_data = sycl::malloc_device<T>(num_entries, *q);
    }
}
template<typename T> void GpuVector<T>::clear(){
    num_entries = 0;
    if (gpu_data != nullptr){
        sycl::queue *q = reinterpret_cast<sycl::queue*>(sycl_queue.get());
        sycl::free(gpu_data, *q);
    }
    gpu_data = nullptr;
}
template<typename T> void GpuVector<T>::load(AccelerationContext const *acc, size_t count, const T* cpu_data){
    resize(acc, count);
    sycl::queue *q = getSyclQueue(acc);
    q->memcpy(gpu_data, cpu_data, count * sizeof(T));
    q->wait();
}
template<typename T> void GpuVector<T>::unload(AccelerationContext const *acc, size_t num, T* cpu_data) const{
    sycl::queue *q = getSyclQueue(acc);
    q->memcpy(cpu_data, gpu_data, num * sizeof(T));
    q->wait();
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

template void GpuVector<std::int64_t>::resize(AccelerationContext const*, size_t);
template void GpuVector<std::int64_t>::clear();
template void GpuVector<std::int64_t>::load(AccelerationContext const*, size_t, const std::int64_t*);
template void GpuVector<std::int64_t>::unload(AccelerationContext const*, size_t, std::int64_t*) const;

GpuEngine::~GpuEngine(){}

void GpuEngine::setSyclQueue(void *queue){
    sycl_gpu_queue = queue;
    own_gpu_queue = false;
}

int AccelerationMeta::getNumGpuDevices(){
    return 1; // fake device for now, will actually use the CPU
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

void transpose_matrix(sycl::queue *q, int m, int n, double const A[], double AT[]){
    q->submit([&](sycl::handler& h) {
            h.parallel_for<class tsg_transpose>(sycl::range<2>{static_cast<size_t>(m), static_cast<size_t>(n)}, [=](sycl::id<2> i){
                AT[i[0] * n + i[1]] = A[i[0] + m * i[1]];
            });
        });
    q->wait();
}

//! \brief Wrapper around rocsolver_dgetrf().
void factorizePLU(AccelerationContext const *acceleration, int n, double A[], int_gpu_lapack ipiv[]){
    sycl::queue *q = getSyclQueue(acceleration);
    size_t size = oneapi::mkl::lapack::getrf_scratchpad_size<double>(*q, n, n, n);
    q->wait();
    GpuVector<double> workspace(acceleration, size);
    oneapi::mkl::lapack::getrf(*q, n, n, A, n, ipiv, workspace.data(), size);
    q->wait();
}

void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int_gpu_lapack const ipiv[], double b[]){
    sycl::queue *q = getSyclQueue(acceleration);
    size_t size = oneapi::mkl::lapack::getrs_scratchpad_size<double>(*q, (trans == 'T') ? oneapi::mkl::transpose::T :oneapi::mkl::transpose::N, n, 1, n, n);
    q->wait();
    GpuVector<double> workspace(acceleration, size);
    oneapi::mkl::lapack::getrs(*q, (trans == 'T') ? oneapi::mkl::transpose::T :oneapi::mkl::transpose::N, n, 1,
                               const_cast<double*>(A), n, const_cast<int_gpu_lapack*>(ipiv), b, n, workspace.data(), size);
    q->wait();
}
void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int_gpu_lapack const ipiv[], int nrhs, double B[]){
    sycl::queue *q = getSyclQueue(acceleration);
    size_t size = oneapi::mkl::lapack::getrs_scratchpad_size<double>(*q, (trans == 'T') ? oneapi::mkl::transpose::T :oneapi::mkl::transpose::N, n, nrhs, n, n);
    q->wait();
    GpuVector<double> workspace(acceleration, size);
    GpuVector<double> BT(acceleration, n, nrhs);
    transpose_matrix(q, nrhs, n, B, BT.data());
    oneapi::mkl::lapack::getrs(*q, (trans == 'T') ? oneapi::mkl::transpose::T :oneapi::mkl::transpose::N, n, nrhs,
                               const_cast<double*>(A), n, const_cast<int_gpu_lapack*>(ipiv), BT.data(), n, workspace.data(), size);
    q->wait();
    transpose_matrix(q, n, nrhs, BT.data(), B);
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
    sycl::queue *q = getSyclQueue(acceleration);
    if (M > 1){
        if (N > 1){ // matrix-matrix mode
            oneapi::mkl::blas::column_major::gemm(*q, oneapi::mkl::transpose::N, oneapi::mkl::transpose::N, M, N, K, alpha, A.data(), M, B.data(), K, beta, C, M);
        }else{ // matrix vector, A * v = C
            oneapi::mkl::blas::column_major::gemv(*q, oneapi::mkl::transpose::N, M, K, alpha, A.data(), M, B.data(), 1, beta, C, 1);
        }
    }else{ // matrix vector B^T * v = C
        oneapi::mkl::blas::column_major::gemv(*q, oneapi::mkl::transpose::T, K, N, alpha, B.data(), K, A.data(), 1, beta, C, 1);
    }
    q->wait();
}

template void denseMultiply<float>(AccelerationContext const*, int, int, int, float,
                                   GpuVector<float> const&, GpuVector<float> const&, float, float[]);
template void denseMultiply<double>(AccelerationContext const*, int, int, int, double,
                                    GpuVector<double> const&, GpuVector<double> const&, double, double[]);

template<typename scalar_type>
void sparseMultiply(AccelerationContext const *acceleration, int M, int N, int K, typename GpuVector<scalar_type>::value_type alpha,
                    GpuVector<scalar_type> const &A, GpuVector<int> const &pntr, GpuVector<int> const &indx,
                    GpuVector<scalar_type> const &vals, scalar_type C[]){
    sycl::queue *q = getSyclQueue(acceleration);
    oneapi::mkl::sparse::matrix_handle_t mat = nullptr;
    oneapi::mkl::sparse::init_matrix_handle(&mat);

    oneapi::mkl::sparse::set_csr_data(mat, N, K, oneapi::mkl::index_base::zero,
                                      const_cast<int*>(pntr.data()), const_cast<int*>(indx.data()), const_cast<scalar_type*>(vals.data()));

    Utils::Wrapper2D<scalar_type> ywrap(M, C);
    Utils::Wrapper2D<scalar_type const> surpluses(M, A.data());
    int const *spntr = pntr.data();
    int const *sindx = indx.data();
    scalar_type const *svals = vals.data();
    #pragma omp parallel for
    for(int i=0; i<N; i++){
        scalar_type *this_y = ywrap.getStrip(i);
        std::fill(this_y, this_y + M, 0.0);
        for(int j=spntr[i]; j<spntr[i+1]; j++){
            double v = svals[j];
            scalar_type const *s = surpluses.getStrip(sindx[j]);
            for(int k=0; k<M; k++) this_y[k] += v * s[k];
        }
    }

    q->wait();
    oneapi::mkl::sparse::release_matrix_handle(&mat);
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
