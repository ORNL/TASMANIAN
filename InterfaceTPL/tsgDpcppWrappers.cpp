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

void InternalSyclQueue::init_testing(){
    use_testing = true;
    test_queue = makeNewQueue(0);
}
InternalSyclQueue test_queue;

template<typename T> void GpuVector<T>::resize(AccelerationContext const *acc, size_t count){
    if (count != num_entries){ // if the current array is not big enoug
        clear(); // resets dynamic_mode
        num_entries = count;
        sycl::queue *q = getSyclQueue(acc);
        sycl_queue = reinterpret_cast<void*>(q);
        gpu_data = sycl::malloc_device<T>(num_entries, *q);
    }
}
template<typename T> void GpuVector<T>::clear(){
    num_entries = 0;
    if (gpu_data != nullptr){
        sycl::queue *q = reinterpret_cast<sycl::queue*>(sycl_queue);
        sycl::free(gpu_data, *q);
        q->wait(); // wait is needed here in case the pointer q gets deleted
    }
    gpu_data = nullptr;
}
template<typename T> void GpuVector<T>::load(AccelerationContext const *acc, size_t count, const T* cpu_data){
    resize(acc, count);
    sycl::queue *q = getSyclQueue(acc);
    q->memcpy(gpu_data, cpu_data, count * sizeof(T)).wait();
}
template<typename T> void GpuVector<T>::unload(AccelerationContext const *acc, size_t num, T* cpu_data) const{
    sycl::queue *q = getSyclQueue(acc);
    q->memcpy(cpu_data, gpu_data, num * sizeof(T)).wait();
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

template<> void deleteHandle<AccHandle::Syclqueue>(int *p){
    sycl::queue *q = reinterpret_cast<sycl::queue*>(p);
    delete q;
}

void GpuEngine::setSyclQueue(void *queue){
    internal_queue = std::unique_ptr<int, HandleDeleter<AccHandle::Syclqueue>>
        (reinterpret_cast<int*>(queue), HandleDeleter<AccHandle::Syclqueue>(false));
}

int AccelerationMeta::getNumGpuDevices(){
    return readSyclDevices().names.size();
}
void AccelerationMeta::setDefaultGpuDevice(int){}
unsigned long long AccelerationMeta::getTotalGPUMemory(int deviceID){ // int deviceID
    return readSyclDevices().memory[deviceID];
}
std::string AccelerationMeta::getGpuDeviceName(int deviceID){ // int deviceID
    return readSyclDevices().names[deviceID];
}
template<typename T> void AccelerationMeta::recvGpuArray(AccelerationContext const *acc, size_t num_entries, const T *gpu_data, std::vector<T> &cpu_data){
    sycl::queue *q = getSyclQueue(acc);
    cpu_data.resize(num_entries);
    q->memcpy(cpu_data.data(), gpu_data, num_entries * sizeof(T)).wait();
}
template<typename T> void AccelerationMeta::delGpuArray(AccelerationContext const *acc, T *x){
    sycl::queue *q = getSyclQueue(acc);
    sycl::free(x, *q);
    q->wait(); // wait is needed here in case the pointer q gets deleted
}

template void AccelerationMeta::recvGpuArray<double>(AccelerationContext const*, size_t num_entries, const double*, std::vector<double>&);
template void AccelerationMeta::recvGpuArray<float>(AccelerationContext const*, size_t num_entries, const float*, std::vector<float>&);
template void AccelerationMeta::recvGpuArray<int>(AccelerationContext const*, size_t num_entries, const int*, std::vector<int>&);

template void AccelerationMeta::delGpuArray<double>(AccelerationContext const*, double*);
template void AccelerationMeta::delGpuArray<float>(AccelerationContext const*, float*);
template void AccelerationMeta::delGpuArray<int>(AccelerationContext const*, int*);

namespace TasGpu{

template<typename scalar_type> struct tsg_transpose{};

template<typename scalar_type>
void transpose_matrix(sycl::queue *q, int m, int n, scalar_type const A[], scalar_type AT[]){
    q->submit([&](sycl::handler& h){
            h.parallel_for<tsg_transpose<scalar_type>>(sycl::range<2>{static_cast<size_t>(m), static_cast<size_t>(n)}, [=](sycl::id<2> i){
                AT[i[0] * n + i[1]] = A[i[0] + m * i[1]];
            });
        }).wait();
}

//! \brief Wrapper around rocsolver_dgetrf().
void factorizePLU(AccelerationContext const *acceleration, int n, double A[], int_gpu_lapack ipiv[]){
    sycl::queue *q = getSyclQueue(acceleration);
    size_t size = oneapi::mkl::lapack::getrf_scratchpad_size<double>(*q, n, n, n);
    GpuVector<double> workspace(acceleration, size);
    oneapi::mkl::lapack::getrf(*q, n, n, A, n, ipiv, workspace.data(), size).wait();
}

void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int_gpu_lapack const ipiv[], double b[]){
    sycl::queue *q = getSyclQueue(acceleration);
    size_t size = oneapi::mkl::lapack::getrs_scratchpad_size<double>(*q, (trans == 'T') ? oneapi::mkl::transpose::T :oneapi::mkl::transpose::N, n, 1, n, n);
    GpuVector<double> workspace(acceleration, size);
    oneapi::mkl::lapack::getrs(*q, (trans == 'T') ? oneapi::mkl::transpose::T :oneapi::mkl::transpose::N, n, 1,
                               const_cast<double*>(A), n, const_cast<int_gpu_lapack*>(ipiv), b, n, workspace.data(), size).wait();
}
void solvePLU(AccelerationContext const *acceleration, char trans, int n, double const A[], int_gpu_lapack const ipiv[], int nrhs, double B[]){
    sycl::queue *q = getSyclQueue(acceleration);
    size_t size = oneapi::mkl::lapack::getrs_scratchpad_size<double>(*q, (trans == 'T') ? oneapi::mkl::transpose::T :oneapi::mkl::transpose::N, n, nrhs, n, n);
    GpuVector<double> workspace(acceleration, size);
    GpuVector<double> BT(acceleration, n, nrhs);
    transpose_matrix(q, nrhs, n, B, BT.data());
    oneapi::mkl::lapack::getrs(*q, (trans == 'T') ? oneapi::mkl::transpose::T :oneapi::mkl::transpose::N, n, nrhs,
                               const_cast<double*>(A), n, const_cast<int_gpu_lapack*>(ipiv), BT.data(), n,
                               workspace.data(), size).wait();
    transpose_matrix(q, n, nrhs, BT.data(), B);
}

template<typename scalar_type>
void dump_data(sycl::queue *q, size_t m, size_t n, scalar_type *data){
    std::vector<scalar_type> cpu_data(m * n);
    q->memcpy(cpu_data.data(), data, Utils::size_mult(m, n) * sizeof(scalar_type)).wait();
    std::cout << std::scientific; std::cout.precision(4);
    for(size_t i=0; i<m; i++){
        for(size_t j=0; j<n; j++){
            std::cout << std::setw(15) << cpu_data[j * m + i];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template<typename scalar_type>
std::int64_t gemqr_workspace(cl::sycl::queue *q, oneapi::mkl::side side, oneapi::mkl::transpose trans, std::int64_t m, std::int64_t n, std::int64_t k, std::int64_t lda, std::int64_t ldc){
    return oneapi::mkl::lapack::ormqr_scratchpad_size<scalar_type>(*q, side, trans, m, n, k, lda, ldc);
}

template<>
std::int64_t gemqr_workspace<std::complex<double>>(cl::sycl::queue *q, oneapi::mkl::side side, oneapi::mkl::transpose trans, std::int64_t m, std::int64_t n, std::int64_t k, std::int64_t lda, std::int64_t ldc){
    return oneapi::mkl::lapack::unmqr_scratchpad_size<std::complex<double>>(*q, side, trans, m, n, k, lda, ldc);
}

template<typename scalar_type>
inline void gemqr(cl::sycl::queue *q, oneapi::mkl::side side, oneapi::mkl::transpose trans, std::int64_t m, std::int64_t n, std::int64_t k, scalar_type* A, std::int64_t lda, scalar_type *T, scalar_type* C, std::int64_t ldc, scalar_type* workspace, std::int64_t worksize){
    oneapi::mkl::lapack::ormqr(*q, side, trans, m, n, k, A, lda, T, C, ldc, workspace, worksize).wait();
}

template<>
void gemqr<std::complex<double>>(cl::sycl::queue *q, oneapi::mkl::side side, oneapi::mkl::transpose trans, std::int64_t m, std::int64_t n, std::int64_t k, std::complex<double>* A, std::int64_t lda, std::complex<double> *T, std::complex<double>* C, std::int64_t ldc, std::complex<double>* workspace, std::int64_t worksize){
    oneapi::mkl::lapack::unmqr(*q, side, trans, m, n, k, A, lda, T, C, ldc, workspace, worksize).wait();
}

/*
 * Algorithm section
 */
template<typename scalar_type>
void solveLSmultiGPU(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]){
    sycl::queue *q = getSyclQueue(acceleration);

    auto side = oneapi::mkl::side::left;
    auto trans = (std::is_same<scalar_type, double>::value) ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::conjtrans;

    std::int64_t worksize = oneapi::mkl::lapack::geqrf_scratchpad_size<scalar_type>(*q, n, m, n);
    worksize = std::max(worksize, gemqr_workspace<scalar_type>(q, side, trans, n, nrhs, m, n, n));

    GpuVector<scalar_type> AT(acceleration, n, m);
    GpuVector<scalar_type> T(acceleration, m); // tau parameter, or the weights of the orthogonal shift
    transpose_matrix(q, m, n, A, AT.data());

    GpuVector<scalar_type> workspace(acceleration, worksize);
    try{
        oneapi::mkl::lapack::geqrf(*q, n, m, AT.data(), n, T.data(), workspace.data(), worksize).wait();
    }catch(oneapi::mkl::lapack::exception &e){
        std::cout << "lapack geqrf() error code:  " << e.info()
                  << "\nlapack geqrf() detail code: " << e.detail() << std::endl;
        throw;
    }

    if (nrhs == 1){
        try{
            gemqr(q, side, trans, n, 1, m, AT.data(), n, T.data(), B, n, workspace.data(), worksize);
        }catch(oneapi::mkl::lapack::exception &e){
            std::cout << "lapack " << ((std::is_same<scalar_type, double>::value) ? "ormqr()" : "unmqr()")
                      << " error code:  " << e.info()
                      << "\nlapack " << ((std::is_same<scalar_type, double>::value) ? "ormqr()" : "unmqr()")
                      << " detail code:  " << e.detail() << std::endl;
            throw;
        }
        oneapi::mkl::blas::column_major::trsv(*q, oneapi::mkl::uplo::U, oneapi::mkl::transpose::N, oneapi::mkl::diag::N, m, AT.data(), n, B, 1).wait();
    }else{
        GpuVector<scalar_type> BT(acceleration, n, nrhs);
        transpose_matrix(q, nrhs, n, B, BT.data());
        gemqr(q, side, trans, n, nrhs, m, AT.data(), n, T.data(), BT.data(), n, workspace.data(), worksize);
        oneapi::mkl::blas::column_major::trsm(*q, oneapi::mkl::side::L, oneapi::mkl::uplo::U, oneapi::mkl::transpose::N, oneapi::mkl::diag::N, m,  nrhs, 1.0, AT.data(), n, BT.data(), n).wait();
        transpose_matrix(q, n, nrhs, BT.data(), B);
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
    sycl::queue *q = getSyclQueue(acceleration);
    if (M > 1){
        if (N > 1){ // matrix-matrix mode
            oneapi::mkl::blas::column_major::gemm(*q, oneapi::mkl::transpose::N, oneapi::mkl::transpose::N, M, N, K, alpha, A.data(), M, B.data(), K, beta, C, M).wait();
        }else{ // matrix vector, A * v = C
            oneapi::mkl::blas::column_major::gemv(*q, oneapi::mkl::transpose::N, M, K, alpha, A.data(), M, B.data(), 1, beta, C, 1).wait();
        }
    }else{ // matrix vector B^T * v = C
        oneapi::mkl::blas::column_major::gemv(*q, oneapi::mkl::transpose::T, K, N, alpha, B.data(), K, A.data(), 1, beta, C, 1).wait();
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
    sycl::queue *q = getSyclQueue(acceleration);
    oneapi::mkl::sparse::matrix_handle_t mat = nullptr;
    oneapi::mkl::sparse::init_matrix_handle(&mat);

    oneapi::mkl::sparse::set_csr_data(mat, N, K, oneapi::mkl::index_base::zero,
                                      const_cast<int*>(pntr.data()), const_cast<int*>(indx.data()), const_cast<scalar_type*>(vals.data()));

    if (M == 1){ // using sparse-blas level 2
        oneapi::mkl::sparse::gemv(*q, oneapi::mkl::transpose::nontrans, alpha, mat,
                                  const_cast<scalar_type*>(A.data()), 0.0, C).wait();
    }else{ // using sparse-blas level 3
        oneapi::mkl::sparse::gemm(*q, oneapi::mkl::layout::row_major, oneapi::mkl::transpose::nontrans, oneapi::mkl::transpose::nontrans,
                                  alpha, mat, const_cast<scalar_type*>(A.data()), M, M, 0.0, C, M).wait();
    }

    oneapi::mkl::sparse::release_matrix_handle(&mat);
}

template void sparseMultiply<float>(AccelerationContext const*, int, int, int, float, GpuVector<float> const &A,
                                    GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<float> const &vals, float C[]);
template void sparseMultiply<double>(AccelerationContext const*, int, int, int, double, GpuVector<double> const &A,
                                     GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<double> const &vals, double C[]);

template<typename T> void load_n(AccelerationContext const *acc, T const *cpu_data, size_t num_entries, T *gpu_data){
    sycl::queue *q = getSyclQueue(acc);
    q->memcpy(gpu_data, cpu_data, num_entries * sizeof(T)).wait();
}

template void load_n<int>(AccelerationContext const*, int const*, size_t, int*);
template void load_n<float>(AccelerationContext const*, float const*, size_t, float*);
template void load_n<double>(AccelerationContext const*, double const*, size_t, double*);
template void load_n<std::complex<double>>(AccelerationContext const*, std::complex<double> const*, size_t, std::complex<double>*);

}
}

#endif
