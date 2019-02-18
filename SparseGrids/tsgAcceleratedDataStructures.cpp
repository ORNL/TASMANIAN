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

#ifndef __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_CPP
#define __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_CPP

#include "tsgAcceleratedDataStructures.hpp"

#include "tsgCudaMacros.hpp"

#ifdef Tasmanian_ENABLE_MAGMA
#include "magma_v2.h"
#include "magmasparse.h"
#endif

namespace TasGrid{

#ifdef Tasmanian_ENABLE_CUDA
template<typename T> void CudaVector<T>::resize(size_t count){
    if (count != num_entries){ // if the current array is not big enough
        clear(); // resets dynamic_mode
        num_entries = count;
        cudaError_t cudaStat = cudaMalloc(((void**) &gpu_data), num_entries * sizeof(T));
        AccelerationMeta::cudaCheckError((void*) &cudaStat, "CudaVector::resize(), call to cudaMalloc()");
    }
}
template<typename T> void CudaVector<T>::clear(){
    num_entries = 0;
    if (dynamic_mode && (gpu_data != nullptr)){ // if I own the data and the data is not null
        cudaError_t cudaStat = cudaFree(gpu_data);
        AccelerationMeta::cudaCheckError((void*) &cudaStat, "CudaVector::clear(), call to cudaFree()");
    }
    gpu_data = nullptr;
    dynamic_mode = true; // The old data is gone and I own the current (null) data
}
template<typename T> void CudaVector<T>::load(size_t count, const T* cpu_data){
    resize(count);
    cudaError_t cudaStat = cudaMemcpy(gpu_data, cpu_data, num_entries * sizeof(T), cudaMemcpyHostToDevice);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "CudaVector::load(), call to cudaMemcpy()");
}
template<typename T> void CudaVector<T>::unload(T* cpu_data) const{
    cudaError_t cudaStat = cudaMemcpy(cpu_data, gpu_data, num_entries * sizeof(T), cudaMemcpyDeviceToHost);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "CudaVector::unload(), call to cudaMemcpy()");
}

template void CudaVector<double>::resize(size_t);
template void CudaVector<double>::clear();
template void CudaVector<double>::load(size_t, const double*);
template void CudaVector<double>::unload(double*) const;

template void CudaVector<int>::resize(size_t);
template void CudaVector<int>::clear();
template void CudaVector<int>::load(size_t, const int*);
template void CudaVector<int>::unload(int*) const;

CudaEngine::~CudaEngine(){
    if (cublasHandle != nullptr){
        cublasDestroy((cublasHandle_t) cublasHandle);
        cublasHandle = nullptr;
    }
    if (cusparseHandle != nullptr){
        cusparseDestroy((cusparseHandle_t) cusparseHandle);
        cusparseHandle = nullptr;
    }
    #ifdef Tasmanian_ENABLE_MAGMA
    if (magmaCudaQueue != nullptr){
        magma_queue_destroy((magma_queue*) magmaCudaQueue);
        magmaCudaQueue = nullptr;
        magma_finalize();
    }
    if (magmaCudaStream != nullptr) cudaStreamDestroy((cudaStream_t) magmaCudaStream);
    magmaCudaStream = nullptr;
    #endif
}
void CudaEngine::cuBlasPrepare(){
    if (cublasHandle == nullptr){
        cublasHandle_t cbh;
        cublasCreate(&cbh);
        cublasHandle = (void*) cbh;
    }
}
void CudaEngine::cuSparsePrepare(){
    if (cusparseHandle == nullptr){
        cusparseHandle_t csh;
        cusparseCreate(&csh);
        cusparseHandle = (void*) csh;
    }
}
void CudaEngine::magmaPrepare(){
    #ifdef Tasmanian_ENABLE_MAGMA
    if (magmaCudaQueue == nullptr){
        magma_init();
        cuBlasPrepare();
        cuSparsePrepare();
        magma_queue_create_from_cuda(gpu, (cudaStream_t) magmaCudaStream, (cublasHandle_t) cublasHandle, (cusparseHandle_t) cusparseHandle, ((magma_queue**) &magmaCudaQueue));
    }
    #endif
}
void CudaEngine::denseMultiply(int M, int N, int K, double alpha, const CudaVector<double> &A, const CudaVector<double> &B, double beta, CudaVector<double> &C){
    #ifdef Tasmanian_ENABLE_MAGMA
    if (magma){
        magmaPrepare();
        if (M > 1){
            if (N > 1){ // matrix mode
                magma_dgemm(MagmaNoTrans, MagmaNoTrans, M, N, K,
                            alpha, A.data(), M, B.data(), K, beta, C.data(), M, (magma_queue_t) magmaCudaQueue);
            }else{ // matrix vector, A * v = C
                magma_dgemv(MagmaNoTrans, M, K,
                            alpha, A.data(), M, B.data(), 1, beta, C.data(), 1, (magma_queue_t) magmaCudaQueue);
            }
        }else{ // matrix vector B^T * v = C
            magma_dgemv(MagmaTrans, N, K,
                        alpha, B.data(), K, A.data(), 1, beta, C.data(), 1, (magma_queue_t) magmaCudaQueue);
        }
        return;
    }
    #endif
    cublasStatus_t stat;
    cuBlasPrepare();
    if (M > 1){
        if (N > 1){ // matrix mode
            stat = cublasDgemm((cublasHandle_t) cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, K,
                                &alpha, A.data(), M, B.data(), K, &beta, C.data(), M);
        }else{ // matrix vector, A * v = C
            stat= cublasDgemv((cublasHandle_t) cublasHandle, CUBLAS_OP_N, M, K,
                            &alpha, A.data(), M, B.data(), 1, &beta, C.data(), 1);
        }
    }else{ // matrix vector B^T * v = C
        stat= cublasDgemv((cublasHandle_t) cublasHandle, CUBLAS_OP_T, N, K,
                            &alpha, B.data(), K, A.data(), 1, &beta, C.data(), 1);
    }
    AccelerationMeta::cublasCheckError((void*) &stat, "while calling CudaEngine::denseMultiply()");
}
void CudaEngine::sparseMultiply(int M, int N, int K, double alpha, const CudaVector<double> &A,
                                const CudaVector<int> &pntr, const CudaVector<int> &indx, const CudaVector<double> &vals, double beta, CudaVector<double> &C){
    #ifdef Tasmanian_ENABLE_MAGMA
    //if (magma){ // TODO: Enable more MAGMA sparse capabilities
    //    return;
    //}
    #endif
    cusparseStatus_t sparse_stat;
    cuSparsePrepare();
    if (N > 1){ // dense matrix has many columns
        if (M > 1){ // dense matrix has many rows, use matrix-matrix algorithm
            cusparseMatDescr_t mat_desc;
            sparse_stat = cusparseCreateMatDescr(&mat_desc);
            AccelerationMeta::cusparseCheckError((void*) &sparse_stat, "cusparseCreateMatDescr() in CudaEngine::sparseMultiply()");
            cusparseSetMatType(mat_desc, CUSPARSE_MATRIX_TYPE_GENERAL);
            cusparseSetMatIndexBase(mat_desc, CUSPARSE_INDEX_BASE_ZERO);
            cusparseSetMatDiagType(mat_desc, CUSPARSE_DIAG_TYPE_NON_UNIT);

            CudaVector<double> tempC(C.size());
            sparse_stat = cusparseDcsrmm2((cusparseHandle_t) cusparseHandle,
                                          CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_TRANSPOSE, N, M, K, (int) indx.size(),
                                          &alpha, mat_desc, vals.data(), pntr.data(), indx.data(), A.data(), M, &beta, tempC.data(), N);
            AccelerationMeta::cusparseCheckError((void*) &sparse_stat, "cusparseDcsrmm2() in CudaEngine::sparseMultiply()");

            cusparseDestroyMatDescr(mat_desc);

            cuBlasPrepare();
            double talpha = 1.0, tbeta = 0.0;
            cublasStatus_t dense_stat;
            dense_stat = cublasDgeam((cublasHandle_t) cublasHandle,
                                     CUBLAS_OP_T, CUBLAS_OP_T, M, N, &talpha, tempC.data(), N, &tbeta, tempC.data(), N, C.data(), M);
            AccelerationMeta::cublasCheckError((void*) &dense_stat, "cublasDgeam() in CudaEngine::sparseMultiply()");
        }else{ // dense matrix has only one row, use sparse matrix times dense vector
            cusparseMatDescr_t mat_desc;
            sparse_stat = cusparseCreateMatDescr(&mat_desc);
            AccelerationMeta::cusparseCheckError((void*) &sparse_stat, "cusparseCreateMatDescr() in CudaEngine::sparseMultiply()");
            cusparseSetMatType(mat_desc, CUSPARSE_MATRIX_TYPE_GENERAL);
            cusparseSetMatIndexBase(mat_desc, CUSPARSE_INDEX_BASE_ZERO);
            cusparseSetMatDiagType(mat_desc, CUSPARSE_DIAG_TYPE_NON_UNIT);

            sparse_stat = cusparseDcsrmv((cusparseHandle_t) cusparseHandle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE, N, K, (int) indx.size(),
                                        &alpha, mat_desc, vals.data(), pntr.data(), indx.data(), A.data(), &beta, C.data());
            AccelerationMeta::cusparseCheckError((void*) &sparse_stat, "cusparseDcsrmv() in CudaEngine::sparseMultiply()");

            cusparseDestroyMatDescr(mat_desc);
        }
    }else{ // sparse matrix has only one column, use dense matrix times sparse vector
        // quote from Nvidia CUDA cusparse manual at https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-gemvi
        // "This function requires no extra storage for the general matrices when operation CUSPARSE_OPERATION_NON_TRANSPOSE is selected."
        // Yet, buffer is required when num_nz exceeds 32 even with CUSPARSE_OPERATION_NON_TRANSPOSE
        int buffer_size;
        sparse_stat = cusparseDgemvi_bufferSize((cusparseHandle_t) cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                M, K, (int) indx.size(), &buffer_size);
        AccelerationMeta::cusparseCheckError((void*) &sparse_stat, "cusparseDgemvi_bufferSize() in CudaEngine::sparseMultiply()");
        CudaVector<double> gpu_buffer((size_t) buffer_size);

        sparse_stat = cusparseDgemvi((cusparseHandle_t) cusparseHandle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE, M, K, &alpha, A.data(), M, (int) indx.size(), vals.data(),
                                        indx.data(), &beta, C.data(), CUSPARSE_INDEX_BASE_ZERO, gpu_buffer.data());
        AccelerationMeta::cusparseCheckError((void*) &sparse_stat, "cusparseDgemvi() in CudaEngine::sparseMultiply()");
    }
}
void CudaEngine::setDevice() const{ cudaSetDevice(gpu); }
#endif


TypeAcceleration AccelerationMeta::getIOAccelerationString(const char * name){
    if (strcmp(name, "cpu-blas") == 0){
        return accel_cpu_blas;
    }else if (strcmp(name, "gpu-default") == 0){
        return accel_gpu_default;
    }else if (strcmp(name, "gpu-cublas") == 0){
        return accel_gpu_cublas;
    }else if (strcmp(name, "gpu-cuda") == 0){
        return accel_gpu_cuda;
    }else if (strcmp(name, "gpu-magma") == 0){
        return accel_gpu_magma;
    }else{
        return accel_none;
    }
}
const char* AccelerationMeta::getIOAccelerationString(TypeAcceleration accel){
    switch (accel){
        case accel_cpu_blas:       return "cpu-blas";
        case accel_gpu_default:    return "gpu-default";
        case accel_gpu_cublas:     return "gpu-cublas";
        case accel_gpu_cuda:       return "gpu-cuda";
        case accel_gpu_magma:      return "gpu-magma";
        default: return "none";
    }
}
int AccelerationMeta::getIOAccelerationInt(TypeAcceleration accel){
    switch (accel){
        case accel_cpu_blas:       return 1;
        case accel_gpu_default:    return 3;
        case accel_gpu_cublas:     return 4;
        case accel_gpu_cuda:       return 5;
        case accel_gpu_magma:      return 6;
        default: return 0;
    }
}
TypeAcceleration AccelerationMeta::getIOIntAcceleration(int accel){
    switch (accel){
        case 1:  return accel_cpu_blas;
        case 3:  return accel_gpu_default;
        case 4:  return accel_gpu_cublas;
        case 5:  return accel_gpu_cuda;
        case 6:  return accel_gpu_magma;
        default: return accel_none;
    }
}
bool AccelerationMeta::isAccTypeGPU(TypeAcceleration accel){
    switch (accel){
        case accel_gpu_default:
        case accel_gpu_cublas:
        case accel_gpu_cuda:
        case accel_gpu_magma: return true;
        default:
            return false;
    }
}

TypeAcceleration AccelerationMeta::getAvailableFallback(TypeAcceleration accel){
    // sparse grids are evaluated in 2 stages:
    // - s1: convert multi-index to matrix B
    // - s2: multiply matrix B by stored matrix A
    // Mode   | Stage 1 device | Stage 2 device | Library for stage 2
    // CUBLAS |      CPU       |     GPU        | Nvidia cuBlas (or cuSparse)
    // CUDA   |      GPU       |     GPU        | Nvidia cuBlas (or cuSparse)
    // MAGMA  |      GPU*      |     GPU        | UTK magma and magma_sparse
    // BLAS   |      CPU       |     CPU        | BLAS
    // none   | all done on CPU, still using OpenMP (if available)
    // *if CUDA is not simultaneously available with MAGMA, then MAGMA will use the CPU for stage 1
    // Note: using CUDA without either cuBlas or MAGMA is a bad idea (it will still work, just slow)

    // accel_gpu_default should always point to the potentially "best" option (currently MAGMA)
    if (accel == accel_gpu_default) accel = accel_gpu_magma;
    #if !defined(Tasmanian_ENABLE_CUDA) || !defined(Tasmanian_ENABLE_MAGMA) || !defined(Tasmanian_ENABLE_BLAS)
    // if any of the 3 acceleration modes is missing, then add a switch statement to guard against setting that mode
    switch(accel){
        #ifndef Tasmanian_ENABLE_CUDA
        // if CUDA is missing: just use the CPU
        case accel_gpu_cublas:
        case accel_gpu_cuda:
            #ifdef Tasmanian_ENABLE_BLAS
            accel = accel_cpu_blas;
            #else
            accel = accel_none;
            #endif
            break;
        #endif // Tasmanian_ENABLE_CUDA
        #ifndef Tasmanian_ENABLE_MAGMA
        // MAGMA tries to use CUDA kernels with magma linear algebra, this CUDA is the next best thing
        case accel_gpu_magma:
            #ifdef Tasmanian_ENABLE_CUDA
            accel = accel_gpu_cuda;
            #elif defined(Tasmanian_ENABLE_BLAS)
            accel = accel_cpu_blas;
            #else
            accel = accel_none;
            #endif
            break;
        #endif // Tasmanian_ENABLE_MAGMA
        #ifndef Tasmanian_ENABLE_BLAS
        // if BLAS is missing, do not attempt to use the GPU but go directly to "none" mode
        case accel_cpu_blas:
            accel = accel_none;
            break;
        #endif // Tasmanian_ENABLE_BLAS
        default: // compiler complains if there is no explicit "default", even if empty
            break;
    }
    #endif
    return accel;
}

#ifdef Tasmanian_ENABLE_CUDA
void AccelerationMeta::cudaCheckError(void *cudaStatus, const char *info){
    if (*((cudaError_t*) cudaStatus) != cudaSuccess){
        std::string message = "ERROR: cuda failed at ";
        message += info;
        message += " with error: ";
        message += cudaGetErrorString(*((cudaError_t*) cudaStatus));
        throw std::runtime_error(message);
    }
}
void AccelerationMeta::cublasCheckError(void *cublasStatus, const char *info){
    if (*((cublasStatus_t*) cublasStatus) != CUBLAS_STATUS_SUCCESS){
        std::string message = "ERROR: cuBlas failed with code: ";
        if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_NOT_INITIALIZED){
            message += "CUBLAS_STATUS_NOT_INITIALIZED";
        }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_ALLOC_FAILED){
            message += "CUBLAS_STATUS_ALLOC_FAILED";
        }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_INVALID_VALUE){
            message += "CUBLAS_STATUS_INVALID_VALUE";
        }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_ARCH_MISMATCH){
            message += "CUBLAS_STATUS_ARCH_MISMATCH";
        }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_MAPPING_ERROR){
            message += "CUBLAS_STATUS_MAPPING_ERROR";
        }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_EXECUTION_FAILED){
            message += "CUBLAS_STATUS_EXECUTION_FAILED";
        }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_INTERNAL_ERROR){
            message += "CUBLAS_STATUS_INTERNAL_ERROR";
        }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_NOT_SUPPORTED){
            message += "CUBLAS_STATUS_NOT_SUPPORTED";
        }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_LICENSE_ERROR){
            message += "CUBLAS_STATUS_LICENSE_ERROR";
        }else{
            message += "UNKNOWN";
        }
        message += " at ";
        message += info;
        throw std::runtime_error(message);
    }
}
void AccelerationMeta::cusparseCheckError(void *cusparseStatus, const char *info){
    if (*((cusparseStatus_t*) cusparseStatus) != CUSPARSE_STATUS_SUCCESS){
        std::string message = "ERROR: cuSparse failed with code: ";
        if (*((cusparseStatus_t*) cusparseStatus) == CUSPARSE_STATUS_NOT_INITIALIZED){
            message += "CUSPARSE_STATUS_NOT_INITIALIZED";
        }else if (*((cusparseStatus_t*) cusparseStatus) == CUSPARSE_STATUS_ALLOC_FAILED){
            message += "CUSPARSE_STATUS_ALLOC_FAILED";
        }else if (*((cusparseStatus_t*) cusparseStatus) == CUSPARSE_STATUS_INVALID_VALUE){
            message += "CUSPARSE_STATUS_INVALID_VALUE";
        }else if (*((cusparseStatus_t*) cusparseStatus) == CUSPARSE_STATUS_ARCH_MISMATCH){
            message += "CUSPARSE_STATUS_ARCH_MISMATCH";
        }else if (*((cusparseStatus_t*) cusparseStatus) == CUSPARSE_STATUS_INTERNAL_ERROR){
            message += "CUSPARSE_STATUS_INTERNAL_ERROR";
        }else if (*((cusparseStatus_t*) cusparseStatus) == CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED){
            message += "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
        }else if (*((cusparseStatus_t*) cusparseStatus) == CUSPARSE_STATUS_EXECUTION_FAILED){
            message += "CUSPARSE_STATUS_EXECUTION_FAILED";
        }else{
            message += "UNKNOWN";
        }
        message += " at ";
        message += info;
        throw std::runtime_error(message);
    }
}
int AccelerationMeta::getNumCudaDevices(){
    int gpu_count = 0;
    cudaGetDeviceCount(&gpu_count);
    return gpu_count;
}
void AccelerationMeta::setDefaultCudaDevice(int deviceID){
    cudaSetDevice(deviceID);
}
unsigned long long AccelerationMeta::getTotalGPUMemory(int deviceID){
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, deviceID);
    return prop.totalGlobalMem;
}
char* AccelerationMeta::getCudaDeviceName(int deviceID){
    char *name = new char[1];
    name[0] = '\0';
    if ((deviceID < 0) || (deviceID >= getNumCudaDevices())) return name;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, deviceID);

    int c = 0; while(prop.name[c] != '\0'){ c++; }
    delete[] name;
    name = new char[c+1];
    for(int i=0; i<c; i++){ name[i] = prop.name[i]; }
    name[c] = '\0';

    return name;
}
template<typename T> void AccelerationMeta::recvCudaArray(size_t num_entries, const T *gpu_data, std::vector<T> &cpu_data){
    cpu_data.resize(num_entries);
    cudaError_t cudaStat = cudaMemcpy(cpu_data.data(), gpu_data, num_entries * sizeof(T), cudaMemcpyDeviceToHost);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "cudaRecv(type, type)");
}
template<typename T> void AccelerationMeta::delCudaArray(T *x){
    TasCUDA::cudaDel<T>(x);
}

template void AccelerationMeta::recvCudaArray<double>(size_t num_entries, const double*, std::vector<double>&);
template void AccelerationMeta::recvCudaArray<int>(size_t num_entries, const int*, std::vector<int>&);

template void AccelerationMeta::delCudaArray<double>(double*);
template void AccelerationMeta::delCudaArray<int>(int*);
#endif // Tasmanian_ENABLE_CUDA


#ifdef Tasmanian_ENABLE_CUDA
AccelerationDomainTransform::AccelerationDomainTransform() : num_dimensions(0), padded_size(0){}
AccelerationDomainTransform::~AccelerationDomainTransform(){}

void AccelerationDomainTransform::clear(){
    num_dimensions = 0;
    padded_size = 0;
    gpu_trans_a.clear();
    gpu_trans_b.clear();
}
bool AccelerationDomainTransform::empty(){ return (num_dimensions == 0); }
void AccelerationDomainTransform::load(const std::vector<double> &transform_a, const std::vector<double> &transform_b){
    // The points are stored contiguously in a vector with stride equal to num_dimensions
    // Using the contiguous memory in a contiguous fashion on the GPU implies that thread 0 works on dimension 0, thread 1 on dim 1 ...
    // But the number of dimensions is often way less than the number of threads
    // Therefore, we lump vectors together into large vectors of sufficient dimension
    // The dimension is least 512, but less than max CUDA threads 1024
    // The domain transforms are padded accordingly
    num_dimensions = (int) transform_a.size();
    padded_size = num_dimensions;
    while(padded_size < 512) padded_size += num_dimensions;

    std::vector<double> rate(padded_size);
    std::vector<double> shift(padded_size);
    int c = 0;
    for(int i=0; i<padded_size; i++){
        // instead of storing upper/lower limits (as in TasmanianSparseGrid) use rate and shift
        double diff = transform_b[c] - transform_a[c];
        rate[i] = 2.0 / diff;
        shift[i] = (transform_b[c] + transform_a[c]) / diff;
        c++;
        c = (c % num_dimensions);
    }

    gpu_trans_a.load(rate);
    gpu_trans_b.load(shift);
}
void AccelerationDomainTransform::getCanonicalPoints(bool use01, const double *gpu_transformed_x, int num_x, CudaVector<double> &gpu_canonical_x){
    gpu_canonical_x.resize(((size_t) num_dimensions) * ((size_t) num_x));
    TasCUDA::dtrans2can(use01, num_dimensions, num_x, padded_size, gpu_trans_a.data(), gpu_trans_b.data(), gpu_transformed_x, gpu_canonical_x.data());
}
#endif // Tasmanian_ENABLE_CUDA


}

#endif
