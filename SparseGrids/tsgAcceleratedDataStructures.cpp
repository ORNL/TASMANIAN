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

BaseAccelerationData::BaseAccelerationData(){}
BaseAccelerationData::~BaseAccelerationData(){}

AccelerationDataGPUFull::AccelerationDataGPUFull() :
    gpu_values(0), gpu_nodes(0), gpu_support(0),
    gpu_hpntr(0), gpu_hindx(0), gpu_roots(0), logstream(0){
#ifdef Tasmanian_ENABLE_CUDA
    cublasHandle = 0;
    cusparseHandle = 0;
#endif
#ifdef Tasmanian_ENABLE_MAGMA
    magma_initialized = false; // call init once per object (must simplify later)
    magmaCudaStream = 0;
    magmaCudaQueue = 0;
#endif
}
AccelerationDataGPUFull::~AccelerationDataGPUFull(){
    #ifdef Tasmanian_ENABLE_CUDA
    resetGPULoadedData();
    if (cublasHandle != 0){
        cublasDestroy((cublasHandle_t) cublasHandle);
        cublasHandle = 0;
    }
    if (cusparseHandle != 0){
        cusparseDestroy((cusparseHandle_t) cusparseHandle);
        cusparseHandle = 0;
    }
    #endif // Tasmanian_ENABLE_CUDA
    #ifdef Tasmanian_ENABLE_MAGMA
    if (magmaCudaQueue != 0) magma_queue_destroy((magma_queue*) magmaCudaQueue);
    magmaCudaQueue = 0;
    if (magma_initialized) magma_finalize();
    if (magmaCudaStream != 0) cudaStreamDestroy((cudaStream_t) magmaCudaStream);
    magmaCudaStream = 0;
    #endif
}
void AccelerationDataGPUFull::makeCuBlasHandle(){
    #ifdef Tasmanian_ENABLE_CUDA
    if (cublasHandle == 0){
        cublasHandle_t cbh;
        cublasCreate(&cbh);
        cublasHandle = (void*) cbh;
    }
    #endif // Tasmanian_ENABLE_CUDA
}
void AccelerationDataGPUFull::makeCuSparseHandle(){
    #ifdef Tasmanian_ENABLE_CUDA
    if (cusparseHandle == 0){
        cusparseHandle_t csh;
        cusparseCreate(&csh);
        cusparseHandle = (void*) csh;
    }
    #endif // Tasmanian_ENABLE_CUDA
}

#ifdef Tasmanian_ENABLE_MAGMA
void AccelerationDataGPUFull::initializeMagma(int gpuID){
    if (!magma_initialized){
        magma_init();
        magma_initialized = true;
    }
    makeCuBlasHandle();
    makeCuSparseHandle();
    magma_queue_create_from_cuda(gpuID, (cudaStream_t) magmaCudaStream, (cublasHandle_t) cublasHandle, (cusparseHandle_t) cusparseHandle, ((magma_queue**) &magmaCudaQueue));
}
#else
void AccelerationDataGPUFull::initializeMagma(int){}
#endif

void AccelerationDataGPUFull::setLogStream(std::ostream *os){ logstream = os; }

bool AccelerationDataGPUFull::isCompatible(TypeAcceleration acc) const{ return AccelerationMeta::isAccTypeFullMemoryGPU(acc); }

#ifdef Tasmanian_ENABLE_CUDA
void AccelerationDataGPUFull::loadGPUValues(size_t total_entries, const double *cpu_values){
    gpu_values = TasCUDA::cudaSend<double>(total_entries, cpu_values, logstream);
}
void AccelerationDataGPUFull::resetGPULoadedData(){
    if (gpu_values != 0){ TasCUDA::cudaDel<double>(gpu_values, logstream); gpu_values = 0; }
    if (gpu_nodes != 0){ TasCUDA::cudaDel<double>(gpu_nodes, logstream); gpu_nodes = 0; }
    if (gpu_support != 0){ TasCUDA::cudaDel<double>(gpu_support, logstream); gpu_support = 0; }
    if (gpu_hpntr != 0){ TasCUDA::cudaDel<int>(gpu_hpntr, logstream); gpu_hpntr = 0; }
    if (gpu_hindx != 0){ TasCUDA::cudaDel<int>(gpu_hindx, logstream); gpu_hindx = 0; }
    if (gpu_roots != 0){ TasCUDA::cudaDel<int>(gpu_roots, logstream); gpu_roots = 0; }
}
void AccelerationDataGPUFull::loadGPUNodesSupport(int total_entries, const double *cpu_nodes, const double *cpu_support){
    gpu_nodes   = TasCUDA::cudaSend<double>(total_entries, cpu_nodes,   logstream);
    gpu_support = TasCUDA::cudaSend<double>(total_entries, cpu_support, logstream);
}
void AccelerationDataGPUFull::loadGPUHierarchy(int num_points, const int *pntr, const int *indx, int num_roots, const int *roots){
    gpu_hpntr = TasCUDA::cudaSend<int>(num_points + 1, pntr, logstream);
    gpu_hindx = TasCUDA::cudaSend<int>(pntr[num_points], indx, logstream);
    gpu_roots = TasCUDA::cudaSend<int>(num_roots, roots, logstream);
}
#else
void AccelerationDataGPUFull::loadGPUValues(size_t, const double *){}
void AccelerationDataGPUFull::resetGPULoadedData(){}
void AccelerationDataGPUFull::loadGPUNodesSupport(int, const double *, const double *){}
void AccelerationDataGPUFull::loadGPUHierarchy(int, const int*, const int*, int, const int*){}
#endif // Tasmanian_ENABLE_CUDA

double* AccelerationDataGPUFull::getGPUValues() const{ return gpu_values; }
double* AccelerationDataGPUFull::getGPUNodes() const{ return gpu_nodes; }
double* AccelerationDataGPUFull::getGPUSupport() const{ return gpu_support; }
int* AccelerationDataGPUFull::getGPUpntr() const{ return gpu_hpntr; }
int* AccelerationDataGPUFull::getGPUindx() const{ return gpu_hindx; }
int* AccelerationDataGPUFull::getGPUroots() const{ return gpu_roots; }

#ifdef Tasmanian_ENABLE_CUDA
void AccelerationDataGPUFull::cublasDGEMM(bool cpu_pointers, int num_outputs, int num_x, int num_points, const double weights[], double *result){
    makeCuBlasHandle(); // creates a handle only if one doesn't exist

    const double *gpu_weights = 0;
    double *gpu_result = 0, *gpu_temp = 0;

    gpu_weights = (cpu_pointers) ? TasCUDA::cudaSendConst<double>(num_x * num_points, weights, gpu_temp, logstream) : weights;
    gpu_result  = (cpu_pointers) ? TasCUDA::cudaNew<double>(num_outputs * num_x, logstream) : result;

    double alpha = 1.0, beta = 0.0;
    if (num_x > 1){ // matrix-matrix mode
        cublasStatus_t stat = cublasDgemm((cublasHandle_t) cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, num_outputs, num_x, num_points,
                                        &alpha, gpu_values, num_outputs, gpu_weights, num_points, &beta, gpu_result, num_outputs);
        AccelerationMeta::cublasCheckError((void*) &stat, "cublasDgemm in DGEMM", logstream);
    }else{ // matrix-vector mode
        cublasStatus_t stat= cublasDgemv((cublasHandle_t) cublasHandle, CUBLAS_OP_N, num_outputs, num_points,
                                        &alpha, gpu_values, num_outputs, gpu_weights, 1, &beta, gpu_result, 1);
        AccelerationMeta::cublasCheckError((void*) &stat, "cublasDgemv in DGEMV", logstream);
    }

    if (cpu_pointers){
        TasCUDA::cudaRecv<double>(num_outputs * num_x, gpu_result, result, logstream);
        TasCUDA::cudaDel<double>(gpu_result, logstream);
        TasCUDA::cudaDel<double>(gpu_temp, logstream);
    }
}
void AccelerationDataGPUFull::cusparseMatmul(bool cpu_pointers, int num_points, int num_outputs, int num_x, const int *spntr, const int *sindx, const double *svals, int num_nz, double *result){
    makeCuSparseHandle(); // creates a handle only if one doesn't exist
    makeCuBlasHandle(); // creates a handle only if one doesn't exist
    int *tempp = 0, *tempi = 0;
    double *tempv = 0;

    if (cpu_pointers) num_nz = spntr[num_x];
    const int *gpu_pntr    = (cpu_pointers) ? TasCUDA::cudaSendConst<int>(num_x + 1, spntr, tempp, logstream) : spntr;
    const int *gpu_indx    = (cpu_pointers) ? TasCUDA::cudaSendConst<int>(num_nz, sindx, tempi, logstream)    : sindx;
    const double *gpu_vals = (cpu_pointers) ? TasCUDA::cudaSendConst<double>(num_nz, svals, tempv, logstream) : svals;

    double *gpu_result = (cpu_pointers) ? TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) num_outputs), logstream) : result;

    double *gpu_result_t = TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) num_outputs), logstream);

    // call cusparse
    cusparseStatus_t stat;
    double alpha = 1.0, beta = 0.0;
    cusparseMatDescr_t mat_desc;
    stat = cusparseCreateMatDescr(&mat_desc);
    AccelerationMeta::cusparseCheckError((void*) &stat, "alloc mat_desc in Matmul", logstream);
    cusparseSetMatType(mat_desc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(mat_desc, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatDiagType(mat_desc, CUSPARSE_DIAG_TYPE_NON_UNIT);

    stat = cusparseDcsrmm2((cusparseHandle_t) cusparseHandle,
            CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_TRANSPOSE, num_x, num_outputs, num_points, num_nz,
            &alpha, mat_desc, gpu_vals, gpu_pntr, gpu_indx, gpu_values, num_outputs, &beta, gpu_result_t, num_x);
    AccelerationMeta::cusparseCheckError((void*) &stat, "cusparseDcsrmm2 in Matmul", logstream);

    cusparseDestroyMatDescr(mat_desc);

    // transpose the result! (incomplete sparse standard)
    cublasStatus_t bstat;
    bstat = cublasDgeam((cublasHandle_t) cublasHandle,
                        CUBLAS_OP_T, CUBLAS_OP_T, num_outputs, num_x,
                        &alpha, gpu_result_t, num_x, &beta, gpu_result_t, num_x, gpu_result, num_outputs);
    AccelerationMeta::cublasCheckError((void*) &bstat, "cublasDgeam in Matmul", logstream);

    TasCUDA::cudaDel<double>(gpu_result_t, logstream);

    if (cpu_pointers){
        TasCUDA::cudaRecv<double>(((size_t) num_x) * ((size_t) num_outputs), gpu_result, result, logstream);

        TasCUDA::cudaDel<double>(gpu_result, logstream);
        TasCUDA::cudaDel<double>(tempv, logstream);
        TasCUDA::cudaDel<int>(tempi, logstream);
        TasCUDA::cudaDel<int>(tempp, logstream);
    }
}
void AccelerationDataGPUFull::cusparseMatvec(int num_points, int num_x, const int *spntr, const int *sindx, const double *svals, int num_nz, double *result){
    makeCuSparseHandle(); // creates a handle only if one doesn't exist
    cusparseStatus_t stat;
    double alpha = 1.0, beta = 0.0;
    cusparseMatDescr_t mat_desc;
    stat = cusparseCreateMatDescr(&mat_desc);
    AccelerationMeta::cusparseCheckError((void*) &stat, "alloc mat_desc in Matvec", logstream);
    cusparseSetMatType(mat_desc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(mat_desc, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatDiagType(mat_desc, CUSPARSE_DIAG_TYPE_NON_UNIT);

    stat = cusparseDcsrmv((cusparseHandle_t) cusparseHandle,
            CUSPARSE_OPERATION_NON_TRANSPOSE, num_x, num_points, num_nz,
            &alpha, mat_desc, svals, spntr, sindx, gpu_values, &beta, result);
    AccelerationMeta::cusparseCheckError((void*) &stat, "cusparseDcsrmv in Matvec", logstream);

    cusparseDestroyMatDescr(mat_desc);
}
void AccelerationDataGPUFull::cusparseMatveci(int num_outputs, int num_points, int num_nz, const int *sindx, const double *svals, double *result){
    makeCuSparseHandle(); // creates a handle only if one doesn't exist
    cusparseStatus_t stat;
    double alpha = 1.0, beta = 0.0;

    // quote from Nvidia CUDA cusparse manual at https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-gemvi
    // "This function requires no extra storage for the general matrices when operation CUSPARSE_OPERATION_NON_TRANSPOSE is selected."
    // Yet, buffer is required when num_nz exceeds 32 even with CUSPARSE_OPERATION_NON_TRANSPOSE
    int buffer_size;
    double *gpu_buffer = 0;
    stat = cusparseDgemvi_bufferSize((cusparseHandle_t) cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     num_outputs, num_points, num_nz, &buffer_size);
    AccelerationMeta::cusparseCheckError((void*) &stat, "cusparseDgemvi_bufferSize in Matveci", logstream);
    if (buffer_size > 0){
        gpu_buffer = TasCUDA::cudaNew<double>(buffer_size, logstream);
    }

    stat = cusparseDgemvi((cusparseHandle_t) cusparseHandle,
            CUSPARSE_OPERATION_NON_TRANSPOSE, num_outputs, num_points, &alpha,
            gpu_values, num_outputs, num_nz, svals, sindx, &beta, result, CUSPARSE_INDEX_BASE_ZERO, gpu_buffer);
    AccelerationMeta::cusparseCheckError((void*) &stat, "cusparseDgemvi in Matveci", logstream);
    if (gpu_buffer != 0) TasCUDA::cudaDel<double>(gpu_buffer, logstream);
}
#else
void AccelerationDataGPUFull::cublasDGEMM(bool, int, int, int, const double *, double *){}
void AccelerationDataGPUFull::cusparseMatmul(bool, int, int, int, const int*, const int*, const double*, int, double*){}
void AccelerationDataGPUFull::cusparseMatvec(int, int, const int*, const int*, const double*, int, double*){}
void AccelerationDataGPUFull::cusparseMatveci(int, int, int, const int*, const double*, double*){}
#endif // Tasmanian_ENABLE_CUDA

#ifdef Tasmanian_ENABLE_MAGMA
void AccelerationDataGPUFull::magmaCudaDGEMM(bool cpu_pointers, int gpuID, int num_outputs, int num_x, int num_points, const double weights[], double *result){
    initializeMagma(gpuID); // calls magma_init(), but only if not called before by this object

    const double *gpu_weights = 0;
    double *gpu_result = 0, *gpu_temp = 0;

    gpu_weights = (cpu_pointers) ? TasCUDA::cudaSendConst<double>(num_x * num_points, weights, gpu_temp, logstream) : weights;
    gpu_result  = (cpu_pointers) ? TasCUDA::cudaNew<double>(num_outputs * num_x, logstream) : result;

    double alpha = 1.0, beta = 0.0;
    magma_trans_t noTranspose = MagmaNoTrans;
    if (num_x > 1){ // matrix-matrix mode
        magma_dgemm(noTranspose, noTranspose, num_outputs, num_x, num_points, alpha, gpu_values, num_outputs,
                    gpu_weights, num_points, beta, gpu_result, num_outputs, (magma_queue_t) magmaCudaQueue);
    }else{ // matrix-vector mode
        magma_dgemv(noTranspose, num_outputs, num_points, alpha, gpu_values, num_outputs,
                    gpu_weights, 1, beta, gpu_result, 1, (magma_queue_t) magmaCudaQueue);
    }

    if (cpu_pointers){
        TasCUDA::cudaRecv<double>(num_outputs * num_x, gpu_result, result, logstream);
        TasCUDA::cudaDel<double>(gpu_result, logstream);
        TasCUDA::cudaDel<double>(gpu_temp, logstream);
    }
}
#else
void AccelerationDataGPUFull::magmaCudaDGEMM(bool, int, int, int, int, const double[], double*){}
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
bool AccelerationMeta::isAccTypeFullMemoryGPU(TypeAcceleration accel){
    switch (accel){
        case accel_gpu_default:
        case accel_gpu_cublas:
        case accel_gpu_cuda:
        case accel_gpu_magma: return true;
        default:
            return false;
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
void AccelerationMeta::cudaCheckError(void *cudaStatus, const char *info, std::ostream *os){
    if (*((cudaError_t*) cudaStatus) != cudaSuccess){
        std::string message = "ERROR: cuda failed at ";
        message += info;
        message += " with error: ";
        message += cudaGetErrorString(*((cudaError_t*) cudaStatus));
        throw std::runtime_error(message);
    }
}
void AccelerationMeta::cublasCheckError(void *cublasStatus, const char *info, std::ostream *os){
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
void AccelerationMeta::cusparseCheckError(void *cusparseStatus, const char *info, std::ostream *os){
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
#else
void AccelerationMeta::cudaCheckError(void *, const char *, std::ostream *){}
void AccelerationMeta::cublasCheckError(void *, const char *, std::ostream *){}
void AccelerationMeta::cusparseCheckError(void *, const char *, std::ostream *){}
#endif // Tasmanian_ENABLE_CUDA


#ifdef Tasmanian_ENABLE_CUDA
AccelerationDomainTransform::AccelerationDomainTransform(int num_dimensions, const double *transform_a, const double *transform_b, std::ostream *os) :
    gpu_trans_a(0), gpu_trans_b(0), padded_size(0), logstream(os)
{
    padded_size = num_dimensions;
    while(padded_size < 512) padded_size += num_dimensions;

    double *rate = new double[padded_size];
    double *shift = new double[padded_size];
    int c = 0;
    for(int i=0; i<padded_size; i++){
        double diff = transform_b[c] - transform_a[c];
        rate[i] = 2.0 / diff;
        shift[i] = (transform_b[c] + transform_a[c]) / diff;
        c++;
        c = (c % num_dimensions);
    }

    gpu_trans_a = TasCUDA::cudaSend<double>(padded_size, rate, logstream);
    gpu_trans_b = TasCUDA::cudaSend<double>(padded_size, shift, logstream);

    delete[] rate;
    delete[] shift;
}
AccelerationDomainTransform::~AccelerationDomainTransform(){
    TasCUDA::cudaDel<double>(gpu_trans_a, logstream);
    TasCUDA::cudaDel<double>(gpu_trans_b, logstream);
}
double* AccelerationDomainTransform::getCanonicalPoints(int num_dimensions, int num_x, const double *gpu_transformed_x){
    double *gpu_x_canonical = TasCUDA::cudaNew<double>(num_dimensions * num_x, logstream);
    TasCUDA::dtrans2can(num_dimensions, num_x, padded_size, gpu_trans_a, gpu_trans_b, gpu_transformed_x, gpu_x_canonical);
    return gpu_x_canonical;
}
#else
AccelerationDomainTransform::AccelerationDomainTransform(int, const double*, const double*, std::ostream *){}
AccelerationDomainTransform::~AccelerationDomainTransform(){}
double* AccelerationDomainTransform::getCanonicalPoints(int, int, const double*){ return 0; }
#endif // Tasmanian_ENABLE_CUDA


}

#endif
