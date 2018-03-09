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

#if defined(TASMANIAN_CUBLAS) || defined(TASMANIAN_CUDA)
#include <cuda_runtime_api.h>
#include <cuda.h>
#endif

#ifdef TASMANIAN_CUBLAS
#include <cublas_v2.h>
#include <cusparse.h>
#endif // TASMANIAN_CUBLAS

#ifdef _TASMANIAN_DEBUG_
#define _IF_DEBUG_MACRO(x) x
#else
#define _IF_DEBUG_MACRO(x)
#endif

namespace TasGrid{

BaseAccelerationData::BaseAccelerationData(){}
BaseAccelerationData::~BaseAccelerationData(){}

AccelerationDataGPUFull::AccelerationDataGPUFull() : gpu_values(0), cublasHandle(0), cusparseHandle(0), logstream(0){}
AccelerationDataGPUFull::~AccelerationDataGPUFull(){
    #if defined(TASMANIAN_CUBLAS) || defined(TASMANIAN_CUDA)
    if (gpu_values != 0){
        cudaFree(gpu_values);
        gpu_values = 0;
    }
    #endif // TASMANIAN_CUBLAS || TASMANIAN_CUDA
    #ifdef TASMANIAN_CUBLAS
    if (cublasHandle != 0){
        cublasDestroy((cublasHandle_t) cublasHandle);
        cublasHandle = 0;
    }
    if (cusparseHandle != 0){
        cusparseDestroy((cusparseHandle_t) cusparseHandle);
        cusparseHandle = 0;
    }
    #endif // TASMANIAN_CUBLAS
}
void AccelerationDataGPUFull::makeCuBlasHandle(){
    #ifdef TASMANIAN_CUBLAS
    if (cublasHandle == 0){
        cublasHandle_t cbh;
        cublasCreate(&cbh);
        cublasHandle = (void*) cbh;
    }
    #endif // TASMANIAN_CUBLAS
}
void AccelerationDataGPUFull::makeCuSparseHandle(){
    #ifdef TASMANIAN_CUBLAS
    if (cusparseHandle == 0){
        cusparseHandle_t csh;
        cusparseCreate(&csh);
        cusparseHandle = (void*) csh;
    }
    #endif // TASMANIAN_CUBLAS
}

void AccelerationDataGPUFull::setLogStream(std::ostream *os){ logstream = os; }

bool AccelerationDataGPUFull::isCompatible(TypeAcceleration acc) const{ return AccelerationMeta::isAccTypeFullMemoryGPU(acc); }

#if defined(TASMANIAN_CUBLAS) || defined(TASMANIAN_CUDA)
void AccelerationDataGPUFull::loadGPUValues(int total_entries, const double *cpu_values){
    cudaError_t cudaStat = cudaMalloc(((void**) &gpu_values), total_entries * sizeof(double));
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_values", logstream);
    cudaStat = cudaMemcpy(gpu_values, cpu_values, total_entries * sizeof(double), cudaMemcpyHostToDevice);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy gpu_values", logstream);
}
#else
void AccelerationDataGPUFull::loadGPUValues(int, const double *){}
#endif // TASMANIAN_CUBLAS or TASMANIAN_CUDA
double* AccelerationDataGPUFull::getGPUValues() const{ return gpu_values; }

#if defined(TASMANIAN_CUBLAS) || defined(TASMANIAN_CUDA)
void AccelerationDataGPUFull::resetValuesAndSurpluses(){
    if (gpu_values != 0){ cudaFree(gpu_values); gpu_values = 0; }
}
#else
void AccelerationDataGPUFull::resetValuesAndSurpluses(){}
#endif

#ifdef TASMANIAN_CUBLAS
void AccelerationDataGPUFull::cublasDGEMV(int num_outputs, int num_points, const double cpu_weights[], double *cpu_result){
    makeCuBlasHandle(); // creates a handle only if one doesn't exist
    double *gpu_weights = 0;
    cudaError_t cudaStat = cudaMalloc(((void**) &gpu_weights), num_points * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_weights in dgemv", logstream);)
    cudaStat = cudaMemcpy(gpu_weights, cpu_weights, num_points * sizeof(double), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy gpu_weights in dgemv", logstream);)

    double *gpu_result = 0;
    cudaStat = cudaMalloc(((void**) &gpu_result), num_outputs * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_result in dgemv", logstream);)

    double alpha = 1.0, beta = 0.0;
    cublasStatus_t stat;
    stat = cublasDgemv((cublasHandle_t) cublasHandle, CUBLAS_OP_N, num_outputs, num_points, &alpha, gpu_values, num_outputs, gpu_weights, 1, &beta, gpu_result, 1);
    AccelerationMeta::cublasCheckError((void*) &stat, "cublasDgemv in dgemv", logstream);

    cudaStat = cudaMemcpy(cpu_result, gpu_result, num_outputs * sizeof(double), cudaMemcpyDeviceToHost);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "retrieve gpu_result in dgemv", logstream);

    cudaStat = cudaFree(gpu_result);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_result in dgemv", logstream);)
    cudaStat = cudaFree(gpu_weights);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_weights in dgemv", logstream);)
}
#else
void AccelerationDataGPUFull::cublasDGEMV(int, int, const double *, double *){}
#endif // TASMANIAN_CUBLAS

#ifdef TASMANIAN_CUBLAS
void AccelerationDataGPUFull::cublasDGEMM(int num_outputs, int num_points, int num_x, const double cpu_weights[], double *cpu_result){
    makeCuBlasHandle(); // creates a handle only if one doesn't exist
    double *gpu_weights = 0;
    cudaError_t cudaStat = cudaMalloc(((void**) &gpu_weights), num_points * num_x * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_weights in dgemm", logstream);)
    cudaStat = cudaMemcpy(gpu_weights, cpu_weights, num_points * num_x * sizeof(double), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy gpu_weights in dgemm", logstream);)

    double *gpu_result = 0;
    cudaStat = cudaMalloc(((void**) &gpu_result), num_outputs * num_x * sizeof(double));
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_result in dgemm", logstream);

    double alpha = 1.0, beta = 0.0;
    cublasStatus_t stat;
    stat = cublasDgemm((cublasHandle_t) cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, num_outputs, num_x, num_points, &alpha, gpu_values, num_outputs, gpu_weights, num_points, &beta, gpu_result, num_outputs);
    AccelerationMeta::cublasCheckError((void*) &stat, "cublasDgemm in dgemm", logstream);

    cudaStat = cudaMemcpy(cpu_result, gpu_result, num_outputs * num_x * sizeof(double), cudaMemcpyDeviceToHost);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "retrieve gpu_result in dgemm", logstream);

    cudaStat = cudaFree(gpu_result);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_result in dgemm", logstream);)
    cudaStat = cudaFree(gpu_weights);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_weights in dgemm", logstream);)
}
void AccelerationDataGPUFull::cusparseDCRMM2(int num_points, int num_outputs, int num_x, const int *cpu_pntr, const int *cpu_indx, const double *cpu_vals, double *cpu_result){
    makeCuSparseHandle(); // creates a handle only if one doesn't exist
    makeCuBlasHandle(); // creates a handle only if one doesn't exist
    cudaError_t cudaStat;
    int *gpu_pntr, *gpu_indx;
    double *gpu_vals;
    double *gpu_result, *gpu_result_t;

    int num_nz = cpu_pntr[num_x];

    cudaStat = cudaMalloc(((void**) &gpu_pntr), (num_x+1) * sizeof(int));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_pntr in DCRMM2", logstream);)
    cudaStat = cudaMemcpy(gpu_pntr, cpu_pntr, (num_x+1) * sizeof(int), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy gpu_pntr in DCRMM2", logstream);)

    cudaStat = cudaMalloc(((void**) &gpu_indx), num_nz * sizeof(int));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_indx in DCRMM2", logstream);)
    cudaStat = cudaMemcpy(gpu_indx, cpu_indx, num_nz * sizeof(int), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy gpu_indx in DCRMM2", logstream);)

    cudaStat = cudaMalloc(((void**) &gpu_vals), num_nz * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_vals in DCRMM2", logstream);)
    cudaStat = cudaMemcpy(gpu_vals, cpu_vals, num_nz * sizeof(double), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy gpu_vals in DCRMM2", logstream);)

    cudaStat = cudaMalloc(((void**) &gpu_result), num_x * num_outputs * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_result in DCRMM2", logstream);)
    cudaStat = cudaMalloc(((void**) &gpu_result_t), num_x * num_outputs * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_result_t in DCRMM2", logstream);)

    // call cusparse
    cusparseStatus_t stat;
    double alpha = 1.0, beta = 0.0;
    cusparseMatDescr_t mat_desc;
    stat = cusparseCreateMatDescr(&mat_desc);
    _IF_DEBUG_MACRO(AccelerationMeta::cusparseCheckError((void*) &stat, "alloc mat_desc in DCRMM2", logstream);)
    cusparseSetMatType(mat_desc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(mat_desc, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatDiagType(mat_desc, CUSPARSE_DIAG_TYPE_NON_UNIT);
    //mat_desc->MatrixType = CUSPARSE_MATRIX_TYPE_GENERAL;

    stat = cusparseDcsrmm2((cusparseHandle_t) cusparseHandle,
            CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_TRANSPOSE, num_x, num_outputs, num_points, num_nz,
            &alpha, mat_desc, gpu_vals, gpu_pntr, gpu_indx, gpu_values, num_outputs, &beta, gpu_result_t, num_x);
    AccelerationMeta::cusparseCheckError((void*) &stat, "cusparseDcsrmm2 in DCRMM2", logstream);

    cusparseDestroyMatDescr(mat_desc);

    // transpose the result! (incomplete sparse standard)
    cublasStatus_t bstat;
    bstat = cublasDgeam((cublasHandle_t) cublasHandle,
                        CUBLAS_OP_T, CUBLAS_OP_T, num_outputs, num_x,
                        &alpha, gpu_result_t, num_x, &beta, gpu_result_t, num_x, gpu_result, num_outputs);
    AccelerationMeta::cublasCheckError((void*) &bstat, "cublasDgeam in DCRMM2", logstream);

    cudaStat = cudaMemcpy(cpu_result, gpu_result, num_x * num_outputs * sizeof(double), cudaMemcpyDeviceToHost);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "retrieve gpu_result in DCRMM2", logstream);

    cudaStat = cudaFree(gpu_result_t); _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_result_t in DCRMM2", logstream);)
    cudaStat = cudaFree(gpu_result);   _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_result in DCRMM2", logstream);)
    cudaStat = cudaFree(gpu_vals);     _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_vals in DCRMM2", logstream);)
    cudaStat = cudaFree(gpu_indx);     _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_indx in DCRMM2", logstream);)
    cudaStat = cudaFree(gpu_pntr);     _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_pntr in DCRMM2", logstream);)
}
void AccelerationDataGPUFull::cusparseDCRSMM(int num_points, int num_outputs, const int *cpu_pntr, const int *cpu_indx, const double *cpu_vals, const double *values, double *surpluses){
    makeCuSparseHandle(); // creates a handle only if one doesn't exist
    makeCuBlasHandle(); // creates a handle only if one doesn't exist
    cudaError_t cudaStat;
    int *gpu_pntr, *gpu_indx;
    double *gpu_vals;
    double *gpu_a, *gpu_b;
    double alpha = 1.0, beta = 0.0;

    int num_nz = cpu_pntr[num_points];

    cudaStat = cudaMalloc(((void**) &gpu_pntr), (num_points+1) * sizeof(int));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_pntr in DCRSMM", logstream);)
    cudaStat = cudaMemcpy(gpu_pntr, cpu_pntr, (num_points+1) * sizeof(int), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy gpu_pntr in DCRSMM", logstream);)

    cudaStat = cudaMalloc(((void**) &gpu_indx), num_nz * sizeof(int));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_indx in DCRSMM", logstream);)
    cudaStat = cudaMemcpy(gpu_indx, cpu_indx, num_nz * sizeof(int), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy gpu_indx in DCRSMM", logstream);)

    cudaStat = cudaMalloc(((void**) &gpu_vals), num_nz * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_vals in DCRSMM", logstream);)
    cudaStat = cudaMemcpy(gpu_vals, cpu_vals, num_nz * sizeof(double), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy gpu_vals in DCRSMM", logstream);)

    cudaStat = cudaMalloc(((void**) &gpu_a), num_points * num_outputs * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_result in DCRSMM", logstream);)
    cudaStat = cudaMalloc(((void**) &gpu_b), num_points * num_outputs * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "alloc gpu_result_t in DCRSMM", logstream);)

    // Load values
    cudaStat = cudaMemcpy(gpu_a, values, num_points * num_outputs * sizeof(double), cudaMemcpyHostToDevice); // gpu_a = values
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "send matrix to gpu in cusparseDCRSMM", logstream);

    // Transpose values
    cublasStatus_t bstat;
    bstat = cublasDgeam((cublasHandle_t) cublasHandle,
                        CUBLAS_OP_T, CUBLAS_OP_T, num_points, num_outputs,
                        &alpha, gpu_a, num_outputs, &beta, gpu_a, num_outputs, gpu_b, num_points); // gpu_b = gpu_a^T = values^T
    AccelerationMeta::cublasCheckError((void*) &bstat, "cublasDgeam first in DCRSMM", logstream);

    // Set matrix properties
    cusparseStatus_t stat;
    cusparseMatDescr_t mat_desc;
    stat = cusparseCreateMatDescr(&mat_desc);
    _IF_DEBUG_MACRO(AccelerationMeta::cusparseCheckError((void*) &stat, "alloc mat_desc in DCRSMM", logstream);)
    cusparseSetMatType(mat_desc, CUSPARSE_MATRIX_TYPE_TRIANGULAR);
    cusparseSetMatIndexBase(mat_desc, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatDiagType(mat_desc, CUSPARSE_DIAG_TYPE_UNIT);
    cusparseSetMatFillMode(mat_desc, CUSPARSE_FILL_MODE_LOWER);

    // Analyze matrix
    cusparseSolveAnalysisInfo_t info;
    cusparseCreateSolveAnalysisInfo(&info);
    stat = cusparseDcsrsm_analysis((cusparseHandle_t) cusparseHandle,
                                   CUSPARSE_OPERATION_NON_TRANSPOSE, num_points, num_nz,
                                   mat_desc, gpu_vals, gpu_pntr, gpu_indx, info);
    AccelerationMeta::cusparseCheckError((void*) &stat, "cusparseDcsrsm_analysis in DCRSMM", logstream);

    // Solve
    stat = cusparseDcsrsm_solve((cusparseHandle_t) cusparseHandle,
                                 CUSPARSE_OPERATION_NON_TRANSPOSE, num_points, num_outputs,
                                 &alpha, mat_desc, gpu_vals, gpu_pntr, gpu_indx, info,
                                 gpu_b, num_points, gpu_a, num_points);
    AccelerationMeta::cusparseCheckError((void*) &stat, "cusparseDcsrsm_solve in DCRSMM", logstream); // gpu_a = S * gpu_b = s^{-T} values^T

    // transpose back
    bstat = cublasDgeam((cublasHandle_t) cublasHandle,
                        CUBLAS_OP_T, CUBLAS_OP_T, num_outputs, num_points,
                        &alpha, gpu_a, num_points, &beta, gpu_a, num_points, gpu_b, num_outputs); // gpu_b = gpu_a^T = values * s^{-1}
    AccelerationMeta::cublasCheckError((void*) &bstat, "cublasDgeam second in DCRSMM", logstream);

    cudaStat = cudaMemcpy(surpluses, gpu_b, num_points * num_outputs * sizeof(double), cudaMemcpyDeviceToHost); // gpu_a = values
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "retrieve answer from gpu in cusparseDCRSMM", logstream);

    cusparseDestroyMatDescr(mat_desc);

    cudaStat = cudaFree(gpu_a);     _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_a in cusparseDCRSMM", logstream);)
    cudaStat = cudaFree(gpu_b);     _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_b in cusparseDCRSMM", logstream);)
    cudaStat = cudaFree(gpu_vals);  _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_vals in cusparseDCRSMM", logstream);)
    cudaStat = cudaFree(gpu_indx);  _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_indx in cusparseDCRSMM", logstream);)
    cudaStat = cudaFree(gpu_pntr);  _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_pntr in cusparseDCRSMM", logstream);)
}
#else
void AccelerationDataGPUFull::cublasDGEMM(int, int, int, const double *, double *){}
void AccelerationDataGPUFull::cusparseDCRMM2(int, int, int, const int *, const int *, const double *, double *){}
void AccelerationDataGPUFull::cusparseDCRSMM(int, int, const int*, const int*, const double*, const double*, double*){}
#endif // TASMANIAN_CUBLAS




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
    }else if (strcmp(name, "gpu-fullmem") == 0){
        return accel_gpu_fullmemory;
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
        case accel_gpu_fullmemory: return "gpu-fullmem";
        default: return "none";
    }
}
int AccelerationMeta::getIOAccelerationInt(TypeAcceleration accel){
    switch (accel){
        case accel_cpu_blas:       return 1;
        case accel_gpu_fullmemory: return 2;
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
        case accel_gpu_magma:
        case accel_gpu_fullmemory: return true;
        default:
            return false;
    }
}
bool AccelerationMeta::isAccTypeGPU(TypeAcceleration accel){
    switch (accel){
        case accel_gpu_default:
        case accel_gpu_cublas:
        case accel_gpu_cuda:
        case accel_gpu_magma:
        case accel_gpu_fullmemory: return true;
        default:
            return false;
    }
}

#if defined(TASMANIAN_CUBLAS) || defined(TASMANIAN_CUDA)
void AccelerationMeta::cudaCheckError(void *cudaStatus, const char *info, std::ostream *os){
    if (*((cudaError_t*) cudaStatus) != cudaSuccess){
        if (os != 0){
            (*os) << "ERROR: cuda failed at " << info << " with error: " << endl;
            (*os) << cudaGetErrorString(*((cudaError_t*) cudaStatus)) << endl;
        }
    }
}
#else
void AccelerationMeta::cudaCheckError(void *, const char *, std::ostream *){}
#endif // TASMANIAN_CUBLAS or TASMANIAN_CUDA

#ifdef TASMANIAN_CUBLAS
void AccelerationMeta::cublasCheckError(void *cublasStatus, const char *info, std::ostream *os){
    if (*((cublasStatus_t*) cublasStatus) != CUBLAS_STATUS_SUCCESS){
        if (os != 0){
            (*os) << "ERROR: cublas failed with code: ";
            if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_NOT_INITIALIZED){
                (*os) << "CUBLAS_STATUS_NOT_INITIALIZED";
            }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_ALLOC_FAILED){
                (*os) << "CUBLAS_STATUS_ALLOC_FAILED";
            }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_INVALID_VALUE){
                (*os) << "CUBLAS_STATUS_INVALID_VALUE";
            }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_ARCH_MISMATCH){
                (*os) << "CUBLAS_STATUS_ARCH_MISMATCH";
            }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_MAPPING_ERROR){
                (*os) << "CUBLAS_STATUS_MAPPING_ERROR";
            }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_EXECUTION_FAILED){
                (*os) << "CUBLAS_STATUS_EXECUTION_FAILED";
            }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_INTERNAL_ERROR){
                (*os) << "CUBLAS_STATUS_INTERNAL_ERROR";
            }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_NOT_SUPPORTED){
                (*os) << "CUBLAS_STATUS_NOT_SUPPORTED";
            }else if (*((cublasStatus_t*) cublasStatus) == CUBLAS_STATUS_LICENSE_ERROR){
                (*os) << "CUBLAS_STATUS_LICENSE_ERROR";
            }else{
                (*os) << "UNKNOWN";
            }
            (*os) << " at " << info << endl;
        }
    }
}
void AccelerationMeta::cusparseCheckError(void *cusparseStatus, const char *info, std::ostream *os){
    if (*((cusparseStatus_t*) cusparseStatus) != CUSPARSE_STATUS_SUCCESS){
        if (os != 0){
            (*os)  << "ERROR: cusparse failed at " << info << endl;
        }
    }
}
#else
void AccelerationMeta::cublasCheckError(void *, const char *, std::ostream *){}
void AccelerationMeta::cusparseCheckError(void *, const char *, std::ostream *){}
#endif // TASMANIAN_CUBLAS
}

#endif
