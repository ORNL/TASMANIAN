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

#ifndef __TASMANIAN_SPARSE_GRID_CUDA_KERNELS_CU
#define __TASMANIAN_SPARSE_GRID_CUDA_KERNELS_CU

//#include <iostream> // for debugging

#include "tsgAcceleratedDataStructures.hpp"
//#include "tasCudaKernels.hpp"

#define TASMANIAN_CUDA_NUM_THREADS 1024

namespace TasGrid{
    
__global__
void tasgpu_d3gecs_1(int N, int M, int K, const double *gpuA, const int* gpuBpntr, const int *gpuBindx, const double *gpuBvals, double *gpuC){
    __shared__ int k[TASMANIAN_CUDA_NUM_THREADS];
    k[threadIdx.x] = blockIdx.x*blockDim.x + threadIdx.x;
    if (k[threadIdx.x] < N){
        __shared__ double c[TASMANIAN_CUDA_NUM_THREADS];
        __shared__ int j[TASMANIAN_CUDA_NUM_THREADS];
        __shared__ int J[TASMANIAN_CUDA_NUM_THREADS];
        __shared__ int i[TASMANIAN_CUDA_NUM_THREADS];
        __shared__ int sM[TASMANIAN_CUDA_NUM_THREADS];
        sM[threadIdx.x] = M;
        i[threadIdx.x] = 0;
        while(i[threadIdx.x] < sM[threadIdx.x]){
            c[threadIdx.x] = 0.0;
            j[threadIdx.x] = gpuBpntr[i[threadIdx.x]];
            J[threadIdx.x] = gpuBpntr[i[threadIdx.x]+1];
            while(j[threadIdx.x] < J[threadIdx.x]){
                c[threadIdx.x] += gpuBvals[j[threadIdx.x]] * gpuA[gpuBindx[j[threadIdx.x]]*N + k[threadIdx.x]];
                j[threadIdx.x]++;
            }
            gpuC[(i[threadIdx.x]++)*N + k[threadIdx.x]] = c[threadIdx.x];
        }
    }
}
__global__
void tasgpu_d3gecs_2(int N, int M, int K, const double *gpuA, const int* gpuBpntr, const int *gpuBindx, const double *gpuBvals, double *gpuC, int stride){
    __shared__ int k[TASMANIAN_CUDA_NUM_THREADS];
    k[threadIdx.x] = blockIdx.x*blockDim.x + threadIdx.x;
    //__shared__ int sN = N;
    while(k[threadIdx.x] < N){
        __shared__ double c[TASMANIAN_CUDA_NUM_THREADS];
        __shared__ int j[TASMANIAN_CUDA_NUM_THREADS];
        __shared__ int J[TASMANIAN_CUDA_NUM_THREADS];
        for(int i=0; i<M; i++){
            c[threadIdx.x] = 0.0;
            j[threadIdx.x] = gpuBpntr[i];
            J[threadIdx.x] = gpuBpntr[i+1];
            while(j[threadIdx.x] < J[threadIdx.x]){
                c[threadIdx.x] += gpuBvals[j[threadIdx.x]] * gpuA[gpuBindx[j[threadIdx.x]]*N + k[threadIdx.x]];
                j[threadIdx.x]++;
            }
            gpuC[i*N + k[threadIdx.x]] = c[threadIdx.x];
        }
        k[threadIdx.x] += stride;
    }
}
 
#ifdef _TASMANIAN_DEBUG_
#define _IF_DEBUG_MACRO(x) x
#else
#define _IF_DEBUG_MACRO(x) 
#endif

TasCUDA::TasCUDA(){}
TasCUDA::~TasCUDA(){}

void TasCUDA::d3gecs(int N, int M, int K, const double *gpuA, const int *cpuBpntr, const int *cpuBindx, const double *cpuBvals, double *cpuC, std::ostream *os){
    cudaError_t cudaStat;
    int *gpuBpntr, *gpuBindx;
    double *gpuBvals;
    double *gpuC;
    
    int num_nz = cpuBpntr[M];
    
    cudaStat = cudaMalloc(((void**) &gpuBpntr), (M+1) * sizeof(int));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs alloc gpuBpntr", os);)
    cudaStat = cudaMemcpy(gpuBpntr, cpuBpntr, (M+1) * sizeof(int), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs copy gpuBpntr", os);)
    
    cudaStat = cudaMalloc(((void**) &gpuBindx), num_nz * sizeof(int));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs alloc gpuBindx", os);)
    cudaStat = cudaMemcpy(gpuBindx, cpuBindx, num_nz * sizeof(int), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs copy gpuBindx", os);)
    
    cudaStat = cudaMalloc(((void**) &gpuBvals), num_nz * sizeof(double));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs alloc gpuBvals", os);)
    cudaStat = cudaMemcpy(gpuBvals, cpuBvals, num_nz * sizeof(double), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs copy gpuBvals", os);)
    
    cudaStat = cudaMalloc(((void**) &gpuC), N * M * sizeof(double));
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs alloc gpuC", os);

    // call kernel (for now assume N < num_threads * num_blocks:  1024 * 65536 = 67108864
    int num_blocks = N / TASMANIAN_CUDA_NUM_THREADS + ((N % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks < 65536){
        tasgpu_d3gecs_1<<< num_blocks, TASMANIAN_CUDA_NUM_THREADS >>>(N, M, K, gpuA, gpuBpntr, gpuBindx, gpuBvals, gpuC);
    }else{
        tasgpu_d3gecs_2<<< num_blocks, TASMANIAN_CUDA_NUM_THREADS >>>(N, M, K, gpuA, gpuBpntr, gpuBindx, gpuBvals, gpuC, num_blocks * TASMANIAN_CUDA_NUM_THREADS);
    }

    cudaStat = cudaMemcpy(cpuC, gpuC, M * N * sizeof(double), cudaMemcpyDeviceToHost);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy back gpuC", os);

    cudaStat = cudaFree(gpuC);     _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuC", os);)
    cudaStat = cudaFree(gpuBvals); _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuBvals", os);)
    cudaStat = cudaFree(gpuBindx); _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuBindx", os);)
    cudaStat = cudaFree(gpuBpntr); _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuBpntr", os);)
}

}

#endif
