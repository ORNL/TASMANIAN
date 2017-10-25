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

#define TASMANIAN_CUDA_NUM_THREADS 1024
//#define TASMANIAN_CUDA_NUM_CACHE 4

namespace TasGrid{

//__global__
//void tasgpu_d3gecs_1(int N, int M, const double *gpuA, const int* gpuBpntr, const int *gpuBindx, const double *gpuBvals, double *gpuC){
//    __shared__ int k[TASMANIAN_CUDA_NUM_THREADS];
//    k[threadIdx.x] = blockIdx.x*blockDim.x + threadIdx.x;
//    if (k[threadIdx.x] < N){
//        __shared__ double c[TASMANIAN_CUDA_NUM_THREADS];
//        for(int i=0; i<M; i++){ // index on num_x
//            c[threadIdx.x] = 0.0;
//            for(int j=gpuBpntr[i]; j<gpuBpntr[i+1]; j++){
//                c[threadIdx.x] += gpuBvals[j] * gpuA[gpuBindx[j]*N + k[threadIdx.x]];
//            }
//            gpuC[i*N + k[threadIdx.x]] = c[threadIdx.x];
//        }
//    }
//}

__global__
void tasgpu_d3gecs_1_v2(int N, int M, int num_nz, const double *gpuA, const int* gpuBpntr, const int *gpuBindx, const double *gpuBvals, double *gpuC, int k_stride){
    __shared__ int k[TASMANIAN_CUDA_NUM_THREADS];
    __shared__ int cache_gpuBindx[TASMANIAN_CUDA_NUM_THREADS];
    __shared__ double cache_gpuBvals[TASMANIAN_CUDA_NUM_THREADS];
    __shared__ int cache_count[TASMANIAN_CUDA_NUM_THREADS];
    __shared__ double c[TASMANIAN_CUDA_NUM_THREADS];
    cache_count[threadIdx.x] = 0;
    if (threadIdx.x < num_nz){
        cache_gpuBindx[threadIdx.x] = gpuBindx[threadIdx.x];
        cache_gpuBvals[threadIdx.x] = gpuBvals[threadIdx.x];
    }
    for(int i=0; i<M; i++){ // index on num_x
        __syncthreads();
        if ((gpuBpntr[i+1] - cache_count[threadIdx.x]) > TASMANIAN_CUDA_NUM_THREADS){
            cache_count[threadIdx.x] = gpuBpntr[i];
            if ((threadIdx.x + cache_count[threadIdx.x]) < num_nz){
                cache_gpuBindx[threadIdx.x] = gpuBindx[threadIdx.x + cache_count[threadIdx.x]];
                cache_gpuBvals[threadIdx.x] = gpuBvals[threadIdx.x + cache_count[threadIdx.x]];
            }
        }
        __syncthreads();
        k[threadIdx.x] = blockIdx.x*blockDim.x + threadIdx.x;
        while (k[threadIdx.x] < N){
            //double c = 0.0;
            c[threadIdx.x] = 0.0;
            for(int j=gpuBpntr[i]-cache_count[threadIdx.x]; j<gpuBpntr[i+1]-cache_count[threadIdx.x]; j++){
                c[threadIdx.x] += cache_gpuBvals[j] * gpuA[cache_gpuBindx[j]*N + k[threadIdx.x]];
            }
            gpuC[i*N + k[threadIdx.x]] = c[threadIdx.x];
            k[threadIdx.x] += k_stride;
        }
    }
}

//__global__
//void tasgpu_d3gecs_numx(int N, int M, const double *gpuA, const int* gpuBpntr, const int *gpuBindx, const double *gpuBvals, double *gpuC){
//    __shared__ int i[TASMANIAN_CUDA_NUM_THREADS];
//    i[threadIdx.x] = blockIdx.x*blockDim.x + threadIdx.x;
//    if (i[threadIdx.x] < M){
//        __shared__ double vals[TASMANIAN_CUDA_NUM_THREADS];
//        __shared__ int idx[TASMANIAN_CUDA_NUM_THREADS];
//        __shared__ int off[TASMANIAN_CUDA_NUM_THREADS];
//        off[threadIdx.x] = i[threadIdx.x] * N;
//        vals[threadIdx.x] = gpuBvals[gpuBpntr[i[threadIdx.x]]];
//        for(int k=0; k<N; k++) gpuC[off[threadIdx.x] + k] = vals[threadIdx.x] * gpuA[gpuBindx[gpuBpntr[i[threadIdx.x]]] * N + k];
//        for(int j=gpuBpntr[i[threadIdx.x]]+1; j<gpuBpntr[i[threadIdx.x]+1]; j++){
//            vals[threadIdx.x] = gpuBvals[j];
//            idx[threadIdx.x] = gpuBindx[j] * N;
//            for(int k=0; k<N; k++)
//                gpuC[off[threadIdx.x] + k] += vals[threadIdx.x] * gpuA[idx[threadIdx.x] + k];
//        }
//    }
//}
//
//__global__
//void tasgpu_d3gecs_2(int N, int M, const double *gpuA, const int* gpuBpntr, const int *gpuBindx, const double *gpuBvals, double *gpuC, int stride){
//    __shared__ int k[TASMANIAN_CUDA_NUM_THREADS];
//    k[threadIdx.x] = blockIdx.x*blockDim.x + threadIdx.x;
//    __shared__ double c[TASMANIAN_CUDA_NUM_THREADS];
//    while(k[threadIdx.x] < N){
//        for(int i = 0; i<M; i++){ // index on num_x
//            c[threadIdx.x] = 0.0;
//            for(int j=gpuBpntr[i]; j<gpuBpntr[j+1]; j++){
//                c[threadIdx.x] += gpuBvals[j] * gpuA[gpuBindx[j]*N + k[threadIdx.x]];
//            }
//            gpuC[i*N + k[threadIdx.x]] = c[threadIdx.x];
//        }
//        k[threadIdx.x] += stride;
//    }
//}


 
#ifdef _TASMANIAN_DEBUG_
#define _IF_DEBUG_MACRO(x) x
#else
#define _IF_DEBUG_MACRO(x) 
#endif

TasCUDA::TasCUDA(){}
TasCUDA::~TasCUDA(){}

void TasCUDA::d3gecs(int N, int M, const double *gpuA, const int *cpuBpntr, const int *cpuBindx, const double *cpuBvals, double *cpuC, std::ostream *os){
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

    int num_blocks = N / TASMANIAN_CUDA_NUM_THREADS + ((N % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_d3gecs_1_v2<<< num_blocks, TASMANIAN_CUDA_NUM_THREADS >>>(N, M, num_nz, gpuA, gpuBpntr, gpuBindx, gpuBvals, gpuC, num_blocks*TASMANIAN_CUDA_NUM_THREADS);

    cudaStat = cudaMemcpy(cpuC, gpuC, M * N * sizeof(double), cudaMemcpyDeviceToHost);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy back gpuC", os);

    cudaStat = cudaFree(gpuC);     _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuC", os);)
    cudaStat = cudaFree(gpuBvals); _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuBvals", os);)
    cudaStat = cudaFree(gpuBindx); _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuBindx", os);)
    cudaStat = cudaFree(gpuBpntr); _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuBpntr", os);)
}

__global__
void tasgpu_d3gecss(int N, int M, int *order, int top_level, const int *gpuBpntr, const int *gpuBindx, const double *gpuBvals, double *gpuX, int k_stride){
    __shared__ int k[TASMANIAN_CUDA_NUM_THREADS];
    k[threadIdx.x] = blockIdx.x*blockDim.x + threadIdx.x;
    if (k[threadIdx.x] < N){
        for(int o=0; o<M; o++){
            int i = order[o];
            double sum = 0.0;
            for(int j=gpuBpntr[i]; j<gpuBpntr[i+1]; j++){
                sum += gpuBvals[j] * gpuX[gpuBindx[j] * N + k[threadIdx.x]];
            }
            gpuX[i*N + k[threadIdx.x]] -= sum;
        }
        k[threadIdx.x] += k_stride;
    }
}

void TasCUDA::d3gecss(int N, int M, int *level, int top_level, const int *cpuBpntr, const int *cpuBindx, const double *cpuBvals, const double *cpuA, double *cpuC, std::ostream *os){
    cudaError_t cudaStat;
    int *gpuBpntr, *gpuBindx, *gpu_order;
    double *gpuBvals;
    double *gpuX;
    
    //cout << "Calling d3gecss" << endl;
    
    int num_nz = cpuBpntr[M];
    
    int *order = new int[M], c = 0;
    for(int l=1; l<=top_level; l++){
        for(int i=0; i<M; i++){
            if (level[i] == l){
                order[c++] = i;
            }
        }
    }
    while(c < M) order[c++] = 0;
    
    cudaStat = cudaMalloc(((void**) &gpu_order), M * sizeof(int));
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs alloc levels", os);)
    cudaStat = cudaMemcpy(gpu_order, order, M * sizeof(int), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs copy levels", os);)
    
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
    
    cudaStat = cudaMalloc(((void**) &gpuX), N * M * sizeof(double));
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs alloc gpuX", os);
    cudaStat = cudaMemcpy(gpuX, cpuA, N * M * sizeof(double), cudaMemcpyHostToDevice);
    _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs copy gpuX", os);)
    
    int num_blocks = N / TASMANIAN_CUDA_NUM_THREADS + ((N % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_d3gecss<<< num_blocks, TASMANIAN_CUDA_NUM_THREADS >>>(N, M, gpu_order, top_level, gpuBpntr, gpuBindx, gpuBvals, gpuX, num_blocks * TASMANIAN_CUDA_NUM_THREADS);
    
    cudaStat = cudaMemcpy(cpuC, gpuX, M * N * sizeof(double), cudaMemcpyDeviceToHost);
    AccelerationMeta::cudaCheckError((void*) &cudaStat, "copy back gpuX", os);

    delete[] order;
    cudaStat = cudaFree(gpuX);      _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuC", os);)
    cudaStat = cudaFree(gpuBvals);  _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuBvals", os);)
    cudaStat = cudaFree(gpuBindx);  _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuBindx", os);)
    cudaStat = cudaFree(gpuBpntr);  _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpuBpntr", os);)
    cudaStat = cudaFree(gpu_order); _IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "free gpu_order", os);)
}
// Double - Sparse - TRiangular - Dense Row format - Solve (DSTRDRS)

        //int maxc = TASMANIAN_CUDA_NUM_THREADS * TASMANIAN_CUDA_NUM_CACHE;
        //int np = 1, c = 0;
        //for(int i=0; i<M; i++) if (cpuBpntr[i+1] > c + maxc){ np++; c = cpuBpntr[i]; }
        //int *pntri = new int[np], *gpu_pntri;
        //pntri[0] = 0;
        //c = 0;
        //np = 0;
        //for(int i=0; i<M; i++) if (cpuBpntr[i+1] > c + maxc){ np++; pntri[np] = i; c = cpuBpntr[i]; }
        //np++;
        //pntri[np] = M;
        ////for(int i=0; i<=M; i++) cout << cpuBpntr[i] << endl;
        ////cout << "maxc = " << maxc << "  " << N << "  " << M << "  " << num_nz << endl;
        ////for(int i=0; i<=np; i++) cout << pntri[i] << endl;
        //cudaStat = cudaMalloc(((void**) &gpu_pntri), (np+1) * sizeof(int));
        //_IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs alloc pntri", os);)
        //cudaStat = cudaMemcpy(gpu_pntri, pntri, (np+1) * sizeof(int), cudaMemcpyHostToDevice);
        //_IF_DEBUG_MACRO(AccelerationMeta::cudaCheckError((void*) &cudaStat, "d3gecs copy pntr", os);)
        //int num_blocks = N / TASMANIAN_CUDA_NUM_THREADS + ((N % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
        //if (num_blocks > 32) num_blocks = 32;
        //tasgpu_d3gecs_cached<<< num_blocks, TASMANIAN_CUDA_NUM_THREADS >>>(N, M, num_nz, gpuA, gpuBpntr, gpuBindx, gpuBvals, gpuC, np, gpu_pntri, num_blocks*TASMANIAN_CUDA_NUM_THREADS);

//__global__
//void tasgpu_d3gecs_cached(int N, int M, int num_nz, const double *gpuA, const int* gpuBpntr, const int *gpuBindx, const double *gpuBvals, double *gpuC, int num_ijumps, int *pntri, int k_stride){
//    //__shared__ int k[TASMANIAN_CUDA_NUM_THREADS];
//    __shared__ int cache_gpuBindx[TASMANIAN_CUDA_NUM_THREADS * TASMANIAN_CUDA_NUM_CACHE];
//    __shared__ double cache_gpuBvals[TASMANIAN_CUDA_NUM_THREADS * TASMANIAN_CUDA_NUM_CACHE];
//    __shared__ int cache_count[TASMANIAN_CUDA_NUM_THREADS];
//    //cache_count[threadIdx.x] = - gpuBpntr[1] - TASMANIAN_CUDA_NUM_THREADS * TASMANIAN_CUDA_NUM_CACHE;
//    //cache_count[threadIdx.x] = 0;
//    for(int c=0; c<num_ijumps; c++){
//        __syncthreads();
//        //int Nc = (c + TASMANIAN_CUDA_NUM_CACHE < M) ? c + TASMANIAN_CUDA_NUM_CACHE : M;
//        //if ((gpuBpntr[Nc] - cache_count[threadIdx.x]) > TASMANIAN_CUDA_NUM_THREADS * TASMANIAN_CUDA_NUM_CACHE){
//        cache_count[threadIdx.x] = gpuBpntr[pntri[c]];
//        //if (threadIdx.x == 0 && blockIdx.x == 0) printf("Caching ccount %d\n",cache_count[threadIdx.x]);
//        //if (threadIdx.x == 0 && blockIdx.x == 1) printf("Caching count %d, numnz = %d  %d\n",cache_count[threadIdx.x], num_nz, num_ijumps);
//        for(int cn=0; cn<TASMANIAN_CUDA_NUM_CACHE; cn++){
//            if ((threadIdx.x + cache_count[threadIdx.x] + cn * TASMANIAN_CUDA_NUM_THREADS) < num_nz){
//                cache_gpuBindx[threadIdx.x + cn * TASMANIAN_CUDA_NUM_THREADS] = gpuBindx[threadIdx.x + cache_count[threadIdx.x] + cn * TASMANIAN_CUDA_NUM_THREADS];
//                cache_gpuBvals[threadIdx.x + cn * TASMANIAN_CUDA_NUM_THREADS] = gpuBvals[threadIdx.x + cache_count[threadIdx.x] + cn * TASMANIAN_CUDA_NUM_THREADS];
//                //if (threadIdx.x == 0 && blockIdx.x == 0) printf("Caching indx %d\n",cache_gpuBindx[threadIdx.x + cn * TASMANIAN_CUDA_NUM_THREADS]);
//            }
//        }
//        //}
//        __syncthreads();
//        //if (threadIdx.x == 0 && blockIdx.x == 1) printf("Print  %d  %d  %d\n", pntri[c], pntri[c+1], c);
//        for(int i=pntri[c]; i<pntri[c+1]; i++){ // index on num_x
//            //if (threadIdx.x == 0 && blockIdx.x == 0) printf("consider i = %d\n",i);
//            //k[threadIdx.x] = blockIdx.x*blockDim.x + threadIdx.x;
//            //while(k[threadIdx.x] < N){
//            for(int kk=blockIdx.x*blockDim.x + threadIdx.x; kk<N; kk += k_stride){
//                double sum = 0.0;
//                for(int j=gpuBpntr[i]; j<gpuBpntr[i+1]; j++){
//                    //sum += cache_gpuBvals[j-cache_count[threadIdx.x]] * gpuA[cache_gpuBindx[j-cache_count[threadIdx.x]]*N + k[threadIdx.x]];
//                    sum += cache_gpuBvals[j-cache_count[threadIdx.x]] * gpuA[cache_gpuBindx[j-cache_count[threadIdx.x]]*N + kk];
//                }
//                //gpuC[i*N + k[threadIdx.x]] = sum;
//                //k[threadIdx.x] += k_stride;
//                gpuC[i*N + kk] = sum;
//            }
//        }
//    }
//}

}

#endif
