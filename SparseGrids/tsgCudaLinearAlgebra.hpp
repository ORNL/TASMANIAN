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

#ifndef __TASMANIAN_SPARSE_GRID_CUDA_LINEAR_ALGEBRA_HPP
#define __TASMANIAN_SPARSE_GRID_CUDA_LINEAR_ALGEBRA_HPP

#include "tasmanianConfig.hpp"

#ifdef TASMANIAN_CUDA

namespace TasGrid{

// sparse triangular solve using cuda kernel, does not out-perform the CPU since the values have to be moved from host to device memory
// sparse solves are also inefficient when using cuSparse
// double precision - triangular - general - column sparse - solve = d3gecss, probably too cryptic
__global__
void tasgpu_d3gecss(int N, int M, int *order, int top_level, const int *gpuBpntr, const int *gpuBindx, const double *gpuBvals, double *gpuX, int k_stride){
    int k = blockIdx.x*blockDim.x + threadIdx.x;
    if (k < N){
        for(int o=0; o<M; o++){
            int i = order[o];
            double sum = 0.0;
            for(int j=gpuBpntr[i]; j<gpuBpntr[i+1]; j++){
                sum += gpuBvals[j] * gpuX[gpuBindx[j] * N + k];
            }
            gpuX[i*N + k] -= sum;
        }
        k += k_stride;
    }
}

// really only works with 64 threads with double, but has to be launched with 2xTHREADS  size(T) * THREADS^2 + 4 * size(T) * THREADS shared memory
template <typename T, int THREADS>
__global__ void tasgpu_cudaTgemm_v3(int M, int N, int K, const T *gpu_a, const T *gpu_b, T *gpu_c){
    __shared__ T cache_c[THREADS * THREADS];
    __shared__ T a_flip[THREADS];
    __shared__ T a_flop[THREADS];
    __shared__ T b_flip[THREADS];
    __shared__ T b_flop[THREADS];

    int num_blocks_m = M / THREADS + ((M % THREADS == 0) ? 0 : 1);
    int num_blocks_n = N / THREADS + ((N % THREADS == 0) ? 0 : 1);

    int block = blockIdx.x;
    int half_thread = threadIdx.x - THREADS;
    int half_threadK = half_thread * K;

    while(block < num_blocks_m * num_blocks_n){

        // find out the current matrix block for this block of threads
        int block_m = block / num_blocks_n;
        int maxM = THREADS;
        if (block_m * THREADS + maxM > M) maxM = M - block_m * THREADS;

        int block_n = block % num_blocks_n;
        int maxN = THREADS;
        if (block_n * THREADS + maxN > N) maxN = N - block_n * THREADS;

        int offsetm = block_m * THREADS; // identify row/column of C
        int offsetn = block_n * THREADS * K;

        // prepare cache
        if (threadIdx.x < THREADS){
            for(int i=0; i<THREADS; i++){
                cache_c[i * THREADS + threadIdx.x] = 0.0;
            }
        }else{
            if (half_thread < maxM)
                a_flip[half_thread] = gpu_a[offsetm + half_thread];
            if (half_thread < maxN)
                b_flip[half_thread] = gpu_b[offsetn + half_threadK];
                //printf(" b_flip[%d] = gpu_b[%d + %d] = %1.6e\n", half_thread, offsetn, half_threadK, b_flip[half_thread]);
        }
        __syncthreads();

        // loop over K
        for(int k=0; k<K; k++){
            if (threadIdx.x < THREADS){ // these threads are computing
                T *cached_a;
                T *cached_b;
                if (k % 2 == 0){
                    cached_a = a_flip;
                    cached_b = b_flip;
                }else{
                    cached_a = a_flop;
                    cached_b = b_flop;
                }
                if (threadIdx.x < maxM){
                    T val = cached_a[threadIdx.x];
                    for(int n = 0; n < maxN; n++){
                        cache_c[n * THREADS + threadIdx.x] += val * cached_b[n];
                        //printf(" m = %d  n = %d  adding = %1.4e * %1.4e\n", threadIdx.x, n, val, cached_b[n]);
                    }
                }
            }else{ // these threads are caching the next row/columns
                offsetm += M; // move to next column of a
                offsetn += 1; // move to next row of b

                if (k+1 < K){
                    T *cached_a;
                    T *cached_b;
                    if (k % 2 == 0){
                        cached_a = a_flop;
                        cached_b = b_flop;
                    }else{
                        cached_a = a_flip;
                        cached_b = b_flip;
                    }
                    if (half_thread < maxM)
                        cached_a[half_thread] = gpu_a[offsetm + half_thread];
                    if (half_thread < maxN)
                        cached_b[half_thread] = gpu_b[offsetn + half_threadK];
                }
            }
            __syncthreads();
        }

        offsetm = block_n * THREADS * M + block_m * THREADS;

        // everyone should write over to c
        if (threadIdx.x < THREADS){
            if (threadIdx.x < maxM){
                for(int n=0; n<maxN; n+=2){
                    gpu_c[offsetm + n * M + threadIdx.x] = cache_c[n * THREADS + threadIdx.x];
                }
            }
        }else{
            if (half_thread < maxM){
                for(int n=1; n<maxN; n+=2){
                    gpu_c[offsetm + n * M + half_thread] = cache_c[n * THREADS + half_thread];
                }
            }
        }
        __syncthreads();

        block += gridDim.x;
    }


}

// simple, but not efficient (can use with 64 threads, CN = 48, and CK = 61, which gives exactly 48K of shared memory for T = double)
template <typename T, int THREADS, int CN, int CK>
__global__ void tasgpu_cudaTgemm_v2(int M, int N, int K, const T *gpu_a, const T *gpu_b, T *gpu_c){
    __shared__ T cache_c[THREADS * CN];
    __shared__ T cache_b[CN * CK];

    int m = blockIdx.x * THREADS + threadIdx.x;

    int num_blockn = N / CN;
    if (N % CN > 0) num_blockn++;
    int num_blockk = K / CK;
    if (K % CK > 0) num_blockk++;

    for(int blockn = 0; blockn < num_blockn; blockn++){
        // cache C
        for(int n=0; n<CN; n++){
            cache_c[n * THREADS + threadIdx.x] = 0.0;
        }

        int maxN = N - blockn * CN;
        if (maxN > CN) maxN = CN;

        for(int blockk = 0; blockk < num_blockk; blockk++){
            // cache B
            int maxK = K - blockk * CK;
            if (maxK > CK) maxK = CK;

            if (threadIdx.x < maxK){
                int offk = blockn * CN * K + blockk * CK + threadIdx.x;
                for(int n=0; n<maxN; n++){
                    cache_b[n * CK + threadIdx.x] = gpu_b[n * K + offk];
                }
            }
            __syncthreads();

            if (m < M){
                int off = blockk * CK * M + m;
                for(int k=0; k<maxK; k++){
                    T a = gpu_a[off + k * M];
                    for(int n=0; n<maxN; n++){
                        cache_c[n * THREADS + threadIdx.x] += a * cache_b[n * CK + k];
                    }
                }

                off = blockn * CN * M + m;
                for(int n=0; n<maxN; n++){
                    gpu_c[off + n * M] = cache_c[n * THREADS + threadIdx.x];
                }
            }
            __syncthreads();
        }
    }
}

// really only works with 32 threads with double since we need size(T) * THREADS^2 * 3 bytes of shared memory
template <typename T, int THREADS>
__global__ void tasgpu_cudaTgemm(int num_cblocks_m, int num_cblocks_n, int M, int N, int K, const T *gpu_a, const T *gpu_b, T *gpu_c){
    __shared__ T cache_a[THREADS * THREADS];
    __shared__ T cache_b[THREADS * THREADS];
    __shared__ T cache_c[THREADS * THREADS];

    int num_cblocks_k = K / THREADS;
    if (num_cblocks_k * THREADS < K) num_cblocks_k++;

    int cblock_m = blockIdx.x / num_cblocks_n;
    int maxM = THREADS;
    if (cblock_m * THREADS + maxM > M) maxM = M - cblock_m * THREADS;

    while(maxM > 0){

        int cblock_n = blockIdx.x % num_cblocks_n;
        int maxN = THREADS;
        if (cblock_n * THREADS + maxN > N) maxN = N - cblock_n * THREADS;

        while(maxN > 0){

            // zero the c-block
            for(int i=0; i<THREADS; i++){
                cache_c[i * THREADS + threadIdx.x] = 0.0;
            }
            __syncthreads();

            for(int cblock_k=0; cblock_k < num_cblocks_k; cblock_k++){
                // cache blocks of A and B
                int maxK = THREADS;
                if (cblock_k * THREADS + maxK > K) maxK = K - cblock_k * THREADS;

                if (threadIdx.x < maxM){
                    int a_row_offset = cblock_m * THREADS + cblock_k * THREADS * M + threadIdx.x;
                    for(int i=0; i<maxK; i++){
                        cache_a[i * THREADS + threadIdx.x] = gpu_a[i * M + a_row_offset];
                    }
                }

                if (threadIdx.x < maxK){
                    int b_row_offset = cblock_n * THREADS * K + cblock_k * THREADS + threadIdx.x;
                    for(int i=0; i<maxN; i++){
                        cache_b[i * THREADS + threadIdx.x] = gpu_b[i * K + b_row_offset];
                    }
                }
                __syncthreads();

                if (threadIdx.x < maxM){
                    for(int i=0; i<maxN; i++){
                        T sum = 0.0;
                        for(int j=0; j<maxK; j++){
                            sum += cache_a[j * THREADS + threadIdx.x] * cache_b[i * THREADS + j];
                        }
                        cache_c[i * THREADS + threadIdx.x] += sum;
                    }
                }

                __syncthreads();
            }

            if (threadIdx.x < maxM){
                int offset = cblock_n * THREADS * M + cblock_m * THREADS +  threadIdx.x;
                for(int i=0; i<maxN; i++){
                    gpu_c[i * M + offset] = cache_c[i * THREADS + threadIdx.x];
                }
            }
            __syncthreads();

            cblock_n += num_cblocks_n;
            maxN = THREADS;
            if (cblock_n * THREADS + maxN > N) maxN = N - cblock_n * THREADS;
        }

        cblock_m += num_cblocks_m;
        maxM = THREADS;
        if (cblock_m * THREADS + maxM > M) maxM = M - cblock_m * THREADS;
    }
}

// convert sparse matrix to dense format
template <typename T, int THREADS>
__global__ void tascuda_fill(int n, T value, T *destination){
    int k = blockIdx.x * THREADS + threadIdx.x;
    T v = value;
    while(k < n){
        destination[k] = v;
        k += gridDim.x * THREADS;
    }
}

// assumes row-compressed and row-major format
// each block processes one row, so use few threads per block
// this ensures contiguous read from indx and vals
template <typename T, int THREADS>
__global__ void tascuda_sparse_to_dense(int num_rows, int num_columns, const int *pntr, const int *indx, const T *vals, T *destination){
    int i = blockIdx.x; // row to process
    while(i < num_rows){ // each block processes one row
        int j = pntr[i];
        int endj = pntr[i+1];
        int offi = i * num_columns;
        while(j < endj){
            if (j + threadIdx.x < endj){
                destination[offi + indx[j + threadIdx.x]] = vals[j + threadIdx.x];
            }
            j += THREADS;
        }
        i += gridDim.x;
    }
}

}

#endif

#endif
