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

#include "TasmanianConfig.hpp"

namespace TasGrid{

// NOTE: the kernels here were used for testing and debugging (and place holders before using the MAGMA library)
//       presently, none of the templates are instantiated unless __TASMANIAN_COMPILE_FALLBACK_CUDA_KERNELS__ is defined in tsgAcceleratedDataStructures.hpp (TasCUDA namespace)

// only works with SHORT = 32 and BLOCK = 96
// using blocks, block in C has dims BLOCK by BLOCK = 3 * SHORT by 3 * SHORT,
// blocks of A and B are both BLOCK by SHORT or 3 * SHORT by SHORT
template <typename T, int SHORT, int BLOCK> // assuming ratio LONG = 2 * SHORT, THREADS = SHORT * SHORT
__global__ void tasgpu_cudaTgemm(int M, int N, int K, const T *gpu_a, const T *gpu_b, T *gpu_c){
    __shared__ T cache_a[SHORT * BLOCK];
    __shared__ T cache_b[SHORT * BLOCK];

    int num_blocks_m = M / BLOCK + ((M % BLOCK == 0) ? 0 : 1); // number of total blocks in M direction
    int num_blocks_n = N / BLOCK + ((N % BLOCK == 0) ? 0 : 1); // number of total blocks in N direction
    num_blocks_m *= num_blocks_n; // num_blocks_m is never used without a product to num_blocks_n

    int block = blockIdx.x; // keep track of the block being processed (each thread block takes on a block of C)

    int height = threadIdx.x % SHORT; // threads are oranized logivally in a square of size SHORT by SHORT (relates to M and K)
    int swidth = threadIdx.x / SHORT; // height and swidth are the indexes of this thread in the block (relates to N and K)
    // sequential threads are adjacent in height (column major logic)

    while(block < num_blocks_m){

        // find out the current matrix block for this block of threads
        int block_m = block / num_blocks_n;
        int maxM = BLOCK;
        if (block_m * BLOCK + maxM > M) maxM = M - block_m * BLOCK;

        int block_n = block % num_blocks_n;
        int maxN = BLOCK;
        if (block_n * BLOCK + maxN > N) maxN = N - block_n * BLOCK;

        int offsetm = block_m * BLOCK; // identify row/column of C

        // prepare cache in c
        T c11 = 0.0;    T c12 = 0.0;    T c13 = 0.0;
        T c21 = 0.0;    T c22 = 0.0;    T c23 = 0.0;
        T c31 = 0.0;    T c32 = 0.0;    T c33 = 0.0;

        for(int k=0; k<K; k += SHORT){
            int maxK = SHORT;
            if (k + maxK > K) maxK = K - k;

            // preload A
            if (swidth < maxK){
                int offseta  = k * M + swidth * M + offsetm + height;
                int offsetl = swidth * SHORT + height;
                if (height < maxM) cache_a[offsetl] = gpu_a[offseta];
                offseta  += SHORT;
                offsetl += SHORT*SHORT;
                if (height + SHORT < maxM) cache_a[offsetl] = gpu_a[offseta];
                offseta  += SHORT;
                offsetl += SHORT*SHORT;
                if (height + 2*SHORT < maxM) cache_a[offsetl] = gpu_a[offseta];
            }

            // preload B
            if (height < maxK){
                // start of the block + k processes so far + local column + the k for this thread
                int offsetb = block_n * BLOCK * K + k + swidth * K + height;
                if (swidth < maxN) cache_b[threadIdx.x] = gpu_b[offsetb];
                offsetb += SHORT * K;
                if (swidth + SHORT < maxN) cache_b[SHORT * SHORT + threadIdx.x] = gpu_b[offsetb];
                offsetb += SHORT * K;
                if (swidth + 2 * SHORT < maxN) cache_b[2 * SHORT * SHORT + threadIdx.x] = gpu_b[offsetb];
            }
            __syncthreads();

            for(int local_k=0; local_k < maxK; local_k++){
                // process the 4 blocks of C
                int offa = height         + local_k * SHORT; // position in a cache
                int offb = swidth * SHORT + local_k;         // position in b cache
                // for memory in cache a threads read adjacent values, in bache b all threads read one valua at a time
                T val = cache_a[offa];
                c11 += val * cache_b[offb];
                c12 += val * cache_b[offb +     SHORT * SHORT];
                c13 += val * cache_b[offb + 2 * SHORT * SHORT];

                val = cache_a[offa + SHORT * SHORT];
                c21 += val * cache_b[offb];
                c22 += val * cache_b[offb +     SHORT * SHORT];
                c23 += val * cache_b[offb + 2 * SHORT * SHORT];

                val = cache_a[offa + 2 * SHORT * SHORT];
                c31 += val * cache_b[offb];
                c32 += val * cache_b[offb +     SHORT * SHORT];
                c33 += val * cache_b[offb + 2 * SHORT * SHORT];
            }
            __syncthreads();
        }

        // write into C
        //        block column          thread column        block row         thread row
        offsetm = block_n * BLOCK * M + swidth * M + block_m * BLOCK + height;
        if (swidth < maxN){
            if (height             < maxM) gpu_c[offsetm            ] = c11;
            if (height +     SHORT < maxM) gpu_c[offsetm +     SHORT] = c21;
            if (height + 2 * SHORT < maxM) gpu_c[offsetm + 2 * SHORT] = c31;
        }
        offsetm += SHORT * M; // jump by short columns
        if (swidth + SHORT < maxN){
            if (height             < maxM) gpu_c[offsetm            ] = c12;
            if (height +     SHORT < maxM) gpu_c[offsetm +     SHORT] = c22;
            if (height + 2 * SHORT < maxM) gpu_c[offsetm + 2 * SHORT] = c32;
        }
        offsetm += SHORT * M; // jump by short columns
        if (swidth + 2 * SHORT < maxN){
            if (height             < maxM) gpu_c[offsetm            ] = c13;
            if (height +     SHORT < maxM) gpu_c[offsetm +     SHORT] = c23;
            if (height + 2 * SHORT < maxM) gpu_c[offsetm + 2 * SHORT] = c33;
        }

        block += gridDim.x;
    }
}

// sparse matrix-matrix multiply
// C = A * B, where A is in sparse row compressed form (pntr, indx, vals, num_nz) and C and B are in row major formats, C is M by N
// threads in a block first cache indx and vals, then use the cached values to march accross N incrementing C contiguously
template<typename T, int THREADS>
__global__ void tasgpu_sparse_matmul(int M, int N, int num_nz, const int *pntr, const int *indx, const T *vals, const T *B, T *C){
    __shared__ int cache_indx[THREADS];
    __shared__ T cache_vals[THREADS];
    int i = blockIdx.x * THREADS; // indexes rows of C
    while(i < M){
        int endi = i + THREADS;
        if (endi > M) endi = M;
        while(i < endi){
            int c = threadIdx.x; // indexes columns of C
            while(c < N){ // intialize C to zero
                C[i * N + c] = 0.0;
                c += THREADS;
            }
            int offj = pntr[i];
            while(offj < pntr[i+1]){
                if (offj + threadIdx.x < num_nz){ // cache a bunch of indx and vals
                    cache_indx[threadIdx.x] = indx[offj + threadIdx.x];
                    cache_vals[threadIdx.x] = vals[offj + threadIdx.x];
                }
                __syncthreads();
                int endj = THREADS;
                if (offj + endj > pntr[i+1]) endj = pntr[i+1] - offj; // stop when reaching the next i
                c = threadIdx.x;
                while(c < N){
                    T sum = 0.0;
                    for(int j=0; j<endj; j++){ // indexes non-zeros
                        sum += cache_vals[j] * B[cache_indx[j] * N + c];
                    }
                    C[i * N + c] += sum;
                    c += THREADS;
                }
                offj += endj;
                __syncthreads();
            }
            i++;
        }
        i += gridDim.x * THREADS; // move to the next set of rows
    }
}

// dense row major matrix times a sparse vector, the result is a row-major matrix
// indx, vals and num_nz descripbe the vector, A is M by N, C is M by 1, C = A * vector (M here is num_outputs)
// each thread processes up to 4 entries of A increasing reuse of the cached vals and indx
// NBLOCKS specifies how many entries
template<typename T, int THREADS, int NBLOCKS>
__global__ void tasgpu_sparse_matveci(int M, int num_nz, const T *A, const int *indx, const T *vals, T *C){
    __shared__ int cache_indx[THREADS];
    __shared__ T cache_vals[THREADS];
    int m = blockIdx.x * THREADS * NBLOCKS ; // m is the starting row of the block
    while(m < M){
        m += threadIdx.x; // the row this thread will process

        T c1 = 0.0;
        T c2 = 0.0;
        T c3 = 0.0;
        T c4 = 0.0;

        int sparse_row = 0;
        while(sparse_row < num_nz){
            int maxVals = THREADS;
            if (sparse_row + maxVals > num_nz) maxVals = num_nz - sparse_row;
            // cache the vector
            if (threadIdx.x < maxVals){
                cache_indx[threadIdx.x] = indx[sparse_row + threadIdx.x];
                cache_vals[threadIdx.x] = vals[sparse_row + threadIdx.x];
            }
            __syncthreads();

            for(int i=0; i<maxVals; i++){
                int off = cache_indx[i] * M + m;
                T val = cache_vals[i];
                c1 += cache_vals[i] * A[off];
                if (NBLOCKS > 1) c2 += val * A[off +     THREADS];
                if (NBLOCKS > 2) c3 += val * A[off + 2 * THREADS];
                if (NBLOCKS > 3) c4 += val * A[off + 3 * THREADS];
            }
            sparse_row += THREADS;
            __syncthreads();
        }

        C[m] = c1;
        if (NBLOCKS > 1) C[m +     THREADS] = c2;
        if (NBLOCKS > 2) C[m + 2 * THREADS] = c3;
        if (NBLOCKS > 3) C[m + 3 * THREADS] = c4;

        m -= threadIdx.x;
        m += THREADS * NBLOCKS; // the starting row of this block
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

// vector fill, fills the data with the given value
template<typename T, long long THREADS>
__global__ void tascuda_vfill(long long n, T data[], T value){
    long long k = blockIdx.x * THREADS + threadIdx.x;
    while(k < n){
        data[k] = value;
        k += gridDim.x * THREADS;
    }
}
// strided fill, fills the data with the given value and stride
template<typename T, long long THREADS>
__global__ void tascuda_sfill(long long n, long long stride, T data[], T value){
    long long k = blockIdx.x * THREADS + threadIdx.x;
    while(k < n){
        data[k * stride] = value;
        k += gridDim.x * THREADS;
    }
}


}

#endif
