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
 
#ifndef __TASMANIAN_SPARSE_GRID_CUDA_BASIS_EVALUATIONS_HPP
#define __TASMANIAN_SPARSE_GRID_CUDA_BASIS_EVALUATIONS_HPP

#include "tasmanianConfig.hpp"

#ifdef TASMANIAN_CUDA

namespace TasGrid{

// convert a transformed domain to a canonical one
// gpu_trans_a and gpu_trans_b are the rate and shifts (precomputed on the cpu)
// size_a is the size of the gpu_trans_a and gpu_trans_b, note that the rates and shifts are repeated in some integer multiple of dims to allow contiguous access
// num_x is the number of points that need converstion
// T is the input type, i.e., the type of the transformed points
// C is the output type, i.e., the type of the canonical points and the rate and shift arrays
template<typename T, typename C, int THREADS> // transformed and canonical types
__global__ void tasgpu_transformed_to_canonical(int dims, int num_x, int size_a, const C *gpu_trans_a, const C *gpu_trans_b, const T *gpu_x_transformed, C *gpu_x_canonical){
    extern __shared__ C rate[];
    C *shift = &(rate[size_a]);

    int i = threadIdx.x;
    while(i < size_a){
        rate[i] = gpu_trans_a[i];
        shift[i] = gpu_trans_b[i];
        i += THREADS;
    }
    __syncthreads();

    int k = blockIdx.x * THREADS + threadIdx.x;
    int k_stride = THREADS * gridDim.x;
    i = k % size_a;
    while(k < (dims * num_x)){
        gpu_x_canonical[k] = ((C) gpu_x_transformed[k]) * rate[i] - shift[i];
        k += k_stride;
        i = k % size_a;
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


// evaluates the sparse grid functions for localp0, local0 and semilocalp rules using a DENSE algorithm
// this kernel will fail if the sparse grid has more than 8388608 basis functions (i.e., num_points = getNumPoints()), at this stage, we should really use the sparse algorithm (see below)
// the algorithm will also fail if the cache size: 3 * dims * num_threads * sizeof(T) exceeds the shared memory capacity (at 48K and 128 threads with 64-bit T, we can go up to 16 dimensions)
// grid blocks form a 2D grid with size: block_height x block_width
// block_height corresponds to the columns of gpu_y, i.e., the number of sparse grids nodes given by num_points = getNumPoints()
// block_width corresponds to the rows of gpu_y, i.e., num_x
// gpu_x is canonical x (same format as CPU, dim-major format)
template <typename T, int order, TypeOneDRule rule, int THREADS> // rule: localp0, localp, semilocalp
__global__ void tasgpu_devalpwpoly(int block_height, int block_width, int dims, int num_x, int num_points, const T *gpu_x, const T *gpu_nodes, const T *gpu_support, T *gpu_y){
    extern __shared__ T cache_nodes[];
    T *cache_support = &(cache_nodes[dims * THREADS]);
    T *cache_x = &(cache_support[dims * THREADS]);

    // find out how many points (nodes this block of threads must process)
    int my_point = (blockIdx.x % block_height) * THREADS; // my_point is block start point
    int local_num_points = THREADS;
    if (my_point + local_num_points > num_points)
        local_num_points = num_points - my_point;

    // cache the nodes and support values for this block of threads (never changes)
    my_point *= dims; // my_point is the offset where block points start
    for(int s=0; s<dims; s++){
        int i = s * THREADS + threadIdx.x;
        if (i < local_num_points * dims){
            cache_nodes[i] = gpu_nodes[my_point + i];
            cache_support[i] = gpu_support[my_point + i];
        }
    }
    my_point /= dims;
    my_point += threadIdx.x; // this is my sg point (or node)
    __syncthreads();

    // find out which x-points will be given to this block
    int block_x = (blockIdx.x / block_height) * THREADS; // index of the x value to be processed on this iteration (using square blocks)
    int local_num_x = THREADS;
    if (block_x + local_num_x > num_x)
        local_num_x = num_x - block_x;

    while(local_num_x > 0){
        T *local_y = &(gpu_y[block_x * num_points]);

        // cache a humber of x-points
        block_x *= dims; // my_point is the offset where block x start
        for(int s=0; s<dims; s++){
            int i = s * THREADS + threadIdx.x;
            if (i < local_num_x * dims){
                cache_x[i] = gpu_x[block_x + i];
            }
        }
        block_x /= dims;
        __syncthreads();

        if (my_point < num_points){
            for(int i=0; i<local_num_x; i++){
                T p = 1.0;
                for(int j=0; j<dims; j++){
                    T v;
                    if (order == 1) v = fabs(cache_x[i*dims + j] - cache_nodes[threadIdx.x*dims + j]);
                    if (order == 2){
                        v = (cache_x[i*dims + j] - cache_nodes[threadIdx.x*dims + j]);
                        v *= v;
                    }
                    v /= cache_support[threadIdx.x*dims + j];
                    if (v > 1.0) p = 0.0; // special points make v negative, p will not be set to 0
                    v = 1.0 - v;
                    if ((rule == rule_localp) || (rule == rule_semilocalp)) if (cache_support[threadIdx.x*dims + j] == -1.0) v = 1.0; // localp and semilocalp, point 0 (constant)
                    if ((rule == rule_localp) || (order == 2)){
                        if (cache_support[threadIdx.x*dims + j] == -2.0) v = -cache_x[i*dims + j]; // localp, point 1 (left end of the domain)
                        if (cache_support[threadIdx.x*dims + j] == -3.0) v =  cache_x[i*dims + j]; // localp, point 2 (right end of the domain)
                        if (v < 0.0) v = 0.0; // set the localp and localp0 functions to 0 outside the support
                    }
                    if ((rule == rule_semilocalp) || (order == 2)){
                        if (cache_support[threadIdx.x*dims + j] == -4.0) v = 0.5 * cache_x[i*dims + j] * (cache_x[i*dims + j] - 1.0); // semilocalp, point 1 (left end of the domain)
                        if (cache_support[threadIdx.x*dims + j] == -5.0) v = 0.5 * cache_x[i*dims + j] * (cache_x[i*dims + j] + 1.0); // semilocalp, point 2 (left end of the domain)
                    }
                    p *= v;
                }
                local_y[my_point + i*num_points] = p;
            }
        }

        block_x += block_width * THREADS + local_num_x;
        local_num_x = THREADS;
        if (block_x + local_num_x > num_x)
            local_num_x = num_x - block_x;
        __syncthreads();
    }
}


// evaluates a single basis function, for use with sparse evaluations
template <typename T, int order, TypeOneDRule rule>
__device__ inline T tasgpu_devalpwpoly_sparse_feval(int dims, int i, int ip, const T *x, const T *nodes, const T *support){
    T p = 1.0;
    for(int j=0; j<dims; j++){
        T v;
        if (rule == rule_localp){
            if (order == 1){
                v = 1.0 - fabs(x[i * dims + j] - nodes[ip * dims + j]) / support[ip * dims + j];
                if (support[ip * dims + j] == -1.0) v = 1.0;
                if (v < 0.0) v = 0.0;
            }else if (order == 2){
                v = x[i * dims + j] - nodes[ip * dims + j];
                v *= v;
                v = 1.0 - v / support[ip * dims + j];
                if (support[ip * dims + j] == -1.0) v = 1.0;
                if (support[ip * dims + j] == -2.0) v = -x[i * dims + j];
                if (support[ip * dims + j] == -3.0) v =  x[i * dims + j];
                if (v < 0.0) v = 0.0;
            }
        }else if (rule == rule_localp0){
            if (order == 1){
                v = 1.0 - fabs(x[i * dims + j] - nodes[ip * dims + j]) / support[ip * dims + j];
                if (v < 0.0) v = 0.0;
            }else if (order == 2){
                v = x[i * dims + j] - nodes[ip * dims + j];
                v *= v;
                v = 1.0 - v / support[ip * dims + j];
                if (v < 0.0) v = 0.0;
            }
        }else{
            if (order == 2){
                v = x[i * dims + j] - nodes[ip * dims + j];
                v *= v;
                v = 1.0 - v / support[ip * dims + j];
                if (v < 0.0) v = 0.0;
                if (support[ip * dims + j] == -1.0) v = 1.0;
                if (support[ip * dims + j] == -4.0) v = 0.5 * x[i * dims + j] * (x[i * dims + j] - 1.0);
                if (support[ip * dims + j] == -5.0) v = 0.5 * x[i * dims + j] * (x[i * dims + j] + 1.0);
            }
        }
        p *= v;
    }
    return p;
}

// essentially the same as local polynomial evaluate, just repeated for many threads
// at 64 threads and 46 TOPLEVEL, the 48K of shared memory is ehxausted with 4 byte integer
// note that the resolution at level 46 is 1 / 2^46 = O(1.E-13), which is comparative to the precision, i.e., never use more than 46 levels, makes no sense
// the kernel is called twice with fill false and then true
// false only counts the non-zeros, sindx and svals are ignored
// true fills the non-zeros, sindx and svals are pre-allocated and will be filled
// dense being true means that the algorithm will fill the dense matrix "dense" and ignore the sparse ones (call only with fill == false)
template <typename T, int THREADS, int TOPLEVEL, int order, TypeOneDRule rule, bool fill, bool dense>
__global__ void tasgpu_devalpwpoly_sparse(int dims, int num_x, int num_points,
                                          const T *x, const T *nodes, const T *support,
                                          const int *hpntr, const int *hindx, int num_roots, const int *roots,
                                          int *spntr, int *sindx, T *svals, T *gpu_dense = 0){
    __shared__ int mcount[TOPLEVEL][THREADS];
    __shared__ int mstop[TOPLEVEL][THREADS];

    int i = blockIdx.x * THREADS + threadIdx.x;

    while(i < num_x){
        int c = 0;

        if (fill) c = spntr[i];
        //printf("i = %d   c = %d\n", i, c);

        for(int r=0; r<num_roots; r++){
            int ip = roots[r];
            T p = tasgpu_devalpwpoly_sparse_feval<T, order, rule>(dims, i, ip, x, nodes, support);

            if (p != 0.0){
                if (fill){
                    sindx[c] = ip;
                    svals[c] = p;
                }
                c++;
                if (dense) gpu_dense[i * num_points + ip] = p;

                int current = 0;
                mstop[0][threadIdx.x] = hpntr[ip + 1];
                mcount[0][threadIdx.x] = hpntr[ip];

                while(mcount[0][threadIdx.x] < mstop[0][threadIdx.x]){
                    if (mcount[current][threadIdx.x] < mstop[current][threadIdx.x]){
                        ip = hindx[mcount[current][threadIdx.x]];
                        p = tasgpu_devalpwpoly_sparse_feval<T, order, rule>(dims, i, ip, x, nodes, support);

                        if (p != 0.0){
                            if (fill){
                                sindx[c] = ip;
                                svals[c] = p;
                            }
                            c++;
                            if (dense) gpu_dense[i * num_points + ip] = p;

                            current++;
                            mstop[current][threadIdx.x] = hpntr[ip + 1];
                            mcount[current][threadIdx.x] = hpntr[ip];
                        }else{
                            mcount[current][threadIdx.x]++;
                        }
                    }else{
                        current--;
                        mcount[current][threadIdx.x]++;
                    }
                }
            }
        }

        if (!fill && !dense){
            spntr[i+1] = c;
        }

        i += gridDim.x * THREADS;
    }
}

}

#endif

#endif
