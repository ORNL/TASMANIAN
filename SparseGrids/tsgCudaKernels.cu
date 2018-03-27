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

#include "tsgCudaMacros.hpp"


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
void TasCUDA::d3gecss(int N, int M, int *level, int top_level, const int *cpuBpntr, const int *cpuBindx, const double *cpuBvals, const double *cpuA, double *cpuC, std::ostream *os){
    int *gpuBpntr, *gpuBindx, *gpu_order;
    double *gpuBvals;
    double *gpuX;
    
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
    
    gpu_order = TasCUDA::cudaSend<int>(M, order, os);
    gpuBpntr = TasCUDA::cudaSend<int>(M+1, cpuBpntr, os);
    gpuBindx = TasCUDA::cudaSend<int>(num_nz, cpuBindx, os);
    gpuBvals = TasCUDA::cudaSend<double>(num_nz, cpuBvals, os);
    gpuX = TasCUDA::cudaSend<double>(M*N, cpuA, os);

    int num_blocks = N / TASMANIAN_CUDA_NUM_THREADS + ((N % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_d3gecss<<< num_blocks, TASMANIAN_CUDA_NUM_THREADS >>>(N, M, gpu_order, top_level, gpuBpntr, gpuBindx, gpuBvals, gpuX, num_blocks * TASMANIAN_CUDA_NUM_THREADS);

    TasCUDA::cudaRecv<double>(M*N, gpuX, cpuC, os);

    delete[] order;
    TasCUDA::cudaDel<double>(gpuX);
    TasCUDA::cudaDel<double>(gpuBvals);
    TasCUDA::cudaDel<int>(gpuBindx);
    TasCUDA::cudaDel<int>(gpuBpntr);
    TasCUDA::cudaDel<int>(gpu_order);
}


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

void TasCUDA::dtrans2can(int dims, int num_x, int pad_size, const double *gpu_trans_a, const double *gpu_trans_b, const double *gpu_x_transformed, double *gpu_x_canonical){
    int num_blocks = (num_x * dims) / TASMANIAN_CUDA_NUM_THREADS + (((num_x * dims) % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_transformed_to_canonical<double, double, TASMANIAN_CUDA_NUM_THREADS><<<num_blocks, TASMANIAN_CUDA_NUM_THREADS, (2*pad_size) * sizeof(double)>>>(dims, num_x, pad_size, gpu_trans_a, gpu_trans_b, gpu_x_transformed, gpu_x_canonical);
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

void TasCUDA::cudaSparseMatmul(int M, int N, int num_nz, const int* gpu_spntr, const int* gpu_sindx, const double* gpu_svals, const double *gpu_B, double *gpu_C){
    int blocks = M / 64 + ((M % 64 == 0) ? 0 : 1);
    tasgpu_sparse_matmul<double, 64><<<blocks, 64>>>(M, N, num_nz, gpu_spntr, gpu_sindx, gpu_svals, gpu_B, gpu_C);
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

template <typename T, int order, int rule>
__device__ inline T tasgpu_devalpwpoly_sparse_feval(int dims, int i, int ip, const T *x, const T *nodes, const T *support){
    T p = 1.0;
    for(int j=0; j<dims; j++){
        T v;
        if (rule == rule_localp){
            if (order == 1){
                v = fabs(x[i * dims + j] - nodes[ip * dims + j]) / support[ip * dims + j];
                if (support[ip * dims + j] == -1.0) v = 1.0;
                if (v < 0.0) v = 0.0;
            }
        }else if (rule == rule_localp0){
        }else{
        }
        p *= v;
    }
    return p;
}

template <typename T, int THREADS, int TOPLEVEL, int order, int rule, bool fill> // rule: 0 -> localp0, 1 - localp, 2 - semi-localp>
__global__ void tasgpu_devalpwpoly_sparse(int dims, int num_x, int num_points, 
                                          const T *x, const T *nodes, const T *support, 
                                          const int *hpntr, const int *hindx, int num_roots, const int *roots,
                                          int *spntr, int *sindx, T *svals){
    __shared__ int mcount[TOPLEVEL][THREADS];
    __shared__ int mstop[TOPLEVEL][THREADS];
    
    int i = blockIdx.x * THREADS + threadIdx.x;
    
    while(i < num_x){
        int c = 0;
        
        if (fill) c = spntr[i];
        
        for(int r=0; r<num_roots; r++){
            int ip = roots[r];
            T p = tasgpu_devalpwpoly_sparse_feval<T, order, rule>(dims, i, ip, x, nodes, support);
            
            if (p != 0.0){
                if (fill){
                    sindx[c] = ip;
                    svals[c] = p;
                }
                c++;
                
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
        
        if (!fill) spntr[i] = c;
    
        i += gridDim.x * THREADS;
    }
}

void TasCUDA::devalpwpoly(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support, double *gpu_y){
    int block_height = num_points / TASMANIAN_CUDA_NUM_THREADS_SHORT + ((num_points % TASMANIAN_CUDA_NUM_THREADS_SHORT == 0) ? 0 : 1);
    // size of __shared__ 3 * dims * num_threads
    // num_blocks has to be a multiple of block_stride
    int block_width = num_x / TASMANIAN_CUDA_NUM_THREADS_SHORT + ((num_x % TASMANIAN_CUDA_NUM_THREADS_SHORT == 0) ? 0 : 1);
    if (block_height * block_width >= 65536) block_width = 65536 / block_height;
    if (rule == rule_localp0){
        switch(order){
            case 2: tasgpu_devalpwpoly<double, 2, rule_localp0, TASMANIAN_CUDA_NUM_THREADS_SHORT><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
                            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    tasgpu_devalpwpoly<double, 1, rule_localp0, TASMANIAN_CUDA_NUM_THREADS_SHORT><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
                            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else if (rule == rule_localp){
        switch(order){
            case 2: tasgpu_devalpwpoly<double, 2, rule_localp, TASMANIAN_CUDA_NUM_THREADS_SHORT><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
                            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    tasgpu_devalpwpoly<double, 1, rule_localp, TASMANIAN_CUDA_NUM_THREADS_SHORT><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
                            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else{
        tasgpu_devalpwpoly<double, 2, rule_semilocalp, TASMANIAN_CUDA_NUM_THREADS_SHORT><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
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

void TasCUDA::devalpwpoly_sparse(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                                 int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nz){
    gpu_spntr = cudaNew<int>(num_x + 1);
    int num_blocks = num_x / 64 + ((num_x % 64 == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    // call with fill == false to count the non-zeros per row of the matrix
    if (rule == rule_localp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp, false, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp, false, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        }
    }else if (rule == rule_semilocalp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_semilocalp, false, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
        break;
        }
    }else{
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp0, false, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp0, false, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        }
    }
    int *cpu_spntr = new int[num_x+1];
    cudaRecv(num_x+1, gpu_spntr, cpu_spntr);
    cpu_spntr[0] = 0;
    for(int i=1; i<=num_x; i++) cpu_spntr[i] += cpu_spntr[i-1];
    num_nz = cpu_spntr[num_x]; // save the number of non-zeros
    cudaSend<int>(num_x + 1, cpu_spntr, gpu_spntr);
    delete[] cpu_spntr;
    gpu_sindx = cudaNew<int>(num_nz);
    gpu_svals = cudaNew<double>(num_nz);
    // call with fill == true to load the non-zeros
    if (rule == rule_localp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp, true, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp, true, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        }
    }else if (rule == rule_semilocalp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_semilocalp, true, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
        break;
        }
    }else{
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp0, true, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp0, true, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        }
    }
}

void TasCUDA::devalpwpoly_sparse_dense(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                                 int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, double *gpu_dense){
    int num_blocks = num_x / 64 + ((num_x % 64 == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    if (rule == rule_localp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp, false, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, 0, 0, 0, gpu_dense);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp, false, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, 0, 0, 0, gpu_dense);
        break;
        }
    }else if (rule == rule_semilocalp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_semilocalp, false, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, 0, 0, 0, gpu_dense);
        break;
        default:
        break;
        }
    }else{
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp0, false, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, 0, 0, 0, gpu_dense);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp0, false, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, 0, 0, 0, gpu_dense);
        break;
        }
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

void TasCUDA::convert_sparse_to_dense(int num_rows, int num_columns, const int *pntr, const int *indx, const double *vals, double *destination){
    int n = num_rows * num_columns;
    int num_blocks = n / TASMANIAN_CUDA_NUM_THREADS + ((n % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tascuda_fill<double, TASMANIAN_CUDA_NUM_THREADS><<<num_blocks, TASMANIAN_CUDA_NUM_THREADS>>>(n, 0.0, destination);
    num_blocks = num_rows;
    if (num_blocks >= 65536) num_blocks = 65536;
    tascuda_sparse_to_dense<double, 64><<<num_blocks, 64>>>(num_rows, num_columns, pntr, indx, vals, destination);
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

void TasCUDA::cudaDgemm(int M, int N, int K, const double *gpu_a, const double *gpu_b, double *gpu_c){ // gpu_c = gpu_a * gpu_b, gpu_c is M by N
    const int num_threads = 32;

    int activeN = N;
    int blocks_n = (activeN / num_threads) + (((activeN % num_threads) == 0) ? 0 : 1);
    while(blocks_n > 65536){
        activeN /= 2;
        blocks_n = (activeN / num_threads) + (((activeN % num_threads) == 0) ? 0 : 1);
    }
    int blocks_m = 1;
    while((blocks_m * num_threads < M) && ((blocks_m + 1) * blocks_n < 65536)) blocks_m++;
    tasgpu_cudaTgemm<double, 32><<<blocks_m * blocks_n, num_threads>>>(blocks_m, blocks_n, M, N, K, gpu_a, gpu_b, gpu_c);

    // check the last error
    //cudaError_t cudaStat = cudaPeekAtLastError();
    //AccelerationMeta::cudaCheckError((void*) &cudaStat, "kernel()", &cerr);

    // pull out gpu_a and gpu_b for debugging purposes
    //double *A = new double[M*K]; cudaRecv<double>(M * K, gpu_a, A);
    //double *B = new double[N*K]; cudaRecv<double>(N * K, gpu_b, B);
    //double *C = new double[N*M];

    // print gpu_a and gpu_b
    //cout << std::scientific; cout.precision(16);
    //for(int i=0; i<1; i++){ for(int j=0; j<K; j++){ cout << A[j*M + i] << "   "; } cout << endl; }
    //cout << endl;
    //for(int i=0; i<K; i++){ for(int j=0; j<1; j++){ cout << B[j*K + i] << "   "; } cout << endl; }
}

}

#endif
