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


//#define TASMANIAN_CUDA_NUM_CACHE 4

namespace TasGrid{

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

void TasCUDA::d3gecs(int N, int M, const double *gpuA, const int *cpuBpntr, const int *cpuBindx, const double *cpuBvals, double *cpuC, std::ostream *os){
    int *gpuBpntr, *gpuBindx;
    double *gpuBvals;
    double *gpuC;

    int num_nz = cpuBpntr[M];

    gpuBpntr = TasCUDA::cudaSend<int>(M+1, cpuBpntr, os);
    gpuBindx = TasCUDA::cudaSend<int>(num_nz, cpuBindx, os);
    gpuBvals = TasCUDA::cudaSend<double>(num_nz, cpuBvals, os);
    gpuC = TasCUDA::cudaNew<double>(M*N, os);

    int num_blocks = N / TASMANIAN_CUDA_NUM_THREADS + ((N % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_d3gecs_1_v2<<< num_blocks, TASMANIAN_CUDA_NUM_THREADS >>>(N, M, num_nz, gpuA, gpuBpntr, gpuBindx, gpuBvals, gpuC, num_blocks*TASMANIAN_CUDA_NUM_THREADS);

    TasCUDA::cudaRecv<double>(M*N, gpuC, cpuC, os);

    TasCUDA::cudaDel<double>(gpuC);
    TasCUDA::cudaDel<double>(gpuBvals);
    TasCUDA::cudaDel<int>(gpuBindx);
    TasCUDA::cudaDel<int>(gpuBpntr);
}

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
    //cudaError_t cudaStat;
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
// Double - Sparse - TRiangular - Dense Row format - Solve (DSTRDRS)

// applies the transform, but also transposes the array, not really needed in hid sight but I want to keep the transposed logic
//__global__
//void tasgpu_transformed_to_canonical(int dims, int num_x, const double *gpu_trans_a, const double *gpu_trans_b, const double *gpu_x_transformed, double *gpu_x_canonical, int k_stride){
//    extern __shared__ double rate[];
//    double *shift = &(rate[dims]);
//    double *tempx = &(shift[dims]);
//    int k = threadIdx.x;
//    while (k < dims){
//        double a = gpu_trans_a[k];
//        double b = gpu_trans_b[k];
//        double c = b - a;
//        rate[k] = 2.0 / c;
//        shift[k] = (b + a) / c;
//
//        k += blockDim.x;
//    }
//    __syncthreads();
//
//    k = blockIdx.x*blockDim.x + threadIdx.x;
//    //printf("%1.4e\n",shift[k]);
//    //printf("%d\n", k);
//
//    int offset = blockIdx.x*blockDim.x * dims;
//    int active_threads = ((k - threadIdx.x + blockDim.x) < num_x) ? blockDim.x : (num_x - k + threadIdx.x);
//    while (k < num_x){
//        // load a block of points
//        int l = threadIdx.x;
//        //printf("l = %d,  off = %d\n", l, offset);
//        for(int s=0; s<dims; s++){
//            int i = l / dims;
//            //printf("l = %d   i = %d  blah = %d  active = %d   s = %d   dims = %d\n", l, i, (l - i*dims) * active_threads + i, active_threads, s, dims);
//            tempx[(l - i*dims) * active_threads + i] = gpu_x_transformed[offset + l];
//            l += active_threads;
//        }
//        __syncthreads();
//
//        for(int s=0; s<dims; s++){
//            double x = tempx[threadIdx.x + active_threads * s] * rate[s];
//            gpu_x_canonical[num_x * s + k] = x - shift[s];
//            //printf(" target = %d,  tmpx idx = %d  temp = %1.3f\n", num_x * s + k, threadIdx.x + active_threads * s, tempx[threadIdx.x + active_threads * s]);
//        }
//
//        // figure out the active_threads after incrementing k
//
//        k += k_stride;
//        offset += blockDim.x * dims;
//        active_threads = ((k - threadIdx.x + blockDim.x) < num_x) ? blockDim.x : (num_x - k + threadIdx.x);
//        __syncthreads();
//    }
//}

__global__
void tasgpu_transformed_to_canonical(int dims, int num_x, int size_a, const double *gpu_trans_a, const double *gpu_trans_b, const double *gpu_x_transformed, double *gpu_x_canonical, int k_stride){
    extern __shared__ double rate[];
    double *shift = &(rate[size_a]);

    int i = threadIdx.x;
    while(i < size_a){
        rate[i] = gpu_trans_a[i];
        shift[i] = gpu_trans_b[i];
        i += blockDim.x;
    }
    __syncthreads();

    int k = blockIdx.x * blockDim.x + threadIdx.x;
    i = k % size_a;
    while(k < (dims * num_x)){
        gpu_x_canonical[k] = gpu_x_transformed[k] * rate[i] - shift[i];
        k += k_stride;
        i = k % size_a;
    }
}

void TasCUDA::dtrans2can(int dims, int num_x, int pad_size, const double *gpu_trans_a, const double *gpu_trans_b, const double *gpu_x_transformed, double *gpu_x_canonical){
    int num_blocks = (num_x * dims) / TASMANIAN_CUDA_NUM_THREADS + (((num_x * dims) % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_transformed_to_canonical<<<num_blocks, TASMANIAN_CUDA_NUM_THREADS, (2*pad_size) * sizeof(double)>>>(dims, num_x, pad_size, gpu_trans_a, gpu_trans_b, gpu_x_transformed, gpu_x_canonical, TASMANIAN_CUDA_NUM_THREADS * num_blocks);
}

template<typename T, int THREADS>
__global__ void tasgpu_sparse_matmul(int num_rows, int num_nz, int num_points, const int *pntr, const int *indx, const T *vals, const T *B, const T *C){
    __shared__ int cache_indx[THREADS];
    __shared__ T cache_vals[THREADS];
    int k = blockIdx.x * THREADS;
    while(k < num_rows){
        int i = k;
        int endi = k + THREADS;
        if (endi > num_rows) endi = num_rows - k;
        int jcached = pntr[i];
        while(i < k + THREADS){
            // stage 1: load into cache
            if (jcached + threadIdx.x < num_nz){
                cache_indx[threadIdx.x] = indx[jcached + threadIdx.x];
                cache_vals[threadIdx.x] = vals[jcached + threadIdx.x];
            }
            __syncthreads();
            // stage 2: loop over cached
            int endj = THREADS;
            if (jcached + endj > num_nz) endj = num_nz - jcached;
            if (jcached + endj > pntr[endi]) endj = pntr[endi] - jcached;
            int c = threadIdx.x;
            while(c < num_points){
                int sum = 0.0;
                for(int j=0; j<endj; j++){
                    if (j >= pntr[i+1]) i++;

                }
                c += THREADS;
            }
            __syncthreads();
        }

        k += gridDim.x * THREADS;
    }
}


template <int order, int rule> // rule: 0 -> localp0, 1 - localp, 2 - semi-localp
__global__ void tasgpu_devalpwpoly(int block_height, int block_width, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support, double *gpu_y){
    extern __shared__ double cache_nodes[];
    double *cache_support = &(cache_nodes[dims * blockDim.x]);
    double *cache_x = &(cache_support[dims * blockDim.x]);

    // find out how many points (nodes this block of threads must process)
    int my_point = (blockIdx.x % block_height) * blockDim.x; // my_point is block start point
    int local_num_points = blockDim.x;
    if (my_point + local_num_points > num_points)
        local_num_points = num_points - my_point;

    // cache the nodes and support values for this block of threads (never changes)
    my_point *= dims; // my_point is the offset where block points start
    for(int s=0; s<dims; s++){
        int i = s * blockDim.x + threadIdx.x;
        if (i < local_num_points * dims){
            cache_nodes[i] = gpu_nodes[my_point + i];
            cache_support[i] = gpu_support[my_point + i];
        }
    }
    my_point /= dims;
    my_point += threadIdx.x; // this is my sg point (or node)
    __syncthreads();

    // find out which x-points will be given to this block
    int block_x = (blockIdx.x / block_height) * blockDim.x; // index of the x value to be processed on this iteration (using square blocks)
    int local_num_x = blockDim.x;
    if (block_x + local_num_x > num_x)
        local_num_x = num_x - block_x;

    while(local_num_x > 0){
        double *local_y = &(gpu_y[block_x * num_points]);

        // cache a humber of x-points
        block_x *= dims; // my_point is the offset where block x start
        for(int s=0; s<dims; s++){
            int i = s * blockDim.x + threadIdx.x;
            if (i < local_num_x * dims){
                cache_x[i] = gpu_x[block_x + i];
            }
        }
        block_x /= dims;
        __syncthreads();

        if (my_point < num_points){
            for(int i=0; i<local_num_x; i++){
                double p = 1.0;
                for(int j=0; j<dims; j++){
                    double v;
                    if (order == 1) v = fabs(cache_x[i*dims + j] - cache_nodes[threadIdx.x*dims + j]);
                    if (order == 2){
                        v = (cache_x[i*dims + j] - cache_nodes[threadIdx.x*dims + j]);
                        v *= v;
                    }
                    v /= cache_support[threadIdx.x*dims + j];
                    if (v > 1.0) p = 0.0; // special points make v negative, p will not be set to 0
                    v = 1.0 - v;
                    if ((rule == 1) || (rule == 2)) if (cache_support[threadIdx.x*dims + j] == -1.0) v = 1.0; // localp and semilocalp, point 0 (constant)
                    if ((rule == 1) || (order == 2)){
                        if (cache_support[threadIdx.x*dims + j] == -2.0) v = -cache_x[i*dims + j]; // localp, point 1 (left end of the domain)
                        if (cache_support[threadIdx.x*dims + j] == -3.0) v =  cache_x[i*dims + j]; // localp, point 2 (right end of the domain)
                        if (v < 0.0) v = 0.0; // set the localp and localp0 functions to 0 outside the support
                    }
                    if ((rule == 2) || (order == 2)){
                        if (cache_support[threadIdx.x*dims + j] == -4.0) v = 0.5 * cache_x[i*dims + j] * (cache_x[i*dims + j] - 1.0); // semilocalp, point 1 (left end of the domain)
                        if (cache_support[threadIdx.x*dims + j] == -5.0) v = 0.5 * cache_x[i*dims + j] * (cache_x[i*dims + j] + 1.0); // semilocalp, point 2 (left end of the domain)
                    }
                    p *= v;
                }
                local_y[my_point + i*num_points] = p;
            }
        }

        block_x += block_width * blockDim.x + local_num_x;
        local_num_x = blockDim.x;
        if (block_x + local_num_x > num_x)
            local_num_x = num_x - block_x;
        __syncthreads();
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
            case 2: tasgpu_devalpwpoly<2, 0><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
                            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    tasgpu_devalpwpoly<1, 0><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
                            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else if (rule == rule_localp){
        switch(order){
            case 2: tasgpu_devalpwpoly<2, 1><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
                            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    tasgpu_devalpwpoly<1, 1><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
                            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else{
        tasgpu_devalpwpoly<2, 2><<<block_width * block_height, TASMANIAN_CUDA_NUM_THREADS_SHORT, 3 * dims * TASMANIAN_CUDA_NUM_THREADS_SHORT * sizeof(double)>>>
            (block_height, block_width, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
    }
}

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
    //printf("Eval device = %1.4e  %1.4e  %1.4e -- %1.4e  %1.4e  %1.4e -- p = %1.4e\n", x[i * dims + 0], nodes[ip * dims + 0], support[ip * dims + 0], x[i * dims + 1], nodes[ip * dims + 1], support[ip * dims + 1], p);
    return p;
}

template <typename T, int THREADS, int TOPLEVEL, int order, TypeOneDRule rule, bool fill> // rule: 0 - localp0, 1 - localp, 2 - semi-localp
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
        //printf("i = %d   c = %d\n", i, c);

        for(int r=0; r<num_roots; r++){
            int ip = roots[r];
            T p = tasgpu_devalpwpoly_sparse_feval<T, order, rule>(dims, i, ip, x, nodes, support);

            if (p != 0.0){
                if (fill){
                    sindx[c] = ip;
                    svals[c] = p;
                    //printf(" c = %d   ip = %d   p = %1.4e\n", c, ip, p);
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
                                //printf(" c = %d   ip = %d   p = %1.4e\n", c, ip, p);
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

        if (!fill){
            spntr[i+1] = c;
            //printf("i = %d  c = %d\n", i, c);
        }

        i += gridDim.x * THREADS;
    }
}

void TasCUDA::devalpwpoly_sparse(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                                 int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nz){
    gpu_spntr = cudaNew<int>(num_x + 1);
    int num_blocks = num_x / 64 + ((num_x % 64 == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    if (rule == rule_localp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        }
    }else if (rule == rule_semilocalp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_semilocalp, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
        break;
        }
    }else{
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp0, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp0, false><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        }
    }
    int *cpu_spntr = new int[num_x+1];
    cudaRecv(num_x+1, gpu_spntr, cpu_spntr);
    cpu_spntr[0] = 0;
    for(int i=1; i<=num_x; i++) cpu_spntr[i] += cpu_spntr[i-1];
    num_nz = cpu_spntr[num_x];
    cudaSend<int>(num_x + 1, cpu_spntr, gpu_spntr);
    gpu_sindx = cudaNew<int>(num_nz);
    gpu_svals = cudaNew<double>(num_nz);
    if (rule == rule_localp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        }
    }else if (rule == rule_semilocalp){
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_semilocalp, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
        break;
        }
    }else{
        switch(order){
        case 2:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 2, rule_localp0, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        default:
            tasgpu_devalpwpoly_sparse<double, 64, 46, 1, rule_localp0, true><<<num_blocks, 64>>>
                (dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals);
        break;
        }
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
