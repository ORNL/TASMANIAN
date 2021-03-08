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

#include "TasmanianConfig.hpp"

#ifdef Tasmanian_ENABLE_CUDA

#include <cuda.h>
#include <cuda_runtime_api.h>
// Global grids require atomic add
#include <sm_20_atomic_functions.h>
#if __CUDA_ARCH__ >= 600
#include <sm_60_atomic_functions.h> // atomic add for doubles was added in arch 6.0
#else
#ifdef __CUDA_ARCH__
// borrowed from https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
// implements less efficient (software) atomic add that works on archs before 6.0
__device__ inline double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif
#endif
#endif

#ifdef Tasmanian_ENABLE_HIP
#ifndef __HIP_PLATFORM_HCC__
#define __HIP_PLATFORM_HCC__
#endif
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#endif


namespace TasGrid{
// convert a transformed domain to a canonical one
// gpu_trans_a and gpu_trans_b are the rate and shifts (precomputed on the cpu)
// size_a is the size of the gpu_trans_a and gpu_trans_b, note that the rates and shifts are repeated in some integer multiple of dims to allow contiguous access
// num_x is the number of points that need converstion
// T is the input type, i.e., the type of the transformed points
// C is the output type, i.e., the type of the canonical points and the rate and shift arrays
// THREADS can be as high as possible (e.g., 1024), but the algorithm will fail if num_dimensions exceeds THREADS (which is highly unlikely)
// also see AccelerationDomainTransform::AccelerationDomainTransform(), where gpu_trans_a/b are encoded and size_a is computed
template<typename T, typename C, int THREADS> // transformed and canonical types
__global__ void tasgpu_transformed_to_canonical(int dims, int num_x, int size_a, const C *gpu_trans_a, const C *gpu_trans_b, const T *gpu_x_transformed, T *gpu_x_canonical){
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
        gpu_x_canonical[k] = ((T) gpu_x_transformed[k]) * rate[i] - shift[i];
        k += k_stride;
        i = k % size_a;
    }
}

// convert (-1, 1) to (0, 1)
template<typename T, int THREADS>
__global__ void tasgpu_m11_to_01(int num_points, T *gpu_x){
    int i = blockIdx.x * THREADS + threadIdx.x;
    while (i < num_points){
        T v = gpu_x[i];
        v += 1.0;
        v /= 2.0;
        gpu_x[i] = v;
        i += THREADS * gridDim.x;
    }
}

__device__ inline double linear_boundary_wavelet(double x){
    if (fabs(x + 0.5) > 0.5) return 0.0;
    if (x <= -0.75){
        return 0.75 * (7.0 * x + 6.0);
    }else if (x <= -0.5){
        return -0.25 * (11.0 * x + 6.0);
    }else{
        return 0.25 * x;
    }
}
__device__ inline double linear_central_wavelet(double x){
    if (fabs(x + 0.25) > 0.75) return 0.0;
    if (x <= -0.5){
        return -0.5 * x - 0.5;
    }else if (x >= 0.0){
        return 0.5 * x - 0.25;
    }else if (x <= -0.25){
        return 4.0 * x + 1.75;
    }else{
        return -4.0 * x - 0.25;
    }
}

// evaluates a single basis function
// syncronize with GridLocalPolynomial::encodeSupportForGPU() for the meaning of negative support
// in essense, the formula is that feval = 1 - |x - node| / support, or feval = 1 - (x - node)^2 / support^2
// there are exceptions for the constant function on level 0, the "half linear" functions on level 1 for localp, and the global functions on level 1 for semi-localp
template <typename T, int order, TypeOneDRule rule>
__device__ inline T tasgpu_devalpwpoly_feval(const double x, const double node, const double support){ // <- arrays are cached
    T v;
    if (rule == rule_localp){
        if (order == 0){
            v = (fabs(x - node) > support) ? 0.0 : 1.0;
        }else if (order == 1){
            v = 1.0 - fabs(x - node) / support;
            if (support == -1.0) v = 1.0;
            if (v < 0.0) v = 0.0;
        }else if (order == 2){
            v = x - node;
            v *= v;
            v = 1.0 - v / support;
            if (support == -1.0) v = 1.0;
            if (support == -2.0) v = -x;
            if (support == -3.0) v =  x;
            if (v < 0.0) v = 0.0;
        }
    }else if (rule == rule_localp0){
        if (order == 1){
            v = 1.0 - fabs(x - node) / support;
            if (v < 0.0) v = 0.0;
        }else if (order == 2){
            v = x - node;
            v *= v;
            v = 1.0 - v / support;
            if (v < 0.0) v = 0.0;
        }
    }else if (rule == rule_localpb){
        if (order == 1){
            v = 1.0 - fabs(x - node) / support;
            if (v < 0.0) v = 0.0;
        }else if (order == 2){
            if (support == -2.0){
                v = 1.0 + fabs(x - node) / support;
            }else{
                v = x - node;
                v *= v;
                v = 1.0 - v / support;
                if (v < 0.0) v = 0.0;
            }
        }
    }else if (rule == rule_semilocalp){
        if (order == 2){
            v = x - node;
            v *= v;
            v = 1.0 - v / support;
            if (v < 0.0) v = 0.0;
            if (support == -1.0) v = 1.0;
            if (support == -4.0) v = 0.5 * x * (x - 1.0);
            if (support == -5.0) v = 0.5 * x * (x + 1.0);
        }
    }else{ // wavelet, node = scale, support = shift and encodes the type of wavelet
        // sync with RuleWavelet::getShiftScale()
        if (support == -1.0){ // level 0, using basic hat-functions
            v = 1.0 - fabs(x - node);
            if (v < 0.0) v = 0.0;
        }else if (support == -2.0){ // left boundary
            v = linear_boundary_wavelet(node * (x + 1.0) - 1.0);
        }else if (support == -3.0){ // right boundary
            v = linear_boundary_wavelet(node * (1.0 - x) - 1.);
        }else{
            v = linear_central_wavelet(node * (x + 1.0) - 1.0 - support);
        }
    }
    //printf(" v = %1.4e   order = %d   x = %1.4e  node = %1.4e  supp = %1.4e \n", v, order, x, node, support);
    return v;
}

// evaluates the sparse grid functions for localp0, local0 and semilocalp rules using a DENSE algorithm
// use with MAX_DIM = 64 and SHORT = 32 or MAX_DIM = 32 and SHORT = 64
// dims is num_dimensions, this kernel cannot handle more than 64-dimensions, which is plenty for a sparse grid
// grid blocks take 32 points each, cache the points, then iteratively take successive blocks of 32 x values (32 x dims)
// each block of 32 gpu_x values is cached and used to form a 32 by 32 block of the values matrix
// thus, the nodes and support are reused and only the smaller array of x vlaues is read repeatedly (each thread-block reads the whole array once)
// gpu_x is canonical x (same format as CPU, dim-major format)
// use 1D thread-blocks with (num_points / 32) blocks
template <typename T, int order, TypeOneDRule rule, int SHORT, int MAX_DIM> // rule: localp0, localp, semilocalp
__global__ void tasgpu_devalpwpoly(int dims, int num_x, int num_points, const T *gpu_x, const T *gpu_nodes, const T *gpu_support, T *gpu_y){
    __shared__ T cache_nodes[SHORT * MAX_DIM];
    __shared__ T cache_support[SHORT * MAX_DIM];
    __shared__ T cache_x[SHORT * MAX_DIM];

    int height = threadIdx.x % SHORT; // threads are oranized logically in a square of size SHORT by SHORT
    int swidth = threadIdx.x / SHORT; // height and swidth are the indexes of this thread in the block (relates to num_points and num_x respectively)
    int max_p = dims * num_points; // maximum number of entries in nodes and support, avoid reading beyond
    int max_x = dims * num_x; // maximum number of entries in gpu_x, avoid reading beyond

    int id_p = SHORT * blockIdx.x;

    while(id_p < num_points){

        // cache support and nodes (reused a lot)
        int i = threadIdx.x;
        while(i < SHORT * dims){
            if (id_p * dims + i < max_p){
                cache_nodes[(i % dims) * SHORT + (i / dims)] = gpu_nodes[id_p * dims + i];
                cache_support[(i % dims) * SHORT + (i / dims)] = gpu_support[id_p * dims + i];
            }
            i += 1024;
        }
        __syncthreads();

        int id_x = 0;
        while(id_x < num_x){

            // cache gpu_x (used only for this block)
            i = threadIdx.x;
            while(i < SHORT * dims){
                if (id_x * dims + i < max_x)
                    cache_x[i] = gpu_x[id_x * dims + i];
                i += 1024;
            }
            __syncthreads();

            T v = 1.0;
            for(int j=0; j<dims; j++){
                v *= tasgpu_devalpwpoly_feval<T, order, rule>(cache_x[swidth * dims + j], cache_nodes[height + SHORT * j], cache_support[height + SHORT * j]);
            }

            if ((id_p + height < num_points) && (id_x + swidth < num_x))
                gpu_y[(id_x + swidth) * num_points + (id_p + height)] = v;
            __syncthreads();

            id_x += (1024 / SHORT);

        }
        id_p += SHORT * gridDim.x;
    }
}


// evaluates a single basis function, for use with sparse evaluations
template <typename T, int order, TypeOneDRule rule>
__device__ inline T tasgpu_devalpwpoly_basis_multid(int dims, int i, int ip, const T *x, const T *nodes, const T *support){
    T p = 1.0;
    for(int j=0; j<dims; j++){
        const T this_x = x[i * dims + j];
        const T this_node = nodes[ip * dims + j];
        const T this_supp = support[ip * dims + j];
        p *= tasgpu_devalpwpoly_feval<T, order, rule>(this_x, this_node, this_supp);
    }
    return p;
}

// the refinement process for piece wise constant rule (pwc) should cover a region larger than the regular support
// this is needed since refinement has to be performed between existing nodes and not just within the support
// see: Miroslav Stoyanov, "Adaptive sparse grid construction in a context of local anisotropy and multiple hierarchical parents"
// Sparse Grids with Applications, 2016 (pre-print)
template<typename T>
__device__ inline T tasgpu_devalpwpoly_support_pwc_multid(int dims, int i, int ip, const T *x, const T *nodes, const T *support){
    T p = 1.0;
    for(int j=0; j<dims; j++){
        p *= ((fabs(x[i * dims + j] - nodes[ip * dims + j]) > 2.0 * support[ip * dims + j]) ? 0.0 : 1.0);
    }
    return p;
}

// essentially the same as local polynomial evaluate, just repeated for many threads, this is SPARSE algorithm
// at 64 threads and 46 TOPLEVEL, the 48K of shared memory is ehxausted with 4 byte integer
// note that the resolution at level 46 is 1 / 2^46 = O(1.E-13), which is comparative to the precision, i.e., never use more than 46 levels, makes no sense
// the kernel is called twice with fill false and then true
// false only counts the non-zeros, sindx and svals are ignored
// true fills the non-zeros in sindx and svals, sindx and svals must be pre-allocated
template <typename T, int THREADS, int TOPLEVEL, int order, TypeOneDRule rule, bool fill>
__global__ void tasgpu_devalpwpoly_sparse(int dims, int num_x, const T *x, const T *nodes, const T *support,
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
            T p;
            if (order > 0){
                p = tasgpu_devalpwpoly_basis_multid<T, order, rule>(dims, i, ip, x, nodes, support);
            }else{
                p = tasgpu_devalpwpoly_support_pwc_multid<T>(dims, i, ip, x, nodes, support);
            }

            if (p != 0.0){
                if (fill){
                    sindx[c] = ip;
                    if (order == 0){
                        p = tasgpu_devalpwpoly_basis_multid<T, order, rule>(dims, i, ip, x, nodes, support);;
                    }
                    svals[c] = p;
                }
                c++;

                int current = 0;
                mstop[0][threadIdx.x] = hpntr[ip + 1];
                mcount[0][threadIdx.x] = hpntr[ip];

                while(mcount[0][threadIdx.x] < mstop[0][threadIdx.x]){
                    if (mcount[current][threadIdx.x] < mstop[current][threadIdx.x]){
                        ip = hindx[mcount[current][threadIdx.x]];
                        if (order > 0){
                            p = tasgpu_devalpwpoly_basis_multid<T, order, rule>(dims, i, ip, x, nodes, support);
                        }else{
                            p = tasgpu_devalpwpoly_support_pwc_multid<T>(dims, i, ip, x, nodes, support);
                        }

                        if (p != 0.0){
                            if (fill){
                                sindx[c] = ip;
                                if (order == 0){
                                    p = tasgpu_devalpwpoly_basis_multid<T, order, rule>(dims, i, ip, x, nodes, support);;
                                }
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

        if (!fill){
            spntr[i+1] = c;
        }

        i += gridDim.x * THREADS;
    }
}

// call with THREADS > dims, THREADS > max_num_nodes
// creates a cache data structure for evaluating basis functions with sequence grids
// dims and num_x are number of dimensions and number of x points, gpuX holds the points in the same format as the CPU version
// nodes and coeffs are the 1D sequence nodes and the normalizing hierarchical Newton coefficients (same as the GPU version)
// max_num_nodes is the maxumum number of nodes in any direction (needed to correctly capture the cache)
// the result (i.e., cache) has block form, each of the offsets give cumulative index of the next block
// each bock has size num_nodes[d] X num_x, holding the values of the Newton polynomials associated with x in direction d
// blocks do not have identical structure, due to anisotropy block sizes can vary
// num_nodes is the maximum number of nodes in the corresponding direction
// cache is contiguous in the direction of x
template <typename T, int THREADS>
__global__ void tasgpu_dseq_build_cache(int dims, int num_x, const T *gpuX, const T *nodes, const T *coeffs, int max_num_nodes, const int *offsets, const int *num_nodes, T *result){
    __shared__ int soffsets[THREADS]; // cache offsets, nodes and coefficients
    __shared__ T snodes[THREADS];
    __shared__ T scoeffs[THREADS];
    if (threadIdx.x < dims)
        soffsets[threadIdx.x] = offsets[threadIdx.x];

    if (threadIdx.x < max_num_nodes){
        snodes[threadIdx.x] = nodes[threadIdx.x];
        scoeffs[threadIdx.x] = coeffs[threadIdx.x];
    }
    __syncthreads();

    int i = blockIdx.x * THREADS + threadIdx.x;

    while(i < num_x){
        for(int d=0; d<dims; d++){
            int offset = soffsets[d] + i;
            T x = gpuX[i * dims + d]; // non-strided read

            T v = 1.0; // v holds newton polynomial corresponding to x in direction d and node j
            result[offset] = 1.0;
            offset += num_x;
            for(int j=1; j<num_nodes[d]; j++){
                v *= (x - snodes[j-1]); // everybody in the block will read the same shared node value
                result[offset] = v / scoeffs[j];
                offset += num_x;
            }
        }

        i += gridDim.x * THREADS;
    }
}

// use with SHORT = 32 or 64
// cache is computed above, dims is the number of dimensions
// num_x and num_points give the dimension of the basis matrix (data is contiguous in index of points)
// points is an array of integers, but it the transpose of the CPU points
template <typename T, int SHORT>
__global__ void tasgpu_dseq_eval_sharedpoints(int dims, int num_x, int num_points, const int *points, const int *offsets, const T *cache, T *result){
    __shared__ int cpoints[1024 / SHORT][SHORT]; // uses only 32 * 32 * 4 = 4K of 48K shared memory

    int height = threadIdx.x % SHORT; // threads are oranized logically in a square of size SHORT by SHORT
    int swidth = threadIdx.x / SHORT; // height and swidth are the indexes of this thread in the block (relates to num_points and num_x respectively)

    int id_p = SHORT * blockIdx.x;

    while(id_p < num_points){
        if ((swidth < dims) && (id_p + height < num_points)){
            cpoints[swidth][height] = points[swidth * num_points + id_p + height]; // load a block of points
        }
        __syncthreads();

        if (id_p + height < num_points){
            int id_x = swidth;
            while(id_x < num_x){ // loop over the x-values

                // read and multiply from cached data strucutre, but uses only L2/1 cache, no shared GPU memory
                double v = cache[num_x * cpoints[0][height] + id_x];
                for(int j=1; j<dims; j++){
                    v *= cache[offsets[j] + num_x * cpoints[j][height] + id_x];
                }

                result[num_points * id_x + id_p + height] = v;

                id_x += (1024 / SHORT);
            }
        }
        id_p += SHORT * gridDim.x;
        __syncthreads();
    }
}

// call with THREADS > dims
// creates a cache data structure for evaluating basis functions with fourier grids
// dims and num_x are number of dimensions and number of x points, gpuX holds the points in the same format as the CPU version
// each bock has size 2 X num_nodes[d] X num_x, holding the values of the Fourier basis associated with x in direction d
// the result (i.e., cache) has block form, each of the offsets give cumulative index of the next block
// blocks do not have identical structure, due to anisotropy block sizes can vary
// num_nodes is the maximum number of nodes in the corresponding direction (i.e., maximum power)
// cache is contiguous in the direction of x
template <typename T, int THREADS>
__global__ void tasgpu_dfor_build_cache(int dims, int num_x, const T *gpuX, const int *offsets, const int *num_nodes, T *result){
    __shared__ int soffsets[THREADS]; // cache offsets, nodes and coefficients
    if (threadIdx.x < dims)
        soffsets[threadIdx.x] = offsets[threadIdx.x];
    __syncthreads();

    int i = blockIdx.x * THREADS + threadIdx.x;

    while(i < num_x){
        for(int d=0; d<dims; d++){
            int offset = soffsets[d] + i;
            T x = gpuX[i * dims + d]; // non-strided read

            T step_real = cos(-2.0 * 3.14159265358979323846 * x);
            T step_imag = sin(-2.0 * 3.14159265358979323846 * x); // start with exp(-i x)

            T vreal = 1.0; // start with 1.0 (holds the current power)
            T vimag = 0.0;

            result[offset] = 1.0;
            offset += num_x;
            result[offset] = 0.0;
            offset += num_x;
            for(int j=1; j<num_nodes[d]; j+=2){
                T v = vreal * step_real - vimag * step_imag;
                vimag = vreal * step_imag + vimag * step_real;
                vreal = v;

                result[offset] = vreal;
                offset += num_x;
                result[offset] = vimag; // exp(-((j+1)/2) * i * x)
                offset += num_x;
                result[offset] = vreal;
                offset += num_x;
                result[offset] = -vimag; // conjugate for exp( ((j+1)/2) * i * x)
                offset += num_x;
            }
        }

        i += gridDim.x * THREADS;
    }
}

// use with SHORT = 32 or 64, Dimensions < 32
// cache is computed above, dims is the number of dimensions
// num_x and num_points give the dimension of the basis matrix (data is contiguous in index of points)
// points is an array of integers, but it the transpose of the CPU points
template <typename T, int SHORT, bool interlace>
__global__ void tasgpu_dfor_eval_sharedpoints(int dims, int num_x, int num_points, const int *points, const int *offsets, const T *cache, T *wreal, T *wimag){
    __shared__ int cpoints[1024 / SHORT][SHORT]; // uses only 32 * 32 * 4 = 4K of 48K shared memory

    int height = threadIdx.x % SHORT; // threads are oranized logically in a square of size SHORT by SHORT
    int swidth = threadIdx.x / SHORT; // height and swidth are the indexes of this thread in the block (relates to num_points and num_x respectively)

    int id_p = SHORT * blockIdx.x;

    while(id_p < num_points){
        if ((swidth < dims) && (id_p + height < num_points)){
            cpoints[swidth][height] = points[swidth * num_points + id_p + height]; // load a block of points
        }
        __syncthreads();

        if (id_p + height < num_points){
            int id_x = swidth;
            while(id_x < num_x){ // loop over the x-values

                // read and multiply from cached data strucutre, but uses only L2/1 cache, no shared GPU memory
                T vreal = cache[2 * num_x * cpoints[0][height] + id_x];
                T vimag = cache[2 * num_x * cpoints[0][height] + num_x + id_x];
                for(int j=1; j<dims; j++){
                    T sreal = cache[offsets[j] + 2 * num_x * cpoints[j][height] + id_x];
                    T simag = cache[offsets[j] + 2 * num_x * cpoints[j][height] + num_x + id_x];

                    T v = vreal * sreal - vimag * simag;
                    vimag = vreal * simag + vimag * sreal;
                    vreal = v;
                }

                if (interlace){
                    wreal[2 * (num_points * id_x + id_p + height)] = vreal;
                    wreal[2 * (num_points * id_x + id_p + height) + 1] = vimag;
                }else{
                    wreal[num_points * id_x + id_p + height] = vreal;
                    wimag[num_points * id_x + id_p + height] = vimag;
                }

                id_x += (1024 / SHORT);
            }
        }
        id_p += SHORT * gridDim.x;
        __syncthreads();
    }
}

// Builds the cache matrix for global basis evals
template<typename T, int NUM_THREADS, bool nested, bool is_cc0>
__global__ void tasgpu_dglo_build_cache(int num_dims, int num_x, int cache_lda, T const *gpu_x, T const *nodes, T const *coeff,
                                        int const *nodes_per_level, int const *offset_per_level, int const *dim_offsets,
                                        int const *map_dimension, int const *map_level, T *cache){
    int blkid = blockIdx.x;

    while(blkid < cache_lda){ // loop over all dim-level pairs
        int idx = threadIdx.x;

        int dim = map_dimension[blkid]; // get the dim for this thread block
        int lvl = map_level[blkid]; // get the level for the thread block
        int num_nodes = nodes_per_level[lvl];
        int opl = offset_per_level[lvl];

        while(idx < num_x){ // loop over all num_x
            T x = gpu_x[idx * num_dims + dim]; // get the local x

            int cache_offset = (dim_offsets[dim] + opl) * num_x;

            T c = 1.0;
            cache[cache_offset + idx] = c;
            for(int j=0; j<num_nodes-1; j++){
                c *= (nested) ? (x - nodes[j]) : (x - nodes[opl + j]);
                cache_offset += num_x;
                cache[cache_offset + idx] = c;
            }

            c = (is_cc0) ? (x * x - 1.0) : 1.0;
            cache[cache_offset + idx] *= c * coeff[opl + num_nodes - 1];
            cache_offset -= num_x;
            for(int j=num_nodes-1; j>0; j--){
                c *= (nested) ? (x - nodes[j]) : (x - nodes[opl + j]);
                cache[cache_offset + idx] *= c * coeff[opl + j - 1];
                cache_offset -= num_x;
            }

            idx += NUM_THREADS;
        }

        blkid += gridDim.x;
    }
}

template <typename T, int THREADS>
__global__ void tasgpu_dglo_eval_zero(int size, T *result){
    int i = blockIdx.x * THREADS + threadIdx.x;
    while(i < size){
        result[i] = 0.0;
        i += gridDim.x * THREADS;
    }
}

template <typename T, int NUM_THREADS>
__global__ void tasgpu_dglo_eval_sharedpoints(int num_dims, int num_x, int loop_lda, int result_lda, T const *cache,
                                              T const *tweights,
                                              int const *offset_per_level, int const *dim_offsets,
                                              int const *active_tensors, int const *active_num_points,
                                              int const *map_tensor, int const *map_index, int const *map_reference,
                                              T *result){
    int blkid = blockIdx.x;

    while(blkid < loop_lda){ // loop over all dim-level pairs
        int idx = threadIdx.x;

        int tensor = map_tensor[blkid]; // get the tensor for this thread block
        int ref    = map_reference[blkid];
        int toff   = tensor * num_dims;

        while(idx < num_x){ // loop over all num_x
            int index  = map_index[blkid]; // get the index for the thread block

            T w = 1.0;
            for(int j=num_dims-1; j>=0; j--){
                w *= cache[(dim_offsets[j] + offset_per_level[active_tensors[toff + j]] + index % active_num_points[toff + j]) * num_x + idx];
                index /= active_num_points[toff + j];
            }

            atomicAdd(&result[idx * result_lda + ref], tweights[tensor] * w);

            idx += NUM_THREADS;
        }

        blkid += gridDim.x;
    }
}

}

#endif
