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
    }else{
        if (order == 2){
            v = x - node;
            v *= v;
            v = 1.0 - v / support;
            if (v < 0.0) v = 0.0;
            if (support == -1.0) v = 1.0;
            if (support == -4.0) v = 0.5 * x * (x - 1.0);
            if (support == -5.0) v = 0.5 * x * (x + 1.0);
        }
    }
    //printf(" v = %1.4e   order = %d   x = %1.4e  node = %1.4e  supp = %1.4e \n", v, order, x, node, support);
    return v;
}

// evaluates the sparse grid functions for localp0, local0 and semilocalp rules using a DENSE algorithm
// use with MAX_DIM = 64 and SHORT = 32
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
            i += SHORT * SHORT;
        }
        __syncthreads();

        int id_x = 0;
        while(id_x < num_x){

            // cache gpu_x (used only for this block)
            i = threadIdx.x;
            while(i < SHORT * dims){
                if (id_x * dims + i < max_x)
                    cache_x[i] = gpu_x[id_x * dims + i];
                i += SHORT * SHORT;
            }
            __syncthreads();

            T v = 1.0;
            for(int j=0; j<dims; j++){
                v *= tasgpu_devalpwpoly_feval<T, order, rule>(cache_x[swidth * dims + j], cache_nodes[height + SHORT * j], cache_support[height + SHORT * j]);
            }

            if ((id_p + height < num_points) && (id_x + swidth < num_x))
                gpu_y[(id_x + swidth) * num_points + (id_p + height)] = v;
            __syncthreads();

            id_x += SHORT;

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
// dense being true means that the algorithm will fill the dense matrix "dense" and ignore the sparse (spntr, sindx, svals) (call dense == true only when fill == false)
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
                if (dense) gpu_dense[i * num_points + ip] = p;

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
