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

#ifndef __TASMANIAN_SPARSE_GRID_DPCPP_BASIS_EVALUATIONS_HPP
#define __TASMANIAN_SPARSE_GRID_DPCPP_BASIS_EVALUATIONS_HPP

#include "TasmanianConfig.hpp"

#ifndef Tasmanian_ENABLE_DPCPP
#error "Cannot use tsgDpcppBasisEvaluations.hpp without Tasmanian_ENABLE_DPCPP"
#endif

#define MKL_INT int
#include <CL/sycl.hpp>
#include <CL/sycl/usm.hpp>
#include "oneapi/mkl.hpp"

namespace TasGrid{

template<typename T, int NUM_THREADS, bool nested, bool is_cc0>
void tasgpu_dglo_build_cache(sycl::queue *q,
                             int num_dims, int num_x, int cache_lda, T const *gpu_x, T const *nodes, T const *coeff,
                             int const *nodes_per_level, int const *offset_per_level, int const *dim_offsets,
                             int const *map_dimension, int const *map_level, T *cache){

//     q->submit([&](sycl::handler& h) {
//             h.parallel_for<class tsg_transpose>(sycl::range<2>{static_cast<size_t>(m), static_cast<size_t>(n)}, [=](sycl::id<2> i){
//                 AT[i[0] * n + i[1]] = A[i[0] + m * i[1]];
//             });
//         });
//     q->wait();

    q->submit([&](sycl::handler& h) {
        h.parallel_for<class tsg_glo_build_cache_kernel>(sycl::range<2>{static_cast<size_t>(cache_lda), static_cast<size_t>(NUM_THREADS)}, [=](sycl::id<2> i){

            int blkid = i[0];

            int idx = i[1];

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
        });
    });
    q->wait();

//     int blkid = blockIdx.x;
//
//     while(blkid < cache_lda){ // loop over all dim-level pairs
//         int idx = threadIdx.x;
//
//         int dim = map_dimension[blkid]; // get the dim for this thread block
//         int lvl = map_level[blkid]; // get the level for the thread block
//         int num_nodes = nodes_per_level[lvl];
//         int opl = offset_per_level[lvl];
//
//         while(idx < num_x){ // loop over all num_x
//             T x = gpu_x[idx * num_dims + dim]; // get the local x
//
//             int cache_offset = (dim_offsets[dim] + opl) * num_x;
//
//             T c = 1.0;
//             cache[cache_offset + idx] = c;
//             for(int j=0; j<num_nodes-1; j++){
//                 c *= (nested) ? (x - nodes[j]) : (x - nodes[opl + j]);
//                 cache_offset += num_x;
//                 cache[cache_offset + idx] = c;
//             }
//
//             c = (is_cc0) ? (x * x - 1.0) : 1.0;
//             cache[cache_offset + idx] *= c * coeff[opl + num_nodes - 1];
//             cache_offset -= num_x;
//             for(int j=num_nodes-1; j>0; j--){
//                 c *= (nested) ? (x - nodes[j]) : (x - nodes[opl + j]);
//                 cache[cache_offset + idx] *= c * coeff[opl + j - 1];
//                 cache_offset -= num_x;
//             }
//
//             idx += NUM_THREADS;
//         }
//
//         blkid += gridDim.x;
//     }
}

template <typename T>
void tasgpu_dglo_eval_zero(sycl::queue *q, int size, T *result){
    q->submit([&](sycl::handler& h) {
        h.parallel_for<class tsg_glo_build_cache_kernel>(sycl::range<1>{static_cast<size_t>(size),}, [=](sycl::id<1> i){
            result[i[0]] = 0.0;
        });
    });
    q->wait();
}

template <typename T, int NUM_THREADS>
void tasgpu_dglo_eval_sharedpoints(sycl::queue *q,
                                   int num_dims, int num_x, int loop_lda, int result_lda, T const *cache,
                                   T const *tweights,
                                   int const *offset_per_level, int const *dim_offsets,
                                   int const *active_tensors, int const *active_num_points,
                                   int const *map_tensor, int const *map_index, int const *map_reference, T *result){
    for(int blkid = 0; blkid < loop_lda; blkid++){
    //while(blkid < loop_lda){ // loop over all dim-level pairs
        //int idx = threadIdx.x;

        int tensor = map_tensor[blkid]; // get the tensor for this thread block
        int ref    = map_reference[blkid];
        int toff   = tensor * num_dims;

        //while(idx < num_x){ // loop over all num_x
        for(int idx = 0; idx < num_x; idx++){
            int index  = map_index[blkid]; // get the index for the thread block

            T w = 1.0;
            for(int j=num_dims-1; j>=0; j--){
                w *= cache[(dim_offsets[j] + offset_per_level[active_tensors[toff + j]] + index % active_num_points[toff + j]) * num_x + idx];
                index /= active_num_points[toff + j];
            }

            //atomicAdd(&result[idx * result_lda + ref], tweights[tensor] * w);
            result[idx * result_lda + ref] += tweights[tensor] * w;

            //idx += NUM_THREADS;
        }

        //blkid += gridDim.x;
    }
}

}

#endif