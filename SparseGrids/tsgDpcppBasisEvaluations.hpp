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

template<typename T, typename C> struct tsg_transformed_to_canonical_kernel{};

template<typename T, typename C> // transformed and canonical types
void tasgpu_transformed_to_canonical(sycl::queue *q, int dims, int num_x, int size_a, const C *gpu_trans_a, const C *gpu_trans_b, const T *gpu_x_transformed, T *gpu_x_canonical){

    q->submit([&](sycl::handler& h) {
        h.parallel_for<tsg_transformed_to_canonical_kernel<T, C>>(sycl::range<1>{static_cast<size_t>(dims * num_x), }, [=](sycl::id<1> threadId){
	    int i = threadId[0] % size_a;
	    gpu_x_canonical[threadId[0]] = static_cast<T>(gpu_x_transformed[threadId[0]] * gpu_trans_a[i] - gpu_trans_b[i]);
        });
    });
    q->wait();

}

template<typename T> struct tsg_m11_to_01_kernel{};

// convert (-1, 1) to (0, 1)
template<typename T>
void tasgpu_m11_to_01(sycl::queue *q, int num_points, T *gpu_x){

    q->submit([&](sycl::handler& h) {
        h.parallel_for<tsg_m11_to_01_kernel<T>>(sycl::range<1>{static_cast<size_t>(num_points), }, [=](sycl::id<1> threadId){
            size_t i = threadId[0];
            gpu_x[i] = ( gpu_x[i] + 1.0 ) / 2.0;
        });
    });
    q->wait();

}

template <typename T> struct tsg_seq_build_cache_kernel{};

template <typename T>
void tasgpu_dseq_build_cache(sycl::queue *q, int dims, int num_x, const T *gpuX, const T *nodes, const T *coeffs, const int *offsets, const int *num_nodes, T *result){

    q->submit([&](sycl::handler& h) {
        h.parallel_for<tsg_seq_build_cache_kernel<T>>(sycl::range<1>{static_cast<size_t>(num_x), }, [=](sycl::id<1> threadId){
            int i = threadId[0];
            for(int d=0; d<dims; d++){
                int offset = offsets[d] + i;
                T x = gpuX[i * dims + d];

                T v = 1.0; // v holds newton polynomial corresponding to x in direction d and node j
                result[offset] = 1.0;
                offset += num_x;
                for(int j=1; j<num_nodes[d]; j++){
                    v *= (x - nodes[j-1]); // everybody in the block will read the same shared node value
                    result[offset] = v / coeffs[j];
                    offset += num_x;
                }
            }
        });
    });
    q->wait();
}

template <typename T> struct tsg_seq_eval_sharedpoints_kernel{};

template <typename T>
void tasgpu_dseq_eval_sharedpoints(sycl::queue *q, int dims, int num_x, int num_points, const int *points, const int *offsets, const T *cache, T *result){

    q->submit([&](sycl::handler& h) {
        h.parallel_for<tsg_seq_eval_sharedpoints_kernel<T>>(sycl::range<2>{static_cast<size_t>(num_points), static_cast<size_t>(num_x)}, [=](sycl::id<2> threadId){
            int id_p = threadId[0];
            int id_x = threadId[1];

            double v = cache[num_x * points[id_p] + id_x];
            for(int j=1; j<dims; j++){
                v *= cache[offsets[j] + num_x * points[j * num_points + id_p] + id_x];
            }

            result[num_points * id_x + id_p] = v;
        });
    });
    q->wait();
}
template<typename T>
inline T linear_boundary_wavelet(T x){
    if (fabs(x + 0.5) > 0.5) return 0.0;
    if (x <= -0.75){
        return 0.75 * (7.0 * x + 6.0);
    }else if (x <= -0.5){
        return -0.25 * (11.0 * x + 6.0);
    }else{
        return 0.25 * x;
    }
}
template<typename T>
inline T linear_central_wavelet(T x){
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
template <typename T, int order, TypeOneDRule rule>
inline T tasgpu_devalpwpoly_feval(const double x, const double node, const double support){ // <- arrays are cached
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
    return v;
}
template<typename T, int order, int rule> struct tasgpu_devalpwpoly_kernel{};

template <typename T, int order, TypeOneDRule rule> // rule: localp0, localp, semilocalp
void tasgpu_devalpwpoly(sycl::queue *q, int dims, int num_x, int num_points, const T *gpu_x, const T *gpu_nodes, const T *gpu_support, T *gpu_y){

    q->submit([&](sycl::handler& h) {
        h.parallel_for<tasgpu_devalpwpoly_kernel<T, order, rule>>(sycl::range<2>{static_cast<size_t>(num_points), static_cast<size_t>(num_x)}, [=](sycl::id<2> threadId){

            int id_p = threadId[0];
            int id_x = threadId[1];

            T v = 1.0;
            for(int j=0; j<dims; j++)
                v *= tasgpu_devalpwpoly_feval<T, order, rule>(gpu_x[id_x * dims + j], gpu_nodes[id_p * dims + j], gpu_support[id_p * dims + j]);

            gpu_y[id_x * num_points + id_p] = v;
        });
    });
    q->wait();
}

template <typename T, int order, TypeOneDRule rule>
inline T tasgpu_devalpwpoly_basis_multid(int dims, int i, int ip, const T *x, const T *nodes, const T *support){
    T p = 1.0;
    for(int j=0; j<dims; j++){
        const T this_x = x[i * dims + j];
        const T this_node = nodes[ip * dims + j];
        const T this_supp = support[ip * dims + j];
        p *= tasgpu_devalpwpoly_feval<T, order, rule>(this_x, this_node, this_supp);
    }
    return p;
}
template<typename T>
inline T tasgpu_devalpwpoly_support_pwc_multid(int dims, int i, int ip, const T *x, const T *nodes, const T *support){
    T p = 1.0;
    for(int j=0; j<dims; j++){
        p *= ((fabs(x[i * dims + j] - nodes[ip * dims + j]) > 2.0 * support[ip * dims + j]) ? 0.0 : 1.0);
    }
    return p;
}

template <typename T, int TOPLEVEL, int order, int rule, bool fill> struct tasgpu_devalpwpoly_sparse_kernel{};

template <typename T, int TOPLEVEL, int order, TypeOneDRule rule, bool fill>
void tasgpu_devalpwpoly_sparse(sycl::queue *q, int dims, int num_x, const T *x, const T *nodes, const T *support,
                               const int *hpntr, const int *hindx, int num_roots, const int *roots,
                               int *spntr, int *sindx, T *svals){

    q->submit([&](sycl::handler& h) {
        h.parallel_for<tasgpu_devalpwpoly_sparse_kernel<T, TOPLEVEL, order, rule, fill>>(sycl::range<1>{static_cast<size_t>(num_x),}, [=](sycl::id<1> threadId){
            int mcount[TOPLEVEL];
            int mstop[TOPLEVEL];

            int c = 0;
            int i = threadId[0];

            if (fill) c = spntr[i];

            for(int r=0; r<num_roots; r++){
                int ip = roots[r];
                T p = (order > 0) ?
                    tasgpu_devalpwpoly_basis_multid<T, order, rule>(dims, i, ip, x, nodes, support) :
                    tasgpu_devalpwpoly_support_pwc_multid<T>(dims, i, ip, x, nodes, support);

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
                    mstop[0] = hpntr[ip + 1];
                    mcount[0] = hpntr[ip];

                    while(mcount[0] < mstop[0]){
                        if (mcount[current] < mstop[current]){
                            ip = hindx[mcount[current]];
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
                                mstop[current] = hpntr[ip + 1];
                                mcount[current] = hpntr[ip];
                            }else{
                                mcount[current]++;
                            }
                        }else{
                            current--;
                            mcount[current]++;
                        }
                    }
                }
            }

            if (not fill) spntr[i+1] = c;
        });
    });
    q->wait();
}

template <typename T> struct tasgpu_dfor_build_cache_kernel{};

template <typename T>
void tasgpu_dfor_build_cache(sycl::queue *q, int dims, int num_x, const T *gpuX, const int *offsets, const int *num_nodes, T *result){

    q->submit([&](sycl::handler& h) {
        h.parallel_for<tasgpu_dfor_build_cache_kernel<T>>(sycl::range<1>{static_cast<size_t>(num_x),}, [=](sycl::id<1> threadId){
            int i = threadId[0];
            for(int d=0; d<dims; d++){
                int offset = offsets[d] + i;
                T x = gpuX[i * dims + d]; // non-strided read

                T step_real = sycl::cos(-2.0 * 3.14159265358979323846 * x);
                T step_imag = sycl::sin(-2.0 * 3.14159265358979323846 * x); // start with exp(-i x)

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
        });
    });
    q->wait();
}

template <typename T, bool interlace> struct tasgpu_dfor_eval_sharedpoints_kernel{};

template <typename T, bool interlace>
void tasgpu_dfor_eval_sharedpoints(sycl::queue *q, int dims, int num_x, int num_points, const int *points, const int *offsets, const T *cache, T *wreal, T *wimag){

    q->submit([&](sycl::handler& h) {
        h.parallel_for<tasgpu_dfor_eval_sharedpoints_kernel<T, interlace>>(sycl::range<2>{static_cast<size_t>(num_points), static_cast<size_t>(num_x)}, [=](sycl::id<2> threadId){
            int id_p = threadId[0];
            int id_x = threadId[1];

            T vreal = cache[2 * num_x * points[id_p] + id_x];
            T vimag = cache[2 * num_x * points[id_p] + num_x + id_x];
            for(int j=1; j<dims; j++){
                T sreal = cache[offsets[j] + 2 * num_x * points[j * num_points + id_p] + id_x];
                T simag = cache[offsets[j] + 2 * num_x * points[j * num_points + id_p] + num_x + id_x];

                T v = vreal * sreal - vimag * simag;
                vimag = vreal * simag + vimag * sreal;
                vreal = v;
            }

            if (interlace){
                wreal[2 * (num_points * id_x + id_p)] = vreal;
                wreal[2 * (num_points * id_x + id_p) + 1] = vimag;
            }else{
                wreal[num_points * id_x + id_p] = vreal;
                wimag[num_points * id_x + id_p] = vimag;
            }
        });
    });
    q->wait();
}

template<typename T, bool nested, bool is_cc0> struct tsg_glo_build_cache_kernel{};

template<typename T, bool nested, bool is_cc0>
void tasgpu_dglo_build_cache(sycl::queue *q,
                             int num_dims, int num_x, int cache_lda, T const *gpu_x, T const *nodes, T const *coeff,
                             int const *nodes_per_level, int const *offset_per_level, int const *dim_offsets,
                             int const *map_dimension, int const *map_level, T *cache){

    q->submit([&](sycl::handler& h) {
        h.parallel_for<tsg_glo_build_cache_kernel<T, nested, is_cc0>>(sycl::range<2>{static_cast<size_t>(cache_lda), static_cast<size_t>(num_x)}, [=](sycl::id<2> threadID){
            int blkid = threadID[0];
            int idx   = threadID[1];

            int dim = map_dimension[blkid]; // get the dim for this thread block
            int lvl = map_level[blkid]; // get the level for the thread block
            int num_nodes = nodes_per_level[lvl];
            int opl = offset_per_level[lvl];

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
        });
    });
    q->wait();
}

template <typename T> struct tsg_glo_eval_zero_kernel{};

template <typename T>
void tasgpu_dglo_eval_zero(sycl::queue *q, size_t size, T *result){
    q->submit([&](sycl::handler& h) {
        h.parallel_for<tsg_glo_eval_zero_kernel<T>>(sycl::range<1>{size,}, [=](sycl::id<1> i){
            result[i[0]] = static_cast<T>(0.0);
        });
    });
    q->wait();
}

template <typename T> struct tsg_glo_eval_sharedpoints_kernel{};

template <typename T>
void tasgpu_dglo_eval_sharedpoints(sycl::queue *q,
                                   int num_dims, int num_x, int loop_lda, int result_lda, T const *cache,
                                   T const *tweights,
                                   int const *offset_per_level, int const *dim_offsets,
                                   int const *active_tensors, int const *active_num_points,
                                   int const *map_tensor, int const *map_index, int const *map_reference, T *result){

    q->submit([&](sycl::handler& h){
        h.parallel_for<tsg_glo_eval_sharedpoints_kernel<T>>(sycl::range<1>{static_cast<size_t>(1),}, [=](sycl::id<1>){
            // this can be called in parallel_for over both bklid and idx (using loop_lda and num_x)
            // running this in parallel requires the atomicAdd()
            for(int blkid = 0; blkid < loop_lda; blkid++){
                int tensor = map_tensor[blkid]; // get the tensor for this thread block
                int ref    = map_reference[blkid];
                int toff   = tensor * num_dims;

                for(int idx = 0; idx < num_x; idx++){
                    int index  = map_index[blkid]; // get the index for the thread block

                    T w = 1.0;
                    for(int j=num_dims-1; j>=0; j--){
                        w *= cache[(dim_offsets[j] + offset_per_level[active_tensors[toff + j]] + index % active_num_points[toff + j]) * num_x + idx];
                        index /= active_num_points[toff + j];
                    }

                    //sycl::atomic_fetch_add<T, sycl::access::address_space::global_space>(result[idx * result_lda + ref], tweights[tensor] * w);
                    result[idx * result_lda + ref] += tweights[tensor] * w;
                }
            }
        });
    });
    q->wait();
}

}

#endif
