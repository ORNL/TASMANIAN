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

#include "tsgAcceleratedDataStructures.hpp"

#include "tsgCudaMacros.hpp"

#include "tsgCudaLinearAlgebra.hpp"

#include "tsgCudaBasisEvaluations.hpp"

// several kernels assume linear distribution the threads and can be executed with "practically unlimited" number of threads
// thus we can set this to the CUDA max number of threads, based on the computer version
#define _MAX_CUDA_THREADS 1024

namespace TasGrid{

void TasCUDA::dtrans2can(bool use01, int dims, int num_x, int pad_size, const double *gpu_trans_a, const double *gpu_trans_b, const double *gpu_x_transformed, double *gpu_x_canonical){
    int num_blocks = (num_x * dims) / _MAX_CUDA_THREADS + (((num_x * dims) % _MAX_CUDA_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_transformed_to_canonical<double, double, _MAX_CUDA_THREADS><<<num_blocks, _MAX_CUDA_THREADS, (2*pad_size) * sizeof(double)>>>(dims, num_x, pad_size, gpu_trans_a, gpu_trans_b, gpu_x_transformed, gpu_x_canonical);
    if (use01) tasgpu_m11_to_01<double, _MAX_CUDA_THREADS><<<num_blocks, _MAX_CUDA_THREADS>>>(dims * num_x, gpu_x_canonical);
}

// local polynomial basis functions, DENSE algorithm
void TasCUDA::devalpwpoly(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support, double *gpu_y){
    // each block thread runs 1024 threads and processes 32 points (or basis functions)
    int num_blocks = (num_points / 32) + ((num_points % 32 == 0) ? 0 : 1);
    // order == 1 is considered "default" so that the compiler doesn't complain about missing default statement
    // semilocalp cannot have order less than 2, only rule_localp can have order 0 (this gets overwrittein in makeLocalPolynomialGrid())
    if (rule == rule_localp){
        switch(order){
            case 0:
                    tasgpu_devalpwpoly<double, 0, rule_localp, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            case 2: tasgpu_devalpwpoly<double, 2, rule_localp, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    tasgpu_devalpwpoly<double, 1, rule_localp, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else if (rule == rule_localp0){
        switch(order){
            case 2: tasgpu_devalpwpoly<double, 2, rule_localp0, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    tasgpu_devalpwpoly<double, 1, rule_localp0, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else if (rule == rule_localpb){
        switch(order){
            case 2: tasgpu_devalpwpoly<double, 2, rule_localpb, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    tasgpu_devalpwpoly<double, 1, rule_localpb, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else{ // rule == rule_semilocalp
        tasgpu_devalpwpoly<double, 2, rule_semilocalp, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
    }
}

// there is a switch statement that realizes templates for each combination of rule/order
// make one function that covers that switch, the rest is passed from devalpwpoly_sparse
template<typename T, int THREADS, int TOPLEVEL, bool fill>
inline void devalpwpoly_sparse_realize_rule_order(int order, TypeOneDRule rule,
                                          int dims, int num_x, int num_points,
                                          const T *x, const T *nodes, const T *support,
                                          const int *hpntr, const int *hindx, int num_roots, const int *roots,
                                          int *spntr, int *sindx, T *svals){
    int num_blocks = num_x / THREADS + ((num_x % THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    if (rule == rule_localp){
        switch(order){
            case 0:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 0, rule_localp, fill><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
                break;
            case 2:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_localp, fill><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
                break;
            default:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 1, rule_localp, fill><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
        }
    }else if (rule == rule_localp0){
        switch(order){
            case 2:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_localp0, fill><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
                break;
            default:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 1, rule_localp0, fill><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
        }
    }else if (rule == rule_localpb){
        switch(order){
            case 2:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_localpb, fill><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
                break;
            default:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 1, rule_localpb, fill><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
        }
    }else{ // rule == rule_semilocalp
        tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_semilocalp, fill><<<num_blocks, THREADS>>>
            (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
    }
}

// local polynomial basis functions, SPARSE algorithm (2 passes, one pass to compue the non-zeros and one pass to evaluate)
void TasCUDA::devalpwpoly_sparse(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const cudaDoubles &gpu_nodes, const cudaDoubles &gpu_support,
                            const cudaInts &gpu_hpntr, const cudaInts &gpu_hindx, const  cudaInts &gpu_roots, cudaInts &gpu_spntr, cudaInts &gpu_sindx, cudaDoubles &gpu_svals){
    gpu_spntr.resize(num_x + 1);
    // call with fill == false to count the non-zeros per row of the matrix
    devalpwpoly_sparse_realize_rule_order<double, 64, 46, false>
        (order, rule, dims, num_x, num_points, gpu_x, gpu_nodes.data(), gpu_support.data(), gpu_hpntr.data(), gpu_hindx.data(), (int) gpu_roots.size(), gpu_roots.data(), gpu_spntr.data(), 0, 0);

    std::vector<int> cpu_spntr;
    gpu_spntr.unload(cpu_spntr);
    cpu_spntr[0] = 0;
    int nz = 0;
    for(auto &i : cpu_spntr){
        i += nz;
        nz = i;
    }
    gpu_spntr.load(cpu_spntr);
    gpu_sindx.resize(nz);
    gpu_svals.resize(nz);
    // call with fill == true to load the non-zeros
    devalpwpoly_sparse_realize_rule_order<double, 64, 46, true>
        (order, rule, dims, num_x, num_points, gpu_x, gpu_nodes.data(), gpu_support.data(), gpu_hpntr.data(), gpu_hindx.data(), (int) gpu_roots.size(), gpu_roots.data(), gpu_spntr.data(), gpu_sindx.data(), gpu_svals.data());
}

// Sequence Grid basis evaluations
void TasCUDA::devalseq(int dims, int num_x, const std::vector<int> &max_levels, const double *gpu_x, const cudaInts &num_nodes, const cudaInts &points, const cudaDoubles &nodes, const cudaDoubles &coeffs, double *gpu_result){
    std::vector<int> offsets(dims);
    offsets[0] = 0;
    for(int d=1; d<dims; d++) offsets[d] = offsets[d-1] + num_x * (max_levels[d-1] + 1);
    size_t num_total = offsets[dims-1] + num_x * (max_levels[dims-1] + 1);

    int maxl = max_levels[0]; for(auto l : max_levels) if (maxl < l) maxl = l;

    cudaInts gpu_offsets(offsets);
    cudaDoubles cache1D(num_total);
    int num_blocks = num_x / _MAX_CUDA_THREADS + ((num_x % _MAX_CUDA_THREADS == 0) ? 0 : 1);

    tasgpu_dseq_build_cache<double, _MAX_CUDA_THREADS><<<num_blocks, _MAX_CUDA_THREADS>>>
        (dims, num_x, gpu_x, nodes.data(), coeffs.data(), maxl+1, gpu_offsets.data(), num_nodes.data(), cache1D.data());

    num_blocks = num_x / 32 + ((num_x % 32 == 0) ? 0 : 1);
    tasgpu_dseq_eval_sharedpoints<double, 32><<<num_blocks, 1024>>>
        (dims, num_x, (int) points.size() / dims, points.data(), gpu_offsets.data(), cache1D.data(), gpu_result);
}

// Fourier Grid basis evaluations
void TasCUDA::devalfor(int dims, int num_x, const std::vector<int> &max_levels, const double *gpu_x, const cudaInts &num_nodes, const cudaInts &points, double *gpu_wreal, double *gpu_wimag){
    std::vector<int> max_nodes(dims);
    for(int j=0; j<dims; j++){
        int n = 1;
        for(int i=0; i<max_levels[j]; i++) n *= 3;
        max_nodes[j] = n;
    }

    std::vector<int> offsets(dims);
    offsets[0] = 0;
    for(int d=1; d<dims; d++) offsets[d] = offsets[d-1] + 2 * num_x * (max_nodes[d-1] + 1);
    size_t num_total = offsets[dims-1] + 2 * num_x * (max_nodes[dims-1] + 1);

    cudaInts gpu_offsets(offsets);
    cudaDoubles cache1D(num_total);
    int num_blocks = num_x / _MAX_CUDA_THREADS + ((num_x % _MAX_CUDA_THREADS == 0) ? 0 : 1);

    tasgpu_dfor_build_cache<double, _MAX_CUDA_THREADS><<<num_blocks, _MAX_CUDA_THREADS>>>
        (dims, num_x, gpu_x, gpu_offsets.data(), num_nodes.data(), cache1D.data());

    num_blocks = num_x / 32 + ((num_x % 32 == 0) ? 0 : 1);
    if (gpu_wimag == 0){
        tasgpu_dfor_eval_sharedpoints<double, 32, true><<<num_blocks, 1024>>>
            (dims, num_x, (int) points.size() / dims, points.data(), gpu_offsets.data(), cache1D.data(), gpu_wreal, 0);
    }else{
        tasgpu_dfor_eval_sharedpoints<double, 32, false><<<num_blocks, 1024>>>
            (dims, num_x, (int) points.size() / dims, points.data(), gpu_offsets.data(), cache1D.data(), gpu_wreal, gpu_wimag);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       Linear Algebra
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __TASMANIAN_COMPILE_FALLBACK_CUDA_KERNELS__
void TasCUDA::cudaDgemm(int M, int N, int K, const double *gpu_a, const double *gpu_b, double *gpu_c){ // gpu_c = gpu_a * gpu_b, gpu_c is M by N
    int blocks = (N / 96) + (((N % 96) == 0) ? 0 : 1);
    blocks *= (M / 96) + (((M % 96) == 0) ? 0 : 1);
    while(blocks > 65536) blocks = 65536;
    tasgpu_cudaTgemm<double, 32, 96><<<blocks, 1024>>>(M, N, K, gpu_a, gpu_b, gpu_c);
}

void TasCUDA::cudaSparseMatmul(int M, int N, int num_nz, const int* gpu_spntr, const int* gpu_sindx, const double* gpu_svals, const double *gpu_B, double *gpu_C){
    int blocks = M / 64 + ((M % 64 == 0) ? 0 : 1);
    tasgpu_sparse_matmul<double, 64><<<blocks, 64>>>(M, N, num_nz, gpu_spntr, gpu_sindx, gpu_svals, gpu_B, gpu_C);
}

void TasCUDA::cudaSparseVecDenseMat(int M, int N, int num_nz, const double *A, const int *indx, const double *vals, double *C){
    int num_blocks = N / _MAX_CUDA_THREADS + ((N % _MAX_CUDA_THREADS == 0) ? 0 : 1);
    if (num_blocks< 65536){
        tasgpu_sparse_matveci<double, _MAX_CUDA_THREADS, 1><<<num_blocks, _MAX_CUDA_THREADS>>>(M, N, num_nz, A, indx, vals, C);
    }else{
        num_blocks = N / (2 * _MAX_CUDA_THREADS) + ((N % (2 * _MAX_CUDA_THREADS) == 0) ? 0 : 1);
        if (num_blocks< 65536){
            tasgpu_sparse_matveci<double, _MAX_CUDA_THREADS, 2><<<num_blocks, _MAX_CUDA_THREADS>>>(M, N, num_nz, A, indx, vals, C);
        }else{
            num_blocks = N / (3 * _MAX_CUDA_THREADS) + ((N % (3 * _MAX_CUDA_THREADS) == 0) ? 0 : 1);
            if (num_blocks< 65536){
                tasgpu_sparse_matveci<double, _MAX_CUDA_THREADS, 3><<<num_blocks, _MAX_CUDA_THREADS>>>(M, N, num_nz, A, indx, vals, C);
            }
        }
    }
}

void TasCUDA::convert_sparse_to_dense(int num_rows, int num_columns, const int *pntr, const int *indx, const double *vals, double *destination){
    int n = num_rows * num_columns;
    int num_blocks = n / _MAX_CUDA_THREADS + ((n % _MAX_CUDA_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tascuda_fill<double, _MAX_CUDA_THREADS><<<num_blocks, _MAX_CUDA_THREADS>>>(n, 0.0, destination);
    num_blocks = num_rows;
    if (num_blocks >= 65536) num_blocks = 65536;
    tascuda_sparse_to_dense<double, 64><<<num_blocks, 64>>>(num_rows, num_columns, pntr, indx, vals, destination);
}
#endif

}

#endif
