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

void TasCUDA::dtrans2can(int dims, int num_x, int pad_size, const double *gpu_trans_a, const double *gpu_trans_b, const double *gpu_x_transformed, double *gpu_x_canonical){
    int num_blocks = (num_x * dims) / _MAX_CUDA_THREADS + (((num_x * dims) % _MAX_CUDA_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_transformed_to_canonical<double, double, _MAX_CUDA_THREADS><<<num_blocks, _MAX_CUDA_THREADS, (2*pad_size) * sizeof(double)>>>(dims, num_x, pad_size, gpu_trans_a, gpu_trans_b, gpu_x_transformed, gpu_x_canonical);
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
    }else{ // rule == rule_semilocalp
        tasgpu_devalpwpoly<double, 2, rule_semilocalp, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
    }
}

// there is a switch statement that realizes templates for each combination of rule/order
// make one function that covers that switch, the rest is passed from devalpwpoly_sparse and devalpwpoly_sparse_dense
template<typename T, int THREADS, int TOPLEVEL, bool fill, bool dense>
inline void devalpwpoly_sparse_realize_rule_order(int order, TypeOneDRule rule,
                                          int dims, int num_x, int num_points,
                                          const T *x, const T *nodes, const T *support,
                                          const int *hpntr, const int *hindx, int num_roots, const int *roots,
                                          int *spntr, int *sindx, T *svals, T *gpu_dense){
    int num_blocks = num_x / THREADS + ((num_x % THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    if (rule == rule_localp){
        switch(order){
            case 0:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 0, rule_localp, fill, dense><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals, gpu_dense);
                break;
            case 2:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_localp, fill, dense><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals, gpu_dense);
                break;
            default:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 1, rule_localp, fill, dense><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals, gpu_dense);
        }
    }else if (rule == rule_localp0){
        switch(order){
            case 2:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_localp0, fill, dense><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals, gpu_dense);
                break;
            default:
                tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 1, rule_localp0, fill, dense><<<num_blocks, THREADS>>>
                    (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals, gpu_dense);
        }
    }else{ // rule == rule_semilocalp
        tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_semilocalp, fill, dense><<<num_blocks, THREADS>>>
            (dims, num_x, num_points, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals, gpu_dense);
    }
}

// local polynomial basis functions, SPARSE algorithm (2 passes, one pass to compue the non-zeros and one pass to evaluate)
void TasCUDA::devalpwpoly_sparse(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                                 int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nz,
                                 std::ostream *os){
    gpu_spntr = cudaNew<int>(num_x + 1, os);
    // call with fill == false to count the non-zeros per row of the matrix
    devalpwpoly_sparse_realize_rule_order<double, 64, 46, false, false>
        (order, rule, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, 0, 0, 0);

    int *cpu_spntr = new int[num_x+1];
    cudaRecv(num_x+1, gpu_spntr, cpu_spntr, os);
    cpu_spntr[0] = 0;
    for(int i=1; i<=num_x; i++) cpu_spntr[i] += cpu_spntr[i-1];
    num_nz = cpu_spntr[num_x]; // save the number of non-zeros
    cudaSend<int>(num_x + 1, cpu_spntr, gpu_spntr, os);
    delete[] cpu_spntr;
    gpu_sindx = cudaNew<int>(num_nz, os);
    gpu_svals = cudaNew<double>(num_nz, os);
    // call with fill == true to load the non-zeros
    devalpwpoly_sparse_realize_rule_order<double, 64, 46, true, false>
        (order, rule, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, gpu_spntr, gpu_sindx, gpu_svals, 0);
}

void TasCUDA::devalpwpoly_sparse_dense(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                                 int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, double *gpu_dense){
    // call with fill = false, dense = true, to load the values in the dense arrays
    devalpwpoly_sparse_realize_rule_order<double, 64, 46, false, true>
        (order, rule, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_hpntr, gpu_hindx, num_roots, gpu_roots, 0, 0, 0, gpu_dense);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       Linear Algebra
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TasCUDA::cudaDgemm(int M, int N, int K, const double *gpu_a, const double *gpu_b, double *gpu_c){ // gpu_c = gpu_a * gpu_b, gpu_c is M by N
    int blocks = (N / 96) + (((N % 96) == 0) ? 0 : 1);
    blocks *= (M / 96) + (((M % 96) == 0) ? 0 : 1);
    while(blocks > 65536) blocks = 65536;
    tasgpu_cudaTgemm<double, 32, 96><<<blocks, 1024>>>(M, N, K, gpu_a, gpu_b, gpu_c);
}

// sparse triangular solve, usnig cuda kernels (as opposed to cusparse)
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

    int num_blocks = N / _MAX_CUDA_THREADS + ((N % _MAX_CUDA_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_d3gecss<<< num_blocks, _MAX_CUDA_THREADS >>>(N, M, gpu_order, top_level, gpuBpntr, gpuBindx, gpuBvals, gpuX, num_blocks * _MAX_CUDA_THREADS);

    TasCUDA::cudaRecv<double>(M*N, gpuX, cpuC, os);

    delete[] order;
    TasCUDA::cudaDel<double>(gpuX, os);
    TasCUDA::cudaDel<double>(gpuBvals, os);
    TasCUDA::cudaDel<int>(gpuBindx, os);
    TasCUDA::cudaDel<int>(gpuBpntr, os);
    TasCUDA::cudaDel<int>(gpu_order, os);
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

}

#endif
