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

namespace TasGrid{

void TasCUDA::dtrans2can(int dims, int num_x, int pad_size, const double *gpu_trans_a, const double *gpu_trans_b, const double *gpu_x_transformed, double *gpu_x_canonical){
    int num_blocks = (num_x * dims) / TASMANIAN_CUDA_NUM_THREADS + (((num_x * dims) % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_transformed_to_canonical<double, double, TASMANIAN_CUDA_NUM_THREADS><<<num_blocks, TASMANIAN_CUDA_NUM_THREADS, (2*pad_size) * sizeof(double)>>>(dims, num_x, pad_size, gpu_trans_a, gpu_trans_b, gpu_x_transformed, gpu_x_canonical);
}

// local polynomial basis functions, DENSE algorithm
void TasCUDA::devalpwpoly(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support, double *gpu_y){
    // each block thread runs 1024 threads and processes 32 points (or basis functions)
    int num_blocks = (num_points / 32) + ((num_points % 32 == 0) ? 0 : 1);
    // order == 1 is considered "default" so that the compiler doesn't complain about missing default statement
    // semilocalp cannot have order less than 2, only rule_localp can have order 0 (this gets overwrittein in makeLocalPolynomialGrid())
    if (rule == rule_localp0){
        switch(order){
            case 2: tasgpu_devalpwpoly<double, 2, rule_localp0, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    tasgpu_devalpwpoly<double, 1, rule_localp0, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else if (rule == rule_localp){
        switch(order){
            case 0:
                    tasgpu_devalpwpoly<double, 0, rule_localp, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            case 2: tasgpu_devalpwpoly<double, 2, rule_localp, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    tasgpu_devalpwpoly<double, 1, rule_localp, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else{ // rule == rule_semilocalp
        tasgpu_devalpwpoly<double, 2, rule_semilocalp, 32, 64><<<num_blocks, 1024>>>(dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
    }
}

// local polynomial basis functions, SPARSE algorithm (2 passes, one pass to compue the non-zeros and one pass to evaluate)
void TasCUDA::devalpwpoly_sparse(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                                 int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nz,
                                 std::ostream *os){
    gpu_spntr = cudaNew<int>(num_x + 1, os);
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
    cudaRecv(num_x+1, gpu_spntr, cpu_spntr, os);
    cpu_spntr[0] = 0;
    for(int i=1; i<=num_x; i++) cpu_spntr[i] += cpu_spntr[i-1];
    num_nz = cpu_spntr[num_x]; // save the number of non-zeros
    cudaSend<int>(num_x + 1, cpu_spntr, gpu_spntr, os);
    delete[] cpu_spntr;
    gpu_sindx = cudaNew<int>(num_nz, os);
    gpu_svals = cudaNew<double>(num_nz, os);
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

void TasCUDA::devalseq(int dims, int num_x, int num_points, int num_nodes, const double *gpu_x, const double *gpu_nodes, const double *gpu_coeff, const int *points, double *gpu_dense){

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

    int num_blocks = N / TASMANIAN_CUDA_NUM_THREADS + ((N % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tasgpu_d3gecss<<< num_blocks, TASMANIAN_CUDA_NUM_THREADS >>>(N, M, gpu_order, top_level, gpuBpntr, gpuBindx, gpuBvals, gpuX, num_blocks * TASMANIAN_CUDA_NUM_THREADS);

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
    int num_blocks = n / TASMANIAN_CUDA_NUM_THREADS + ((n % TASMANIAN_CUDA_NUM_THREADS == 0) ? 0 : 1);
    if (num_blocks >= 65536) num_blocks = 65536;
    tascuda_fill<double, TASMANIAN_CUDA_NUM_THREADS><<<num_blocks, TASMANIAN_CUDA_NUM_THREADS>>>(n, 0.0, destination);
    num_blocks = num_rows;
    if (num_blocks >= 65536) num_blocks = 65536;
    tascuda_sparse_to_dense<double, 64><<<num_blocks, 64>>>(num_rows, num_columns, pntr, indx, vals, destination);
}

void TasCUDA::cudaSparseMatmul(int M, int N, int num_nz, const int* gpu_spntr, const int* gpu_sindx, const double* gpu_svals, const double *gpu_B, double *gpu_C){
    int blocks = M / 64 + ((M % 64 == 0) ? 0 : 1);
    tasgpu_sparse_matmul<double, 64><<<blocks, 64>>>(M, N, num_nz, gpu_spntr, gpu_sindx, gpu_svals, gpu_B, gpu_C);
}

}

#endif
