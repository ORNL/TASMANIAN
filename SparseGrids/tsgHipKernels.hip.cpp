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

#ifndef __TASMANIAN_SPARSE_GRID_HIP_KERNELS_HIP
#define __TASMANIAN_SPARSE_GRID_HIP_KERNELS_HIP

#include "tsgAcceleratedDataStructures.hpp"
#include "tsgCudaBasisEvaluations.hpp"
#include "tsgCudaLinearAlgebra.hpp"

// several kernels assume a linear distribution of the threads and can be executed with "practically unlimited" number of threads
// thus we can set this to the HIP max number of threads, based on the current cuda version
constexpr int _MAX_HIP_THREADS  = 1024;

// max number of blocks per grid direction
constexpr int _MAX_HIP_BLOCKS = 65535;

/*
 * Create a 1-D HIP thread grid using the total_threads and number of threads per block.
 * Basically, computes the number of blocks but no more than _MAX_CUDA_BLOCKS.
 */
struct ThreadGrid1d{
    // Compute the threads and blocks.
    ThreadGrid1d(long long total_threads, long long num_per_block) :
        threads(static_cast<int>(num_per_block)),
        blocks(static_cast<int>(std::min(total_threads / threads + ((total_threads % threads == 0) ? 0 : 1), static_cast<long long>(_MAX_HIP_BLOCKS))))
    {}
    // number of threads
    int const threads;
    // number of blocks
    int const blocks;
};

namespace TasGrid{

template<typename T>
void TasGpu::dtrans2can(AccelerationContext const*, bool use01, int dims, int num_x, int pad_size, double const *gpu_trans_a, double const *gpu_trans_b, T const *gpu_x_transformed, T *gpu_x_canonical){
    ThreadGrid1d task(num_x * dims, _MAX_HIP_THREADS);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_transformed_to_canonical<T, double, _MAX_HIP_THREADS>),
                       task.blocks, task.threads, (2*pad_size) * sizeof(double), 0,
                       dims, num_x, pad_size, gpu_trans_a, gpu_trans_b, gpu_x_transformed, gpu_x_canonical);
    if (use01) hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_m11_to_01<T, _MAX_HIP_THREADS>),
                                  task.blocks, task.threads, 0, 0, dims * num_x, gpu_x_canonical);
}

template void TasGpu::dtrans2can<double>(AccelerationContext const*, bool, int, int, int, double const*, double const*, double const*, double*);
template void TasGpu::dtrans2can<float>(AccelerationContext const*, bool, int, int, int, double const*, double const*, float const*, float*);

// local polynomial basis functions, DENSE algorithm
template<typename T>
void TasGpu::devalpwpoly(AccelerationContext const*, int order, TypeOneDRule rule, int dims, int num_x, int num_points, const T *gpu_x, const T *gpu_nodes, const T *gpu_support, T *gpu_y){
    // each block thread runs 1024 threads and processes 32 points (or basis functions)
    int num_blocks = (num_points / 32) + ((num_points % 32 == 0) ? 0 : 1);
    // order == 1 is considered "default" so that the compiler doesn't complain about missing default statement
    // semilocalp cannot have order less than 2, only rule_localp can have order 0 (this gets overwrittein in makeLocalPolynomialGrid())
    if (rule == rule_localp){
        switch(order){
            case 0:
                    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly<T, 0, rule_localp, 32, 64>), num_blocks, 1024, 0, 0, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            case 2: hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly<T, 2, rule_localp, 32, 64>), num_blocks, 1024, 0, 0, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly<T, 1, rule_localp, 32, 64>), num_blocks, 1024, 0, 0, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else if (rule == rule_localp0){
        switch(order){
            case 2: hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly<T, 2, rule_localp0, 32, 64>), num_blocks, 1024, 0, 0, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly<T, 1, rule_localp0, 32, 64>), num_blocks, 1024, 0, 0, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else if (rule == rule_localpb){
        switch(order){
            case 2: hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly<T, 2, rule_localpb, 32, 64>), num_blocks, 1024, 0, 0, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
                    break;
            default:
                    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly<T, 1, rule_localpb, 32, 64>), num_blocks, 1024, 0, 0, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
        }
    }else if (rule == rule_semilocalp){
        hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly<T, 2, rule_semilocalp, 32, 64>), num_blocks, 1024, 0, 0, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
    }else{ // rule == wavelet
        hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly<T, 1, rule_wavelet, 32, 64>), num_blocks, 1024, 0, 0, dims, num_x, num_points, gpu_x, gpu_nodes, gpu_support, gpu_y);
    }
}

template void TasGpu::devalpwpoly<double>(AccelerationContext const*, int, TypeOneDRule, int, int, int, const double*, const double*, const double*, double*);
template void TasGpu::devalpwpoly<float>(AccelerationContext const*, int, TypeOneDRule, int, int, int, const float*, const float*, const float*, float*);

// there is a switch statement that realizes templates for each combination of rule/order
// make one function that covers that switch, the rest is passed from devalpwpoly_sparse
template<typename T, int THREADS, int TOPLEVEL, bool fill>
inline void devalpwpoly_sparse_realize_rule_order(int order, TypeOneDRule rule, int dims, int num_x,
                                          const T *x, const T *nodes, const T *support,
                                          const int *hpntr, const int *hindx, int num_roots, const int *roots,
                                          int *spntr, int *sindx, T *svals){
    int num_blocks = num_x / THREADS + ((num_x % THREADS == 0) ? 0 : 1);
    if (num_blocks >= _MAX_HIP_BLOCKS) num_blocks = _MAX_HIP_BLOCKS;
    if (rule == rule_localp){
        switch(order){
            case 0:
                hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 0, rule_localp, fill>), num_blocks, THREADS, 0, 0,
                    dims, num_x, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
                break;
            case 2:
                hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_localp, fill>), num_blocks, THREADS, 0, 0,
                    dims, num_x, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
                break;
            default:
                hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 1, rule_localp, fill>), num_blocks, THREADS, 0, 0,
                    dims, num_x, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
        }
    }else if (rule == rule_localp0){
        switch(order){
            case 2:
                hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_localp0, fill>), num_blocks, THREADS, 0, 0,
                    dims, num_x, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
                break;
            default:
                hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 1, rule_localp0, fill>), num_blocks, THREADS, 0, 0,
                    dims, num_x, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
        }
    }else if (rule == rule_localpb){
        switch(order){
            case 2:
                hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_localpb, fill>), num_blocks, THREADS, 0, 0,
                    dims, num_x, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
                break;
            default:
                hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 1, rule_localpb, fill>), num_blocks, THREADS, 0, 0,
                    dims, num_x, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
        }
    }else{ // rule == rule_semilocalp
        hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_devalpwpoly_sparse<T, THREADS, TOPLEVEL, 2, rule_semilocalp, fill>), num_blocks, THREADS, 0, 0,
            dims, num_x, x, nodes, support, hpntr, hindx, num_roots, roots, spntr, sindx, svals);
    }
}

// local polynomial basis functions, SPARSE algorithm (2 passes, one pass to compue the non-zeros and one pass to evaluate)
template<typename T>
void TasGpu::devalpwpoly_sparse(AccelerationContext const*, int order, TypeOneDRule rule, int dims, int num_x, const T *gpu_x,
                                const GpuVector<T> &gpu_nodes, const GpuVector<T> &gpu_support,
                                const GpuVector<int> &gpu_hpntr, const GpuVector<int> &gpu_hindx, const GpuVector<int> &gpu_hroots,
                                GpuVector<int> &gpu_spntr, GpuVector<int> &gpu_sindx, GpuVector<T> &gpu_svals){
    gpu_spntr.resize(nullptr, num_x + 1);
    // call with fill == false to count the non-zeros per row of the matrix
    devalpwpoly_sparse_realize_rule_order<T, 64, 46, false>
        (order, rule, dims, num_x, gpu_x, gpu_nodes.data(), gpu_support.data(),
        gpu_hpntr.data(), gpu_hindx.data(), (int) gpu_hroots.size(), gpu_hroots.data(), gpu_spntr.data(), 0, 0);

    std::vector<int> cpu_spntr;
    gpu_spntr.unload(nullptr, cpu_spntr);
    cpu_spntr[0] = 0;
    int nz = 0;
    for(auto &i : cpu_spntr){
        i += nz;
        nz = i;
    }
    gpu_spntr.load(nullptr, cpu_spntr);
    gpu_sindx.resize(nullptr, nz);
    gpu_svals.resize(nullptr, nz);
    // call with fill == true to load the non-zeros
    devalpwpoly_sparse_realize_rule_order<T, 64, 46, true>
        (order, rule, dims, num_x, gpu_x, gpu_nodes.data(), gpu_support.data(),
        gpu_hpntr.data(), gpu_hindx.data(), (int) gpu_hroots.size(), gpu_hroots.data(), gpu_spntr.data(), gpu_sindx.data(), gpu_svals.data());
}
template void TasGpu::devalpwpoly_sparse<double>(AccelerationContext const*, int, TypeOneDRule, int, int,
                                                 const double*, const GpuVector<double>&, const GpuVector<double>&,
                                                 const GpuVector<int>&, const GpuVector<int>&, const GpuVector<int>&,
                                                 GpuVector<int>&, GpuVector<int>&, GpuVector<double>&);
template void TasGpu::devalpwpoly_sparse<float>(AccelerationContext const*, int, TypeOneDRule, int, int,
                                                const float*, const GpuVector<float>&, const GpuVector<float>&,
                                                const GpuVector<int>&, const GpuVector<int>&, const GpuVector<int>&,
                                                GpuVector<int>&, GpuVector<int>&, GpuVector<float>&);

// Sequence Grid basis evaluations
template<typename T>
void TasGpu::devalseq(AccelerationContext const*, int dims, int num_x,
                      const std::vector<int> &max_levels, const T *gpu_x, const GpuVector<int> &num_nodes,
                      const GpuVector<int> &points, const GpuVector<T> &nodes, const GpuVector<T> &coeffs, T *gpu_result){
    std::vector<int> offsets(dims);
    offsets[0] = 0;
    for(int d=1; d<dims; d++) offsets[d] = offsets[d-1] + num_x * (max_levels[d-1] + 1);
    size_t num_total = offsets[dims-1] + num_x * (max_levels[dims-1] + 1);

    int maxl = max_levels[0]; for(auto l : max_levels) if (maxl < l) maxl = l;

    GpuVector<int> gpu_offsets(nullptr, offsets);
    GpuVector<T> cache1D(nullptr, num_total);

    ThreadGrid1d task(num_x, _MAX_HIP_THREADS);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dseq_build_cache<T, _MAX_HIP_THREADS>), task.blocks, task.threads, 0, 0,
            dims, num_x, gpu_x, nodes.data(), coeffs.data(), maxl+1, gpu_offsets.data(), num_nodes.data(), cache1D.data());

    ThreadGrid1d task_compute(num_x, 32);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dseq_eval_sharedpoints<T, 32>), task_compute.blocks, 1024, 0, 0,
            dims, num_x, (int) points.size() / dims, points.data(), gpu_offsets.data(), cache1D.data(), gpu_result);
}

template void TasGpu::devalseq<double>(AccelerationContext const*, int dims, int num_x,
                                       const std::vector<int> &max_levels, const double *gpu_x, const GpuVector<int> &num_nodes,
                                       const GpuVector<int> &points, const GpuVector<double> &nodes, const GpuVector<double> &coeffs, double *gpu_result);
template void TasGpu::devalseq<float>(AccelerationContext const*, int dims, int num_x,
                                      const std::vector<int> &max_levels, const float *gpu_x, const GpuVector<int> &num_nodes,
                                      const GpuVector<int> &points, const GpuVector<float> &nodes, const GpuVector<float> &coeffs, float *gpu_result);

// Fourier Grid basis evaluations
template<typename T>
void TasGpu::devalfor(AccelerationContext const*, int dims, int num_x, const std::vector<int> &max_levels, const T *gpu_x,
                      const GpuVector<int> &num_nodes, const GpuVector<int> &points, T *gpu_wreal, typename GpuVector<T>::value_type *gpu_wimag){
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

    GpuVector<int> gpu_offsets(nullptr, offsets);
    GpuVector<T> cache1D(nullptr, num_total);
    ThreadGrid1d task(num_x, _MAX_HIP_THREADS);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dfor_build_cache<T, _MAX_HIP_THREADS>), task.blocks, task.threads, 0, 0,
            dims, num_x, gpu_x, gpu_offsets.data(), num_nodes.data(), cache1D.data());

    ThreadGrid1d task_compute(num_x, 32);
    if (gpu_wimag == nullptr){
        hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dfor_eval_sharedpoints<T, 32, true>), task_compute.blocks, 1024, 0, 0,
            dims, num_x, (int) points.size() / dims, points.data(), gpu_offsets.data(), cache1D.data(), gpu_wreal, nullptr);
    }else{
        hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dfor_eval_sharedpoints<T, 32, false>), task_compute.blocks, 1024, 0, 0,
            dims, num_x, (int) points.size() / dims, points.data(), gpu_offsets.data(), cache1D.data(), gpu_wreal, gpu_wimag);
    }
}

template void TasGpu::devalfor<double>(AccelerationContext const*, int, int, const std::vector<int>&, const double*, const GpuVector<int>&, const GpuVector<int>&, double*, double*);
template void TasGpu::devalfor<float>(AccelerationContext const*, int, int, const std::vector<int>&, const float*, const GpuVector<int>&, const GpuVector<int>&, float*, float*);

template<typename T>
void TasGpu::devalglo(AccelerationContext const*, bool is_nested, bool is_clenshawcurtis0, int dims, int num_x, int num_p, int num_basis,
                      T const *gpu_x, GpuVector<T> const &nodes, GpuVector<T> const &coeff, GpuVector<T> const &tensor_weights,
                      GpuVector<int> const &nodes_per_level, GpuVector<int> const &offset_per_level, GpuVector<int> const &map_dimension, GpuVector<int> const &map_level,
                      GpuVector<int> const &active_tensors, GpuVector<int> const &active_num_points, GpuVector<int> const &dim_offsets,
                      GpuVector<int> const &map_tensor, GpuVector<int> const &map_index, GpuVector<int> const &map_reference, T *gpu_result){
    GpuVector<T> cache(nullptr, num_x, num_basis);
    int num_blocks = (int) map_dimension.size();
    if (num_blocks >= _MAX_HIP_BLOCKS) num_blocks = _MAX_HIP_BLOCKS;

    if (is_nested){
        if (is_clenshawcurtis0){
            hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dglo_build_cache<T, _MAX_HIP_THREADS, true, true>), num_blocks, _MAX_HIP_THREADS, 0, 0,
                 dims, num_x, (int) map_dimension.size(), gpu_x, nodes.data(), coeff.data(),
                                        nodes_per_level.data(), offset_per_level.data(), dim_offsets.data(),
                                        map_dimension.data(), map_level.data(), cache.data());
        }else{
            hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dglo_build_cache<T, _MAX_HIP_THREADS, true, false>), num_blocks, _MAX_HIP_THREADS, 0, 0,
                 dims, num_x, (int) map_dimension.size(), gpu_x, nodes.data(), coeff.data(),
                                        nodes_per_level.data(), offset_per_level.data(), dim_offsets.data(),
                                        map_dimension.data(), map_level.data(), cache.data());
        }
    }else{
        hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dglo_build_cache<T, _MAX_HIP_THREADS, false, false>), num_blocks, _MAX_HIP_THREADS, 0, 0,
             dims, num_x, (int) map_dimension.size(), gpu_x, nodes.data(), coeff.data(),
                                    nodes_per_level.data(), offset_per_level.data(), dim_offsets.data(),
                                    map_dimension.data(), map_level.data(), cache.data());
    }

    int mat_size = num_x * num_p;
    num_blocks = num_x / _MAX_HIP_THREADS + ((mat_size % _MAX_HIP_THREADS == 0) ? 0 : 1);
    if (num_blocks >= _MAX_HIP_THREADS) num_blocks = _MAX_HIP_THREADS;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dglo_eval_zero<T, _MAX_HIP_THREADS>), num_blocks, _MAX_HIP_THREADS, 0, 0, mat_size, gpu_result);

    num_blocks = (int) map_tensor.size();
    if (num_blocks >= _MAX_HIP_BLOCKS) num_blocks = _MAX_HIP_BLOCKS;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(tasgpu_dglo_eval_sharedpoints<T, _MAX_HIP_THREADS>), num_blocks, _MAX_HIP_THREADS, 0, 0,
        dims, num_x, (int) map_tensor.size(), num_p, cache.data(),
        tensor_weights.data(), offset_per_level.data(), dim_offsets.data(), active_tensors.data(), active_num_points.data(),
        map_tensor.data(), map_index.data(), map_reference.data(), gpu_result);
}

template void TasGpu::devalglo<double>(AccelerationContext const*, bool, bool, int, int, int, int,
                                       double const*, GpuVector<double> const&, GpuVector<double> const&, GpuVector<double> const&,
                                       GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                       GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                       GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, double*);
template void TasGpu::devalglo<float>(AccelerationContext const*, bool, bool, int, int, int, int,
                                      float const*, GpuVector<float> const&, GpuVector<float> const&, GpuVector<float> const&,
                                      GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                      GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                      GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, float*);

void TasGpu::fillDataGPU(AccelerationContext const*, double value, long long n, long long stride, double data[]){
    if (stride == 1){
        ThreadGrid1d tgrid(n, _MAX_HIP_THREADS);
        hipLaunchKernelGGL(HIP_KERNEL_NAME(tascuda_vfill<double, _MAX_HIP_THREADS>), tgrid.blocks, tgrid.threads, 0, 0, n, data, value);
    }else{
        ThreadGrid1d tgrid(n, 32);
        hipLaunchKernelGGL(HIP_KERNEL_NAME(tascuda_sfill<double, 32>), tgrid.blocks, tgrid.threads, 0, 0, n, stride, data, value);
    }
}

}

#endif
