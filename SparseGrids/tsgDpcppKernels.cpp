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

namespace TasGrid{

template<typename T>
void TasGpu::dtrans2can(bool use01, int dims, int num_x, int pad_size, double const *gpu_trans_a, double const *gpu_trans_b, T const *gpu_x_transformed, T *gpu_x_canonical){

}

template void TasGpu::dtrans2can<double>(bool, int, int, int, double const*, double const*, double const*, double*);
template void TasGpu::dtrans2can<float>(bool, int, int, int, double const*, double const*, float const*, float*);

// local polynomial basis functions, DENSE algorithm
template<typename T>
void TasGpu::devalpwpoly(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const T *gpu_x, const T *gpu_nodes, const T *gpu_support, T *gpu_y){

}

template void TasGpu::devalpwpoly<double>(int, TypeOneDRule, int, int, int, const double*, const double*, const double*, double*);
template void TasGpu::devalpwpoly<float>(int, TypeOneDRule, int, int, int, const float*, const float*, const float*, float*);

// there is a switch statement that realizes templates for each combination of rule/order
// make one function that covers that switch, the rest is passed from devalpwpoly_sparse
template<typename T, int THREADS, int TOPLEVEL, bool fill>
inline void devalpwpoly_sparse_realize_rule_order(int order, TypeOneDRule rule, int dims, int num_x,
                                          const T *x, const T *nodes, const T *support,
                                          const int *hpntr, const int *hindx, int num_roots, const int *roots,
                                          int *spntr, int *sindx, T *svals){

}

// local polynomial basis functions, SPARSE algorithm (2 passes, one pass to compue the non-zeros and one pass to evaluate)
template<typename T>
void TasGpu::devalpwpoly_sparse(int order, TypeOneDRule rule, int dims, int num_x, const T *gpu_x,
                                const GpuVector<T> &gpu_nodes, const GpuVector<T> &gpu_support,
                                const GpuVector<int> &gpu_hpntr, const GpuVector<int> &gpu_hindx, const GpuVector<int> &gpu_hroots,
                                GpuVector<int> &gpu_spntr, GpuVector<int> &gpu_sindx, GpuVector<T> &gpu_svals){

}
template void TasGpu::devalpwpoly_sparse<double>(int, TypeOneDRule, int, int, const double*, const GpuVector<double>&, const GpuVector<double>&,
                                                 const GpuVector<int>&, const GpuVector<int>&, const GpuVector<int>&,
                                                 GpuVector<int>&, GpuVector<int>&, GpuVector<double>&);
template void TasGpu::devalpwpoly_sparse<float>(int, TypeOneDRule, int, int, const float*, const GpuVector<float>&, const GpuVector<float>&,
                                                const GpuVector<int>&, const GpuVector<int>&, const GpuVector<int>&,
                                                GpuVector<int>&, GpuVector<int>&, GpuVector<float>&);

// Sequence Grid basis evaluations
template<typename T>
void TasGpu::devalseq(int dims, int num_x, const std::vector<int> &max_levels, const T *gpu_x, const GpuVector<int> &num_nodes,
                      const GpuVector<int> &points, const GpuVector<T> &nodes, const GpuVector<T> &coeffs, T *gpu_result){

}

template void TasGpu::devalseq<double>(int dims, int num_x, const std::vector<int> &max_levels, const double *gpu_x, const GpuVector<int> &num_nodes,
                                      const GpuVector<int> &points, const GpuVector<double> &nodes, const GpuVector<double> &coeffs, double *gpu_result);
template void TasGpu::devalseq<float>(int dims, int num_x, const std::vector<int> &max_levels, const float *gpu_x, const GpuVector<int> &num_nodes,
                                      const GpuVector<int> &points, const GpuVector<float> &nodes, const GpuVector<float> &coeffs, float *gpu_result);

// Fourier Grid basis evaluations
template<typename T>
void TasGpu::devalfor(int dims, int num_x, const std::vector<int> &max_levels, const T *gpu_x,
                      const GpuVector<int> &num_nodes, const GpuVector<int> &points, T *gpu_wreal, typename GpuVector<T>::value_type *gpu_wimag){

}

template void TasGpu::devalfor<double>(int, int, const std::vector<int>&, const double*, const GpuVector<int>&, const GpuVector<int>&, double*, double*);
template void TasGpu::devalfor<float>(int, int, const std::vector<int>&, const float*, const GpuVector<int>&, const GpuVector<int>&, float*, float*);

template<typename T>
void TasGpu::devalglo(bool is_nested, bool is_clenshawcurtis0, int dims, int num_x, int num_p, int num_basis,
                      T const *gpu_x, GpuVector<T> const &nodes, GpuVector<T> const &coeff, GpuVector<T> const &tensor_weights,
                      GpuVector<int> const &nodes_per_level, GpuVector<int> const &offset_per_level, GpuVector<int> const &map_dimension, GpuVector<int> const &map_level,
                      GpuVector<int> const &active_tensors, GpuVector<int> const &active_num_points, GpuVector<int> const &dim_offsets,
                      GpuVector<int> const &map_tensor, GpuVector<int> const &map_index, GpuVector<int> const &map_reference, T *gpu_result){

}

template void TasGpu::devalglo<double>(bool, bool, int, int, int, int,
                                       double const*, GpuVector<double> const&, GpuVector<double> const&, GpuVector<double> const&,
                                       GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                       GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                       GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, double*);
template void TasGpu::devalglo<float>(bool, bool, int, int, int, int,
                                      float const*, GpuVector<float> const&, GpuVector<float> const&, GpuVector<float> const&,
                                      GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                      GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                      GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, float*);

void TasGpu::fillDataGPU(double value, long long n, long long stride, double data[]){

}

}

#endif
