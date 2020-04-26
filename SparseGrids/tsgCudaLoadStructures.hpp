/*
 * Copyright (c) 2018, Miroslav Stoyanov
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

#ifndef __TASMANIAN_SPARSE_CUDA_LOAD_STRUCTS_HPP
#define __TASMANIAN_SPARSE_CUDA_LOAD_STRUCTS_HPP

#include "tsgAcceleratedDataStructures.hpp"

namespace TasGrid{


#ifndef __TASMANIAN_DOXYGEN_SKIP

/*!
 * \internal
 * \brief Wrapper structure for the vectors needed by Global grid CUDA methods.
 *
 * \endinternal
 */
template<typename FP>
struct CudaGlobalData{
    GpuVector<FP> values;
    // cache stage
    int num_basis; // number of 1d basis functions
    GpuVector<FP> nodes;
    GpuVector<FP> coeff;
    GpuVector<int> nodes_per_level;
    GpuVector<int> offset_per_level; // non-nested case nodes, always used for coefficients
    GpuVector<int> map_dimension;
    GpuVector<int> map_level;
    // compute stage
    GpuVector<FP> tensor_weights;
    GpuVector<int> active_tensors;
    GpuVector<int> active_num_points;
    GpuVector<int> dim_offsets; // relates to the cache
    GpuVector<int> map_tensor;
    GpuVector<int> map_index;
    GpuVector<int> map_reference;
    void clearNodes(){
        num_basis = 0;
        nodes.clear();
        coeff.clear();
        nodes_per_level.clear();
        offset_per_level.clear();
        map_dimension.clear();
        map_level.clear();
        tensor_weights.clear();
        active_tensors.clear();
        active_num_points.clear();
        dim_offsets.clear();
        map_tensor.clear();
        map_index.clear();
        map_reference.clear();
    }
};

//! \internal
//! \brief Wrapper structure for the vectors needed by Sequence grid CUDA methods.

//! The \b surplusses is a copy of the Sequence grid surpluses and is used in the linear algebra (stage 2) of the evaluations.
//! The \b nodes and \b coeffs are the 1D nodes and coefficients, needed for the basis evaluations (stage 1),
//! more specifically used in the caching stage of the evaluation call.
//! The \b points is a transposed copy of the nodes in the MultiIndexSet and \b num_nodes is one more than the number of nodes in each dimension.
//! The \b points and \b num_nodes are needed for the accumulate stage of the basis evaluations.
template<typename FP>
struct CudaSequenceData{
    GpuVector<FP> surpluses, nodes, coeff;
    GpuVector<int> num_nodes, points;
    void clearNodes(){
        nodes.clear();
        coeff.clear();
        num_nodes.clear();
        points.clear();
    }
};

//! \internal
//! \brief Wrapper structure for the vectors needed by Fourier grid CUDA methods.

//! The \b real and \b imag values correspond to the real and complex part of the Fourier coefficients.
//! The \b points is a transposed copy of the nodes in the MultiIndexSet and \b num_nodes is the number of nodes in each dimension.
template<typename FP>
struct CudaFourierData{
    GpuVector<FP> real, imag;
    GpuVector<int> num_nodes, points;
};

//! \internal
//! \brief Wrapper structure for the vectors needed by Local Polynomial grid CUDA methods.

//! The \b surpluses are the hierarchical surpluses and used in the linear algebra (stage 2) of the evaluations.
//! The \b nodes and \b support describe the basis functions, negative support is used to distinguish
//! between local and global support and different order of the first few basis functions.
//! The \b hpntr, \b hindx, and \b hroots, describe the hierarchy and are used in the sparse matrix algorithm.
template<typename FP>
struct CudaLocalPolynomialData{
    GpuVector<FP> surpluses, nodes, support;
    GpuVector<int> hpntr, hindx, hroots;
    void clearHierarchy(){
        hpntr.clear();
        hindx.clear();
        hroots.clear();
    }
    void clearBasisHierarchy(){
        nodes.clear();
        support.clear();
        clearHierarchy();
    }
};

//! \internal
//! \brief Wrapper structure for the vectors needed by Wavelet grid CUDA methods.

//! The \b coefficients are the hierarchical coefficients and used in the linear algebra (stage 2) of the evaluations.
template<typename FP>
struct CudaWaveletData{
    GpuVector<FP> coefficients, nodes, support;
    void clearNodes(){
        nodes.clear();
        support.clear();
    }
};
#endif // __TASMANIAN_DOXYGEN_SKIP

}

#endif
