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

#ifndef __TASMANIAN_SPARSE_GRID_CUDA_MACROS_HPP
#define __TASMANIAN_SPARSE_GRID_CUDA_MACROS_HPP

#include "TasmanianConfig.hpp"

#include <iostream>

#if defined(Tasmanian_ENABLE_CUBLAS) || defined(Tasmanian_ENABLE_CUDA)
#include <cuda_runtime_api.h>
#include <cuda.h>
#endif

#ifdef Tasmanian_ENABLE_CUBLAS
#include <cublas_v2.h>
#include <cusparse.h>
#endif // Tasmanian_ENABLE_CUBLAS

namespace TasGrid{

#if defined(Tasmanian_ENABLE_CUBLAS) || defined(Tasmanian_ENABLE_CUDA)
namespace TasCUDA{
    // general GPU/CPU I/O, new/delete vectors, send/recv data
    template <typename T>
    inline T* cudaNew(size_t num_entries, std::ostream *os){
        T* x = 0;
        cudaError_t cudaStat = cudaMalloc(((void**) &x), num_entries * sizeof(T));
        AccelerationMeta::cudaCheckError((void*) &cudaStat, "cudaNew()", os);
        return x;
    }

    template <typename T>
    inline T* cudaSend(size_t num_entries, const T *cpu_array, std::ostream *os){
        T *x = cudaNew<T>(num_entries, os);
        cudaError_t cudaStat = cudaMemcpy(x, cpu_array, num_entries * sizeof(T), cudaMemcpyHostToDevice);
        AccelerationMeta::cudaCheckError((void*) &cudaStat, "cudaSend(type)", os);
        return x;
    }

    template <typename T>
    inline void cudaSend(size_t num_entries, const T *cpu_array, T *gpu_array, std::ostream *os){
        cudaError_t cudaStat = cudaMemcpy(gpu_array, cpu_array, num_entries * sizeof(T), cudaMemcpyHostToDevice);
        AccelerationMeta::cudaCheckError((void*) &cudaStat, "cudaSend(type, type)", os);
    }

    template <typename T>
    inline void cudaRecv(size_t num_entries, const T *gpu_array, T *cpu_array, std::ostream *os){
        cudaError_t cudaStat = cudaMemcpy(cpu_array, gpu_array, num_entries * sizeof(T), cudaMemcpyDeviceToHost);
        AccelerationMeta::cudaCheckError((void*) &cudaStat, "cudaRecv(type, type)", os);
    }

    template <typename T>
    inline void cudaDel(T *gpu_array, std::ostream *os){
        cudaError_t cudaStat = cudaFree(gpu_array);
        AccelerationMeta::cudaCheckError((void*) &cudaStat, "cudaDel(type)", os);
    }

}
#endif

}

#endif
