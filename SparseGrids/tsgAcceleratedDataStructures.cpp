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

#ifndef __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_CPP
#define __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_CPP

#include "tsgAcceleratedDataStructures.hpp"

namespace TasGrid{

std::map<std::string, TypeAcceleration> AccelerationMeta::getStringToAccelerationMap(){
    return {
        {"none",        accel_none},
        {"cpu-blas",    accel_cpu_blas},
        {"gpu-default", accel_gpu_default},
        {"gpu-cublas",  accel_gpu_cublas},
        {"gpu-cuda",    accel_gpu_cuda},
        {"gpu-rocblas", accel_gpu_rocblas},
        {"gpu-hip",     accel_gpu_hip},
        {"gpu-magma",   accel_gpu_magma}};
}

TypeAcceleration AccelerationMeta::getIOAccelerationString(const char * name){
    try{
        return getStringToAccelerationMap().at(name);
    }catch(std::out_of_range &){
        return accel_none;
    }
}
const char* AccelerationMeta::getIOAccelerationString(TypeAcceleration accel){
    switch (accel){
        case accel_cpu_blas:       return "cpu-blas";
        case accel_gpu_default:    return "gpu-default";
        case accel_gpu_cublas:     return "gpu-cublas";
        case accel_gpu_cuda:       return "gpu-cuda";
        case accel_gpu_magma:      return "gpu-magma";
        default: return "none";
    }
}
int AccelerationMeta::getIOAccelerationInt(TypeAcceleration accel){
    switch (accel){
        case accel_cpu_blas:       return 1;
        case accel_gpu_default:    return 3;
        case accel_gpu_cublas:     return 4;
        case accel_gpu_cuda:       return 5;
        case accel_gpu_magma:      return 6;
        default: return 0;
    }
}
TypeAcceleration AccelerationMeta::getIOIntAcceleration(int accel){
    switch (accel){
        case 1:  return accel_cpu_blas;
        case 3:  return accel_gpu_default;
        case 4:  return accel_gpu_cublas;
        case 5:  return accel_gpu_cuda;
        case 6:  return accel_gpu_magma;
        default: return accel_none;
    }
}
bool AccelerationMeta::isAccTypeGPU(TypeAcceleration accel){
    switch (accel){
        case accel_gpu_default:
        case accel_gpu_cublas:
        case accel_gpu_cuda:
        case accel_gpu_magma: return true;
        default:
            return false;
    }
}

TypeAcceleration AccelerationMeta::getAvailableFallback(TypeAcceleration accel){
    // sparse grids are evaluated in 2 stages:
    // - s1: convert multi-index to matrix B
    // - s2: multiply matrix B by stored matrix A
    // Mode   | Stage 1 device | Stage 2 device | Library for stage 2
    // CUBLAS |      CPU       |     GPU        | Nvidia cuBlas (or cuSparse)
    // CUDA   |      GPU       |     GPU        | Nvidia cuBlas (or cuSparse)
    // MAGMA  |      GPU*      |     GPU        | UTK magma and magma_sparse
    // BLAS   |      CPU       |     CPU        | BLAS
    // none   | all done on CPU, still using OpenMP (if available)
    // note that CUDA, HIP and DPCPP are interchangeable based on the selected backend at compiler time

    #ifdef Tasmanian_ENABLE_DPCPP
    // temporary workaround
//     if (accel == accel_gpu_magma or accel == accel_gpu_cuda)
//         return accel_gpu_cublas;
//     return accel;
    #endif

    // accel_gpu_default should always point to the potentially "best" option (currently MAGMA)
    if (accel == accel_gpu_default) accel = accel_gpu_magma;
    #if !defined(Tasmanian_ENABLE_GPU) || !defined(Tasmanian_ENABLE_MAGMA) || !defined(Tasmanian_ENABLE_BLAS)
    // if any of the 3 acceleration modes is missing, then add a switch statement to guard against setting that mode
    switch(accel){
        #ifndef Tasmanian_ENABLE_GPU
        // if CUDA is missing: just use the CPU
        case accel_gpu_cublas:
        case accel_gpu_cuda:
            #ifdef Tasmanian_ENABLE_BLAS
            accel = accel_cpu_blas;
            #else
            accel = accel_none;
            #endif
            break;
        #endif // Tasmanian_ENABLE_GPU
        #ifndef Tasmanian_ENABLE_MAGMA
        // MAGMA tries to use CUDA kernels with magma linear algebra, this CUDA is the next best thing
        case accel_gpu_magma:
            #ifdef Tasmanian_ENABLE_GPU
            accel = accel_gpu_cuda;
            #elif defined(Tasmanian_ENABLE_BLAS)
            accel = accel_cpu_blas;
            #else
            accel = accel_none;
            #endif
            break;
        #endif // Tasmanian_ENABLE_MAGMA
        #ifndef Tasmanian_ENABLE_BLAS
        // if BLAS is missing, do not attempt to use the GPU but go directly to "none" mode
        case accel_cpu_blas:
            accel = accel_none;
            break;
        #endif // Tasmanian_ENABLE_BLAS
        default: // compiler complains if there is no explicit "default", even if empty
            break;
    }
    #endif
    return accel;
}

AccelerationDomainTransform::AccelerationDomainTransform(AccelerationContext const *acc, std::vector<double> const &transform_a, std::vector<double> const &transform_b) : acceleration(acc){
    // The points are stored contiguously in a vector with stride equal to num_dimensions
    // Using the contiguous memory in a contiguous fashion on the GPU implies that thread 0 works on dimension 0, thread 1 on dim 1 ...
    // But the number of dimensions is often way less than the number of threads
    // Therefore, we lump vectors together into large vectors of sufficient dimension
    // The dimension is least 512, but less than max CUDA threads 1024
    // The domain transforms are padded accordingly
    num_dimensions = (int) transform_a.size();
    padded_size = num_dimensions;
    while(padded_size < 512) padded_size += num_dimensions;

    std::vector<double> rate(padded_size);
    std::vector<double> shift(padded_size);
    int c = 0;
    for(int i=0; i<padded_size; i++){
        // instead of storing upper/lower limits (as in TasmanianSparseGrid) use rate and shift
        double diff = transform_b[c] - transform_a[c];
        rate[i] = 2.0 / diff;
        shift[i] = (transform_b[c] + transform_a[c]) / diff;
        c++;
        c = (c % num_dimensions);
    }

    gpu_trans_a.load(acc, rate);
    gpu_trans_b.load(acc, shift);
}

template<typename T>
void AccelerationDomainTransform::getCanonicalPoints(bool use01, const T *gpu_transformed_x, int num_x, GpuVector<T> &gpu_canonical_x){
    gpu_canonical_x.resize(acceleration, ((size_t) num_dimensions) * ((size_t) num_x));
    TasGpu::dtrans2can(acceleration, use01, num_dimensions, num_x, padded_size, gpu_trans_a.data(), gpu_trans_b.data(), gpu_transformed_x, gpu_canonical_x.data());
}

template void AccelerationDomainTransform::getCanonicalPoints<float>(bool, float const[], int, GpuVector<float>&);
template void AccelerationDomainTransform::getCanonicalPoints<double>(bool, double const[], int, GpuVector<double>&);

}

#endif
