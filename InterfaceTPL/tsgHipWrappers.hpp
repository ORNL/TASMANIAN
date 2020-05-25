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

#ifndef __TASMANIAN_HIP_WRAPPERS_HPP
#define __TASMANIAN_HIP_WRAPPERS_HPP

#include "tsgGpuWrappers.hpp"

#ifndef Tasmanian_ENABLE_HIP
#error "Cannot use tsgHipWrappers.cpp without Tasmanian_ENABLE_HIP"
#endif

#define __HIP_PLATFORM_HCC__
#include <hip/hip_runtime.h>
#include <rocblas.h>

/*!
 * \file tsgHipWrappers.hpp
 * \brief Wrappers to HIP functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * Helper methods for the HIP backend.
 */

namespace TasGrid{
namespace TasGpu{
/*
 * Common methods
 */
inline void hipcheck(hipError_t status, std::string info){
   if (status != hipSuccess)
       throw std::runtime_error(std::string("rocm-hip encountered an error with code: ") + std::to_string(status) + " at " + info);
}

inline std::string hip_info_from(rocblas_status status){
    switch(status){
        case rocblas_status_success : return "success";
        case rocblas_status_invalid_handle : return "invlida handle";
        case rocblas_status_not_implemented : return "not implemented";
        case rocblas_status_invalid_pointer : return "invalid pointer";
        case rocblas_status_invalid_size : return "invalid size";
        case rocblas_status_memory_error : return "memory error";
        case rocblas_status_internal_error : return "internal error";
        case rocblas_status_perf_degraded : return "performance degraded due to low memory";
        case rocblas_status_size_query_mismatch : return "unmatched start/stop size query";
        case rocblas_status_size_increased : return "queried device memory size increased";
        case rocblas_status_size_unchanged : return "queried device memory size unchanged";
        case rocblas_status_invalid_value : return "passed argument not valid";
        case rocblas_status_continue : return "nothing preventing function to proceed";
        default:
            return "unknown";
    }
}
inline void hipcheck(rocblas_status status, std::string info){
    if (status != rocblas_status_success)
        throw std::runtime_error(std::string("rocblas encountered an error with code: ") + hip_info_from(status) + ", at " + info);
}

//! \brief Make cuBlas handle.
inline rocblas_handle getRocBlasHandle(AccelerationContext const *acceleration){
    if (acceleration->engine->rocblasHandle == nullptr){
        hipcheck( rocblas_create_handle(reinterpret_cast<rocblas_handle*>(&acceleration->engine->rocblasHandle)), "create handle" );
        acceleration->engine->own_rocblas_handle = true;
    }
    return reinterpret_cast<rocblas_handle>(acceleration->engine->rocblasHandle);
}

}
}

#endif
