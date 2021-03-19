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

#ifndef __TASMANIAN_DPCPP_WRAPPERS_HPP
#define __TASMANIAN_DPCPP_WRAPPERS_HPP

#include "tsgGpuWrappers.hpp"

#ifndef Tasmanian_ENABLE_DPCPP
#error "Cannot use tsgDpcppWrappers.cpp without Tasmanian_ENABLE_DPCPP"
#endif

#define MKL_INT int
#include <CL/sycl.hpp>
#include <CL/sycl/usm.hpp>
#include "oneapi/mkl.hpp"

/*!
 * \file tsgDpcppWrappers.hpp
 * \brief Wrappers to DPC++ functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * Helper methods for the DPC++ backend.
 */

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianTPLWrappers
 * \brief Returns the SYCL queue associated with the given AccelerationContext, creates a new queue if needed.
 *
 * If there is a current active sycl::queue the method will return a pointer to the queue,
 * otherwise tries to create a new sycl::queue using the sycl::gpu_selector.
 * If creating the GPU queue fails, e.g., throws sycl::exception,
 * then a queue will be created using the sycl::cpu_selector.
 *
 * \endinternal
 */
inline sycl::queue* getSyclQueue(AccelerationContext const *acceleration){
    if (acceleration->engine->sycl_gpu_queue == nullptr){
        sycl::queue *qq = nullptr;
        try{
            sycl::gpu_selector g_selector;
            qq = new sycl::queue(g_selector);
        }catch(sycl::exception const&){
            sycl::cpu_selector c_selector;
            qq = new sycl::queue(c_selector);
        }

        acceleration->engine->internal_queue = std::shared_ptr<int>(
            reinterpret_cast<int*>(qq),
            [](int* q_int){
                sycl::queue *q = reinterpret_cast<sycl::queue*>(q_int);
                delete q;
            }
        );
        acceleration->engine->sycl_gpu_queue = acceleration->engine->internal_queue.get();
        acceleration->engine->own_gpu_queue = true;
    }
    return reinterpret_cast<sycl::queue*>(acceleration->engine->sycl_gpu_queue);
}

}

#endif
