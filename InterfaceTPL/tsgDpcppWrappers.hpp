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
 * \brief Derives from sycl::device_selector and creates a list of device names to be used in Tasmanian.
 *
 * \endinternal
 */
class tsg_gpu_selector : public sycl::device_selector{
public:
    //! \brief Default constructor.
    tsg_gpu_selector() : deviceID(-1), has_gpu(false){}
    //! \brief Constructor that selects a specific GPU device.
    tsg_gpu_selector(int device) : deviceID(device), has_gpu(false){}
    //! \brief Used during the device selection, also populates the two lists.
    int operator()(const sycl::device &device) const override{
        if (device.is_gpu()){
            names.push_back(device.get_info<sycl::info::device::name>());
            memory.push_back(device.get_info<sycl::info::device::global_mem_size>());
            has_gpu = true;
        }
        if (deviceID > -1 and static_cast<size_t>(deviceID + 1) == names.size())
            return 100; // if a specific device is desired, mark it high
        else
            return 1;
    }
    //! \brief The ID of the GPU device to be selected.
    int const deviceID;
    //! \brief Returns true if a GPU device has been found.
    mutable bool has_gpu;
    //! \brief Holds a list of the device names.
    mutable std::vector<std::string> names;
    //! \brief Holds a list of the device memory.
    mutable std::vector<unsigned long long> memory;
};
/*!
 * \internal
 * \ingroup TasmanianTPLWrappers
 * \brief Creates a new tsg_gpu_selector and populates it with the names and memory for each device.
 *
 * \endinternal
 */
inline tsg_gpu_selector readSyclDevices(){
    tsg_gpu_selector selector;
    sycl::queue q(selector);
    if (not selector.has_gpu){ // add the default CPU device
        q = sycl::queue();
        selector.names.push_back(q.get_device().get_info<sycl::info::device::name>());
        selector.memory.push_back(q.get_device().get_info<sycl::info::device::global_mem_size>());
    }
    return selector;
}

/*!
 * \internal
 * \ingroup TasmanianTPLWrappers
 * \brief Return a new sycl::queue wrapped in a unique_ptr object in owning mode and associated with the given deviceID.
 *
 * \endinternal
 */
inline std::unique_ptr<int, HandleDeleter<AccHandle::Syclqueue>> makeNewQueue(int deviceID){
    sycl::queue *qq = nullptr;
    if (readSyclDevices().has_gpu){
        tsg_gpu_selector selector(deviceID);
        qq = new sycl::queue(selector);
    }else{
        qq = new sycl::queue();
    }
    return std::unique_ptr<int, HandleDeleter<AccHandle::Syclqueue>>(
                                    reinterpret_cast<int*>(qq),
                                    HandleDeleter<AccHandle::Syclqueue>()
                                );
}

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
    if (not acceleration->engine->internal_queue){
        if (test_queue.use_testing){
            acceleration->engine->internal_queue = test_queue; // take non-owning copy of the pointer
        }else{
            acceleration->engine->internal_queue = makeNewQueue(acceleration->device);
        }
    }
    return reinterpret_cast<sycl::queue*>(acceleration->engine->internal_queue.get());
}

}

#endif
