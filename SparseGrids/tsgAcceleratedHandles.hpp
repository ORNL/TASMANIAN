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

#ifndef __TASMANIAN_SPARSE_GRID_ACCELERATED_HANDLES
#define __TASMANIAN_SPARSE_GRID_ACCELERATED_HANDLES

#include "tsgEnumerates.hpp"

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianAcceleration
 * \brief Type tags for the different handle types.
 *
 * \endinternal
 */
namespace AccHandle{
    //! \brief cuBlas handle.
    struct Cublas{};
    //! \brief cuSparse handle.
    struct Cusparse{};
    //! \brief cuSolver handle.
    struct Cusolver{};
    //! \brief rocBlas handle.
    struct Rocblas{};
    //! \brief rocSparse handle.
    struct Rocsparse{};
    //! \brief SYCL queue handle.
    struct Syclqueue{};
}

/*!
 * \internal
 * \ingroup TasmanianAcceleration
 * \brief Deletes the handle, specialized for each TPL backend and tag in TasGrid::AccHandle namepace.
 *
 * \endinternal
 */
template<typename> void deleteHandle(int*);

/*!
 * \internal
 * \ingroup TasmanianAcceleration
 * \brief Deleter template for the GPU handles, e.g., cuBlas and rocBlas
 *
 * The deleter is used in conjunction with std::unique_ptr to provide both type-safety
 * and RAII resource management for the opaque pointer handles.
 * The delete operation is defined in the corresponding TPL file.
 * \endinternal
 */
template<typename ehandle>
struct HandleDeleter{
    //! \brief Constructor indicating the ownership of the handle.
    HandleDeleter(bool init_own = true) : own(init_own){}
    //! \brief Deletes the handle, if owned; nothing otherwise.
    void operator()(int *p) const { if (own) deleteHandle<ehandle>(p); }
    //! \brief Save the own/not own state.
    bool own;
};

}

#endif
