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
 * and RAII resource management for the opaque pointer handles for TPL such as cuBlas and rocBlas.
 * The delete operation is defined in the corresponding TPL file.
 * The deleter is initialized with ownership flag which is indicated whether the associated
 * handle should or should not be deleted, e.g., whether the handle was created internally
 * or handed down by the user.
 *
 * Effectively, the deleter with init_own set to false will turn the associated std::unique_ptr
 * into a non-owning container. The positive side of this approach include:
 * - the rest of the code can refer to a single unique_ptr, e.g., this bundles the ownership flag
 *   into the object
 * - the rest of the code can utilize all properties of std::unique_ptr, e.g., automatically
 *   initializing and creating/deleting the proper copy and move constructor and assignment operators
 * - the C++ code outside of the TasGrid::deleteHandle() specializations looks like modern C++
 *
 * The negative side of the above is that a non-owning std::unique_ptr is unexpected and violates
 * the spirit of the container, and possibly the whole point of the container.
 * However, this hack is necessitated by the mix of C-style of memory management for cuBlas/rocBlas
 * handles together with the C++ code. The alternative to this hack is to follow a C-style of code
 * throughout and keeping separate variables to indicate ownership, manually deleting constructors,
 * manually implementing move semantics, and losing C++ abstraction and automation.
 * The loss of abstraction would be obvious to a causal reader, but the code will be prone to all
 * problems that come without RAII regardless of how diligent or how familiar with Tasmanian
 * a contributor may be.
 * Therefore, this deleter is implemented to allow the use of std::unique_ptr in either owning or non-owning mode.
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
