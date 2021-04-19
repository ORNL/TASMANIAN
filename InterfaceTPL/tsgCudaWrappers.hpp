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

#ifndef __TASMANIAN_CUDA_WRAPPERS_HPP
#define __TASMANIAN_CUDA_WRAPPERS_HPP

#include "tsgGpuWrappers.hpp"

#ifdef Tasmanian_ENABLE_CUDA
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cusolverDn.h>
#else
#error "Cannot use tsgCudaWrappers.cpp without Tasmanian_ENABLE_CUDA"
#endif

#ifdef Tasmanian_ENABLE_MAGMA
#include "tsgMagmaWrappers.hpp"
#endif

/*!
 * \file tsgCudaWrappers.hpp
 * \brief Wrappers to CUDA functionality.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * Helper methods for the CUDA backend.
 */

namespace TasGrid{
namespace TasGpu{

/*
 * Common methods
 */
//! \brief Common error checking, prints a message based on the status and includes the \b info.
inline void cucheck(cudaError_t status, std::string info){
    if (status != cudaSuccess){
        std::string message = "ERROR: cuda failed at ";
        message += info + " with error: ";
        message += cudaGetErrorString(status);
        throw std::runtime_error(message);
    }
}

//! \brief Common error checking, prints a message based on the status and includes the \b info.
inline void cucheck(cublasStatus_t status, std::string info){
    if (status != CUBLAS_STATUS_SUCCESS){
        std::string message = "ERROR: cuBlas failed with code: ";
        if (status == CUBLAS_STATUS_NOT_INITIALIZED){
            message += "CUBLAS_STATUS_NOT_INITIALIZED";
        }else if (status == CUBLAS_STATUS_ALLOC_FAILED){
            message += "CUBLAS_STATUS_ALLOC_FAILED";
        }else if (status == CUBLAS_STATUS_INVALID_VALUE){
            message += "CUBLAS_STATUS_INVALID_VALUE";
        }else if (status == CUBLAS_STATUS_ARCH_MISMATCH){
            message += "CUBLAS_STATUS_ARCH_MISMATCH";
        }else if (status == CUBLAS_STATUS_MAPPING_ERROR){
            message += "CUBLAS_STATUS_MAPPING_ERROR";
        }else if (status == CUBLAS_STATUS_EXECUTION_FAILED){
            message += "CUBLAS_STATUS_EXECUTION_FAILED";
        }else if (status == CUBLAS_STATUS_INTERNAL_ERROR){
            message += "CUBLAS_STATUS_INTERNAL_ERROR";
        }else if (status == CUBLAS_STATUS_NOT_SUPPORTED){
            message += "CUBLAS_STATUS_NOT_SUPPORTED";
        }else if (status == CUBLAS_STATUS_LICENSE_ERROR){
            message += "CUBLAS_STATUS_LICENSE_ERROR";
        }else{
            message += "UNKNOWN";
        }
        message += " at ";
        message += info;
        throw std::runtime_error(message);
    }
}
//! \brief Common error checking, prints a message based on the status and includes the \b info.
inline void cucheck(cusparseStatus_t status, std::string info){
    if (status != CUSPARSE_STATUS_SUCCESS){
        std::string message = "ERROR: cuSparse failed with code: ";
        if (status == CUSPARSE_STATUS_NOT_INITIALIZED){
            message += "CUSPARSE_STATUS_NOT_INITIALIZED";
        }else if (status == CUSPARSE_STATUS_ALLOC_FAILED){
            message += "CUSPARSE_STATUS_ALLOC_FAILED";
        }else if (status == CUSPARSE_STATUS_INVALID_VALUE){
            message += "CUSPARSE_STATUS_INVALID_VALUE";
        }else if (status == CUSPARSE_STATUS_ARCH_MISMATCH){
            message += "CUSPARSE_STATUS_ARCH_MISMATCH";
        }else if (status == CUSPARSE_STATUS_INTERNAL_ERROR){
            message += "CUSPARSE_STATUS_INTERNAL_ERROR";
        }else if (status == CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED){
            message += "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
        }else if (status == CUSPARSE_STATUS_EXECUTION_FAILED){
            message += "CUSPARSE_STATUS_EXECUTION_FAILED";
        }else{
            message += "UNKNOWN";
        }
        message += " at ";
        message += info;
        throw std::runtime_error(message);
    }
}
//! \brief Common error checking, prints a message based on the status and includes the \b info.
inline void cucheck(cusolverStatus_t status, std::string info){
    if (status != CUSOLVER_STATUS_SUCCESS){
        std::string message = "ERROR: cuSolver failed with error code " + std::to_string(status) + " at " + info;
        throw std::runtime_error(message);
    }
}
}

namespace TasGpu{
//! \brief Make cuBlas handle.
inline cublasHandle_t getCuBlasHandle(AccelerationContext const *acceleration){
    if (not acceleration->engine->cublas_handle){
        cublasHandle_t handle;
        cucheck( cublasCreate(&handle), "create handle" );
        acceleration->engine->cublas_handle = std::unique_ptr<int, HandleDeleter<AccHandle::Cublas>>
            (reinterpret_cast<int*>(handle), HandleDeleter<AccHandle::Cublas>());
    }
    return reinterpret_cast<cublasHandle_t>(acceleration->engine->cublas_handle.get());
}
//! \brief Make cuSparse handle.
inline cusparseHandle_t getCuSparseHandle(AccelerationContext const *acceleration){
    if (not acceleration->engine->cusparse_handle){
        cusparseHandle_t handle;
        cucheck( cusparseCreate(&handle), "create handle" );
        acceleration->engine->cusparse_handle = std::unique_ptr<int, HandleDeleter<AccHandle::Cusparse>>
            (reinterpret_cast<int*>(handle), HandleDeleter<AccHandle::Cusparse>());
    }
    return reinterpret_cast<cusparseHandle_t>(acceleration->engine->cusparse_handle.get());
}

//! \brief Wrapper around cusparseMatDescr_t
struct cusparseMatDesc{
    //! \brief An instance of cusparseMatDescr_t
    cusparseMatDescr_t mat_desc;
    //! \brief Create a general matrix description.
    cusparseMatDesc(){
        cucheck( cusparseCreateMatDescr(&mat_desc), "cusparseCreateMatDescr()");
        cusparseSetMatType(mat_desc, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(mat_desc, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatDiagType(mat_desc, CUSPARSE_DIAG_TYPE_NON_UNIT);
    }
    //! \brief Use locally and internally, do not copy or move.
    cusparseMatDesc(cusparseMatDesc const&) = delete;
    //! \brief Release all used memory.
    ~cusparseMatDesc(){
        cusparseDestroyMatDescr(mat_desc);
    }
    //! \brief Automatically convert to the wrapper description.
    operator cusparseMatDescr_t (){ return mat_desc; }
};

#if (CUDART_VERSION >= 11000)
template<typename T>
struct cusparseStructDescription{
    //! \brief Construct a null description.
    cusparseStructDescription() : desc(nullptr){}
    //! \brief Destructor, deletes the description.
    ~cusparseStructDescription(){
        if (std::is_same<T, cusparseSpMatDescr_t>::value){
            if (desc != nullptr) cusparseDestroySpMat(reinterpret_cast<cusparseSpMatDescr_t>(desc));
        }else if (std::is_same<T, cusparseDnVecDescr_t>::value){
            if (desc != nullptr) cusparseDestroyDnVec(reinterpret_cast<cusparseDnVecDescr_t>(desc));
        }else if (std::is_same<T, cusparseDnMatDescr_t>::value){
            if (desc != nullptr) cusparseDestroyDnMat(reinterpret_cast<cusparseDnMatDescr_t>(desc));
        }
    }
    //! \brief Move constructor.
    cusparseStructDescription(cusparseStructDescription &&other) : desc(other.desc){
        other.desc = nullptr;
    }

    //! \brief Returns the description, automatic conversion.
    operator T () const{ return desc; }
    //! \brief Holds the description pointer.
    T desc;
};
template<typename scalar_type>
cusparseStructDescription<cusparseSpMatDescr_t>
makeSparseMatDesc(int rows, int cols, int nnz, int const pntr[], int const indx[], scalar_type const vals[]){
    cusparseStructDescription<cusparseSpMatDescr_t> result;
    cudaDataType_t cuda_type = (std::is_same<scalar_type, float>::value) ? CUDA_R_32F : CUDA_R_64F;
    cucheck(
        cusparseCreateCsr(&result.desc, (int64_t) rows, (int64_t) cols, (int64_t) nnz,
                          const_cast<void*>(reinterpret_cast<void const*>(pntr)),
                          const_cast<void*>(reinterpret_cast<void const*>(indx)),
                          const_cast<void*>(reinterpret_cast<void const*>(vals)),
                          CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, cuda_type),
        "cusparseCreateCsr()"
    );
    return result;
}
template<typename scalar_type>
cusparseStructDescription<cusparseDnVecDescr_t>
makeSparseDenseVecDesc(int size, scalar_type const vals[]){
    cusparseStructDescription<cusparseDnVecDescr_t> result;
    cudaDataType_t cuda_type = (std::is_same<scalar_type, float>::value) ? CUDA_R_32F : CUDA_R_64F;
    cucheck(
        cusparseCreateDnVec(&result.desc, static_cast<int64_t>(size),
                            const_cast<void*>(reinterpret_cast<void const*>(vals)), cuda_type),
        "cusparseCreateDnVec()"
    );
    return result;
}
template<typename scalar_type>
cusparseStructDescription<cusparseDnMatDescr_t>
makeSparseDenseMatDesc(int rows, int cols, int lda, scalar_type const vals[]){
    cusparseStructDescription<cusparseDnMatDescr_t> result;
    cudaDataType_t cuda_type = (std::is_same<scalar_type, float>::value) ? CUDA_R_32F : CUDA_R_64F;
    cucheck(
        cusparseCreateDnMat(&result.desc, static_cast<int64_t>(rows), static_cast<int64_t>(cols), static_cast<int64_t>(lda),
                            const_cast<void*>(reinterpret_cast<void const*>(vals)), cuda_type, CUSPARSE_ORDER_COL),
        "cusparseCreateDnMat()"
    );
    return result;
}
#endif

//! \brief Make cuSolver handle.
inline cusolverDnHandle_t getCuSolverDnHandle(AccelerationContext const *acceleration){
    if (not acceleration->engine->cusolver_handle){
        cusolverDnHandle_t handle;
        cucheck( cusolverDnCreate(&handle), "create handle" );
        acceleration->engine->cusolver_handle = std::unique_ptr<int, HandleDeleter<AccHandle::Cusolver>>
            (reinterpret_cast<int*>(handle), HandleDeleter<AccHandle::Cusolver>());
    }
    return reinterpret_cast<cusolverDnHandle_t>(acceleration->engine->cusolver_handle.get());
}

}
}

#endif
