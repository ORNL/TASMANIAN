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

#ifndef __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_HPP
#define __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_HPP

#include <stdexcept>
#include <string>
#include <vector>

#include "tsgEnumerates.hpp"

//! \internal
//! \file tsgAcceleratedDataStructures.hpp
//! \brief Data structures for interacting with CUDA and MAGMA environments.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianAcceleration
//!
//! Classes and namespaces that wrap around basic CUDA and MAGMA functionality,
//! the classes allow RAII style of memory management for CUDA GPU arrays,
//! as well as handy encapsulation of the cuBlas/cuSparse/MAGMA handles and streams.

//! \internal
//! \defgroup TasmanianAcceleration Classes and functions used for acceleration methods.
//!
//! \par RAII Memory Management
//! CUDA uses C-style of memory management with cudaMalloc(), cudaMemcopy(), cudaFree(),
//! but templated C++ std::vector-style class is far more handy.
//! The \b CudaVector template class guards against memory leaks and offers more seamless
//! integration between CPU and GPU data structures.
//! See the \b CudaVector documentation for details.
//!
//! \par Streams and Handles Encapsulation
//! CUDA linear algebra libraries (as well as MAGAM), use streams and handles for all their calls.
//! The handles have to be allocated, deleted, and passed around which causes unnecessary code clutter.
//! Encapsulating the handles in a single \b CudaEngine class greatly simplifies the work-flow.
//! Furthermore, some (sparse) linear operations require multiple calls to CUDA/MAGMA libraries,
//! and it is easier to combine those into a single call to a \b CudaEngine method.
//!
//! \par Acceleration Metadata
//! The \b AccelerationMeta namespace offers several methods used throughout the library and in the testing:
//! - Tasmanian specific acceleration fallback logic
//! - Reading CUDA device properties, e.g., number of devices or total memory
//! - Error handling for common CUDA/cuBlas/cuSparse calls

namespace TasGrid{

#ifdef Tasmanian_ENABLE_CUDA
//! \defgroup TasmanianCudaVector CUDA Vector class for RAII-style of memory management.
//!
//! A template class that allows for RAII style of memory management for CUDA arrays.

#ifndef __TASMANIAN_DOXYGEN_SKIP
// The classes below will be removed later
class cudaInts{
public:
    cudaInts();
    cudaInts(const std::vector<int> &cpu_data);
    ~cudaInts();

    size_t size() const;
    int* data();
    const int* data() const;
    void resize(size_t cnum);
    void clear();

    void load(const std::vector<int> &cpu_data);
    void unload(std::vector<int> &cpu_data) const;

    void eject(int* &destination); // moves the data to the external pointer

private:
    size_t num;
    int *gpu_data;
};

class cudaDoubles{
public:
    cudaDoubles();
    cudaDoubles(size_t cnum);
    cudaDoubles(int a, int b);
    cudaDoubles(int a, int b, const double *cpu_data);
    cudaDoubles(const std::vector<double> &cpu_data);
    ~cudaDoubles();

    size_t size() const;
    double* data();
    const double* data() const;
    void resize(size_t cnum);
    void clear();

    void load(size_t cnum, const double *cpu_data);
    void load(const std::vector<double> &cpu_data);
    void unload(double *cpu_data) const;

    void eject(double* &destination); // moves the data to the external pointer

private:
    size_t num;
    double *gpu_data;
};
#endif

//! \brief Template class that wraps around a single CUDA array, providing functionality that mimics std::vector
//! \ingroup TasmanianCudaVector

//! \b NOTE: This class is not intended to replace CUDA Thrust or CUDA arrays as a whole;
//! the class is primarily used for internal Tasmanian API and can be safely ignored by the majority Tasmanian users.
//! This documentation is provided for developers and can be useful to users that want to do simple data movement between
//! different CUDA enabled libraries without "getting dirty" with non-RAII and (sometimes) non-intuitive CUDA API.
//!
//! \par Wraps Around a CUDA Array
//! The class can be instantiated with either \b int or \b double (other types are not currently available through the external API).
//! The class can either allocate (and deallocate with the descructor) a CUDA array of desired size,
//! and load/unload data to a CPU std::vector with the same type.
//!
//! \par Deliberately Omitted Functionality
//! The class is not intended to do any of the following:
//! - No random memory access to ranges and elements
//! - No algorithm operations (although \b CudaVector.data() can be passes to kernels and functions in place of an array)
//! - No reserve/insert functionality, the size is always the size of the allocated memory
template<typename T>
class CudaVector{
public:
    //! \brief Default constructor, creates an empty (null) array.
    CudaVector() : num_entries(0), dynamic_mode(true), gpu_data(nullptr){}
    //! \brief Construct a vector with \b count number of entries.

    //! Allocates an array that will be automatically deleted
    //! if \b clear() is called, or \b resize() or \b load() functions
    //! are used with size that is different from \b count.
    CudaVector(size_t count) : num_entries(0), dynamic_mode(true), gpu_data(nullptr){ resize(count); }

    //! \brief Same as \b CudaVector(dim1 * dim2), but guards against overflow.

    //! Many of the Tasmanian data-structures are inherently two-dimensional,
    //! for example, passing number of points and number of dimensions separately makes the code more readable,
    //! and both integers are converted to size_t before multiplication which prevents overflow.
    //! Note: the dimensions \b will \b not be stored, the underlying data is still one dimensional.
    CudaVector(int dim1, int dim2) : num_entries(0), dynamic_mode(true), gpu_data(nullptr){ resize(((size_t) dim1) * ((size_t) dim2)); }
    //! \brief Create a vector with size that matches \b cpu_data and copy the data to the CUDA device.
    CudaVector(const std::vector<T> &cpu_data) : num_entries(0), dynamic_mode(true), gpu_data(nullptr){ load(cpu_data); }
    //! \brief Destructor, release all allocated memory.
    ~CudaVector(){ clear(); }

    //! \brief Return the current size of the CUDA array.
    size_t size() const{ return num_entries; }
    //! \brief Get a reference to the CUDA array, which an be used as input to CUDA libraries and kernels.
    T* data(){ return gpu_data; }
    //! \brief Get a const-reference to the CUDA array, which an be used as input to CUDA libraries and kernels.
    const T* data() const{ return gpu_data; }

    //! \brief Clear all data currently stored in the vector and allocate a new array (unlike std::vector this does not copy the data).
    void resize(size_t count);
    //! \brief Delete all allocated memory and reset the array to empty.
    void clear();

    //! \brief Copy the content of \b cpu_data to the CUDA device, all pre-existing data is deleted and the vector is resized to match \b cpu_data.

    //! Makes a call to the other overload with \b cpu_data.size() as \b count. See the overload.
    void load(const std::vector<T> &cpu_data){ load(cpu_data.size(), cpu_data.data()); }
    //! \brief Copy the first \b count entries of \b cpu_data to the CUDA device, all pre-existing data is deleted and the vector is resized to match \b count.

    //! If \b count does not match the current size, the current array will be deleted and new array will be allocated (even if \b size() exceeds \b count).
    //! If the currently stored array has been assigned with \b wrap(), the alias is released.
    //! However, if \b count matches \b size(), the data is overwritten but no arrays will be allocated, deleted, or de-aliased.
    void load(size_t count, const T* cpu_data);
    //! \brief Copy the data from the CUDA array to \b cpu_data, the \b cpu_data will be resized and overwritten, the new size will match the \b size().
    void unload(std::vector<T> &cpu_data) const{
        cpu_data.resize(num_entries);
        unload(cpu_data.data());
    }
    //! \brief Copy the data from the CUDA array to the \b cpu_data buffer, assumes that the buffer is sufficiently large.
    void unload(T* cpu_data) const;

    //! \brief Move the data to the \b external array, the vector is set to empty (unlike move command on std::vector).
    T* eject(){
        T* external = gpu_data;
        gpu_data = nullptr;
        num_entries = 0;
        dynamic_mode = true;
        return external;
    }

    //! \brief Create internal alias to the \b external buffer, assume the buffer has size \b count; the alias will \b not be deleted.

    //! The purpose of the \b wrap() method is to assume control a user allocated CUDA array passed as an input.
    //! That allows to have a consistent internal API using only \b CudaVector classes,
    //! while the external API does not force the user into a Tasmanian specific data-structure (only native CUDA arrays are required).
    //!
    //! Using \b wrap() suppresses the deletion of the data from the destructor, \b clear(), or \b load() (from a vector/array of different size).
    //! The wrapped array has to be deleted through \b external or with \b AccelerationMeta::delCudaArray(\b vector.eject()).
    void wrap(size_t count, T* external){
        gpu_data = external;
        num_entries = count;
        dynamic_mode = false;
    }

private:
    size_t num_entries; // keep track of the size, update on every call that changes the gpu_data
    bool dynamic_mode; // was the memory allocated or wrapped around an exiting object
    T *gpu_data; // the CUDA array
};

#ifndef __TASMANIAN_DOXYGEN_SKIP
// The class below will be removed later
class LinearAlgebraEngineGPU{
public:
    LinearAlgebraEngineGPU();
    ~LinearAlgebraEngineGPU();

    void reset();

    void cublasDGEMM(int M, int N, int K, double alpha, const cudaDoubles &A, const cudaDoubles &B, double beta, cudaDoubles &C);
    void cublasDGEMM(int M, int N, int K, double alpha, const cudaDoubles &A, const std::vector<double> &B, double beta, cudaDoubles &C);
    void cublasDGEMM(int M, int N, int K, double alpha, const cudaDoubles &A, const std::vector<double> &B, double beta, double C[]);
    // dense matrix-matrix (dgemm) or matrix-vector (dgemv for N == 1) product using Nvidai cuBlas

    void cusparseMatmul(int M, int N, int K, double alpha, const cudaDoubles &A, const std::vector<int> &spntr, const std::vector<int> &sindx, const std::vector<double> &svals, double beta, double C[]);
    void cusparseMatmul(int M, int N, int K, double alpha, const cudaDoubles &A, const cudaInts &spntr, const cudaInts &sindx, const cudaDoubles &svals, double beta, cudaDoubles &C);
    // sparse matrix times dense matrix using Nvidia cuSparse (C = alpha * A * (spntr, sindx, svals) + beta *C), the sparse matrix is in column compressed form

    void cusparseMatvec(int M, int N, double alpha, const cudaInts &spntr, const cudaInts &sindx, const cudaDoubles &svals, const cudaDoubles &x, double beta, double y[]);
    // sparse matrix times a dense vector, makes sense only if the matrix already sits on the gpu (y = alpha * A * x + beta * y)

    void cusparseMatveci(int M, int K, double alpha, const cudaDoubles &A, const std::vector<int> &sindx, const std::vector<double> &svals, double beta, double C[]);
    // dense matrix times a sparse vector defined by sindx and svals (C = beta * C + alpha * A * b), currently the sparse vector can only be computed on the cpu

    #ifdef Tasmanian_ENABLE_MAGMA
    void magmaCudaDGEMM(int gpuID, int M, int N, int K, double alpha, const cudaDoubles &A, const cudaDoubles &B, double beta, cudaDoubles &C);
    void magmaCudaDGEMM(int gpuID, int M, int N, int K, double alpha, const cudaDoubles &A, const std::vector<double> &B, double beta, cudaDoubles &C);
    void magmaCudaDGEMM(int gpuID, int M, int N, int K, double alpha, const cudaDoubles &A, const std::vector<double> &B, double beta, double C[]);
    // dense matrix-matrix (dgemm) or matrix-vector (dgemv for N == 1) product using UTK MAGMA
    #endif

protected:
    void makeCuBlasHandle();
    void makeCuSparseHandle();

    #ifdef Tasmanian_ENABLE_MAGMA
    void initializeMagma(int gpuID);
    #endif

private:
    #ifdef Tasmanian_ENABLE_CUDA
    void *cublasHandle;
    void *cusparseHandle;
    #endif

    #ifdef Tasmanian_ENABLE_MAGMA
    bool magma_initialized; // call init once per object (must simplify later)
    void *magmaCudaStream;
    void *magmaCudaQueue;
    #endif
};
#endif

// stores domain transforms, to be used by the top class TasmanianSparseGrid
class AccelerationDomainTransform{
public:
    AccelerationDomainTransform();
    ~AccelerationDomainTransform();

    void clear();
    bool empty();
    void load(const std::vector<double> &transform_a, const std::vector<double> &transform_b);
    void getCanonicalPoints(bool use01, const double *gpu_transformed_x, int num_x, cudaDoubles &gpu_canonical_x);

private:
    // these actually store the rate and shift and not the hard upper/lower limits
    cudaDoubles gpu_trans_a, gpu_trans_b;
    int num_dimensions, padded_size;
};
#endif


#ifdef Tasmanian_ENABLE_CUDA
// namespace realized in tsgCudaKernels.cu, each function corresponds to a CUDA kernel for evaluations of basis matrix, domain transform, or fallback linear algebra
namespace TasCUDA{
    // convert transformed points to the canonical domain, all inputs live on the GPU
    void dtrans2can(bool use01, int dims, int num_x, int pad_size, const double *gpu_trans_a, const double *gpu_trans_b, const double *gpu_x_transformed, double *gpu_x_canonical);

    // evaluate local polynomial rules
    void devalpwpoly(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support, double *gpu_y);
    void devalpwpoly_sparse(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                            int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nzm);
    void devalpwpoly_sparse(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const cudaDoubles &gpu_nodes, const cudaDoubles &gpu_support,
                            const cudaInts &gpu_hpntr, const cudaInts &gpu_hindx, const  cudaInts &gpu_roots, cudaInts &gpu_spntr, cudaInts &gpu_sindx, cudaDoubles &gpu_svals);
    void devalpwpoly_sparse_dense(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                                 int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, double *gpu_dense);

    // evaluata basis functions for sequence rules
    // most data structures are identical to the CPU version, except num_nodes = max_levels + 1, and points is transposed from the IndexSet data
    void devalseq(int dims, int num_x, const std::vector<int> &max_levels, const double *gpu_x, const cudaInts &num_nodes, const cudaInts &points, const cudaDoubles &nodes, const cudaDoubles &coeffs, double *gpu_result);

    // evaluate basis functions for Fourier grids
    // most data structures are identical to the CPU version, except num_nodes = max_levels + 1, and points is transposed from the IndexSet data
    void devalfor(int dims, int num_x, const std::vector<int> &max_levels, const double *gpu_x, const cudaInts &num_nodes, const cudaInts &points, double *gpu_wreal, double *gpu_wimag);


    // #define __TASMANIAN_COMPILE_FALLBACK_CUDA_KERNELS__ // uncomment to compile a bunch of custom CUDA kernels that provide some functionality similar to cuBlas
    #ifdef __TASMANIAN_COMPILE_FALLBACK_CUDA_KERNELS__
    // CUDA kernels that provide essentially the same functionality as cuBlas and MAGMA, but nowhere near as optimal
    // those functions should not be used in a Release or production builds
    // the kernels are useful because they are simple and do not depend on potentially poorly documented 3d party library
    // since the kernels are useful for testing and some debugging, the code should not be deleted (for now), but also don't waste time compiling in most cases

    void cudaDgemm(int M, int N, int K, const double *gpu_a, const double *gpu_b, double *gpu_c);
    // lazy cuda dgemm, nowhere near as powerful as cuBlas, but does not depend on cuBlas
    // gpu_a is M by K, gpu_b is K by N, gpu_c is M by N, all in column-major format
    // on exit gpu_c = gpu_a * gpu_b

    void cudaSparseMatmul(int M, int N, int num_nz, const int* gpu_spntr, const int* gpu_sindx, const double* gpu_svals, const double *gpu_B, double *gpu_C);
    // lazy cuda sparse dgemm, less efficient (especially for large N), but more memory conservative then cusparse as there is no need for a transpose
    // C is M x N, B is K x N (K is max(gpu_sindx)), both are given in row-major format, num_nz/spntr/sindx/svals describe row compressed A which is M by K
    // on exit C = A * B

    void cudaSparseVecDenseMat(int M, int N, int num_nz, const double *A, const int *indx, const double *vals, double *C);
    // dense matrix A (column major) times a sparse vector defiend by num_nz, indx, and vals
    // A is M by N, C is M by 1,
    // on exit C = A * (indx, vals)

    void convert_sparse_to_dense(int num_rows, int num_columns, const int *gpu_pntr, const int *gpu_indx, const double *gpu_vals, double *gpu_destination);
    // converts a sparse matrix to a dense representation (all data sits on the gpu and is pre-allocated)
    #endif
}
#endif

// generic error checking function and types I/O
namespace AccelerationMeta{
    TypeAcceleration getIOAccelerationString(const char * name);
    const char* getIOAccelerationString(TypeAcceleration accel);
    int getIOAccelerationInt(TypeAcceleration accel);
    TypeAcceleration getIOIntAcceleration(int accel);
    bool isAccTypeFullMemoryGPU(TypeAcceleration accel);
    bool isAccTypeGPU(TypeAcceleration accel);

    TypeAcceleration getAvailableFallback(TypeAcceleration accel);

    int getNumCudaDevices();
    void setDefaultCudaDevice(int deviceID);
    unsigned long long getTotalGPUMemory(int deviceID);
    char* getCudaDeviceName(int deviceID);
    template<typename T> void recvCudaArray(size_t num_entries, const T *gpu_data, std::vector<T> &cpu_data);
    template<typename T> void delCudaArray(T *x);
}

}

#endif // __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_HPP
