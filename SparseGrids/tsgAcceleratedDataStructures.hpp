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

#include "tsgAcceleratedHandles.hpp"

//! \internal
//! \file tsgAcceleratedDataStructures.hpp
//! \brief Data structures for interacting with CUDA and MAGMA environments.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianAcceleration
//!
//! Classes and namespaces that wrap around basic CUDA and MAGMA functionality,
//! the classes allow RAII style of memory management for CUDA GPU arrays,
//! as well as handy encapsulation of the cuBlas/cuSparse/MAGMA handles and streams.

/*!
 * \internal
 * \ingroup TasmanianSG
 * \addtogroup TasmanianAcceleration Classes and functions used for acceleration methods
 *
 * \par RAII Memory Management
 * CUDA uses C-style of memory management with cudaMalloc(), cudaMemcopy(), cudaFree(),
 * but templated C++ std::vector-style class is far more handy and more fail-safe.
 * The \b GpuVector template class guards against memory leaks and offers more seamless
 * integration between CPU and GPU data structures.
 * See the \b GpuVector documentation for details.
 *
 * \par Streams and Handles Encapsulation
 * CUDA linear algebra libraries (as well as MAGAM), use streams and handles for all their calls.
 * The handles have to be allocated, deleted, and passed around which causes unnecessary code clutter.
 * Encapsulating the handles in a single \b GpuEngine class greatly simplifies the work-flow.
 * Furthermore, some (sparse) linear operations require multiple calls to CUDA/MAGMA libraries,
 * and it is easier to combine those into a single call to a \b GpuEngine method.
 *
 * \par Acceleration Metadata
 * The \b AccelerationMeta namespace offers several methods used throughout the library and in the testing:
 * - Tasmanian specific acceleration fallback logic
 * - Reading CUDA device properties, e.g., number of devices or total memory
 * - Error handling for common CUDA/cuBlas/cuSparse calls
 *
 * \par C++ Wrappers to Fortran BLAS API
 * The standard BLAS API follows Fortran calling conventions,
 * e.g., call by value and underscore at the end of function names.
 * A C++ wrapper is provided that handles Tasmanian specific cases of
 * dense matrix-matrix and matrix-vector multiplication using C++ compatible API.
 * \endinternal
 */

namespace TasGrid{

struct AccelerationContext; // forward declaration, CUDA and HIP GpuVector do not use the context, but DPC++ needs it

/*!
 * \ingroup TasmanianAcceleration
 * \brief Template class that wraps around a single GPU array, providing functionality that mimics std::vector.
 *
 * \par Wraps Around a GPU Array
 * The class can be instantiated with either \b int, \b float or \b double (other types are not currently available).
 * The class can either allocate (and deallocate with the descructor) an array of desired size that resided in the GPU device memory,
 * and load/unload data to a CPU std::vector or array with the same type.
 *
 * Note that the class does not provide a single entry access and it is not copyable, only movable in assignment and constructor.
 */
template<typename T>
class GpuVector{
public:
    //! \brief Delete the copy-constructor.
    GpuVector(GpuVector<T> const &) = delete;
    //! \brief Delete the copy-assignment.
    GpuVector<T>& operator =(GpuVector<T> const &) = delete;

    //! \brief Allow for move-construction.
    GpuVector(GpuVector<T> &&other) : num_entries(Utils::exchange(other.num_entries, 0)), gpu_data(Utils::exchange(other.gpu_data, nullptr))
    #ifdef Tasmanian_ENABLE_DPCPP
        , sycl_queue(other.sycl_queue)
    #endif
    {}
    //! \brief Allow for move-assignment.
    GpuVector<T>& operator =(GpuVector<T> &&other){
        GpuVector<T> temp(std::move(other));
        std::swap(num_entries, temp.num_entries);
        std::swap(gpu_data, temp.gpu_data);
        #ifdef Tasmanian_ENABLE_DPCPP
        std::swap(sycl_queue, temp.sycl_queue);
        #endif
        return *this;
    }

    //! \brief Default constructor, creates an empty (null) array.
    GpuVector() : num_entries(0), gpu_data(nullptr){}
    //! \brief Construct a vector with \b count number of entries.
    GpuVector(AccelerationContext const *acc, size_t count) : num_entries(0), gpu_data(nullptr){ resize(acc, count); }

    /*!
     * \brief Same as \b GpuVector(dim1 * dim2), but guards against overflow.
     *
     * Many of the Tasmanian data-structures are inherently two-dimensional,
     * for example, passing number of points and number of dimensions separately makes the code more readable,
     * and both integers are converted to size_t before multiplication which prevents overflow.
     * Note: the dimensions \b will \b not be stored, the underlying data is still one dimensional.
     */
    GpuVector(AccelerationContext const *acc, int dim1, int dim2) : num_entries(0), gpu_data(nullptr){ resize(acc, Utils::size_mult(dim1, dim2)); }
    //! \brief Create a vector with size that matches \b cpu_data and copy the data to the GPU device.
    GpuVector(AccelerationContext const *acc, const std::vector<T> &cpu_data) : num_entries(0), gpu_data(nullptr){ load(acc, cpu_data); }
    //! \brief Construct a vector and load with date provided on to the cpu.
    GpuVector(AccelerationContext const *acc, int dim1, int dim2, T const *cpu_data) : num_entries(0), gpu_data(nullptr){ load(acc, Utils::size_mult(dim1, dim2), cpu_data); }
    //! \brief Construct a vector by loading from a given range.
    template<typename IteratorLike> GpuVector(AccelerationContext const *acc, IteratorLike ibegin, IteratorLike iend) : GpuVector(){ load(acc, ibegin, iend); }
    //! \brief Destructor, release all allocated memory.
    ~GpuVector(){ clear(); }

    //! \brief Return the current size of the GPU array.
    size_t size() const{ return num_entries; }
    //! \brief Get a reference to the GPU array, which an be used as input to GPU libraries and kernels.
    T* data(){ return gpu_data; }
    //! \brief Get a const-reference to the GPU array, which an be used as input to GPU libraries and kernels.
    const T* data() const{ return gpu_data; }

    //! \brief Clear all data currently stored in the vector and allocate a new array (unlike std::vector this does not relocate the old data).
    void resize(AccelerationContext const *acc, size_t count);
    //! \brief Delete all allocated memory and reset the array to empty.
    void clear();
    //! \brief Return \b true if the \b size() is zero.
    bool empty() const{ return (num_entries == 0); }

    //! \brief Copy the content of \b cpu_data to the GPU device, all pre-existing data is deleted and the vector is resized to match \b cpu_data.
    void load(AccelerationContext const *acc, const std::vector<T> &cpu_data){ load(acc, cpu_data.size(), cpu_data.data()); }

    //! \brief Load from a range defined by the begin and end, converts if necessary.
    template<typename IteratorLike>
    void load(AccelerationContext const *acc, IteratorLike ibegin, IteratorLike iend){
        load(acc, std::distance(ibegin, iend), &*ibegin);
    }

    /*!
     * \brief Takes a vector with entries of different precision, converts and loads.
     *
     * Used when the CPU vectors are stored in double-precision format while the GPU entries are prepared to work with single-precision.
     */
    template<typename U>
    Utils::use_if<!std::is_same<U, T>::value> load(AccelerationContext const *acc, const std::vector<U> &cpu_data){
        load(acc, cpu_data.size(), cpu_data.data());
    }

    /*!
     * \brief Copy the first \b count entries of \b cpu_data to the GPU device.
     *
     * If \b count does not match the current size, the current array will be deleted and new array will be allocated
     * (even if \b size() exceeds \b count).
     * The final vector size will match \b count.
     * However, if \b count matches \b size(), the data is overwritten but no arrays will be allocated, deleted, or de-aliased.
     */
    void load(AccelerationContext const *acc, size_t count, const T* cpu_data);
    /*!
     * \brief Takes a vector with entries of different precision, converts and loads.
     *
     * Used when the CPU data is stored in double-precision format while the GPU entries are prepared to work with single-precision.
     */
    template<typename U>
    Utils::use_if<!std::is_same<U, T>::value> load(AccelerationContext const *acc, size_t count, const U* cpu_data){
        std::vector<T> converted(count);
        std::transform(cpu_data, cpu_data + count, converted.begin(), [](U const &x)->T{ return static_cast<T>(x); });
        load(acc, converted);
    }
    //! \brief Copy the data from the GPU array to \b cpu_data, the \b cpu_data will be resized and overwritten.
    void unload(AccelerationContext const *acc, std::vector<T> &cpu_data) const{
        cpu_data.resize(num_entries);
        unload(acc, cpu_data.data());
    }
    //! \brief Return a CPU vector holding the data of the GPU.
    std::vector<T> unload(AccelerationContext const *acc) const{
        std::vector<T> y;
        unload(acc, y);
        return y;
    }
    //! \brief Copy the first \b num entries to the \b cpu_data buffer, assumes that the buffer is sufficiently large.
    void unload(AccelerationContext const *acc, size_t num, T* cpu_data) const;
    //! \brief Copy the data from the GPU array to the \b cpu_data buffer, assumes that the buffer is sufficiently large.
    void unload(AccelerationContext const *acc, T* cpu_data) const{ unload(acc, num_entries, cpu_data); }

    //! \brief Move the data to the \b external array, the vector is set to empty (unlike move command on std::vector).
    T* eject(){
        T* external = gpu_data;
        gpu_data = nullptr;
        num_entries = 0;
        return external;
    }

    //! \brief The data-type of the vector entries.
    using value_type = T;

private:
    size_t num_entries; // keep track of the size, update on every call that changes the gpu_data
    T *gpu_data; // the GPU array
    #ifdef Tasmanian_ENABLE_DPCPP
    void* sycl_queue;
    #endif
};

/*!
 * \ingroup TasmanianAcceleration
 * \brief Wrapper class around calls GPU accelerated linear algebra libraries.
 *
 * The class also manages the required handles and queues and holds the context of the active GPU device.
 */
struct GpuEngine{
    #ifdef Tasmanian_ENABLE_CUDA
    //! \brief Manually sets the cuBlas handle, handle must be a valid cublasHandle_t associated with this CUDA device.
    void setCuBlasHandle(void *handle);
    //! \brief Manually sets the cuSparse handle, handle must be a valid cusparseHandle_t associated with this CUDA device.
    void setCuSparseHandle(void *handle);
    //! \brief Manually sets the cuSparse handle, handle must be a valid cusolverDnHandle_t associated with this CUDA device.
    void setCuSolverDnHandle(void *handle);

    //! \brief Holds the cuBlas handle.
    std::unique_ptr<int, HandleDeleter<AccHandle::Cublas>> cublas_handle;
    //! \brief Holds the cuSparse handle.
    std::unique_ptr<int, HandleDeleter<AccHandle::Cusparse>> cusparse_handle;
    //! \brief Holds the cuSolver handle.
    std::unique_ptr<int, HandleDeleter<AccHandle::Cusolver>> cusolver_handle;
    #endif

    #ifdef Tasmanian_ENABLE_HIP
    //! \brief Set the rocBlas handle, handle must be a valid rocblas_handle associated with this ROCm device.
    void setRocBlasHandle(void *handle);
    //! \brief Set the rocSparse handle, handle must be a valid rocsparse_handle associated with this ROCm device.
    void setRocSparseHandle(void *handle);

    //! \brief Holds the rocBlas handle.
    std::unique_ptr<int, HandleDeleter<AccHandle::Rocblas>> rblas_handle;
    //! \brief Holds the rocSparse handle.
    std::unique_ptr<int, HandleDeleter<AccHandle::Rocsparse>> rsparse_handle;
    #endif

    #ifdef Tasmanian_ENABLE_DPCPP
    //! \brief Set a user provided sycl::queue.
    void setSyclQueue(void *queue);
    //! \brief Holds the actual queue.
    std::unique_ptr<int, HandleDeleter<AccHandle::Syclqueue>> internal_queue;
    #endif

    //! \brief Avoids an empty engine when no acceleration is enabled, allows for default constructor/move/copy, skips extraneous calls to MAGMA init.
    std::unique_ptr<int> called_magma_init;
};

//! \internal
//! \brief Implements the domain transform algorithms in case the user data is provided on the GPU.
//! \ingroup TasmanianAcceleration

//! Takes the upper and lower bounds of a hypercube and transforms the user provided points to the canonical domain (-1, 1) or (0, 1).
//! The transformation is done on the GPU to avoid extraneous data movement.
//!
//! \b Note: Conformal mapping and the non-linear Gauss-Hermite and Gauss-Laguerre transforms are not supported.
class AccelerationDomainTransform{
public:
    //! \brief Constructor, load the transform data to the GPU, the vectors are the same as used in the \b TasmanianSparseGrid class.
    AccelerationDomainTransform(AccelerationContext const *, std::vector<double> const &transform_a, std::vector<double> const &transform_b);

    /*!
     * \brief Transform a set of points, used in the calls to \b evaluateHierarchicalFunctionsGPU()
     *  Takes the user provided \b gpu_transformed_x points of dimension matching the grid num_dimensions
     * and total number \b num_x.
     * The \b gpu_canonical_x is resized to match \b gpu_transformed_x and it loaded with the corresponding canonical points.
     * The \b use01 flag indicates whether to use canonical domain (0, 1) (Fourier grids), or (-1, 1) (almost everything else).
     */
    template<typename T>
    void getCanonicalPoints(bool use01, T const gpu_transformed_x[], int num_x, GpuVector<T> &gpu_canonical_x);

private:
    // these actually store the rate and shift and not the hard upper/lower limits
    GpuVector<double> gpu_trans_a, gpu_trans_b;
    int num_dimensions, padded_size;
    AccelerationContext const *acceleration;
};

//! \internal
//! \brief Wrappers around custom CUDA kernels to handle domain transforms and basis evaluations, the kernels are instantiated in tsgCudaKernels.cu
//! \ingroup TasmanianAcceleration
namespace TasGpu{
    //! \internal
    //! \brief Uses custom kernel to convert \b transformed points to \b canonical points, all arrays live on the CUDA device.
    //! \ingroup TasmanianAcceleration

    //! The array \b gpu_x_transformed is provided by the user and must be of size \b num_x times \b dims.
    //! The points are stored contiguously in strides of \b dim (identical to all other calls).
    //! In order to facilitate contiguous memory access, it is most efficient to assign each thread to different dimension,
    //! but the dimensions are much less than the threads.
    //! Thus, we pad the transforms \b gpu_trans_a and \b gpu_trans_b as if they are applied to much larger vectors (of dimension \b pad_size).
    //! The \b pad_size must be a multiple of \b dims, but \b num_x doesn't have to divide into \b pad_size.
    //! The \b gpu_trans_a and \b gpu_trans_b must have length \b pad_size and must contain the rate and shift of the transform (\b not the upper/lower limits).
    //!
    //! The \b use01 indicates whether to use canonical interval (0, 1) or (-1, 1).
    //! This is called from \b AccelerationDomainTransform::getCanonicalPoints().
    //!
    //! See the implementation of \b AccelerationDomainTransform::load() for more details.
    template<typename T>
    void dtrans2can(AccelerationContext const *acc, bool use01, int dims, int num_x, int pad_size,
                    const double *gpu_trans_a, const double *gpu_trans_b,
                    const T *gpu_x_transformed, T *gpu_x_canonical);

    //! \internal
    //! \brief Evaluate the basis functions for a local polynomial grid using the \b DENSE algorithm.
    //! \ingroup TasmanianAcceleration

    //! Used for basis evaluations in \b GridLocalPolynomial with \b order (0, 1, and 2) and \b rule.
    //! The grid has \b num_dimensions and \b num_basis functions (equal to the number of points in the grid).
    //! The number of locations to evaluate is \b num_x, and \b gpu_x must have size \b num_x times \b num_dimensions, and \b gpu_x must be on a canonical interval.
    //! The output matrix \b gpu_y has dimension \b num_basis (contiguous dimension) times \b num_x.
    //!
    //! The grid nodes and the associated support are "encoded" in \b gpu_nodes and \b gpu_support,
    //! where negative support is used to handle special cases such as global support on level 1 (rule semi-localp).
    //! The interpretation of the negative support must be synchronized between the kernel \b tasgpu_devalpwpoly_feval() and \b GridLocalPolynomial::encodeSupportForGPU().
    template<typename T>
    void devalpwpoly(AccelerationContext const *acc, int order, TypeOneDRule rule, int num_dimensions, int num_x, int num_basis, const T *gpu_x, const T *gpu_nodes, const T *gpu_support, T *gpu_y);

    //! \internal
    //! \brief Evaluate the basis functions for a local polynomial grid using the \b SPARSE algorithm.
    //! \ingroup TasmanianAcceleration

    //! The inputs are identical to \b devalpwpoly() with the addition of hierarchy vectors \b gpu_hpntr, \b gpu_hindx, and \b gpu_roots.
    //! The hierarchy vectors define a series of trees, one for each entry in \b gpu_roots, and all three nodes combined are indexed from 0 to gpu_spntr.size().
    //! The \b gpu_hpntr holds the offsets of the children of each node and the indexes of the children are stored in \b gpu_sindx.
    //! The format is identical to row compressed sparse matrix.
    //!
    //! The output vectors \b gpu_spntr, \b gpu_sindx and \b gpu_svals form a row compressed matrix,
    //! e.g., in format that can directly interface with Nvidia cusparseDcsrmm2().
    template<typename T>
    void devalpwpoly_sparse(AccelerationContext const *acc, int order, TypeOneDRule rule, int dims, int num_x, const T *gpu_x,
                            const GpuVector<T> &gpu_nodes, const GpuVector<T> &gpu_support,
                            const GpuVector<int> &gpu_hpntr, const GpuVector<int> &gpu_hindx, const GpuVector<int> &gpu_hroots,
                            GpuVector<int> &gpu_spntr, GpuVector<int> &gpu_sindx, GpuVector<T> &gpu_svals);

    //! \internal
    //! \brief Evaluate the basis for a Sequence grid.
    //! \ingroup TasmanianAcceleration

    //! The evaluation points are defined on a canonical interval and given in \b gpu_x with size \b dims  times \b num_x.
    //! The \b max_levels indicates the maximum level in each direction (one more than the vector stored in the GridSequence data structure),
    //! the vector is required to compute the offsets of the intermediate cache for the Newton polynomials.
    //! The \b points holds the same information as the \b MultiIndexSet, but in transposed order, i.e., the dimensions of the multi-indexes are contiguous.
    //! The kernel handles one dimension at a time, hence the switched order compared to the CPU which handles one multi-index at a time.
    //! The \b ndoes vector has the cached 1D nodes and \b coeffs holds the cached coefficients of the Newton polynomials.
    //!
    //! The output is \b gpu_result which must have dimension \b num_x by \b num_nodes.size() / \b dims.
    template<typename T>
    void devalseq(AccelerationContext const *acc, int dims, int num_x, const std::vector<int> &max_levels, const T *gpu_x,
                  const GpuVector<int> &num_nodes,
                  const GpuVector<int> &points, const GpuVector<T> &nodes, const GpuVector<T> &coeffs, T *gpu_result);

    //! \internal
    //! \brief Evaluate the basis for a Fourier grid.
    //! \ingroup TasmanianAcceleration

    //! The logic is identical to \b devalseq(), except the Fourier polynomials do not require nodes or coefficients.
    //! The output is two real arrays of size \b num_x by \b num_nodes.size() / \b dims corresponding to the real and complex parts of the basis.
    template<typename T>
    void devalfor(AccelerationContext const *acc, int dims, int num_x, const std::vector<int> &max_levels, const T *gpu_x, const GpuVector<int> &num_nodes, const GpuVector<int> &points, T *gpu_wreal, typename GpuVector<T>::value_type *gpu_wimag);

    /*!
     * \internal
     * \ingroup TasmanianAcceleration
     * \brief Evaluate the basis for Global grid.
     *
     * The logic is more complicated due to the more general nature of the grid.
     * \endinternal
     */
    template<typename T>
    void devalglo(AccelerationContext const *acc, bool is_nested, bool is_clenshawcurtis0, int dims, int num_x, int num_p, int num_basis,
                  T const *gpu_x, GpuVector<T> const &nodes, GpuVector<T> const &coeff, GpuVector<T> const &tensor_weights,
                  GpuVector<int> const &nodes_per_level, GpuVector<int> const &offset_per_level, GpuVector<int> const &map_dimension, GpuVector<int> const &map_level,
                  GpuVector<int> const &active_tensors, GpuVector<int> const &active_num_points, GpuVector<int> const &dim_offsets,
                  GpuVector<int> const &map_tensor, GpuVector<int> const &map_index, GpuVector<int> const &map_reference, T *gpu_result);


    /*!
     * \ingroup TasmanianAcceleration
     * \brief Fills the \b data with the provided real number at the given stride.
     */
    void fillDataGPU(AccelerationContext const *acc, double value, long long N, long long stride, double data[]);

    /*!
     * \ingroup TasmanianAcceleration
     * \brief Similar to copy_n, copies the data from the CPU to the GPU.
     */
    template<typename T> void load_n(AccelerationContext const *acc, T const *cpu_data, size_t num_entries, T *gpu_data);

    /*!
     * \ingroup TasmanianAcceleration
     * \brief Similar to copy_n, copies the data from the CPU to the GPU.
     */
    template<typename T, typename U>
    Utils::use_if<!std::is_same<U, T>::value> load_n(AccelerationContext const *acc, U const *cpu_data, size_t num_entries, T *gpu_data){
        std::vector<T> converted(num_entries);
        std::transform(cpu_data, cpu_data + num_entries, converted.begin(), [](U const &x)->T{ return static_cast<T>(x); });
        load_n(acc, converted.data(), num_entries, gpu_data);
    }

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

//! \internal
//! \brief Common methods for manipulating acceleration options and reading CUDA environment properties.
//! \ingroup TasmanianAcceleration
namespace AccelerationMeta{
    //! \internal
    //! \brief Convert the string (coming from C or Python) into an enumerated type.
    //! \ingroup TasmanianAcceleration
    TypeAcceleration getIOAccelerationString(const char * name);
    //! \brief Creates a map with \b std::string rule names (used by C/Python/CLI) mapped to \b TypeAcceleration enums.
    std::map<std::string, TypeAcceleration> getStringToAccelerationMap();
    //! \internal
    //! \brief Convert the enumerated type to a string, the inverse of \b getIOAccelerationString()
    //! \ingroup TasmanianAcceleration
    const char* getIOAccelerationString(TypeAcceleration accel);
    //! \internal
    //! \brief Convert the integer (coming from Fortran) into an enumerated type.
    //! \ingroup TasmanianAcceleration
    int getIOAccelerationInt(TypeAcceleration accel);
    //! \internal
    //! \brief Convert the enumerated type to an integer, the inverse of \b getIOAccelerationInt()
    //! \ingroup TasmanianAcceleration
    TypeAcceleration getIOIntAcceleration(int accel);
    //! \internal
    //! \brief Returns \b true if \b accele is cuda, cublas or magma.
    //! \ingroup TasmanianAcceleration
    bool isAccTypeGPU(TypeAcceleration accel);

    //! \brief Identifies whether the acceleration mode is available.
    inline bool isAvailable(TypeAcceleration accel){
        switch(accel){
            #ifdef Tasmanian_ENABLE_MAGMA
            case accel_gpu_magma: return true;
            #endif
            #ifdef Tasmanian_ENABLE_CUDA
            case accel_gpu_cuda: return true;
            case accel_gpu_cublas: return true;
            #endif
            #ifdef Tasmanian_ENABLE_HIP
            case accel_gpu_hip: return true;
            case accel_gpu_rocblas: return true;
            #endif
            #ifdef Tasmanian_ENABLE_DPCPP
            case accel_gpu_cuda: return true;
            case accel_gpu_cublas: return true;
            #endif
            #ifdef Tasmanian_ENABLE_BLAS
            case accel_cpu_blas: return true;
            #endif
            case accel_none:
                return true;
            default:
                return false;
        }
    }

    //! \internal
    //! \brief Implements fallback logic, if \b accel has been enabled through CMake then this returns \b accel, otherwise it returns the "next-best-thing".
    //! \ingroup TasmanianAcceleration

    //! This function always returns a valid acceleration type.
    //! The fallback logic is documented with the enumerated type \b TasGrid::TypeAcceleration.
    TypeAcceleration getAvailableFallback(TypeAcceleration accel);

    //! \internal
    //! \brief Return the number of visible GPU devices.
    //! \ingroup TasmanianAcceleration
    int getNumGpuDevices();

    //! \internal
    //! \brief Selects the active device for this CPU thread, not supported for DPC++.
    //! \ingroup TasmanianAcceleration

    //! The \b deviceID must be a valid ID (between 0 and getNumGpuDevices() -1).
    void setDefaultGpuDevice(int deviceID);

    //! \internal
    //! \brief Return the memory available in the device (in units of bytes).
    //! \ingroup TasmanianAcceleration

    //! The \b deviceID must be a valid ID (between 0 and getNumGpuDevices() -1).
    unsigned long long getTotalGPUMemory(int deviceID);

    //! \internal
    //! \brief Returns the name of the selected GPU device, empty string if no device is available or the index is out of bounds.
    //! \ingroup TasmanianAcceleration

    //! The \b deviceID must be a valid ID (between 0 and getNumGpuDevices() -1).
    std::string getGpuDeviceName(int deviceID);

    //! \internal
    //! \brief Copy a device array to the main memory, used for testing only, always favor using \b GpuVector (if possible).
    //! \ingroup TasmanianAcceleration
    template<typename T> void recvGpuArray(AccelerationContext const*, size_t num_entries, const T *gpu_data, std::vector<T> &cpu_data);

    //! \internal
    //! \brief Deallocate  device array, used primarily for testing, always favor using \b GpuVector (if possible).
    //! \ingroup TasmanianAcceleration
    template<typename T> void delGpuArray(AccelerationContext const*, T *x);

    /*!
     * \ingroup TasmanianAcceleration
     * \brief Creates a new cuBlas handle, used in unit-testing only.
     */
    void *createCublasHandle();
    /*!
     * \ingroup TasmanianAcceleration
     * \brief Destroys the cuBlas handle, used in unit-testing only.
     */
    void deleteCublasHandle(void *);
}

/*!
 * \internal
 * \ingroup TasmanianAcceleration
 * \brief Wrapper class around GPU device ID, acceleration type and GpuEngine.
 *
 * Single acceleration context held by the TasmanianSparseGrid class and aliased into each of the sub-classes.
 * The context is modified through the main class and is persistent on move, copy and read.
 * The sub-classes hold only a const pointer and can read and use the associated variables and the engine.
 * \endinternal
 */
struct AccelerationContext{
    //! \brief Defines the sparse-dense algorithm flavors, whenever applicable.
    enum AlgorithmPreference{
        //! \brief Use dense algorithm.
        algorithm_dense,
        //! \brief Use sparse algorithm.
        algorithm_sparse,
        //! \brief Use automatically select based on heuristics.
        algorithm_autoselect
    };

    /*!
     * \brief Defines the types of acceleration context updates so they can be linked to acceleration cache updates.
     *
     * The update methods will generate a change event that may require an update of the acceleration cache.
     * For example, switching to a different GPU device means that the cache from the old device must be cleared.
     */
    enum ChangeType{
        //! \brief No change, do nothing.
        change_none,
        //! \brief Change the associated GPU device.
        change_gpu_device,
        //! \brief Change from BLAS or none to a GPU acceleration mode.
        change_gpu_enabled,
        //! \brief Change BLAS to none or none to BLAS.
        change_cpu_blas,
        //! \brief Change the sparse-dense AlgorithmPreference.
        change_sparse_dense
    };

    //! \brief The current active acceleration mode.
    TypeAcceleration mode;
    //! \brief The preference to use dense or sparse algorithms.
    AlgorithmPreference algorithm_select;
    //! \brief If using a GPU acceleration mode, holds the active device.
    int device;

    //! \brief Holds the context to the GPU TPL handles, e.g., MAGMA queue.
    mutable std::unique_ptr<GpuEngine> engine;

    #ifdef Tasmanian_ENABLE_BLAS
    //! \brief Creates a default context, the device id is set to 0 and acceleration is BLAS (if available) or none.
    AccelerationContext() : mode(accel_cpu_blas), algorithm_select(algorithm_autoselect), device(0){}
    #else
    AccelerationContext() : mode(accel_none), algorithm_select(algorithm_autoselect), device(0){}
    #endif

    ~AccelerationContext(){}

    //! \brief Sets algorithm affinity in the direction of sparse.
    ChangeType favorSparse(bool favor){
        AlgorithmPreference new_preference = [=]()->AlgorithmPreference{
            if (favor){
                return (algorithm_select == algorithm_dense) ? algorithm_autoselect : algorithm_sparse;
            }else{
                return (algorithm_select == algorithm_sparse) ? algorithm_autoselect : algorithm_dense;
            }
        }();
        if (new_preference != algorithm_select){
            algorithm_select = new_preference;
            return change_sparse_dense;
        }else{
            return change_none;
        }
    }

    //! \brief Returns true if BLAS is enabled and the current mode is not none.
    bool blasCompatible() const{
        #ifdef Tasmanian_ENABLE_BLAS
        return (mode != accel_none);
        #else
        return false;
        #endif
    }

    //!  \brief Returns true if the current mode implies the use of custom GPU kernels.
    bool useKernels() const{
        #if defined(Tasmanian_ENABLE_CUDA) || defined(Tasmanian_ENABLE_HIP)
        return ((mode == accel_gpu_cuda) or (mode == accel_gpu_magma));
        #else
        return false;
        #endif
    }

    //! \brief Returns the ChangeType if enable() is called, but does not change the acceleration.
    ChangeType testEnable(TypeAcceleration acc, int new_gpu_id) const{
        TypeAcceleration effective_acc = AccelerationMeta::getAvailableFallback(acc);
        if (AccelerationMeta::isAccTypeGPU(effective_acc) and ((new_gpu_id < 0 or new_gpu_id >= AccelerationMeta::getNumGpuDevices())))
            throw std::runtime_error("Invalid GPU device ID, see ./tasgrid -v for list of detected devices.");
        ChangeType mode_change = (effective_acc == mode) ? change_none : change_cpu_blas;
        ChangeType device_change = (device == new_gpu_id) ? change_none : change_gpu_device;

        if (on_gpu()){
            return (AccelerationMeta::isAccTypeGPU(effective_acc)) ? device_change : change_gpu_device;
        }else{
            return (AccelerationMeta::isAccTypeGPU(effective_acc)) ? change_gpu_enabled : mode_change;
        }
    }

    //! \brief Accepts parameters directly from TasmanianSparseGrid::enableAcceleration()
    void enable(TypeAcceleration acc, int new_gpu_id){
        // get the effective new acceleration mode (use the fallback if acc is not enabled)
        TypeAcceleration effective_acc = AccelerationMeta::getAvailableFallback(acc);
        // if switching to a GPU mode, check if the device id is valid
        if (AccelerationMeta::isAccTypeGPU(effective_acc) and ((new_gpu_id < 0 or new_gpu_id >= AccelerationMeta::getNumGpuDevices())))
            throw std::runtime_error("Invalid GPU device ID, see ./tasgrid -v for list of detected devices.");

        if (AccelerationMeta::isAccTypeGPU(effective_acc)){
            // if the new mode is GPU-based, make an engine or reset the engine if the device has changed
            // if the engine exists and the device is not changed, then keep the existing engine
            if (!engine or new_gpu_id != device)
                engine = Utils::make_unique<GpuEngine>();
        }else{
            engine.reset();
        }

        // assign the new values for the mode and device
        mode = effective_acc;
        device = new_gpu_id;
    }
    //! \brief Set default device.
    void setDevice() const{ AccelerationMeta::setDefaultGpuDevice(device); }
    //! \brief Custom convert to \b GpuEngine
    operator GpuEngine* () const{ return engine.get(); }
    //! \brief Returns true if any of the GPU-based acceleration modes have been enabled.
    bool on_gpu() const{ return !!engine; }
};

#ifdef Tasmanian_ENABLE_DPCPP
/*!
 * \internal
 * \ingroup TasmanianAcceleration
 * \brief Class holding the sycl queue used for internal testing.
 *
 * Testing uses multiple Tasmanian objects created and destroyed a lot of those.
 * Creating lots of sycl queues causes the jit system to keep recompiling the same kernels,
 * testing goes better when using a single queue.
 * The queue is created internally and reused for all tests, if InternalSyclQueue::init_testing()
 * has not been called, then regular internal queue would be used (separate for each object).
 * Note that in all cases, the user can provide their own sycl queue.
 * \endinternal
 */
struct InternalSyclQueue{
    //! \brief Default constructor, assume we are not in a testing mode.
    InternalSyclQueue() : use_testing(false){}
    //! \brief Initialize the testing, in which case the internal queue would be used in place of a new queue.
    void init_testing();
    //! \brief Auto-converts to a non-owning std::unique_ptr.
    operator std::unique_ptr<int, HandleDeleter<AccHandle::Syclqueue>> (){
        return std::unique_ptr<int, HandleDeleter<AccHandle::Syclqueue>>(test_queue.get(),
                                                                         HandleDeleter<AccHandle::Syclqueue>(false));
    }
    //! \brief Indicates whether this is a testing run.
    bool use_testing;
    //! \brief Holds the internal sycl::queue for testing.
    std::unique_ptr<int, HandleDeleter<AccHandle::Syclqueue>> test_queue;
};
/*!
 * \internal
 * \ingroup TasmanianAcceleration
 * \brief Holds an internal queue for testing purposes.
 *
 * \endinternal
 */
extern InternalSyclQueue test_queue;
#endif

}

#endif // __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_HPP
