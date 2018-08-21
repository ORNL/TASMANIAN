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

namespace TasGrid{

#ifdef Tasmanian_ENABLE_CUDA
class cudaInts{
public:
    cudaInts();
    cudaInts(size_t cnum);
    cudaInts(int a, int b);
    cudaInts(size_t cnum, const int *cpu_data);
    cudaInts(int a, int b, const int *cpu_data);
    cudaInts(const std::vector<int> &cpu_data);
    ~cudaInts();

    size_t size() const;
    int* data();
    const int* data() const;
    void resize(size_t cnum);
    void clear();

    void load(size_t cnum, const int *cpu_data);
    void load(const std::vector<int> &cpu_data);
    void unload(int *cpu_data) const;
    void unload(std::vector<int> &cpu_data) const;

private:
    size_t num;
    int *gpu_data;
};

class cudaDoubles{
public:
    cudaDoubles();
    cudaDoubles(size_t cnum);
    cudaDoubles(int a, int b);
    cudaDoubles(size_t cnum, const double *cpu_data);
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
    void unload(std::vector<double> &cpu_data) const;

private:
    size_t num;
    double *gpu_data;
};


class LinearAlgebraEngineGPU{
public:
    LinearAlgebraEngineGPU();
    ~LinearAlgebraEngineGPU();

    void reset();

    void cublasDGEMM(int M, int N, int K, double alpha, const cudaDoubles &A, const cudaDoubles &B, double beta, cudaDoubles &C);
    void cublasDGEMM(int M, int N, int K, double alpha, const cudaDoubles &A, const std::vector<double> &B, double beta, double C[]);
    // dense matrix-matrix (dgemm) or matrix-vector (dgemv for N == 1) product using Nvidai cuBlas

    #ifdef Tasmanian_ENABLE_MAGMA
    void magmaCudaDGEMM(int gpuID, int M, int N, int K, double alpha, const cudaDoubles &A, const cudaDoubles &B, double beta, cudaDoubles &C);
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

class BaseAccelerationData{
public:
    BaseAccelerationData();
    virtual ~BaseAccelerationData();

    virtual bool isCompatible(TypeAcceleration acc) const = 0;

    virtual void resetGPULoadedData() = 0;
};

// wrapper around cublas and cusparse handles, and magma queue
// it also stores the values and data needed for cuda evaluations
// (points, support, etc.)
class AccelerationDataGPUFull : public BaseAccelerationData{
public:
    AccelerationDataGPUFull();
    ~AccelerationDataGPUFull();

    bool isCompatible(TypeAcceleration acc) const;
    // when setting a new acceleration, check if we need to also set a new object or can reuse this one

    double* getGPUValues() const;
    void loadGPUValues(size_t total_entries, const double *cpu_values);
    // either values (global grid) or coefficients (sequence, local polynomial, etc.)

    double* getGPUNodes() const;
    double* getGPUSupport() const;
    void loadGPUNodesSupport(int total_entries, const double *cpu_nodes, const double *cpu_support);
    // nodes and support for local polynomial evaluations (support also encodes order for special cases)

    int* getGPUpntr() const;
    int* getGPUindx() const;
    int* getGPUroots() const;
    void loadGPUHierarchy(int num_points, const int *pntr, const int *indx, int num_roots, const int *roots);
    // collection of trees to be used in CUDA evaluations of local polynomial grids

    void resetGPULoadedData();
    // deletes nodes, values, and hierarchy

    void cublasDGEMM(bool cpu_pointers, int num_outputs, int num_x, int num_points, const double weights[], double *result);
    // dense matrix-matrix (dgemm) or matrix-vector (dgemv for num_x == 1) product using cublas, the matrix is the getGPUValues()

    void cusparseMatmul(bool cpu_pointers, int num_points, int num_outputs, int num_x, const int *spntr, const int *sindx, const double *svals, int num_nz, double *result);
    // sparse matrix times dense matrix, dense matrix is getGPUValues(), sparse matrix can be given on either cpu or gpu

    void cusparseMatvec(int num_points, int num_x, const int *spntr, const int *sindx, const double *svals, int num_nz, double *result);
    // sparse matrix times a dense vector, makes sense only if the matrix already sits on the gpu (hence no cpu_pointers flag)
    // e.g., assumes num_points == 1, the vectors is getGPUValues(), and the sparse matrix was computed on the GPU

    void cusparseMatveci(int num_outputs, int num_points, int num_nz, const int *sindx, const double *svals, double *result);
    // dense matrix times a sparse vector defined by sindx and svals, currently the sparse vector can only be computed on the cpu
    // the dense matrix is getGPUValues()

    void magmaCudaDGEMM(bool cpu_pointers, int gpuID, int num_outputs, int num_x, int num_points, const double weights[], double *result);
    // dense matrix-matrix (dgemm) or matrix-vector (dgemv for num_x == 1) product using magma, the matrix is the getGPUValues()

protected:
    void makeCuBlasHandle();
    void makeCuSparseHandle();

    void initializeMagma(int gpuID);

private:
    double *gpu_values;
    double *gpu_nodes, *gpu_support; // support is half support and encodes the level of the function, i.e., constant, linear, or quadratic
    int *gpu_hpntr, *gpu_hindx, *gpu_roots;

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

// namespace realized in tsgCudaKernels.cu, each function corresponds to a CUDA kernel for evaluations of basis matrix, domain transform, or fallback linear algebra
namespace TasCUDA{
    // convert transformed points to the canonical domain, all inputs live on the GPU
    void dtrans2can(int dims, int num_x, int pad_size, const double *gpu_trans_a, const double *gpu_trans_b, const double *gpu_x_transformed, double *gpu_x_canonical);

    // evaluate local polynomial rules
    void devalpwpoly(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support, double *gpu_y);
    void devalpwpoly_sparse(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                            int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nzm);
    void devalpwpoly_sparse_dense(int order, TypeOneDRule rule, int dims, int num_x, int num_points, const double *gpu_x, const double *gpu_nodes, const double *gpu_support,
                                 int *gpu_hpntr, int *gpu_hindx, int num_roots, int *gpu_roots, double *gpu_dense);

    // evaluate sequence grids (not done yet)
    //void devalseq(int dims, int num_x, int num_points, int num_nodes, const double *gpu_x, const double *gpu_nodes, const double *gpu_coeff, const int *points, double *gpu_dense);

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

// generic error checking function and types I/O
namespace AccelerationMeta{
    TypeAcceleration getIOAccelerationString(const char * name);
    const char* getIOAccelerationString(TypeAcceleration accel);
    int getIOAccelerationInt(TypeAcceleration accel);
    bool isAccTypeFullMemoryGPU(TypeAcceleration accel);
    bool isAccTypeGPU(TypeAcceleration accel);

    TypeAcceleration getAvailableFallback(TypeAcceleration accel);

    void cudaCheckError(void *cudaStatus, const char *info);
    void cublasCheckError(void *cublasStatus, const char *info);
    void cusparseCheckError(void *cusparseStatus, const char *info);
}

// stores domain transforms, to be used by the top class TasmanianSparseGrid
class AccelerationDomainTransform{
public:
    AccelerationDomainTransform(int num_dimensions, const double *transform_a, const double *transform_b);
    ~AccelerationDomainTransform();

    double* getCanonicalPoints(int num_dimensions, int num_x, const double *gpu_transformed_x);

private:
    // these actually store the rate and shift and not the hard upper/lower limits
    #ifdef Tasmanian_ENABLE_CUDA
    double *gpu_trans_a, *gpu_trans_b;
    int padded_size;
    #endif // Tasmanian_ENABLE_CUDA
};

}

#endif // __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_HPP
