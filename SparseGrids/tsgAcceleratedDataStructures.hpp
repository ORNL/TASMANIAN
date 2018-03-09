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

#include "tsgEnumerates.hpp"

namespace TasGrid{

class BaseAccelerationData{
public:
    BaseAccelerationData();
    virtual ~BaseAccelerationData();

    virtual bool isCompatible(TypeAcceleration acc) const = 0;

    virtual void resetValuesAndSurpluses() = 0;
};

class AccelerationDataGPUFull : public BaseAccelerationData{
public:
    AccelerationDataGPUFull();
    ~AccelerationDataGPUFull();

    void setLogStream(std::ostream *os);

    bool isCompatible(TypeAcceleration acc) const;

    double* getGPUValues() const;
    void loadGPUValues(int total_entries, const double *cpu_values);
    void resetValuesAndSurpluses();

    void cublasDGEMV(int num_outputs, int num_points, const double cpu_weights[], double *cpu_result);
    void cublasDGEMM(int num_outputs, int num_points, int num_x, const double cpu_weights[], double *cpu_result);

    void cusparseDCRMM2(int num_points, int num_outputs, int num_x, const int *cpu_pntr, const int *cpu_indx, const double *cpu_vals, double *cpu_result);
    void cusparseDCRSMM(int num_points, int num_outputs, const int *cpu_pntr, const int *cpu_indx, const double *cpu_vals, const double *values, double *surpluses);

protected:
    void makeCuBlasHandle();
    void makeCuSparseHandle();

private:
    double *gpu_values;
    void *cublasHandle;
    void *cusparseHandle;

    std::ostream *logstream;
};

namespace TasCUDA{
    // matrix multiply double-precision (d), level 3, general (ge), column compressed sparse (cs): C = A * B
    // A is N by K, B is K by M, C is N by M
    // A and C are stored in column format, B is sparse column compresses
    void d3gecs(int N, int M, const double *gpuA, const int *cpuBpntr, const int *cpuBindx, const double *cpuBvals, double *cpuC, std::ostream *os = 0);
    // matrix solve double-precision (d), level 3, general (ge), column compressed sparse (cs) (solve): C = A * B
    // A is N by M, B is M by M, C is N by M
    // A and C are stored in column format, B is sparse column compresses and unit triangular (the diagonal entry is no include in the pattern)
    void d3gecss(int N, int M, int *levels, int top_level, const int *cpuBpntr, const int *cpuBindx, const double *cpuBvals, const double *cpuA, double *cpuC, std::ostream *os);
}

namespace AccelerationMeta{
    TypeAcceleration getIOAccelerationString(const char * name);
    const char* getIOAccelerationString(TypeAcceleration accel);
    int getIOAccelerationInt(TypeAcceleration accel);
    bool isAccTypeFullMemoryGPU(TypeAcceleration accel);
    bool isAccTypeGPU(TypeAcceleration accel);

    void cudaCheckError(void *cudaStatus, const char *info, std::ostream *os);
    void cublasCheckError(void *cublasStatus, const char *info, std::ostream *os);
    void cusparseCheckError(void *cusparseStatus, const char *info, std::ostream *os);
}

}

#endif // __TASMANIAN_SPARSE_GRID_ACCELERATED_DATA_STRUCTURES_HPP
