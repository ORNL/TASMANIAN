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

#ifndef __TASMANIAN_CUDA_WRAPPERS_CPP
#define __TASMANIAN_CUDA_WRAPPERS_CPP

#include "tsgGpuWrappers.hpp"

/*!
 * \file tsgGpuNull.cpp
 * \brief Wrappers to no-op functions with GPU signature.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * Reduces the need for preprocessor macros by providing methods with GPU signature but no-op implementation.
 */

namespace TasGrid{
/*
 * Meta methods
 */
template<typename T> void GpuVector<T>::resize(AccelerationContext const*, size_t){}
template<typename T> void GpuVector<T>::clear(){}
template<typename T> void GpuVector<T>::load(AccelerationContext const*, size_t, const T[]){}
template<typename T> void GpuVector<T>::unload(AccelerationContext const*, size_t, T[]) const{}

template void GpuVector<double>::resize(AccelerationContext const*, size_t);
template void GpuVector<double>::clear();
template void GpuVector<double>::load(AccelerationContext const*, size_t, const double*);
template void GpuVector<double>::unload(AccelerationContext const*, size_t, double*) const;

template void GpuVector<std::complex<double>>::resize(AccelerationContext const*, size_t);
template void GpuVector<std::complex<double>>::clear();
template void GpuVector<std::complex<double>>::load(AccelerationContext const*, size_t, const std::complex<double>*);
template void GpuVector<std::complex<double>>::unload(AccelerationContext const*, size_t, std::complex<double>*) const;

template void GpuVector<float>::resize(AccelerationContext const*, size_t);
template void GpuVector<float>::clear();
template void GpuVector<float>::load(AccelerationContext const*, size_t, const float*);
template void GpuVector<float>::unload(AccelerationContext const*, size_t, float*) const;

template void GpuVector<int>::resize(AccelerationContext const*, size_t);
template void GpuVector<int>::clear();
template void GpuVector<int>::load(AccelerationContext const*, size_t, const int*);
template void GpuVector<int>::unload(AccelerationContext const*, size_t, int*) const;

int AccelerationMeta::getNumGpuDevices(){ return 0; }
void AccelerationMeta::setDefaultGpuDevice(int){}
unsigned long long AccelerationMeta::getTotalGPUMemory(int){ return 0; }
std::string AccelerationMeta::getGpuDeviceName(int){ return std::string(); }
template<typename T> void AccelerationMeta::recvGpuArray(AccelerationContext const*, size_t, const T[], std::vector<T>&){}
template<typename T> void AccelerationMeta::delGpuArray(AccelerationContext const*, T*){}

void* AccelerationMeta::createCublasHandle(){ return nullptr; }
void AccelerationMeta::deleteCublasHandle(void*){}

template void AccelerationMeta::recvGpuArray<double>(AccelerationContext const*, size_t num_entries, const double*, std::vector<double>&);
template void AccelerationMeta::recvGpuArray<float>(AccelerationContext const*, size_t num_entries, const float*, std::vector<float>&);
template void AccelerationMeta::recvGpuArray<int>(AccelerationContext const*, size_t num_entries, const int*, std::vector<int>&);

template void AccelerationMeta::delGpuArray<double>(AccelerationContext const*, double*);
template void AccelerationMeta::delGpuArray<float>(AccelerationContext const*, float*);
template void AccelerationMeta::delGpuArray<int>(AccelerationContext const*, int*);

namespace TasGpu{
/*
 * Algorithm section
 */
void factorizePLU(AccelerationContext const*, int, double[], int[]){}
void solvePLU(AccelerationContext const*, char, int, double const[], int const[], double[]){}
void solvePLU(AccelerationContext const*, char, int, double const[], int const[], int, double[]){}

template<typename scalar_type>
void solveLSmultiGPU(AccelerationContext const*, int, int, scalar_type[], int, scalar_type[]){}

template void solveLSmultiGPU<double>(AccelerationContext const*, int, int, double[], int, double[]);
template void solveLSmultiGPU<std::complex<double>>(AccelerationContext const*, int, int, std::complex<double>[], int, std::complex<double>[]);

template<typename scalar_type>
void solveLSmultiOOC(AccelerationContext const*, int, int, scalar_type[], int, scalar_type[]){}

template void solveLSmultiOOC<double>(AccelerationContext const*, int, int, double[], int, double[]);
template void solveLSmultiOOC<std::complex<double>>(AccelerationContext const*, int, int, std::complex<double>[], int, std::complex<double>[]);

template<typename scalar_type>
void denseMultiply(AccelerationContext const*, int, int, int, typename GpuVector<scalar_type>::value_type, GpuVector<scalar_type> const&,
                   GpuVector<scalar_type> const&, typename GpuVector<scalar_type>::value_type, scalar_type[]){}

template void denseMultiply<float>(AccelerationContext const*, int, int, int, float,
                                   GpuVector<float> const&, GpuVector<float> const&, float, float[]);
template void denseMultiply<double>(AccelerationContext const*, int, int, int, double,
                                    GpuVector<double> const&, GpuVector<double> const&, double, double[]);

template<typename scalar_type>
void sparseMultiply(AccelerationContext const*, int, int, int, typename GpuVector<scalar_type>::value_type,
                    GpuVector<scalar_type> const&, GpuVector<int> const&, GpuVector<int> const&,
                    GpuVector<scalar_type> const&, scalar_type[]){}

template void sparseMultiply<float>(AccelerationContext const*, int, int, int, float, GpuVector<float> const &A,
                                    GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<float> const &vals, float C[]);
template void sparseMultiply<double>(AccelerationContext const*, int, int, int, double, GpuVector<double> const &A,
                                     GpuVector<int> const &pntr, GpuVector<int> const &indx, GpuVector<double> const &vals, double C[]);

template<typename T> void dtrans2can(AccelerationContext const*, bool, int, int, int, double const[], double const[], T const[], T[]){}
template<typename T> void devalpwpoly(AccelerationContext const*, int, TypeOneDRule, int, int, int, const T[], const T[], const T[], T[]){}

template<typename T>
void devalpwpoly_sparse(AccelerationContext const*, int, TypeOneDRule, int, int, const T[], GpuVector<T> const&, GpuVector<T> const&,
                                GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int>&, GpuVector<int>&, GpuVector<T>&){}

template<typename T>
void devalseq(AccelerationContext const*, int, int, std::vector<int> const&, T const[], GpuVector<int> const&, GpuVector<int> const&, GpuVector<T> const&, GpuVector<T> const&, T[]){}
template<typename T>
void devalfor(AccelerationContext const*, int, int, std::vector<int> const&, const T[], GpuVector<int> const&, GpuVector<int> const&, T[], typename GpuVector<T>::value_type[]){}

template<typename T>
void devalglo(AccelerationContext const*, bool, bool, int, int, int, int, T const[],
              GpuVector<T> const&, GpuVector<T> const&, GpuVector<T> const&,
              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, T[]){}

void fillDataGPU(AccelerationContext const*, double, long long, long long, double[]){}

template void dtrans2can<float>(AccelerationContext const*, bool, int, int, int, double const*, double const*, float const*, float*);
template void dtrans2can<double>(AccelerationContext const*, bool, int, int, int, double const*, double const*, double const*, double*);

template void devalpwpoly<float>(AccelerationContext const*, int, TypeOneDRule, int, int, int, float const*, float const*, float const*, float*);
template void devalpwpoly<double>(AccelerationContext const*, int, TypeOneDRule, int, int, int, double const*, double const*, double const*, double*);

template void devalpwpoly_sparse<float>(AccelerationContext const*, int, TypeOneDRule, int, int, float const*,
                                        GpuVector<float> const&, GpuVector<float> const&,
                                        GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                        GpuVector<int>&, GpuVector<int>&, GpuVector<float>&);
template void devalpwpoly_sparse<double>(AccelerationContext const*, int, TypeOneDRule, int, int, double const*,
                                         GpuVector<double> const&, GpuVector<double> const&,
                                         GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                                         GpuVector<int>&, GpuVector<int>&, GpuVector<double>&);

template void devalfor<float>(AccelerationContext const*, int, int, std::vector<int> const&, float const*, GpuVector<int> const&, GpuVector<int> const&, float*, float*);
template void devalfor<double>(AccelerationContext const*, int, int, std::vector<int> const&, double const*, GpuVector<int> const&, GpuVector<int> const&, double*, double*);

template void devalseq<float>(AccelerationContext const*, int, int, std::vector<int> const&, float const*, GpuVector<int> const&,
                              GpuVector<int> const&, GpuVector<float> const&, GpuVector<float> const&, float*);
template void devalseq<double>(AccelerationContext const*, int, int, std::vector<int> const&, double const*, GpuVector<int> const&,
                               GpuVector<int> const&, GpuVector<double> const&, GpuVector<double> const&, double*);

template void devalglo<float>(AccelerationContext const*, bool, bool, int, int, int, int,
                              float const*, GpuVector<float> const&, GpuVector<float> const&, GpuVector<float> const&,
                              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                              GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, float*);
template void devalglo<double>(AccelerationContext const*, bool, bool, int, int, int, int,
                               double const*, GpuVector<double> const&, GpuVector<double> const&, GpuVector<double> const&,
                               GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                               GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&,
                               GpuVector<int> const&, GpuVector<int> const&, GpuVector<int> const&, double*);

template<typename T> void load_n(AccelerationContext const*, T const*, size_t, T*){}

template void load_n<int>(AccelerationContext const*, int const*, size_t, int*);
template void load_n<float>(AccelerationContext const*, float const*, size_t, float*);
template void load_n<double>(AccelerationContext const*, double const*, size_t, double*);
template void load_n<std::complex<double>>(AccelerationContext const*, std::complex<double> const*, size_t, std::complex<double>*);

}
}

#endif
