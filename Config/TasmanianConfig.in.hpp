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

#ifndef __TASMANIAN_CONFIG_HPP
#define __TASMANIAN_CONFIG_HPP

#define TASMANIAN_VERSION_MAJOR @Tasmanian_VERSION_MAJOR@
#define TASMANIAN_VERSION_MINOR @Tasmanian_VERSION_MINOR@
#define TASMANIAN_VERSION_STRING "@Tasmanian_VERSION_MAJOR@.@Tasmanian_VERSION_MINOR@@Tasmanian_version_comment@"
#define TASMANIAN_LICENSE "@Tasmanian_license@"

#define TASMANIAN_GIT_COMMIT_HASH "@Tasmanian_git_hash@"
#define TASMANIAN_CXX_FLAGS "@Tasmanian_cxx_flags@"

// cmake options propagated to the source code
#cmakedefine Tasmanian_ENABLE_BLAS
#cmakedefine Tasmanian_ENABLE_CUDA
#cmakedefine Tasmanian_ENABLE_HIP
#cmakedefine Tasmanian_ENABLE_DPCPP
#cmakedefine Tasmanian_ENABLE_MAGMA
#cmakedefine Tasmanian_ENABLE_MPI

#if defined(Tasmanian_ENABLE_CUDA) || defined(Tasmanian_ENABLE_HIP) || defined(Tasmanian_ENABLE_DPCPP)
#define Tasmanian_ENABLE_GPU // One GPU definition
#endif

// handle variations in library standards
#cmakedefine Tasmanian_BLAS_HAS_ZGELQ

#endif
