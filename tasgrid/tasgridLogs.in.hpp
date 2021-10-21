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

#include <iostream>
#include <fstream>

#define TasmanianGPTableBuild "@CMAKE_CURRENT_BINARY_DIR@/SparseGrids/GaussPattersonRule.table"
#define TasmanianGPTableInstall "@Tasmanian_final_install_path@/share/Tasmanian/GaussPattersonRule.table"

inline void show_log(){
    std::ifstream logfile("@Tasmanian_final_install_path@/share/Tasmanian/Tasmanian.log");
    std::cout << "\n" << logfile.rdbuf() << std::endl;
}

inline void show_cmake_log(){
    std::cout << "\nCMake parameters:\n--------------------------------------------------------------------------------"
    << R"CMAKEVARS(
CMAKE_COMMAND                      @CMAKE_COMMAND@
CMAKE_BUILD_TYPE                   @CMAKE_BUILD_TYPE@
CMAKE_INSTALL_PREFIX               @CMAKE_INSTALL_PREFIX@
CMAKE_CXX_FLAGS                    @CMAKE_CXX_FLAGS@
CMAKE_CXX_FLAGS_DEBUG              @CMAKE_CXX_FLAGS_DEBUG@
CMAKE_CXX_FLAGS_RELEASE            @CMAKE_CXX_FLAGS_RELEASE@
CMAKE_CXX_COMPILER                 @CMAKE_CXX_COMPILER@
CMAKE_CUDA_COMPILER                @CMAKE_CUDA_COMPILER@
CMAKE_CUDA_FLAGS                   @CMAKE_CUDA_FLAGS@
CMAKE_CUDA_FLAGS_DEBUG             @CMAKE_CUDA_FLAGS_DEBUG@
CMAKE_CUDA_FLAGS_RELEASE           @CMAKE_CUDA_FLAGS_RELEASE@
Tasmanian_cudamathlibs             @Tasmanian_cudamathlibs@
Tasmanian_cudaruntime              @Tasmanian_cudaruntime@
PYTHON_EXECUTABLE                  @PYTHON_EXECUTABLE@
BUILD_SHARED_LIBS                  @BUILD_SHARED_LIBS@
BLAS_LIBRARIES                     @BLAS_LIBRARIES@
CMAKE_THREAD_LIBS_INIT             @CMAKE_THREAD_LIBS_INIT@
OpenMP_CXX_LIBRARIES               @OpenMP_CXX_LIBRARIES@
OpenMP_CXX_FLAGS                   @OpenMP_CXX_FLAGS@
MPI_CXX_LIBRARIES                  @MPI_CXX_LIBRARIES@
MPI_CXX_INCLUDE_PATH               @MPI_CXX_INCLUDE_PATH@
MPI_CXX_COMPILE_FLAGS              @MPI_CXX_COMPILE_FLAGS@
MPI_CXX_LINK_FLAGS                 @MPI_CXX_LINK_FLAGS@
MPIEXEC_EXECUTABLE                 @MPIEXEC_EXECUTABLE@
Tasmanian_ENABLE_RECOMMENDED       @Tasmanian_ENABLE_RECOMMENDED@
Tasmanian_ENABLE_OPENMP            @Tasmanian_ENABLE_OPENMP@
Tasmanian_ENABLE_BLAS              @Tasmanian_ENABLE_BLAS@
Tasmanian_ENABLE_PYTHON            @Tasmanian_ENABLE_PYTHON@
Tasmanian_ENABLE_CUDA              @Tasmanian_ENABLE_CUDA@
Tasmanian_ENABLE_HIP               @Tasmanian_ENABLE_HIP@
Tasmanian_ENABLE_DPCPP             @Tasmanian_ENABLE_DPCPP@
Tasmanian_ENABLE_FORTRAN           @Tasmanian_ENABLE_FORTRAN@
Tasmanian_ENABLE_MPI               @Tasmanian_ENABLE_MPI@
Tasmanian_ENABLE_MAGMA             @Tasmanian_ENABLE_MAGMA@
Tasmanian_MATLAB_WORK_FOLDER       @Tasmanian_MATLAB_WORK_FOLDER@
Tasmanian_ENABLE_DOXYGEN           @Tasmanian_ENABLE_DOXYGEN@
DOXYGEN_INTERNAL_DOCS              @DOXYGEN_INTERNAL_DOCS@
Tasmanian_PYTHONPATH               @Tasmanian_PYTHONPATH@
Tasmanian_final_install_path       @Tasmanian_final_install_path@)CMAKEVARS"
    << "\n--------------------------------------------------------------------------------\n"
    << std::endl;
}
