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

#ifndef __TASMANIAN_BLAS_NULL_HPP
#define __TASMANIAN_BLAS_NULL_HPP

/*!
 * \file tsgBlasNull.hpp
 * \brief Wrappers to no-op functions with BLAS signature.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianTPLWrappers
 *
 * Removed the need for preprocessor directives by providing methods
 * with the correct signature but are  no-op.
 */

namespace TasBLAS{

    inline double norm2(int, double const[], int){ return 0.0; }
    inline void vswap(int, double[], int, double[], int){}
    inline void scal(int, double, double[], int){}

    inline void gemv(char, int, int, double, double const[], int, double const[], int, double, double[], int){}
    inline void trsv(char, char, char, int, double const[], int, double[], int){}

    template<typename T>
    inline void trsm(char, char, char, char, int, int, T, T const[], int, T[], int){}

    inline void getrf(int, int, double[], int, int[]){}
    inline void getrs(char, int, int, double const[], int, int const[], double[], int){}

    inline void sterf(int, double[], double[]){}

    template<typename T>
    inline void denseMultiply(int, int, int, T, const T[], const T[], T, T[]){}

    template<typename scalar_type>
    inline void solveLS(char, int, int, scalar_type[], scalar_type[], int = 1){}

    template<typename scalar_type>
    inline void factorizeLQ(int, int, scalar_type[], std::vector<scalar_type> &){}

    template<typename scalar_type>
    inline void multiplyQ(int, int, int, scalar_type const[], std::vector<scalar_type> const&, scalar_type[]){}

    template<typename scalar_type>
    void solveLSmulti(int, int, scalar_type[], int, scalar_type[]){}
}

#endif
