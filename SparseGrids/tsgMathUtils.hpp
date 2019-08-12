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

#ifndef __TASMANIAN_SPARSE_GRID_MATHUTILS_HPP
#define __TASMANIAN_SPARSE_GRID_MATHUTILS_HPP

/*!
 * \internal
 * \file tsgMathUtils.hpp
 * \brief Math functions and constants.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianMaths
 *
 * Simple functions and constnats used throughout the internal algorithms of Tasmanian.
 * \endinternal
 */

/*!
 * \internal
 * \ingroup TasmanianSG
 * \addtogroup TasmanianMaths Math functions and constants
 *
 * \endinternal
 */

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianMaths
 * \brief Math functions and constants.
 */
namespace Maths{
/*!
 * \ingroup TasmanianMaths
 * \brief Computes std::floor(std::log2(i)), but uses only integer operations.
 */
inline int intlog2(int i){
    int result = 0;
    while (i >>= 1){ result++; }
    return result;
}
/*!
 * \ingroup TasmanianMaths
 * \brief Computes std::pow(2, std::floor(std::log2(i))), but uses only integer operations.
 */
inline int int2log2(int i){ // this is effectively: 2^(floor(log_2(i))), when i == 0, this returns 1
    int result = 1;
    while (i >>= 1){ result <<= 1; }
    return result;
}
/*!
 * \ingroup TasmanianMaths
 * \brief Computes std::pow(3, std::floor(std::log(i) / std::log(2))), but uses only integer operations.
 */
inline int int3log3(int i){
    int result = 1;
    while(i >= 1){ i /= 3; result *= 3; }
    return result;
}
/*!
 * \ingroup TasmanianMaths
 * \brief Computes std::pow(2, p), but uses only integer operations.
 */
inline int pow2(int p){ return (1 << p); }
/*!
 * \ingroup TasmanianMaths
 * \brief Computes std::pow(3, p), but uses only integer operations.
 */
inline int pow3(int p){
    int result = 1;
    for(int i=0; i<p; i++) result *= 3;
    return result;
}

/*!
 * \ingroup TasmanianMaths
 * \brief Returns +1 or -1 with the same sign as \b p.
 */
inline double sign(double p){ return (p >= 0.0) ? 1.0 : -1.0; }

/*!
 * \internal
 * \ingroup TasmanianMaths
 * \brief Half-period of the \b std::sin() and \b std::cos() functions.
 *
 * Borrowed from "math.h", let's hope they don't change the value in a future standard.
 * \endinternal
 */
constexpr double pi = 3.14159265358979323846;

/*!
 * \internal
 * \ingroup TasmanianMaths
 * \brief Numerical tolerance for various algorithms.
 *
 *  num_tol is used in many places:
 *  - as a stopping criteria for various iterative schemes (e.g., finding leja points)
 *  - drop criteria for eigenvalue solver related to Gauss rules
 *  - comparison between nodes to detect repeated points in non-nested rules (e.g., all odd Chebyshev rules include zero)
 *  - determining sparse matrix pattern, entries smaller than \b num_tol will be ignored (for wavelet grids)
 *  - drop criteria in estimating anisotropic coefficients (refinement or just the coefficients) surpluses or Legendre coefficients below 10^3 times \b num_tol will be ignored
 * \endinternal
 */
constexpr double num_tol = 1.E-12;

} // namespace Maths

} // namespace TasGrid

#endif
