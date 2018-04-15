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

#ifndef __TASMANIAN_SPARSE_GRID_ENUMERATES_HPP
#define __TASMANIAN_SPARSE_GRID_ENUMERATES_HPP

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <string.h>

#include "TasmanianConfig.hpp"

namespace TasGrid{

using std::cout; // for debugging purposes
using std::endl; // for debugging purposes

using std::cerr; // for error messages

enum TypeIndexRelation{ // internal for IndexSets, unambiguous set comparison
    type_abeforeb, type_bbeforea, type_asameb
};

enum TypeDepth{
    type_none,
    type_level, type_curved,
    type_iptotal, type_ipcurved,
    type_qptotal, type_qpcurved,
    type_hyperbolic, type_iphyperbolic, type_qphyperbolic,
    type_tensor, type_iptensor, type_qptensor
};

enum TypeOneDRule{ // list of one-d rules
    rule_none, // encountering is an indication of error
    rule_clenshawcurtis,
    rule_clenshawcurtis0, // assumes zero boundary
    rule_chebyshev,
    rule_chebyshevodd,
    rule_gausslegendre,
    rule_gausslegendreodd,
    rule_gausspatterson,
    rule_leja,
    rule_lejaodd,
    rule_rleja,
    rule_rlejadouble2,
    rule_rlejadouble4,
    rule_rlejaodd,
    rule_rlejashifted,
    rule_rlejashiftedeven,
    rule_rlejashifteddouble,
    rule_maxlebesgue,
    rule_maxlebesgueodd,
    rule_minlebesgue,
    rule_minlebesgueodd,
    rule_mindelta,
    rule_mindeltaodd,
    rule_gausschebyshev1,
    rule_gausschebyshev1odd,
    rule_gausschebyshev2,
    rule_gausschebyshev2odd,
    rule_fejer2,
    rule_gaussgegenbauer,
    rule_gaussgegenbauerodd,
    rule_gaussjacobi,
    rule_gaussjacobiodd,
    rule_gausslaguerre,
    rule_gausslaguerreodd,
    rule_gausshermite,
    rule_gausshermiteodd,
    rule_customtabulated,
    // Piece-Wise rules
    rule_localp,
    rule_localp0,
    rule_semilocalp,
    // Wavelet rules
    rule_wavelet
};

enum TypeRefinement{
    refine_classic, refine_parents_first, refine_direction_selective, refine_fds, refine_none /* FDS = parents_first + direction_selective */
};

enum TypeAcceleration{
    accel_none,
    accel_cpu_blas, // default if compiled with CPU-BLAS flag
    accel_gpu_fullmemory, // load all values entirely into the GPU
    accel_gpu_default,
    accel_gpu_cublas,
    accel_gpu_cuda,
    accel_gpu_magma
};

enum TypeLocalPolynomialBackendFlavor{
    flavor_auto,
    flavor_sparse_sparse,
    flavor_sparse_dense,
    flavor_dense_sparse,
    flavor_dense_dense,
    flavor_cuda
};



////////////////////////////////////////////////////////////
//                                                        //
//     ---  Moved from tsgHardcodedConstants.hpp  ---     //
//                                                        //
////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////
//
//  On the purpose of this section:
//
//
//  Hardcoding constants is a bad coding practice and should be avoided.
//  On the other hand, numerical algorithms depend on many small tweaking parameters.
//  For example, this code computes the nodes and abscissas of Gauss-Legendre rule on the fly,
//  this is done with an iterative eigenvalue decomposition method that needs a stopping criteria,
//  we can add another user specified parameter to the corresponding "makeGlobalGrid()" rule,
//  however, this is an additional tweaking variable that the user has to understand and control.
//  The goal of our code is to provide a most seamless experience to the user,
//  one should think about the desired properties of a quadrature rule or interpolant and not
//  about convergence of an iterative scheme (especially one hidden from top view)
//  Therefore, the tolerance for such convergence criteria needs to be set inside the code to
//  a "reasonable" value, which is a value that would work well for the overwhelming majority
//  of use cases.
//
//  On the other hand, every application is different and the usage of the code will change
//  over time. It is unreasonable to believe that one single value would work for absolutely
//  everyone. Instead of "hiding" hardcoded constant throughout the code, all such constants
//  will be exposed here so the user can make adjustments in compile time.
//
//  Long story short: do not adjust those variables unless you have a good reason.
//
///////////////////////////////////////////////////////////////////////////////////////////////

//  NUM_TOL is used in many places:
// - as a stopping criteria for various iterative schemes (e.g., finding leja points)
// - drop criteria for eigenvalue solver related to Gauss rules
// - comparison between nodes to detect repeated points in non-nested rules (e.g., all odd Chebyshev rules include zero)
// - determining sparse matrix pattern, entries smaller than NUM_TOL will be ignored (for wavelet grids)
// - drop criteria in estimating anisotropic coefficients (refinement or just the coefficients) surpluses or Legendre coefficients below 10^3 times NUM_TOL will be ignored
#define TSG_NUM_TOL 1.E-12 // to move in hardcoded constants

// this defines the maximum number of secant method iterations to be used for finding Leja, Lebesgue, and Delta points
// this is a safeguard criteria to prevent "hanging" in a loop
#define TSG_MAX_SECANT_ITERATIONS 1000 // to move in hardcoded constants

// this defines the threshold for switching between sparse and dense version of batch evaluate for Local Polynomial Grids
// generally, because of caching and reuse of data, dense operations have 10x more flops per second than sparse ops
// unless the sparse version of the algorithm can cut total flops by 10x (i.e., fill is less than 10%), then it
// is faster to convert the sparse matrix into a dense one. This assumes matrix-matrix product, which happens only
// if the number of outputs is sufficiently large, if the outputs are few, there is very little caching anyway
// hence use the sparse version regardless
#define TSG_LOCALP_BLAS_NUM_OUTPUTS 2048

}

#endif
