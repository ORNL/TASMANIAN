/*
 * Copyright (c) 2022, Miroslav Stoyanov & Weiwei Kong
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
 * conditions are met:
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
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND
 * IMPLIED. THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
 * THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL
 * ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE. THE USER ASSUMES
 * RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING
 * FROM OR ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TASMANIAN_OPTIMIZATION_HPP
#define __TASMANIAN_OPTIMIZATION_HPP

#include "tsgParticleSwarm.hpp"
#include "tsgGradientDescent.hpp"

/*!
 * \internal
 * \file TasmanianOptimization.hpp
 * \brief Optimization states and methods.
 * \author Weiwei Kong & Miroslav Stoyanov
 * \ingroup TasmanianOptimization
 *
 * The main header required to gain access to the Optimization capabilities of Tasmanian.
 * The header will include all files needed by the Optimization module including the TasmanianSparseGrid.hpp and
 * TasmanianDREAM headers.
 * \endinternal
 */

/*!
 * \defgroup TasmanianOptimization Optimization
 *
 * \par Optimization
 * A collection of optimization algorithms for minimizing multivariate real-valued functions.
 * The algorithms can be applied to both surrogates constructed with sparse grids and
 * user-provided lambdas.
 */

/*!
 * \ingroup TasmanianOptimization
 * \addtogroup OptimizationState Optimization States
 *
 * The Tasmanian framework uses \b states objects to encapsulate meta-data related to the optimization algorithms.
 * Each algorithm is associated with a separate state class containing specific parameters.
 */

/*!
 * \ingroup TasmanianOptimization
 * \addtogroup OptimizationAlgorithm Optimization Algorithms
 *
 * The optimization algorithms are written in functional programming applied to combinations
 * of objective functionals and optimization states.
 * The state and the algorithm are split so that different functionals can be used
 * with a single state in a multi-fidelity paradigm.
 * An example would be the use of a sparse grid surrogate for the first few steps
 * of the process and switching to the full-model for the last few iterations.
 */

/*!
 * \ingroup TasmanianOptimization
 * \brief Encapsulates the Tasmanian Optimization module.
 */
namespace TasOptimization {}

#endif
