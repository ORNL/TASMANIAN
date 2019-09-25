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

#ifndef __TASMANIAN_DREAM_HPP
#define __TASMANIAN_DREAM_HPP

#include "tsgDreamSample.hpp"
#include "tsgDreamLikelyGaussian.hpp"

/*!
 * \internal
 * \file TasmanianDREAM.hpp
 * \brief DiffeRential Evolution Adaptive Metropolis methods.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAM
 *
 * The main header required to gain access to the DREAM capabilities of Tasmanian.
 * The header will include all files needed by the DREAM module including
 * the TasmanianSparseGrid.hpp header.
 * \endinternal
 */

/*!
 * \defgroup TasmanianDREAM DREAM: DiffeRential Evolution Adaptive Metropolis
 *
 * \par DREAM
 * The DiffeRential Evolution Adaptive Metropolis is a method to draw samples
 * from an arbitrary probability distribution defined by an arbitrary non-negative function
 * (not necessarily normalized to integrate to 1).
 * In the Tasmanian DREAM module, the samples (and the history) are stored in
 * a TasDREAM::TasmanianDREAM state object which also defines the number of samples
 * and the number of dimensions of the input space.
 * The sampling is performed by the TasDREAM::SampleDREAM() template
 * that takes an initialized state and several callable objects that describe
 * the geometry of the domain and the parameters of the sampling, e.g., number of iterations.
 *
 * \par Bayesian Inference
 * One of the most common applications for DREAM is in the context of Bayesian inference,
 * where the probability distribution is comprised of a model, likelihood and prior.
 * The three components can be combined together with the TasDREAM::posterior()
 * template, which returns a callable object that represents the probability distribution.
 *
 * \par Examples
 * See the included examples.
 */

/*!
 * \ingroup TasmanianDREAM
 * \brief Encapsulates the Tasmanian DREAM module.
 *
 * DREAM related classes and methods sit under the TasDREAM namespace.
 */
namespace TasDREAM{}

#endif
