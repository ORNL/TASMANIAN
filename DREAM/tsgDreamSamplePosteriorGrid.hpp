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

#ifndef __TASMANIAN_DREAM_SAMPLE_POSTERIOR_GRID_HPP
#define __TASMANIAN_DREAM_SAMPLE_POSTERIOR_GRID_HPP

#include "tsgDreamSampleGrid.hpp"
#include "tsgDreamSamplePosterior.hpp"

//! \file tsgDreamSamplePosteriorGrid.hpp
//! \brief Templates to sample from Bayesian posterior with sparse grid model.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAM
//!
//! Defines the DREAM template for sampling from a posterior distribution with sparse grid model.

/*!
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMGridModel Level 4 DREAM templates, sampling a posterior from a sparse grids model
 *
 * Level 4 templates assume sampling from a posterior distribution composed of
 * a \b TasDREAM::TasmanianLikelihood, a \b TasGrid::TasmanianSparseGrid approximation to the model,
 * and a \b prior distribution. The likelihood is either one of the Gaussian variations
 * provided by the user, or a use defined class that inherits from the \b TasDREAM::TasmanianLikelihood.
 */

namespace TasDREAM{

/*!
 * \internal
 * \brief Macro to create a model lambda from sparse grid.
 * \ingroup DREAMAux
 * \endinternal
 */
inline std::function<void(std::vector<double> const &candidates, std::vector<double> &values)> makeGridModel(TasGrid::TasmanianSparseGrid const &grid){
    return [&](const std::vector<double> &candidates, std::vector<double> &values)->void{ grid.evaluateBatch(candidates, values); };
}

//! \brief Overloads of \b SampleDREAMPost() which uses a sparse grid as model.
//! \ingroup DREAMGridModel

//! The \b model variable is replaced by the \b grid, see the other overloads for the rest of the variables.
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::function<bool(const std::vector<double> &x)> inside,
                     std::function<void(std::vector<double> &x)> independent_update,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    checkGridSTate(grid, state);
    SampleDREAMPost<form>(num_burnup, num_collect, likelihood, makeGridModel(grid), prior, inside, independent_update, state, differential_update, get_random01);
}


//! \brief Overloads of \b SampleDREAMPost() which uses a sparse grid as model, domain is a hypercube.
//! \ingroup DREAMGridModel

//! The \b model variable is replaced by the \b grid, see the other overloads for the rest of the variables.
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     const std::vector<double> &lower, const std::vector<double> &upper,
                     std::function<void(std::vector<double> &x)> independent_update,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    checkGridSTate(grid, state);
    SampleDREAMPost<form>(num_burnup, num_collect, likelihood, makeGridModel(grid), prior, lower, upper, independent_update, state, differential_update, get_random01);
}


//! \brief Overloads of \b SampleDREAMPost() which uses a sparse grid as model with domain matching the grid.
//! \ingroup DREAMGridModel

//! The \b model variable is replaced by the \b grid, see the other overloads for the rest of the variables.
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::function<void(std::vector<double> &x)> independent_update,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    TasGrid::TypeOneDRule rule = grid.getRule();
    auto domain = extractDomain(grid);
    if ((rule == TasGrid::rule_gausshermite) || (rule == TasGrid::rule_gausshermiteodd)){ // unbounded domain
        SampleDREAMPost<form>(num_burnup, num_collect, likelihood, makeGridModel(grid), prior, domainGaussHermite, independent_update, state, differential_update, get_random01);
    }else if ((rule == TasGrid::rule_gausslaguerre) || (rule == TasGrid::rule_gausslaguerreodd)){ // bounded from below
        SampleDREAMPost<form>(num_burnup, num_collect, likelihood,
                              makeGridModel(grid), prior, makeDomainGaussLaguerre(domain.first), independent_update, state, differential_update, get_random01);
    }else{
        SampleDREAMPost<form>(num_burnup, num_collect, likelihood,
                              makeGridModel(grid), prior, domain.first, domain.second, independent_update, state, differential_update, get_random01);
    }
}


//! \brief Overloads of \b SampleDREAMPost() which uses a sparse grid as model with independent update from known distribution.
//! \ingroup DREAMGridModel

//! The \b model variable is replaced by the \b grid, see the other overloads for the rest of the variables.
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::function<bool(const std::vector<double> &x)> inside,
                     TypeDistribution independent_dist, double independent_magnitude,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    SampleDREAMPost<form>(num_burnup, num_collect, likelihood,
                          makeGridModel(grid), prior, inside, independent_dist, independent_magnitude, state, differential_update, get_random01);
}


//! \brief Overloads of \b SampleDREAMPost() which uses a sparse grid as model with independent update from known distribution and a hypercube.
//! \ingroup DREAMGridModel

//! The \b model variable is replaced by the \b grid, see the other overloads for the rest of the variables.
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     const std::vector<double> &lower, const std::vector<double> &upper,
                     TypeDistribution independent_dist, double independent_magnitude,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    SampleDREAMPost<form>(num_burnup, num_collect, likelihood,
                          makeGridModel(grid), prior, lower, upper, independent_dist, independent_magnitude, state, differential_update, get_random01);
}


//! \brief Overload of \b SampleDREAMGrid() assuming independent update from a list of internals and domain matching the grid.
//! \ingroup DREAMGridModel

//! See other overloads for details.
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     TypeDistribution independent_dist, double independent_magnitude,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    TasGrid::TypeOneDRule rule = grid.getRule();
    auto domain = extractDomain(grid);
    if ((rule == TasGrid::rule_gausshermite) || (rule == TasGrid::rule_gausshermiteodd)){ // unbounded domain
        SampleDREAMPost<form>(num_burnup, num_collect, likelihood, makeGridModel(grid), prior, domainGaussHermite,
                              independent_dist, independent_magnitude, state, differential_update, get_random01);
    }else if ((rule == TasGrid::rule_gausslaguerre) || (rule == TasGrid::rule_gausslaguerreodd)){ // bounded from below
        SampleDREAMPost<form>(num_burnup, num_collect, likelihood, makeGridModel(grid), prior, makeDomainGaussLaguerre(domain.first),
                              independent_dist, independent_magnitude, state, differential_update, get_random01);
    }else{
        SampleDREAMPost<form>(num_burnup, num_collect, likelihood, makeGridModel(grid), prior, domain.first, domain.second,
                              independent_dist, independent_magnitude, state, differential_update, get_random01);
    }
}

}

#endif
