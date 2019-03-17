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

#ifndef __TASMANIAN_DREAM_SAMPLE_GRID_HPP
#define __TASMANIAN_DREAM_SAMPLE_GRID_HPP

#include "tsgDreamSample.hpp"

//! \file tsgDreamSampleGrid.hpp
//! \brief DREAM methods using likelihood approximated by a sparse grid.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAM
//!
//! Defines the templates and overloads for sampling a posterior distribution using a likelihood defined by a sparse grid approximation.

/*!
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMSampleGrid Level 2 DREAM templates, sampling using a sparse grids likelihood
 *
 * Level 2 templates assume sampling from a posterior distribution (i.e., Bayesian inference),
 * where the likelihood is approximated using a Sparse Grid.
 * Note that using \b uniform_prior() effectively samples from a sparse grids approximated
 * probability distribution, and hence the templates have a whider potential usage.
 */

namespace TasDREAM{

/*!
 * \internal
 * \brief Macro creating a \b lambda for the probability distribution combining a sparse grid and a prior.
 * \ingroup DREAMAux
 *
 * The same \b lambda function is used for multiple overloads, so long as the variable names match, this will work.
 * The variable names must match, in order to have consistency anyway. Specifically,
 * - \b grid is the TasmanianSparseGrid variable
 * - \b prior is the \b lambda for the prior
 * - \b form is the template parameter that chooses between regular and logarithmic form
 * \endinternal
 */
template<TypeSamplingForm form> std::function<void(const std::vector<double> &candidates, std::vector<double> &values)>
makePDFGridPrior(TasGrid::TasmanianSparseGrid const &grid, std::function<void(std::vector<double> const &candidates, std::vector<double> &values)> prior){
    return [&](std::vector<double> const &candidates, std::vector<double> &values) -> void{
        grid.evaluateBatch(candidates, values);

        std::vector<double> priors(values.size());
        prior(candidates, priors);

        auto iv = values.begin();
        if (form == regform){
            for(auto p : priors) *iv++ *= p;
        }else{
            for(auto p : priors) *iv++ += p;
        }
    };
}

/*!
 * \internal
 * \brief Macro that checks if \b grid and \b state have the same number of dimensions, throws runtime_error otherwise.
 * \ingroup DREAMAux
 * \endinternal
 */
inline void checkGridSTate(TasGrid::TasmanianSparseGrid const &grid, TasmanianDREAM const &state){
    if (grid.getNumDimensions() != state.getNumDimensions()) throw std::runtime_error("ERROR: mismatch between the dimensions of the grid and the DREAM state.");
}

//! \internal
//! \brief Extract the \b grid rule and domain.
//! \ingroup DREAMAux

//! Extract the \b grid rule and domain.
#define __TASDREAM_GRID_EXTRACT_RULE \
    TasGrid::TypeOneDRule rule = grid.getRule(); \
    std::vector<double> transform_a, transform_b; \
    if (grid.isSetDomainTransfrom()) grid.getDomainTransform(transform_a, transform_b); \

//! \internal
//! \brief Get the Gauss-Hermite lambda.
//! \ingroup DREAMAux

//! Get the Gauss-Hermite lambda.
#define __TASDREAM_GRID_DOMAIN_GHLAMBDA [&](const std::vector<double> &)->bool{ return true; }

//! \internal
//! \brief Get the Gauss-Laguerre lambda.
//! \ingroup DREAMAux

//! Get the Gauss-Laguerre lambda.
#define __TASDREAM_GRID_DOMAIN_GLLAMBDA [&](const std::vector<double> &x)->bool{ \
    auto ix = x.begin(); \
    for(auto bound : transform_a) if (*ix++ < bound) return false; \
    return true; }

//! \internal
//! \brief Get the default transform weights.
//! \ingroup DREAMAux

//! Get the default transform weights.
#define __TASDREAM_GRID_DOMAIN_DEFAULTS \
    if (!grid.isSetDomainTransfrom()){\
        transform_a.resize(grid.getNumDimensions(), ((rule == TasGrid::rule_fourier) ? 0.0 : -1.0));\
        transform_b.resize(grid.getNumDimensions(), 1.0);\
    }


//! \brief Variation of \b SampleDREAM() which assumes a Bayesian inference problem with likelihood approximated by a sparse grid.
//! \ingroup DREAMSampleGrid

//! The inputs are identical to the ones defined in \b SampleDREAM(), with \b probability_distribution() replaced by the combination of sparse grid and \b prior.
//! The values obtained from the grid are multiplied (or incremented in \b logform) with the values from the prior.
//!
//! Note that \b SampleDREAMGrid() with \b uniform_prior() can also be used to generate samples from a probability distribution that has been
//! approximated using a sparse grid.
template<TypeSamplingForm form = regform>
void SampleDREAMGrid(int num_burnup, int num_collect,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::function<bool(const std::vector<double> &x)> inside,
                     std::function<void(std::vector<double> &x)> independent_update,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    checkGridSTate(grid, state);
    SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), inside, independent_update, state, differential_update, get_random01);
}


//! \brief Overload of \b SampleDREAMGrid() assuming the domain is a hyperbube with min/max values given by \b lower and \b upper.
//! \ingroup DREAMSampleGrid

//! See the overloads of \b SampleDREAM() regarding the vectors \b lower and \b upper.
template<TypeSamplingForm form = regform>
void SampleDREAMGrid(int num_burnup, int num_collect,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     const std::vector<double> &lower, const std::vector<double> &upper,
                     std::function<void(std::vector<double> &x)> independent_update,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    checkGridSTate(grid, state);
    SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), lower, upper, independent_update, state, differential_update, get_random01);
}


//! \brief Overload of \b SampleDREAMGrid() assuming the domain matches the sparse grid domain.
//! \ingroup DREAMSampleGrid

//! Grid constructed using a Gauss-Hermite rule will use the entire real plane as domain, the \b inside() function is always true.
//! The Gauss-Laguerre rule is bounded from below, but not above.
//! The canonical (no domain transform) Fourier grids are supported on hypercube [0,1]
//! In all other cases, the domain is set to the canonical hypercube [-1, 1] or the hypercube specified by the domain transform loaded in the grid.
//! See \b TasGrid::TypeOneDRule for details regarding the rules and their support.
template<TypeSamplingForm form = regform>
void SampleDREAMGrid(int num_burnup, int num_collect,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::function<void(std::vector<double> &x)> independent_update,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    __TASDREAM_GRID_EXTRACT_RULE
    if ((rule == TasGrid::rule_gausshermite) || (rule == TasGrid::rule_gausshermiteodd)){ // unbounded domain
        SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), __TASDREAM_GRID_DOMAIN_GHLAMBDA, independent_update, state, differential_update, get_random01);
    }else if ((rule == TasGrid::rule_gausslaguerre) || (rule == TasGrid::rule_gausslaguerreodd)){ // bounded from below
        SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), __TASDREAM_GRID_DOMAIN_GHLAMBDA, independent_update, state, differential_update, get_random01);
    }else{
        __TASDREAM_GRID_DOMAIN_DEFAULTS
        SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), transform_a, transform_b, independent_update, state, differential_update, get_random01);
    }
}


//! \brief Overload of \b SampleDREAMGrid() assuming independent update from a list of internals.
//! \ingroup DREAMSampleGrid

//! See the overloads of \b SampleDREAM() regarding the known updates.
template<TypeSamplingForm form = regform>
void SampleDREAMGrid(int num_burnup, int num_collect,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::function<bool(const std::vector<double> &x)> inside,
                     TypeDistribution independent_dist, double independent_magnitude,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), inside, independent_dist, independent_magnitude, state, differential_update, get_random01);
}


//! \brief Overload of \b SampleDREAMGrid() assuming hypercube domain and independent update from a list of internals.
//! \ingroup DREAMSampleGrid

//! See the overloads of \b SampleDREAM() regarding the known updates and the \b lower and \b upper vectors.
template<TypeSamplingForm form = regform>
void SampleDREAMGrid(int num_burnup, int num_collect,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     const std::vector<double> &lower, const std::vector<double> &upper,
                     TypeDistribution independent_dist, double independent_magnitude,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), lower, upper, independent_dist, independent_magnitude, state, differential_update, get_random01);
}


//! \brief Overload of \b SampleDREAMGrid() assuming independent update from a list of internals and domain matching the grid.
//! \ingroup DREAMSampleGrid

//! See other overloads for details.
template<TypeSamplingForm form = regform>
void SampleDREAMGrid(int num_burnup, int num_collect,
                     const TasGrid::TasmanianSparseGrid &grid,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     TypeDistribution independent_dist, double independent_magnitude,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    __TASDREAM_GRID_EXTRACT_RULE
    if ((rule == TasGrid::rule_gausshermite) || (rule == TasGrid::rule_gausshermiteodd)){ // unbounded domain
        SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), __TASDREAM_GRID_DOMAIN_GHLAMBDA,
                          independent_dist, independent_magnitude, state, differential_update, get_random01);
    }else if ((rule == TasGrid::rule_gausslaguerre) || (rule == TasGrid::rule_gausslaguerreodd)){ // bounded from below
        SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), __TASDREAM_GRID_DOMAIN_GHLAMBDA,
                          independent_dist, independent_magnitude, state, differential_update, get_random01);
    }else{
        __TASDREAM_GRID_DOMAIN_DEFAULTS
        SampleDREAM<form>(num_burnup, num_collect, makePDFGridPrior<form>(grid, prior), transform_a, transform_b, independent_dist, independent_magnitude, state, differential_update, get_random01);
    }
}

}

#endif
