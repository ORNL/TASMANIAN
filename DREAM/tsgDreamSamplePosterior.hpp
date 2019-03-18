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

#ifndef __TASMANIAN_DREAM_SAMPLE_POSTERIOR_HPP
#define __TASMANIAN_DREAM_SAMPLE_POSTERIOR_HPP

#include "tsgDreamSample.hpp"
#include "tsgDreamLikelyGaussian.hpp"

//! \file tsgDreamSamplePosterior.hpp
//! \brief Templates to sample from Bayesian posterior.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAM
//!
//! Defines the DREAM template for sampling from a posterior distribution.

/*!
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMSampleModel Level 3 DREAM templates, sampling a posterior from a user provided model
 *
 * Level 3 templates assume sampling from a posterior distribution composed of
 * a \b TasDREAM::TasmanianLikelihood, a \b model defined by a user provided \b lambda,
 * and a \b prior distribution. The likelihood is either one of the Gaussian variations
 * provided by the user, or a use defined class that inherits from the \b TasDREAM::TasmanianLikelihood.
 */

namespace TasDREAM{

/*!
 * \internal
 * \brief Create a probability distribution lambda from the model, likelihood and prior.
 * \ingroup DREAMAux
 *
 * The same \b lambda function is used for multiple overloads.
 * - \b likelihood is the TasmanianLikelihood
 * - \b model is the \b lambda for the model
 * - \b prior is the \b lambda for the prior
 * - \b form is the template parameter that chooses between regular and logarithmic form
 * \endinternal
 */
template<TypeSamplingForm form> std::function<void(const std::vector<double> &candidates, std::vector<double> &values)>
makePDFPosterior(TasmanianLikelihood const &likelihood,
                 std::function<void(std::vector<double> const &candidates, std::vector<double> &values)> model,
                 std::function<void(std::vector<double> const &candidates, std::vector<double> &values)> prior){
    return [&](const std::vector<double> &candidates, std::vector<double> &values)->void{
        std::vector<double> model_outs(likelihood.getNumOuputs() * values.size());
        model(candidates, model_outs);
        likelihood.getLikelihood(form, model_outs, values);

        std::vector<double> prior_vals(values.size());
        prior(candidates, prior_vals);

        auto iv = values.begin();
        if (form == regform){
            for(auto p : prior_vals) *iv++ *= p;
        }else{
            for(auto p : prior_vals) *iv++ += p;
        }
    };
}

//! \brief Variation of \b SampleDREAM() which assumes a Bayesian inference problem with likelihood, model and prior.
//! \ingroup DREAMSampleModel

//! The inputs are identical to the ones defined in \b SampleDREAM(), with \b probability_distribution() replaced by the combination of
//! \b TasmanianLikelihood, a \b model lambda and \b prior.
//!
//! - The \b likelihood must inherit from the \b TasmanianLikelihood virtual class, e.g., \b LikelihoodGaussIsotropic
//! - The \b model operates similar to the \b probability_distribution(), but could return multiple outputs for each candidate;
//!   the number of outputs of the \b model must match the dimensions of the \b TasmanianLikelihood (no explicit check performed).
//! - The prior is identical to the one in \b SampleDREAMGrid().
//!
//! Note that \b SampleDREAMGrid() with \b uniform_prior() can also be used to generate samples from a probability distribution that has been
//! approximated using a sparse grid.
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> model,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::function<bool(const std::vector<double> &x)> inside,
                     std::function<void(std::vector<double> &x)> independent_update,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    SampleDREAM<form>(num_burnup, num_collect, makePDFPosterior<form>(likelihood, model, prior),
                      inside, independent_update, state, differential_update, get_random01);
}


//! \brief Overload of \b SampleDREAMPost() that works on a hypercube domain, see the overloads of \b SampleDREAM() for the \b lower and \b upper definition.
//! \ingroup DREAMSampleModel
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> model,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::vector<double> &lower, std::vector<double> &upper,
                     std::function<void(std::vector<double> &x)> independent_update,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    SampleDREAM<form>(num_burnup, num_collect, makePDFPosterior<form>(likelihood, model, prior),
                      lower, upper, independent_update, state, differential_update, get_random01);
}


//! \brief Overload of \b SampleDREAMPost() that uses independent update from a list of known ones, see the overloads of \b SampleDREAM() for details.
//! \ingroup DREAMSampleModel
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> model,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::function<bool(const std::vector<double> &x)> inside,
                     TypeDistribution independent_dist, double independent_magnitude,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    SampleDREAM<form>(num_burnup, num_collect, makePDFPosterior<form>(likelihood, model, prior),
                      inside, independent_dist, independent_magnitude, state, differential_update, get_random01);
}


//! \brief Overload of \b SampleDREAMPost() that uses independent update from a list of known ones and operates on a hypercube domain, see the other overloads for details.
//! \ingroup DREAMSampleModel
template<TypeSamplingForm form = regform>
void SampleDREAMPost(int num_burnup, int num_collect,
                     const TasmanianLikelihood &likelihood,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> model,
                     std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> prior,
                     std::vector<double> &lower, std::vector<double> &upper,
                     TypeDistribution independent_dist, double independent_magnitude,
                     TasmanianDREAM &state,
                     std::function<double(void)> differential_update = const_one,
                     std::function<double(void)> get_random01 = tsgCoreUniform01){
    SampleDREAM<form>(num_burnup, num_collect, makePDFPosterior<form>(likelihood, model, prior),
                      lower, upper, independent_dist, independent_magnitude, state, differential_update, get_random01);
}

}

#endif
