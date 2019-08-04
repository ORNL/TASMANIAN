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

#ifndef __TASMANIAN_DREAM_SAMPLE_HPP
#define __TASMANIAN_DREAM_SAMPLE_HPP

#include "tsgDreamState.hpp"
#include "tsgDreamCoreRandom.hpp"
#include "tsgDreamCorePDF.hpp"

//! \file tsgDreamSample.hpp
//! \brief Core sampling templates.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAM
//!
//! Defines the core MCMC template for sampling from an arbitrarily defined unscaled probability density.

/*!
 * \internal
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMAux Macros and Simall Auxilary functions
 *
 * Several mactos and one-two line function to simplify the work-flow.
 * \endinternal
 */

/*!
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMSampleCore Level 1 DREAM templates, sampling from a user defined distribution
 *
 * Level 1 templates use a custom propability distribution defined wither with a \b lambda.
 * There is no assumption whether the distribution is associated with a Bayesian
 * inference problem or a completely different context.
 */

namespace TasDREAM{

/*!
 * \internal
 * \brief Checks if vectors with names \b lower and \b upper have the same size as the dimensions in TasmanianDREAM \b state.
 * \ingroup DREAMAux
 *
 * Throws \b runtime_error if the size of vectors \b lower and \b upper does not match \b state.getNumDimensions().
 * \endinternal
 */
inline void checkLowerUpper(std::vector<double> const &lower, std::vector<double> const &upper, TasmanianDREAM const &state){
    if (lower.size() != (size_t) state.getNumDimensions()) throw std::runtime_error("ERROR: the size of lower does not match the dimension in state.");
    if (upper.size() != (size_t) state.getNumDimensions()) throw std::runtime_error("ERROR: the size of upper does not match the dimension in state.");
}

//! \internal
//! \brief Returns \b true if the entries in \b x obey the \b lower and \b upper values (sizes must match, does not check).
//! \ingroup DREAMAux

//! The three vectors \b lower, \b upper and \b x must have the same size,
//! returns \b true if every entry of \b x lies between the \b lower and \b upper boundaries.
inline bool inHypercube(const std::vector<double> &lower, const std::vector<double> &upper, const std::vector<double> &x){
    auto il = lower.begin(), iu = upper.begin();
    for(auto v : x) if ((v < *il++) || (v > *iu++)) return false;
    return true;
}

/*!
 * \internal
 * \ingroup DREAMAux
 * \brief Make a lambda that matches the \b inside signature in \b SampleDREAM() and the vector x is in the hyperbube described by \b lower and \b upper.
 * \endinternal
 */
inline std::function<bool(std::vector<double> const &x)> makeHypercudabeLambda(std::vector<double> const &lower, std::vector<double> const &upper){
    return [&](const std::vector<double> &x)->bool{ return inHypercube(lower, upper, x); };
}

/*!
 * \ingroup DREAMSampleCore
 * \brief Make a lambda that matches the \b inside signature in \b SampleDREAM() and the vector x is in the hyperbube described by \b lower and \b upper.
 */
inline std::function<bool(std::vector<double> const &x)> hypercube(std::vector<double> const &lower, std::vector<double> const &upper){
    return [=](const std::vector<double> &x)->bool{
        auto il = lower.begin(), iu = upper.begin();
        for(auto v : x) if ((v < *il++) || (v > *iu++)) return false;
        return true;
    };
}

//! \internal
//! \brief Dummy function that returns 1.0, used as default for the \b differential_update() in \b SampleDREAM().
//! \ingroup DREAMAux

//! Just an inline function that returns 1.0.
inline double const_one(){ return 1.0; }

//! \brief Template that returns a constant based on the percentage, i.e., \b weight_percent / 100.0
//! \ingroup DREAMSampleCore

//! The template simplifies the syntax when calling \b SampleDREAM() with a constant differential update.
//! For example, setting the update to 0.5 can be done with
//! \code
//! TasmanianDREAM state(...);
//! SampleDREAM(..., state, const_percent<50>);
//! \endcode
template<int weight_percent>
double const_percent(){ return ((double) weight_percent) / 100.0; }

//! \brief Uniform prior distribution for both regular and log form.
//! \ingroup DREAMSampleCore

//! Applies uniform (non-informative) prior, can be used with any of the Bayesian inference methods.
//! In practice, this is no-op, since adding zero (in \b logform) or mulitplying by 1 (in \b regform) amounts to nothing.
inline void uniform_prior(const std::vector<double> &, std::vector<double> &values){ values.clear(); }

/*!
 * \ingroup DREAMSampleCore
 * \brief
 */
template<TypeSamplingForm form = regform>
std::function<void(const std::vector<double> &candidates, std::vector<double> &values)>
posterior(std::function<void(TypeSamplingForm, const std::vector<double> &model, std::vector<double> &likely)> &likelihood,
          std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> &model,
          std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> &prior){
    return [&](const std::vector<double> &candidates, std::vector<double> &values)->void{
        std::vector<double> model_outs;
        model(candidates, model_outs);
        likelihood(form, model_outs, values);

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


//! \brief Core template for the sampling algorithm, usually called through any of the overloads.
//! \ingroup DREAMSampleCore

//! Evolves the chains of the \b state using the Metropolis algorithm where the updates are comprised of two components, independent and differential.
//! The independent component relies on independently sampled (pesudo)-random numbers with zero mean and could also be constant zero.
//! The differential component is based on the difference between two randomly chosen chains, which effectively exchanges information between the chains.
//!
//! This is the core algorithm that is oblivious to the source of the probability distribution (i.e., posterior for Bayesian inference based on custom model or sparse grid).
//! - \b form indicates whether the \b probability_distribution() function return the regular form or the logarithm of the desired pdf.
//! - \b num_burnup is the number of initial iterations that will not be saved in the history.
//! - \b num_collect is the number of iterations that will be saved in the \b state history, the total number of collected samples is \b num_collect times \b state.getNumChains().
//! - \b probability_distribution(candidates, values) accepts \b candidates state that consists of multiple points in num_dimensions that correspond to the candidate samples
//!   and must return in \b values the probability density at the points;
//!   the number of points is \b candidates.size() / num_dimensions and will divide evenly, also the the number of points will not exceed the number of chains;
//!   the structure of \b candidates is the same as the input to \b TasmanianSparseGrid::evaluateBatch();
//!   the \b values vector will have size matching the number of samples.
//! - \b inside() accepts a vector of size num_dimensions and returns \b false if the vector is outside of the sampling domain;
//!   \b inside() is called for every candidate point and failing the test will result in automatic rejection;
//!   \b probability_distribution() will never be called with a vector for which \b inside() fails.
//! - \b independent_update() accepts a vector \b x of size num_dimensions and adds a random perturbation to each entry;
//!   the random perturbation must have zero mean (e.g., standard normal distribution);
//!   the perturbation could be deterministic zero, but that may lock the chains into a fixed pattern based on the initial state and thus fail to consider the entire domain.
//! - \b state has initialized state which will be used as the initial iteration, then the state will be evolved num_burnup + num_collect times
//!   and the history will be saved for the last \b num_collect iterations.
//! - \b differential_update() returns a random scaling factor for the differential correction, could be a random or deterministic number, e.g., uniform [0,1] or just 1;
//!   the differential update magnitude defaults to 1.0 and can be set to a constant with the \b const_percent() template;
//!   negative values are statistically equivalent to the positive ones by symmetry of the selection of the chains;
//!   values of more than 1.0 may lead to decrease in acceptance rate and poor mixing (i.e., slow statistical convergence);
//!   setting the update to 0 will result in evolving the chains independently each following the classical Metropolis-Hastings algorithm.
//! - \b get_random01() is the random number generation engine, must return pseudo-random numbers uniformly distributed over [0, 1];
//!   by default, Tasmanian will use the C++ \b rand() function (divided by RAND_MAX).
template<TypeSamplingForm form = regform>
void SampleDREAM(int num_burnup, int num_collect,
                 std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> probability_distribution,
                 std::function<bool(const std::vector<double> &x)> inside,
                 std::function<void(std::vector<double> &x)> independent_update,
                 TasmanianDREAM &state,
                 std::function<double(void)> differential_update = const_one,
                 std::function<double(void)> get_random01 = tsgCoreUniform01){

    if (!state.isStateReady()) throw std::runtime_error("ERROR: DREAM sampling requires that the setState() has been called first on the TasmanianDREAM.");

    if (!state.isPDFReady()) // initialize probability density (if not initialized already)
        state.setPDFvalues(probability_distribution);

    if (num_collect > 0) // pre-allocate memory for the new history
        state.expandHistory(num_collect);

    size_t num_chains = (size_t) state.getNumChains(), num_dimensions = (size_t) state.getNumDimensions();
    double unitlength = (double) num_chains;

    int total_iterations = std::max(num_burnup, 0) + std::max(num_collect, 0);
    for(int t = 0; t < total_iterations; t++){
        std::vector<double> candidates, values;
        candidates.reserve(num_chains * num_dimensions);
        values.reserve(num_chains);

        std::vector<bool> valid(num_chains, true); // keep track whether the samples need to be evaluated

        for(size_t i=0; i<num_chains; i++){
            std::vector<double> propose(num_dimensions);

            size_t jindex = (size_t) (get_random01() * unitlength);
            size_t kindex = (size_t) (get_random01() * unitlength);
            if (jindex >= num_chains) jindex = num_chains - 1; // this is needed in case get_random01() returns 1
            if (kindex >= num_chains) jindex = num_chains - 1;

            state.getIJKdelta(i, jindex, kindex, differential_update(), propose); // propose = s_i + w ( s_k - s_j)
            independent_update(propose); // propose += correction

            if (inside(propose)){
                candidates.insert(candidates.end(), propose.begin(), propose.end());
                values.resize(values.size() + 1);
            }else{
                valid[i] = false;
            }
        }

        probability_distribution(candidates, values);

        std::vector<double> new_state(num_chains * num_dimensions), new_values(num_chains);

        auto icand = candidates.begin(); // loop over all candidates and values, accept or reject
        auto ival = values.begin();

        size_t accepted = 0;

        for(size_t i=0; i<num_chains; i++){
            bool keep_new = valid[i]; // if not valid, automatically reject
            if (valid[i]){ // apply random test
                if (*ival > state.getPDFvalue(i)){ // if the new value has higher probability, automatically accept
                    keep_new = true;
                }else{
                    if (form == regform){
                        keep_new = (*ival / state.getPDFvalue(i) >= get_random01()); // keep if the new value has higher probability
                    }else{
                        keep_new = (*ival - state.getPDFvalue(i) >= log(get_random01()));
                    }
                    //std::cout << "Trsh = " << *ival / state.getPDFvalue(i) << "   " << ((keep_new) ? "Accept" : "Reject") << std:: endl;
                }
            }

            if (keep_new){
                std::copy_n(icand, num_dimensions, new_state.begin() + i * num_dimensions);
                new_values[i] = *ival;
                accepted++; // accepted one more proposal
            }else{ // reject and reuse the old state
                state.getChainState((int) i, &*(new_state.begin() + i * num_dimensions));
                new_values[i] = state.getPDFvalue(i);
            }

            if (valid[i]){ // kept or rejected, if this sample was valid then move to the next sample in the list
                std::advance(icand, num_dimensions);
                ival++;
            }
        }

        state.setState(new_state);
        state.setPDFvalues(new_values);

        if (t >= num_burnup)
            state.saveStateHistory(accepted);
    }
}


//! \brief Overload of \b SampleDREAM() assuming independent update from a list of internals.
//! \ingroup DREAMSampleCore

//! Uses independent update is applied uniformly to all dimensions and comes from a list of internal functions, e.g., uniform or Gaussian.
//! This overload wraps around functions such as \b applyUniformCorrection() and \b applyGaussianCorrection().
template<TypeSamplingForm form = regform>
void SampleDREAM(int num_burnup, int num_collect,
                 std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> probability_distribution,
                 std::function<bool(const std::vector<double> &x)> inside,
                 TypeDistribution dist, double magnitude,
                 TasmanianDREAM &state,
                 std::function<double(void)> differential_update = const_one,
                 std::function<double(void)> get_random01 = tsgCoreUniform01){
    if (dist == dist_uniform){
        SampleDREAM<form>(num_burnup, num_collect, probability_distribution, inside,
                         [&](std::vector<double> &x)->void{ applyUniformUpdate(x, magnitude, get_random01); }, state, differential_update, get_random01);
    }else{ // assuming Gaussian
        SampleDREAM<form>(num_burnup, num_collect, probability_distribution, inside,
                         [&](std::vector<double> &x)->void{ applyGaussianUpdate(x, magnitude, get_random01); }, state, differential_update, get_random01);
    }
}


//! \brief Overload of \b SampleDREAM() assuming the domain is a hyperbube with min/max values given by \b lower and \b upper.
//! \ingroup DREAMSampleCore

//! The two vectors \b lower and \b upper must have the same size as the dimensions of the state
//! and the lower/upper limits of the internals are stores in the corresponding entries.
template<TypeSamplingForm form = regform>
void SampleDREAM(int num_burnup, int num_collect,
                 std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> probability_distribution,
                 const std::vector<double> &lower, const std::vector<double> &upper,
                 std::function<void(std::vector<double> &x)> independent_update,
                 TasmanianDREAM &state,
                 std::function<double(void)> differential_update = const_one,
                 std::function<double(void)> get_random01 = tsgCoreUniform01){
    checkLowerUpper(lower, upper, state);
    SampleDREAM<form>(num_burnup, num_collect, probability_distribution, makeHypercudabeLambda(lower, upper), independent_update, state, differential_update, get_random01);
}


//! \brief Overload of \b SampleDREAM() assuming the domain is a hyperbube and independent update is from a known list.
//! \ingroup DREAMSampleCore

//! The two vectors \b lower and \b upper must have the same size as the dimensions of the state
//! and the lower/upper limits of the internals are stores in the corresponding entries.
//!
//! See other overloads for the \b independent_dist and \b independent_magnitude.
template<TypeSamplingForm form = regform>
void SampleDREAM(int num_burnup, int num_collect,
                 std::function<void(const std::vector<double> &candidates, std::vector<double> &values)> probability_distribution,
                 const std::vector<double> &lower, const std::vector<double> &upper,
                 TypeDistribution independent_dist, double independent_magnitude,
                 TasmanianDREAM &state,
                 std::function<double(void)> differential_update = const_one,
                 std::function<double(void)> get_random01 = tsgCoreUniform01){
    checkLowerUpper(lower, upper, state);
    SampleDREAM<form>(num_burnup, num_collect, probability_distribution, makeHypercudabeLambda(lower, upper), independent_dist, independent_magnitude, state, differential_update, get_random01);
}

}

#endif
