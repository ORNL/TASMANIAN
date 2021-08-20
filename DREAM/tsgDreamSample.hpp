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

/*!
 * \internal
 * \file tsgDreamSample.hpp
 * \brief Core sampling templates.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAM
 *
 * Defines the core MCMC template for sampling from an arbitrarily defined unscaled probability density.
 * \endinternal
 */

/*!
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMSampleCore DREAM Sampling Templates
 *
 * Templates and auxiliary methods for the DREAM sampling.
 * The main template is TasDREAM::SampleDREAM() with one overload
 * and several helper functions.
 * The helpers provide ways to define the probability distribution:
 * either custom defined, interpolated with a sparse grid, or product of Bayesian inference problem.
 */

namespace TasDREAM{

/*!
 * \ingroup DREAMSampleCore
 * \brief Generic test function whether a sample belongs in the domain.
 *
 * The function accepts a single vector with size equal to the number of domain dimensions
 * and must return \b true if the vector belongs to the domain or \b false otherwise.
 */
using DreamDomain = std::function<bool(std::vector<double> const &x)>;

/*!
 * \ingroup DREAMSampleCore
 * \brief Make a lambda that matches the \b inside signature in \b SampleDREAM(), test if the vector x is in the hyperbube described by \b lower and \b upper.
 */
inline DreamDomain hypercube(std::vector<double> const &lower, std::vector<double> const &upper){
    return [=](const std::vector<double> &x)->bool{
        auto il = lower.begin(), iu = upper.begin();
        for(auto v : x) if ((v < *il++) || (v > *iu++)) return false;
        return true;
    };
}

/*!
 * \ingroup DREAMSampleCore
 * \brief Dummy function that does not make any changes to the vector as default for the \b independent_update() in \b SampleDREAM().
 *
 * The function is no-op.
 */
inline void no_update(std::vector<double> &){}

/*!
 * \brief Dummy function that returns 1.0, used as default for the \b differential_update() in \b SampleDREAM().
 * \ingroup DREAMSampleCore
 *
 * Just an inline function that returns 1.0.
 */
inline double const_one(){ return 1.0; }

/*!
 * \brief Template that returns a constant based on the percentage, i.e., \b weight_percent / 100.0
 * \ingroup DREAMSampleCore
 *
 * The template simplifies the syntax when calling \b SampleDREAM() with a constant differential update.
 * For example, setting the update to 0.5 can be done with
 * \code
 * TasmanianDREAM state(...);
 * SampleDREAM(..., independent_update, const_percent<50>, state);
 * \endcode
 */
template<int weight_percent>
double const_percent(){ return ((double) weight_percent) / 100.0; }

/*!
 * \brief Uniform prior distribution for both regular and log form.
 * \ingroup DREAMSampleCore
 *
 * Applies uniform (non-informative) prior, can be used with any of the Bayesian inference methods.
 * In practice, this actually does nothing, since adding zero (in \b logform) or mulitplying by 1 (in \b regform) amounts to nothing.
 */
inline void uniform_prior(TypeSamplingForm, const std::vector<double> &, std::vector<double> &values){ values.clear(); }

/*!
 * \ingroup DREAMSampleCore
 * \brief Generic probability distribution used by Tasmanian.
 *
 * The probability distribution must be set to accept multiple candidates and return
 * the value of the unscaled probability distribution for each point.
 *
 * \param candidates is a vector with size that is an even multiple of the dimensions,
 *      the vector is organized logically into strips of size num-dimensions.
 *      Each strip represents a test point that is guaranteed to be inside the domain.
 * \param values is a vector with size equal to the number of strips (samples) in \b candidates.
 *      The vector should not be resized, instead each value has to be overwritten
 *      with the corresponding unscaled probability distribution.
 *
 * \b Note: the generic probability distribution does not accept a parameter to specify the sampling
 *      form, e.g., TasDREAM::logform. The returned values \b must corresponding to the sampling
 *      form set in the TasDREAM::SampleDREAM() template. If logarithm form is used, the values
 *      can be negative, in regular form the values must be positive, Tasmanian will not throw
 *      an exception but using negative values with TasDREAM::regform leads to undefined behavior.
 */
using DreamPDF = std::function<void(const std::vector<double> &candidates, std::vector<double> &values)>;

/*!
 * \ingroup DREAMSampleCore
 * \brief Generic model signature used by Tasmanian.
 *
 * The model is very similar to the TasDREAM::DreamPDF, in fact the input \b candidates is the same.
 * The differences are two:
 * - the model may have multiple outputs and the number of returned outputs must match
 *   the number of outputs used by the likelihood.
 * - the \b outputs vector will not be set to the correct size, it \b must be resized,
 *   because Tasmanian does not keep track of the number of outputs.
 *
 * \b Note: the \b outputs vector will be fed as input to the TasDREAM::DreamLikelihood.
 */
using DreamModel = std::function<void(const std::vector<double> &candidates, std::vector<double> &outputs)>;

/*!
 * \ingroup DREAMSampleCore
 * \brief Generic likelihood signature used by Tasmanian.
 *
 * The likelihood assigns a value to how likely the model outputs are given data (i.e., measurements).
 * Classes that inherit TasmanianLikelihood will automatically convert to a lambda object with
 * this signature.
 *
 * \param form is the sampling type used in the call to TasDREAM::SampleDREAM(), in a custom object
 *      it is sufficient to implement only one form, Tasmanian likelihood classes implement both
 *      and thus the parameter specifies which implementation to use.
 * \param model_outputs is identical to \b outputs parameter in TasDREAM::DreamModel
 *      except here it is used as an input to compute the likelihood.
 * \param likely is a vector with size equal to the entries with \b outputs,
 *      the entries must be overwritten with the corresponding likelihood.
 */
using DreamLikelihood = std::function<void(TypeSamplingForm form, const std::vector<double> &model_outputs, std::vector<double> &likely)>;

/*!
 * \ingroup DREAMSampleCore
 * \brief Generic signature for the prior distributions used by Tasmanian.
 *
 * Specifies the prior distribution, the TasDREAM::uniform_prior() satisfies this signature.
 * \param form is the same as in TasDREAM::DreamLikelihood.
 * \param candidates is the same as in TasDREAM::DreamModel and TasDREAM::DreamPDF.
 * \param values are similar to the \b likely in TasDREAM::DreamLikelihood, but instead of the likelihood
 *      the values define the prior distribution in the specified \b form.
 */
using DreamPrior = std::function<void(TypeSamplingForm form, const std::vector<double> &candidates, std::vector<double> &values)>;

/*!
 * \ingroup DREAMSampleCore
 * \brief Generic signature for a combination of a likelihood and a model.
 *
 * The likelihood and the model are not always separated, e.g., a sparse grid
 * approximation can be used to interpolated the likelihood which has a single output
 * and is therefore cheaper than interpolating multi-output model.
 * The implementation should be equivalent to:
 * \param candidates is same as in TasDREAM::DreamModel.
 * \param values is same as in TasDREAM::DreamLikelihood \b likely parameter.
 */
using DreamMergedLikelyModel = std::function<void(const std::vector<double> &candidates, std::vector<double> &values)>;

/*!
 * \ingroup DREAMSampleCore
 * \brief Combines the three components of a Bayesian posterior into a single distribution.
 *
 * The Bayesian posterior has a model, likelihood function and a prior distribution.
 * This function combines three function objects into a single probability distribution
 * that be passed to SampleDREAM, e.g.,
 * \code
 * SampleDREAM(num_burnup, num_collect, posterior(likely, model, uniform_prior), ...);
 * \endcode
 *
 * The generality of the approach used here comes at the price of volatility.
 * There is no builtin error-checking and error-detection on vector sizes.
 * Specifically, the number of inputs provided by the model must match
 * the outputs accepted by the likelihood, and the number of dimensions
 * accepted by the model and prior must be the same.
 *
 * \tparam form indicates whether to use the regular or logarithm of the sampling problem;
 *      Gaussian-types of likelihood functions are often used where the regular
 *      form can be \f$ \exp( -0.5 x^2 ) \f$ while the log-form is \f$ -0.5 x^2 \f$,
 *      working with a simple quadratic function can be more stable with respect
 *      to rounding error.
 *      The template parameter \b must be the same as in the call to TasDREAM::SampleDREAM().
 *
 * \param model accepts a set of model inputs and will return the corresponding model values.
 *      - \b candidates is the same as in the input of \b probability_distribution() in
 *        the call to TasDREAM::SampleDREAM().
 *        Logically the candidates will be arranged in strips of size equal to the problem
 *        dimensions, the vector size will divide evenly by the dimensions and the factor
 *        is the number of candidates.
 *      - \b outputs must be resized to match the number of candidates times the number of outputs,
 *        the behavior must match that of TasmanianSparseGrid::evaluateBatch().
 *
 * \param likelihood accepts a set of model outputs and provides a measure of
 *      how likely those outputs are given some observed data with some noise.
 *      Tasmanian provides likelihood functions that can be used here, e.g.,
 *      TasDREAM::LikelihoodGaussIsotropic and TasDREAM::LikelihoodGaussAnisotropic.
 *      - The TasDREAM::TypeSamplingForm will always match the template parameter \b form,
 *        thus, it is sufficient to implement only one sampling from. The Tasmanian likelihood
 *        classes implement both forms, hence the added flexibility.
 *      - The \b model_outputs is a vector with size equal to the number of candidates times
 *        the number of outputs, i.e., must match the output of the \b model.
 *      - The \b likely will have size equal to the number of candidates and must
 *        be filled (without resize) with the likelihood for each set of model outputs.
 *
 * \param prior provides the values of the prior distribution in either regular or logarithm form.
 *      The prior will take in the same \b candidates as the model and a vector of the same size as \b likely,
 *      and must return the values of the corresponding prior distribution in either regular or logarithm form.
 *
 * Example usage:
 * \code
 *  auto model = TasGrid::read("foo");
 *  TasDREAM::LikelihoodGaussIsotropic likely(0.1, data);
 *  TasDREAM::SampleDREAM(..., TasDREAM::posterior(likely, model, TasDREAM::uniform_prior), ...);
 * \endcode
 *
 */
template<TypeSamplingForm form = regform>
DreamPDF posterior(DreamModel model,
          DreamLikelihood likelihood,
          DreamPrior prior){
    return [=](const std::vector<double> &candidates, std::vector<double> &values)->void{
        std::vector<double> model_outs;
        model(candidates, model_outs);
        likelihood(form, model_outs, values);

        std::vector<double> prior_vals(values.size());
        prior(form, candidates, prior_vals);

        auto iv = values.begin();
        if (form == regform){
            for(auto p : prior_vals) *iv++ *= p;
        }else{
            for(auto p : prior_vals) *iv++ += p;
        }
    };
}

/*!
 * \ingroup DREAMSampleCore
 * \brief Overload where the model and likelihood are combined into a single call.
 *
 * There are situations where splitting the model and likelihood is undesirable,
 * e.g., if the model has a huge number of outputs it may be easier to construct
 * a sparse grid surrogate to the single-output combined model and likelihood.
 * This is a short hand-template that uses such model-likelihood and combines
 * it with a prior distribution.
 */
template<TypeSamplingForm form = regform>
DreamPDF posterior(DreamMergedLikelyModel likelihood_model,
          DreamPrior prior){
    return [=](const std::vector<double> &candidates, std::vector<double> &values)->void{
        likelihood_model(candidates, values);

        std::vector<double> prior_vals(values.size());
        prior(form, candidates, prior_vals);

        auto iv = values.begin();
        if (form == regform){
            for(auto p : prior_vals) *iv++ *= p;
        }else{
            for(auto p : prior_vals) *iv++ += p;
        }
    };
}


/*!
 * \brief Core template for the sampling algorithm.
 * \ingroup DREAMSampleCore
 *
 * Evolves the chains of the \b state using the Metropolis algorithm where the updates are comprised of two components, independent and differential.
 * The independent component relies on independently sampled (pesudo)-random numbers with zero mean (could be deterministic constant zero).
 * The differential component is based on the difference between two randomly chosen chains, which effectively exchanges information between the chains.
 *
 * The implementation is very generic using lambda objects to describe most aspects of the problem,
 * including the probability distribution, domain, etc.
 * However, the generality comes with some sacrifice in resilience, i.e.,
 * - each lambda object must respect the problem dimensions, the template will populate the inputs with the correct
 *   number of entries, but the lambdas have to properly utilize the entries.
 * - each lambda will have to capture external objects and variables and should avoid capture by value for large
 *   vectors, e.g., a sparse grid or a large data vector, while objects captured by reference must remain alive
 *   during the call to the template.
 * See couple of examples after the variable listing.
 *
 * \tparam form indicates whether the \b probability_distribution() function return the regular form or the logarithm of the desired pdf;
 *      in some cases, taking the exponential of a very large negative number can lead to problems with rounding comparison
 *      between numbers very close to zero, hence the logarithm from could be more stable from the round-off error perspective.
 *
 * \param num_burnup is the number of initial iterations that will not be saved in the history;
 *      the Metropolis algorithm guarantees convergence to the desired distribution but only in the limit;
 *      when working with a finite number of samples the results can be contaminated by the initial state
 *      which can be significantly different from the limit.
 *      A common practice is to take \b num_burnup to be roughly equal to \b num_collect, but that is only
 *      a heuristic suggestion.
 *
 * \param num_collect is the number of iterations that will be saved in the \b state history,
 *      the total number of collected samples is \b num_collect times \b state.getNumChains().
 *
 * \param probability_distribution must accept a vector of \b candidates locations and return the \b values
 *      of the unscaled probability distribution at those points.
 *      The input-output relation is similar to working with TasmanianSparseGrid::evaluateBatch()
 *      when the grid has a a single output.
 *      - The \b candidates input will be logically divided into strips of size \b state.getNumDimensions(),
 *        the total number of strips will be between 1 and \b state.getNumChains() depending
 *        on the number of proposals that fall within the domain, i.e.,
 *        the inputs will always satisfy the \b inside() condition.
 *      - The \b values will have size equal to the number of strips and each entry should be filled
 *        the corresponding value of the probability density function.
 *        The \b values vector should not be resized.
 *      - Any TasGrid::TasmanianSparseGrid object with a single output and non-zero loaded points
 *        can be passed in as a \b probability_distribution; when using TasDREAM::regform
 *        the output of the grid must be always non-negative, using TasDREAM::logform has no lower bound requirements.
 *      - In the context of Bayesian inference, there probability distribution is comprised
 *        of three components: model, likelihood, and prior. See TasDREAM::posterior().
 *
 * \param inside must accept a vector of size equal tot he number of dimensions and return \b true
 *      if the vector describes a point in the domain and \b false otherwise.
 *      The function will be called for every proposed sample and unlike the \b probability_distribution
 *      only one vector will be given as input at a time.
 *      Each TasGrid::TasmanianSparseGrid object some with a canonical domain that depends on the rule and
 *      can be optionally transformed, the TasmanianSparseGrid::getDomainInside() method will produce
 *      a lambda object describing the sparse grid domain.
 *      See also TasDREAM::hyperbube().
 *
 * \param state must be an initialized instance of TasmanianDREAM. The number of chains and dimensions must be set
 *      to match the lambda objects and TasmanianDREAM::setState() must be called to load the state with
 *      a valid initial set of vectors.
 *      The \b state will be evolved in total of \b num_burnup + \b num_collect iterations
 *      and the last set of \b num_collect iterations will be appended to the state history.
 *
 * \param independent_update must accept a vector of size \b state.getNumDimensions() and (without resizing)
 *      perturb each entry by adding an independent zero-mean random number.
 *      It is possible to omit the \b independent_update, e.g., pass an object that will not modify
 *      the input vector x; however, this can lock the chains to a set of values not dense in the domain
 *      which in turn will break convergence. For example, if the domain is the real line, the initial
 *      state has integer entries, and the \b differential_update is set to 100%, then all proposals
 *      will be integers, non-integer values will never be considered.
 *      TasDREAM::SampleDREAM provides an overload where the independent update is set to
 *      a distribution from an included list with known magnitude.
 *
 * \param differential_update is a random or deterministic number between zero and one,
 *      indicating the magnitude of the difference between two randomly chosen chains that will be
 *      added to the next proposal.
 *      Using deterministic zero will decouple the chains, i.e., the state of one chain will
 *      never affect any other chain and the algorithm will reduce to the propagation of
 *      multiple independent chains of the standard Metropolis-Hastings method with
 *      the selected \b independent_update.
 *      Using negative values (by symmetry) is equivalent to using a positive value with the same magnitude,
 *      and values larger than 1 are likely to result in poor mixing and bad convergence rate.
 *      The default differential update is deterministic one and deterministic percentage
 *      can be specified with the \b const_percent() template.
 *
 * \param get_random01 is the pseudo-random number generator to be used in the sampling procedure.
 *      By default, Tasmanian will use \b rand() divided by \b RAND_MAX, but this is implementation
 *      dependent and not always optimal.
 *
 * Correct call using a sparse grid object as input:
 * \code
 *   auto grid = TasGrid::read("foo"); // create a grid object
 *   SampleDREAM(..., grid, ...);
 *   // above, the grid will create a lambda object that will capture grid by reference
 *   // the lambda object is destroyed at the end of the call and before grid
 * \endcode
 * Incorrect call, \b never \b do \b that:
 * \code
 * SampleDREAM(..., TasGrid::read("foo"), ...); // <- Incorrect
 * // above, read() will create a TasmanianSparseGrid object which will create a lambda
 * // but the grid will be destroyed before entering the SampleDREAM() and cause a segfault.
 * \endcode
 */
template<TypeSamplingForm form = regform>
void SampleDREAM(int num_burnup, int num_collect,
                 DreamPDF probability_distribution,
                 DreamDomain inside,
                 TasmanianDREAM &state,
                 std::function<void(std::vector<double> &x)> independent_update = no_update,
                 std::function<double(void)> differential_update = const_one,
                 std::function<double(void)> get_random01 = tsgCoreUniform01){

    size_t num_chains = (size_t) state.getNumChains(), num_dimensions = (size_t) state.getNumDimensions();
    double unitlength = (double) num_chains;

    if (num_chains == 0) return; // no sampling with a null state

    if (!state.isStateReady()) throw std::runtime_error("ERROR: DREAM sampling requires that the setState() has been called first on the TasmanianDREAM.");

    if (!state.isPDFReady()) // initialize probability density (if not initialized already)
        state.setPDFvalues(probability_distribution);

    if (num_collect > 0) // pre-allocate memory for the new history
        state.expandHistory(num_collect);

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

        if (!candidates.empty()) // block the pathological case of all proposals leaving the domain
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


/*!
 * \ingroup DREAMSampleCore
 * \brief Overload of \b SampleDREAM() assuming independent update from a list of internally implemented options.
 *
 * Uses independent update is applied uniformly to all dimensions
 * and comes from a list of internal functions, e.g., uniform or Gaussian.
 * This overload wraps around functions such as
 * \b TasDREAM::applyUniformCorrection() and \b TasDREAM::applyGaussianCorrection().
 */
template<TypeSamplingForm form = regform>
void SampleDREAM(int num_burnup, int num_collect,
                 DreamPDF probability_distribution,
                 DreamDomain inside,
                 TasmanianDREAM &state,
                 TypeDistribution dist, double magnitude,
                 std::function<double(void)> differential_update = const_one,
                 std::function<double(void)> get_random01 = tsgCoreUniform01){
    if (dist == dist_uniform){
        SampleDREAM<form>(num_burnup, num_collect, probability_distribution, inside, state,
                         [&](std::vector<double> &x)->void{ applyUniformUpdate(x, magnitude, get_random01); }, differential_update, get_random01);
    }else if (dist == dist_gaussian){
        SampleDREAM<form>(num_burnup, num_collect, probability_distribution, inside, state,
                         [&](std::vector<double> &x)->void{ applyGaussianUpdate(x, magnitude, get_random01); }, differential_update, get_random01);
    }else{ // assuming none
        SampleDREAM<form>(num_burnup, num_collect, probability_distribution, inside, state, no_update, differential_update, get_random01);
    }
}

}

#endif
