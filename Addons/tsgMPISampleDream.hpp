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

#ifndef __TASMANIAN_ADDONS_MPISAMPLEDREAM_HPP
#define __TASMANIAN_ADDONS_MPISAMPLEDREAM_HPP

/*!
 * \internal
 * \file tsgMPISampleDream.hpp
 * \brief DREAM posterior sampling using MPI.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianAddonsCommon
 *
 * Allows calling some of DREAM templates in a distributed MPI context.
 * \endinternal
 */

#include "tsgMPIConstructGrid.hpp"

/*!
 * \ingroup TasmanianAddons
 * \addtogroup TasmanianAddonsMPIDream MPI Sampling from Posterior
 *
 * Procedure to collect samples from a posterior distribution using DREAM.
 */

#ifdef Tasmanian_ENABLE_MPI

namespace TasDREAM{

/*!
 * \ingroup TasmanianAddonsMPIDream
 * \brief Class that enables distributed DREAM sampling with MPI.
 *
 * Models with many outputs are computationally expensive to sample when used in a Bayesian posterior.
 * Such models can be distributed across ranks of an MPI communicator and the DREAM candidates
 * can be computed in parallel, i.e., distributing the work.
 * Of specific interest are posteriors with likelihood functions that can be expressed as a product (or sum)
 * of one term per model outputs, e.g., constant or diagonal Gaussian likelihoods.
 * The actual sampling is performed on a single rank while the values of the probability
 * distribution at the candidate points is computed by all MPI ranks together (each ranks handling a separate set of inputs).
 * The high-level algorithm can be expressed with two MPI calls, MPI_Bcast() all candidates
 * to all ranks, then MPI_Reduce() the result.
 * See TasGrid::MPIGridScatterOutputs() and TasDREAM::MPILikelihoodScatter()
 * for ways to distribute a sparse grid model and a likelihood function.
 *
 * The usage of this class is very similar to the calls to TasDREAM::posterior() template,
 * the first two inputs are the local portions of the model and likelihood,
 * followed by the prior which will only be used on the root rank.
 * The other inputs define parameters of the MPI side of the algorithm.
 * The constructor can be inlined in a call to TasDREAM::SampleDREAM().
 *
 * \par MPI Synchronization
 * The call to the constructor will block all ranks except the root,
 * all non-root ranks will enter into a worker cycle waiting for candidates from root.
 * The root process will send-out the unblock signal to all workers whenever the root object is
 * destroyed or the clear() method is called.
 * The class can be passed as the probability distribution to TasDREAM::SampleDREAM(),
 * but only the root rank can compute valid samples, the other ranks will use a no-op function.
 * If SampleDREAM() is called with an empty TasmanianDREAM state then sampling will not be performed,
 * thus the class should be coupled with an empty state on each of the non-root ranks;
 * which helps mirror the code across all ranks and avoid extraneous if-statements.
 * Here are some valid calls:
 *
 * \code
 *   int me;
 *   MPI_Comm_rank(MPI_COMM_WORLD, &me);
 *   auto full_grid = TasGrid::readGrid("foo");
 *   TasGrid::TasmanianSparseGrid grid;
 *   MPIGridScatterOutputs(full_grid, grid, root, 1, MPI_COMM_WORLD); // distribute the grid
 *   TasDREAM::LikelihoodGaussIsotropic full_likely(...);
 *   TasDREAM::LikelihoodGaussIsotropic likely;
 *   MPILikelihoodScatter(full_likely, likely, root, 1, MPI_COMM_WORLD); // distribute the likelihood
 *   TasDREAM::TasmanianDREAM state = (me == root) ? TasmanianDREAM(...) : TasmanianDREAM();
 *   TasDREAM::SampleDREAM(..., DistributedPosterior(likely, grid, ..., root, ...), ..., state, ...);
 * \endcode
 * Alternatively, the object can also be assigned to a variable:
 * \code
 *   TasDREAM::DistributedPosterior post(likely, grid, ..., root, ...);
 *   if (me == root) TasDREAM::SampleDREAM(..., post, ...); // the state on non-root ranks is irrelevant
 *   post.clear(); // unblock the non-root ranks
 * \endcode
 */
template<TypeSamplingForm form = regform>
class DistributedPosterior{
public:
    /*!
     * \brief Constructor that sets the parameters for the distribued posterior.
     *
     * Constructs a distributed posterior object from local model and likelihood objects.
     * See the class description for example usage, details on the parameters follow here:
     *
     * \tparam form is the same as in the call to TasDREAM::SampleDREAM() and both \b must match.
     *
     * \param distributed_model is the local portion of the model, the sub-set of the outputs
     *      must match the set used by the \b likelihood.
     * \param likelihood is the local portion of the likelihood. Same as in the call to TasDREAM::posterior().
     * \param prior will be used only by the root rank, same as in the call to TasDREAM::posterior().
     * \param num_inputs must match the dimensions of the state on the root rank,
     *      since the non-root ranks do not need a valid state this parameter is used to synchronize
     *      the dimensions across all ranks.
     * \param num_chains same as the number set by the state, see the \b num_inputs.
     * \param mpi_root is the root process that will perform the actual sampling.
     * \param communicator is the communicator where all ranks reside.
     */
    DistributedPosterior(DreamModel distributed_model,
                         DreamLikelihood likelihood,
                         DreamPrior prior,
                         int num_inputs, int num_chains, int mpi_root, MPI_Comm communicator)
    : model(distributed_model), likely(likelihood), dist_prior(prior),
      num_dimensions(num_inputs), num_batch(num_chains), root(mpi_root), me(TasGrid::getMPIRank(communicator)), comm(communicator),
      x(Utils::size_mult(num_dimensions, num_batch) + 1), y((size_t) num_batch){

          if (me != root){ // enter work loop
              int num_candidates = 1;
              do{
                MPI_Bcast(x.data(), num_dimensions*num_batch+1, MPI_DOUBLE, root, comm);
                num_candidates = (int) x.back(); // the last entry holds the effective number of candidates
                if (num_candidates > 0){ // if not the shutdown signal
                    x.resize(Utils::size_mult(num_dimensions, num_candidates)); // set x to the correct size
                    y.resize((size_t) num_candidates);

                    std::vector<double> model_outs;
                    model(x, model_outs);
                    likelihood(form, model_outs, y);

                    MPI_Reduce(y.data(), nullptr, num_candidates, MPI_DOUBLE, ((form == regform) ? MPI_PROD : MPI_SUM), root, comm);

                    x.resize(Utils::size_mult(num_dimensions, num_batch) + 1);
                }
              }while(num_candidates > 0);
          }
    }

    //! \brief Destructor, unblocks the non-root ranks (if still blocked).
    ~DistributedPosterior(){ clear(); }

    //! \brief Unblocks the non-root ranks, the object cannot be used after this calls (can be destroyed only).
    void clear(){
        if ((me == root) && (!x.empty())){ // send out the shutdown signal
            x.back() = 0.0;
            MPI_Bcast(x.data(), num_dimensions*num_batch+1, MPI_DOUBLE, root, comm);
        }
    }

    //! \brief Allows passing the object as an input to TasDREAM::SampleDREAM().
    operator DreamPDF(){
        if (me == root){
            return [&](const std::vector<double> &candidates, std::vector<double> &values)->void{
                std::copy_n(candidates.begin(), candidates.size(), x.begin());
                int num_candidates = (int) candidates.size() / num_dimensions;
                x.back() = (double) num_candidates;
                MPI_Bcast(x.data(), num_dimensions*num_batch+1, MPI_DOUBLE, root, comm);

                y.resize((size_t) num_candidates);
                std::vector<double> model_outs;
                model(candidates, model_outs);
                likely(form, model_outs, y);

                MPI_Reduce(y.data(), values.data(), num_candidates, MPI_DOUBLE, ((form == regform) ? MPI_PROD : MPI_SUM), root, comm);

                std::vector<double> prior_vals(values.size());
                dist_prior(form, candidates, prior_vals);

                auto iv = values.begin();
                if (form == regform){
                    for(auto p : prior_vals) *iv++ *= p;
                }else{
                    for(auto p : prior_vals) *iv++ += p;
                }
            };
        }else{ // no-op distribution
            return [](const std::vector<double> &, std::vector<double> &)->void{};
        }
    }

private:
    std::function<void(std::vector<double> const &x, std::vector<double> &y)> model;
    std::function<void(TypeSamplingForm, const std::vector<double> &model_outputs, std::vector<double> &likely)> likely;
    std::function<void(TypeSamplingForm, const std::vector<double> &candidates, std::vector<double> &values)> dist_prior;
    int num_dimensions, num_batch, root, me;
    MPI_Comm comm;
    std::vector<double> x, y;
};

}

#endif // Tasmanian_ENABLE_MPI

#endif
