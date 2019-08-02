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

template<TypeSamplingForm form = regform>
class DistributedPosterior{
public:
    DistributedPosterior(std::function<void(TypeSamplingForm, const std::vector<double> &model_outputs, std::vector<double> &likely)> likelihood,
                         std::function<void(const std::vector<double> &candidates, std::vector<double> &outputs)> distributed_model,
                         std::function<void(TypeSamplingForm, const std::vector<double> &candidates, std::vector<double> &values)> prior,
                         int num_inputs, int num_chains, int mpi_root, MPI_Comm communicator)
    : likely(likelihood), model(distributed_model), dist_prior(prior),
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

    ~DistributedPosterior(){ clear(); }

    void clear(){
        if ((me == root) && (!x.empty())){ // send out the shutdown signal
            x.back() = 0.0;
            MPI_Bcast(x.data(), num_dimensions*num_batch+1, MPI_DOUBLE, root, comm);
        }
    }

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
    std::function<void(TypeSamplingForm, const std::vector<double> &model_outputs, std::vector<double> &likely)> likely;
    std::function<void(std::vector<double> const &x, std::vector<double> &y)> model;
    std::function<void(TypeSamplingForm, const std::vector<double> &candidates, std::vector<double> &values)> dist_prior;
    int num_dimensions, num_batch, root, me;
    MPI_Comm comm;
    std::vector<double> x, y;
};

}

#endif // Tasmanian_ENABLE_MPI

#endif
