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

#ifndef __TASMANIAN_ADDONS_CONSTRUCT_SURROGATE_HPP
#define __TASMANIAN_ADDONS_CONSTRUCT_SURROGATE_HPP

/*!
 * \internal
 * \file tsgConstructSurrogate.hpp
 * \brief Automated parallel construction procedure.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianAddonsCommon
 *
 * Templates for automated surrogate construction.
 * \endinternal
 */

#include "tsgCandidateManager.hpp"

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianAddonsConstruct
 * \brief Construction algorithm using generic candidates generation.
 *
 * \endinternal
 */
template<bool parallel_construction>
void constructCommon(std::function<void(std::vector<double> const &x, std::vector<double> &y)> model,
                     size_t max_num_samples, size_t num_parallel_jobs,
                     TasmanianSparseGrid &grid,
                     std::function<std::vector<double>(TasmanianSparseGrid &)> candidates){
    CandidateManager manager(grid.getNumDimensions());
    CompleteStorage complete(grid.getNumDimensions());

    if (!grid.isUsingConstruction())
        grid.beginConstruction();

    manager = candidates(grid);

    std::vector<double> x(grid.getNumDimensions()), y(grid.getNumOutputs());

    auto has_budget = [&]()->bool{
        return (max_num_samples == -1) || (grid.getNumLoaded() + complete.getNumStored() + manager.getNumRunning() < max_num_samples);
    };

    while(has_budget() && (manager.getNumCandidates() > 0)){
        x = manager.next();
        if (x.empty()){ // need more candidates
            complete.load(grid);
            manager = candidates(grid);
            x = manager.next(); // if this is empty, then we have exhausted the candidates
        }
        if (!x.empty()){ // could be empty if there are no more candidates
            model(x, y); // compute a sample
            complete.add(x, y);
            manager.complete(x);

            // the fist thousand points can be loaded one at a time, then add when % increase of the grid is achieved
            if ((grid.getNumLoaded() < 1000) || (double(complete.getNumStored()) / double(grid.getNumLoaded()) > 0.2))
                complete.load(grid);
            // if done with the top % of the grid points, recompute the candidates
            if (double(manager.getNumDone()) / double(manager.getNumCandidates()) > 0.2){
                complete.load(grid); // always load before computing new candidates
                manager = candidates(grid);
            }
        }
    }
    complete.load(grid);

    grid.finishConstruction();
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the lambda model.
 */
template<bool parallel_construction = true>
void constructSurrogate(std::function<void(std::vector<double> const &x, std::vector<double> &y)> model,
                        size_t max_num_samples, size_t num_parallel_jobs,
                        TasmanianSparseGrid &grid,
                        double tolerance, TypeRefinement criteria, int output = -1,
                        std::vector<int> const &level_limits = std::vector<int>(),
                        std::vector<double> const &scale_correction = std::vector<double>()){
    constructCommon<parallel_construction>(model, max_num_samples, num_parallel_jobs, grid,
                                           [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(tolerance, criteria, output, level_limits, scale_correction);
                                           });
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the lambda model.
 */
template<bool parallel_construction = true>
void constructSurrogate(std::function<void(std::vector<double> const &x, std::vector<double> &y)> model,
                        size_t max_num_samples, size_t num_parallel_jobs,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, std::vector<int> const &anisotropic_weights = std::vector<int>(),
                        std::vector<int> const &level_limits = std::vector<int>()){
    constructCommon<parallel_construction>(model, max_num_samples, num_parallel_jobs, grid,
                                           [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(type, anisotropic_weights, level_limits);
                                           });
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the lambda model.
 */
template<bool parallel_construction = true>
void constructSurrogate(std::function<void(std::vector<double> const &x, std::vector<double> &y)> model,
                        size_t max_num_samples, size_t num_parallel_jobs,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, int output, std::vector<int> const &level_limits = std::vector<int>()){
    constructCommon<parallel_construction>(model, max_num_samples, num_parallel_jobs, grid,
                                           [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(type, output, level_limits);
                                           });
}


}

#endif
