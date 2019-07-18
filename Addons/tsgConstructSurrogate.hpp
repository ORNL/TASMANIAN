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
void constructCommon(std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t thread_id)> model,
                     size_t max_num_samples, size_t num_parallel_jobs,
                     TasmanianSparseGrid &grid,
                     std::function<std::vector<double>(TasmanianSparseGrid &)> candidates,
                     std::string const &checkpoint_filename){
    CandidateManager manager(grid.getNumDimensions()); // keeps track of started and ordered samples
    CompleteStorage complete(grid.getNumDimensions()); // temporarily stores complete samples (batch loading is faster)

    std::string filename = checkpoint_filename;
    std::string filename_old = checkpoint_filename + "_old";
    if (!filename.empty())
        grid.write(filename.c_str(), true); // initial checkpoint

    // prepare several commonly used steps
    auto checkpoint = [&]()->void{ // keeps two saved states for the constructed grid
        if (!filename.empty()){
            { // copy current into old and write to current
                std::ifstream current_state(filename, std::ios::binary);
                std::ofstream previous_state(filename, std::ios::binary);
                previous_state << current_state.rdbuf();
            }
            grid.write(filename.c_str(), true); // write grid to current
        }
    };

    auto load_complete = [&]()->void{ // loads any complete points, does nothing if getNumStored() is zero
        if (complete.getNumStored() > 0){
            complete.load(grid);
            checkpoint();
        }
    };

    auto refresh_candidates = [&]()->void{ // loads complete and asks for new candidates
        load_complete(); // always load before computing new candidates
        manager = candidates(grid); // get new candidates
    };

    auto checkout_sample = [&]()->std::vector<double>{ // get the "most-important" point that has not started yet
        auto x = manager.next();
        if (x.empty()){ // did not find a job, maybe we need to refresh the candidates
            refresh_candidates();
            x = manager.next(); // if this is empty, then we have exhausted all possible candidates
        }
        return x;
    };

    if (!grid.isUsingConstruction())
        grid.beginConstruction();

    refresh_candidates();

    if (parallel_construction){
        // allocate space for all x and y pairs, will be filled by workers and processed by main
        std::vector<std::vector<double>> x(num_parallel_jobs), y(num_parallel_jobs, std::vector<double>(grid.getNumOutputs()));
        size_t total_num_launched = 0; // count all launched jobs

        std::vector<size_t> finished_jobs;
        std::mutex finished_jobs_access;

        // lambda that will handle the work
        auto do_work = [&](size_t thread_id)->void{

            model(x[thread_id], y[thread_id], thread_id); // does the model evaluations

            {
                std::lock_guard<std::mutex> lock(finished_jobs_access);
                finished_jobs.push_back(thread_id);
            }
        };

        // launch initial set of jobs
        std::vector<std::thread> workers(num_parallel_jobs);
        for(size_t id=0; id<num_parallel_jobs; id++){
            x[id] = manager.next();
            if (!x[id].empty()){
                total_num_launched++;
                workers[id] = std::thread(do_work, id);
            }
        }

        while(manager.getNumRunning() > 0){
            // collect finished jobs
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            std::vector<size_t> done; // jobs complete in this iteration
            {
                std::lock_guard<std::mutex> lock(finished_jobs_access);
                std::swap(done, finished_jobs);
            }

            if (done.empty()){ // nothing to do, keep waiting
                std::this_thread::yield();
            }else{ // some threads have finished, process the result
                for(auto id : done){
                    workers[id].join(); // thread has completed
                    if (!x.empty()){
                        complete.add(x[id], y[id]);
                        manager.complete(x[id]);
                    }
                }

                if ((grid.getNumLoaded() < 1000) || (double(complete.getNumStored()) / double(grid.getNumLoaded()) > 0.2))
                    load_complete(); // also does checkpoint save

                // relaunch the threads, if we still need to keep working
                for(auto id : done){
                    if (total_num_launched < max_num_samples){
                        // refresh the candidates if enough of the current candidates have completed
                        if (double(manager.getNumDone()) / double(manager.getNumCandidates()) > 0.2)
                            refresh_candidates();

                        x[id] = checkout_sample(); // if necessary this will call refresh_candidates()

                        if (!x[id].empty()){ // if empty, then we have finished all possible candidates (reached tolerance)
                            total_num_launched++;
                            workers[id] = std::thread(do_work, id);
                        }
                    }
                }

            }
        }

        complete.load(grid);

    }else{
        std::vector<double> x(grid.getNumDimensions()), y(grid.getNumOutputs());

        while((grid.getNumLoaded() + complete.getNumStored() + manager.getNumRunning() < max_num_samples) && (manager.getNumCandidates() > 0)){
            x = manager.next();
            if (x.empty()){ // need more candidates
                refresh_candidates();
                x = manager.next(); // if this is empty, then we have exhausted the candidates
            }
            if (!x.empty()){ // could be empty if there are no more candidates
                model(x, y, 0); // compute a sample
                complete.add(x, y);
                manager.complete(x);

                // the fist thousand points can be loaded one at a time, then add when % increase of the grid is achieved
                if ((grid.getNumLoaded() < 1000) || (double(complete.getNumStored()) / double(grid.getNumLoaded()) > 0.2))
                    load_complete(); // also does checkpoint save
                // if done with the top % of the grid points, recompute the candidates
                if (double(manager.getNumDone()) / double(manager.getNumCandidates()) > 0.2)
                    refresh_candidates();
            }
        }
        complete.load(grid);
    }

    grid.finishConstruction();
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the lambda model.
 */
template<bool parallel_construction = true>
void constructSurrogate(std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t thread_id)> model,
                        size_t max_num_samples, size_t num_parallel_jobs,
                        TasmanianSparseGrid &grid,
                        double tolerance, TypeRefinement criteria, int output = -1,
                        std::vector<int> const &level_limits = std::vector<int>(),
                        std::vector<double> const &scale_correction = std::vector<double>(),
                        std::string const &checkpoint_filename = std::string()){
    constructCommon<parallel_construction>(model, max_num_samples, num_parallel_jobs, grid,
                                           [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(tolerance, criteria, output, level_limits, scale_correction);
                                           }, checkpoint_filename);
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the lambda model.
 */
template<bool parallel_construction = true>
void constructSurrogate(std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t thread_id)> model,
                        size_t max_num_samples, size_t num_parallel_jobs,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, std::vector<int> const &anisotropic_weights = std::vector<int>(),
                        std::vector<int> const &level_limits = std::vector<int>(),
                        std::string const &checkpoint_filename = std::string()){
    constructCommon<parallel_construction>(model, max_num_samples, num_parallel_jobs, grid,
                                           [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(type, anisotropic_weights, level_limits);
                                           }, checkpoint_filename);
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the lambda model.
 */
template<bool parallel_construction = true>
void constructSurrogate(std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t thread_id)> model,
                        size_t max_num_samples, size_t num_parallel_jobs,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, int output, std::vector<int> const &level_limits = std::vector<int>(),
                        std::string const &checkpoint_filename = std::string()){
    constructCommon<parallel_construction>(model, max_num_samples, num_parallel_jobs, grid,
                                           [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(type, output, level_limits);
                                           }, checkpoint_filename);
}

}

#endif
