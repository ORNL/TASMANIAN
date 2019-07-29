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
 * \ingroup TasmanianAddonsConstruct
 * \brief Allows for expressive calls to TasGrid::constructSurrogate().
 */
constexpr bool mode_parallel = true;

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Allows for expressive calls to TasGrid::constructSurrogate().
 */
constexpr bool mode_sequential = false;

/*!
 * \internal
 * \ingroup TasmanianAddonsConstruct
 * \brief Used for expressive template calls, indicates model that accepts an initial guess.
 *
 * \endinternal
 */
struct AcceleratedMode{};
/*!
 * \internal
 * \ingroup TasmanianAddonsConstruct
 * \brief Used for expressive template calls, indicates model that does not accept initial guess.
 *
 * \endinternal
 */
struct StaticMode{};

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Model described by inputs \b x, outputs \b y, initial guess \b y0, and running on thread with \b thread_id, see TasGrid::constructSurrogate().
 */
using AcceleratedModelLambda = std::function<void(std::vector<double> const &x, std::vector<double> &y, std::vector<double> const &y0, size_t thread_id)>;
/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Model described by inputs \b x, outputs \b y, no initial guess, and running on thread with \b thread_id, see TasGrid::constructSurrogate().
 */
using SimpleModelLambda = std::function<void(std::vector<double> const &, std::vector<double> &, size_t)>;

/*!
 * \internal
 * \ingroup TasmanianAddonsConstruct
 * \brief Takes a lambda of a model without an initial guess and converts to model that ignored the initial guess.
 *
 * \endinternal
 */
AcceleratedModelLambda modelNoInit(SimpleModelLambda &model){
    return [&](std::vector<double> const &x, std::vector<double> &y, std::vector<double> const &, size_t thread_id)->void{ model(x, y, thread_id); };
}

/*!
 * \internal
 * \ingroup TasmanianAddonsConstruct
 * \brief Construction algorithm using generic candidates procedure.
 *
 * Implements the parallel and sequential algorithms into one template that uses
 * the \b candidates lambda that wraps around the appropriate
 * TasmanianSparseGrid::getCandidateConstructionPoints() overload.
 * \endinternal
 */
template<bool parallel_construction, class InitialGuessMode>
void constructCommon(AcceleratedModelLambda model,
                     size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                     TasmanianSparseGrid &grid,
                     std::function<std::vector<double>(TasmanianSparseGrid &)> candidates,
                     std::string const &checkpoint_filename){
    num_parallel_jobs   = std::max(size_t(1), num_parallel_jobs);
    max_samples_per_job = std::max(size_t(1), max_samples_per_job);
    size_t num_dimensions = (size_t) grid.getNumDimensions();
    CandidateManager manager(num_dimensions, max_samples_per_job); // keeps track of started and ordered samples
    CompleteStorage complete(num_dimensions); // temporarily stores complete samples (batch loading is faster)

    std::string filename = checkpoint_filename;
    std::string filename_old = checkpoint_filename + "_old";

    if (!filename.empty()){ // recover from an existing checkpoint
        auto save_grid = std::move(grid);
        std::ifstream infile(filename, std::ios::binary);
        try{ // attempt to recover from filename
            if (!infile.good()) throw std::runtime_error("missing main checkpoint");
            grid.read(infile, mode_binary);
            complete.read(infile);
        }catch(std::runtime_error &){
            // main file is missing or is corrupt, try the older version
            std::ifstream oldfile(filename_old, std::ios::binary);
            try{
                if (!oldfile.good()) throw std::runtime_error("missing main checkpoint");
                grid.read(oldfile, mode_binary);
                complete.read(oldfile);
            }catch(std::runtime_error &){
                // nothing could be recovered, start over
                grid = std::move(save_grid);
            }
        }
    }

    if (!filename.empty()){ // initial checkpoint
        std::ofstream ofs(filename, std::ios::binary);
        grid.write(ofs, mode_binary); // write grid to current
        complete.write(ofs);
    }

    // prepare several commonly used steps
    auto checkpoint = [&]()->void{ // keeps two saved states for the constructed grid
        if (!filename.empty()){
            { // copy current into old and write to current
                std::ifstream current_state(filename, std::ios::binary);
                std::ofstream previous_state(filename, std::ios::binary);
                previous_state << current_state.rdbuf();
            }
            std::ofstream ofs(filename, std::ios::binary);
            grid.write(ofs, mode_binary); // write grid to current
            complete.write(ofs);
        }
    };

    auto load_complete = [&]()->void{ // loads any complete points, does nothing if getNumStored() is zero
        if (complete.getNumStored() > 0)
            complete.load(grid);
    };

    auto refresh_candidates = [&]()->void{ // loads complete and asks for new candidates
        load_complete(); // always load before computing new candidates
        manager = candidates(grid); // get new candidates
    };

    size_t total_num_launched = complete.getNumStored() + grid.getNumLoaded(); // count all launched jobs, including the ones already complete
    auto checkout_sample = [&]()->std::vector<double>{ // get the "most-important" point that has not started yet
        auto x = manager.next(max_num_points - total_num_launched);
        if (x.empty()){ // did not find a job, maybe we need to refresh the candidates
            refresh_candidates();
            x = manager.next(max_num_points - total_num_launched); // if this is empty, then we have exhausted all possible candidates
        }
        return x;
    };

    if (!grid.isUsingConstruction())
        grid.beginConstruction();

    refresh_candidates();

    if (parallel_construction == mode_parallel){
        // allocate space for all x and y pairs, will be filled by workers and processed by main
        std::vector<std::vector<double>> x(num_parallel_jobs),
                                         y(num_parallel_jobs, std::vector<double>(grid.getNumOutputs())),
                                         y0(num_parallel_jobs);

        std::vector<size_t> finished_jobs;
        std::mutex finished_jobs_access;

        // lambda that will handle the work
        auto do_work = [&](size_t thread_id)->void{

            model(x[thread_id], y[thread_id], y0[thread_id], thread_id); // does the model evaluations

            {
                std::lock_guard<std::mutex> lock(finished_jobs_access);
                finished_jobs.push_back(thread_id);
            }
        };

        // launch initial set of jobs
        std::vector<std::thread> workers(num_parallel_jobs);
        for(size_t id=0; id<num_parallel_jobs; id++){
            x[id] = manager.next(max_num_points - total_num_launched);
            if (!x[id].empty()){
                total_num_launched += x[id].size() / num_dimensions;
                if (std::is_same<InitialGuessMode, AcceleratedMode>::value && (grid.getNumLoaded() > 0))
                    grid.evaluateBatch(x[id], y0[id]);
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
                    if (total_num_launched < max_num_points){
                        // refresh the candidates if enough of the current candidates have completed
                        if (double(manager.getNumDone()) / double(manager.getNumCandidates()) > 0.2)
                            refresh_candidates();

                        x[id] = checkout_sample(); // if necessary this will call refresh_candidates()

                        if (!x[id].empty()){ // if empty, then we have finished all possible candidates (reached tolerance)
                            total_num_launched += x[id].size() / num_dimensions;
                            if (std::is_same<InitialGuessMode, AcceleratedMode>::value && (grid.getNumLoaded() > 0))
                                grid.evaluateBatch(x[id], y0[id]);
                            workers[id] = std::thread(do_work, id);
                        }
                    }
                    checkpoint();
                }

            }
        }

        complete.load(grid);

    }else{
        std::vector<double> x(grid.getNumDimensions()), y(grid.getNumOutputs()), y0(grid.getNumOutputs());

        while((total_num_launched < max_num_points) && (manager.getNumCandidates() > 0)){
            x = manager.next(max_num_points - total_num_launched);
            if (x.empty()){ // need more candidates
                refresh_candidates();
                x = manager.next(max_num_points - total_num_launched); // if this is empty, then we have exhausted the candidates
            }
            if (!x.empty()){ // could be empty if there are no more candidates
                total_num_launched += x.size() / num_dimensions;
                if (std::is_same<InitialGuessMode, AcceleratedMode>::value && (grid.getNumLoaded() > 0))
                    grid.evaluateBatch(x, y0);
                model(x, y, y0, 0); // compute a sample
                complete.add(x, y);
                manager.complete(x);

                // the fist thousand points can be loaded one at a time, then add when % increase of the grid is achieved
                if ((grid.getNumLoaded() < 1000) || (double(complete.getNumStored()) / double(grid.getNumLoaded()) > 0.2))
                    load_complete(); // also does checkpoint save
                // if done with the top % of the grid points, recompute the candidates
                if (double(manager.getNumDone()) / double(manager.getNumCandidates()) > 0.2)
                    refresh_candidates();
            }
            checkpoint();
        }
        complete.load(grid);
    }

    grid.finishConstruction();
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the model defined by the lambda.
 *
 * Creates a two way feedback between the model and the grid, samples from the model are collected
 * and loaded into the grid, while algorithms from the grid propose new samples according to
 * estimated importance.
 * The procedure is carried out until either tolerance is reached, the budget is exhausted,
 * or no more samples satisfy the level limits.
 * If set to parallel mode, the lambda will be called in separate threads concurrently
 * up to the specified maximum number.
 *
 * \tparam parallel_construction defines the use of parallel or sequential mode.
 *      The variable is of type \b bool but the constexpr constants
 *      TasGrid::mode_parallel and TasGrid::mode_sequential can be used to for more
 *      expressive calls.
 *
 * \param model defines the input-output relation to be approximated by the surrogate.
 *      In each call, \b x will have size equal to an even multiple of the dimension of
 *      the gird and will hold the required sample inputs for a set of points;
 *      the number of points is controlled by \b max_num_points
 *      The \b y will have size equal to the number of samples time
 *      the number of outputs and must be loaded the with corresponding outputs.
 *      If using the parallel mode, \b thread_id will be a number between 0 and
 *      \b max_num_samples \b -1, all threads running simultaneously will be
 *      given a different thread id.
 *      The \b y0 input will be the best guess for \b y given the current state of the grid,
 *      it is useful in cases where the model can benefit from a good initial guess,
 *      e.g., when using an iterative linear solver.
 *      If the model cannot use such information, it is safe to ignore the parameter
 *      or better yet, call the overload that does not accept \b y0 which will skip
 *      the calls to \b grid.evaluateBatch().
 *      If the grid does not have any loaded points (i.e., in the beginning of the
 *      construction process), the initial guess will be an empty vector.
 * \param max_num_points defines the computational budget for the surrogate construction.
 *      The construction procedure will terminate after the grid has reached
 *      the maximum number of points.
 * \param num_parallel_jobs defined the maximum number of concurent thread
 *      executing different calls to the \b model.
 *      In sequential mode, i.e., when \b parallel_construction is \b false,
 *      this number will loosely control the frequency of recalculating
 *      the list of candidate "most important" samples.
 *      If set to 0, it will be used as if set to 1.
 * \param max_samples_per_job defines the largest number of samples per call to \b model.
 *      In some cases, outputs can be more efficiently computed when the samples are
 *      lumped together. If the model evaluations are not oprimized for batching
 *      then this can be simply set to one.
 *      If set to 0, it will be used as if set to 1.
 * \param grid is the resulting surrogate model.
 *      The grid must be initialized with the appropriate type, number of dimensions,
 *      number of outputs, and sufficiently large number of initial points
 *      (e.g., there are enough candidates to load all threads).
 *      This template works with local polynomial grids.
 * \param tolerance defines the target tolerance, identical to
 *      TasmanianSparseGrid::setSurplusRefinement().
 * \param criteria defines the refinement algorithm, identical to
 *      TasmanianSparseGrid::setSurplusRefinement().
 * \param output defines the output to use in the algorithm, identical to
 *      TasmanianSparseGrid::setSurplusRefinement().
 * \param level_limits defines the maximum level to use in each direction,
 *      identical to TasmanianSparseGrid::setSurplusRefinement().
 *      If level limits are already set in the construction and/or refinement
 *      those weights will be used, this parameter can overwrite them.
 * \param scale_correction defines multiplicative correction term to the
 *      hierarchical coefficients used in the candidate selection procedure.
 *      Identical to TasmanianSparseGrid::setSurplusRefinement().
 * \param checkpoint_filename defines the two filenames to be used in to store the
 *      intermediate constructed grids so that the procedure can be restarted
 *      in case of a system crash.
 *      If the filename is "foo" then the files will be called "foo" and "foo_old".
 *      No intermediate saves will be made if the string is empty.
 *      If the string is not empty, the procedure will first attempt to recover
 *      from "foo" and "foo_old".
 *
 * \b WARNING: if the checkpoint files contain data from an older runs, the files must be deleted
 *             to avoid recovering from the old executing.
 *
 * \throws std::runtime_error if called for a grid that is not local polynomial.
 *
 *
 * \par Checkpoint Restart
 * Large scale simulations that take long time to complete, run a significant risk of
 * system failure which can waste significant computing resources.
 * Checkpoint-restart is a technique to periodically save the computed result
 * and, in case of a crash, restart from the last saved point.
 * For example, suppose the following call crashes mid-way:
 * \code
 * // original call
 * TasGrid::TasmanianSparseGrid grid = TasGrid::makeLocalPolynomialGrid(...);
 * constructSurrogate(model, budget, num_parallel, grid, 1.E-6, TasGrid::refine_classic, -1,
 *                    std::vector<int>(), std::vector<double>(), "foo");
 * \endcode
 * After the crash, Tasmanian will attempt to recover the computed samples from
 * files "foo" and "foo_old", if the files don't exist or do not contain
 * valid recovery data, the procedure will restart from scratch.
 * However, if the files were saved successfully the procedure will restart
 * mid-way and samples will not have to be recomputed.
 */
template<bool parallel_construction = TasGrid::mode_parallel>
void constructSurrogate(AcceleratedModelLambda model,
                        size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                        TasmanianSparseGrid &grid,
                        double tolerance, TypeRefinement criteria, int output = -1,
                        std::vector<int> const &level_limits = std::vector<int>(),
                        std::vector<double> const &scale_correction = std::vector<double>(),
                        std::string const &checkpoint_filename = std::string()){
    if (!grid.isLocalPolynomial()) throw std::runtime_error("ERROR: construction (with tolerance and criteria) called for a grid that is not local polynomial.");
    constructCommon<parallel_construction, AcceleratedMode>(model, max_num_points, num_parallel_jobs, max_samples_per_job, grid,
                                                            [&](TasmanianSparseGrid &g)->std::vector<double>{
                                                                return g.getCandidateConstructionPoints(tolerance, criteria, output, level_limits, scale_correction);
                                                            }, checkpoint_filename);
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Overload with model that cannot utilize an initial guess.
 */
template<bool parallel_construction = TasGrid::mode_parallel>
void constructSurrogate(SimpleModelLambda model,
                        size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                        TasmanianSparseGrid &grid,
                        double tolerance, TypeRefinement criteria, int output = -1,
                        std::vector<int> const &level_limits = std::vector<int>(),
                        std::vector<double> const &scale_correction = std::vector<double>(),
                        std::string const &checkpoint_filename = std::string()){
    if (!grid.isLocalPolynomial()) throw std::runtime_error("ERROR: construction (with tolerance and criteria) called for a grid that is not local polynomial.");
    constructCommon<parallel_construction, StaticMode>(modelNoInit(model), max_num_points, num_parallel_jobs, max_samples_per_job, grid,
                                                       [&](TasmanianSparseGrid &g)->std::vector<double>{
                                                           return g.getCandidateConstructionPoints(tolerance, criteria, output, level_limits, scale_correction);
                                                       }, checkpoint_filename);
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the model defined by the lambda.
 *
 * Uses the user provided \b anisotropic_weights to order the samples by importance
 * and calles the anisotropic overload of TasmanianSparseGrid::getCandidateConstructionPoints().
 * Otherwise the function is identical to TasGrid::constructSurrogate().
 *
 * \b WARNING: anisotropic refinement does not target a tolerance,
 * thus sampling will continue until the budget is exhausted
 * or the level limits are reached (which will produce a full tensor grid).
 */
template<bool parallel_construction = TasGrid::mode_parallel>
void constructSurrogate(AcceleratedModelLambda model,
                        size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, std::vector<int> const &anisotropic_weights = std::vector<int>(),
                        std::vector<int> const &level_limits = std::vector<int>(),
                        std::string const &checkpoint_filename = std::string()){
    constructCommon<parallel_construction, AcceleratedMode>(model, max_num_points, num_parallel_jobs, max_samples_per_job, grid,
                                                            [&](TasmanianSparseGrid &g)->std::vector<double>{
                                                                return g.getCandidateConstructionPoints(type, anisotropic_weights, level_limits);
                                                            }, checkpoint_filename);
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Overload with model that cannot utilize an initial guess.
 */
template<bool parallel_construction = TasGrid::mode_parallel>
void constructSurrogate(SimpleModelLambda model,
                        size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, std::vector<int> const &anisotropic_weights = std::vector<int>(),
                        std::vector<int> const &level_limits = std::vector<int>(),
                        std::string const &checkpoint_filename = std::string()){
    constructCommon<parallel_construction, StaticMode>(modelNoInit(model), max_num_points, num_parallel_jobs, max_samples_per_job, grid,
                                                       [&](TasmanianSparseGrid &g)->std::vector<double>{
                                                           return g.getCandidateConstructionPoints(type, anisotropic_weights, level_limits);
                                                       }, checkpoint_filename);
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the model defined by the lambda.
 *
 * Uses anisotropic weights to order the samples by importance,
 * starts with a fully isotropic grid until enough points are loaded to allow to estimate the weightss.
 * The procedure uses the anisotropic overload of TasmanianSparseGrid::getCandidateConstructionPoints(),
 * otherwise the function is identical to TasGrid::constructSurrogate().
 *
 * \b WARNING: anisotropic refinement does not target a tolerance,
 * thus sampling will continue until the budget is exhausted
 * or the level limits are reached (which will produce a full tensor grid).
 */
template<bool parallel_construction = TasGrid::mode_parallel>
void constructSurrogate(AcceleratedModelLambda model,
                        size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, int output, std::vector<int> const &level_limits = std::vector<int>(),
                        std::string const &checkpoint_filename = std::string()){
    constructCommon<parallel_construction, AcceleratedMode>(model, max_num_points, num_parallel_jobs, max_samples_per_job, grid,
                                                            [&](TasmanianSparseGrid &g)->std::vector<double>{
                                                                return g.getCandidateConstructionPoints(type, output, level_limits);
                                                            }, checkpoint_filename);
}



/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Overload with model that cannot utilize an initial guess.
 */
template<bool parallel_construction = TasGrid::mode_parallel>
void constructSurrogate(SimpleModelLambda model,
                        size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, int output, std::vector<int> const &level_limits = std::vector<int>(),
                        std::string const &checkpoint_filename = std::string()){
    constructCommon<parallel_construction, StaticMode>(modelNoInit(model), max_num_points, num_parallel_jobs, max_samples_per_job, grid,
                                                       [&](TasmanianSparseGrid &g)->std::vector<double>{
                                                           return g.getCandidateConstructionPoints(type, output, level_limits);
                                                       }, checkpoint_filename);
}

}

#endif
