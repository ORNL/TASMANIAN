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
 * \brief Signature of a model function to be used in the construction procedures.
 *
 * This types describes an abstract model with inputs, outputs and a specific thread-id.
 *
 * \param x is the model inputs, the vector is logically organized in strips each containing a set of model inputs.
 * \param y is the model outputs, the vector must be logically organized in strips each containing a set of model
 *      outputs corresponding to a set of inputs (the order must match, first set of inputs to first set outputs and so on).
 *      In most cases, \b y will have the correct size before the call, the exception is when there is insufficient data
 *      to compute an initial guess in a context that asks for an initial guess, then \b y will be empty and must be resized.
 * \param thread_id is a unique identifier of the thread associated with this model call, useful in multi-threaded context
 *      e.g., to assign a GPU accelerator to the current call.
 *
 * The call to TasGrid::constructSurrogate() will control the maximum number of inputs in a single model call,
 * the total number of threads, and whether \b y will have the correct size and/or contain an initial guess.
 */
using ModelSignature = std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t thread_id)>;

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
 * \ingroup TasmanianAddonsConstruct
 * \brief Allows for expressive calls to TasGrid::constructSurrogate().
 */
constexpr bool with_initial_guess = true;

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Allows for expressive calls to TasGrid::constructSurrogate().
 */
constexpr bool no_initial_guess = false;

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
template<bool parallel_construction, bool use_initial_guess>
void constructCommon(ModelSignature model,
                     size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                     TasmanianSparseGrid &grid,
                     std::function<std::vector<double>(TasmanianSparseGrid &)> candidates,
                     std::string const &checkpoint_filename){

    num_parallel_jobs   = std::max(size_t(1), num_parallel_jobs);
    max_samples_per_job = std::max(size_t(1), max_samples_per_job);
    size_t num_dimensions = (size_t) grid.getNumDimensions();
    size_t num_outputs    = (size_t) grid.getNumOutputs();
    CandidateManager manager(num_dimensions, max_samples_per_job); // keeps track of started and ordered samples
    CompleteStorage complete(num_dimensions); // temporarily stores complete samples (batch loading is faster)

    std::string filename = checkpoint_filename;
    std::string filename_old = checkpoint_filename + "_old";

    if (!filename.empty()){ // recover from an existing checkpoint
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
                // nothing could be recovered, start over from the current grid
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

    // load the initial guess into y (is using initial guess), otherwise set y to the correct size
    auto set_initial_guess = [&](std::vector<double> const &x, std::vector<double> &y)->void{
        if (use_initial_guess){
            if (grid.getNumLoaded()) grid.evaluateBatch(x, y);
            else y.clear();
        }else{
            y.resize(num_outputs * (x.size() / num_dimensions));
        }
    };

    if (!grid.isUsingConstruction()) // the procedure assumes dynamic construction
        grid.beginConstruction();

    refresh_candidates();

    if (parallel_construction == mode_parallel){ // parallel version
        // allocate space for all x and y pairs, will be filled by workers and processed by main
        std::vector<std::vector<double>> x(num_parallel_jobs),
                                         y(num_parallel_jobs, std::vector<double>(max_samples_per_job * num_outputs));

        std::vector<int> work_flag(num_parallel_jobs);
        constexpr int flag_done = 0;
        constexpr int flag_computing = 1;
        constexpr int flag_shutdown = 2;

        std::condition_variable until_someone_done;
        std::condition_variable until_new_job;
        std::mutex access_count_done;
        int count_done = 0;

        // lambda that will handle the work
        auto do_work = [&](size_t thread_id)->void{

            int my_flag = flag_computing;
            while(my_flag == flag_computing){
                model(x[thread_id], y[thread_id], thread_id); // does the model evaluations

                { // must guarantee sync between work_flag and count_done, use a lock
                    std::lock_guard<std::mutex> lock(access_count_done);
                    work_flag[thread_id] = flag_done;
                    count_done++;
                }
                until_someone_done.notify_one(); // just finished some work, notify the main thread

                { // wait till the main thread gives us an new piece of work
                    std::unique_lock<std::mutex> lock(access_count_done);
                    until_new_job.wait(lock, [&]()->bool{ return (work_flag[thread_id] != flag_done); });
                    my_flag = work_flag[thread_id];
                }
            }
        };

        // launch initial set of jobs
        std::vector<std::thread> workers(num_parallel_jobs);
        for(size_t id=0; id<num_parallel_jobs; id++){
            x[id] = manager.next(max_num_points - total_num_launched);
            if (!x[id].empty()){
                total_num_launched += x[id].size() / num_dimensions;
                set_initial_guess(x[id], y[id]);
                work_flag[id] = flag_computing;
                workers[id] = std::thread(do_work, id);
            }else{
                work_flag[id] = flag_shutdown; // not enough samples, cancel the thread
            }
        }

        auto collect_finished = [&]()->bool{
            bool any_done = false;
            for(size_t id=0; id<num_parallel_jobs; id++){
                if (work_flag[id] == flag_done){
                    if (!x.empty()){ // shouldn't be empty
                        complete.add(x[id], y[id]);
                        manager.complete(x[id]);
                        any_done = true;
                    }
                    if ((grid.getNumLoaded() < 1000) || (double(complete.getNumStored()) / double(grid.getNumLoaded()) > 0.2))
                        load_complete(); // move from complete into the grid

                    if (total_num_launched < max_num_points){
                        // refresh the candidates if enough of the current candidates have completed
                        if (double(manager.getNumDone()) / double(manager.getNumCandidates()) > 0.2)
                            refresh_candidates();

                        x[id] = checkout_sample(); // if necessary this will call refresh_candidates()

                        if (!x[id].empty()){ // if empty, then we have finished all possible candidates (reached tolerance)
                            total_num_launched += x[id].size() / num_dimensions;
                            set_initial_guess(x[id], y[id]);
                            work_flag[id] = flag_computing;
                        }else{
                            work_flag[id] = flag_shutdown; // not enough samples, cancel the thread
                        }
                    }else{
                        work_flag[id] = flag_shutdown; // reached the budget, shutdown the thread
                    }
                }
            }
            return any_done;
        };

        while(manager.getNumRunning() > 0){ // main loop
            {   // lock access to the count_done variable
                std::unique_lock<std::mutex> lock(access_count_done);
                // unlock and wait until some else increments the "done" count
                until_someone_done.wait(lock, [&]()->bool{ return (count_done > 0); });
                // the lock is back on at this point, process the completed samples, reset the count and go back to waiting
                count_done = 0;
                if (collect_finished()) checkpoint(); // if new samples were computed, save the state
            } // unlock the access_count_done and notify that we have loaded new jobs
            // without the unlock, the threads will wake up but will not be able to read the worker flags
            until_new_job.notify_all();
        }

        load_complete(); // flush completed jobs

        for(auto &w : workers) if (w.joinable()) w.join(); // join all threads

    }else{
        std::vector<double> x(grid.getNumDimensions()), y( grid.getNumOutputs());

        while((total_num_launched < max_num_points) && (manager.getNumCandidates() > 0)){
            x = manager.next(max_num_points - total_num_launched);
            if (x.empty()){ // need more candidates
                refresh_candidates();
                x = manager.next(max_num_points - total_num_launched); // if this is empty, then we have exhausted the candidates
            }
            if (!x.empty()){ // could be empty if there are no more candidates
                total_num_launched += x.size() / num_dimensions;
                set_initial_guess(x, y);
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
            checkpoint();
        }

        load_complete(); // flush completed jobs
    }
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
 * Two notable options are the ability to call the model with batches of samples and
 * the ability to assign an initial guess for each sample.
 * - In some cases, model evaluations can be performed much more efficiently using batches of
 * samples, e.g., by allowing the usage of BLAS level 3 methods (as opposed to level 2),
 * better occupancy on GPU accelerators, or better memory cohesion in sparse linear algebra.
 * For more details, see the work by:\n
 * E. Phipps, M. D'Elia, H. Edwards, M. Hoemmen, J. Hu,
 * <a style="font-weight:bold" href="https://epubs.siam.org/doi/abs/10.1137/15M1044679">Embedded Ensemble Propagation for Improving
 * Performance, Portability, and Scalability of Uncertainty Quantification on Emerging Computational Architectures</a>,
 * SIAM Journal on Scientific Computing, vol. 39, num. 2, pp. C162--C193, 2017.\n
 * M. D'Elia, E. Phipps, A. Rushdi, M. Ebeida,
 * <a style="font-weight:bold" href="https://arxiv.org/abs/1705.02003">Surrogate-based Ensemble Grouping Strategies for
 * Embedded Sampling-based Uncertainty Quantification</a>, arXiv:1705.02003.
 *
 * - Some models can utilize (or even require) a good initial guess to perform a simulation,
 * e.g., when using linear or non-linear iterative solvers.
 * As the sparse grid surrogate becomes closer and closer to the model, increasingly better guess can be
 * computes, effectively relying on the model to solve only for the correction between the current
 * sparse grid approximation and the actual model. For more details and examples see the work of:\n
 * D. Galindo, P. Jantsch, C. G. Webster, and G. Zhang,
 * <a style="font-weight:bold" href="https://doi.org/10.1137/15M1019568">Accelerating Stochastic Collocation Methods for Partial Differential
 * Equations with Random Input Data</a>, SIAM/ASA Journal on Uncertainty Quantification, vol. 4, num. 1, pp. 1111--1137, 2016.
 *
 * The template can be instantiated in either parallel or single threaded mode, with or without an initial guess feedback,
 * and sampling can be performed in batches or single point.
 * The rest of the parameters control the computational budget and the specifics of the refinement scheme.
 *
 * \tparam parallel_construction defines the use of parallel or sequential mode.
 *      The variable is of type \b bool but the constexpr constants
 *      TasGrid::mode_parallel and TasGrid::mode_sequential can be used to for more
 *      expressive calls.
 * \tparam initial_guess defines the state of the output vector \b y when the model is called.
 *      If the initial guess is set to \b true (or TasGrid::with_initial_guess),
 *      then the input will be filled with the best guess for \b y based on the current grid,
 *      i.e., the previously computed samples.
 *      This mode is useful when the model can benefit from a good initial guess,
 *      e.g., when using an iterative solver.
 *      If the model cannot use an initial guess, then use \b false
 *      or TasGrid::no_initial_guess to avoid extraneous calls to \b grid.evaluateBatch().
 *
 *
 * \param model defines the input-output relation to be approximated by the surrogate.
 *      - Entering each call, \b x will have size equal to an even multiple of the dimension of
 *      the gird and will hold the required sample inputs for a set of points;
 *      the number of points is controlled by \b max_samples_per_job.
 *      - Exiting the call, \b y must have size equal to the number of samples time
 *      the number of outputs and must be loaded the outputs corresponding to \b x.
 *      If running with \b initial_guess, on entry, \b y will be either loaded with
 *      the best guess based on the previously computed samples or empty which
 *      will indicate that there are not enough samples to make even a coarse grid.
 *      If running without \b initial_guess, the vector will always have the correct
 *      size, but the values will be unspecified and should be overwritten.
 *      - If using the parallel mode, \b thread_id will be a number between 0 and
 *      \b max_num_samples \b -1, all threads running simultaneously will be
 *      given a different thread id.
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
 * \par Example
 * Constructing a surrogate to a simple exponential function with budget of 100 points
 * and running on 4 threads with one sample per thread:
 * \code
 *   auto grid = TasGrid::makeLocalPolynomialGrid(2, ..); // sync the dimensions with number of model inputs
 *   TasGrid::constructSurrogate([&](std::vector<double> const &x, std::vector<double> &y)
 *                               ->void{ y[0] = std::exp(x[0] + x[1]); },
 *                               100, 4, 1, grid, 1.E-6, TasGrid::refine_classic);
 * \endcode
 * The procedure will terminate after 100 samples or after the tolerance of 1.E-6 has been reached.
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
template<bool parallel_construction = TasGrid::mode_parallel, bool initial_guess = no_initial_guess>
void constructSurrogate(ModelSignature model,
                        size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                        TasmanianSparseGrid &grid,
                        double tolerance, TypeRefinement criteria, int output = -1,
                        std::vector<int> const &level_limits = std::vector<int>(),
                        std::string const &checkpoint_filename = std::string()){
    if (!grid.isLocalPolynomial() && !grid.isWavelet()) throw std::runtime_error("ERROR: construction (with tolerance and criteria) called for a grid that is not local polynomial or wavelet.");
    constructCommon<parallel_construction, initial_guess>
                                          (model, max_num_points, num_parallel_jobs, max_samples_per_job, grid,
                                           [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(tolerance, criteria, output, level_limits);
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
template<bool parallel_construction = TasGrid::mode_parallel, bool initial_guess = no_initial_guess>
void constructSurrogate(ModelSignature model,
                        size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, std::vector<int> const &anisotropic_weights = std::vector<int>(),
                        std::vector<int> const &level_limits = std::vector<int>(),
                        std::string const &checkpoint_filename = std::string()){
    constructCommon<parallel_construction, initial_guess>
                                           (model, max_num_points, num_parallel_jobs, max_samples_per_job, grid,
                                           [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(type, anisotropic_weights, level_limits);
                                           }, checkpoint_filename);
}

/*!
 * \ingroup TasmanianAddonsConstruct
 * \brief Construct a sparse grid surrogate to the model defined by the lambda.
 *
 * Uses anisotropic weights to order the samples by importance,
 * starts with a fully isotropic grid until enough points are loaded to allow to estimate the weights.
 * The procedure uses the anisotropic overload of TasmanianSparseGrid::getCandidateConstructionPoints(),
 * otherwise the function is identical to TasGrid::constructSurrogate().
 *
 * \b WARNING: anisotropic refinement does not target a tolerance,
 * thus sampling will continue until the budget is exhausted
 * or the level limits are reached (which will produce a full tensor grid).
 */
template<bool parallel_construction = TasGrid::mode_parallel, bool initial_guess = no_initial_guess>
void constructSurrogate(ModelSignature model,
                        size_t max_num_points, size_t num_parallel_jobs, size_t max_samples_per_job,
                        TasmanianSparseGrid &grid,
                        TypeDepth type, int output, std::vector<int> const &level_limits = std::vector<int>(),
                        std::string const &checkpoint_filename = std::string()){
    constructCommon<parallel_construction, initial_guess>
                                           (model, max_num_points, num_parallel_jobs, max_samples_per_job, grid,
                                           [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(type, output, level_limits);
                                           }, checkpoint_filename);
}

}

#endif
