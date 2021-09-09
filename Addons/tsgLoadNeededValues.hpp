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

#ifndef __TASMANIAN_ADDONS_LOADNEEDEDVALS_HPP
#define __TASMANIAN_ADDONS_LOADNEEDEDVALS_HPP

/*!
 * \internal
 * \file tsgLoadNeededValues.hpp
 * \brief Templates for non-adaptive sampling from lambda models.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianAddonsCommon
 *
 * Templates that call call a model at the needed points and
 * load the corresponding model values.
 * \endinternal
 */

#include "tsgAddonsCommon.hpp"

/*!
 * \ingroup TasmanianAddons
 * \addtogroup TasmanianAddonsLoadNeededVals Static Load Model Values
 *
 * Handy templates to automate the common use case of evaluating
 * a model at the needed (or loaded) points and feeding the values to the grid.
 */

namespace TasGrid{

/*!
 * \ingroup TasmanianAddonsLoadNeededVals
 * \brief Loads the current grid with model values, does not perform any refinement.
 *
 * Loads model values into the grid using the provided lambda function,
 * the model can be computed either sequentially or in parallel.
 * This is a non-adaptive procedure, i.e., the points will not be changes
 * only the associated model values will be modified.
 *
 * \tparam parallel_construction defines whether to run in parallel or sequential mode.
 * \tparam overwrite_loaded defines whether to overwrite the currently loaded model
 *                          values or to assign values to the needed points.
 *
 * \param model is the lambda representation of a model with inputs \b x and outputs \b y,
 *              the lambda must accept arrays with size equal to the grid dimensions and
 *              outputs. For each \b x, the \b y must be overwritten with the corresponding
 *              model values. If parallel sampling is uses, then \b thread_id will be
 *              a unique number between 0 and \b num_threads-1 associated with the running
 *              thread. The id can help associate each thread with additional external resources,
 *              e.g., a separate CUDA device can be associated with each thread.
 *
 * \param grid is the sparse grid that will be loaded. The grid must not be set for construction
 *             and the number of inputs must be positive. The method grid.loadNeededValues()
 *             will be called with values corresponding to the model outputs at either
 *             the loaded or the needed points.
 *
 * \param num_threads is the number of parallel calls to the \b model lambda.
 *                    If set to zero, sequential mode will be used without launching any threads,
 *                    the number is ignored in sequential mode.
 *                    Note that this is using the C++ native std::thread as opposed to OpenMP.
 *
 * \throws std::runtime_error if grid.isUsingConstruction() is true or if grid.getNumOutputs() is zero.
 *
 * Example:
 * \code
 * // construct a grid for 10-th order global polynomials
 * auto grid = TasGrid::makeGlobalGrid(4, 1, 10, TasGrid::type_iptotal, TasGrid::rule_leja);
 * // note that the grid is set to 4 inputs and 1 output, the model function should be the same
 * auto model = [](double const x[], double y[], size_t)->void{
 *      y[0] = std::exp(x[0] + x[1] + x[2] + x[3]);
 * }
 * loadNeededValues(model, grid, 4); // using parallel sampling with 4 threads
 * // at this point, grid is a 10-th order polynomial approximation to exp(x0 + x2 + x3 + x4)
 * \endcode
 */
template<bool parallel_construction = true, bool overwrite_loaded = false>
void loadNeededValues(std::function<void(double const x[], double y[], size_t thread_id)> model, TasmanianSparseGrid &grid, size_t num_threads){
    int num_points = (overwrite_loaded) ? grid.getNumLoaded() : grid.getNumNeeded();
    int num_outputs = grid.getNumOutputs();
    if (grid.isUsingConstruction()) throw std::runtime_error("ERROR: cannot call loadNeededPoints() addon when isUsingConstruction() is true");
    if (num_outputs == 0) throw std::runtime_error("ERROR: cannot call loadNeededPoints() addon when the grid has no outputs");
    if (num_points == 0) return; // nothing to do here
    if (overwrite_loaded && (grid.getNumNeeded() != 0)) grid.clearRefinement(); // using loaded points only, clear the refinement

    // get the points and values
    auto points = (overwrite_loaded) ? grid.getLoadedPoints() : grid.getNeededPoints();
    std::vector<double> values(Utils::size_mult(num_points, num_outputs));

    // divide logically into strips
    Utils::Wrapper2D<double> xwrap(grid.getNumDimensions(), points.data());
    Utils::Wrapper2D<double> ywrap(num_outputs, values.data());

    if (parallel_construction && (num_threads > 0)){
        std::vector<bool> checked_out(num_points, false);
        std::mutex checked_out_lock;

        std::vector<std::thread> workers;
        workers.reserve(num_threads);
        for(size_t thread_id=0; thread_id<num_threads; thread_id++){
            workers.emplace_back( // create a new worker thread
                [&, thread_id](void)->void{
                    int sample = 0;
                    do{
                        { // find the next sample
                            std::lock_guard<std::mutex> lock(checked_out_lock);
                            while ((sample < num_points) && checked_out[sample]) sample++;
                            if (sample < num_points) checked_out[sample] = true;
                        }
                        if (sample < num_points) // if found, compute the next sample
                            model(xwrap.getStrip(sample), ywrap.getStrip(sample), thread_id);
                    }while(sample < num_points);
                }
            );
        }

        for(auto &w : workers) w.join(); // wait till finished

    }else{
        for(int i=0; i<num_points; i++)
            model(xwrap.getStrip(i), ywrap.getStrip(i), 0);
    }

    grid.loadNeededPoints(values);
}

/*!
 * \ingroup TasmanianAddonsLoadNeededVals
 * \brief Overload that uses vectors for the model inputs and outputs.
 *
 * The template will accept the array inputs from the main implementation
 * and copy those to the vectors.
 * Thus, there is an extra copy of data, which is unavoidable since
 * TasGrid::TasmanianSparseGrid returns contiguous vectors that cannot be split
 * without copy.
 */
template<bool parallel_construction = true, bool overwrite_loaded = false>
void loadNeededValues(std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t thread_id)> model,
                      TasmanianSparseGrid &grid, size_t num_threads){
    int num_dimensions = grid.getNumDimensions();
    int num_outputs = grid.getNumOutputs();
    loadNeededValues<parallel_construction, overwrite_loaded>(
        [&](double const x[], double y[], size_t thread_id)->void{
            std::vector<double> vecy(num_outputs);
            model(std::vector<double>(x, x + num_dimensions), vecy, thread_id);
            std::copy(vecy.begin(), vecy.end(), y);
        }, grid, num_threads);
}
/*!
 * \ingroup TasmanianAddonsLoadNeededVals
 * \brief Alias to loadNeededValues(), array variant.
 */
template<bool parallel_construction = true, bool overwrite_loaded = false>
void loadNeededPoints(std::function<void(double const x[], double y[], size_t thread_id)> model, TasmanianSparseGrid &grid, size_t num_threads){
    loadNeededValues<parallel_construction, overwrite_loaded>(model, grid, num_threads);
}
/*!
 * \ingroup TasmanianAddonsLoadNeededVals
 * \brief Alias to loadNeededValues(), vector variant.
 */
template<bool parallel_construction = true, bool overwrite_loaded = false>
void loadNeededPoints(std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t thread_id)> model,
                      TasmanianSparseGrid &grid, size_t num_threads){
    loadNeededValues<parallel_construction, overwrite_loaded>(model, grid, num_threads);
}

}

#endif
