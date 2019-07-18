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

#include "TasmanianAddons.hpp"
#include "tasgridCLICommon.hpp"

using std::cout;
using std::setw;

//! \brief Runs tests and throws if the surrogates do not match to the given tolerance.
inline void compareGrids(double tolerance, TasGrid::TasmanianSparseGrid const &a, TasGrid::TasmanianSparseGrid const &b){
    if (a.getNumDimensions() != b.getNumDimensions())
        throw std::runtime_error("grids have wrong dimensions");
    if (a.getNumOutputs() != b.getNumOutputs())
        throw std::runtime_error("grids have wrong outputs");
    if (a.getNumPoints() != b.getNumPoints())
        throw std::runtime_error("grids have wrong number of points");

    std::minstd_rand park_miller(77);
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::vector<double> test_points(2 * 1000);
    for(auto &t : test_points) t = unif(park_miller);

    std::vector<double> resa, resb; // reference and actual result
    a.evaluateBatch(test_points, resa);
    b.evaluateBatch(test_points, resb);
    double err = 0.0;
    for(auto ia = resa.begin(), ib = resb.begin(); ia != resa.end(); ia++, ib++)
        err = std::max(err, std::abs(*ia - *ib));

    if (err > tolerance)
        throw std::runtime_error("grids have outputs that differ by more than the tolerance");
}

//! \brief Simple sequential construction for reference purposes.
void simpleSequentialConstruction(std::function<void(std::vector<double> const &x, std::vector<double> &y)> model,
                                  size_t max_num_samples, size_t num_parallel_jobs,
                                  TasmanianSparseGrid &grid,
                                  double tolerance, TypeRefinement criteria, int output = -1,
                                  std::vector<int> const &level_limits = std::vector<int>(),
                                  std::vector<double> const &scale_correction = std::vector<double>()){
    if (!grid.isUsingConstruction())
        grid.beginConstruction();

    auto candidates = grid.getCandidateConstructionPoints(tolerance, criteria, output, level_limits, scale_correction);
    while(!candidates.empty() && (grid.getNumLoaded() < (int) max_num_samples)){
        CompleteStorage storage(grid.getNumDimensions());
        size_t num_samples = std::min(num_parallel_jobs, max_num_samples - grid.getNumLoaded());
        for(size_t i=0; i<num_samples; i++){
            std::vector<double> x(&candidates[i*grid.getNumDimensions()], &candidates[i*grid.getNumDimensions()] + grid.getNumDimensions());
            std::vector<double> y(grid.getNumOutputs());
            model(x, y);
            storage.add(x, y);
        }
        storage.load(grid);
        auto pnts = grid.getLoadedPoints();
        candidates = grid.getCandidateConstructionPoints(tolerance, criteria, output, level_limits, scale_correction);
    }

    grid.finishConstruction();
}

//! \brief Tests the templates for automated construction.
bool testConstructSurrogate(){
    auto model_exp = [&](std::vector<double> const &x, std::vector<double> &y)->void{ y = {std::exp(x[0] + x[1])}; };
    auto grid = TasGrid::makeLocalPolynomialGrid(2, 1, 3, 2);
    auto reference_grid = grid;

    TasGrid::constructSurrogate(model_exp, -1, 2, grid, 1.E-4, TasGrid::refine_classic);
    simpleSequentialConstruction(model_exp, 700, 2, reference_grid, 1.E-4, TasGrid::refine_classic);
    compareGrids(1.E-9, grid, reference_grid);

    return true;
}
