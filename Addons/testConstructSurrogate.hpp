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
#include <limits>

using std::cout;
using std::setw;

//! \brief Runs tests and throws if the surrogates do not match to the given tolerance.
inline void compareGrids(double tolerance, TasGrid::TasmanianSparseGrid const &a, TasGrid::TasmanianSparseGrid const &b, bool check_num_points){
    if (a.getNumDimensions() != b.getNumDimensions())
        throw std::runtime_error("grids have wrong dimensions");
    if (a.getNumOutputs() != b.getNumOutputs())
        throw std::runtime_error("grids have wrong outputs");
    if (check_num_points){ // some tests have a random number of points
        if (a.getNumPoints() != b.getNumPoints()){
            cout << "number of points = " << a.getNumPoints() << "  " << b.getNumPoints() << endl;
            throw std::runtime_error("grids have wrong number of points");
        }
    }

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

    if (err > tolerance){
        cout << "difference = " << err << " expected = " << tolerance << endl;
        throw std::runtime_error("grids have outputs that differ by more than the tolerance");
    }
}

//! \brief Simple sequential construction for reference purposes.
void simpleSequentialConstruction(std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t)> model,
                                  size_t max_num_samples, size_t num_parallel_jobs,
                                  TasmanianSparseGrid &grid,
                                  std::function<std::vector<double>(TasGrid::TasmanianSparseGrid&)> get_candidates){
    if (!grid.isUsingConstruction())
        grid.beginConstruction();

    auto candidates = get_candidates(grid);

    while(!candidates.empty() && (grid.getNumLoaded() < (int) max_num_samples)){
        CompleteStorage storage(grid.getNumDimensions());
        size_t num_samples = std::min(std::min(num_parallel_jobs, max_num_samples - grid.getNumLoaded()), candidates.size() / grid.getNumDimensions());
        for(size_t i=0; i<num_samples; i++){
            std::vector<double> x(&candidates[i*grid.getNumDimensions()], &candidates[i*grid.getNumDimensions()] + grid.getNumDimensions());
            std::vector<double> y(grid.getNumOutputs());
            model(x, y, i);
            storage.add(x, y);
        }
        storage.load(grid);
        auto pnts = grid.getLoadedPoints();
        candidates = get_candidates(grid);
    }

    grid.finishConstruction();
}

//! \brief Reference construction for localp grids.
void simpleSequentialConstruction(std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t)> model,
                                  size_t max_num_samples, size_t num_parallel_jobs,
                                  TasmanianSparseGrid &grid,
                                  double tolerance, TypeRefinement criteria, int output = -1,
                                  std::vector<int> const &level_limits = std::vector<int>(),
                                  std::vector<double> const &scale_correction = std::vector<double>()){
    simpleSequentialConstruction(model, max_num_samples, num_parallel_jobs, grid,
                                 [&](TasGrid::TasmanianSparseGrid& g)->std::vector<double>{
                                     return g.getCandidateConstructionPoints(tolerance, criteria, output, level_limits, scale_correction);
                                });
}

//! \brief Reference construction for anisotropic refinement.
void simpleSequentialConstruction(std::function<void(std::vector<double> const &x, std::vector<double> &y, size_t)> model,
                                  size_t max_num_samples, size_t num_parallel_jobs,
                                  TasmanianSparseGrid &grid,
                                  TypeDepth type, int output, std::vector<int> const &level_limits = std::vector<int>()){
    simpleSequentialConstruction(model, max_num_samples, num_parallel_jobs, grid,
                                 [&](TasGrid::TasmanianSparseGrid& g)->std::vector<double>{
                                     return g.getCandidateConstructionPoints(type, output, level_limits);
                                });
}

//! \brief Tests the templates for automated construction.
bool testConstructSurrogate(bool verbose){
    std::atomic_int last;
    last = -1;
    constexpr unsigned int delay_on_lock = 2;
    // test the simple loadNeededPoints()
    auto model_trig = [&](double const x[], double y[], size_t)->void{ y[0] = std::sin(x[0]) * std::cos(x[1]); };
    auto model_trig_vec = [&](std::vector<double> const &x, std::vector<double> &y, size_t id)->void{ model_trig(x.data(), y.data(), id); };
    auto grid = TasGrid::makeSequenceGrid(2, 1, 9, TasGrid::type_level, TasGrid::rule_leja);
    auto reference_grid = grid;
    auto points = reference_grid.getPoints();
    std::vector<double> values(reference_grid.getNumPoints());
    for(size_t i=0; i<values.size(); i++)
        model_trig(&points[2*i], &values[i], 0);
    reference_grid.loadNeededPoints(values);

    loadNeededPoints<false, false>(model_trig, grid, 128); // sequential version ignores the num_threads
    compareGrids(1.E-10, grid, reference_grid, true);

    grid.loadNeededPoints(std::vector<double>(reference_grid.getNumPoints(), -1.0)); // set wrong values
    loadNeededPoints<true, true>(model_trig, grid, 4); // reload the needed points using 4 threads
    compareGrids(1.E-10, grid, reference_grid, true);

    grid = TasGrid::makeSequenceGrid(2, 1, 9, TasGrid::type_level, TasGrid::rule_leja);
    loadNeededPoints(model_trig_vec, grid, 4);
    compareGrids(1.E-10, grid, reference_grid, true);

    if (verbose) cout << std::setw(40) << "simple load values" << std::setw(10) << "Pass" << endl;

    // parallel construction is susceptible to order of execution, number of points and which points may change from one run to the next
    auto model_exp_seq = [&](std::vector<double> const &x, std::vector<double> &y, size_t)->void{ y = {std::exp(x[0] + x[1])}; };
    auto model_exp = [&](std::vector<double> const &x, std::vector<double> &y, size_t id)->void{
        model_exp_seq(x, y, id);
        if (last == (int) id) std::this_thread::sleep_for(std::chrono::milliseconds(delay_on_lock));
        else last = (int) id;
    };
    auto model_exp2 = [&](std::vector<double> const &x, std::vector<double> &y, size_t id)->void{
        if (x.size() == 2) y = {std::exp(x[0] + x[1])}; // one sample
        else y = {std::exp(x[0] + x[1]), std::exp(x[2] + x[3])}; // two samples
        if (last == (int) id) std::this_thread::sleep_for(std::chrono::milliseconds(delay_on_lock));
        else last = (int) id;
    }; // two samples
    grid = TasGrid::makeLocalPolynomialGrid(2, 1, 3, 2);
    reference_grid = grid;

    // Basic test, run until grid points are exhausted, number of points can still vary
    // due to child surplus being computed before or after the parent point has completed
    TasGrid::constructSurrogate<TasGrid::mode_parallel, no_initial_guess>
                               (model_exp, std::numeric_limits<size_t>::max(), 2, 1, grid, 1.E-4, TasGrid::refine_classic); // parallel
    simpleSequentialConstruction(model_exp, 300, 2, reference_grid, 1.E-4, TasGrid::refine_classic);
    compareGrids(5.E-4, grid, reference_grid, false);
    if (verbose) cout << std::setw(40) << "parallel localp unlimited budget" << std::setw(10) << "Pass" << endl;

    grid = TasGrid::makeLocalPolynomialGrid(2, 1, 3, 2);
    TasGrid::constructSurrogate<TasGrid::mode_parallel, no_initial_guess>
                               (model_exp2, std::numeric_limits<size_t>::max(), 2, 2, grid, 1.E-4, TasGrid::refine_classic); // parallel, batch 2
    compareGrids(5.E-4, grid, reference_grid, false);
    if (verbose) cout << std::setw(40) << "parallel localp batch" << std::setw(10) << "Pass" << endl;

    // Similar test, the number of points must match, but the actual points can be different
    // compare the grids by how similar they are to each other
    grid = TasGrid::makeLocalPolynomialGrid(2, 1, 3, 1);
    reference_grid = grid;
    TasGrid::constructSurrogate<TasGrid::mode_parallel, no_initial_guess>
                               (model_exp, 500, 4, 1, grid, 1.E-4, TasGrid::refine_classic); // parallel, limit points
    simpleSequentialConstruction(model_exp, 500, 4, reference_grid, 1.E-4, TasGrid::refine_classic);
    compareGrids(5.E-3, grid, reference_grid, true);
    if (verbose) cout << std::setw(40) << "parallel localp limited budget" << std::setw(10) << "Pass" << endl;

    // Sequential test, checks the simpler algorithm, this should be deterministic
    grid = TasGrid::makeLocalPolynomialGrid(2, 1, 3, 2);
    reference_grid = grid;
    TasGrid::constructSurrogate<TasGrid::mode_sequential, no_initial_guess>
                               (model_exp_seq, std::numeric_limits<size_t>::max(), 2, 1, grid, 1.E-4, TasGrid::refine_classic); // sequential
    simpleSequentialConstruction(model_exp, 300, 2, reference_grid, 1.E-4, TasGrid::refine_classic);
    compareGrids(1.E-9, grid, reference_grid, true);
    if (verbose) cout << std::setw(40) << "sequential localp limited budget" << std::setw(10) << "Pass" << endl;

    // The construction algorithm is the same, but check if the getCandidateConstructionPoints() lambda works right
    auto model_aniso_seq = [&](std::vector<double> const &x, std::vector<double> &y, size_t)->void{
        y = {std::exp(x[0] + 0.1 * x[1])};
    };
    auto model_aniso = [&](std::vector<double> const &x, std::vector<double> &y, size_t id)->void{
        model_aniso_seq(x, y, id);
        if (last == (int) id) std::this_thread::sleep_for(std::chrono::milliseconds(delay_on_lock));
        else last = (int) id;
    };
    grid = TasGrid::makeGlobalGrid(2, 1, 3, TasGrid::type_level, TasGrid::rule_rleja);
    reference_grid = grid;
    TasGrid::constructSurrogate<TasGrid::mode_parallel, no_initial_guess>
                               (model_aniso, 200, 3, 1, grid, TasGrid::type_iptotal, 0); // parallel
    simpleSequentialConstruction(model_aniso, 200, 2, reference_grid, TasGrid::type_iptotal, 0);
    compareGrids(1.E-8, grid, reference_grid, true);
    if (verbose) cout << std::setw(40) << "parallel anisotropic rleja" << std::setw(10) << "Pass" << endl;

    // additional fluctuation of number of points can happen due to not enough points to complete a tensor
    // nevertheless the grid accuracy should match reasonably well
    grid = TasGrid::makeGlobalGrid(2, 1, 3, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
    reference_grid = grid;
    TasGrid::constructSurrogate<TasGrid::mode_parallel, no_initial_guess>
                               (model_aniso, 400, 3, 1, grid, TasGrid::type_iptotal, 0); // parallel
    simpleSequentialConstruction(model_aniso, 400, 2, reference_grid, TasGrid::type_iptotal, 0);
    compareGrids(5.E-4, grid, reference_grid, false);
    if (verbose) cout << std::setw(40) << "parallel anisotropic global" << std::setw(10) << "Pass" << endl;

    // fix the weights, the computed grid must be very similar to the direct anisotropic make grid
    // i.e., the "most important" points are defined the same way through the user provided anisotropic weights
    // the atomic trick is used to ensure that no thread lags and holds important samples that affect the estimate for the anisotropy
    last = -1;
    std::vector<int> aweights = {1, 2};
    reference_grid = TasGrid::makeSequenceGrid(2, 1, 9, TasGrid::type_level, TasGrid::rule_leja, aweights);
    loadNeededPoints<mode_sequential>(model_aniso_seq, reference_grid, 0);

    grid = TasGrid::makeSequenceGrid(2, 1, 1, TasGrid::type_level, TasGrid::rule_leja);
    TasGrid::constructSurrogate<TasGrid::mode_parallel, no_initial_guess>
                               (model_aniso, reference_grid.getNumLoaded(), 4, 1, grid, TasGrid::type_iptotal, aweights); // parallel
    compareGrids(5.E-7, grid, reference_grid, true);
    if (verbose) cout << std::setw(40) << "parallel weighted sequence" << std::setw(10) << "Pass" << endl;

    // test checkpoint-restart
    std::atomic_int num_good;
    num_good = 0;
    auto model_crash = [&](std::vector<double> const &x, std::vector<double> &y, size_t)->void{
        y = {std::exp(x[0] + 0.1 * x[1])};
        num_good++;
        if (num_good == 10) throw std::runtime_error("test fail"); // simulate a crash on iteration 10
    };

    grid = TasGrid::makeSequenceGrid(2, 1, 1, TasGrid::type_level, TasGrid::rule_leja);
    std::remove("checkpoint"); // make sure the checkpoint filename is not present (e.g., left after an earlier crash)
    try{
        TasGrid::constructSurrogate<TasGrid::mode_sequential, no_initial_guess>
                                   (model_crash, reference_grid.getNumLoaded(),
                                    2, 1, grid, TasGrid::type_iptotal, aweights, {}, "checkpoint");
        std::cout << "ERROR: testConstructSurrogate() could not simulate a crash." << std::endl;
        return false;
    }catch(std::runtime_error &){
        //cout << "Error caught!" << endl;
        grid = TasGrid::makeSequenceGrid(2, 1, 1, TasGrid::type_level, TasGrid::rule_rleja); // recovery should change the rule
        TasGrid::constructSurrogate<TasGrid::mode_sequential, no_initial_guess>
                                   (model_crash, reference_grid.getNumLoaded(),
                                    2, 1, grid, TasGrid::type_iptotal, aweights, {}, "checkpoint");
        if (grid.getRule() != rule_leja) throw std::runtime_error("Failed recovery from crash, wrong rule.");
        compareGrids(1.E-9, grid, reference_grid, true);
        if (verbose) cout << std::setw(40) << "parallel resilient sequence" << std::setw(10) << "Pass" << endl;
    }
    if (std::remove("checkpoint") != 0) throw std::runtime_error("Could not delete the 'checkpoint' file, the file must exists after the test.");

    return true;
}
