
#include "Tasmanian.hpp"
#include <random>

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_06.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 6.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples6 Tasmanian Sparse Grids module, example 6
 *
 * \par Example 6
 * Build a surrogate model using automatic construction.
 */

/*!
 * \ingroup TasmanianSGExamples6
 * \brief Sparse Grids Example 6: adaptive surrogate modeling
 *
 * Example 5 demonstrates how to construct a surrogate by processing
 * the samples in batches. The batch construction can facilitate
 * parallel sampling, but the exact number of samples is hard to control
 * since the points from the entire batch are needed at once.
 * Furthermore, if a fixed number of computing nodes is available
 * (e.g., processor cores or MPI ranks) it may be hard to keep all
 * nodes loaded if the batch size does not divide evenly.
 *
 * The construction procedure offers a more flexible refinement alternative.
 * The batch size can be controlled at a much finer scale to accommodate
 * a fixed budget and occupancy across all parallel resources can be guaranteed
 * (unless the level limits or tolerance are about to be reached,
 * or the initial grid is too coarse).
 * See the documentation for TasGrid::mpiConstructSurrogate() for
 * the MPI variant of this method that differs only in some MPI related
 * parameters.
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_06.cpp SG_Example_06 example
 */
void sparse_grids_example_06(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_06 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Example 6: interpolate f(x,y) = exp(-x^2) * exp(y) * cos(z), using the rleja rule\n"
         << "           employs adaptive construction\n";

    // define the model as a C++ lambda expression, the advantage of lambdas
    // is that they can wrap around very complex models
    // the construction signature uses std::vector<double>
    int const num_inputs  = 3;
    int const num_outputs = 1;
    auto model = [&](std::vector<double> const &x, std::vector<double> &y, size_t)->
        void{
            // if no initial guess is used, then y will be guaranteed to have
            // the correct size, i.e., num_outputs
            // otherwise, y must be correctly resized
            y.resize(num_outputs);
            y[0] = std::exp(-x[0] * x[0]) * std::exp(x[1]) * std::cos(x[2]);
        };

    // the accuracy of the surrogate models is measure from 1000 random reference points
    int const num_test_points = 1000;
    std::vector<double> test_points(num_test_points * num_inputs);
    std::minstd_rand park_miller(42);
    std::uniform_real_distribution<double> domain(-1.0, 1.0);
    for(auto &t : test_points) t = domain(park_miller);

    std::vector<double> reference_values(num_test_points * num_outputs);
    for(int i=0; i<num_test_points; i++){
        std::vector<double> y(1);
        model(std::vector<double>(test_points.begin() + i * num_inputs,
                                  test_points.begin() + (i+1) * num_inputs),
              y, 0);
        std::copy(y.begin(), y.end(), reference_values.begin() + i * num_outputs);
    }

    // computes the maximum error between the reference values and the grid interpolant
    auto test_grid = [&](TasGrid::TasmanianSparseGrid const &grid)->
        double{
            std::vector<double> result;
            grid.evaluateBatch(test_points, result);
            double diff = 0.0;
            for(int i=0; i<num_test_points; i++)
                diff = std::max(diff, std::abs(result[i] - reference_values[i]));
            return diff;
        };

    // when the model is cheap to compute, running multiple-threads results in a large
    // difference in the time-per-sample due to the dominating cost of thread
    // synchronization
    // unfortunately, this can have adverse effect when important samples are delayed
    // and less important samples are computed thus wasting the budget
    // for that reason, the example runs sequentially
    constexpr auto mode = TasGrid::mode_sequential;

    auto grid = TasGrid::makeSequenceGrid(num_inputs, num_outputs, 2,
                                          TasGrid::type_level, TasGrid::rule_rleja);

    cout << setw(10) << "points" << setw(20) << "error\n";

    int budget = 25; // try budget 50, 100, 200, 400, etc.

    // running examples with ever increasing budget
    // the budget is defined as the sum of loaded + new points
    for(int itr=0; itr<6; itr++){ // run 6 example iterations

        budget *= 2; // try budget 100, 200, 400, ...

        TasGrid::constructSurrogate<mode>(model, budget, 12, 1, grid,
                                          TasGrid::type_iptotal, 0);

        cout << setw(10) << grid.getNumPoints() << setw(20) << test_grid(grid) << "\n";
    }

    // Note: when using a sequence rule, the actual number of points always matches the budget
    // but this cannot be guaranteed for rules that grow by more than one node per level
    // because a tensor can be added only when all points are present and each tensor requires
    // more than one point
    // for fast growing rules, sometime it is simply not possible to build a grid
    // with exactly the given number of points

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_06 example]
#endif
}
