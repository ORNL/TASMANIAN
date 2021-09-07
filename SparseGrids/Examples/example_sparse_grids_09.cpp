
#include "Tasmanian.hpp"
#include <random>

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_09.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 9.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples9 Tasmanian Sparse Grids module, example 9
 *
 * \par Example 9
 * Different local polynomial refinement.
 */

/*!
 * \ingroup TasmanianSGExamples9
 * \brief Sparse Grids Example 9: local polynomial refinement
 *
 * Example 8 showed the difference between the local polynomial rules and how
 * the basis can be tuned to a specific model. This example shows refinement
 * in the local polynomial context using hierarchical coefficients (surpluses).
 * A sharp model is selected so that TasGrid::rule_localp is the natural choice,
 * and two different refinement strategies are tested.
 *
 * \b Note: the refinement process can be used in the context of construction,
 * e.g., similar to Example 6 and TasGrid::constructSurrogate().
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_09.cpp SG_Example_09 example
 */
void sparse_grids_example_09(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_09 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Example 9: comparison between local polynomial refinement strategies\n\n";

    int const num_inputs = 2;

    // using random points to test the error
    int const num_test_points = 1000;
    std::vector<double> test_points(num_test_points * num_inputs);
    std::minstd_rand park_miller(42);
    std::uniform_real_distribution<double> domain(-1.0, 1.0);
    for(auto &t : test_points) t = domain(park_miller);

    // computes the error between the gird surrogate model and the actual model
    // using the test points, finds the largest absolute error
    auto get_error = [&](TasGrid::TasmanianSparseGrid const &grid,
                         std::function<void(double const x[], double y[], size_t)> model)->
        double{
            std::vector<double> grid_result;
            grid.evaluateBatch(test_points, grid_result);
            double err = 0.0;
            for(int i=0; i<num_test_points; i++){
                double model_result; // using only one output
                model(&test_points[i*num_inputs], &model_result, 0);
                err = std::max(err, std::abs(grid_result[i] - model_result));
            }
            return err;
        };

    // using a model with sharp transition, the modified rules are not
    // suitable for this model, using rule_localp
    auto sharp_model = [](double const x[], double y[], size_t)->
        void{ y[0] = exp(-x[0]) / (1.0 + 100.0 * exp(-10.0 * x[1])); };

    // using maximal order
    int const order = -1;
    auto grid_classic = TasGrid::makeLocalPolynomialGrid(num_inputs, 1, 2, order,
                                                         TasGrid::rule_localp);
    auto grid_fds = TasGrid::copyGrid(grid_classic);

    // setup the top rows of the output table
    cout << "Using batch refinement:\n"
         << setw(22) << "classic" << setw(22) << "fds\n"
         << setw(8) << "points" << setw(14) << "error"
         << setw(8) << "points" << setw(14) << "error\n";

    double const tolerance = 1.E-5;
    while((grid_classic.getNumNeeded() > 0) || (grid_fds.getNumNeeded() > 0)){
        TasGrid::loadNeededValues(sharp_model, grid_classic, 4);
        TasGrid::loadNeededValues(sharp_model, grid_fds, 4);

        // output the grids and results at the current stage
        cout << setw(8) << grid_classic.getNumLoaded()
             << setw(14) << get_error(grid_classic, sharp_model)
             << setw(8) << grid_fds.getNumLoaded()
             << setw(14) << get_error(grid_fds, sharp_model) << "\n";

        // setting refinement for each grid
        grid_classic.setSurplusRefinement(tolerance, TasGrid::refine_classic);
        grid_fds.setSurplusRefinement(tolerance, TasGrid::refine_fds);
    }

    // Repeat the example using construction
    auto vector_model = [=](std::vector<double> const &x,
                           std::vector<double> &y,
                           size_t tid)->
        void{
            y.resize(1);
            sharp_model(x.data(), y.data(), tid);
        };

    // reset the grids
    grid_classic = TasGrid::makeLocalPolynomialGrid(num_inputs, 1, 2, order,
                                                         TasGrid::rule_localp);
    grid_fds = TasGrid::copyGrid(grid_classic);

    // setup the top rows of the output table
    cout << "\nUsing construction:\n"
         << setw(22) << "classic" << setw(22) << "fds\n"
         << setw(8) << "points" << setw(14) << "error"
         << setw(8) << "points" << setw(14) << "error\n";

    for(size_t budget = 50; budget < 800; budget += 100){
        TasGrid::constructSurrogate(vector_model, budget, 1, 1, grid_classic,
                                    tolerance, TasGrid::refine_classic);
        TasGrid::constructSurrogate(vector_model, budget, 1, 1, grid_fds,
                                    tolerance, TasGrid::refine_fds);

        // output the grids and results at the current stage
        cout << setw(8) << grid_classic.getNumLoaded()
             << setw(14) << get_error(grid_classic, sharp_model)
             << setw(8) << grid_fds.getNumLoaded()
             << setw(14) << get_error(grid_fds, sharp_model) << "\n";
    }

    // The FDS strategy does not win at every step of the way, but at the final
    // stage the FDS approach results in fewer nodes while keeping the error
    // at near the same magnitude

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_09 example]
#endif
}
