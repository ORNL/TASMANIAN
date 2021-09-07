
#include "Tasmanian.hpp"
#include <random>

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_08.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 8.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples8 Tasmanian Sparse Grids module, example 8
 *
 * \par Example 8
 * Different local polynomial rules.
 */

/*!
 * \ingroup TasmanianSGExamples8
 * \brief Sparse Grids Example 8: local polynomial rules
 *
 * Local polynomial grids use a hierarchy of basis functions with decreasing support,
 * Tasmanian offers several different types of basis each tuned to different types
 * of models. This example demonstrates the efficiency (i.e., error per number of points)
 * of different local polynomial grids when applied to different models.
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_08.cpp SG_Example_08 example
 */
void sparse_grids_example_08(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_08 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Example 8: interpolate different functions demonstrating the different\n"
         << "           local polynomial rules\n\n";

    int const num_inputs = 2; // using two inputs for all tests

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

    // test 1: using smooth function
    auto smooth_model = [](double const x[], double y[], size_t)->
        void{ y[0] = std::exp(-x[0] * x[0]) * std::cos(x[1]); };

    int const order = 2; // at order 1 localp and semilocalp are identical
    auto grid_localp = TasGrid::makeLocalPolynomialGrid(num_inputs, 1, 7, order,
                                                        TasGrid::rule_localp);
    auto grid_semilocalp = TasGrid::makeLocalPolynomialGrid(num_inputs, 1, 7, order,
                                                            TasGrid::rule_semilocalp);

    TasGrid::loadNeededValues(smooth_model, grid_localp, 4);
    TasGrid::loadNeededValues(smooth_model, grid_semilocalp, 4);

    cout << "Using smooth model: f(x, y) = exp(-x*x) * cos(y)\n"
         << " rule_localp,     points = " << grid_localp.getNumPoints()
         << "   error = " << get_error(grid_localp, smooth_model) << "\n"
         << " rule_semilocalp, points = " << grid_semilocalp.getNumPoints()
         << "   error = " << get_error(grid_semilocalp, smooth_model) << "\n"
         << " If the model is smooth, rule_semilocalp has an advantage.\n\n";

    // test 1: using model with zero-boundary conditions
    constexpr double pi = 3.14159265358979323846;
    auto zero_model = [=](double const x[], double y[], size_t)->
        void{ y[0] = std::cos(0.5 * pi * x[0]) * std::cos(0.5 * pi * x[1]); };

    auto grid_localp0 = TasGrid::makeLocalPolynomialGrid(num_inputs, 1, 6, order,
                                                        TasGrid::rule_localp0);

    // the true value indicates overwrite of the currently loaded model
    TasGrid::loadNeededValues<TasGrid::mode_parallel, true>(zero_model, grid_localp, 4);
    TasGrid::loadNeededValues(zero_model, grid_localp0, 4);

    cout << "Using homogeneous model: f(x, y) = cos(pi * x / 2) * cos(pi * y / 2)\n"
         << " rule_localp,   points = " << grid_localp.getNumPoints()
         << "   error = " << get_error(grid_localp, zero_model) << "\n"
         << " rule_localp0,  points = " << grid_localp0.getNumPoints()
         << "   error = " << get_error(grid_localp0, zero_model) << "\n"
         << " The rule_localp0 uses basis tuned for models with zero boundary.\n";

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_08 example]
#endif
}
