#include "Tasmanian.hpp"

using namespace std;

/*!
 * \internal
 * \file example_optimization_02.cpp
 * \brief Examples for the Tasmanian Optimization module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianOPTExamples
 *
 * Tasmanian Optimization Example 2
 * \endinternal
 */

/*!
 * \ingroup TasmanianOPTExamples
 * \addtogroup TasmanianOPTExamples2 Tasmanian Optimization module, example 2
 *
 * Example 2: Gradient Descent method
 */

/*!
 * \ingroup TasmanianOPTExamples2
 * \brief Optimization Example 2: Demonstrates the use of the Gradient Descent method
 *
 * The Gradient Descent method is both fast converging and very stable,
 * provided that the objective function has a single relative minimum (e.g., convex),
 * or the initial starting point is close to the global minimum.
 * The method is a good complement to the Particle Swarm, where the more slowly
 * converging particle method can identify the global extrema for a non-convex problem,
 * and few gradient iterations can quickly zoom into the solution.
 *
 * This example uses a very simple quadratic function, the primary goal here is to
 * demonstrate the syntax rather than compare, contrast or combine different optimization methods.
 *
 * \snippet DREAM/Examples/example_optimization_02.cpp OPT_Example_02 example
 */

void optimizaiton_example_02(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [OPT_Example_02 example]
#endif
    // Example 2:
    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 2: use Gradient Descent algorithms on a simple quadratic\n"
         << "           f(x,y) =  2.0 * (x - 1.0) * (x - 1.0) + (y - 2.0) * (y - 2.0) / 2.0\n"
         << "    See the comments in example_optimization_02.cpp\n\n";

     cout << "Exact solution is at (1.0, 2.0)\n\n";

     auto objective = [](std::vector<double> const &x)
          ->double{
               return 2.0 * (x[0] - 1.0) * (x[0] - 1.0) + (x[1] - 2.0) * (x[1] - 2.0) / 2.0;
          };
     // the gradient descent methods require the gradient of the objective
     auto gradient = [](std::vector<double> const &x, std::vector<double> &grad)
          ->void{
               grad[0] = 4.0 * (x[0] - 1.0);
               grad[1] = x[1] - 2.0;
          };

     int num_dimensions = 2;
     double initial_stepsize = 0.0;

     // the {0.0, 0.0} is a vector with initial state, e.g., initial guess
     // coupling particle swarm and gradient descent can happen with
     // auto gradient_state =
     //              TasOptimization::GradientDescentState(particle_state.getBestPosition(), 0.0);
     TasOptimization::GradientDescentState state({0.0, 0.0}, initial_stepsize);

     int max_iterations = 200;
     double tolerance = 1.E-3;
     TasOptimization::OptimizationStatus status =
          TasOptimization::GradientDescent(gradient, 1.0/8.0, max_iterations, tolerance, state);

     std::vector<double> opt_x = state.getX();

     cout << "Using the objective function\n"
          << "aiming at tolerance 1.e-3\n"
          << "performed iterations: " << status.performed_iterations << "\n"
          << "best value for x = " << opt_x[0] << "\n"
          << "best value for y = " << opt_x[1] << "\n\n";

     // same problem, but using a surrogate model
     // construct the surrogate
     auto grid = TasGrid::makeSequenceGrid(num_dimensions, 1, 2,
                                           TasGrid::type_iptotal, TasGrid::rule_leja);
     std::vector<double> points = grid.getNeededPoints();
     std::vector<double> values(grid.getNumNeeded());
     for(int i=0; i<grid.getNumNeeded(); i++){
          values[i] = objective(std::vector<double>(points.begin() + 2 * i,
                                                    points.begin() + 2 * i + 2));
     }
     grid.loadNeededValues(values);

     // reset the state to zero
     state = TasOptimization::GradientDescentState({0.0, 0.0}, initial_stepsize);
     status = TasOptimization::GradientDescent(
          [&](std::vector<double> const &x, std::vector<double> &grad)
          ->void{
               grid.differentiate(x, grad);
          },
          1.0/8.0, max_iterations, tolerance, state);

     opt_x = state.getX();

     // actually, since the objective is quadratic, the surrogate is an exact match
     // the example demonstrates syntax more than anything else
     cout << "Using the surrogate function\n"
          << "aiming at tolerance 1.e-3\n"
          << "performed iterations: " << status.performed_iterations << "\n"
          << "best value for x = " << opt_x[0] << "\n"
          << "best value for y = " << opt_x[1] << "\n\n";

     // finding optimum inside a box
     // define the projection operator of a vector into a box
     // restriction function is probably better name here, since the projection
     // is not orthogonal
     // Note: that such restriction is necessary when working with a sparse grid
     //       since the surrogates are usually restricted to a hypercube box.
     auto projection = [](const std::vector<double> &x, std::vector<double> &p) {
          p[0] = std::min(std::max(x[0], 0.0), 0.5);
          p[1] = std::min(std::max(x[1], 0.0), 1.5);
     };

     // reset the state and take the initial_stepsize to 1.0
     state = TasOptimization::GradientDescentState({0.0, 0.0}, 1.0);

     double increase_coeff = 1.25;
     double decrease_coeff = 1.25;
     status = TasOptimization::GradientDescent(objective, gradient, projection,
                                               increase_coeff, decrease_coeff,
                                               max_iterations, tolerance, state);

     opt_x = state.getX();

     // actually, since the objective is quadratic, the surrogate is an exact match
     // the example demonstrates syntax more than anything else
     cout << "Using the projection inside of a box\n"
          << "the solution is at the edge of the box at (0.5, 1.5)\n"
          << "performed iterations: " << status.performed_iterations << "\n"
          << "best value for x = " << opt_x[0] << "\n"
          << "best value for y = " << opt_x[1] << "\n";

    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [OPT_Example_02 example]
#endif
}
