#include "Tasmanian.hpp"

using namespace std;

/*!
 * \internal
 * \file example_optimization_01.cpp
 * \brief Examples for the Tasmanian Optimization module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianOPTExamples
 *
 * Tasmanian Optimization Example 1
 * \endinternal
 */

/*!
 * \ingroup TasmanianOPTExamples
 * \addtogroup TasmanianOPTExamples1 Tasmanian Optimization module, example 1
 *
 * Example 1: Particle Swarm method
 */

/*!
 * \ingroup TasmanianOPTExamples1
 * \brief Optimization Example 1: Demonstrates the use of the Particle Swarm method.
 *
 * Find the minimum of the six-hump camel function
 * \f$  f(x,y) = ( 4 - 2.1 x^2 + x^4 / 3) x^2 + x y + ( - 4 + 4 y^2) y^2 \f$
 * the problem is challenging due to the multiple relative and global extrema.
 * Classic gradient based methods often stagnate and fail when applied to this benchmark problem.
 * In contrast, the Particle Swarm method uses multiple "particles" that move around the domain
 * in search for an optimal position.
 * The "swarm" shares global information and is fairly insensitive to local extrema.
 * The algorithm shares many similarities with the DREAM sampling procedure
 * and is a good fit for the Tasmanian framework.
 * While the method is probabilistic, identifying a correct global minimum
 * comes with a high probability of success.
 */

//! \snippet DREAM/Examples/example_optimization_01.cpp OPT_Example_01 example
void optimizaiton_example_01(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [OPT_Example_01 example]
#endif
    // using the default random engine, but must reset the random number generator
    std::srand(std::time(nullptr));

    // Example 1:
    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 1: use the Particle Swarm algorithm to minimize the six-hump camel function\n"
         << "           f(x,y) = (4-2.1 *x^2 + x^4 / 3) * x^2 + x * y + (-4 + 4 * y^2) * y^2\n"
         << "           domain is x in (-3, +3), y in (-2, +2)\n"
         << "           using 50 particles and 200 iterations\n"
         << "    See the comments in example_optimization_01.cpp\n\n";

     cout << "the problem has two solutions at (-8.98420e-02, 7.12656e-01) and (8.98420e-02, -7.12656e-01)\n"
          << "the objective at the solutions is -1.03163.\n\n";

     int num_dimensions = 2;
     int num_particles = 50;

     TasOptimization::ParticleSwarmState state(num_dimensions, num_particles);
     state.initializeParticlesInsideBox({-3.0, -2.0}, {3.0, 2.0});

     auto objective = TasOptimization::makeObjectiveFunction(
                    num_dimensions,
                    [](const std::vector<double> &x)->double {
                         return (4.0 - 2.1 * x[0]*x[0] + x[0]*x[0]*x[0]*x[0] / 3.0) * x[0]*x[0] +
                              x[0] * x[1] +
                              (-4.0 + 4.0 * x[1]*x[1]) * x[1]*x[1];
                    });

     int num_iterations = 200;
     double inertia_weight = 0.5, cognitive_coeff = 2.0, social_coeff = 2.0;
     TasOptimization::ParticleSwarm(objective,
                                    TasDREAM::hypercube({-3.0, -2.0}, {3.0, 2.0}),
                                    inertia_weight, cognitive_coeff, social_coeff,
                                    num_iterations, state);

     std::vector<double> best_position = state.getBestPosition();
     std::vector<double> best_objective(1);
     objective(best_position, best_objective);

     cout << "Using the objective function\n"
          << "found best position after 200 iteration:\n"
          << "best value for x = " << best_position[0] << "\n"
          << "best value for y = " << best_position[1] << "\n"
          << "value of the objective = " << best_objective[0] << "\n\n";


     // solve the same optimization problem, but use a sparse grid surrogate
     // first we create the surrogate, then we optimize

     auto grid = TasGrid::makeLocalPolynomialGrid(2, 1, 10);
     grid.setDomainTransform({-3.0, -2.0}, {3.0, 2.0});

     std::vector<double> points = grid.getNeededPoints();
     std::vector<double> values(grid.getNumNeeded());
     objective(points, values);
     grid.loadNeededValues(values); // at this point, we have the surrogate

     // reset the state
     state = TasOptimization::ParticleSwarmState(num_dimensions, num_particles);
     state.initializeParticlesInsideBox({-3.0, -2.0}, {3.0, 2.0});

     // the values of the objective function for all particle positions are computed in batch
     // thus, the problem is amenable to GPU acceleration
     // Note: if GPU acceleration is not available, Tasmanian will automatically fallback to CPU
     grid.enableAcceleration(TasGrid::accel_gpu_cuda);
     TasOptimization::ParticleSwarm(
                              [&](std::vector<double> const &x, std::vector<double> &y)->void{
                                   grid.evaluateBatch(x, y);
                              },
                              TasDREAM::hypercube({-3.0, -2.0}, {3.0, 2.0}),
                              inertia_weight, cognitive_coeff, social_coeff,
                              num_iterations, state);

     best_position = state.getBestPosition();
     objective(best_position, best_objective);

     cout << "Using the surrogate to the objective function\n"
          << "found best position after 200 iteration:\n"
          << "best value for x = " << best_position[0] << "\n"
          << "best value for y = " << best_position[1] << "\n"
          << "value of the objective = " << best_objective[0] << "\n\n";

     cout << "Note: the method will find only one of the minimums\n"
          << "      the random seed will determine which one.\n";


    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [OPT_Example_01 example]
#endif
}
