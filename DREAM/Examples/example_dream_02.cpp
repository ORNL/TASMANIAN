#include "Tasmanian.hpp"

using namespace std;

/*!
 * \internal
 * \file example_dream_02.cpp
 * \brief Examples for the Tasmanian DREAM module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAMExamples
 *
 * Tasmanian DREAM Example 2
 * \endinternal
 */

/*!
 * \ingroup TasmanianDREAMExamples
 * \addtogroup TasmanianDREAMExamples2 Tasmanian DREAM module, example 2
 *
 * Example 2:
 * Demonstrates how to solve a Bayesian inference problem for parameter calibration
 * using a sparse grid approximation to the Bayesian likelihood function.
 */

//! \brief DREAM Example 2: perform parameter calibration using likelihood that is approximated with a sparse grid.
//! \ingroup TasmanianDREAMExamples2

//! \snippet DREAM/Examples/example_dream_02.cpp DREAM_Example_02 example
void dream_example_02(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_02 example]
#endif
    // using the default random engine, but must reset the random number generator
    srand((int) time(nullptr));

    // Example 2
    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 2: set the inference problem: identify model parameters x_0 and x_1\n"
         << "           from data (noise free example)\n"
         << "           model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(M_PI*t + 0.3*M_PI)\n"
         << "           t in [0,1], t is discretized with 32 equidistant nodes\n"
         << "           the likelihood is exp(- 16 * (f(x) - d)^2)\n"
         << "           using a sparse grid to interpolate the likelihood\n"
         << "     NOTE: 16 = 32/2 corresponds to the discretization error in t\n\n";

    constexpr double pi = 3.14159265358979323846;

    int num_dimensions = 2;
    int num_chains = 100;
    int num_burnup_iterations = 3000;
    int num_sample_iterations = 100;
    // the total number of samples is num_chains * num_iterations = 100 * 100 = 10,000
    int num_discrete_nodes = 32;

    // create a lambda function that represents the model,
    // normally this would be a call to an external code
    auto model = [&](double x0, double x1, std::vector<double> &data)->
        void{
            double dt = 1.0 / ((double) data.size());
            double t = 0.5 * dt;
            for(auto &d : data){
                d = std::sin(x0 * pi * t + x1);
                t += dt;
            }
        };

    // create a lambda function that represents the likelihood
    // the Tasmanian LikelihoodGaussIsotropic class can be used here,
    // but the objective of the example is to use sparse grid interpolant for the likelihood
    // thus we assume that the likelihood formula is potentially much more complex
    auto likelihood = [&](double scale,
                          std::vector<double> &model_values,
                          std::vector<double> &data)->
        double{
            double likelihood_value = 0.0;
            for(size_t j=0; j<model_values.size(); j++){
                likelihood_value += (model_values[j] - data[j]) * (model_values[j] - data[j]);
            }
            return - 0.5 * scale * likelihood_value;
        };

    // generate the data, normally this will come from an external source
    std::vector<double> data(num_discrete_nodes);
    model(1.0, 0.3 * pi, data); // "true" values: x0 = 1.0, x1 = 0.3 * pi

    // construct a sparse grid approximation for the log of the likelihood
    // the sampling will be done in logform
    // for details consult the documentation of the sparse grids module

    // define the domain over which the sparse grid will be constructed
    std::vector<double> domain_a = {0.5, -0.1}, domain_b = {8.0,  1.7};

    auto grid = TasGrid::makeGlobalGrid(num_dimensions, 1, 10, // using 10-th order polynomial
                                        TasGrid::type_iptotal, TasGrid::rule_clenshawcurtis);
    grid.setDomainTransform(domain_a, domain_b);

    TasGrid::loadNeededPoints<TasGrid::mode_sequential>(
        [&](std::vector<double> const &x, std::vector<double> &y, size_t)->void{
            std::vector<double> model_at_point(num_discrete_nodes);
            model(x[0], x[1], model_at_point);
            y[0] = likelihood((double) num_discrete_nodes, model_at_point, data);
        },
        grid, 1);

    TasDREAM::TasmanianDREAM state(num_chains, grid); // get the dimensions from the sparse grid
    std::vector<double> init_chains = TasDREAM::genUniformSamples(domain_a, domain_b, num_chains);
    state.setState(init_chains); // use chains distributed uniformly over the domain

    // Call to Tasmanian DREAM Sampling algorithm
    TasDREAM::SampleDREAM<TasDREAM::logform>
                          (num_burnup_iterations, num_sample_iterations,
                           TasDREAM::posterior(grid, TasDREAM::uniform_prior),
                           grid.getDomainInside(),
                           state,
                           TasDREAM::dist_uniform, 0.5, // independent update of magnitude 0.5
                           TasDREAM::const_percent<90> // use 90% differential update
                           );

    std::vector<double> expectation, variance;
    state.getHistoryMeanVariance(expectation, variance);

    cout << "Inferred values (using 10th order polynomial sparse grid):\n";
    cout << "   frequency:" << setw(12) << std::fixed << expectation[0]
         << "   error:" << setw(12) << std::scientific << std::abs(expectation[0] - 1.0) << "\n";
    cout << "  correction:" << setw(12) << std::fixed << expectation[1]
         << "   error:" << setw(12) << std::scientific << std::abs(expectation[1] - 0.3 * pi)
         << "\n\n";

    // Solve the same example, but switch to 30-th order polynomial with a Sequence grid

    // rebuild the sparse grid interpolant
    // use 30-th order polynomial and the faster Sequence grid
    grid.makeSequenceGrid(2, 1, 30, TasGrid::type_iptotal, TasGrid::rule_leja);
    grid.setDomainTransform(domain_a, domain_b);

    TasGrid::loadNeededPoints<TasGrid::mode_sequential>(
        [&](std::vector<double> const &x, std::vector<double> &y, size_t)->
        void{
            std::vector<double> model_at_point(num_discrete_nodes);
            model(x[0], x[1], model_at_point);
            y[0] = likelihood((double) num_discrete_nodes, model_at_point, data);
        },
        grid, 1);

    state.clearPDFvalues(); // forget the previous history
    state.setState(init_chains); // reset the chains

    // Call to Tasmanian DREAM Sampling algorithm
    TasDREAM::SampleDREAM<TasDREAM::logform>
                         (num_burnup_iterations, num_sample_iterations,
                          TasDREAM::posterior<TasDREAM::logform>(grid, TasDREAM::uniform_prior),
                          grid.getDomainInside(),
                          state,
                          TasDREAM::dist_uniform, 0.5, // independent update of magnitude 0.5
                          TasDREAM::const_percent<90> // use 90% differential update
                          );

    state.getHistoryMeanVariance(expectation, variance); // recompute the statistics

    cout << "Inferred values (using 30th order polynomial sparse grid):\n";
    cout << "   frequency:" << setw(12) << std::fixed << expectation[0]
         << "   error:" << setw(12) << std::scientific << std::abs(expectation[0] - 1.0) << "\n";
    cout << "  correction:" << setw(12) << std::fixed << expectation[1]
         << "   error:" << setw(12) << std::scientific << std::abs(expectation[1] - 0.3 * pi)
         << "\n\n";

    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_02 example]
#endif
}
