#include <iostream>
#include <iomanip>
#include <ctime>
#include "math.h"

#include "TasmanianDREAM.hpp"

using namespace std;

//! \file example_dream_02.cpp
//! \brief Examples for the Tasmanian DREAM module.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAMExamples
//!
//! Tasmanian DREAM Example 2

//! \defgroup TasmanianDREAMExamples2 Tasmanian DREAM module, example 2.
//!
//! Example 2:
//! Demonstrates how to solve a Bayesian inference problem for parameter calibration
//! using a sparse grid approximation to the Bayesian likelihood function.

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
    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 2: set inference problem, identify x_0 and x_1 model parameters from data (noise free example)" << endl;
    cout << "           model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(M_PI*t + 0.3*M_PI)" << endl;
    cout << "           t in [0,1], discretized with 32 equidistant nodes" << endl;
    cout << "           likelihood is exp(- 16 * (f(x) - d)^2)" << endl;
    cout << "           use sparse grid to interpolate the likelihood" << endl;
    cout << "     NOTE: 32 corresponds to discretization error in t" << endl << endl;

    int num_dimensions = 2;
    int num_chains = 100;
    int num_burnup_iterations = 3000;
    int num_sample_iterations = 100;
    // the total number of samples is num_chains * num_iterations = 100 * 100 = 10,000
    int num_discrete_nodes = 32;

    // create a lambda function that represents the model, normally this would be a call to an external code
    auto model = [&](double x0, double x1, std::vector<double> &data)->
        void{
            double dt = 1.0 / ((double) data.size());
            double t = 0.5 * dt;
            for(auto &d : data){
                d = sin(x0 * M_PI * t + x1);
                t += dt;
            }
        };

    // create a lambda function that represents the likelihood
    // the Tasmanian LikelihoodGaussIsotropic class can be used here,
    // but the objective of the example is to use sparse grid interpolant for the likelihood
    // thus we assume that the likelihood formula is potentially much more complex
    auto likelihood = [&](double scale, std::vector<double> &model_values, std::vector<double> &data)->
        double{
            double likelihood_value = 0.0;
            for(size_t j=0; j<model_values.size(); j++){
                likelihood_value += (model_values[j] - data[j]) * (model_values[j] - data[j]);
            }
            return - 0.5 * scale * likelihood_value;
        };

    // generate the data, normally this will come from an external source
    std::vector<double> data(num_discrete_nodes);
    model(1.0, 0.3 * M_PI, data); // "true" values: x0 = 1.0, x1 = 0.3 * M_PI

    // construct a sparse grid approximation for the log of the likelihood (i.e., will use logform later)
    // for details consult the documentation of the sparse grids module

    // define the domain over which the sparse grid will be constructed
    std::vector<double> domain_a = {0.5, -0.1}, domain_b = {8.0,  1.7};

    TasGrid::TasmanianSparseGrid grid;
    grid.makeGlobalGrid(num_dimensions, 1, 10, TasGrid::type_iptotal, TasGrid::rule_clenshawcurtis); // using 10-th order polynomial
    grid.setDomainTransform(domain_a, domain_b);

    std::vector<double> points; // get the sparse grid points
    grid.getNeededPoints(points);
    std::vector<double> values(grid.getNumNeeded()); // allocate memory for the values

    for(int i=0; i<grid.getNumNeeded(); i++){ // compte the likelihood for each sparse grid point (log-form)
        std::vector<double> model_at_point(num_discrete_nodes);
        model(points[2*i], points[2*i+1], model_at_point);
        values[i] = likelihood((double) num_discrete_nodes, model_at_point, data);
    }
    grid.loadNeededPoints(values); // load the values into the grid

    TasDREAM::TasmanianDREAM state(num_chains, grid); // assume the dimensions from the sparse grid
    std::vector<double> initial_chains;
    TasDREAM::genUniformSamples(domain_a, domain_b, num_chains, initial_chains);
    state.setState(initial_chains); // use chains distributed uniformly over the domain

    // Ccall to Tasmanian DREAM Sampling algorithm
    TasDREAM::SampleDREAMGrid<TasDREAM::logform>
                              (num_burnup_iterations, num_sample_iterations,
                               grid, // provide the likelihood
                               TasDREAM::uniform_prior, // assume non-informative prior
                               TasDREAM::dist_uniform, 0.5, // uniform independent update of magnitude 0.5
                               state,
                               TasDREAM::const_percent<90> // use 90% of differential update
                              );

    std::vector<double> expectation, variance;
    state.getHistoryMeanVariance(expectation, variance);

    cout << "Inferred values (using 10th order polynomial sparse grid):" << endl;
    cout << "   frequency:" << setw(12) << std::fixed << expectation[0]
         << "   error:" << setw(12) << std::scientific << fabs(expectation[0] - 1.0) << endl;
    cout << "  correction:" << setw(12) << std::fixed << expectation[1]
         << "   error:" << setw(12) << std::scientific << fabs(expectation[1] - 0.3 * M_PI) << endl << endl;

    // Solve the same example, but switch to 30-th order polynomial with Sequence grid

    // rebuild the sparse grid interpolant
    grid.makeSequenceGrid(2, 1, 30, TasGrid::type_iptotal, TasGrid::rule_leja); // 30-th order polynomial (using faster Sequence grid)
    grid.setDomainTransform(domain_a, domain_b);

    grid.getNeededPoints(points); // get the sparse grid points
    values.resize(grid.getNumNeeded()); // allocate memory for the values

    for(int i=0; i<grid.getNumNeeded(); i++){ // compte the likelihood for each sparse grid point (log-form)
        std::vector<double> model_at_point(num_discrete_nodes);
        model(points[2*i], points[2*i+1], model_at_point);
        values[i] = likelihood((double) num_discrete_nodes, model_at_point, data);
    }
    grid.loadNeededPoints(values); // load the values into the grid

    state.clearPDFvalues(); // forget the previous history
    state.setState(initial_chains); // reset the chains

    // Ccall to Tasmanian DREAM Sampling algorithm
    TasDREAM::SampleDREAMGrid<TasDREAM::logform>
                              (num_burnup_iterations, num_sample_iterations,
                               grid, // provide the likelihood
                               TasDREAM::uniform_prior, // assume non-informative prior
                               TasDREAM::dist_uniform, 0.5, // uniform independent update of magnitude 0.5
                               state,
                               TasDREAM::const_percent<90> // use 90% of differential update
                              );

    state.getHistoryMeanVariance(expectation, variance); // recompute the statistics

    cout << "Inferred values (using 30th order polynomial sparse grid):" << endl;
    cout << "   frequency:" << setw(12) << std::fixed << expectation[0]
         << "   error:" << setw(12) << std::scientific << fabs(expectation[0] - 1.0) << endl;
    cout << "  correction:" << setw(12) << std::fixed << expectation[1]
         << "   error:" << setw(12) << std::scientific << fabs(expectation[1] - 0.3 * M_PI) << endl << endl;

    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_02 example]
#endif
}
