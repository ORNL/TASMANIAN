#include <iostream>
#include <iomanip>
#include <ctime>

#include "TasmanianDREAM.hpp"

using namespace std;

//! \file example_dream_04.cpp
//! \brief Examples for the Tasmanian DREAM module.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAMExamples
//!
//! Tasmanian DREAM Example 4

/*!
 * \ingroup TasmanianDREAMExamples
 * \addtogroup TasmanianDREAMExamples4 Tasmanian DREAM module, example 4
 *
 * Example 4:
 * Given data that is the superposition of two exponential curves,
 * use Sparse Grid and Bayesian inference to identify the scale and growth of each curve.
 * Higher dimensions and bi-modal posterior will decrease the acceptance rate,
 * thus we need more samples than the previous example.
 */

//! \brief DREAM Example 4: signal decomposition, using bi-modal posterior distribution.
//! \ingroup TasmanianDREAMExamples4

//! \snippet DREAM/Examples/example_dream_04.cpp DREAM_Example_04 example
void dream_example_04(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_04 example]
#endif
    // using the default random engine, but must reset the random number generator
    srand((int) time(nullptr));

    // EXAMPLE 4:
    cout << "\n" << "-------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 4: similar to Example 3, but data is noisy and the multi-modal posterior comes from symmetry in the model\n";
    cout << "           set inference problem, identify x_0, x_1, x_2 and x_3 (four) model parameters from noisy data\n";
    cout << "           model: f(x) = x_0*exp(x_1*t) + x_2*exp(x_3*t),  data: d = exp(t) + 0.4*exp(3*t)\n";
    cout << "           note the symmetry between x_1 and x_3 which causes bi-modal posterior\n";
    cout << "           t in [0,1], discretized with 64 equidistant nodes\n";
    cout << "           likelihood is exp(-64 * (f(x) - d)^2)\n";
    cout << "           using sparse grid to interpolate the model\n";
    cout << "     NOTE: 64 corresponds to discretization error in t\n" << endl;

    // multi-modal distributions require a lot more chains and samples
    int num_chains = 50;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 1000;
    // the total number of samples is num_chains * num_iterations
    int num_discrete_nodes = 64;

    // create a lambda function that represents the model, normally this would be a call to an external code
    auto model = [&](double x0, double x1, double x2, double x3, std::vector<double> &data)->
        void{
            double dt = 1.0 / ((double) data.size());
            double t = 0.5 * dt;
            for(auto &d : data){
                d = x0 * std::exp(x1 * t) + x2 * std::exp(x3 * t);
                t += dt;
            }
        };

    // create data corresponding to 1.0, 1.0, 0.4, 3.0
    std::vector<double> data(num_discrete_nodes);
    model(1.0, 1.0, 0.4, 3.0, data);

    // add noise to the data, use magnitude (standard deviation) of 1 / num_discrete_nodes
    // you can adjust the example to consider more/less noise
    TasDREAM::applyGaussianUpdate(data, 1.0 / ((double) num_discrete_nodes));

    TasGrid::TasmanianSparseGrid grid;
    grid.makeSequenceGrid(4, num_discrete_nodes, 15, TasGrid::type_iptotal, TasGrid::rule_leja); // 15-th order polynomial
    std::vector<double> domain_a = {0.2, 0.5, 0.2, 0.5}; // set search interval x_0, x_2 in [0.2, 1.2], x_1, x_3 in [0.5, 4.0]
    std::vector<double> domain_b = {1.2, 4.0, 1.2, 4.0};
    grid.setDomainTransform(domain_a, domain_b);

    std::vector<double> points; // get the sparse grid points
    grid.getNeededPoints(points);
    std::vector<double> values; // stores the model values

    for(int i=0; i<grid.getNumNeeded(); i++){ // compute the model for each sparse grid point
        std::vector<double> model_at_point(num_discrete_nodes);
        model(points[4*i], points[4*i+1], points[4*i+2], points[4*i+3], model_at_point);
        values.insert(values.end(), model_at_point.begin(), model_at_point.end()); // append to the total vector
    }
    grid.loadNeededPoints(values); // load the values into the grid

    // when working with grids with large N, acceleration becomes very useful
    // if BLAS is enabled on compile time, the grid will use BLAS by default
    // GPU acceleration can be enabled here using:
    //   grid.enableAcceleration(TasGrid::accel_gpu_cuda);
    //   grid.setGPUID(0);

    // define the likelihood function and load the data
    // even though the example is noise free, we assume that the "noise" is due to discretization error
    // using piece-wise constant approximation the discretization error is 1.0 / num_discrete_nodes
    TasDREAM::LikelihoodGaussIsotropic likely(1.0 / ((double) num_discrete_nodes), data);

    TasDREAM::TasmanianDREAM state(num_chains, grid); // assume the dimensions from the sparse grid
    std::vector<double> initial_chains;
    TasDREAM::genUniformSamples(domain_a, domain_b, num_chains, initial_chains);
    state.setState(initial_chains); // use chains distributed uniformly over the domain

    // Call to Tasmanian DREAM Sampling algorithm
    TasDREAM::SampleDREAM<TasDREAM::logform>
                         (num_burnup_iterations, num_sample_iterations,
                          TasDREAM::posterior<TasDREAM::logform>(likely, // provide the likelihood
                                                                 grid, // provide the model
                                                                 TasDREAM::uniform_prior), // assume non-informative prior
                          grid.getDomainInside(),
                          TasDREAM::dist_gaussian, 0.01, // Gaussian independent update of magnitude 0.01
                          state,
                          TasDREAM::const_percent<100> // use 100% of differential update
                          );

    // get the vector containing the sampling history, could also use "auto history = state.getHistory();"
    const std::vector<double> &history = state.getHistory();

    // splitting the bimodal data is tricky,
    // add constraint that orders the rates
    double rate_low = 0.0,  scale_low = 0.0;
    double rate_high = 0.0, scale_high = 0.0;
    for(size_t i=0; i<history.size(); i+=4){
        // group the samples by looking at the second and fourth input
        if (history[i+1] < history[i+3]){
            scale_low  += history[i];
            rate_low   += history[i+1];
            scale_high += history[i+2];
            rate_high  += history[i+3];
        }else{
            scale_low  += history[i+2];
            rate_low   += history[i+3];
            scale_high += history[i];
            rate_high  += history[i+1];
        }
    }
    double num_samples = (double) (history.size() / 4);
    rate_low   /= num_samples;
    scale_low  /= num_samples;
    rate_high  /= num_samples;
    scale_high /= num_samples;

    // High dimensions and multiple modes reduce the acceptance rate,
    // and tuning sampling parameters such as the differential update magnitude could have both positive and negative effects.
    // Low acceptance rate is undesirable (called poor mixing in literature)
    // and extremely low rate indicates ill-posed or incorrectly implemented problem.
    cout << "Acceptance rate: " << std::fixed << state.getAcceptanceRate() << "\n\n";
    cout << "High dimensions and multiple modes reduce the acceptance rate, and sampling parameters (e.g., differential update magnitude), aff.\n\n";
    cout << "Inferred values (noise free case):" << endl;
    cout << " low   rate:" << setw(12) << std::fixed << rate_low   << "   error:" << setw(12) << std::scientific << std::abs(rate_low - 1.0)   << endl;
    cout << " low  scale:" << setw(12) << std::fixed << scale_low  << "   error:" << setw(12) << std::scientific << std::abs(scale_low - 1.0)  << "\n" << endl;
    cout << " high  rate:" << setw(12) << std::fixed << rate_high  << "   error:" << setw(12) << std::scientific << std::abs(rate_high - 3.0)  << endl;
    cout << " high scale:" << setw(12) << std::fixed << scale_high << "   error:" << setw(12) << std::scientific << std::abs(scale_high - 0.4) << "\n" << endl;

    cout << "\n" << "-------------------------------------------------------------------------------------------------" << endl;

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_04 example]
#endif
}
