#include <iostream>
#include <iomanip>
#include <ctime>
#include "math.h"

#include "TasmanianDREAM.hpp"

using namespace std;

//! \file example_dream_03.cpp
//! \brief Examples for the Tasmanian DREAM module.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAMExamples
//!
//! Tasmanian DREAM Example 3

/*!
 * \ingroup TasmanianDREAMExamples
 * \addtogroup TasmanianDREAMExamples3 Tasmanian DREAM module, example 3
 *
 * Example 3:
 * Given data that is the superposition of two sin-waves,
 * use Sparse Grid and Bayesian inference to identify the frequencies and shifts of each wave.
 * Higher dimensions and bi-modal posterior will decrease the acceptance rate,
 * thus we need more samples than the previous example.
 */

//! \brief DREAM Example 3: signal decomposition, using bi-modal posterior distribution.
//! \ingroup TasmanianDREAMExamples3

//! \snippet DREAM/Examples/example_dream_03.cpp DREAM_Example_03 example
void dream_example_03(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_03 example]
#endif
    // using the default random engine, but must reset the random number generator
    srand((int) time(nullptr));

    // EXAMPLE 3:
    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
    cout << "EXAMPLE 3: set inference problem, identify x_0 and x_1 model parameters from data (noise free example)" << endl;
    cout << "           model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(5*M_PI*t + 0.3*M_PI) + sin(10*M_PI*t + 0.1*M_PI)" << endl;
    cout << "           compared to Example 2, the data is superposition of two signals and the posterior is multi-modal" << endl;
    cout << "             -- problem parameters --" << endl;
    cout << "           t in [0,1], discretized with 32 equidistant nodes" << endl;
    cout << "           likelihood is exp(- 16 * (f(x) - d)^2)" << endl;
    cout << "           use sparse grid to interpolate the model" << endl;
    cout << "     NOTE: 16 corresponds to discretization error in t (divided by 2)" << endl << endl;

    // multi-modal distributions require a lot more chains and samples
    int num_chains = 500;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 300;
    // the total number of samples is num_chains * num_iterations
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

    // create data which is a superposition of two signals
    std::vector<double> signal1(num_discrete_nodes), signal2(num_discrete_nodes), data(num_discrete_nodes);
    model( 5.0, 0.3 * M_PI, signal1);
    model(10.0, 0.1 * M_PI, signal2);
    std::transform(signal1.begin(), signal1.end(), signal2.begin(), data.begin(), std::plus<double>());

    TasGrid::TasmanianSparseGrid grid;
    grid.makeSequenceGrid(2, num_discrete_nodes, 30, TasGrid::type_iptotal, TasGrid::rule_leja); // 30-th order polynomial
    std::vector<double> domain_a = { 1.0, -0.1}; // set search interval x_0 in [1, 12], x_1 in [-0.1, 1.7]
    std::vector<double> domain_b = {12.0,  1.7};
    grid.setDomainTransform(domain_a, domain_b);

    std::vector<double> points; // get the sparse grid points
    grid.getNeededPoints(points);
    std::vector<double> values; // stores the model values

    for(int i=0; i<grid.getNumNeeded(); i++){ // compute the model for each sparse grid point
        std::vector<double> model_at_point(num_discrete_nodes);
        model(points[2*i], points[2*i+1], model_at_point);
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
    TasDREAM::SampleDREAMPost<TasDREAM::logform>
                              (num_burnup_iterations, num_sample_iterations,
                               likely, // provide the likelihood
                               grid, // provide the model
                               TasDREAM::uniform_prior, // assume non-informative prior
                               TasDREAM::dist_gaussian, 0.01, // Gaussian independent update of magnitude 0.01
                               state,
                               TasDREAM::const_percent<90> // use 90% of differential update
                              );

    // get the vector containing the sampling history (could also use "auto history = state.getHistory();"
    const std::vector<double> &history = state.getHistory();

    // splitting the bi-modal history is tricky, use the mid-point of the domain
    // frequencies below 6.5 will be added to low frequency signal
    // frequencies above 6.5 will be added to high frequency signal
    double frequency_low = 0.0, correction_low = 0.0;
    double frequency_high = 0.0, correction_high = 0.0;
    int num_low = 0, num_high = 0;
    for(size_t i=0; i<history.size(); i+=2){ // 2 parameters, the samples in a stride of 2
        if (history[i] < 6.5){
            frequency_low += history[i];
            correction_low += history[i+1];
            num_low++;
        }else{
            frequency_high += history[i];
            correction_high += history[i+1];
            num_high++;
        }
    }
    // optional report the number of samples used to compute the low and high frequencies and the MCMC acceptance rate.
    //cout << "low samples = " << num_low << "  high samples = " << num_high << ", acceptance rate = " << state.getAcceptanceRate() << endl;
    frequency_low /= (double)(num_low);
    correction_low /= (double)(num_low);
    frequency_high /= (double)(num_high);
    correction_high /= (double)(num_high);

    cout.precision(5);
    cout << "Inferred values:" << endl;
    cout << " low   frequency:" << setw(12) << std::fixed << frequency_low
         << "   error:" << setw(12) << std::scientific << fabs(frequency_low - 5.0) << endl;
    cout << " low  correction:" << setw(12) << std::fixed << correction_low
         << "   error:" << setw(12) << std::scientific << fabs(correction_low - 0.3 * M_PI) << endl << endl;
    cout << " high  frequency:" << setw(12) << std::fixed << frequency_high
         << "   error:" << setw(12) << std::scientific << fabs(frequency_high - 10.0) << endl;
    cout << " high correction:" << setw(12) << std::fixed << correction_high
         << "   error:" << setw(12) << std::scientific << fabs(correction_high - 0.1 * M_PI) << endl << endl;

    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_03 example]
#endif
}
