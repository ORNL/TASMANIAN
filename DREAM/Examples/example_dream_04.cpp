#include "Tasmanian.hpp"

using namespace std;

/*!
 * \internal
 * \file example_dream_04.cpp
 * \brief Examples for the Tasmanian DREAM module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAMExamples
 *
 * Tasmanian DREAM Example 4
 * \endinternal
 */

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
    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 4: similar to Example 3, but the data is noisy and the multi-modal posterior comes\n"
         << "           from a symmetry in the model\n"
         << "           set inference problem: identify x_0, x_1, x_2 and x_3 (four) model parameters\n"
         << "                                  from noisy data\n"
         << "           model: f(x) = x_0*exp(x_1*t) + x_2*exp(x_3*t),  data: d = exp(t) + 0.4*exp(3*t)\n"
         << "           note the symmetry between x_1 and x_3 which results in a bi-modal posterior\n"
         << "           t in [0,1], t is discretized with 64 equidistant nodes\n"
         << "           likelihood is exp(-64 * (f(x) - d)^2)\n"
         << "           using a sparse grid to interpolate the model\n"
         << "     NOTE: 64 corresponds to the discretization error and the added noise\n" << endl;

    // multi-modal distributions require a lot more chains and samples
    int num_chains = 50;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 1000;
    // the total number of samples is num_chains * num_iterations
    int num_discrete_nodes = 64;

    // create a lambda function that represents the model
    // normally this would be a call to an external code
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

    auto grid = TasGrid::makeSequenceGrid(4, num_discrete_nodes, 15, // 15-th order polynomial
                                          TasGrid::type_iptotal, TasGrid::rule_leja);
    // set search interval x_0, x_2 in [0.2, 1.2], x_1, x_3 in [0.5, 4.0]
    std::vector<double> domain_a = {0.2, 0.5, 0.2, 0.5};
    std::vector<double> domain_b = {1.2, 4.0, 1.2, 4.0};
    grid.setDomainTransform(domain_a, domain_b);

    TasGrid::loadNeededValues<TasGrid::mode_sequential>(
        [&](std::vector<double> const &x, std::vector<double> &y, size_t)->
        void{
            model(x[0], x[1], x[2], x[3], y);
        },
        grid, 1);

    // when working with grids with large N, acceleration becomes very useful
    // if BLAS is enabled on compile time, the grid will use BLAS by default
    // GPU acceleration can be enabled here using:
    //   grid.enableAcceleration(TasGrid::accel_gpu_cuda);
    //   grid.setGPUID(0);

    // define the likelihood function and load the data
    // the discretization error and noise are both is 1.0 / num_discrete_nodes
    // the Gaussian formula adds another factor of 0.5 which cancels one
    TasDREAM::LikelihoodGaussIsotropic likely(1.0 / ((double) num_discrete_nodes), data);

    TasDREAM::TasmanianDREAM state(num_chains, grid); // get the dimensions from the sparse grid
    std::vector<double> initial_chains;
    TasDREAM::genUniformSamples(domain_a, domain_b, num_chains, initial_chains);
    state.setState(initial_chains); // use chains distributed uniformly over the domain

    // Call to Tasmanian DREAM Sampling algorithm
    TasDREAM::SampleDREAM<TasDREAM::logform>
                         (num_burnup_iterations, num_sample_iterations,
                          TasDREAM::posterior<TasDREAM::logform>
                                (grid,   // provide the model
                                 likely, // provide the likelihood
                                 TasDREAM::uniform_prior), // assume non-informative prior
                          grid.getDomainInside(),
                          state,
                          TasDREAM::dist_gaussian, 0.01, // independent update of magnitude 0.01
                          TasDREAM::const_percent<100> // use 100% of differential update
                          );

    // get the vector containing the sampling history
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
    // and tuning sampling parameters such as the differential update magnitude
    // could have either positive or negative effects.
    // Low acceptance rate is undesirable (called poor mixing in literature)
    // and extremely low rate indicates ill-posed or incorrectly implemented problem.
    cout << "Acceptance rate: " << std::fixed << state.getAcceptanceRate() << "\n\n";
    cout << "High dimensions and multiple modes reduce the acceptance rate,\n"
         << "and sampling parameters (e.g., differential update magnitude) can affect it either way.\n\n";
    cout << "Inferred values (noise free case):\n";
    cout << " low   rate:" << setw(12) << std::fixed << rate_low
         << "   error:" << setw(12) << std::scientific << std::abs(rate_low - 1.0)   << "\n";
    cout << " low  scale:" << setw(12) << std::fixed << scale_low
         << "   error:" << setw(12) << std::scientific << std::abs(scale_low - 1.0)  << "\n\n";
    cout << " high  rate:" << setw(12) << std::fixed << rate_high
         << "   error:" << setw(12) << std::scientific << std::abs(rate_high - 3.0)  << "\n";
    cout << " high scale:" << setw(12) << std::fixed << scale_high
         << "   error:" << setw(12) << std::scientific << std::abs(scale_high - 0.4) << "\n\n";

    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_04 example]
#endif
}
