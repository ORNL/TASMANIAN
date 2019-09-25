#include "Tasmanian.hpp"

using namespace std;

/*!
 * \internal
 * \file example_dream_03.cpp
 * \brief Examples for the Tasmanian DREAM module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAMExamples
 *
 * Tasmanian DREAM Example 3
 * \endinternal
 */

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
    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
    cout << "EXAMPLE 3: set the inference problem: identify x_0 and x_1 model parameters\n"
         << "           from data (noise free example)\n"
         << "           model: f(x) = sin(x_0*M_PI*t + x_1),\n"
         << "           data: d = sin(5*M_PI*t + 0.3*M_PI) + sin(10*M_PI*t + 0.1*M_PI)\n"
         << "           compared to Example 2, the data is a superposition of two signals\n"
         << "           and the posterior is multi-modal\n"
         << "             -- problem setup --\n"
         << "           t in [0,1], t is discretized with 32 equidistant nodes\n"
         << "           the likelihood is exp(- 16 * (f(x) - d)^2)\n"
         << "           using a sparse grid to interpolate the model\n"
         << "     NOTE: 16 = 32/2 corresponds to the discretization error in t\n\n";

    constexpr double pi = 3.14159265358979323846;

    // multi-modal distributions require a lot more chains and samples
    int num_chains = 500;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 300;
    // the total number of samples is num_chains * num_iterations
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

    // create data which is a superposition of two signals
    std::vector<double> signal1(num_discrete_nodes),
                        signal2(num_discrete_nodes),
                        data(num_discrete_nodes);
    model( 5.0, 0.3 * pi, signal1);
    model(10.0, 0.1 * pi, signal2);
    std::transform(signal1.begin(), signal1.end(), signal2.begin(),
                   data.begin(), std::plus<double>());

    auto grid = TasGrid::makeSequenceGrid(2, num_discrete_nodes, 30, // 30-th order polynomial
                                          TasGrid::type_iptotal, TasGrid::rule_leja);
    // set the search interval x_0 in [1, 12], x_1 in [-0.1, 1.7]
    std::vector<double> domain_a = { 1.0, -0.1};
    std::vector<double> domain_b = {12.0,  1.7};
    grid.setDomainTransform(domain_a, domain_b);

    TasGrid::loadNeededPoints<TasGrid::mode_sequential>(
        [&](std::vector<double> const &x, std::vector<double> &y, size_t)->
        void{
            model(x[0], x[1], y);
        },
        grid, 1);

    // when working with grids with large N, acceleration becomes very useful
    // if BLAS is enabled on compile time, the grid will use BLAS by default
    // GPU acceleration can be enabled here using:
    //   grid.enableAcceleration(TasGrid::accel_gpu_cuda);
    //   grid.setGPUID(0);

    // define the likelihood function and load the data
    // even though the example is noise free,
    // we assume that the "noise" is due to discretization error
    // using piece-wise constant approximation in t, the error is 1.0 / num_discrete_nodes
    TasDREAM::LikelihoodGaussIsotropic likely(1.0 / ((double) num_discrete_nodes), data);

    TasDREAM::TasmanianDREAM state(num_chains, grid); // get the dimensions from the sparse grid
    // use initial chains that are distributed uniformly over the domain
    state.setState(TasDREAM::genUniformSamples(domain_a, domain_b, num_chains));

    // Call to Tasmanian DREAM Sampling algorithm
    TasDREAM::SampleDREAM<TasDREAM::logform>
                         (num_burnup_iterations, num_sample_iterations,
                          TasDREAM::posterior<TasDREAM::logform>
                                (grid,   // provide the surrogate model
                                 likely, // provide the likelihood
                                 TasDREAM::uniform_prior), // assume non-informative prior
                          grid.getDomainInside(),
                          state,
                          TasDREAM::dist_gaussian, 0.01, // independent update of magnitude 0.01
                          TasDREAM::const_percent<90> // use 90% differential update
                          );

    // get the vector containing the sampling history
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
    // optional report the number of samples used to compute the low and high frequencies
    // and the MCMC acceptance rate.
    //cout << "low samples = " << num_low << "  high samples = " << num_high
    //     << ", acceptance rate = " << state.getAcceptanceRate() << "\n";
    frequency_low /= (double)(num_low);
    correction_low /= (double)(num_low);
    frequency_high /= (double)(num_high);
    correction_high /= (double)(num_high);

    cout.precision(5);
    cout << "Inferred values:\n"
         << " low   frequency:" << setw(12) << std::fixed << frequency_low
         << "   error:" << setw(12) << std::scientific << std::abs(frequency_low - 5.0) << "\n"
         << " low  correction:" << setw(12) << std::fixed << correction_low
         << "   error:" << setw(12) << std::scientific << std::abs(correction_low - 0.3 * pi)
         << "\n\n"
         << " high  frequency:" << setw(12) << std::fixed << frequency_high
         << "   error:" << setw(12) << std::scientific << std::abs(frequency_high - 10.0) << "\n"
         << " high correction:" << setw(12) << std::fixed << correction_high
         << "   error:" << setw(12) << std::scientific << std::abs(correction_high - 0.1 * pi)
         << "\n\n";

    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_03 example]
#endif
}
