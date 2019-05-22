#include <iostream>
#include <iomanip>
#include <ctime>

#include "TasmanianDREAM.hpp"

using namespace std;

//! \file example_dream_01.cpp
//! \brief Examples for the Tasmanian DREAM module.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAMExamples
//!
//! Tasmanian DREAM Example 1

/*!
 * \ingroup TasmanianDREAMExamples
 * \addtogroup TasmanianDREAMExamples1 Tasmanian DREAM module, example 1
 *
 * Example 1:
 * Demonstrates how to make a custom probability distribution and use Tasmanian DREAM
 * to generate random samples.
 */

//! \brief DREAM Example 1: sample from a custom defined probability distribution.
//! \ingroup TasmanianDREAMExamples1

//! \snippet DREAM/Examples/example_dream_01.cpp DREAM_Example_01 example
void dream_example_01(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_01 example]
#endif
    // using the default random engine, but must reset the random number generator
    srand((int) time(nullptr));

    // Example 1:
    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 1: make your own probability distribution" << endl;
    cout << "           sample from Gaussian distribution: f(x) = exp(-x^2)" << endl;
    cout << "           ingnoring scaling constants, using 3000 smaples" << endl;
    cout << "    See the comments in example_dream_01.cpp" << endl << endl;

    int num_dimensions = 1, num_chains = 30;
    int num_burnup_iterations = 200;
    int num_collect_iterations = 1000; // 1000 iterations with 30 chains gives 3000 samples

    TasDREAM::TasmanianDREAM state(num_chains, num_dimensions);

    // need to initialize the chains
    std::vector<double> initial_chains;
    TasDREAM::genUniformSamples({-2.0}, {2.0}, num_chains, initial_chains); // create uniform samples in (-2, 2)
    state.setState(initial_chains);

    // Main call to Tasmanian DREAM Sampling algorithm
    TasDREAM::SampleDREAM(num_burnup_iterations, num_collect_iterations,
                          [&](const std::vector<double> &candidates, std::vector<double> &values){ // use lambda to implement the formula
                              std::transform(candidates.begin(), candidates.end(), values.begin(),
                                             [&](double x)->double{ return std::exp(-x*x); }); // implement the formula
                          },
                          [&](const std::vector<double>&)->bool{ return true; }, // unbound domain
                          TasDREAM::dist_uniform, 0.5, // uniform independent update of magnitude 0.5
                          state,
                          TasDREAM::const_percent<90> // use 90% of differential update
                         );

    // compute the mean and variance of the samples
    std::vector<double> mean, variance;
    state.getHistoryMeanVariance(mean, variance);

    cout << "Using regular form:" << endl;
    cout << "       mean:" << setw(13) << std::fixed << mean[0]
         << "   error:" << setw(12) << std::scientific << std::abs(mean[0]) << endl;
    cout << "   variance:" << setw(13) << std::fixed << variance[0]
         << "   error:" << setw(12) << std::scientific << std::abs(variance[0] - 0.5) << endl;


    // Repeat the same experiment, but using log-form
    state = TasDREAM::TasmanianDREAM(num_chains, num_dimensions); // reset the state
    state.setState(initial_chains); // set the initial state

    // sample again, but use the logarithm form of the formula
    TasDREAM::SampleDREAM<TasDREAM::logform>
                         (num_burnup_iterations, num_collect_iterations,
                          [&](const std::vector<double> &candidates, std::vector<double> &values){
                              std::transform(candidates.begin(), candidates.end(), values.begin(),
                                             [&](double x)->double{ return -x*x; }); // implement the logarithm of the formula
                          },
                          [&](const std::vector<double>&)->bool{ return true; }, // unbound domain
                          TasDREAM::dist_uniform, 0.5, // uniform independent update of magnitude 0.5
                          state,
                          TasDREAM::const_percent<90> // use 90% of differential update
                         );

    // get the mean and variance for the logform samples
    state.getHistoryMeanVariance(mean, variance);

    cout << "Using regular form:" << endl;
    cout << "       mean:" << setw(13) << std::fixed << mean[0]
         << "   error:" << setw(12) << std::scientific << std::abs(mean[0]) << endl;
    cout << "   variance:" << setw(13) << std::fixed << variance[0]
         << "   error:" << setw(12) << std::scientific << std::abs(variance[0] - 0.5) << endl;

    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_01 example]
#endif
}
