/*
 * Copyright (c) 2017, Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
 *    and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
 * THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
 * THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
 * IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TASMANIAN_TASDREAM_EXTERNAL_TESTS_HPP
#define __TASMANIAN_TASDREAM_EXTERNAL_TESTS_HPP

#include "TasmanianDREAM.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

using namespace TasDREAM;

//! \internal
//! \defgroup TasDREAMTesting: Testing for the Tasmanian DREAM module.
//!
//! \par Testing
//! A series of tests covering sampling under different circumstances,
//! sample from arbitrary distribution, perform Bayesian inference on
//! a model (custom or defined by a sparse grid), or solve an optimization
//! problem.

//! \internal
//! \brief Allows to select a specific category for testing.
//! \ingroup TasDREAMTesting
enum TypeDREAMTest{
    //! \brief Execute all tests.
    test_all,
    //! \brief Tests for sampling from known probability density.
    test_analytic,
    //! \brief Tests for sampling from posterior distributions.
    test_posterior
};

//! \internal
//! \brief Report Pass/FAIL bases on pass, format is "name Pass/FAIL".
//! \ingroup TasDREAMTesting
inline void reportPassFail(bool pass, const char *name){
    cout << setw(45) << name << setw(15) << ((pass) ? "Pass" : "FAIL") << endl;
}

//! \internal
//! \brief Report Pass/FAIL bases on pass, format is "name variant Pass/FAIL"
//! \ingroup TasDREAMTesting
inline void reportPassFail(bool pass, const char *name, const char *variant){
    cout << setw(15) << name << setw(30) << variant << setw(15) << ((pass) ? "Pass" : "FAIL") << endl;
}

//! \internal
//! \brief Get a random number to be used for a random seed, uses \b std::time(\b nullptr).
//! \ingroup TasDREAMTesting
inline long unsigned getRandomRandomSeed(){ return static_cast<long unsigned>(std::time(nullptr)); }

//! \internal
//! \brief Dump the history of the state to \b cout using one point per line (debug use only).
//! \ingroup TasDREAMTesting
inline void printHistory(const TasmanianDREAM &state){
    const std::vector<double> &hist = state.getHistory();
    const std::vector<double> &vals = state.getHistoryPDF();
    int num_dimensions = state.getNumDimensions();
    auto ih = hist.begin();
    cout << std::scientific; cout.precision(6);
    for(auto v : vals){
        for(int i=0; i<num_dimensions; i++) cout << *ih++ << "  ";
        cout << v << endl;
    }
}

//! \internal
//! \brief Print the mean and standard deviation of the history (debug use only).
//! \ingroup TasDREAMTesting
inline void printStats(const TasmanianDREAM &state, const char *message = nullptr){
    std::vector<double> mean, variance;
    state.getHistoryMeanVariance(mean, variance);
    if (message != nullptr) cout << message << endl;
    cout << std::scientific; cout.precision(6);
    for(auto m : mean) cout << m << "  ";
    cout << endl;
    for(auto v : variance) cout << std::sqrt(v) << "  ";
    cout << endl;
}

//! \internal
//! \brief Print the approximate mode from the history.
//! \ingroup TasDREAMTesting
inline void printMode(const TasmanianDREAM &state, const char *message = nullptr){
    std::vector<double> mode;
    state.getApproximateMode(mode);
    if (message != nullptr) cout << message << endl;
    cout << std::scientific; cout.precision(6);
    for(auto m : mode) cout << m << "  ";
    cout << endl;
}

//! \internal
//! \brief Simple model of a signal with two overlapping frequencies.
//! \ingroup TasDREAMTesting

//! The model is \f$ f(t) = \sin(\pi t) + M \sin( F \pi t) \f$ where M is the \b magnitude and F is the \b frequency.
//! The result is sampled for t in (\b time_step, \b num_steps * \b time_step) resulting in a vector of length \b num_steps.
//! The entries are written to the array \b y.
inline void getSinSinModel(double magnitude, double frequency, double time_step, int num_steps, double *y){
    double t = time_step;
    for(int i=0; i<num_steps; i++){
        y[i] = std::sin(DreamMaths::pi * t) + magnitude * std::sin(frequency * DreamMaths::pi * t);
        t += time_step;
    }
}

//! \internal
//! \brief Tester class, manages general test parameters (e.g., verbose mode) and runs the individual tests.
//! \ingroup TasDREAMTesting

//! The tester class that calls individual tests and manages common testing parameters.
class DreamExternalTester{
public:
    //! \brief Default constructor.
    DreamExternalTester() : verbose(false), showvalues(false), usetimeseed(false){}
    //! \brief Default destructor.
    ~DreamExternalTester(){}

    //! \brief Enable verbose mode, show information for each individual test.
    void showVerbose(){ verbose = true; }
    //! \brief Show the statistics and critical values for each test (more verbose).
    void showStatsValues(){ showvalues = true; }
    //! \brief Use \b srand() with \b time(\b nullptr) and a random number for the random seed (as opposed to the hardcoded random value).
    void useRandomRandomSeed(){ usetimeseed = true; }

    //! \brief Main test driver, selects the test to run and writes begin and end information, returns \b false if any test fails.
    bool performTests(TypeDREAMTest test);

protected:
    //! \brief Perform test for sampling from known distributions.
    bool testKnownDistributions();

    //! \brief Generate 3D Gaussian samples using DREAM.
    bool testGaussian3D();

    //! \brief Generate 2D Gaussian samples using DREAM and Sparse Grids.
    bool testGaussian2D();

    //! \brief Perform test for sampling from inferred posterior distributions.
    bool testPosteriorDistributions();

    //! \brief Generate samples from custom models.
    bool testCustomModel();

    //! \brief Generate samples from sparse grid model.
    bool testGridModel();

    //! \brief Hardcoded table with chi-squared values.
    double getChiValue(size_t num_degrees);

    //! \brief Performs Pearson's chi-squared test whether the two cell counts describe the same distribution (returns \b false if they differ above the 99-th percentile).

    //! Note: theoretically, the test will fail in 1 of 100 cases due to sheer randomness, it happens more often if weak pseudo-random generator is used.
    //! Performing multiple tests across different compilers will sooner or later lead to a false-positive.
    //! Hardcoding the random seeds and pseudo-random distributions is not enough due to round-off error generated by compiler optimizations.
    //!
    //! The actual test statistics is: \f$ X = \sum_{i} \frac{A_i \sqrt{B/A} - B_i \sqrt{A/B}}{A_i + B_i} \f$, where \f$ A = \sum_{i} A_i \f$ and \f$ B = \sum_{i} B_i \f$.
    //! The test statistics has chi-squared distribution with degrees of freedom one less than the total number of cells.
    bool testFit(const std::vector<int> &cell_count_a, const std::vector<int> &cell_count_b);

    //! \brief Bin hypercube data.

    //! Takes \b data that represents vectors of dimensions \b lower.size() that lay in hypercube defined by \b lower and \b upper limits,
    //! the data is binned into uniform cell grid with 1-D size \b num_bins1D.
    void binHypercubeSamples(const std::vector<double> &lower, const std::vector<double> &upper, int num_bins1D, const std::vector<double> &data, std::vector<int> &bin_count);

    //! \brief Test if sample follow the same distribution, calls \b binHypercubeSamples() and \b testFit().
    bool compareSamples(const std::vector<double> &lower, const std::vector<double> &upper, int num_bins1D, const std::vector<double> &data1, const std::vector<double> &data2);

private:
    bool verbose, showvalues, usetimeseed;
};


//! \internal
//! \brief Normally an empty function, but the user can put code inside for quick testing.

//! Performance, debug, and other tests of a library require that client code is written and linked to the library,
//! which require a separate project with source files, build system, etc.
//! This function offers an easy way to compile additional code using the default Tasmanian build engine,
//! the code can be invoked from the root of the build folder using:
//! \code ./DREAM/dreamtest -debug \endcode
//! The function is implemented at the very bottom of \b tasdreamExternalTests.cpp.
void testDebug();


#endif
