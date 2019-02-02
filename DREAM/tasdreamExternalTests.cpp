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

#ifndef __TASMANIAN_TASDREAM_EXTERNAL_TESTS_CPP
#define __TASMANIAN_TASDREAM_EXTERNAL_TESTS_CPP

#include "tasdreamExternalTests.hpp"

double DreamExternalTester::getChiValue(size_t num_degrees){
    switch(num_degrees){
        case   9: return  21.666;
        case  15: return  30.578;
        case  19: return  36.191;
        case  20: return  37.566;
        case  24: return  42.980;
        case  49: return  74.919;
        case  99: return 134.642;
        case 124: return 163.546;
        default:
            throw std::runtime_error("ERROR: Unknown degrees of freedom for the Chi-squared test.");
    }
}

bool DreamExternalTester::testFit(const std::vector<int> &cell_count_a, const std::vector<int> &cell_count_b){
    double suma = (double) std::accumulate(cell_count_a.begin(), cell_count_a.end(), 0);
    double sumb = (double) std::accumulate(cell_count_b.begin(), cell_count_b.end(), 0);

    double scale = sqrt(sumb / suma);

    double test_value = 0.0;

    auto ia = cell_count_a.begin(), ib = cell_count_b.begin();
    while(ia != cell_count_a.end()){
        double diff = ((double) *ia) * scale - ((double) *ib) / scale;
        double sum = (double) (*ia++ + *ib++);
        if (sum > 0.0) test_value += diff * diff / sum;
    }

    bool pass = (test_value < getChiValue(cell_count_a.size() - 1));
    if (!pass || showvalues){
        if (!pass) cout << "Chi-Squared test FAILED" << endl;
        cout << "Totals: " << suma << "  " << sumb << endl;
        cout << "Chi-Squared test value = " << test_value << " num cells: " << cell_count_a.size() << endl;
        cout << "Critical Chi-Squared value = " << getChiValue(cell_count_a.size() - 1) << endl;
    }

    return pass;
}

void DreamExternalTester::binHypercubeSamples(const std::vector<double> &lower, const std::vector<double> &upper, int num_bins1D, const std::vector<double> &data, std::vector<int> &bin_count){
    size_t num_dimensions = lower.size();
    if (upper.size() != num_dimensions) throw std::runtime_error("ERROR: upper and lower must have the same size in binHypercubeSamples() DREAM testing");

    std::vector<double> dx(num_dimensions);
    auto il = lower.begin(), iu = upper.begin();
    for(auto &d : dx) d = (*iu++ - *il++) / ((double) num_bins1D);

    size_t num_bins = 1;
    for(size_t i=0; i<num_dimensions; i++) num_bins *= num_bins1D;
    bin_count = std::vector<int>(num_bins, 0);

    auto id = data.begin();
    while(id != data.end()){
        std::vector<size_t> binid(num_dimensions);
        il = lower.begin();
        iu = dx.begin();
        for(auto &i : binid){
            i = (size_t) ((*id++ - *il++) / *iu++);
            if (i >= (size_t) num_bins1D) i = num_bins1D-1;
        }
        size_t bin_index = 0;
        for(auto i : binid) bin_index = num_dimensions * bin_index + i;
        bin_count[bin_index]++;
    }
}

bool DreamExternalTester::compareSamples(const std::vector<double> &lower, const std::vector<double> &upper, int num_bins1D,
                                         const std::vector<double> &data1, const std::vector<double> &data2){
    std::vector<int> count1, count2;
    binHypercubeSamples(lower, upper, num_bins1D, data1, count1);
    binHypercubeSamples(lower, upper, num_bins1D, data2, count2);
    return testFit(count1, count2);
}

bool DreamExternalTester::testGaussian3D(){
    bool passAll = true;
    int num_dimensions = 3;
    int rseed = 42, num_samples = 1000, num_chains = 20;
    int num_iterations = num_samples / num_chains + 2;
    int num_burnup = 20 * num_iterations;
    if (usetimeseed) rseed = getRandomRandomSeed();

    std::minstd_rand park_miller;
    park_miller.seed(rseed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    
    // compute reference samples, mean 2.0, std 3.0
    std::vector<double> means(num_dimensions, 2.0), deviations(num_dimensions, 3.0);
    std::vector<double> tresult;
    genGaussianSamples(means, deviations, num_samples, tresult, [&]()->double{ return unif(park_miller); });

    // Use DREAM with zero-weight (i.e., standard Metropolis-Hastings)
    TasmanianDREAM state(num_chains, num_dimensions);
    std::vector<double> initial_state;
    genGaussianSamples(means, deviations, num_chains, initial_state, [&]()->double{ return unif(park_miller); }); // initialize with correct mean 2.0, std 3.0
    state.setState(initial_state);

    SampleDREAM(num_burnup, num_iterations,
        [&](const std::vector<double> &candidates, std::vector<double> &values){
            // 3D Gaussian PDF with standard deviation of 3.0
            auto ix = candidates.begin();
            for(auto &v : values)
                v = getDensity<dist_gaussian>(*ix++, 2.0, 9.0) * getDensity<dist_gaussian>(*ix++, 2.0, 9.0) * getDensity<dist_gaussian>(*ix++, 2.0, 9.0);
        },
        [&](const std::vector<double>&)->bool{ return true; }, // unbounded domain
        [&](std::vector<double> &x){
            applyGaussianUpdate(x, 3.0, [&]()->double{ return unif(park_miller); });
        },
        state,
        const_percent<0>, // independent chains, no differential proposal
        [&]()->double{ return unif(park_miller); }
    );

    std::vector<double> upper(num_dimensions, 11.0), lower(num_dimensions, -7.0); // compute over a box of 3 standard deviations

    bool pass = compareSamples(lower, upper, 5, tresult, state.getHistory());
    passAll = passAll && pass;
    if (verbose || !pass) reportPassFail(pass, "Gaussian 3D", "with independent chains");

    state = TasmanianDREAM(num_chains, num_dimensions); // reinitialize
    state.setState(initial_state);

    SampleDREAM(num_burnup, 2*num_iterations,
        [&](const std::vector<double> &candidates, std::vector<double> &values){
            // 3D Gaussian PDF with standard deviation of 3.0
            auto ix = candidates.begin();
            for(auto &v : values)
                v = getDensity<dist_gaussian>(*ix++, 2.0, 9.0) * getDensity<dist_gaussian>(*ix++, 2.0, 9.0) * getDensity<dist_gaussian>(*ix++, 2.0, 9.0);
        },
        lower, upper, // large domain
        dist_uniform, 0.2, // uniform proposal
        state,
        const_percent<50>, // differential proposal is weighted by 50%
        [&]()->double{ return unif(park_miller); }
    );

    pass = compareSamples(lower, upper, 5, tresult, state.getHistory());
    passAll = passAll && pass;

    if (verbose || !pass) reportPassFail(pass, "Gaussian 3D", "with correlated chains");

    reportPassFail(passAll, "Gaussian 3D", "DREAM vs Box-Muller");

    return passAll;
}

bool DreamExternalTester::testGaussian2D(){
    bool passAll = true;
    int num_dimensions = 2;
    int rseed = 42, num_samples = 1000, num_chains = 20;
    int num_iterations = num_samples / num_chains + 2;
    int num_burnup = 20 * num_iterations;
    if (usetimeseed) rseed = getRandomRandomSeed();

    std::minstd_rand park_miller;
    park_miller.seed(rseed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    // compute reference samples, mean 0.3, std 0.15 (3 deviations fit in [-1, 1]^2)
    std::vector<double> tresult(num_dimensions * num_samples, 0.3);
    applyGaussianUpdate(tresult, 0.15, [&]()->double{ return unif(park_miller); });

    // approximate the pdf in log-form, log-form of the Gaussian pdf is quadratic, the grid gives exact match
    TasGrid::TasmanianSparseGrid grid;
    grid.makeSequenceGrid(2, 1, 2, TasGrid::type_iptotal, TasGrid::rule_rleja); // interpolates exactly all quadratic polynomials
    std::vector<double> grid_points, values;
    grid.getNeededPoints(grid_points);
    values.resize(grid_points.size() / 2);
    auto ip = grid_points.begin();
    for(auto &v : values)
        v = getDensity<dist_gaussian, logform>(*ip++, 0.3, 0.0225) + getDensity<dist_gaussian, logform>(*ip++, 0.3, 0.0225);
    grid.loadNeededPoints(values);

    // initialize the DREAM state
    TasmanianDREAM state(num_chains, num_dimensions);
    std::vector<double> initial_set(num_chains * num_dimensions, 0.0); // initialize with uniform samples
    applyUniformUpdate(initial_set, 1.0, [&]()->double{ return unif(park_miller); });
    state.setState(initial_set);

    SampleDREAMGrid<logform>(num_burnup, num_iterations, grid, uniform_prior,
        dist_gaussian, 0.1,
        state,
        const_percent<98>, // correlated chains
        [&]()->double{ return unif(park_miller); }
    );

    std::vector<double> upper(num_dimensions, 1.0), lower(num_dimensions, -1.0); // compute over a box of over 3 standard deviations
    bool pass = compareSamples(lower, upper, 10, tresult, state.getHistory());

    passAll = passAll && pass;
    if (verbose || !pass) reportPassFail(pass, "Gaussian 2D", "with inferred domain");


    // ------------------------------------------------------------ //
    // next test uses a sub-domain of the first quadrant, the standard deviation is smaller
    std::fill(tresult.begin(), tresult.end(), 0.3);
    applyGaussianUpdate(tresult, 0.1, [&]()->double{ return unif(park_miller); });

    // approximate the pdf in regular form, true approximation
    grid.makeSequenceGrid(2, 1, 24, TasGrid::type_iptotal, TasGrid::rule_rleja); // interpolates exactly all quadratic polynomials
    grid.getNeededPoints(grid_points);
    values.resize(grid_points.size() / 2);
    ip = grid_points.begin();
    for(auto &v : values) // using tighter variance of 0.01
        v = getDensity<dist_gaussian>(*ip++, 0.3, 0.01) * getDensity<dist_gaussian>(*ip++, 0.3, 0.01);
    grid.loadNeededPoints(values);

    // re-initialize the DREAM state
    state = TasmanianDREAM(num_chains, num_dimensions);
    initial_set = std::vector<double>(tresult.begin(), tresult.begin() + num_chains * num_dimensions);
    state.setState(initial_set);

    lower = std::vector<double>(num_dimensions, 0.0); // consider only the first quadrant
    upper = std::vector<double>(num_dimensions, 1.0);

    SampleDREAMGrid<regform>(num_burnup, num_iterations, grid, uniform_prior,
        lower, upper,
        dist_uniform, 0.1,
        state,
        const_percent<98>, // correlated chains
        [&]()->double{ return unif(park_miller); }
    );

    // check if any of the samples fall outside of the domain
    pass = compareSamples(lower, upper, 10, tresult, state.getHistory()) &&
           std::none_of(state.getHistory().begin(), state.getHistory().end(), [&](double x)->bool{ return ((x < 0.0) || (x>1.0)); });

    passAll = passAll && pass;
    if (verbose || !pass) reportPassFail(pass, "Gaussian 2D", "with custom domain");


    // ------------------------------------------------------------ //
    // next test uses the same sub-domain of the first quadrant, but the grid and prior each define different dimensions

    // approximate the pdf in regular form, true approximation
    grid.makeSequenceGrid(2, 1, 24, TasGrid::type_iptotal, TasGrid::rule_rleja); // interpolates exactly all quadratic polynomials
    grid.getNeededPoints(grid_points);
    values.resize(grid_points.size() / 2);
    ip = grid_points.begin();
    for(auto &v : values){ // using tighter variance of 0.01
        v = getDensity<dist_gaussian>(*ip++, 0.3, 0.01);
        ip++; // skip the second dimension in the likelihood
    }
    grid.loadNeededPoints(values);

    // re-initialize the DREAM state
    state = TasmanianDREAM(num_chains, num_dimensions);
    initial_set = std::vector<double>(tresult.begin(), tresult.begin() + num_chains * num_dimensions);
    state.setState(initial_set);

    SampleDREAMGrid<regform>(num_burnup, num_iterations, grid,
        [&](const std::vector<double> &candidates, std::vector<double> &vals)->void{
            auto ic = candidates.begin();
            for(auto &v : vals){ // using tighter variance of 0.01
                ic++; // skip the first dimension in the likelihood
                v = getDensity<dist_gaussian>(*ic++, 0.3, 0.01);
            }
        },
        lower, upper,
        dist_uniform, 0.1,
        state,
        const_percent<98>, // correlated chains
        [&]()->double{ return unif(park_miller); }
    );

    // check if any of the samples fall outside of the domain
    pass = compareSamples(lower, upper, 10, tresult, state.getHistory()) &&
           std::none_of(state.getHistory().begin(), state.getHistory().end(), [&](double x)->bool{ return ((x < 0.0) || (x>1.0)); });

    passAll = passAll && pass;
    if (verbose || !pass) reportPassFail(pass, "Gaussian 2D", "with custom prior");

    reportPassFail(passAll, "Gaussian 2D", "DREAM-Grid vs Box-Muller");

    return passAll;
}

bool DreamExternalTester::testKnownDistributions(){
    // Test Gaussian distribution

    bool pass1 = testGaussian3D();
    bool pass2 = testGaussian2D();

    return pass1 && pass2;
}

bool DreamExternalTester::testCustomModel(){
    bool passAll = true;
    int num_dimensions = 3;
    int rseed = 42, num_samples = 1000, num_chains = 40;
    int num_iterations = num_samples / num_chains + 2;
    int num_burnup = 20 * num_iterations;
    if (usetimeseed) rseed = getRandomRandomSeed();

    std::minstd_rand park_miller;
    park_miller.seed(rseed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    // compute reference samples, means 1.5, 2.0 and 2.5, variance 4.0, 9.0, 4.0
    std::vector<double> tresult;
    genGaussianSamples({1.5, 2.0, 2.5}, {2.0, 3.0, 2.0}, num_samples, tresult, [&]()->double{ return unif(park_miller); });

    // Use DREAM with custom model of identity (all information comes form the prior and likelihood)
    TasmanianDREAM state(num_chains, num_dimensions);
    std::vector<double> initial_state(num_chains * num_dimensions, 2.0); // initialize with random samples
    applyGaussianUpdate(initial_state, 3.0, [&]()->double{ return unif(park_miller); });
    state.setState(initial_state);

    LikelihoodGaussIsotropic likely(4.0, {1.5, 2.5});
    SampleDREAMPost(num_burnup, num_iterations, likely,
                    [&](const std::vector<double> &candidates, std::vector<double> &values)->void{ // model
                        auto ic = candidates.begin();
                        auto iv = values.begin();
                        while(iv != values.end()){ // takes the first and last parameters
                            *iv++ = *ic++;
                            ic++;
                            *iv++ = *ic++;
                        }
                    },
                    [&](const std::vector<double> &candidates, std::vector<double> &values)->void{ // prior
                        auto ic = candidates.begin() + 1; // uses the second input entries only
                        for(auto &v : values){
                            v = getDensity<dist_gaussian>(*ic, 2.0, 9.0);
                            std::advance(ic, num_dimensions);
                        }
                    },
                    [&](const std::vector<double>&)->bool{ return true; }, // unbounded domain
                    [&](std::vector<double> &x){
                        applyGaussianUpdate(x, 0.5, [&]()->double{ return unif(park_miller); });
                    },
                    state,
                    const_percent<65>,
                    [&]()->double{ return unif(park_miller); }
                );

    std::vector<double> upper(num_dimensions, 11.0), lower(num_dimensions, -7.0); // compute over a box of more than 3 standard deviations

    bool pass = compareSamples(lower, upper, 5, tresult, state.getHistory());
    passAll = passAll && pass;
    if (verbose || !pass) reportPassFail(pass, "Inference 3D", "with custom model");
    
    state = TasmanianDREAM(num_chains, num_dimensions); // reinitialize
    genUniformSamples({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, num_chains, initial_state, [&]()->double{ return unif(park_miller); });
    state.setState(initial_state);

    lower = std::vector<double>(num_dimensions, 0.0);
    upper = std::vector<double>(num_dimensions, 1.0);

    likely = LikelihoodGaussIsotropic(0.01, {0.0, 0.0});

    SampleDREAMPost<logform>(num_burnup, num_iterations, likely,
                             [&](const std::vector<double> &candidates, std::vector<double> &values)->void{ // model
                                 auto ic = candidates.begin();
                                 auto iv = values.begin();
                                 while(iv != values.end()){ // takes the first and last parameters
                                     *iv++ = 1.0 - sin(M_PI * *ic++);
                                     ic++;
                                     *iv++ = 1.0 - sin(M_PI * *ic++);
                                 }
                             },
                             uniform_prior,
                             lower, upper,
                             dist_gaussian, 0.01,
                             state,
                             const_percent<50>,
                             [&]()->double{ return unif(park_miller); }
                        );

    std::vector<double> mode;
    state.getApproximateMode(mode);
    //cout << mode[0] << "  " << mode[1] << "  " << mode[2] << endl;
    pass = ((mode[0] > 0.4) && (mode[0] < 0.6) && (mode[2] > 0.4) && (mode[2] < 0.6));
    passAll = passAll && pass;
    if (verbose || !pass) reportPassFail(pass, "Inference 3D", "optimization objective");
    
    reportPassFail(pass, "Inference 3D", "DREAM Bayesian inference");

    return passAll;
}

bool DreamExternalTester::testGridModel(){
    bool passAll = true;
    int num_dimensions = 2, num_outputs = 64;
    int rseed = 42, num_samples = 1000, num_chains = 40;
    int num_iterations = num_samples / num_chains + 2;
    int num_burnup = 20 * num_iterations;
    if (usetimeseed) rseed = getRandomRandomSeed();

    std::minstd_rand park_miller;
    park_miller.seed(rseed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    
    // Construct sparse grid approximation to the SinSin model
    std::vector<double> lower = {0.0, 2.0}, upper = {4.0, 6.0};
    TasGrid::TasmanianSparseGrid grid;
    grid.makeLocalPolynomialGrid(num_dimensions, num_outputs, 8, 2); // using quadratic basis of level 4
    grid.setDomainTransform(lower, upper); // magnitude is in range (0, 4), frequency in range (2.0, 6.0)
    std::vector<double> points, values(num_outputs * grid.getNumPoints());
    grid.getNeededPoints(points);
    
    auto ip = points.begin(), iv = values.begin();
    while(ip != points.end()){
        getSinSinModel(*ip, *(ip+1), 1.0 / ((double) num_outputs), num_outputs, &*iv);
        std::advance(ip, num_dimensions);
        std::advance(iv, num_outputs);
    }
    grid.loadNeededPoints(values); // surrogate constructed
    
    // initialize the state
    TasmanianDREAM state(num_chains, grid);
    std::vector<double> initial_state;
    genUniformSamples(lower, upper, num_chains, initial_state, [&]()->double{ return unif(park_miller); });
    state.setState(initial_state);
    
    // initialize the likelihood
    std::vector<double> data(num_outputs);
    getSinSinModel(2.0, 5.0, 1.0 / ((double) num_outputs), num_outputs, data.data()); // true magnitude 2.0, frequency 5.0
    LikelihoodGaussIsotropic likely(0.01, data);
    
    // sample using uniform prior
    SampleDREAMPost<logform>(num_burnup, num_chains, likely, grid, uniform_prior, dist_gaussian, 0.1, state, const_percent<50>, [&]()->double{ return unif(park_miller); });
    
    //printMode(state, "mode");
    std::vector<double> mode;
    state.getApproximateMode(mode);
    bool pass = ((mode[0] > 1.0) && (mode[0] < 3.0) && (mode[1] > 4.5) && (mode[1] < 5.5));
    passAll = passAll && pass;
    if (verbose || !pass) reportPassFail(pass, "Inference 2D", "grid frequency model");
    
    reportPassFail(pass, "Inference 2D", "DREAM Bayesian grid model");
    
    return passAll;
}

bool DreamExternalTester::testPosteriorDistributions(){
    // Tests using posteriors constructed from model and prior distributions

    bool pass1 = testCustomModel();
    bool pass2 = testGridModel();

    return pass1 && pass2;
}

bool DreamExternalTester::performTests(TypeDREAMTest test){
    cout << endl << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "           Tasmanian DREAM Module: Functionality Test" << endl;
    cout << "---------------------------------------------------------------------" << endl << endl;

    bool pass = true;

    std::vector<int> results(10, 1); // results for all possible tests

    if ((test == test_all) || (test == test_analytic))  results[0] = (testKnownDistributions()) ? 1 : 0;
    if ((test == test_all) || (test == test_posterior)) results[1] = (testPosteriorDistributions()) ? 1 : 0;

    pass = std::all_of(results.begin(), results.end(), [&](int i)->bool{ return (i == 1); });

    cout << endl;
    if (pass){
        cout << "---------------------------------------------------------------------" << endl;
        cout << "           All Tests Completed Successfully" << endl;
        cout << "---------------------------------------------------------------------" << endl << endl;
    }else{
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl;
        cout << "         Some Tests Have Failed" << endl;
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl << endl;
    }
    return pass;
}

void testDebug(){
    cout << "Debug Test" << endl;
    cout << "Put here testing code and call this with ./dreamtest debug" << endl;
}

/*
TestRNG::TestRNG(int seed): s(seed % 2097165){}
TestRNG::~TestRNG(){}

double TestRNG::getSample01() const{
    s *= 511;
    s += 127;
    s %= 2097165;
    return (double)(((double) s) / 2097165.0);
}

ExternalTester::ExternalTester(int num_monte_carlo) : num_mc(num_monte_carlo), rngseed(-1){
    // rngseed = -1 indicates to use the seed coded for each test
    wfirst = 15; wsecond = 30; wthird = 15;
    verbose = false;
}
ExternalTester::~ExternalTester(){}
void ExternalTester::resetRandomSeed(){ rngseed = (int) time(0); }

void ExternalTester::setVerbose(bool new_verbose){ verbose = new_verbose; }

// test functions to add, sample from existing pdfs vs. dream vs. interpolated
// compare two data sets, make a grid of cells
// chi_square: sum_i ( sqrt(S/R) R_i - sqrt(R/S) S_i )^2 / ( R_i + S_i ), R = sum_i R_i, S = sum_i S_i
// here we have cells indexed by i
// R_i are the hits in cell i by first distribution
// S_i are the hits in cell i by second distribution

// test uniform, Beta and Gamma with different parameters, 1-D, 2-D, 3-D ...
// test truncated Gaussian, too (harder test)
// test vs interpolated pdf
bool ExternalTester::Test(TestList test){
    cout << endl << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "          Tasmanian DREAM Module: Functionality Test" << endl;
    cout << "---------------------------------------------------------------------" << endl << endl;

    bool pass = true;

    bool passUniform1D = true;
    bool passBeta1D = true;
    bool passGamma1D = true;
    bool passGauss2D = true;
    bool passInferAlpha = true;

    if ((test == test_all) || (test == test_analytic)){
        passUniform1D = testUniform1D();
        passBeta1D = testBeta1D();
        passGamma1D = testGamma1D();
        passGauss2D = testGaussian2D();
    }
    if ((test == test_all) || (test == test_model)){
        passInferAlpha = testModelLikelihoodAlpha();
    }


    pass = passUniform1D && passBeta1D && passGamma1D && passGauss2D && passInferAlpha;

    cout << endl;
    if (pass){
        cout << "---------------------------------------------------------------------" << endl;
        cout << "           All Tests Completed Successfully*" << endl;
        cout << "           * this module is still in experimental phase" << endl;
        cout << "---------------------------------------------------------------------" << endl << endl;
    }else{
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl;
        cout << "         Some Tests Have Failed" << endl;
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl << endl;

        cout << "The random nature of probability distributions and Monte Carlo sampling" << endl;
        cout << "makes testing very difficult. While we work on this problem, you can" << endl;
        cout << "try running the tests again to see if the fail is persistent or just random." << endl << endl;
    }

    return pass;
}

bool ExternalTester::testUniform1D(){
    int s = (rngseed == -1) ? 12 : rngseed;
    TestRNG rng(s);

    int num_cells = 16; double delta = 2.0 / ((double) num_cells);
    int num_chains = 50;
    std::vector<double> samples_true(num_mc);
    for(int i=0; i<num_mc; i++) samples_true[i] = -1.0 + 2.0 * rng.getSample01();

    TasDREAM::TasmanianDREAM dream;
    dream.overwriteBaseUnifrom(&rng);
    dream.setProbabilityWeightFunction(&uu1D);
    dream.setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.01);
    gauss.overwriteBaseUnifrom(&rng);
    dream.setCorrectionAll(&gauss);

    std::vector<double> samples_dream;
    dream.collectSamples(3*num_mc, num_mc / num_chains, samples_dream, false);
    //for(int i=0; i<num_mc; i++) samples_dream[i] = -1.0 + 2.0 * u.getSample01();

    std::vector<int> cells_a(num_cells, 0);
    std::vector<int> cells_b(num_cells, 0);
    for(int i=0; i<num_mc; i++){
        int c = (int) floor((samples_true[i] + 1.0) / delta);
        if (c < num_cells) cells_a[c]++;
        c = (int) floor((samples_dream[i] + 1.0) / delta);
        if (c < num_cells) cells_b[c]++;
        //if ( c >= num_cells ) cout << samples_dream[i] << endl;
    }
    //for(int i=0; i<num_cells; i++) cout << cells_a[i] << "  " << cells_b[i] << endl;

    bool pass = true;

    if (testKS(num_cells, cells_a.data(), cells_b.data())){
        cout << setw(wfirst) << "Distribution" << setw(wsecond) << "Uniform 1D" << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Distribution" << setw(wsecond) << "Uniform 1D" << setw(wthird) << "FAIL" << endl;
        pass = false;
    }

    return pass;
}

bool ExternalTester::testBeta1D(){
    int s = (rngseed == -1) ? 12 : rngseed;
    TestRNG rng(s);

    int num_cells = 10; double delta = 2.0 / ((double) num_cells);
    int num_chains = 50;
    std::vector<double> samples_true(num_mc);
    TasDREAM::BetaPDF Btrue(-1.0, 1.0, 2.0, 5.0);
    Btrue.overwriteBaseUnifrom(&rng);

    for(int i=0; i<num_mc; i++) samples_true[i] = Btrue.getSample();

    TasDREAM::TasmanianDREAM dream;
    dream.overwriteBaseUnifrom(&rng);
    dream.setProbabilityWeightFunction(&distBeta1D);
    dream.setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.0009);
    gauss.overwriteBaseUnifrom(&rng);
    dream.setCorrectionAll(&gauss);

    std::vector<double> samples_dream;
    dream.collectSamples(3*num_mc, num_mc / num_chains, samples_dream, false);
    //for(int i=0; i<num_mc; i++) samples_dream[i] = -1.0 + 2.0 * u.getSample01();
    //for(int i=0; i<num_mc; i++) cout << samples_dream[i] << endl;

    std::vector<int> cells_a(num_cells, 0);
    std::vector<int> cells_b(num_cells, 0);
    for(int i=0; i<num_mc; i++){
        int c = (int) floor((samples_true[i] + 1.0) / delta);
        if (c < num_cells) cells_a[c]++;
        c = (int) floor((samples_dream[i] + 1.0) / delta);
        if (c < num_cells) cells_b[c]++;
        //if ( c >= num_cells ) cout << samples_dream[i] << endl;
    }
    //for(int i=0; i<num_cells; i++){ cout << cells_a[i] << "  " << cells_b[i] << endl; }

    bool pass = true;

    if (testKS(num_cells, cells_a.data(), cells_b.data())){
        cout << setw(wfirst) << "Distribution" << setw(wsecond) << "Beta 1D" << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Distribution" << setw(wsecond) << "Beta 1D" << setw(wthird) << "FAIL" << endl;
        pass = false;
    }

    return pass;
}

bool ExternalTester::testGamma1D(){
    int s = (rngseed == -1) ? 12 : rngseed;
    TestRNG rng(s);

    int num_cells = 20; double delta = 10.0 / ((double) num_cells);
    int num_chains = 50;
    std::vector<double> samples_true(num_mc);
    TasDREAM::GammaPDF Gtrue(-2.0, 9.0, 2.0);
    Gtrue.overwriteBaseUnifrom(&rng);

    for(int i=0; i<num_mc; i++) samples_true[i] = Gtrue.getSample();

    TasDREAM::TasmanianDREAM dream;
    dream.overwriteBaseUnifrom(&rng);
    dream.setProbabilityWeightFunction(&distGamma1D);
    dream.setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.04);
    gauss.overwriteBaseUnifrom(&rng);
    dream.setCorrectionAll(&gauss);

    double *samples_dream = dream.collectSamples(3*num_mc, num_mc / num_chains, false);
    //for(int i=0; i<num_mc; i++) samples_dream[i] = -1.0 + 2.0 * u.getSample01();
    //for(int i=0; i<num_mc; i++) cout << samples_dream[i] << endl;

    std::vector<int> cells_a(num_cells+1, 0);
    std::vector<int> cells_b(num_cells+1, 0);
    for(int i=0; i<num_mc; i++){
        int c = (int) floor((samples_true[i] + 2.0) / delta);
        if (c < num_cells) cells_a[c]++;
        if (c >= num_cells) cells_a[num_cells]++;
        c = (int) floor((samples_dream[i] + 2.0) / delta);
        if (c < num_cells) cells_b[c]++;
        if (c >= num_cells) cells_b[num_cells]++;
        if (samples_dream[i] < -2.0 ) cout << "ERROR: bad value: " << samples_dream[i] << endl;
        //if ( c >= num_cells ) cout << samples_dream[i] << endl;
    }
    //for(int i=0; i<num_cells; i++){ cout << cells_a[i] << "  " << cells_b[i] << endl; }

    bool pass = true;

    //if (testChi(num_cells+1, cells_a, cells_b)){
    if (testKS(num_cells+1, cells_a.data(), cells_b.data())){
        cout << setw(wfirst) << "Distribution" << setw(wsecond) << "Gamma 1D" << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Distribution" << setw(wsecond) << "Gamma 1D" << setw(wthird) << "FAIL" << endl;
        pass = false;
    }

    delete[] samples_dream;

    return pass;
}

bool ExternalTester::testGaussian2D(){
    int s = (rngseed == -1) ? 12 : rngseed;
    #ifdef _MSC_VER // using MS Visual Studio
    #if _MSC_VER > 1909 // using MSVC 2017 or newer
    s = 23;
    #endif
    #endif
    TestRNG rng(s);

    int num_cells1d = 4; double delta = 2.0 / ((double) num_cells1d);
    int num_cells = num_cells1d * num_cells1d;
    int num_chains = 100;
    std::vector<double> samples_true(2*num_mc);
    TasDREAM::TruncatedGaussianPDF Gtrue(-0.5, 0.1, -1.0, 1.0);
    Gtrue.overwriteBaseUnifrom(&rng);

    for(int i=0; i<2*num_mc; i++) samples_true[i] = Gtrue.getSample();

    TasDREAM::TasmanianDREAM dream;
    dream.overwriteBaseUnifrom(&rng);
    dream.setProbabilityWeightFunction(&distGauss2D);
    dream.setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.0001);
    gauss.overwriteBaseUnifrom(&rng);
    dream.setCorrectionAll(&gauss);

    std::vector<double> samples_dream;
    dream.collectSamples(3*num_mc, num_mc / num_chains, samples_dream, false);
    //for(int i=0; i<num_mc; i++) samples_dream[i] = -1.0 + 2.0 * u.getSample01();
    //for(int i=0; i<num_mc; i++) cout << samples_dream[i] << endl;

    std::vector<int> cells_a(num_cells, 0);
    std::vector<int> cells_b(num_cells, 0);
    for(int i=0; i<num_mc; i++){
        if ((fabs(samples_true[2*i]) > 1.0) || (fabs(samples_true[2*i+1]) > 1.0)) cout << "ERROR: bad value: " << samples_true[2*i] << "  " << samples_true[2*i + 1] << endl;
        int cx = (int) floor((samples_true[2*i  ] + 1.0) / delta);
        int cy = (int) floor((samples_true[2*i+1] + 1.0) / delta);
        cells_a[cx*num_cells1d + cy]++;

        if ((fabs(samples_dream[2*i]) > 1.0) || (fabs(samples_dream[2*i+1]) > 1.0)) cout << "ERROR: bad value: " << samples_dream[2*i] << "  " << samples_dream[2*i + 1] << endl;
        cx = (int) floor((samples_dream[2*i  ] + 1.0) / delta);
        cy = (int) floor((samples_dream[2*i+1] + 1.0) / delta);
        cells_b[cx*num_cells1d + cy]++;
    }
    //for(int i=0; i<num_cells; i++){ cout << cells_a[i] << "  " << cells_b[i] << endl; }

    bool pass = true;
    if (!testChi(num_cells, cells_a.data(), cells_b.data())){
        cout << "FAIL: mismatch between true Gaussian PDF and dream samples from Gaussian PFD" << endl;
        pass = false;
    }
    //if (pass) cout << "Pass first" << endl;

    TasGrid::TasmanianSparseGrid grid;
    grid.makeGlobalGrid(2, 1, 15, TasGrid::type_iptotal, TasGrid::rule_clenshawcurtis);
    int num_grid_points = grid.getNumNeeded();
    std::vector<double> points;
    grid.getNeededPoints(points);
    std::vector<double> values(num_grid_points);
    for(int i=0; i<num_grid_points; i++){
        values[i] = Gtrue.getDensityLog(points[2*i]) + Gtrue.getDensityLog(points[2*i+1]);
    }
    grid.loadNeededPoints(values);

    TasDREAM::LikelihoodTSG grid_pdf(&grid, true);

    dream.setProbabilityWeightFunction(&grid_pdf);
    dream.setNumChains(num_chains);
    dream.setCorrectionAll(&gauss);
    dream.collectSamples(3*num_mc, num_mc / num_chains, samples_dream, true);

    std::fill(cells_b.begin(), cells_b.end(), 0);
    for(int i=0; i<num_mc; i++){
        if ((fabs(samples_dream[2*i]) > 1.0) || (fabs(samples_dream[2*i+1]) > 1.0)) cout << "ERROR: bad value: " << samples_dream[2*i] << "  " << samples_dream[2*i + 1] << endl;
        int cx = (int) floor((samples_dream[2*i  ] + 1.0) / delta);
        int cy = (int) floor((samples_dream[2*i+1] + 1.0) / delta);
        cells_b[cx*num_cells1d + cy]++;
        //cout << samples_dream[2*i] << "  " << samples_dream[2*i+1] << endl;
    }
    //for(int i=0; i<num_cells; i++){ cout << cells_a[i] << "  " << cells_b[i] << endl; }

    if (!testChi(num_cells, cells_a.data(), cells_b.data())){
        cout << "FAIL: mismatch between true Gaussian PDF and dream samples from interpolated Gaussian PFD" << endl;
        pass = false;
    }

    if (pass){
        cout << setw(wfirst) << "Distribution" << setw(wsecond) << "Truncated Gaussian 2D" << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Distribution" << setw(wsecond) << "Truncated Gaussian 2D" << setw(wthird) << "FAIL" << endl;
    }

    return pass;
}

bool ExternalTester::testModelLikelihoodAlpha(){
    int s = (rngseed == -1) ? 12 : rngseed;
    TestRNG rng(s);

    // model sin(p_0 t M_PI + p_1), data cos(M_PI t) / cos(M_PI t) + cos(5 M_PI y)
    int N = 32; // sample points
    double dt = 1.0 / ((double) N), dt2 = 0.5 * dt;
    int num_chains = 100;

    // set the model
    TasGrid::TasmanianSparseGrid grid;
    grid.makeSequenceGrid(2, N, 36, TasGrid::type_iptotal, TasGrid::rule_leja);
    double domain_a[2] = {0.5, -0.1};
    double domain_b[2] = {8.0,  1.7};
    grid.setDomainTransform(domain_a, domain_b);
    int num_grid_points = grid.getNumNeeded();
    std::vector<double> points;
    grid.getNeededPoints(points);
    std::vector<double> values(num_grid_points * N);
    for(int i=0; i<num_grid_points; i++){
        for(int j=0; j<N; j++){
            values[i*N + j] = sin(points[2*i] * M_PI * (dt2 + j*dt) + points[2*i+1]);
            //values[i*N + j] = sin(points[2*i] * M_PI * (dt2 + j*dt));
        }
    }
    grid.loadNeededPoints(values);

    // set the data
    std::vector<double> data(N);
    for(int j=0; j<N; j++){
        data[j] = sin(M_PI * (dt2 + j*dt) + 0.3 * M_PI);
    }

    TasDREAM::PosteriorFromModel post(&grid);

    double scale = 1.0 / ((double) (N));
    TasDREAM::GaussianLikelihood likely(N, TasDREAM::likely_gauss_scale, &scale, 1, data.data());
    post.setLikelihood(&likely);

    TasDREAM::TasmanianDREAM dream;
    dream.overwriteBaseUnifrom(&rng);
    dream.setProbabilityWeightFunction(&post);

    dream.setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.0004);
    gauss.overwriteBaseUnifrom(&rng);
    dream.setCorrectionAll(&gauss);

    std::vector<double> mcmc;
    dream.collectSamples(num_mc, num_mc / 10, mcmc);
    double frequency = 0.0, correction = 0.0;
    for(int i=0; i<num_chains * num_mc / 10; i++){
        //cout << mcmc[2*i] << "   " << mcmc[2*i+1] << endl;
        frequency += mcmc[2*i];
        correction += mcmc[2*i+1];
    }
    frequency /= (double)(num_chains * num_mc / 10);
    correction /= (double)(num_chains * num_mc / 10);

    if (verbose){
        cout << "inferred frequency: " << frequency << "   actual: 1.0" << endl;
        cout << "inferred correction: " << correction << "   actual: " << 0.3 * M_PI << endl;
    }

    bool pass = true;
    if ((fabs(frequency - 1.0) > 0.015) || (fabs(correction - 0.3 * M_PI) > 0.03)){
        cout << "FAIL: could not identify the frequency and correction for the inference test" << endl;
        pass = false;
    }

    if (pass){
        cout << setw(wfirst) << "Inference" << setw(wsecond) << "Signal Frequency 2D" << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Inference" << setw(wsecond) << "Signal Frequency 2D" << setw(wthird) << "FAIL" << endl;
    }

    return true;
}

bool ExternalTester::testKS(int num_cells, const int count_a[], const int count_b[]){
    int total_a = count_a[0];
    for(int i=1; i<num_cells; i++) total_a += count_a[i];
    int total_b = count_b[0];
    for(int i=1; i<num_cells; i++) total_b += count_b[i];
    if (verbose) cout << "Totals: " << total_a << "  " << total_b << endl;

    int sum_a = 0, sum_b = 0;
    double sum_max = 0.0;

    for(int i=0; i<num_cells; i++){
        sum_a += count_a[i];
        sum_b += count_b[i];

        double diff = fabs(((double) sum_a) / ((double) total_a) - ((double) sum_b) / ((double) total_b));
        if (diff > sum_max) sum_max = diff;
    }

    bool pass = (sum_max < 1.63 * sqrt(((double)(total_a + total_b)) / ((double)(total_a * total_b))));
    if (!pass || verbose){
        cout << "Komogorov-Smirnov test value = " << sum_max << " num cells: " << num_cells << endl;
        cout << "Critical K-S value = " << 1.63 * sqrt(((double)(total_a + total_b)) / ((double)(total_a * total_b))) << endl;
    }
    return pass;
}

bool ExternalTester::testChi(int num_cells, const int count_a[], const int count_b[]){
    double test_value = 0.0;

    int total_a = count_a[0];
    for(int i=1; i<num_cells; i++) total_a += count_a[i];
    int total_b = count_b[0];
    for(int i=1; i<num_cells; i++) total_b += count_b[i];
    if (verbose) cout << "Totals: " << total_a << "  " << total_b << endl;

    double scale_a = sqrt(((double) total_b) / ((double) total_a));
    double scale_b = sqrt(((double) total_a) / ((double) total_b));
    for(int i=0; i<num_cells; i++){
        double diff = scale_a * ((double) count_a[i]) - scale_b * ((double) count_b[i]);
        int sum_count = count_a[i] + count_b[i];
        if (sum_count > 0) test_value += diff * diff / ((double) sum_count);
    }

    // degrees of freedom is one less than num_cells
    bool pass = (test_value < getChiValue(num_cells - 1));

    if (!pass || verbose){
        cout << "Chi-Squared test value = " << test_value << " num cells: " << num_cells << endl;
        cout << "Critical Chi-Squared value = " << getChiValue(num_cells - 1) << endl;
    }

    return pass;
}
double ExternalTester::getChiValue(int num_degrees){
    switch(num_degrees){
        case   9: return  21.666;
        case  15: return  30.578;
        case  19: return  36.191;
        case  20: return  37.566;
        case  24: return  42.980;
        case  49: return  74.919;
        case  99: return 134.642;
        case 124: return 163.546;
        default:
            cerr << "ERROR: unhanded number of degrees of freedom!" << endl;
            return -1.0;
    }
}
*/
#endif
