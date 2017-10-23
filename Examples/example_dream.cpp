#include <iostream>
#include <iomanip>
#include <ctime>
#include "math.h"

//#include "TasmanianSparseGrid.hpp" // TasmanianDREAM.hpp includes TasmanianSparseGrid.hpp
#include "TasmanianDREAM.hpp"

using namespace std;

int main(int argc, const char**){
/*
 * The purpose of this file is to demonstrate the proper way to call
 * functions from the TASMANIAN DREAM Library.
 * The 6 examples are for illustrative purposes only, those are not
 * necessarily the best way to solve the corresponding problem.
 *
 */

    srand(time(0));

    bool limit_examples = false;
    if (argc > 1){
        // some of the exmaples take long time, for post-install testing purposes
        // optionally limit examples to the first 2 (covers dream and sparse grids)
        limit_examples = true;
    }

{
// EXAMPLE 1: make your own probability distribution
//            sample from Gaussian distribution: f(x) = C*exp(-x^2/2) over (-inf, +inf)
    TasGrid::TasmanianSparseGrid grid;

    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 1: sample from Gaussian distribution: f(x) = exp(-x^2/2)" << endl;
    cout << "           ingnoring scaling constants, using 3000 smaples" << endl << endl;

    // decalre a class that describes the PDF
    class UnscaledPDF : public TasDREAM::ProbabilityWeightFunction{
    public:
        UnscaledPDF(){}
        ~UnscaledPDF(){}
        int getNumDimensions() const{ return 1; }
        void evaluate(int num_points, const double x[], double y[], bool useLogForm){
            for(int i=0; i<num_points; i++){ // set the pdf values
                y[i] = -0.5 * x[i] * x[i];
            }
            if (!useLogForm){
                for(int i=0; i<num_points; i++){ // take exponential (if not working with logForm)
                    y[i] = exp(y[i]);
                }
            }
        }
        void getDomainBounds(bool* lower_bound, bool* upper_bound){
            lower_bound[0] = false; // Gaussian is unbounded
            upper_bound[0] = false;
        }
        void getDomainBounds(double* lower_bound, double* upper_bound){
            lower_bound[0] = 0.0; // since bounds above give false,
            upper_bound[0] = 0.0; // those here are dummy values
        }
        void getInitialSample(double x[]){
            // initialize with samples unifromly in [-1,1]
            x[0] = -1.0 + 2.0 * u.getSample01();
        }
    private:
        TasDREAM::CppUniformSampler u;
    };

    UnscaledPDF *exmaple_pdf = new UnscaledPDF();

    int num_chains = 30;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 100;
    // the total number of samples is num_chains * num_iterations
    int num_samples = num_chains * num_sample_iterations;

    TasDREAM::TasmanianDREAM *dream = new TasDREAM::TasmanianDREAM();
    dream->setProbabilityWeightFunction(exmaple_pdf);
    dream->setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.0001);
    dream->setCorrectionAll(&gauss);

    double *samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, false);

    double mean, variance;
    // compute the mean and variance
    mean = 0.0;
    variance = 0.0;
    for(int i=0; i<num_samples; i++){
        mean += samples[i];
        variance += samples[i] * samples[i];
    }
    mean /= (double) num_samples;
    variance /= (double) num_samples;
    variance -= mean * mean;

    cout << "Using regular form:" << endl;
    cout << "       mean:" << setw(13) << mean     << "   error:" << setw(12) << fabs(mean) << endl;
    cout << "   variance:" << setw(13) << variance << "   error:" << setw(12) << fabs(variance - 1.0) << endl;

    delete[] samples;

    // reset the dream samples
    dream->setProbabilityWeightFunction(exmaple_pdf);
    dream->setNumChains(num_chains);
    dream->setCorrectionAll(&gauss);

    samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

    // compute the mean and variance
    mean = 0.0;
    variance = 0.0;
    for(int i=0; i<num_samples; i++){
        mean += samples[i];
        variance += samples[i] * samples[i];
    }
    mean /= (double) num_samples;
    variance /= (double) num_samples;
    variance -= mean * mean;

    cout << "Using log form:" << endl;
    cout << "       mean:" << setw(13) << mean     << "   error:" << setw(12) << fabs(mean) << endl;
    cout << "   variance:" << setw(13) << variance << "   error:" << setw(12) << fabs(variance - 1.0) << endl << endl;

    delete[] samples;
    delete dream;
    delete exmaple_pdf; // must delete the example_pdf, dream will not delete that pointer
}

{
// EXAMPLE 2: use sparse grids to interpolate a likelihood function
//            set inference problem, identify x_0 and x_1 model parameters from data (noise free example)
//            model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(M_PI*t + 0.3*M_PI)
//            t in [0,1], discretized with 32 equidistant nodes
//            likelihood is exp(- 32 * (f(x) - d)^2),
//  NOTE: the likelihood scale must correspond to the noise in the data and d is deteministic
//        however, our focus is the L-2 norm of || f(x) - d ||_2 and discretization introduces error 1/32
    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 2: set inference problem, identify x_0 and x_1 model parameters from data (noise free example)" << endl;
    cout << "           model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(M_PI*t + 0.3*M_PI)" << endl;
    cout << "           t in [0,1], discretized with 32 equidistant nodes" << endl;
    cout << "           likelihood is exp(- 32 * (f(x) - d)^2)" << endl;
    cout << "           use sparse grid to interpolate the likelihood" << endl;
    cout << "     NOTE: 32 corresponds to discretization error in t" << endl << endl;

    int num_chains = 100;
    int num_burnup_iterations = 3000;
    int num_sample_iterations = 100;
    // the total number of samples is num_chains * num_iterations
    int num_samples = num_chains * num_sample_iterations;

    int N = 32;
    double dt = 1.0 / ((double) N), dt2 = 0.5 * dt; // discretization in t

    double *data = new double[N];
    for(int j=0; j<N; j++){
        data[j] = sin(M_PI * (dt2 + j*dt) + 0.3 * M_PI);
    }

    TasGrid::TasmanianSparseGrid grid;
    grid.makeGlobalGrid(2, 1, 15, TasGrid::type_iptotal, TasGrid::rule_clenshawcurtis); // 15-th order polynomial
    double domain_a[2] = {0.5, -0.1}; // set search interval x_0 in [0.5, 8], x_1 in [-0.1, 1.7]
    double domain_b[2] = {8.0,  1.7};
    grid.setDomainTransform(domain_a, domain_b);
    int num_grid_points = grid.getNumNeeded();
    double *points = grid.getNeededPoints(); // sparse grids nodes
    double *values = new double[num_grid_points];
    for(int i=0; i<num_grid_points; i++){
        double exponent = 0.0;
        for(int j=0; j<N; j++){
            double difference = sin(points[2*i] * M_PI * (dt2 + j*dt) + points[2*i+1]) - data[j];
            exponent += difference * difference;
        }
        values[i] = - ((double) N) * exponent; // interpolate the log form of the likelihood
    }
    grid.loadNeededPoints(values);
    delete[] values; // release sparse grids nodes and loaded values
    delete[] points;

    // if not using the log form above, i.e., values[i] = exp(- ((double) N) * exponent)
    // then on the line below will have to say "false"
    // LikelihoodTSG will automatically assume prior Uniform over the sparse grid domain
    TasDREAM::LikelihoodTSG *interpolated_likelihood = new TasDREAM::LikelihoodTSG(&grid, true);

    TasDREAM::TasmanianDREAM *dream = new TasDREAM::TasmanianDREAM();
    dream->setProbabilityWeightFunction(interpolated_likelihood);
    dream->setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.0004);
    dream->setCorrectionAll(&gauss);

    double *samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

    // compute expectation for the frequency x_0 and correction x_1
    double frequency = 0.0, correction = 0.0;
    for(int i=0; i<num_samples; i++){
        frequency += samples[2*i];
        correction += samples[2*i+1];
    }
    frequency /= (double)(num_samples);
    correction /= (double)(num_samples);

    cout << "Inferred values (using 15th order polynomial sparse grid):" << endl;
    cout << "   frequency:" << setw(12) << frequency  << "   error:" << setw(12) << fabs(frequency - 1.0) << endl;
    cout << "  correction:" << setw(12) << correction << "   error:" << setw(12) << fabs(correction - 0.3 * M_PI) << endl << endl;

    grid.makeSequenceGrid(2, 1, 30, TasGrid::type_iptotal, TasGrid::rule_leja); // 30-th order polynomial (using faster Sequence grid)
    grid.setDomainTransform(domain_a, domain_b);
    num_grid_points = grid.getNumNeeded();
    points = grid.getNeededPoints(); // sparse grids nodes
    values = new double[num_grid_points];
    for(int i=0; i<num_grid_points; i++){
        double exponent = 0.0;
        for(int j=0; j<N; j++){
            double difference = sin(points[2*i] * M_PI * (dt2 + j*dt) + points[2*i+1]) - data[j];
            exponent += difference * difference;
        }
        values[i] = - ((double) N) * exponent; // interpolate the log form of the likelihood
    }
    grid.loadNeededPoints(values);
    delete[] values; // release sparse grids nodes and loaded values
    delete[] points;

    // interpolated_likelihood still holds the pointer to grid
    // but we need to reset the DREAM sampler
    dream->setProbabilityWeightFunction(interpolated_likelihood); // dream will automatically assume prior Uniform over the sparse grid domain
    dream->setNumChains(num_chains);
    dream->setCorrectionAll(&gauss);

    delete[] samples;
    samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

    // compute expectation for the frequency x_0 and correction x_1
    frequency = 0.0;
    correction = 0.0;
    for(int i=0; i<num_samples; i++){
        frequency += samples[2*i];
        correction += samples[2*i+1];
    }
    frequency /= (double)(num_samples);
    correction /= (double)(num_samples);

    cout << "Inferred values (using 30th order polynomial sparse grid):" << endl;
    cout << "   frequency:" << setw(12) << frequency  << "   error:" << setw(12) << fabs(frequency - 1.0) << endl;
    cout << "  correction:" << setw(12) << correction << "   error:" << setw(12) << fabs(correction - 0.3 * M_PI) << endl << endl;

    delete[] samples;
    delete dream;
    delete interpolated_likelihood;
    delete[] data;
}

    // some of the exmaples take too long for post-install testing purposes
    // optionally limit examples to the first 2 (covers dream and sparse grids)
    if (limit_examples){
        return 0;
    }

{
// EXAMPLE 3: use sparse grids to interpolate a model and use build-in Gaussian likelilhood, posterior is multimodel due to data
//            set inference problem, identify x_0 and x_1 model parameters from data (noise free example)
//            model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(5*M_PI*t + 0.3*M_PI) + sin(10*M_PI*t + 0.1*M_PI)
//            compared to Example 2, the data is superposition of two signals and the posterior is multi-modal
//            t in [0,1], discretized with 32 equidistant nodes
//            likelihood is exp(- 32 * (f(x) - d)^2),
//  NOTE: the likelihood scale must correspond to the noise in the data and d is deteministic
//        however, our focus is the L-2 norm of || f(x) - d ||_2 and discretization introduces error 1/32
    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 3: set inference problem, identify x_0 and x_1 model parameters from data (noise free example)" << endl;
    cout << "           model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(5*M_PI*t + 0.3*M_PI) + sin(10*M_PI*t + 0.1*M_PI)" << endl;
    cout << "           compared to Example 3, the data is superposition of two signals and the posterior is multi-modal" << endl;
    cout << "           t in [0,1], discretized with 32 equidistant nodes" << endl;
    cout << "           likelihood is exp(- 32 * (f(x) - d)^2)" << endl;
    cout << "           use sparse grid to interpolate the model" << endl;
    cout << "     NOTE: 32 corresponds to discretization error in t" << endl << endl;

    // multi-modal distributions require a lot more chains and samples
    int num_chains = 500;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 300;
    // the total number of samples is num_chains * num_iterations
    int num_samples = num_chains * num_sample_iterations;

    int N = 32;
    double dt = 1.0 / ((double) N), dt2 = 0.5 * dt; // discretization in t

    double *data = new double[N];
    for(int j=0; j<N; j++){
        data[j] = sin(5.0 * M_PI * (dt2 + j*dt) + 0.3 * M_PI) + sin(10.0 * M_PI * (dt2 + j*dt) + 0.1 * M_PI);
    }

    TasGrid::TasmanianSparseGrid grid;
    grid.makeSequenceGrid(2, N, 30, TasGrid::type_iptotal, TasGrid::rule_leja); // 30-th order polynomial
    double domain_a[2] = { 1.0, -0.1}; // set search interval x_0 in [1, 12], x_1 in [-0.1, 1.7]
    double domain_b[2] = {12.0,  1.7};
    grid.setDomainTransform(domain_a, domain_b);
    int num_grid_points = grid.getNumNeeded();
    double *points = grid.getNeededPoints(); // sparse grids nodes
    double *values = new double[num_grid_points * N];
    for(int i=0; i<num_grid_points; i++){
        for(int j=0; j<N; j++){
            values[i*N + j] = sin(points[2*i] * M_PI * (dt2 + j*dt) + points[2*i+1]);
        }
    }
    grid.loadNeededPoints(values);
    delete[] values; // release sparse grids nodes and loaded values
    delete[] points;

    // when working with grids with large N, acceleration becomes very useful
    // if BLAS is enabled on compile time, the grid will use BLAS by default
    // GPU acceleration can be enabled here using
    // if (grid.isCudaEnabled()){
    //     grid.enableAcceleration(TasGrid::accel_gpu_cublas);
    //     grid.setGPUID(0);
    // }

    // PosteriorFromModel will automatically assume prior Uniform over the sparse grid domain
    TasDREAM::PosteriorFromModel *post = new TasDREAM::PosteriorFromModel(&grid);

    // set Gaussian likelihood with diagonal covariance and constant diagonal
    double scale = 1.0 / ((double) (N));
    TasDREAM::GaussianLikelihood *likely = new TasDREAM::GaussianLikelihood(N, TasDREAM::likely_gauss_scale, &scale, 1, data);
    post->setLikelihood(likely);
    //poset->setData(2, data); // the GaussianLikelihood stores the data in the constructor, there is no need to use data again here

    TasDREAM::TasmanianDREAM *dream = new TasDREAM::TasmanianDREAM();
    dream->setProbabilityWeightFunction(post);
    dream->setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.0001);
    dream->setCorrectionAll(&gauss);

    double *samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

    // splitting the bimodal data is tricky, use the mid-point of the domain
    // frequencies below 6.5 will be added to low frequency signal
    // frequencies above 6.5 will be added to high frequency signal
    double frequency_low = 0.0, correction_low = 0.0;
    double frequency_high = 0.0, correction_high = 0.0;
    int num_low = 0, num_high = 0;
    for(int i=0; i<num_samples; i++){
        if (samples[2*i] < 6.5){
            frequency_low += samples[2*i];
            correction_low += samples[2*i+1];
            num_low++;
        }else{
            frequency_high += samples[2*i];
            correction_high += samples[2*i+1];
            num_high++;
        }
    }
    //cout << num_low << "  " << num_high << endl;
    frequency_low /= (double)(num_low);
    correction_low /= (double)(num_low);
    frequency_high /= (double)(num_high);
    correction_high /= (double)(num_high);

    cout << "Inferred values:" << endl;
    cout << " low   frequency:" << setw(12) << frequency_low   << "   error:" << setw(12) << fabs(frequency_low - 5.0) << endl;
    cout << " low  correction:" << setw(12) << correction_low  << "   error:" << setw(12) << fabs(correction_low - 0.3 * M_PI) << endl << endl;
    cout << " high  frequency:" << setw(12) << frequency_high  << "   error:" << setw(12) << fabs(frequency_high - 10.0) << endl;
    cout << " high correction:" << setw(12) << correction_high << "   error:" << setw(12) << fabs(correction_high - 0.1 * M_PI) << endl << endl;

    delete[] samples;
    delete dream;
    delete post;
    delete likely;
    delete[] data;
}

{
// EXAMPLE 4: same as Example 3, but the multi-modal posterior comes from symmetry in the model
//            set inference problem, identify x_0, x_1, x_2 and x_3 (four) model parameters from noisy data
//            model: f(x) = x_0*exp(x_1*t) + x_2*exp(x_3*t),  data: d = exp(t) + 0.4*exp(3*t)
//            compared to Example 2, the data is superposition of two signals and the posterior is multi-modal
//            t in [0,1], discretized with 32 equidistant nodes
//            likelihood is exp(- 64 * (f(x) - d)^2),
//  NOTE: the likelihood scale must correspond to the noise in the data and d is deteministic
//        however, our focus is the L-2 norm of || f(x) - d ||_2 and discretization introduces error 1/100
    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 4: set inference problem, identify x_0, x_1, x_2 and x_3 (four) model parameters from noisy data" << endl;
    cout << "           model: f(x) = x_0*exp(x_1*t) + x_2*exp(x_3*t),  data: d = exp(t) + 0.4*exp(3*t)" << endl;
    cout << "           compared to Example 3, the model is symmetric in x_1 and x_3 which gives multi-modal posterior" << endl;
    cout << "           t in [0,1], discretized with 64 equidistant nodes" << endl;
    cout << "           likelihood is exp(-64 * (f(x) - d)^2)" << endl;
    cout << "           use sparse grid to interpolate the model" << endl;
    cout << "     NOTE: 64 corresponds to discretization error in t" << endl << endl;

    // multi-modal distributions require a lot more chains and samples
    int num_chains = 1000;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 100;
    // the total number of samples is num_chains * num_iterations
    int num_samples = num_chains * num_sample_iterations;

    int N = 64;
    double dt = 1.0 / ((double) N), dt2 = 0.5 * dt; // discretization in t

    double *data = new double[N];
    for(int j=0; j<N; j++){
        data[j] = exp(dt2 + j*dt) + 0.4 * exp(3.0 * (dt2 + j*dt));
    }

    TasGrid::TasmanianSparseGrid grid;
    grid.makeSequenceGrid(4, N, 15, TasGrid::type_iptotal, TasGrid::rule_leja); // 15-th order polynomial
    double domain_a[4] = {0.2, 0.5, 0.2, 0.5}; // set search interval x_0 in [1, 12], x_1 in [-0.1, 1.7]
    double domain_b[4] = {1.2, 4.0, 1.2, 4.0};
    grid.setDomainTransform(domain_a, domain_b);
    int num_grid_points = grid.getNumNeeded();
    double *points = grid.getNeededPoints(); // sparse grids nodes
    double *values = new double[num_grid_points * N];
    for(int i=0; i<num_grid_points; i++){
        for(int j=0; j<N; j++){
            values[i*N + j] = points[4*i] * exp(points[4*i+1] * (dt2 + j*dt)) + points[4*i+2] * exp(points[4*i+3] * (dt2 + j*dt));
        }
    }
    grid.loadNeededPoints(values);
    delete[] values; // release sparse grids nodes and loaded values
    delete[] points;

    // when working with grids with large N, acceleration becomes very useful
    // if BLAS is enabled on compile time, the grid will use BLAS
    // GPU acceleration can be enabled here using
    // if (grid.isCudaEnabled()){
    //     grid.enableAcceleration(TasGrid::accel_gpu_fullmemory);
    //     grid.setGPUID(1);
    // }

    // PosteriorFromModel will automatically assume prior Uniform over the sparse grid domain
    TasDREAM::PosteriorFromModel *post = new TasDREAM::PosteriorFromModel(&grid);

    // set Gaussian likelihood with diagonal covariance and constant diagonal
    double scale = 1.0 / ((double) (N));
    TasDREAM::GaussianLikelihood *likely = new TasDREAM::GaussianLikelihood(N, TasDREAM::likely_gauss_scale, &scale, 1, data);
    post->setLikelihood(likely);
    //poset->setData(1, data); // the GaussianLikelihood stores the data in the constructor, there is no need to use data again here

    TasDREAM::TasmanianDREAM *dream = new TasDREAM::TasmanianDREAM();
    dream->setProbabilityWeightFunction(post);
    dream->setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.0001);
    dream->setCorrectionAll(&gauss);

    double *samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

    // splitting the bimodal data is tricky,
    // add constraint that orders the rates
    double rate_low = 0.0,  scale_low = 0.0;
    double rate_high = 0.0, scale_high = 0.0;
    for(int i=0; i<num_samples; i++){
        if (samples[4*i+1] < samples[4*i+3]){
            scale_low  += samples[4*i];
            rate_low   += samples[4*i+1];
            scale_high += samples[4*i+2];
            rate_high  += samples[4*i+3];
        }else{
            scale_low  += samples[4*i+2];
            rate_low   += samples[4*i+3];
            scale_high += samples[4*i];
            rate_high  += samples[4*i+1];
        }
    }
    rate_low   /= (double)(num_samples);
    scale_low  /= (double)(num_samples);
    rate_high  /= (double)(num_samples);
    scale_high /= (double)(num_samples);

    cout << "Inferred values (noise free case):" << endl;
    cout << " low   rate:" << setw(12) << rate_low   << "   error:" << setw(12) << fabs(rate_low - 1.0)   << endl;
    cout << " low  scale:" << setw(12) << scale_low  << "   error:" << setw(12) << fabs(scale_low - 1.0)  << endl << endl;
    cout << " high  rate:" << setw(12) << rate_high  << "   error:" << setw(12) << fabs(rate_high - 3.0)  << endl;
    cout << " high scale:" << setw(12) << scale_high << "   error:" << setw(12) << fabs(scale_high - 0.4) << endl << endl;

    // repeat the example adding noise to the data
    TasDREAM::GaussianPDF noise(0.0, 1.0 / ((double) N));
    for(int j=0; j<N; j++){
        data[j] = exp(dt2 + j*dt) + 0.4 * exp(3.0 * (dt2 + j*dt)) + noise.getSample();
    }

    // reset the likelihood and update with the new data
    delete likely;
    scale = 1.1 / ((double) (N)); // increase the covariance to accout for the extra noise
    likely = new TasDREAM::GaussianLikelihood(N, TasDREAM::likely_gauss_scale, &scale, 1, data);
    post->setLikelihood(likely);

    dream->setProbabilityWeightFunction(post);
    dream->setNumChains(num_chains);
    dream->setCorrectionAll(&gauss);

    delete[] samples;
    samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

    // splitting the bimodal data is tricky,
    // add constraint that orders the rates
    rate_low = 0.0,  scale_low = 0.0;
    rate_high = 0.0, scale_high = 0.0;
    for(int i=0; i<num_samples; i++){
        if (samples[4*i+1] < samples[4*i+3]){
            scale_low  += samples[4*i];
            rate_low   += samples[4*i+1];
            scale_high += samples[4*i+2];
            rate_high  += samples[4*i+3];
        }else{
            scale_low  += samples[4*i+2];
            rate_low   += samples[4*i+3];
            scale_high += samples[4*i];
            rate_high  += samples[4*i+1];
        }
    }
    rate_low   /= (double)(num_samples);
    scale_low  /= (double)(num_samples);
    rate_high  /= (double)(num_samples);
    scale_high /= (double)(num_samples);

    cout << "Inferred values (noisy case):" << endl;
    cout << " low   rate:" << setw(12) << rate_low   << "   error: " << setw(12) << fabs(rate_low - 1.0)   << endl;
    cout << " low  scale:" << setw(12) << scale_low  << "   error: " << setw(12) << fabs(scale_low - 1.0)  << endl << endl;
    cout << " high  rate:" << setw(12) << rate_high  << "   error: " << setw(12) << fabs(rate_high - 3.0)  << endl;
    cout << " high scale:" << setw(12) << scale_high << "   error: " << setw(12) << fabs(scale_high - 0.4) << endl << endl;

    delete[] samples;
    delete dream;
    delete post;
    delete likely;
    delete[] data;
}

{
// EXAMPLE 5: use inference for a custom model, not sparse grid
//            set inference problem, identify x_0 and x_1 model parameters from data (noise free example)
//            model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(5*M_PI*t + 0.3*M_PI) + sin(10*M_PI*t + 0.1*M_PI)
//            compared to Example 2, the data is superposition of two signals and the posterior is multi-modal
//            t in [0,1], discretized with 32 equidistant nodes
//            likelihood is exp(- 32 * (f(x) - d)^2),
//  NOTE: the likelihood scale must correspond to the noise in the data and d is deteministic
//        however, our focus is the L-2 norm of || f(x) - d ||_2 and discretization introduces error 1/32
    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 5: set inference problem, identify x_0 and x_1 model parameters from data (noise free example)" << endl;
    cout << "           model: f(x) = sin(x_0*M_PI*t + x_1),  data: d = sin(5*M_PI*t + 0.3*M_PI) + sin(10*M_PI*t + 0.1*M_PI)" << endl;
    cout << "           compared to Example 2, the data is superposition of two signals and the posterior is multi-modal" << endl;
    cout << "           t in [0,1], discretized with 32 equidistant nodes" << endl;
    cout << "           likelihood is exp(- 32 * (f(x) - d)^2)" << endl;
    cout << "           compared to Example 3, we use a custom model class as opposed to a sparse grid" << endl;
    cout << "     NOTE: 32 corresponds to discretization error in t" << endl << endl;

    // multi-modal distributions require a lot more chains and samples
    int num_chains = 500;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 300;
    // the total number of samples is num_chains * num_iterations
    int num_samples = num_chains * num_sample_iterations;

    int N = 32;
    double dt = 1.0 / ((double) N), dt2 = 0.5 * dt; // discretization in t

    double *data = new double[N];
    for(int j=0; j<N; j++){
        data[j] = sin(5.0 * M_PI * (dt2 + j*dt) + 0.3 * M_PI) + sin(10.0 * M_PI * (dt2 + j*dt) + 0.1 * M_PI);
    }

    class Exmaple5SignalModel : public TasDREAM::CustomModelWrapper{
    public:
        Exmaple5SignalModel(int num_discrete_t) : N(num_discrete_t){
            dt = 1.0 / ((double) N);
            dt2 = 0.5 * dt; // discretization in t
        }
        ~Exmaple5SignalModel(){}

        int getNumDimensions() const{ return 2; }
        int getNumOutputs() const{ return N; }
        void evaluate(const double x[], int num_points, double y[]) const{
            for(int i=0; i<num_points; i++){
                for(int j=0; j<N; j++){
                    y[i*N + j] = sin(x[2*i] * M_PI * (dt2 + j*dt) + x[2*i+1]);
                }
            }
        }
    private:
        int N;
        double dt, dt2;
    };

    Exmaple5SignalModel *model = new Exmaple5SignalModel(N);

    // PosteriorFromModel cannot guess priors and ranges for the model, must specify manually
    TasDREAM::PosteriorFromModel *post = new TasDREAM::PosteriorFromModel(model);

    TasDREAM::BetaPDF *prior_frequency = new TasDREAM::BetaPDF(1.0, 12.0, 4.0, 2.0);
    post->overwritePDF(0, prior_frequency);
    TasDREAM::UniformPDF *prior_correction = new TasDREAM::UniformPDF(-0.1, 1.7);
    post->overwritePDF(1, prior_correction);
    // make sure to delete the pointers in the end

    // set Gaussian likelihood with diagonal covariance and constant diagonal
    double scale = 1.0 / ((double) (N));
    TasDREAM::GaussianLikelihood *likely = new TasDREAM::GaussianLikelihood(N, TasDREAM::likely_gauss_scale, &scale, 1, data);
    post->setLikelihood(likely);
    //poset->setData(1, data); // the GaussianLikelihood stores the data in the constructor, there is no need to use data again here

    TasDREAM::TasmanianDREAM *dream = new TasDREAM::TasmanianDREAM();
    dream->setProbabilityWeightFunction(post);
    dream->setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.0001);
    dream->setCorrectionAll(&gauss);

    double *samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

    // splitting the bimodal data is tricky, use the mid-point of the domain
    // frequencies below 6.5 will be added to low frequency signal
    // frequencies above 6.5 will be added to high frequency signal
    double frequency_low = 0.0, correction_low = 0.0;
    double frequency_high = 0.0, correction_high = 0.0;
    int num_low = 0, num_high = 0;
    for(int i=0; i<num_samples; i++){
        if (samples[2*i] < 6.5){
            frequency_low += samples[2*i];
            correction_low += samples[2*i+1];
            num_low++;
        }else{
            frequency_high += samples[2*i];
            correction_high += samples[2*i+1];
            num_high++;
        }
    }
    //cout << num_low << "  " << num_high << endl;
    frequency_low /= (double)(num_low);
    correction_low /= (double)(num_low);
    frequency_high /= (double)(num_high);
    correction_high /= (double)(num_high);

    cout << "Inferred values:" << endl;
    cout << " low   frequency:" << setw(12) << frequency_low   << "   error:" << setw(12) << fabs(frequency_low - 5.0) << endl;
    cout << " low  correction:" << setw(12) << correction_low  << "   error:" << setw(12) << fabs(correction_low - 0.3 * M_PI) << endl << endl;
    cout << " high  frequency:" << setw(12) << frequency_high  << "   error:" << setw(12) << fabs(frequency_high - 10.0) << endl;
    cout << " high correction:" << setw(12) << correction_high << "   error:" << setw(12) << fabs(correction_high - 0.1 * M_PI) << endl << endl;

    delete prior_correction;
    delete prior_frequency;
    delete model;
    delete[] samples;
    delete dream;
    delete post;
    delete likely;
    delete[] data;
}

{
// EXAMPLE 6: use inference for a custom model (not sparse grid) and a custom defined likelihood (not Gaussian)
//            shows how to set custom likelihood function and combine that with a model
//            custom model is used here, but sparse grids model can be used too
//            essentially shows how TasmanianDREAM can give approximate solution to an optimization problem
//            assuming here that the optimum corresponds to the pdf mode in a high-probability region
//
    cout << endl << "-------------------------------------------------------------------------------------------------" << endl;
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 6: infer six coefficients of a model f(x) = x_0 + sum_k x_k sin(k*M_PI*t) for k = 1, ..., 5" << endl;
    cout << "           data: d = 2.0 * sin(2*M_PI*t) + sin(4*M_PI*t) + noise" << endl;
    cout << "           t in [0,1], discretized N randomly chosen locations" << endl;
    cout << "           use 2 cases N=36 samples with noise and N=12 samples with no noise" << endl;
    cout << "           testing 2 likelihood choises, gaussian: exp(- N * (f(x) - d)^2)" << endl;
    cout << "           and custom implemented l-1: exp(- N * |f(x) - d|)" << endl;
    cout << "           compared to other exmaples, we are looking for optimum, i.e., best x" << endl;
    cout <<"            not statistical mean of the samples" << endl;
    cout <<"      NOTE: the output here can be used as initial step to a more sophisticated optimization loop" << endl << endl;

    // multi-modal distributions require a lot more chains and samples
    int num_chains = 1000;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 300;
    // the total number of samples is num_chains * num_iterations
    int num_samples = num_chains * num_sample_iterations;

    int N = 40;
    TasDREAM::UniformPDF unif(0.0, 1.0);
    TasDREAM::UniformPDF noise(-0.2, 0.2); // this is strong noise

    double *data = new double[N];
    double *snap = new double[N];
    for(int j=0; j<N; j++){
        double t = unif.getSample();
        snap[j] = t;
        data[j] = 2.0 * sin(2.0 * M_PI * t) + sin(4.0 * M_PI * t) + noise.getSample();
    }

    class Exmaple6ApproximationModel : public TasDREAM::CustomModelWrapper{
    public:
        Exmaple6ApproximationModel(int num_discrete_t, const double measurement_snap[]) :
            N(num_discrete_t), snap(measurement_snap){
        }
        ~Exmaple6ApproximationModel(){}

        int getNumDimensions() const{ return 6; }
        int getNumOutputs() const{ return N; }
        void evaluate(const double x[], int num_points, double y[]) const{
            for(int i=0; i<num_points; i++){
                for(int j=0; j<N; j++){
                    double t = snap[j];
                    y[i*N + j] = x[6*i];
                    for(int k=1; k<6; k++){
                        y[i*N + j] += x[6*i + k] * sin(((double) k) * M_PI * t);
                    }
                }
            }
        }
    private:
        int N;
        const double *snap;
    };

    Exmaple6ApproximationModel *custom_model = new Exmaple6ApproximationModel(N, snap);

    // PosteriorFromModel cannot guess priors and ranges for the model, must specify manually
    TasDREAM::PosteriorFromModel *post = new TasDREAM::PosteriorFromModel(custom_model);

    TasDREAM::GaussianPDF *prior_wide_gaussian = new TasDREAM::GaussianPDF(0.0, 10.0);
    for(int k=0; k<6; k++) post->overwritePDF(k, prior_wide_gaussian);
    // make sure to delete the pointers in the end

    // set Gaussian likelihood with diagonal covariance and constant diagonal, basically L-2 projection
    double scale = 1.0 / ((double) (N));
    TasDREAM::GaussianLikelihood *likely = new TasDREAM::GaussianLikelihood(N, TasDREAM::likely_gauss_scale, &scale, 1, data);
    post->setLikelihood(likely);
    //poset->setData(1, data); // the GaussianLikelihood stores the data in the constructor, there is no need to use data again here

    TasDREAM::TasmanianDREAM *dream = new TasDREAM::TasmanianDREAM();
    dream->setProbabilityWeightFunction(post);
    dream->setNumChains(num_chains);
    TasDREAM::GaussianPDF gauss(0.0, 0.0001);
    dream->setCorrectionAll(&gauss);

    double *samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

	// using the sample mean is not proper when working with optimiation problem
	// since the samples were clustered around the high probability regions
	// assuming high-probability equates a mode (extrema), look for the sample
	// with largest probability
	double true_coefficients[6] = {0.0, 0.0, 2.0, 0.0, 1.0, 0.0};
    double *values = dream->getPDFHistory();
    int best_index = 0;
    double best_value = values[0];
    for(int i=1; i<num_samples; i++){
		if (values[i] > best_value){
			best_value = values[i];
			best_index = i;
		}
	}
	double best_sample[6];
	for(int k=0; k<6; k++) best_sample[k] = samples[6*best_index + k];

	cout << "Large sample size results, N = 36" << endl;
	cout << "Using least-squares mean fit:" << endl;
	cout << "coefficient: ";
    for(int k=0; k<6; k++) cout << setw(13) << best_sample[k];
    cout << endl;
    cout << "      error: ";
    for(int k=0; k<6; k++) cout << setw(13) << fabs(best_sample[k] - true_coefficients[k]);
    cout << endl << endl;

    // define custom likelihood and use L1 minimization
    class L1Likelihood : public TasDREAM::BaseLikelihood{
    public:
        L1Likelihood(int num_outputs, double difference_scale) : N(num_outputs), scale(difference_scale){}
        ~L1Likelihood(){}

        double* getLikelihood(int num_model, const double *model, int num_data, const double *ddata, double *likelihood = 0, bool useLogForm = true){
			double *result = (likelihood != 0) ? likelihood : (new double[num_model]);
			for(int i=0; i<num_model; i++){
				result[i] = 0.0;
				for(int j=0; j<num_data; j++){
					for(int k=0; k<N; k++){
						result[i] -= fabs(model[i*N + k] - ddata[j*N + k]);
					}
				}
			}
			if (!useLogForm){
				for(int i=0; i<num_model; i++) result[i] = exp(scale * result[i]);
			}else{
				for(int i=0; i<num_model; i++) result[i] *= scale;
			}
			return (likelihood != 0) ? 0 : result;
		}
    private:
		int N;
		double scale;
    };

    L1Likelihood *likely1 = new L1Likelihood(N, ((double) N));
    // the model is already set in the posterior, set the new likelihood and data
    post->setLikelihood(likely1);
    // the l1 problem does not rely on precomputing values for speedup
    // we will feed raw data to each call to L1Likelihood::getLikelihood()
    post->setData(1, data);

    dream->setProbabilityWeightFunction(post);
    dream->setNumChains(num_chains);
    dream->setCorrectionAll(&gauss);

	delete[] samples;
    samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

    delete[] values;
    values = dream->getPDFHistory();
    best_index = 0;
    best_value = values[0];
    for(int i=1; i<num_samples; i++){
		if (values[i] > best_value){
			best_value = values[i];
			best_index = i;
		}
	}
	for(int k=0; k<6; k++) best_sample[k] = samples[6*best_index + k];

	cout << "Using l1 minimization fit:" << endl;
	cout << "coefficient: ";
    for(int k=0; k<6; k++) cout << setw(13) << best_sample[k];
    cout << endl;
    cout << "      error: ";
    for(int k=0; k<6; k++) cout << setw(13) << fabs(best_sample[k] - true_coefficients[k]);
    cout << endl << endl;

    //////////////  using small data fit and reusing already defined variables //////////////
    N = 12;

	delete[] data; delete[] snap;
    data = new double[N];
    snap = new double[N];
    for(int j=0; j<N; j++){
        double t = unif.getSample();
        snap[j] = t;
        data[j] = 2.0 * sin(2.0 * M_PI * t) + sin(4.0 * M_PI * t);
    }

	// reset the model
    delete custom_model;
    custom_model = new Exmaple6ApproximationModel(N, snap);

    // PosteriorFromModel cannot guess priors and ranges for the model, must specify manually
    delete post;
    post = new TasDREAM::PosteriorFromModel(custom_model);

    for(int k=0; k<6; k++) post->overwritePDF(k, prior_wide_gaussian);

    // set Gaussian likelihood with diagonal covariance and constant diagonal, basically L-2 projection
    scale = 1.0 / ((double) (2*N));
    delete likely;
    likely = new TasDREAM::GaussianLikelihood(N, TasDREAM::likely_gauss_scale, &scale, 1, data);
    post->setLikelihood(likely);
    //poset->setData(1, data); // the GaussianLikelihood stores the data in the constructor, there is no need to use data again here

    dream->setProbabilityWeightFunction(post); // resets the sampler
    dream->setNumChains(num_chains);
    dream->setCorrectionAll(&gauss);

    delete[] samples;
    samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

	delete[] values;
    values = dream->getPDFHistory();
    best_index = 0;
    best_value = values[0];
    for(int i=1; i<num_samples; i++){
		if (values[i] > best_value){
			best_value = values[i];
			best_index = i;
		}
	}
	for(int k=0; k<6; k++) best_sample[k] = samples[6*best_index + k];

	cout << "Small sample size results, N = 6" << endl;
	cout << "Using least-squares mean fit:" << endl;
	cout << "coefficient: ";
    for(int k=0; k<6; k++) cout << setw(13) << best_sample[k];
    cout << endl;
    cout << "      error: ";
    for(int k=0; k<6; k++) cout << setw(13) << fabs(best_sample[k] - true_coefficients[k]);
    cout << endl << endl;


    delete likely1;
    likely1 = new L1Likelihood(N, ((double)(2*N)));
    // the model is already set in the posterior, set the new likelihood and data
    post->setLikelihood(likely1);
    // the l1 problem does not rely on precomputing values for speedup
    // we will feed raw data to each call to L1Likelihood::getLikelihood()
    post->setData(1, data);

    dream->setProbabilityWeightFunction(post);
    dream->setNumChains(num_chains);
    dream->setCorrectionAll(&gauss);

	delete[] samples;
    samples = dream->collectSamples(num_burnup_iterations, num_sample_iterations, true);

    // could adjust the jump scale here, but did not see much difference
    //delete[] samples;
    //dream->setJumpScale(0.5);
    //samples = dream->collectSamples(0, num_sample_iterations, true);

    delete[] values;
    values = dream->getPDFHistory();
    best_index = 0;
    best_value = values[0];
    for(int i=1; i<num_samples; i++){
		if (values[i] > best_value){
			best_value = values[i];
			best_index = i;
		}
	}
	for(int k=0; k<6; k++) best_sample[k] = samples[6*best_index + k];

	cout << "Using l1 minimization fit:" << endl;
	cout << "coefficient: ";
    for(int k=0; k<6; k++) cout << setw(13) << best_sample[k];
    cout << endl;
    cout << "      error: ";
    for(int k=0; k<6; k++) cout << setw(13) << fabs(best_sample[k] - true_coefficients[k]);
    cout << endl << endl;

    cout << "  Note: the l1 method is expected to be a bit better when the sample size is small" << endl;
    cout << "        and the solution x is sparse, i.e., most of the components are zeros." << endl;
    cout << "        This example demonstrates how to set a custom likelihood, not make l1 vs. l2 comparison." << endl << endl;


    delete prior_wide_gaussian;
    delete custom_model;
    delete[] samples;
    delete[] values;
    delete dream;
    delete post;
    delete likely;
    delete likely1;
    delete[] data;
    delete[] snap;
}

    cout << endl << "-------------------------------------------------------------------------------------------------" << endl << endl;

    return 0;
}
