#include "Tasmanian.hpp"

using namespace std;

/*!
 * \internal
 * \file example_dream_05.cpp
 * \brief Examples for the Tasmanian DREAM module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAMExamples
 *
 * Tasmanian DREAM Example 5
 * \endinternal
 */

/*!
 * \ingroup TasmanianDREAMExamples
 * \addtogroup TasmanianDREAMExamples5 Tasmanian DREAM module, example 5
 *
 * Example 5: Given noisy data that is the superposition of two sin-waves,
 * find the frequency and magnitude of the waves, the frequency is
 * on of 5 possible integer values.
 * Unlike the previous cases where we were looking at the mean and variance
 * of the solution, here we consider the mode of the posterior distribution,
 * i.e., we want the values that give us best-fit in deterministic sense.
 * The problem can be viewed as a deterministic optimization problem,
 * but the noise in the data would offer significant challenges
 * for any global optimization scheme; even if the solution computed
 * by the DREAM algorithm is not sufficiently "close", it is still
 * a very good initial guess for local optimization methods.
 *
 * This example shows how to use DREAM with a custom model
 * (without Tasmanian sparse grids) and a custom likelihood
 * (not implemented in Tasmanian).
 * The example also shows how to use DREAM sampling to search
 * for the (approximate) solution to an optimization problem.
 */

//! \brief DREAM Example 5: signal decomposition, finding the best fit
//! \ingroup TasmanianDREAMExamples5

//! \snippet DREAM/Examples/example_dream_05.cpp DREAM_Example_05 example
void dream_example_05(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_05 example]
#endif
    // using the default random engine, but must reset the random number generator
    srand((int) time(nullptr));

    // EXAMPLE 5:
    cout << "\n" << "---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(5);
    cout << "EXAMPLE 5: infer the frequency and magnitude of two signals from noisy data\n"
         << "           the model has 5 parameters: f(x_1 ... x_5) = sum x_k sin(k * pi * t)\n"
         << "           data = 2.0 * sin(2 * pi * t) + sin(4 * pi * t) + noise\n"
         << "           t in [0, 1], t is discretized using 64 equidistant nodes\n"
         << "           we use two different likelihood functions,"
         <<            "corresponding to l-2 and l-1 norms\n"
         << "           we are looking for the mode of the posterior,"
         <<            "i.e., the optimal fit to the data\n\n";

    constexpr double pi = 3.14159265358979323846; // half-period of the std::sin() function

    // higher dimensions require more samples
    int num_dimensions = 5;
    int num_chains = 50;
    int num_burnup_iterations = 1000;
    int num_sample_iterations = 1000;
    // the total number of samples is num_chains * num_iterations
    int num_discrete_nodes = 64;

    // create a lambda function that represents the model
    // normally this would be a call to an external code
    auto model = [&](std::vector<double> const &x, std::vector<double> &y)->
        void{
            double dt = 1.0 / ((double) y.size());
            double t = 0.5 * dt;
            for(auto &output : y){
                output = 0.0;
                int frequency = 1;
                for(auto const &weight : x)
                    output += weight * std::sin(double(frequency++) * t * pi);
                t += dt;
            }
        };

    // the model defined above deals with single set of inputs
    // sampling requires that models are computed in batch
    // the size of x will always be a multiple of the number of inputs
    // the size of y MUST be set to the corresponding multiple of the number of outputs
    // batch evaluations are best done in parallel
    auto batch_model = [&](std::vector<double> const &x, std::vector<double> &y)->
        void{
            int num_samples = (int) x.size() / num_dimensions;
            y.resize(num_samples * num_discrete_nodes);

            // use up to 4 threads (if available)
            int num_threads = std::min((int) std::thread::hardware_concurrency(), 4);

            std::vector<std::thread> workers(num_threads);
            for(int start = 0; start<num_threads; start++){
                workers[start] = std::thread([&, start]()->
                    void{
                        for(int i=start; i<num_samples; i+=num_threads){
                            std::vector<double> single_input(&x[i*num_dimensions],
                                                            &x[i*num_dimensions] + num_dimensions);
                            std::vector<double> single_output(num_discrete_nodes);

                            model(single_input, single_output);

                            std::copy(single_output.begin(), single_output.end(),
                                      &y[i*num_discrete_nodes]);
                        }
                    });
            }

            for(auto &w : workers) w.join();
        };

    std::vector<double> signal = {0.0, 2.0, 0.0, 1.0, 0.0};
    std::vector<double> data(num_discrete_nodes);
    model(signal, data);

    // add noise to the data, use magnitude 1 / num_discrete_nodes
    // you can adjust the example to consider more/less noise
    TasDREAM::applyUniformUpdate(data, 1.0 / ((double) num_discrete_nodes));

    // first use Gaussian likelihood: exp( - sigma * (f(x) - d)^2 )
    // note that the numerator uses the l-2 norm of the difference between model and data
    // using smaller variance since we are looking for best-fit (even with the noise)
    TasDREAM::LikelihoodGaussIsotropic likely(1.0 / ((double) num_discrete_nodes/2), data);

    // Define the search domain, each parameter is assumed to be in [0, 3]
    std::vector<double> lower(num_dimensions, 0.0);
    std::vector<double> upper(num_dimensions, 3.0);

    TasDREAM::TasmanianDREAM state(num_chains, num_dimensions);
    auto initial_chains = TasDREAM::genUniformSamples(lower, upper, num_chains); // uniform initial state
    state.setState(initial_chains);

    constexpr auto sampling_form = TasDREAM::logform; // ensure uniform sampling form
    TasDREAM::SampleDREAM<sampling_form>
                         (num_burnup_iterations, num_sample_iterations,
                          TasDREAM::posterior<sampling_form>
                                (batch_model, // must use the batch model
                                 likely,      // provide the likelihood
                                 TasDREAM::uniform_prior), // assume non-informative prior
                          TasDREAM::hypercube(lower, upper),
                          state,
                          TasDREAM::no_update,
                          TasDREAM::const_percent<100> // use only the differential update
                          );

    std::vector<double> solution = state.getApproximateMode();

    //cout << " l-2 acceptance rate: " << state.getAcceptanceRate() << endl;
    cout << "Using Gaussian likelihood, the computed solution is:\n"
         << " computed: " << std::fixed;
    for(auto x : solution) cout << setw(13) << x;
    cout << "\n    error: " << std::scientific;
    for(int i=0; i<num_dimensions; i++) cout << setw(13) << std::abs(solution[i] - signal[i]);
    cout << "\n\n";


    // Change the likelihood to use the l-1 norm
    // Combine the model and the likelihood in a single function similar to batch_model()
    // The main difference is that y doesn't need to be resized
    // and the likelihood is applied right after the model
    // Note that the sampling form has to be hard-coded or captured
    auto model_likelihood = [&](std::vector<double> const &x,
                                std::vector<double> &y)->
        void{
            int num_samples = (int) x.size() / num_dimensions;

            // use up to 4 threads (if available)
            int num_threads = std::min((int) std::thread::hardware_concurrency(), 4);

            std::vector<std::thread> workers(num_threads);
            for(int start = 0; start<num_threads; start++){
                workers[start] = std::thread([&, start]()->
                    void{
                        for(int i=start; i<num_samples; i+=num_threads){
                            std::vector<double> single_input(&x[i*num_dimensions],
                                                            &x[i*num_dimensions] + num_dimensions);
                            std::vector<double> single_output(num_discrete_nodes);

                            model(single_input, single_output);

                            y[i] = 0.0; // compute the l-1 norm of the difference
                            for(int j=0; j<num_discrete_nodes; j++)
                                y[i] += std::abs(single_output[j] - data[j]);
                            y[i] = - double(num_discrete_nodes/2) * y[i]; // apply the scale

                            // added for completeness, only logform is used in this example
                            if (sampling_form == TasDREAM::regform) y[i] = std::exp(y[i]);

                        }
                    });
            }

            for(auto &w : workers) w.join();
        };

    // reset the state
    state = TasDREAM::TasmanianDREAM(num_chains, num_dimensions);
    state.setState(initial_chains);

    TasDREAM::SampleDREAM<sampling_form>
                         (num_burnup_iterations, num_sample_iterations,
                          TasDREAM::posterior<sampling_form>
                                (model_likelihood, // provide both the likelihood and the model
                                 TasDREAM::uniform_prior), // assume non-informative prior
                          TasDREAM::hypercube(lower, upper),
                          state,
                          TasDREAM::no_update,
                          TasDREAM::const_percent<100> // use only the differential update
                          );

    solution = state.getApproximateMode();

    //cout << " l-1 acceptance rate: " << state.getAcceptanceRate() << endl;
    cout << "Using l-1 likelihood, the computed solution is:\n"
         << " computed: " << std::fixed;
    for(auto x : solution) cout << setw(13) << x;
    cout << "\n    error: " << std::scientific;
    for(int i=0; i<num_dimensions; i++) cout << setw(13) << std::abs(solution[i] - signal[i]);
    cout << "\n\n";

    // Note: the l-1 likelihood is expected to produce a slightly more accurate solution
    //       but the search is based on random numbers, hence this cannot be guaranteed.
    //       The l-1 likelihood has a much sharper mode, but that reduces the acceptance rate.
    //       Both likelihood methods should identify the signal components to the precision
    //       allowed by the relatively high signal-to-noise ratio.

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [DREAM_Example_05 example]
#endif
}
