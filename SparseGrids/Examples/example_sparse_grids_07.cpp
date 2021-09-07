
#include "Tasmanian.hpp"
#include <chrono>
#include <random>

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_07.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 7.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples7 Tasmanian Sparse Grids module, example 7
 *
 * \par Example 7
 * Benchmark Global vs. Sequence grids.
 */

/*!
 * \ingroup TasmanianSGExamples7
 * \brief Sparse Grids Example 7: Global vs. Sequence grids
 *
 * Sequence rules are those that form the next level by adding one point to the previous one,
 * in Tasmanian those rules are: TasGrid::rule_leja, TasGrid::rule_rleja, TasGrid::rule_rlejashifted,
 * TasGrid::rule_maxlebesgue, TasGrid::rule_minlebesgue, TasGrid::rule_mindelta.
 * Two implementations are provided that deal with such rules, Global grids that use
 * the standard Lagrange polynomial interpolation and Sequence grids that use Newton polynomials.
 * Mathematically the two implementations yield the same result (to within rounding error),
 * but the two processes can have very different computational overhead.
 * Sequence grids offer much faster TasmanianSparseGrid::evaluate() and TasmanianSparseGrid::evaluateBatch()
 * algorithms, at the cost of nearly double the storage and more than double the cost
 * of TasmanianSparseGrid::loadNeededValues().
 *
 * This example serves as a simple demonstration of the difference.
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_07.cpp SG_Example_07 example
 */
void sparse_grids_example_07(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_07 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Example 7: interpolate f(x_1, x_2, x_3, x_4) = exp(-x_1^2 - x_3^2) * exp(x_2) * cos(x_4)\n"
         << "           using rleja rule and comparing Global and Sequence grids\n";

    auto time_start = std::chrono::system_clock::now();
    auto global = TasGrid::makeGlobalGrid(4, 1, 15, TasGrid::type_iptotal, TasGrid::rule_leja);
    auto time_end = std::chrono::system_clock::now();
    long long make_global =
        std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

    time_start = std::chrono::system_clock::now();
    auto sequence = TasGrid::makeSequenceGrid(4, 1, 15, TasGrid::type_iptotal, TasGrid::rule_leja);
    time_end = std::chrono::system_clock::now();
    long long make_sequence =
        std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

    // define batch model
    auto model = [](std::vector<double> const &points)->
        std::vector<double>{
            size_t num_points = points.size() / 4;
            std::vector<double> result(num_points);
            for(size_t i=0; i<num_points; i++){
                double x1 = points[4*i];
                double x2 = points[4*i + 1];
                double x3 = points[4*i + 2];
                double x4 = points[4*i + 3];
                result[i] = std::exp(-x1*x1) * std::cos(x2) * std::exp(-x3*x3) * std::cos(x4);
            }
            return result;
        };

    // load the model values into the grids
    time_start = std::chrono::system_clock::now();
    global.loadNeededValues(model(global.getNeededPoints()));
    time_end = std::chrono::system_clock::now();
    long long load_global =
        std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

    time_start = std::chrono::system_clock::now();
    sequence.loadNeededValues(model(sequence.getNeededPoints()));
    time_end = std::chrono::system_clock::now();
    long long load_sequence =
        std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

    // the benchmark surrogate models is measure from 1000 random reference points
    int const num_test_points = 1000;
    std::vector<double> test_points(4 * num_test_points);
    std::minstd_rand park_miller(42);
    std::uniform_real_distribution<double> domain(-1.0, 1.0);
    for(auto &t : test_points) t = domain(park_miller);

    // reference points
    std::vector<double> reference_result = model(test_points);

    // get the surrogate values at the test points
    time_start = std::chrono::system_clock::now();
    std::vector<double> global_result;
    global.evaluateBatch(test_points, global_result);
    time_end = std::chrono::system_clock::now();
    long long eval_global =
        std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

    time_start = std::chrono::system_clock::now();
    std::vector<double> sequence_result;
    sequence.evaluateBatch(test_points, sequence_result);
    time_end = std::chrono::system_clock::now();
    long long eval_sequence =
        std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

    double global_error = 0.0;
    for(int i=0; i<num_test_points; i++)
        global_error = std::max(global_error,
                                std::abs(global_result[i] - reference_result[i]));

    double sequence_error = 0.0;
    for(int i=0; i<num_test_points; i++)
        sequence_error = std::max(sequence_error,
                                  std::abs(sequence_result[i] - reference_result[i]));

    cout.precision(4);
    cout << std::scientific;
    cout << " Using " << global.getNumPoints() << " points,  "
         << " Global error: " << global_error << ",  Sequence error: " << sequence_error << "\n";

    cout << setw(15) << " " << setw(20) << "Global" << setw(20) << "Sequence" << "\n";
    cout << setw(15) << "make grid" << setw(20) << make_global << setw(20) << make_sequence
         << "  microseconds\n";
    cout << setw(15) << "load values" << setw(20) << load_global << setw(20) << load_sequence
         << "  microseconds\n";
    cout << setw(15) << "evaluate" << setw(20) << eval_global << setw(20) << eval_sequence
         << "  microseconds\n";

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_07 example]
#endif
}
