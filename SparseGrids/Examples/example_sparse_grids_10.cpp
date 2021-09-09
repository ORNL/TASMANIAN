
#include "Tasmanian.hpp"
#include <random>

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_10.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 10.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples10 Tasmanian Sparse Grids module, example 10
 *
 * \par Example 10
 * Local polynomial vs. Wavelet grids.
 */

/*!
 * \ingroup TasmanianSGExamples10
 * \brief Sparse Grids Example 10: local polynomial vs. wavelet grids
 *
 * Wavelet basis has a number of desirable properties compared to the simple local polynomials,
 * most notably, the basis coefficients are much sharper estimates of the local approximation
 * error. In an adaptive refinement context, this often leads to a decrease in the total
 * number of nodes. However, the wavelets also have large Lebesgue constant especially
 * around the boundary, which can have the converse effect.
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_10.cpp SG_Example_10 example
 */
void sparse_grids_example_10(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_10 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Example 10: comparison between local polynomial and wavelet grids\n\n";

    int const num_inputs = 2;

    // using random points to test the error
    int const num_test_points = 1000;
    std::vector<double> test_points(num_test_points * num_inputs);
    std::minstd_rand park_miller(42);
    std::uniform_real_distribution<double> domain(-1.0, 1.0);
    for(auto &t : test_points) t = domain(park_miller);

    // computes the error between the gird surrogate model and the actual model
    // using the test points, finds the largest absolute error
    auto get_error = [&](TasGrid::TasmanianSparseGrid const &grid,
                         std::function<void(double const x[], double y[], size_t)> model)->
        double{
            std::vector<double> grid_result;
            grid.evaluateBatch(test_points, grid_result);
            double err = 0.0;
            for(int i=0; i<num_test_points; i++){
                double model_result; // using only one output
                model(&test_points[i*num_inputs], &model_result, 0);
                err = std::max(err, std::abs(grid_result[i] - model_result));
            }
            return err;
        };

    // using a model with sharp transition, the modified rules are not
    // suitable for this model, using rule_localp
    auto sharp_model = [](double const x[], double y[], size_t)->
        void{ y[0] = x[0] / (1.0 + 100.0 * std::exp(-10.0 * x[1])); };

    // using order 1
    int const order = 1;
    auto grid_poly = TasGrid::makeLocalPolynomialGrid(num_inputs, 1, 3, order,
                                                         TasGrid::rule_localp);
    auto grid_wavelet = TasGrid::makeWaveletGrid(num_inputs, 1, 1, order);

    // setup the top rows of the output table
    cout << "Using batch refinement:\n"
         << setw(22) << "polynomial" << setw(22) << "wavelet\n"
         << setw(8) << "points" << setw(14) << "error"
         << setw(8) << "points" << setw(14) << "error\n";

    double const tolerance = 1.E-5;
    while((grid_poly.getNumNeeded() > 0) || (grid_wavelet.getNumNeeded() > 0)){
        TasGrid::loadNeededValues(sharp_model, grid_poly, 4);
        TasGrid::loadNeededValues(sharp_model, grid_wavelet, 4);

        // output the grids and results at the current stage
        cout << setw(8) << grid_poly.getNumLoaded()
             << setw(14) << get_error(grid_poly, sharp_model)
             << setw(8) << grid_wavelet.getNumLoaded()
             << setw(14) << get_error(grid_wavelet, sharp_model) << "\n";

        // setting refinement for each grid
        grid_poly.setSurplusRefinement(tolerance, TasGrid::refine_fds);
        grid_wavelet.setSurplusRefinement(tolerance, TasGrid::refine_fds);
    }

    // Compared to local polynomial coefficients, the coefficients of the Wavelet basis
    // are a much sharper estimate of the local error, hence, wavelet based refinement
    // can result in approximation with significantly fewer nodes.
    // However, wavelets have larger Lebesgue constant (especially around the boundary)
    // and thus do not always outperform local polynomials.
    cout << "\n";

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_10 example]
#endif
}
