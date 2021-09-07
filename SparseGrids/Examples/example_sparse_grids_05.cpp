
#include "Tasmanian.hpp"
#include <random>

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_05.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 5.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples5 Tasmanian Sparse Grids module, example 5
 *
 * \par Example 5
 * Build a surrogate model using different adaptive schemes.
 */

/*!
 * \ingroup TasmanianSGExamples5
 * \brief Sparse Grids Example 5: adaptive surrogate modeling
 *
 * Example 4 demonstrates how to build a surrogate model in a direct "one-shot" process.
 * However, it is seldom the case that all model inputs have equal effect on the
 * model outputs, adaptive methods can construct much more reliable surrogates
 * with much fewer model evaluations.
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_05.cpp SG_Example_05 example
 */
void sparse_grids_example_05(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_05 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Example 5: interpolate f(x,y) = exp(-x^2) * cos(y), using leja rule\n"
         << "           employ adaptive refinement to increase accuracy per samples\n";

    // define the model as a C++ lambda expression, the advantage of lambdas
    // is that they can wrap around much more complex models
    // see the documentation for TasGrid::loadNeededValues()
    int const num_inputs  = 2;
    int const num_outputs = 1;
    auto model = [](double const x[], double y[], size_t)->
        void{
            y[0] = std::exp(-x[0] * x[0]) * std::cos(x[1]);
        };

    // the accuracy of the surrogate models is measure from 1000 random reference points
    int const num_test_points = 1000;
    std::vector<double> test_points(num_test_points * num_inputs);
    std::minstd_rand park_miller(42);
    std::uniform_real_distribution<double> domain(-1.0, 1.0);
    for(auto &t : test_points) t = domain(park_miller);

    std::vector<double> reference_values(num_test_points * num_outputs);
    for(int i=0; i<num_test_points; i++)
        model(&test_points[num_inputs * i], &reference_values[num_outputs * i], 0);

    // computes the maximum error between the reference values and the grid interpolant
    auto test_grid = [&](TasGrid::TasmanianSparseGrid const &grid)->
        double{
            std::vector<double> result;
            grid.evaluateBatch(test_points, result);
            double diff = 0.0;
            for(int i=0; i<num_test_points; i++)
                diff = std::max(diff, std::abs(result[i] - reference_values[i]));
            return diff;
        };

    // we compare isotropic grid with 3 different refinement strategies
    int initial_level = 4;
    auto grid_isotropic = TasGrid::makeGlobalGrid(num_inputs, num_outputs, initial_level,
                                                  TasGrid::type_level, TasGrid::rule_leja);
    auto grid_iptotal = grid_isotropic; // using the same initial grid
    auto grid_icurved = grid_isotropic;
    auto grid_surplus = grid_isotropic;

    int budget = 100; // define budget for the number of points for the interpolant

    cout << setw(22) << "isotropic" << setw(22) << "iptotal"
         << setw(22) << "ipcurved" << setw(22) << "surplus\n";
    cout << setw(8) << "points" << setw(14) << "error"
         << setw(8) << "points" << setw(14) << "error"
         << setw(8) << "points" << setw(14) << "error"
         << setw(8) << "points" << setw(14) << "error\n";

    // loop while at least one grid hasn't exhausted the budget
    do{
        if (grid_isotropic.getNumLoaded() < budget){
            TasGrid::loadNeededValues(model, grid_isotropic, 0);
            cout << setw(8) << grid_isotropic.getNumLoaded()
                 << setw(14) << test_grid(grid_isotropic);

            // add an isotropic update until we increase the number of points
            int level = 0;
            while(grid_isotropic.getNumNeeded() == 0)
                grid_isotropic.updateGlobalGrid(++level, TasGrid::type_level);
        }else{
            cout << setw(22) << " ";
        }

        if (grid_iptotal.getNumLoaded() < budget){
            TasGrid::loadNeededValues(model, grid_iptotal, 0);
            cout << setw(8) << grid_iptotal.getNumLoaded() << setw(14) << test_grid(grid_iptotal);

           // set anisotropic total degree update using at least 10 new points
            grid_iptotal.setAnisotropicRefinement(TasGrid::type_iptotal, 10, 0);
        }else{
            cout << setw(22) << " ";
        }

        if (grid_icurved.getNumLoaded() < budget){
            TasGrid::loadNeededValues(model, grid_icurved, 0);
            auto w = grid_icurved.estimateAnisotropicCoefficients(TasGrid::type_ipcurved, 0);
            cout << setw(8) << grid_icurved.getNumLoaded() << setw(14) << test_grid(grid_icurved);

            // set anisotropic curved update using at least 10 new points
            grid_icurved.setAnisotropicRefinement(TasGrid::type_ipcurved, 10, 0);
        }else{
            cout << setw(22) << " ";
        }

        if (grid_surplus.getNumLoaded() < budget){
            TasGrid::loadNeededValues(model, grid_surplus, 0);
            cout << setw(8) << grid_surplus.getNumLoaded() << setw(14) << test_grid(grid_surplus);

            // set surplus based update using tolerance of 1.E-8
            grid_surplus.setSurplusRefinement(1.E-8, 0);
        }else{
            cout << setw(22) << " ";
        }
        cout << "\n";

    }while((grid_isotropic.getNumLoaded() < budget) ||
           (grid_iptotal.getNumLoaded() < budget) ||
           (grid_icurved.getNumLoaded() < budget) ||
           (grid_surplus.getNumLoaded() < budget));

    cout << "\nNote: the surplus refinement is sensitive to the non-monotonic coefficient decay\n"
         << "        and the choice of the tolerance.\n";

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_05 example]
#endif
}
