
#include "Tasmanian.hpp"

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_04.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 4.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples4 Tasmanian Sparse Grids module, example 4
 *
 * \par Example 4
 * Build a surrogate model to an simple function,
 * i.e., construct an interpolant that approximates the function.
 */

/*!
 * \ingroup TasmanianSGExamples4
 * \brief Sparse Grids Example 4: basic surrogate modeling
 *
 * Problems associated with statistical analysis and/or optimization often require
 * a huger number of model realizations (i.e., samples) to compute reliable results.
 * However, when the models are computationally expensive,
 * collecting sufficient number of samples becomes infeasible.
 * Surrogate modeling aims at constructing a cheap to evaluate approximate model
 * from a limited number of samples.
 * The model function used in this example is very simple, but it serves as
 * a demonstration on how to construct interpolatory approximations to
 * much harder models.
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_04.cpp SG_Example_04 example
 */
void sparse_grids_example_04(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_04 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Example 4: interpolate f(x,y) = exp(-x^2) * cos(y),"
         << "using clenshaw-curtis iptotal rule\n";

    const int num_inputs  = 2; // just x and y as inputs
    const int num_outputs = 1; // using single-value output

    const int num_samples = 1; // this is a huge number in an actual application
    std::vector<double> pnts = { 0.3, 0.7 }; // points of interest where the model is not known
    double reference_solution = std::exp(-0.3 * 0.3) * std::cos(0.7);

    std::vector<int> polynomial_order = {6, 12}; // use two grids with different polynomial order

    for(auto prec : polynomial_order){
        // the first 3 inputs to make grid are always: inputs, outputs, depth
        // TasGrid::type_iptotal means interpret the precision as "total degree polynomial space"
        auto grid = TasGrid::makeGlobalGrid(num_inputs, num_outputs, prec,
                                            TasGrid::type_iptotal,
                                            TasGrid::rule_clenshawcurtis);

        // having constructed the initial grid
        // to complete the interpolant, Tasmanian needs the model values at the needed points
        auto points = grid.getNeededPoints();
        // allocate a vector for the associated model values
        std::vector<double> model_values(grid.getNumNeeded() * num_outputs);

        for(int i=0; i<grid.getNumNeeded(); i++){
            double x = points[i*num_inputs];
            double y = points[i*num_inputs+1];
            model_values[i] = std::exp(-x*x) * std::cos(y);
        }

        // feed the data to the grid
        grid.loadNeededValues(model_values);

        // interpolation can be performed after loadNeededPoints() is called
        // note that we don't need to set the size, Tasmanian will automatically set the right size
        std::vector<double> result(num_samples * num_outputs);

        // evaluate() deals with a single point
        // if num_samples is large, evaluateBatch() is much more efficient
        grid.evaluate(pnts, result);

        cout << "\n  using polynomials of total degree up to: " << prec << "\n"
             << "                             the grid has: " << grid.getNumPoints() << " poitns\n"
             << "                 interpolant at (0.3,0.7): " << result[0] << "\n"
             << "                                    error: "
             << std::abs(result[0] - reference_solution) << "\n";
    }

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_04 example]
#endif
}
