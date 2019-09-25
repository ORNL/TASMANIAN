
#include "Tasmanian.hpp"

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_01.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 1.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples1 Tasmanian Sparse Grids module, example 1
 *
 * \par Example 1
 * Compute \f$ \int_{-1}^{+1} \int_{-1}^{+1} \exp(-x^2) \cos(y) dy dx \f$
 * using a sparse grid with Clenshaw-Curtis nodes and weights.
 */

/*!
 * \ingroup TasmanianSGExamples1
 * \brief Sparse Grids Example 1: integrate a simple function over a canonical domain.
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_01.cpp SG_Example_01 example
 */
void sparse_grids_example_01(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_01 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(17);
    cout << "Example 1:  integrate f(x,y) = exp(-x^2) * cos(y),\n"
         << "            using clenshaw-curtis nodes and grid of type level\n";

    int dimension = 2;
    int level = 6;

    TasGrid::TasmanianSparseGrid grid;
    grid.makeGlobalGrid(dimension, 0, level, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
    std::vector<double> points  = grid.getPoints();
    std::vector<double> weights = grid.getQuadratureWeights();
    int num_points = grid.getNumPoints(); // also equal to weights.size()

    double I = 0.0;
    for(int i=0; i<num_points; i++){
        double x = points[i*dimension];
        double y = points[i*dimension+1];
        I += weights[i] * std::exp(-x*x) * std::cos(y);
    }

    double exact = 2.513723354063905e+00;
    double E = std::abs(exact - I);
    cout <<   "      at level: " << level
         << "\n      the grid has: " << num_points
         << "\n      integral: " << I
         << "\n         error: " << E << "\n\n";

    level = 7;

    grid.makeGlobalGrid(dimension, 0, level, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
    points = grid.getPoints();
    weights = grid.getQuadratureWeights();
    num_points = grid.getNumPoints();

    I = 0.0;
    for(int i=0; i<num_points; i++){
        double x = points[i*dimension];
        double y = points[i*dimension+1];
        I += weights[i] * std::exp(-x*x) * std::cos(y);
    }

    E = std::abs(exact - I);
    cout <<   "      at level: " << level
         << "\n      the grid has: " << num_points
         << "\n      integral: " << I
         << "\n         error: " << E << endl;

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_01 example]
#endif
}
