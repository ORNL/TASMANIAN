
#include "Tasmanian.hpp"

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_02.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 2.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples2 Tasmanian Sparse Grids module, example 2
 *
 * \par Example 2
 * Compute \f$ \int_{-5}^{+5} \int_{-2}^{+3} \exp(-x^2) \cos(y) dy dx \f$
 * using a sparse grid with Gauss-Patterson nodes and weights.
 * The grid is constructed over a non-canononical domain and the points
 * are selected based on a total degree polynomial space.
 */

/*!
 * \ingroup TasmanianSGExamples2
 * \brief Sparse Grids Example 2: integrate a simple function over a canonical domain.
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_02.cpp SG_Example_02 example
 */
void sparse_grids_example_02(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_02 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(17);
    cout << "Example 2: integrate f(x,y) = exp(-x^2) * cos(y) over [-5,5] x [-2,3]\n"
         << "           using  Gauss-Patterson nodes and total degree polynomial space\n";

    int dimension = 2;
    int exactness = 20;

    // the type_qptotal will guarantee exact integral for all polynomials with degree 20 or less
    auto grid = TasGrid::makeGlobalGrid(dimension, 0, exactness,
                                        TasGrid::type_qptotal, TasGrid::rule_gausspatterson);
    grid.setDomainTransform({-5.0, -2.0}, {5.0, 3.0}); // set the non-canonical domain
    auto points  = grid.getPoints();
    auto weights = grid.getQuadratureWeights();

    double I = 0.0;
    for(int i=0; i<grid.getNumPoints(); i++){
        double x = points[i*dimension];
        double y = points[i*dimension+1];
        I += weights[i] * std::exp(-x*x) * std::cos(y);
    }

    double exact = 1.861816427518323e+00;
    double E = std::abs(exact - I);
    cout <<   "      at polynomial exactness: " << exactness
         << "\n      the grid has: " << weights.size()
         << "\n      integral: " << I
         << "\n         error: " << E << "\n\n";


    exactness = 40;

    grid = TasGrid::makeGlobalGrid(dimension, 0, exactness,
                                   TasGrid::type_qptotal, TasGrid::rule_gausspatterson);
    grid.setDomainTransform({-5.0, -2.0}, {5.0, 3.0}); // set the non-canonical domain
    points  = grid.getPoints();
    weights = grid.getQuadratureWeights();

    I = 0.0;
    for(int i=0; i<grid.getNumPoints(); i++){
        double x = points[i*dimension];
        double y = points[i*dimension+1];
        I += weights[i] * std::exp(-x*x) * std::cos(y);
    }

    E = std::abs(exact - I);
    cout <<   "      at polynomial exactness: " << exactness
         << "\n      the grid has: " << weights.size()
         << "\n      integral: " << I
         << "\n         error: " << E << endl;

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_02 example]
#endif
}
