
#include "Tasmanian.hpp"

using namespace std;

/*!
 * \file example_sparse_grids_05.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 5.
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



#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_05 example]
#endif
}
