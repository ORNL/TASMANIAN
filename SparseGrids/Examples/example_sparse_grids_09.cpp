
#include "Tasmanian.hpp"
#include <random>

using namespace std;

/*!
 * \file example_sparse_grids_09.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 9.
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples9 Tasmanian Sparse Grids module, example 9
 *
 * \par Example 9
 * Different local polynomial refinement.
 */

/*!
 * \ingroup TasmanianSGExamples9
 * \brief Sparse Grids Example 9: local polynomial refinement
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_09.cpp SG_Example_09 example
 */
void sparse_grids_example_09(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_09 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Example 9: comparison between local polynomial refinement strategies\n\n";

    int const num_inputs = 2;

    // using random points to test the error
    int const num_test_points = 1000;
    std::vector<double> test_points(num_test_points * num_inputs);
    std::minstd_rand park_miller(42);
    std::uniform_real_distribution<double> domain(-1.0, 1.0);
    for(auto &t : test_points) t = domain(park_miller);

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_09 example]
#endif
}
