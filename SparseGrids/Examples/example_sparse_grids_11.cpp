
#include "Tasmanian.hpp"
#include <random>

using namespace std;

/*!
 * \internal
 * \file example_sparse_grids_11.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * Tasmanian Sparse Grids Example 11.
 * \endinternal
 */

/*!
 * \ingroup TasmanianSGExamples
 * \addtogroup TasmanianSGExamples11 Tasmanian Sparse Grids module, example 11
 *
 * \par Example 11
 * Constructing a grid from unstructured data.
 */

/*!
 * \ingroup TasmanianSGExamples11
 * \brief Sparse Grids Example 11: Unstructured data
 *
 * Sparse grid approximation (or surrogates) can be constructed from a set of points
 * not necessarily aligned to the points on the grid, using a least-squares fit.
 * - the fit will not interpolate the data, i.e., there will be a difference between
 *   the surrogate and the data at each point
 * - the overall accuracy will be lower and more points are needed to reach the same
 *   accuracy as in the structured case, the fit should be used when the model
 *   cannot be sampled exactly at the sparse grid points
 * - solving the fit requires the solution to a system of linear equations and
 *   uses significant amount of flops and memory, at the minimum BLAS must be enabled
 *   but GPU acceleration is preferable with the MAGMA library which provides
 *   advanced out-of-core methods
 *
 * \snippet SparseGrids/Examples/example_sparse_grids_11.cpp SG_Example_11 example
 */
void sparse_grids_example_11(){
#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_11 example]
#endif

    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Example 11: construction using unstructured data\n\n";

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

    // using a simple model
    auto model = [](double const x[], double y[], size_t)->
        void{ y[0] = std::exp(-x[0]*x[0] -x[1]-x[1]); };

    auto grid = TasGrid::makeGlobalGrid(num_inputs, 1, 4, TasGrid::type_level, TasGrid::rule_clenshawcurtis);

    // generate random data for the inputs, and compute the corresponding outputs
    int const num_data_points = 2000;
    std::vector<double> data_input(num_inputs * num_data_points);
    std::vector<double> data_output(num_data_points);

    for(auto &d : data_input) d = domain(park_miller);
    for(int i=0; i<num_data_points; i++) model(&data_input[i * num_inputs], &data_output[i], 0);

    // check if capability is available
    if (not grid.isAccelerationAvailable(TasGrid::accel_cpu_blas) and
        not grid.isAccelerationAvailable(TasGrid::accel_gpu_cuda) and
        not grid.isAccelerationAvailable(TasGrid::accel_gpu_magma)){
        cout << "Skipping example 11, BLAS, CUDA, or MAGMA acceleration required.\n";
        return;
    }

    // accel_cpu_blas: works on the CPU and can utilize all available RAM
    // accel_gpu_cuda: works on the GPU but it is restricted to the case
    //                 where the data can fit in GPU memory
    // accel_gpu_magma: works out-of-core, the data is stored in CPU RAM
    //                  while computations are still done on the GPU
    grid.enableAcceleration(TasGrid::accel_gpu_magma);

    // constructs the grid, depending on the amount of data data,
    // the side of the grid and the enabled acceleration
    // this can take significant amount of time
    TasGrid::loadUnstructuredDataL2(data_input, data_output, 1.E-4, grid);

    cout << "Using construction from unstructured (random) data\n";
    cout << "   approximatino error = " << get_error(grid, model) << "\n\n";

#ifndef __TASMANIAN_DOXYGEN_SKIP
//! [SG_Example_11 example]
#endif
}
