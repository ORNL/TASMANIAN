/*
 * Copyright (c) 2017, Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
 *    and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
 * THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
 * THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
 * IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#include "TasmanianAddons.hpp"
#include "tasgridCLICommon.hpp"

/*!
 * \brief Returns the model values at the grid points.
 *
 * Useful for verification purposes, the unstructured grid algorithm should be equivalent to the standard loadNeededPoints()
 * if applied to the same data, i.e., points aligned to the grid and the corresponding model values.
 *
 * \param grid is the sparse grid that will be used to load the point.
 * \param model is a callable method with the same signature as the array version of loadNeededPoints()
 *
 * \returns the points (first) and model data (second) in a pair of vectors
 */
template<typename CallableModel>
std::pair<std::vector<double>, std::vector<double>> generateExactData(TasmanianSparseGrid const &grid, CallableModel model){
    int num_dimensions = grid.getNumDimensions();
    int num_outputs = grid.getNumOutputs();
    int num_samples = grid.getNumPoints();
    auto x = grid.getPoints();
    std::vector<double> y(Utils::size_mult(grid.getNumOutputs(), grid.getNumPoints()));

    for(int i=0; i<num_samples; i++)
        model(&x[Utils::size_mult(i, num_dimensions)], &y[Utils::size_mult(i, num_outputs)], 0);

    return std::make_pair(std::move(x), std::move(y));
}

/*!
 * \brief Returns the model values at a set of randomly generated points.
 *
 * Similar to generateExactData() but creates a set of samples given by \b num_samples using uniform distribution.
 */
template<typename CallableModel>
std::pair<std::vector<double>, std::vector<double>> generateRandomData(TasmanianSparseGrid const &grid, CallableModel model, int num_samples){
    int num_dimensions = grid.getNumDimensions();
    int num_outputs = grid.getNumOutputs();

    std::vector<double> transform_a, transform_b;
    grid.getDomainTransform(transform_a, transform_b);

    std::minstd_rand park_miller(42);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    auto x = TasDREAM::genUniformSamples(transform_a, transform_b, num_samples, [&]()->double{ return unif(park_miller); });

    std::vector<double> y(Utils::size_mult(num_outputs, num_samples));

    for(int i=0; i<num_samples; i++)
        model(&x[Utils::size_mult(i, num_dimensions)], &y[Utils::size_mult(i, num_outputs)], 0);

    return std::make_pair(std::move(x), std::move(y));
}

/*!
 * \brief Returns the l-infinity difference in the hierarchical coefficients of the two grid.
 */
inline double coefficientDifference(TasmanianSparseGrid const &gA, TasmanianSparseGrid const &gB){
    size_t num_coeffs = Utils::size_mult(gA.getNumOutputs(), gA.getNumPoints());
    double const *c1 = gA.getHierarchicalCoefficients();
    double const *c2 = gB.getHierarchicalCoefficients();
    return std::inner_product(c1, c1 + num_coeffs, c2, 0.0,
                              [](double a, double b)->double{ return std::max(a, b); },
                              [](double a, double b)->double{ return std::abs(a - b); });
}

/*!
 * \brief Returns the l-infinity difference in the outputs of evaluateBatch() of the two grids.
 */
inline double evalDifference(std::vector<double> const&x, TasmanianSparseGrid const &gA, TasmanianSparseGrid const &gB){
    std::vector<double> yA, yB;

    gA.evaluateBatch(x, yA);
    gB.evaluateBatch(x, yB);

    return std::inner_product(yA.begin(), yA.end(), yB.begin(), 0.0,
                              [](double a, double b)->double{ return std::max(a, b); },
                              [](double a, double b)->double{ return std::abs(a - b); });
}

inline double vecDifference(std::vector<double> const &x, std::vector<double> const &y){
    if (x.size() != y.size()) return 1.E+20;
    return std::inner_product(x.begin(), x.end(), y.begin(), 0.0,
                              [](double a, double b)->double{ return std::max(a, b); },
                              [](double a, double b)->double{ return std::abs(a - b); });
}

/*!
 * \brief Performs an exact test on the specified grid and with the given model.
 *
 * The domain of the grid and the model should be [-0.5, 0.7] in three dimensions, but otherwise
 * the inputs are the same as in generateExactData().
 *
 * \throws std::runtime_error if there is a discrepancy between the computed and expected coefficients
 */
template<typename CallableModel>
void testExactL2(TasmanianSparseGrid &&grid, CallableModel model){
    grid.setDomainTransform({-0.5, -0.5, -0.5,}, {0.7, 0.7, 0.7});
    auto ref_grid = grid;

    auto data = generateExactData(grid, model);

    loadNeededPoints(model, ref_grid, 1);
    if (OneDimensionalMeta::isNonNested(grid.getRule()))
        loadUnstructuredDataL2(data.first, data.second, 1.E-3, grid);
    else
        loadUnstructuredDataL2(data.first, data.second, 0.0, grid);

    if ((not OneDimensionalMeta::isNonNested(grid.getRule()) and coefficientDifference(ref_grid, grid) > Maths::num_tol)
        or (OneDimensionalMeta::isNonNested(grid.getRule()) and evalDifference(data.first, ref_grid, grid) > 1.E-2))
    {
        cout << "Failed exact unstructured construction\n";
        cout << "Observed coefficients error: " << coefficientDifference(ref_grid, grid) << "\n";
        cout << "Observed evaluate error: " << evalDifference(data.first, ref_grid, grid) << "\n";
        grid.printStats();
        throw std::runtime_error("test failed");
    }
}

/*!
 * \brief Similar to testExactL2() but the test is done using randomly selected samples.
 */
template<typename CallableModel>
void testApproxL2(TasmanianSparseGrid &&grid, CallableModel model, int num_samples, double tolerance){
    grid.setDomainTransform({-0.5, -0.5, -0.5,}, {0.7, 0.7, 0.7});

    auto data = generateRandomData(grid, model, num_samples);
    loadUnstructuredDataL2(data.first, data.second, 1.E-8, grid);

    std::vector<double> surrogate;
    grid.evaluateBatch(data.first, surrogate);

    double err = vecDifference(data.second, surrogate);
    if (err > tolerance){
        std::cout << std::scientific;
        std::cout.precision(10);
        cout << "Failed approximate unstructured construction\n";
        cout << "Observed error: " << err << "  tolerance: " << tolerance << "\n";
        grid.printStats();
        throw std::runtime_error("test failed");
    }
}

/*!
 * \brief Prepare a list of tests.
 */
std::vector<std::function<void(void)>> makeTests(TypeAcceleration acc, int gpu_id){
    auto model31 = [](double const x[], double y[], int)->
                   void{
                       y[0] = std::exp(-(x[0] - 0.1) * (x[1] - 0.3) * (x[1] - 0.3) * (x[2] + 0.4));
                   };
    auto model33 = [](double const x[], double y[], int)->
                   void{
                       y[0] = std::exp(-(x[0] - 0.1) * (x[1] - 0.3) * (x[1] - 0.3) * (x[2] + 0.4));
                       y[1] = std::exp(-(x[0] - 0.1) * (x[1] - 0.2) * (x[2] + 0.3));
                       y[2] = std::exp(-(x[0] + 0.3) * (x[1] + 0.1) * (x[2] - 0.1));
                   };
    auto set_acc = [=](TasmanianSparseGrid &&grid)->
                   TasmanianSparseGrid{
                       grid.enableAcceleration(acc, gpu_id);
                       return std::move(grid);
                   };
    return std::vector<std::function<void(void)>>{
        [=](void)->void{
            testExactL2(set_acc(makeGlobalGrid(3, 1, 4, type_level, rule_fejer2)), model31);
        },
        [=](void)->void{
            testExactL2(set_acc(makeGlobalGrid(3, 1, 7, type_level, rule_chebyshev)), model31);
        },
        [=](void)->void{
            testExactL2(set_acc(makeSequenceGrid(3, 1, 12, type_ipcurved, rule_mindelta)), model31);
        },
        [=](void)->void{
            testExactL2(set_acc(makeLocalPolynomialGrid(3, 1, 5, 2)), model31);
        },
        [=](void)->void{
            testExactL2(set_acc(makeWaveletGrid(3, 1, 3, 1)), model31);
        },
        [=](void)->void{
            testExactL2(set_acc(makeFourierGrid(3, 1, 6, type_iphyperbolic)), model31);
        },
        [=](void)->void{
            testExactL2(set_acc(makeGlobalGrid(3, 3, 6, type_level, rule_rlejaodd)), model33);
        },
        [=](void)->void{
            testExactL2(set_acc(makeGlobalGrid(3, 3, 7, type_level, rule_chebyshev)), model33);
        },
        [=](void)->void{
            testExactL2(set_acc(makeSequenceGrid(3, 3, 6, type_ipcurved, rule_rlejashifted)), model33);
        },
        [=](void)->void{
            testExactL2(set_acc(makeLocalPolynomialGrid(3, 3, 4, 3)), model33);
        },
        [=](void)->void{
            testExactL2(set_acc(makeWaveletGrid(3, 3, 1, 3)), model33);
        },
        [=](void)->void{
            testExactL2(set_acc(makeFourierGrid(3, 3, 7, type_iphyperbolic)), model33);
        },
        [=](void)->void{
            testApproxL2(set_acc(makeGlobalGrid(3, 3, 5, type_iptotal, rule_fejer2)), model33, 200, 1.E-4);
        },
        [=](void)->void{
            testApproxL2(set_acc(makeSequenceGrid(3, 3, 6, type_ipcurved, rule_rlejashifted)), model33, 100, 1.E-4);
        },
        [=](void)->void{
            testApproxL2(set_acc(makeLocalPolynomialGrid(3, 1, 3, 2, rule_semilocalp)), model31, 120, 1.E-3);
        },
    };
}

/*!
 * \brief Execute the tests and returns the result.
 */
bool runTests(TypeAcceleration acc, bool verbose, int gpu_id){
    bool pass = true;

    auto tests = makeTests(acc, gpu_id);
    int count = 1;
    for(auto &t : tests){
        try{
            t(); // run the test
        }catch(std::runtime_error &){
            pass = false;
        }
        // this one shows which test has failed (for debugging purposes)
        if (verbose) std::cout << " done test: " << std::setw(2) << count++ << "/" << tests.size() << std::endl;
    }

    return pass;
}

/*!
 * \brief All tests for the loadUnstructuredDataL2() method for all acceleration modes.
 */
bool testLoadUnstructuredL2(bool verbose, int gpu_id){
    #ifdef Tasmanian_ENABLE_DPCPP
    test_queue.init_testing();
    #endif
    bool pass = true;

    bool blas_pass = (AccelerationMeta::isAvailable(accel_cpu_blas)) ? runTests(accel_cpu_blas, verbose, 0) : true;
    if (blas_pass){
        if (verbose) cout << std::setw(10) << "blas   " << std::setw(30) << "unstructured construction" << std::setw(10) << "Pass" << endl;
    }else{
        cout << "Failed testLoadUnstructuredL2() blas case.\n";
    }
    pass = pass and blas_pass;

    if (AccelerationMeta::isAvailable(accel_gpu_cuda)){
        int gpu_begin = (gpu_id == -1) ? 0 : gpu_id;
        int gpu_end   = (gpu_id == -1) ? TasmanianSparseGrid::getNumGPUs() : gpu_id + 1;
        for(int gpu = gpu_begin; gpu < gpu_end; gpu++){
            bool cublas_pass = runTests(accel_gpu_cublas, verbose, gpu);
            if (cublas_pass){
                if (verbose) cout << std::setw(7) << "cublas" << std::setw(3) << gpu << std::setw(30) << "unstructured construction" << std::setw(10) << "Pass" << endl;
            }else{
                cout << "Failed testLoadUnstructuredL2() cublas case on device " << gpu << "\n";
            }
            bool cuda_pass = runTests(accel_gpu_cuda, verbose, gpu);
            if (cuda_pass){
                if (verbose) cout << std::setw(7) << "cuda" << std::setw(3) << gpu << std::setw(30) << "unstructured construction" << std::setw(10) << "Pass" << endl;
            }else{
                cout << "Failed testLoadUnstructuredL2() cuda case on device " << gpu << "\n";
            }
            pass = pass and cuda_pass;
            if (AccelerationMeta::isAvailable(accel_gpu_magma)){
                bool magma_pass = runTests(accel_gpu_magma, verbose, gpu);
                if (magma_pass){
                    if (verbose) cout << std::setw(7) << "magma" << std::setw(3) << gpu
                                      << std::setw(30) << "unstructured construction" << std::setw(10) << "Pass" << endl;
                }else{
                    cout << "Failed testLoadUnstructuredL2() magma case on device " << gpu << "\n";
                }
                pass = pass and magma_pass;
            }
        }
    }

    if (verbose and not AccelerationMeta::isAvailable(accel_cpu_blas) and not AccelerationMeta::isAvailable(accel_gpu_cuda))
        cout << "testLoadUnstructuredL2() requires BLAS or CUDA\n";

    return pass;
}
