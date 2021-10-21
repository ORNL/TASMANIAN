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

#ifndef __TASGRID_TESTER_HPP
#define __TASGRID_TESTER_HPP

#include "tasgridLogs.hpp"
#include "tasgridCLICommon.hpp"
#include "tasgridTestFunctions.hpp"

struct TestResults{
    double error;
    int num_points;
};

enum TestType{
    type_integration, type_nodal_interpolation, type_internal_interpolation
};

/*!
 * \internal
 * \brief Enumerated list of all individual sparse grid ctests.
 *
 * \endinternal
 */
enum TestList{
    //! \brief Perform all tests.
    test_all,
    //! \brief Test consistency in the results of all acceleration types.
    test_acceleration,
    //! \brief Test correctness of domain transforms and quadrature weight scaling.
    test_domain,
    //! \brief Test various refinement strategies.
    test_refinement,
    //! \brief Test correctness of global grid algorithms.
    test_global,
    //! \brief Test correctness of local polynomial grid algorithms.
    test_local,
    //! \brief Test correctness of wavelet grids.
    test_wavelet,
    //! \brief Test correctness of Fourier grids.
    test_fourier,
    //! \brief This is not a test, indicates an error in parsing CLI arguments.
    test_none
};

class ExternalTester{
public:
    ExternalTester(int in_num_mc = 1);
    ~ExternalTester();
    void resetRandomSeed();

    void setVerbose(bool new_verbose);
    void setGPUID(int gpu_id);

    //! \brief Converts the string to a TestList, returns \b test_none if not compatible.
    static TestList hasTest(std::string const &s);

    bool Test(TestList test) const;

    bool testGlobalRule(const BaseFunction *f, TasGrid::TypeOneDRule rule, const int *anisotropic, double alpha, double beta, bool interpolation, const int depths[], const double tols[]) const;
    bool performGlobalTest(TasGrid::TypeOneDRule rule) const;
    bool performGaussTransfromTest(TasGrid::TypeOneDRule rule) const;

    bool testLocalPolynomialRule(const BaseFunction *f, TasGrid::TypeOneDRule rule, const int depths[], const double tols[]) const;
    bool testLocalWaveletRule(const BaseFunction *f, const int depths[], const double tols[], bool flavor) const;
    bool testSurplusRefinement(const BaseFunction *f, TasmanianSparseGrid &grid, double tol, TypeRefinement rtype, const int np[], const double errs[], int max_iter ) const;
    bool testAnisotropicRefinement(const BaseFunction *f, TasmanianSparseGrid &grid, TypeDepth type, int min_growth, const int np[], const double errs[], int max_iter ) const;
    bool testDynamicRefinement(const BaseFunction *f, TasmanianSparseGrid &grid, TypeDepth type, double tolerance, TypeRefinement reftype,
                               const std::vector<int> &np, const std::vector<double> &errs) const;
    bool testAcceleration(const BaseFunction *f, TasmanianSparseGrid &grid) const;
    bool testGpuCaching() const;
    bool testGPU2GPUevaluations() const;
    bool testAcceleratedLoadValues(TasGrid::TypeOneDRule rule) const;

    TestResults getError(const BaseFunction *f, TasGrid::TasmanianSparseGrid &grid, TestType type, std::vector<double> const &x = std::vector<double>()) const;

    bool testAllGlobal() const;
    bool testAllPWLocal() const;
    bool testAllWavelet() const;
    bool testAllFourier() const;
    bool testAllRefinement() const;
    bool testAllDomain() const;
    bool testAllAcceleration() const;

    void benchmark(int argc, const char **argv);

    void debugTest(); // call this with -test debug
    void debugTestII(); // call this with -test debug

    static const char* findGaussPattersonTable();

    static const char* testName(TestType type);

private:
    int num_mc;
    bool verbose;
    int gpuid;

    std::vector<TypeAcceleration> available_acc;

    OneOneP3 f11p3;

    TwoOneExpNX2 f21nx2;
    TwoOneCos f21cos;
    TwoOneSinSin f21sinsin;
    TwoOneCosCos f21coscos;
    TwoOneExpSinCos f21expsincos;
    TwoOneSinCosAxis f21sincosaxis;
    TwoTwoSinCos f22sincos;
    TwoOneDivisionAnisotropic f21aniso;
    TwoOne1DCurved f21curved;
    TwoOneExpm40 f21sharp;

    TwoOneConstGC1 f21constGC1;
    TwoOneConstGC2 f21constGC2;
    TwoOneConstGG f21constGG;
    TwoOneConstGJ f21constGJ;
    TwoOneConstGGL f21constGGL;
    TwoOneConstGH f21constGH;

    TwoOneENX2aniso f21nx2aniso;
    TwoTwoExpAsym f22asym;
    TwoOneExpShiftedDomain f21expDomain;
    TwoOneConformalOne f21conformal;
    SixteenOneActive3 f16active3;
    Two3KExpSinCos f23Kexpsincos;
    TwoOneC1C2Periodic f21c1c2periodic;
};

#endif
