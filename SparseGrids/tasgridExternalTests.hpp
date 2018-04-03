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

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <string.h>
#include <math.h>

#include "TasmanianSparseGrid.hpp"

#include "tasgridTestFunctions.hpp"

using std::cout;
using std::endl;
using std::setw;

using namespace TasGrid;

struct TestResults{
    double error;
    int num_points;
};

enum TestType{
    type_integration, type_nodal_interpolation, type_internal_interpolation
};

enum TestList{
    test_all, test_acceleration, test_domain, test_refinement, test_global, test_local, test_wavelet
};

class ExternalTester{
public:
    ExternalTester(int in_num_mc = 1);
    ~ExternalTester();
    void resetRandomSeed();

    void setVerbose(bool new_verbose);
    void setGPUID(int gpu_id);

    bool Test(TestList test) const;

    bool testGlobalRule(const BaseFunction *f, TasGrid::TypeOneDRule rule, const int *anisotropic, double alpha, double beta, bool interpolation, const int depths[], const double tols[]) const;
    bool performGLobalTest(const TasGrid::TypeOneDRule rule) const;

    bool testLocalPolynomialRule(const BaseFunction *f, TasGrid::TypeOneDRule rule, const int depths[], const double tols[]) const;
    bool testLocalWaveletRule(const BaseFunction *f, const int depths[], const double tols[]) const;
    bool testSurplusRefinement(const BaseFunction *f, TasmanianSparseGrid *grid, double tol, TypeRefinement rtype, const int np[], const double errs[], int max_iter ) const;
    bool testAnisotropicRefinement(const BaseFunction *f, TasmanianSparseGrid *grid, TypeDepth type, int min_growth, const int np[], const double errs[], int max_iter ) const;
    bool testAcceleration(const BaseFunction *f, TasmanianSparseGrid *grid) const;
    bool testGPU2GPUevaluations() const;

    TestResults getError(const BaseFunction *f, TasGrid::TasmanianSparseGrid *grid, TestType type, const double *x = 0) const;

    bool testAllGlobal() const;
    bool testAllPWLocal() const;
    bool testAllWavelet() const;
    bool testAllRefinement() const;
    bool testAllDomain() const;
    bool testAllAcceleration() const;

    void benchmark(int argc, const char **argv);

    void debugTest(); // call this with -test debug
    void debugTestII(); // call this with -test debug

protected:

    void setRandomX(int n, double x[]) const;

private:
    int num_mc;
    bool verbose;
    int gpuid;

    OneOneP3 f11p3;

    TwoOneExpNX2 f21nx2;
    TwoOneCos f21cos;
    TwoOneSinSin f21sinsin;
    TwoOneCosCos f21coscos;
    TwoOneDivisionAnisotropic f21aniso;
    TwoOne1DCurved f21curved;

    TwoOneConstGC1 f21constGC1;
    TwoOneConstGC2 f21constGC2;
    TwoOneConstGG f21constGG;
    TwoOneConstGJ f21constGJ;
    TwoOneConstGGL f21constGGL;
    TwoOneConstGH f21constGH;

    TwoOneENX2aniso f21nx2aniso;
    TwoOneExpShiftedDomain f21expDomain;
    TwoOneConformalOne f21conformal;
    SixteenOneActive3 f16active3;
    Two3KExpSinCos f23Kexpsincos;
};

#endif
