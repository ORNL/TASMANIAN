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

#ifndef __TPL_UNIT_TESTS_CPP
#define __TPL_UNIT_TESTS_CPP

#include "tplUnitTests.hpp"
#include "tsgBlasWrappers.hpp"
#include <string>
#include <math.h>
#include <cmath>

bool testTridiagonal();

TplUnitTester::TplUnitTester() : verbose(false){}
TplUnitTester::~TplUnitTester() {}

void TplUnitTester::setVerbose(bool new_verbose){ verbose = new_verbose; }

bool TplUnitTester::Test(UnitTests test){
    if (test == unit_all) {
        cout << endl << endl;
        cout << "---------------------------------------------------------------------" << endl;
        cout << "       Tasmanian Third Party Libraries (TPL): Unit Tests" << endl;
        cout << "---------------------------------------------------------------------" << endl << endl;
        bool pass = true;
        pass = testTriToeplitzEigen();
        cout << endl;
        if (pass){
            cout << "---------------------------------------------------------------------" << endl;
            cout << "           All Unit Tests Completed Successfully" << endl;
            cout << "---------------------------------------------------------------------" << endl << endl;
        } else {
            cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl;
            cout << "         Some Unit Tests Have Failed" << endl;
            cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl << endl;
        }
        return pass;
    } else if (test == unit_none) {
        cout << "No tests to run.\n";
        return true;
    } else {
        cout << "Unknown test type!\n";
        return false;
    }
}

bool testTriToeplitzEigen() {
    // Initialize.
    const int N = 100;
    bool is_matched;
    double exact_eigs[N], D[N], E[N-1], W[N], Z[N*N], WORK1[2 * N - 2], WORK2[4 * N];
    int IBLOCK1[N], ISPLIT1[N], IWORK1[3 * N];
    double a = M_PI_2;
    double b = M_PI_4;
    ISPLIT1[0] = N;
    IBLOCK1[0] = N;
    for (int i=0; i<N; i++) {
        exact_eigs[i] = a + 2.0 * b * cos((i+1) * M_PI / (N+1));
    }
    std::sort(exact_eigs, exact_eigs + N);

    // Test LAPACK's dsterf function.
    std::fill_n(D, N, a);
    std::fill_n(E, N-1, b);
    TasBLAS::sterf(N, D, E);
    std::sort(D, D + N);
    is_matched = doesMatch(N, D, exact_eigs);

    // Test LAPACK's dsteqr function.
    std::fill_n(D, N, a);
    std::fill_n(E, N-1, b);
    TasBLAS::steqr('N', N, D, E, Z, 1, WORK1);
    std::sort(D, D + N);
    is_matched = doesMatch(N, D, exact_eigs);

    // Test LAPACK's dstebz function.
    std::fill_n(D, N, a);
    std::fill_n(E, N-1, b);
    TasBLAS::stebz('A', 'E', N, 0.0, 0.0, 1, N, 1e-13, D, E, N, 1, W, IBLOCK1, ISPLIT1, WORK2, IWORK1);
    is_matched = doesMatch(N, W, exact_eigs);

    // For debugging only.
    // for (int i=0; i<N; i++) {
    //     cout << std::to_string(D[i]) << " " <<  std::to_string(W[i]) << " " << std::to_string(exact_eigs[i])<< endl;
    // }

    return is_matched;
}

bool doesMatch(const std::vector<double> &a, const std::vector<double> &b, double prec) {
    if (a.size() != b.size()) return false;
    auto ib = b.begin();
    for(auto x : a) if (std::abs(x - *ib++) > prec) return false;
    return true;
}
bool doesMatch(const std::vector<double> &a, const double b[], double prec) {
    auto ib = b;
    for(auto x : a) if (std::abs(x - *ib++) > prec) return false;
    return true;
}
bool doesMatch(const std::vector<int> &a, const int b[]) {
    auto ib = b;
    for(auto x : a) if (x != *ib++) return false;
    return true;
}
bool doesMatch(size_t n, const double a[], const double b[], double prec) {
    for(size_t i=0; i<n; i++) if (std::abs(a[i] - b[i]) > prec) return false;
    return true;
}

#endif
