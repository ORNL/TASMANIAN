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

#ifndef __TASGRID_UNIT_TESTS_CPP
#define __TASGRID_UNIT_TESTS_CPP

#include "tasgridUnitTests.hpp"

GridUnitTester::GridUnitTester() : verbose(false){}
GridUnitTester::~GridUnitTester(){}

void GridUnitTester::setVerbose(bool new_verbose){ verbose = new_verbose; }

bool GridUnitTester::Test(UnitTests test){
    cout << endl << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "       Tasmanian Sparse Grids Module: Unit Tests" << endl;
    cout << "---------------------------------------------------------------------" << endl << endl;

    bool testExceptions = true;

    if ((test == unit_all) || (test == unit_except)) testExceptions = testAllException();

    bool pass = testExceptions;
    //bool pass = true;

    cout << endl;
    if (pass){
        cout << "---------------------------------------------------------------------" << endl;
        cout << "           All Unit Tests Completed Successfully" << endl;
        cout << "---------------------------------------------------------------------" << endl << endl;
    }else{
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl;
        cout << "         Some Unit Tests Have Failed" << endl;
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl << endl;
    }
    return pass;
}

bool GridUnitTester::testAllException(){
    bool pass = true;
    bool passAll = true;
    int wfirst = 15, wsecond = 30, wthird = 15;

    // perform std::invalid_argument tests
    for(int i=0; i<15; i++){
        try{
            invalidArgumentCall(i);
            cout << "Missed arg exception i = " << i << " see GridUnitTester::invalidArgumentCall()" << endl;
            pass = false;
            break;
        }catch(std::invalid_argument &e){
            //cout << "Got arg exception i = " << i << endl;
        }
    }

    if (verbose){
        cout << setw(wfirst) << "Exception" << setw(wsecond) << "std::invalid_argument" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    }

    passAll = passAll && pass;
    pass = true;

    // perform std::runtime_error tests
    for(int i=0; i<2; i++){
        try{
            runtimeErrorCall(i);
            cout << "Missed run exception i = " << i << " see GridUnitTester::runtimeErrorCall()" << endl;
            pass = false;
            break;
        }catch(std::runtime_error &e){
            //cout << "Got arg exception i = " << i << endl;
        }
    }

    if (verbose){
        cout << setw(wfirst) << "Exception" << setw(wsecond) << "std::runtime_error" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    }
    passAll = passAll && pass;

    cout << setw(wfirst+1) << "Exceptions" << setw(wsecond-1) << "" << setw(wthird) << ((passAll) ? "Pass" : "FAIL") << endl;

    return pass;
}

void GridUnitTester::invalidArgumentCall(int i){
    TasmanianSparseGrid grid;
    switch(i){
    case  0: grid.makeGlobalGrid(0, 1, 3, type_level, rule_gausslegendre); break; // dimension is 0
    case  1: grid.makeGlobalGrid(2, -1, 3, type_level, rule_gausslegendre); break; // output is -1
    case  2: grid.makeGlobalGrid(2, 2, -1, type_level, rule_rleja); break; // depth is -1
    case  3: grid.makeGlobalGrid(2, 2, 1, type_level, rule_localp); break; // rule is localp
    case  4: grid.makeGlobalGrid(2, 2, 2, type_level, rule_rleja, std::vector<int>()={3}); break; // aw is too short
    case  5: grid.makeGlobalGrid(2, 2, 2, type_level, rule_customtabulated, std::vector<int>(), 0.0, 0.0, 0); break; // custom filename is empty
    case  6: grid.makeGlobalGrid(2, 2, 2, type_level, rule_chebyshev, std::vector<int>(), 0.0, 0.0, 0, std::vector<int>()={3}); break; // level limits is too short
    case  7: grid.makeSequenceGrid(0, 1, 3, type_level, rule_rleja); break; // dimension is 0
    case  8: grid.makeSequenceGrid(2, -1, 3, type_level, rule_minlebesgue); break; // output is -1
    case  9: grid.makeSequenceGrid(2, 2, -1, type_level, rule_rleja); break; // depth is -1
    case 10: grid.makeSequenceGrid(2, 1, 3, type_level, rule_localp); break; // localp is not a sequence rule
    case 11: grid.makeSequenceGrid(2, 2, 2, type_level, rule_rleja, std::vector<int>()={3}); break; // aw is too short
    case 12: grid.makeSequenceGrid(2, 2, 2, type_level, rule_chebyshev, std::vector<int>(), std::vector<int>()={3}); break; // level limits is too short
    case 13: grid.makeLocalPolynomialGrid(2, 1, 3, -2, rule_localp); break; // -2 is not a valid order
    case 14: grid.makeWaveletGrid(2, 1, 3, 2, 0); break; // 2 is not a valid order (for wavelets)
    default: break;
    }
}

void GridUnitTester::runtimeErrorCall(int i){
    TasmanianSparseGrid grid;
    switch(i){
    case 0: grid.updateGlobalGrid(2, type_level); break; // grid not initialized
    case 1: grid.updateSequenceGrid(2, type_level); break; // grid not initialized
    default: break;
    }
}

#endif
