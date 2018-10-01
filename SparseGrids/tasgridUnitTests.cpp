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

    bool testCover = true;
    bool testExceptions = true;
    bool testAPI = true;
    bool testC = true;

    if ((test == unit_all) || (test == unit_cover)) testCover = testCoverUnimportant();
    if ((test == unit_all) || (test == unit_except)) testExceptions = testAllException();
    if ((test == unit_all) || (test == unit_api)) testAPI = testAPIconsistency();
    if ((test == unit_all) || (test == unit_c)) testC = testCInterface();

    bool pass = testCover && testExceptions && testAPI && testC;
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
    for(int i=0; i<48; i++){
        try{
            invalidArgumentCall(i);
            cout << "Missed arg exception i = " << i << " see GridUnitTester::invalidArgumentCall()" << endl;
            pass = false;
            break;
        }catch(std::invalid_argument &){
            //cout << "Got argument error exception i = " << i << " with message: " << e.what() << endl;
        }
    }

    if (verbose){
        cout << setw(wfirst) << "Exception" << setw(wsecond) << "std::invalid_argument" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    }

    passAll = passAll && pass;
    pass = true;

    // perform std::runtime_error tests
    for(int i=0; i<36; i++){
        try{
            runtimeErrorCall(i);
            cout << "Missed run exception i = " << i << " see GridUnitTester::runtimeErrorCall()" << endl;
            pass = false;
            break;
        }catch(std::runtime_error &){
            //cout << "Got runtime error exception i = " << i << " with message: " << e.what() << endl;
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
    std::vector<int> w;
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
    case 13: grid.makeLocalPolynomialGrid(0,  1,  3,  3, rule_localp); break; // 0 is not valid dimensions
    case 14: grid.makeLocalPolynomialGrid(2, -1,  3,  2, rule_localp); break; // -1 is not valid outputs
    case 15: grid.makeLocalPolynomialGrid(2,  1, -1,  2, rule_localp); break; // -1 is not valid depth
    case 16: grid.makeLocalPolynomialGrid(2,  1,  3, -2, rule_localp); break; // -2 is not a valid order
    case 17: grid.makeLocalPolynomialGrid(2,  1,  3,  2, rule_mindelta); break; // mindelta is not a local rule
    case 18: grid.makeLocalPolynomialGrid(2,  1,  3,  1, rule_localp, std::vector<int>()={3}); break; // level limits is too short
    case 19: grid.makeWaveletGrid(0,  1,  3,  1,  0); break; // 0 is not a valid dimensions
    case 20: grid.makeWaveletGrid(2, -1,  3,  1,  0); break; // -1 is not a valid outputs
    case 21: grid.makeWaveletGrid(2,  1, -3,  1,  0); break; // -3 is not a valid depth
    case 22: grid.makeWaveletGrid(2,  1,  3,  2,  0); break; // 2 is not a valid order (for wavelets)
    case 23: grid.makeWaveletGrid(2,  1,  3,  1,  std::vector<int>()={3}); break; // level limits is too short
    case 24: grid.makeFourierGrid(0, 1, 3, type_level); break; // dimension is 0
    case 25: grid.makeFourierGrid(2, -1, 3, type_level); break; // output is -1
    case 26: grid.makeFourierGrid(2, 2, -1, type_level); break; // depth is -1
    case 27: grid.makeFourierGrid(2, 2, 2, type_level, std::vector<int>()={3}); break; // aw is too short
    case 28: grid.makeFourierGrid(2, 2, 2, type_level, std::vector<int>(), std::vector<int>()={3}); break; // level limits is too short

    case 29: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); grid.updateGlobalGrid(-1, type_level); break; // depth is negative
    case 30: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); grid.updateGlobalGrid(3, type_level, std::vector<int>()={3}); break; // aw is too small
    case 31: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); grid.updateGlobalGrid(3, type_level, std::vector<int>(), std::vector<int>()={3}); break; // ll is too small
    case 32: grid.makeSequenceGrid(2, 1, 3, type_level, rule_rleja); grid.updateSequenceGrid(-1, type_level); break; // depth is negative
    case 33: grid.makeSequenceGrid(2, 1, 3, type_level, rule_rleja); grid.updateSequenceGrid(3, type_level, std::vector<int>()={3}); break; // aw is too small
    case 34: grid.makeSequenceGrid(2, 1, 3, type_level, rule_rleja); grid.updateSequenceGrid(3, type_level, std::vector<int>(), std::vector<int>()={3}); break; // ll is too small

    case 35: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); grid.setAnisotropicRefinement(type_iptotal, -1, 0, std::vector<int>()); break; // min_growth is negative
    case 36: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); gridLoadEN2(&grid); grid.setAnisotropicRefinement(type_iptotal, 1, 2, std::vector<int>()); break; // output out of range
    case 37: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); gridLoadEN2(&grid); grid.setAnisotropicRefinement(type_iptotal, 1, 0, std::vector<int>()={3}); break; // ll is too small

    case 38: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); gridLoadEN2(&grid); grid.estimateAnisotropicCoefficients(type_iptotal, 2, w); break; // output out of range

    case 39: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); gridLoadEN2(&grid); grid.setSurplusRefinement(0.01, 2); break; // output out of range
    case 40: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); gridLoadEN2(&grid); grid.setSurplusRefinement(0.01, 0, std::vector<int>()={3}); break; // ll is too small
    case 41: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); gridLoadEN2(&grid); grid.setSurplusRefinement(-0.1, 0); break; // tolerance is negative

    case 42: grid.makeLocalPolynomialGrid(2, 1, 3); gridLoadEN2(&grid); grid.setSurplusRefinement(0.01, refine_classic, 2); break; // output out of range
    case 43: grid.makeLocalPolynomialGrid(2, 1, 3); gridLoadEN2(&grid); grid.setSurplusRefinement(0.01, refine_classic, 0, std::vector<int>()={3}); break; // ll is too small
    case 44: grid.makeLocalPolynomialGrid(2, 1, 3); gridLoadEN2(&grid); grid.setSurplusRefinement(-0.1, refine_classic, 0); break; // tolerance is negative
    case 45: grid.makeLocalPolynomialGrid(2, 1, 3); gridLoadEN2(&grid); grid.setSurplusRefinement(-0.1, refine_classic, 0, std::vector<int>()={3, 2}, std::vector<double>() = {3.0, 3.0}); break; // scale too small

    case 46: grid.makeLocalPolynomialGrid(2, 1, 3); grid.setDomainTransform(std::vector<double>() = {1.0}, std::vector<double>() = {3.0, 4.0}); break; // a is too small
    case 47: grid.makeLocalPolynomialGrid(2, 1, 3); grid.setDomainTransform(std::vector<double>() = {1.0, 2.0}, std::vector<double>() = {4.0}); break; // b is too small

    default: break;
    }
}

void GridUnitTester::runtimeErrorCall(int i){
    std::vector<double> v, u;
    std::vector<int> w;
    std::vector<int> transformAsin = {4, 4};
    double a[2], b[2];
    TasmanianSparseGrid grid;
    switch(i){
    case  0: grid.updateGlobalGrid(2, type_level); break; // grid not initialized
    case  1: grid.makeSequenceGrid(2, 1, 3, type_level, rule_rleja); grid.updateGlobalGrid(2, type_level); break; // grid not global
    case  2: grid.updateSequenceGrid(2, type_level); break; // grid not initialized
    case  3: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); grid.updateSequenceGrid(2, type_level); break; // grid not sequence
    case  4: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); grid.getInterpolationWeights(std::vector<double>()={0.33}, v); break; // wrong size of x
    case  5: grid.makeGlobalGrid(2, 1, 2, type_level, rule_rleja); grid.loadNeededPoints(std::vector<double>()={0.33, 0.22}); break; // wrong size of loaded data
    case 35: {
             grid.makeGlobalGrid(2, 1, 1, type_level, rule_clenshawcurtis);
             grid.loadNeededPoints(std::vector<double>()={0.33, 0.22, 0.22, 0.22, 0.33});
             grid.loadNeededPoints(std::vector<double>()={0.33, 0.22}); break; // wrong size of loaded data (when overwriting)
             }

    case  6: grid.setAnisotropicRefinement(type_iptotal, 1, 0, 0); break; // grid not init
    case  7: grid.setAnisotropicRefinement(type_iptotal, 1, 0, std::vector<int>()); break; // grid not init
    case  8: grid.makeGlobalGrid(2, 0, 3, type_level, rule_rleja); grid.setAnisotropicRefinement(type_iptotal,  1, 0, std::vector<int>()); break; // no outputs
    case  9: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); grid.setAnisotropicRefinement(type_iptotal,  1, 0, std::vector<int>()); break; // no loaded values
    case 10: grid.makeGlobalGrid(2, 1, 3, type_level, rule_chebyshev); gridLoadEN2(&grid); grid.setAnisotropicRefinement(type_iptotal,  1, 0, std::vector<int>()); break; // rule non-nested
    case 11: grid.makeLocalPolynomialGrid(2, 1, 3); gridLoadEN2(&grid); grid.setAnisotropicRefinement(type_iptotal,  1, 0, std::vector<int>()); break; // grid is localp

    case 12: grid.estimateAnisotropicCoefficients(type_iptotal, 1, w); break; // grid not init
    case 13: grid.makeGlobalGrid(2, 0, 3, type_level, rule_rleja); grid.estimateAnisotropicCoefficients(type_iptotal, 0, w); break; // no outputs
    case 14: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); grid.estimateAnisotropicCoefficients(type_iptotal, 0, w); break; // no loaded values
    case 15: grid.makeGlobalGrid(2, 1, 3, type_level, rule_chebyshev); gridLoadEN2(&grid); grid.estimateAnisotropicCoefficients(type_iptotal, 0, w); break; // rule non-nested
    case 16: grid.makeLocalPolynomialGrid(2, 1, 3); gridLoadEN2(&grid); grid.estimateAnisotropicCoefficients(type_iptotal, 0, w); break; // grid is localp

    case 17: grid.setSurplusRefinement(0.01, 0, 0); break; // grid not init
    case 18: grid.setSurplusRefinement(0.01, 0, std::vector<int>()); break; // grid not init
    case 19: grid.makeGlobalGrid(2, 0, 3, type_level, rule_rleja); grid.setSurplusRefinement(0.01, 0, std::vector<int>()); break; // no outputs
    case 20: grid.makeGlobalGrid(2, 1, 3, type_level, rule_rleja); grid.setSurplusRefinement(0.01, 0, std::vector<int>()); break; // no loaded values
    case 21: grid.makeGlobalGrid(2, 1, 3, type_level, rule_chebyshev); gridLoadEN2(&grid); grid.setSurplusRefinement(0.01, 0, std::vector<int>()); break; // rule non-nested
    case 22: grid.makeLocalPolynomialGrid(2, 1, 3); gridLoadEN2(&grid); grid.setSurplusRefinement(0.01, 0, std::vector<int>()); break; // grid is localp

    case 23: grid.setSurplusRefinement(0.01, refine_classic, 0, 0); break; // grid not init
    case 24: grid.setSurplusRefinement(0.01, refine_classic, 0, std::vector<int>()); break; // grid not init
    case 25: grid.makeLocalPolynomialGrid(2, 0, 3); grid.setSurplusRefinement(0.01, refine_classic, 0); break; // no outputs
    case 26: grid.makeLocalPolynomialGrid(2, 1, 3); grid.setSurplusRefinement(0.01, refine_classic, 0); break; // no loaded values
    case 27: grid.makeGlobalGrid(2, 1, 3, type_level, rule_chebyshev); gridLoadEN2(&grid); grid.setSurplusRefinement(0.01, refine_classic, 0, std::vector<int>()); break; // rule non-local
    case 28: grid.makeGlobalGrid(2, 1, 3, type_level, rule_chebyshev); gridLoadEN2(&grid); grid.setSurplusRefinement(0.01, refine_classic, 0, 0); break; // rule non-local

    case 29: grid.setDomainTransform(a, b); break; // grid is not initialized
    case 30: grid.getDomainTransform(a, b); break; // grid is not initialized
    case 31: grid.setDomainTransform(u, v); break; // grid is not initialized

    case 32: grid.setConformalTransformASIN(transformAsin.data()); break; // grid is not initialized
    case 33: grid.makeGlobalGrid(2, 1, 3, type_level, rule_chebyshev); grid.getConformalTransformASIN(transformAsin.data()); break; // transform not initialized

    case 34: {
             grid.makeGlobalGrid(2, 1, 3, type_level, rule_chebyshev);
             std::vector<int> pntr, indx;
             std::vector<double> vals;
             std::vector<double> x = {-0.33, 0.33};
             grid.evaluateSparseHierarchicalFunctions(x, pntr, indx, vals);
             }

    default: break;
    }
}

void GridUnitTester::gridLoadEN2(TasmanianSparseGrid *grid) const{
    std::vector<double> points;
    grid->getNeededPoints(points);
    int dims = grid->getNumDimensions();
    int outs = grid->getNumOutputs();
    int nump = grid->getNumNeeded();
    std::vector<double> vals(((size_t) nump) * ((size_t) outs));
    auto iter_x = points.begin();
    auto iter_y = vals.begin();
    while(iter_x < points.end()){
        double nrm = 0.0;
        for(int i=0; i<dims; i++){
            nrm += *iter_x * *iter_x;
            iter_x++;
        }
        nrm = exp(-nrm);
        for(int i=0; i<outs; i++) *iter_y++ = nrm;
    }
    grid->loadNeededPoints(vals);
}

bool GridUnitTester::doesMatch(const std::vector<double> &a, const std::vector<double> &b, double prec) const{
    if (a.size() != b.size()) return false;
    auto ib = b.begin();
    for(auto x : a) if (fabs(x - *ib++) > prec) return false;
    return true;
}
bool GridUnitTester::doesMatch(const std::vector<double> &a, const double b[], double prec) const{
    auto ib = b;
    for(auto x : a) if (fabs(x - *ib++) > prec) return false;
    return true;
}
bool GridUnitTester::doesMatch(const std::vector<int> &a, const int b[]) const{
    auto ib = b;
    for(auto x : a) if (x != *ib++) return false;
    return true;
}
bool GridUnitTester::doesMatch(size_t n, double a[], const double b[], double prec) const{
    for(size_t i=0; i<n; i++) if (fabs(a[i] - b[i]) > prec) return false;
    return true;
}

bool GridUnitTester::testAPIconsistency(){
    bool passAll = true;
    int wfirst = 15, wsecond = 30, wthird = 15;

    // test array and vector consistency between the two versions of the API
    bool pass = true;
    double *apoints;
    std::vector<double> vpoints;

    TasmanianSparseGrid grid;
    grid.makeGlobalGrid(2, 1, 4, type_iptotal, rule_clenshawcurtis);
    gridLoadEN2(&grid);
    grid.setAnisotropicRefinement(type_iptotal, 10, 0);

    apoints = grid.getPoints();
    grid.getPoints(vpoints);
    pass = pass && doesMatch(vpoints, apoints);
    delete[] apoints;
    vpoints.clear();

    apoints = grid.getLoadedPoints();
    grid.getLoadedPoints(vpoints);
    pass = pass && doesMatch(vpoints, apoints);
    delete[] apoints;
    vpoints.clear();

    apoints = grid.getNeededPoints();
    grid.getNeededPoints(vpoints);
    pass = pass && doesMatch(vpoints, apoints);
    delete[] apoints;
    vpoints.clear();

    if (verbose) cout << setw(wfirst) << "API variation" << setw(wsecond) << "getPoints()" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    passAll = pass && passAll;

    pass = true;
    std::vector<double> vy, x = {0.333, -0.333};
    double *ay = new double[2];

    grid.evaluate(x, vy);
    grid.evaluate(x.data(), ay);
    pass = pass && doesMatch(vy, ay);
    vy.clear();

    grid.integrate(vy);
    grid.integrate(ay);
    pass = pass && doesMatch(vy, ay);
    vy.clear();
    delete[] ay;

    std::vector<double> vf, vx = {0.333, 0.44, -0.1333, 0.2223};
    double *af = new double[grid.getNumPoints() * 2];
    grid.evaluateHierarchicalFunctions(vx.data(), 2, af);
    grid.evaluateHierarchicalFunctions(vx, vf);
    pass = pass && doesMatch(vf, af);
    delete[] af;

    if (verbose) cout << setw(wfirst) << "API variation" << setw(wsecond) << "evaluate/integrate" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    passAll = pass && passAll;

    std::vector<double> vtransa = {-2.0, 1.0}, vtransb = {1.0, 2.0};
    double atransa[2], atransb[2];

    grid.setDomainTransform(vtransa, vtransb);
    grid.getDomainTransform(atransa, atransb);
    pass = pass && doesMatch(vtransa, atransa) && doesMatch(vtransb, atransb);

    grid.clearDomainTransform();
    grid.getDomainTransform(vtransa, vtransb);
    if (vtransa.size() + vtransb.size() != 0) pass = false;

    if (verbose) cout << setw(wfirst) << "API variation" << setw(wsecond) << "domain transform" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    passAll = pass && passAll;

    std::vector<int> llimits;
    int allimits[3] = {1, 2, 3};
    grid.makeGlobalGrid(3, 2, 5, type_iptotal, rule_fejer2, 0, 0.0, 0.0, 0, allimits);
    grid.getLevelLimits(llimits);
    pass = pass && doesMatch(llimits, allimits) && (llimits.size() == 3);
    grid.clearLevelLimits();
    grid.getLevelLimits(llimits);
    if (llimits.size() != 0) pass = false;

    if (verbose) cout << setw(wfirst) << "API variation" << setw(wsecond) << "level limits" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    passAll = pass && passAll;

    pass = true;
    grid.makeLocalPolynomialGrid(3, 2, 4, rule_localp);
    int *apntr, *aindx;
    double *avals;
    std::vector<int> vpntr, vindx;
    std::vector<double> vvals;
    x = {0.33, 0.33, 0.33, 0.0, 0.44, 0.66, 0.1, -0.33, -0.66};
    grid.evaluateSparseHierarchicalFunctions(x.data(), 3, apntr, aindx, avals);
    grid.evaluateSparseHierarchicalFunctions(x, vpntr, vindx, vvals);
    pass = doesMatch(vpntr, apntr) && doesMatch(vindx, aindx) && doesMatch(vvals, avals);
    delete[] apntr;
    delete[] aindx;
    delete[] avals;
    vpntr.clear();
    vindx.clear();
    vvals.clear();

    if (verbose) cout << setw(wfirst) << "API variation" << setw(wsecond) << "localp sparse basis" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    passAll = pass && passAll;

    pass = true;
    grid.makeWaveletGrid(3, 3, 2, 3);
    x = {0.33, 0.33, 0.33, 0.0, 0.44, 0.66, 0.1, -0.33, -0.66};
    grid.evaluateSparseHierarchicalFunctions(x.data(), 3, apntr, aindx, avals);
    grid.evaluateSparseHierarchicalFunctions(x, vpntr, vindx, vvals);
    pass = doesMatch(vpntr, apntr) && doesMatch(vindx, aindx) && doesMatch(vvals, avals);
    delete[] apntr;
    delete[] aindx;
    delete[] avals;
    vpntr.clear();
    vindx.clear();
    vvals.clear();

    if (verbose) cout << setw(wfirst) << "API variation" << setw(wsecond) << "wavelet sparse basis" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    passAll = pass && passAll;

    // test integer-to-enumerate and string-to-enumerate conversion
    pass = true;
    std::vector<TypeAcceleration> allacc = {accel_none, accel_cpu_blas, accel_gpu_default, accel_gpu_cublas, accel_gpu_cuda, accel_gpu_magma};
    for(auto acc : allacc){
        if (acc != AccelerationMeta::getIOAccelerationString(AccelerationMeta::getIOAccelerationString(acc))){
            cout << "ERROR: mismatch in string to accel conversion: " << AccelerationMeta::getIOAccelerationString(acc) << endl;
            pass = false;
        }
        if (acc != AccelerationMeta::getIOIntAcceleration(AccelerationMeta::getIOAccelerationInt(acc))){
            cout << "ERROR: mismatch in string to accel conversion: " << AccelerationMeta::getIOIntAcceleration(acc) << endl;
            pass = false;
        }
    }

    cout << setw(wfirst+1) << "API variations" << setw(wsecond-1) << "" << setw(wthird) << ((passAll) ? "Pass" : "FAIL") << endl;
    return passAll;
}

bool GridUnitTester::testCInterface(){
    bool pass = (testInterfaceC() != 0);
    int wfirst = 15, wsecond = 30, wthird = 15;
    cout << setw(wfirst+1) << "C interface" << setw(wsecond-1) << "" << setw(wthird) << ((pass) ? "Pass" : "FAIL") << endl;
    return pass;
}

bool GridUnitTester::testCoverUnimportant(){
    // some code is hard/impractical to test automatically, but untested code shows in coverage reports
    // this function gives coverage to such special cases to avoid confusion in the report

    const char *str = TasmanianSparseGrid::getGitCommitHash();
    const char *str2 = TasmanianSparseGrid::getCmakeCxxFlags();
    str = TasmanianSparseGrid::getCmakeCxxFlags();
    if (str[0] != str2[0]){
        cout << "ERROR: mismatch in strings in testCoverUnimportant()" << endl;
        return false;
    }

    #ifndef Tasmanian_ENABLE_CUDA
    AccelerationMeta::cudaCheckError(0, 0);
    AccelerationMeta::cublasCheckError(0, 0);
    AccelerationMeta::cusparseCheckError(0, 0);
    #endif

    std::vector<TypeOneDRule> rules = {rule_none, rule_clenshawcurtis, rule_clenshawcurtis0, rule_chebyshev, rule_chebyshevodd, rule_gausslegendre, rule_gausslegendreodd, rule_gausspatterson, rule_leja, rule_lejaodd, rule_rleja, rule_rlejadouble2, rule_rlejadouble4, rule_rlejaodd, rule_rlejashifted, rule_rlejashiftedeven, rule_rlejashifteddouble, rule_maxlebesgue, rule_maxlebesgueodd, rule_minlebesgue, rule_minlebesgueodd, rule_mindelta, rule_mindeltaodd, rule_gausschebyshev1, rule_gausschebyshev1odd, rule_gausschebyshev2, rule_gausschebyshev2odd, rule_fejer2, rule_gaussgegenbauer, rule_gaussgegenbauerodd, rule_gaussjacobi, rule_gaussjacobiodd, rule_gausslaguerre, rule_gausslaguerreodd, rule_gausshermite, rule_gausshermiteodd, rule_customtabulated, rule_localp, rule_localp0, rule_semilocalp, rule_localpb, rule_wavelet, rule_fourier};
    for(auto r : rules) str = OneDimensionalMeta::getHumanString(r);

    if (!AccelerationMeta::isAccTypeFullMemoryGPU(accel_gpu_default)){
        cout << "ERROR: mismatch in isAccTypeFullMemoryGPU(accel_gpu_default)" << endl;
        return false;
    }
    if (AccelerationMeta::isAccTypeFullMemoryGPU(accel_cpu_blas)){
        cout << "ERROR: mismatch in isAccTypeFullMemoryGPU(accel_cpu_blas)" << endl;
        return false;
    }
    if (!AccelerationMeta::isAccTypeGPU(accel_gpu_default)){
        cout << "ERROR: mismatch in isAccTypeFullMemoryGPU()" << endl;
        return false;
    }

    RuleWavelet rule(1, 10);
    str = rule.getDescription();
    rule.updateOrder(3);
    str = rule.getDescription();

    return true;
}

#endif
