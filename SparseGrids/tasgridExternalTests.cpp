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

#ifndef __TASGRID_TESTER_CPP
#define __TASGRID_TESTER_CPP

#include "tasgridExternalTests.hpp"

#include "tsgCudaMacros.hpp"

ExternalTester::ExternalTester(int in_num_mc) : num_mc(in_num_mc), verbose(false), gpuid(-1){ srand(10); }
ExternalTester::~ExternalTester(){}
void ExternalTester::resetRandomSeed(){ srand((int) time(0)); }

void ExternalTester::setVerbose(bool new_verbose){ verbose = new_verbose; }
void ExternalTester::setGPUID(int gpu_id){ gpuid = gpu_id; }

void ExternalTester::setRandomX(int n, double x[]) const{
    for(int i=0; i<n; i++)
        x[i] = 2.0 * ((double) rand()) / ((double) RAND_MAX) -1.0;
}

bool ExternalTester::Test(TestList test) const{
    cout << endl << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "       Tasmanian Sparse Grids Module: Functionality Test" << endl;
    cout << "---------------------------------------------------------------------" << endl << endl;

    bool passAccel   = true;
    bool passDomain  = true;
    bool passRefine  = true;
    bool passGlobal  = true;
    bool passLocal   = true;
    bool passWavelet = true;
    bool passFourier = true;

    if ((test == test_all) || (test == test_acceleration)) passAccel   = testAllAcceleration();
    if ((test == test_all) || (test == test_domain))       passDomain  = testAllDomain();
    if ((test == test_all) || (test == test_refinement))   passRefine  = testAllRefinement();
    if ((test == test_all) || (test == test_global))       passGlobal  = testAllGlobal();
    if ((test == test_all) || (test == test_local))        passLocal   = testAllPWLocal();
    if ((test == test_all) || (test == test_wavelet))      passWavelet = testAllWavelet();
    if ((test == test_all) || (test == test_fourier))      passFourier = testAllFourier();

    bool pass = passGlobal && passLocal && passWavelet && passFourier && passRefine && passDomain && passAccel;
    //bool pass = true;

    cout << endl;
    if (pass){
        cout << "---------------------------------------------------------------------" << endl;
        cout << "           All Tests Completed Successfully" << endl;
        cout << "---------------------------------------------------------------------" << endl << endl;
    }else{
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl;
        cout << "         Some Tests Have Failed" << endl;
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl << endl;
    }
    return pass;
}

TestResults ExternalTester::getError(const BaseFunction *f, TasGrid::TasmanianSparseGrid *grid, TestType type, const double *x) const{
    TestResults R;
    int num_dimensions = f->getNumInputs();
    int num_outputs = f->getNumOutputs();
    int num_points = grid->getNumPoints();
    if ((type == type_integration) || (type == type_nodal_interpolation)){
        std::vector<double> points, weights;
        grid->getPoints(points);
        if (type == type_integration){
            grid->getQuadratureWeights(weights);
        }else{
            std::vector<double> vx(f->getNumInputs());
            std::copy(x, x + f->getNumInputs(), vx.data());
            grid->getInterpolationWeights(vx, weights);
        }

        std::vector<double> y(num_outputs);
        std::vector<double> r(num_outputs, 0.0);
//      Sequential version: integration
        for(int i=0; i<num_points; i++){
            f->eval(&(points[i*num_dimensions]), y.data());
            for(int k=0; k<num_outputs; k++) r[k] += weights[i] * y[k];
        }

        double err = 0.0;
        if (type == type_integration){
            f->getIntegral(y.data());
        }else{
            f->eval(x, y.data());
        }
        for(int j=0; j<num_outputs; j++){
            err += fabs(y[j] - r[j]);
        };
        R.error = err;
    }else if (type == type_internal_interpolation){
        // load needed points
        int num_needed_points = grid->getNumNeeded();
        if (num_needed_points > 0){
            std::vector<double> values(num_outputs * num_needed_points), needed_points;
            grid->getNeededPoints(needed_points);

            for(int i=0; i<num_needed_points; i++){
                f->eval(&(needed_points[i*num_dimensions]), &(values[i*num_outputs]));
            }

            grid->loadNeededPoints(values);
        }

        std::vector<double> err(num_outputs, 0.0); // absolute error
        std::vector<double> nrm(num_outputs, 0.0); // norm, needed to compute relative error

        std::vector<double> test_x(num_mc * num_dimensions);
        std::vector<double> result_tasm(num_mc * num_outputs);
        std::vector<double> result_true(num_mc * num_outputs);
        setRandomX(num_dimensions * num_mc, test_x.data());

		#pragma omp parallel for // note that iterators do not work with OpenMP, direct indexing does
		for(int i=0; i<num_mc; i++){
			grid->evaluate(&(test_x[i * num_dimensions]), &(result_tasm[i * num_outputs]));
			f->eval(&(test_x[i * num_dimensions]), &(result_true[i * num_outputs]));
		}

		for(int i=0; i<num_mc; i++){
			for(int k=0; k<num_outputs; k++){
                double nrmik = fabs(result_true[i * num_outputs + k]);
                double errik = fabs(result_true[i * num_outputs + k] - result_tasm[i * num_outputs + k]);
                if (nrm[k] < nrmik) nrm[k] = nrmik;
                if (err[k] < errik) err[k] = errik;
			}
		}

		double rel_err = 0.0; // relative error
		for(int k=0; k<num_outputs; k++){
            double relative_errork = err[k] / nrm[k];
			if (rel_err < relative_errork) rel_err = relative_errork;
		}

        R.error = rel_err;
    }
    R.num_points = grid->getNumPoints();
    return R;
}

bool ExternalTester::testGlobalRule(const BaseFunction *f, TasGrid::TypeOneDRule rule, const int *anisotropic, double alpha, double beta, const bool interpolation, const int depths[], const double tols[]) const{
    TasGrid::TasmanianSparseGrid grid;
    TestResults R;
    int num_global_tests = (interpolation) ? 3 : 1;
    TestType tests[3] = { type_integration, type_nodal_interpolation, type_internal_interpolation };
    TasGrid::TypeDepth type = (rule == rule_fourier ? TasGrid::type_level : TasGrid::type_iptotal);
    double *x = new double[f->getNumInputs()]; setRandomX(f->getNumInputs(),x);
    if (rule == rule_fourier){ for(int i=0; i<f->getNumInputs(); i++) x[i] = 0.5*(x[i]+1.0); }    // map to canonical [0,1]^d
    bool bPass = true;
    const char *custom_filename = (rule == rule_customtabulated) ? "SparseGrids/GaussPattersonRule.table" : 0;
    for(int i=0; i<num_global_tests; i++){
        if (rule == rule_fourier){
            grid.makeFourierGrid(f->getNumInputs(), ((interpolation) ? f->getNumOutputs() : 0), depths[i], type, anisotropic);
        }else{
            grid.makeGlobalGrid(f->getNumInputs(), ((interpolation) ? f->getNumOutputs() : 0), depths[i], type, rule, anisotropic, alpha, beta, custom_filename);
        }
        R = getError(f, &grid, ((interpolation) ? tests[i] : type_integration), x);
        if (R.error > tols[i]){
            bPass = false;
            cout << setw(18) << "ERROR: FAILED " << (rule == rule_fourier ? "fourier" : "global") << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString(rule);
            if (interpolation){
                if (tests[i%3] == type_integration){
                    cout << setw(25) << "integration test";
                }else if (tests[i%3] == type_nodal_interpolation){
                    cout << setw(25) << "w-interpolation";
                }else{
                    cout << setw(25) << "interpolation";
                }
            }else{
                cout << setw(25) << "integration test";
            }
            cout << "   failed function: " << f->getDescription();
            cout << setw(10) << "observed: " << R.error << "  expected: " << tols[i] << endl;
        }
    }
    if (rule == rule_customtabulated){
        TasGrid::TasmanianSparseGrid *grid_copy;
        for(int i=0; i<num_global_tests; i++){
            grid.makeGlobalGrid(f->getNumInputs(), ((interpolation) ? f->getNumOutputs() : 0), depths[i], type, rule, anisotropic, alpha, beta, custom_filename);
            grid_copy = new TasGrid::TasmanianSparseGrid(grid);
            R = getError(f, grid_copy, ((interpolation) ? tests[i] : type_integration), x);
            if (R.error > tols[i]){
                bPass = false;
                cout << setw(18) << "ERROR: FAILED global" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString(rule);
                if (interpolation){
                    if (tests[i%3] == type_integration){
                        cout << setw(25) << "integration test";
                    }else if (tests[i%3] == type_nodal_interpolation){
                        cout << setw(25) << "w-interpolation";
                    }else{
                        cout << setw(25) << "interpolation";
                    }
                }else{
                    cout << setw(25) << "integration test";
                }
                cout << "   failed function: " << f->getDescription();
                cout << setw(10) << "observed: " << R.error << "  expected: " << tols[i] << endl;
            }
            delete grid_copy;
        }
    }

    if (TasGrid::OneDimensionalMeta::isSequence(rule)){
        for(int i=0; i<num_global_tests; i++){
            if (interpolation){
                grid.makeSequenceGrid(f->getNumInputs(), f->getNumOutputs(), depths[i], type, rule, anisotropic);
                R = getError(f, &grid, tests[i], x);
            }else{
                grid.makeSequenceGrid(f->getNumInputs(), 0, depths[i], type, rule, anisotropic);
                R = getError(f, &grid, type_integration);
            }
            if (R.error > tols[i]){
                bPass = false;
                cout << setw(18) << "ERROR: FAILED sequence" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString(rule);
                if (interpolation){
                    if (tests[i%3] == type_integration){
                        cout << setw(25) << "integration test";
                    }else if (tests[i%3] == type_nodal_interpolation){
                        cout << setw(25) << "w-interpolation";
                    }else{
                        cout << setw(25) << "interpolation";
                    }
                }else{
                    cout << setw(25) << "integration test";
                }
                cout << "   failed function: " << f->getDescription();
                cout << setw(10) << "observed: " << R.error << "  expected: " << tols[i] << endl;
            }
        }
    }
    delete[] x;
    return bPass;
}

bool ExternalTester::performGLobalTest(TasGrid::TypeOneDRule rule) const{
    double alpha = 0.3, beta = 0.7;
    bool pass = true;
    int wfirst = 10, wsecond = 35, wthird = 15;
    if (rule == TasGrid::rule_clenshawcurtis){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_clenshawcurtis;
        const int depths1[3] = { 25, 25, 25 };
        const double tols1[3] = { 1.E-12, 1.E-12, 1.E-11 };
        const int depths2[3] = { 25, 27, 27 };
        const double tols2[3] = { 1.E-12, 1.E-10, 1.E-11 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_clenshawcurtis0){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_clenshawcurtis0;
        const int depths1[3] = { 25, 25, 25 };
        const double tols1[3] = { 1.E-12, 1.E-12, 1.E-11 };
        if (testGlobalRule(&f21sinsin, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_chebyshev) || (rule == TasGrid::rule_chebyshevodd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[3] = { 22, 22, 22 };
        const double tols1[3] = { 1.E-12, 1.E-10, 1.E-10 };
        const int depths2[3] = { 22, 22, 22 };
        const double tols2[3] = { 1.E-12, 1.E-09, 1.E-09 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_leja) || (rule == TasGrid::rule_lejaodd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_rleja){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_rleja;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 1.E-08, 1.E-08 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_rlejadouble2) || (rule == TasGrid::rule_rlejadouble4)){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejadouble2;
        const int depths1[3] = { 25, 25, 25 };
        const double tols1[3] = { 1.E-12, 1.E-11, 1.E-11 };
        const int depths2[3] = { 25, 27, 27 };
        const double tols2[3] = { 1.E-12, 1.E-10, 1.E-10 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_rlejaodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejaodd;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_rlejashifted){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejashifted;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 1.E-08, 1.E-08 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_rlejashiftedeven){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejashiftedeven;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 6.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_rlejashifteddouble){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejashifteddouble;
        const int depths1[3] = { 25, 25, 25 };
        const double tols1[3] = { 1.E-12, 1.E-12, 1.E-11 };
        const int depths2[3] = { 25, 27, 27 };
        const double tols2[3] = { 1.E-12, 1.E-10, 1.E-11 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_mindelta) || (rule == TasGrid::rule_mindeltaodd) ||
        (rule == TasGrid::rule_minlebesgue) || (rule == TasGrid::rule_minlebesgueodd) ||
        (rule == TasGrid::rule_maxlebesgue) || (rule == TasGrid::rule_maxlebesgueodd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
        // test the hard-coded sequence values vs the optimizer
        if (rule == rule_minlebesgue){
            int n = 22;
            std::vector<double> minleb, precomputed;
            Optimizer::getGreedyNodes<rule_minlebesgue>(n, minleb);
            Optimizer::getPrecomputedMinLebesgueNodes(precomputed);

            TasGrid::Optimizer::tempFunctional<rule_minlebesgue> g(minleb);
            TasGrid::Optimizer::OptimizerResult R = Optimizer::argMaxGlobal(g);
            if (fabs(R.xmax - precomputed[n]) > 1.E-8){
                pass = false;
                cout << "ERROR: mismatch in stored vs computed nodes for rule_minlebesgue rule" << endl;
            }
        }else if (rule == rule_mindelta){
            int n = 22;
            std::vector<double> mindel, precomputed;
            Optimizer::getGreedyNodes<rule_mindelta>(n, mindel);
            Optimizer::getPrecomputedMinDeltaNodes(precomputed);

            TasGrid::Optimizer::tempFunctional<rule_mindelta> d(mindel);
            TasGrid::Optimizer::OptimizerResult R = Optimizer::argMaxGlobal(d);
            if (fabs(R.xmax - precomputed[n]) > 1.E-9){ // this seems large, double-check
                pass = false;
                cout << "ERROR: mismatch in stored vs computed nodes for rule_mindelta rule" << endl;
            }
        }
    }else if ((rule == TasGrid::rule_gausslegendre) || (rule == TasGrid::rule_gausslegendreodd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[3] = { 20, 36, 38 };
        const double tols1[3] = { 1.E-10, 1.E-07, 1.E-07 };
        const int depths2[3] = { 24, 36, 36 };
        const double tols2[3] = { 1.E-10, 1.E-07, 1.E-07 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausspatterson){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausspatterson;
        const int depths1[3] = { 20, 36, 38 };
        const double tols1[3] = { 1.E-10, 1.E-07, 1.E-07 };
        const int depths2[3] = { 24, 36, 36 };
        const double tols2[3] = { 1.E-10, 1.E-07, 1.E-07 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_customtabulated){
        {
        std::ifstream ftest("SparseGrids/GaussPattersonRule.table");
        if (!ftest.good()){
            ftest.close();
            cout << "WARNING: cannot find GaussPattersonRule.table file and cannot test the custom rule!" << endl;
        }else{
            ftest.close();
            TasGrid::TypeOneDRule oned = TasGrid::rule_customtabulated;
            const int depths1[3] = { 20, 36, 38 };
            const double tols1[3] = { 1.E-10, 1.E-07, 1.E-07 };
            const int depths2[3] = { 24, 36, 36 };
            const double tols2[3] = { 1.E-10, 1.E-07, 1.E-07 };
            if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
                if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
            }else{
                cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
            }
        }}
    }else if (rule == TasGrid::rule_fejer2){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_fejer2;
        const int depths1[3] = { 20, 40, 40 };
        const double tols1[3] = { 1.E-14, 1.E-12, 1.E-12 };
        if (testGlobalRule(&f21coscos, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_gausschebyshev1) || (rule == TasGrid::rule_gausschebyshev1odd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 5.E-14, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGC1, oned, 0, alpha, beta, true, depths1, tols1) && performGaussTransfromTest(oned)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_gausschebyshev2) || (rule == TasGrid::rule_gausschebyshev2odd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 1.1E-14, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGC2, oned, 0, alpha, beta, true, depths1, tols1) && performGaussTransfromTest(oned)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_gaussgegenbauer) || (rule == TasGrid::rule_gaussgegenbauerodd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 1.E-11, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGG, oned, 0, alpha, beta, true, depths1, tols1) && performGaussTransfromTest(oned)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_gaussjacobi) || (rule == TasGrid::rule_gaussjacobiodd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 1.E-08, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGJ, oned, 0, alpha, beta, true, depths1, tols1) && performGaussTransfromTest(oned)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_gausslaguerre) || (rule == TasGrid::rule_gausslaguerreodd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[1] = { 20 };
        const double tols1[1] = { 1.E-08 };
        if (testGlobalRule(&f21constGGL, oned, 0, alpha, beta, false, depths1, tols1) && performGaussTransfromTest(oned)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if ((rule == TasGrid::rule_gausshermite) || (rule == TasGrid::rule_gausshermiteodd)){
        { TasGrid::TypeOneDRule oned = rule;
        const int depths1[1] = { 20 };
        const double tols1[1] = { 1.E-09 };
        if (testGlobalRule(&f21constGH, oned, 0, alpha, beta, false, depths1, tols1) && performGaussTransfromTest(oned)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }
    return pass;
}

bool ExternalTester::performGaussTransfromTest(TasGrid::TypeOneDRule oned) const{
    //double alpha = 0.3, beta = 0.7;
    bool pass = true;
    int wfirst = 10, wsecond = 35, wthird = 15;
    if ((oned == TasGrid::rule_gausschebyshev1) || (oned == TasGrid::rule_gausschebyshev1odd)){
        // Gauss-Chebyshev-1 translated to [4, 7], area = M_PI, integral of f(x) = 1 / x is M_PI * sqrt(7.0) / 14.0
        TasGrid::TasmanianSparseGrid grid;
        grid.makeGlobalGrid(1, 1, 6, type_level, oned);
        double transa = 4.0, transb = 7.0;
        grid.setDomainTransform(&transa, &transb);
        double *w = grid.getQuadratureWeights();
        double *p = grid.getNeededPoints();
        int num_p = grid.getNumNeeded();
        double sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i];
        if (fabs(sum - M_PI) > TSG_NUM_TOL){
            cout << sum << "     " << M_PI << endl;
            cout << "ERROR: sum of weight in transformed gauss-chebyshev-1 rule is off by: " << fabs(sum - M_PI) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i] / p[i];
        //cout << "error in integral of 1/x is = " << fabs(sum - M_PI * sqrt(7.0) / 14.0) << endl;
        if (fabs(sum - M_PI * sqrt(7.0) / 14.0) > 1.E-11){
            cout << "ERROR: disrepancy in transformed gauss-chebyshev-1 rule is: " << fabs(sum - M_PI * sqrt(7.0) / 14.0) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        delete[] w;
        delete[] p;
    }else if ((oned == TasGrid::rule_gausschebyshev2) || (oned == TasGrid::rule_gausschebyshev2odd)){
        // Gauss-Chebyshev-2 translated to [4, 7], area = 9.0 * M_PI / 2.0, integral of f(x) = (7 - x)^0.5 (x - 4)^0.5 is 9.0 / 2.0
        TasGrid::TasmanianSparseGrid grid;
        grid.makeGlobalGrid(1, 1, 10, type_level, oned);
        double transa = 4.0, transb = 7.0;
        grid.setDomainTransform(&transa, &transb);
        double *w = grid.getQuadratureWeights();
        double *p = grid.getNeededPoints();
        int num_p = grid.getNumNeeded();
        double sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i];
        if (fabs(sum - 9.0 * M_PI / 8.0) > TSG_NUM_TOL){
            cout << sum << "     " << 9.0 * M_PI / 8.0 << endl;
            cout << "ERROR: sum of weight in transformed gauss-chebyshev-2 rule is off by: " << fabs(sum - 9.0 * M_PI / 8.0) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i] * sqrt(7.0 - p[i]) * sqrt(p[i] - 4.0);
        //cout << "error in integral of (7 - x)^0.5 (x - 4)^0.5 is = " << fabs(sum - 4.5) << endl;
        if (fabs(sum - 4.5) > 1.E-3){
            cout << "ERROR: disrepancy in transformed gauss-chebyshev-2 rule is: " << fabs(sum - 4.5) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        delete[] w;
        delete[] p;
    }else if ((oned == TasGrid::rule_gaussgegenbauer) || (oned == TasGrid::rule_gaussgegenbauerodd)){
        // Gauss-Gegenbauer translated to [4, 7], area = 8.1, integral of f(x) = x^3 is 389367.0 / 280.0
        TasGrid::TasmanianSparseGrid grid;
        grid.makeGlobalGrid(1, 1, 10, type_level, oned, 0, 2.0);
        double transa = 4.0, transb = 7.0;
        grid.setDomainTransform(&transa, &transb);
        double *w = grid.getQuadratureWeights();
        double *p = grid.getNeededPoints();
        int num_p = grid.getNumNeeded();
        double sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i];
        if (fabs(sum - 8.1) > TSG_NUM_TOL){
            cout << sum << "     " << 8.1 << endl;
            cout << "ERROR: sum of weight in transformed gauss-genebauer rule is off by: " << fabs(sum - 8.1) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i] * p[i] * p[i] * p[i];
        if (fabs(sum - 389367.0 / 280.0) > 1.E-10){
            cout << "ERROR: disrepancy in transformed gauss-gegenbauer rule is: " << fabs(sum - 389367.0 / 280.0) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        delete[] w;
        delete[] p;
    }else if ((oned == TasGrid::rule_gaussjacobi) || (oned == TasGrid::rule_gaussjacobiodd)){
        // Gauss-Jacobi translated to [4, 7], area = 12.15, integral of f(x) = x^3 is 389367.0 / 280.0
        TasGrid::TasmanianSparseGrid grid;
        grid.makeGlobalGrid(1, 1, 10, type_level, oned, 0, 3.0, 2.0);
        double transa = 4.0, transb = 7.0;
        grid.setDomainTransform(&transa, &transb);
        double *w = grid.getQuadratureWeights();
        double *p = grid.getNeededPoints();
        int num_p = grid.getNumNeeded();
        double sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i];
        if (fabs(sum - 12.15) > TSG_NUM_TOL){
            cout << sum << "     " << 12.15 << endl;
            cout << "ERROR: sum of weight in transformed gauss-jacobi rule is off by: " << fabs(sum - 12.15) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i] * sin(M_PI * p[i]);
        if (fabs(sum + 18.0 * (3.0 * M_PI * M_PI - 4.0) / pow(M_PI, 5.0)) > 1.E-11){
            cout << "ERROR: disrepancy in transformed gauss-jacobi rule is: " << fabs(sum + 18.0 * (3.0 * M_PI * M_PI - 4.0) / pow(M_PI, 5.0)) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        delete[] w;
        delete[] p;
    }else if ((oned == TasGrid::rule_gausslaguerre) || (oned == TasGrid::rule_gausslaguerreodd)){
        // Gauss-Laguerre, unbounded domain
        TasGrid::TasmanianSparseGrid grid;
        grid.makeGlobalGrid(2, 1, 6, type_level, oned, 0, 3.0);
        double transa[2] = {4.0, 3.0}, transb[2] = {0.5, 0.75};
        grid.setDomainTransform(transa, transb);
        double *w = grid.getQuadratureWeights();
        double *p = grid.getNeededPoints();
        int num_p = grid.getNumNeeded();
        double sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i];
        if (fabs(sum - 96.0 * 512.0 / 27.0) > TSG_NUM_TOL){
            cout << sum << "     " << 96.0 * 512.0 / 27.0 << endl;
            cout << "ERROR: sum of weight in transformed gauss-laguerre rule is off by: " << fabs(sum - 96.0 * 512.0 / 27.0) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i] * (p[2*i]*p[2*i] * p[2*i+1]*p[2*i+1]*p[2*i+1]);
        if (fabs(sum - 15360.0 * 3573248.0 / 243.0) > 6.E-7){
            cout << "ERROR: disrepancy in transformed gauss-laguerre rule is: " << fabs(sum - 15360.0 * 3573248.0 / 243.0) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        double test_x[2] = {3.0 + sqrt(2.0), 2.0 + sqrt(2.0)};
        double *iw = grid.getInterpolationWeights(test_x);
        sum = 0.0; for(int i=0; i<num_p; i++) sum += iw[i] * (p[2*i]*p[2*i] * p[2*i+1]*p[2*i+1]*p[2*i+1]);
        if (fabs(sum - test_x[0] * test_x[0] * test_x[1] * test_x[1] * test_x[1]) > 2.E-9){
            cout << "ERROR: nodal interpolation using gauss-laguerre: " << fabs(sum - test_x[0] * test_x[0] * test_x[1] * test_x[1] * test_x[1]) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        delete[] p;
        delete[] w;
        delete[] iw;
    }else if ((oned == TasGrid::rule_gausshermite) || (oned == TasGrid::rule_gausshermiteodd)){
        // Gauss-Hermite, unbounded domain
        TasGrid::TasmanianSparseGrid grid;
        grid.makeGlobalGrid(2, 1, 6, type_level, oned, 0, 4.0);
        double transa[2] = {4.0, 3.0}, transb[2] = {0.5, 0.75};
        grid.setDomainTransform(transa, transb);
        double *w = grid.getQuadratureWeights();
        double *p = grid.getNeededPoints();
        int num_p = grid.getNumNeeded();
        double sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i];
        if (fabs(sum - (8.0 * M_PI / 3.0) * sqrt(6.0)) > TSG_NUM_TOL){
            cout << sum << "     " << 96.0 * 512.0 / 27.0 << endl;
            cout << "ERROR: sum of weight in transformed gauss-hermite rule is off by: " << fabs(sum - (8.0 * M_PI / 3.0) * sqrt(6.0)) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        sum = 0.0; for(int i=0; i<num_p; i++) sum += w[i] * (p[2*i]*p[2*i] * p[2*i+1]*p[2*i+1]*p[2*i+1]*p[2*i+1]);
        if (fabs(sum - (63.0 * 19912.0 * M_PI / 81.0) * sqrt(6.0)) > 4.E-8){
            cout << "ERROR: disrepancy in transformed gauss-hermite rule is: " << fabs(sum - 15360.0 * 3573248.0 / 243.0) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        double test_x[2] = {3.0 + sqrt(2.0), 2.0 + sqrt(2.0)};
        double *iw = grid.getInterpolationWeights(test_x);
        sum = 0.0; for(int i=0; i<num_p; i++) sum += iw[i] * (p[2*i]*p[2*i] * p[2*i+1]*p[2*i+1]*p[2*i+1]*p[2*i+1]);
        if (fabs(sum - test_x[0] * test_x[0] * test_x[1] * test_x[1] * test_x[1] * test_x[1]) > 1.E-9){
            cout << "ERROR: nodal interpolation using gauss-hermite: " << fabs(sum - test_x[0] * test_x[0] * test_x[1] * test_x[1] * test_x[1]) << endl;
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }
        delete[] p;
        delete[] w;
        delete[] iw;
    }
    return pass;
}

bool ExternalTester::testAllGlobal() const{
    bool pass = true;
    const int nrules = 36; // sync with below
    TasGrid::TypeOneDRule rules[nrules] = {
                    TasGrid::rule_chebyshev,
                    TasGrid::rule_chebyshevodd,
                    TasGrid::rule_clenshawcurtis,
                    TasGrid::rule_clenshawcurtis0,
                    TasGrid::rule_fejer2,
                    TasGrid::rule_leja,
                    TasGrid::rule_lejaodd,
                    TasGrid::rule_rleja,
                    TasGrid::rule_rlejadouble2,
                    TasGrid::rule_rlejadouble4,
                    TasGrid::rule_rlejaodd,
                    TasGrid::rule_rlejashifted,
                    TasGrid::rule_rlejashiftedeven,
                    TasGrid::rule_rlejashifteddouble,
                    TasGrid::rule_maxlebesgue,
                    TasGrid::rule_maxlebesgueodd,
                    TasGrid::rule_minlebesgue,
                    TasGrid::rule_minlebesgueodd,
                    TasGrid::rule_mindelta,
                    TasGrid::rule_mindeltaodd,
                    TasGrid::rule_gausslegendre,
                    TasGrid::rule_gausslegendreodd,
                    TasGrid::rule_gausspatterson,
                    TasGrid::rule_gausschebyshev1,
                    TasGrid::rule_gausschebyshev1odd,
                    TasGrid::rule_gausschebyshev2,
                    TasGrid::rule_gausschebyshev2odd,
                    TasGrid::rule_gaussgegenbauer,
                    TasGrid::rule_gaussgegenbauerodd,
                    TasGrid::rule_gaussjacobi,
                    TasGrid::rule_gaussjacobiodd,
                    TasGrid::rule_gausslaguerre,
                    TasGrid::rule_gausslaguerreodd,
                    TasGrid::rule_gausshermite,
                    TasGrid::rule_gausshermiteodd,
                    TasGrid::rule_customtabulated   };

    for(int i=0; i<nrules; i++){
        if (!performGLobalTest(rules[i])){
            pass = false;
        }
    }
    { TasGrid::TasmanianSparseGrid grid; grid.makeGlobalGrid(2, 1, 3, TasGrid::type_level, TasGrid::rule_fejer2);
        int *indx = 0, *pntr = 0;
        double *vals = 0;
        double *pnts = new double[20]; setRandomX(20, pnts);
        grid.evaluateSparseHierarchicalFunctions(pnts, 10, pntr, indx, vals);
        getError(&f21nx2, &grid, type_internal_interpolation); // this is done to load the values
        const double *coeff = grid.getHierarchicalCoefficients();
        double *y = new double[10];
        grid.evaluateBatch(pnts, 10, y);
        for(int i=0; i<10; i++){
            for(int j=pntr[i]; j<pntr[i+1]; j++){
                y[i] -= coeff[indx[j]] * vals[j];
            }
        }
        for(int i=0; i<10; i++){
            if (fabs(y[i]) > TSG_NUM_TOL){
                cout << "Error in evaluateSparseHierarchicalFunctions() (global)" << endl;
                pass = false;
            }
        }
        delete[] indx; delete[] pntr; delete[] vals; delete[] pnts; delete[] y;
    }
    int wfirst = 11, wsecond = 34, wthird = 15;
    if (pass){
        cout << setw(wfirst) << "Rules" << setw(wsecond) << "global/sequence" << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Rules" << setw(wsecond) << "global/sequence" << setw(wthird) << "FAIL" << endl;
    }
    return pass;
}

bool ExternalTester::testLocalPolynomialRule(const BaseFunction *f, TasGrid::TypeOneDRule rule, const int depths[], const double tols[]) const{
    TasGrid::TasmanianSparseGrid grid;
    TestResults R;
    TestType tests[3] = { type_integration, type_nodal_interpolation, type_internal_interpolation };
    int orders[6] = { 0, 1, 2, 3, 4, -1 };
    double *x = new double[f->getNumInputs()]; setRandomX(f->getNumInputs(),x);
    bool bPass = true;
    for(int i=0; i<18; i++){
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), depths[i], orders[i/3], rule);
        R = getError(f, &grid, tests[i%3], x);
        if (R.error > tols[i]){
            bPass = false;
            cout << setw(18) << "ERROR: FAILED ";
            cout << setw(6) << TasGrid::OneDimensionalMeta::getIORuleString(rule);
            cout << " order: " << orders[i/3];

            if (tests[i%3] == type_integration){
                cout << setw(25) << "integration test";
            }else if (tests[i%3] == type_nodal_interpolation){
                cout << setw(25) << "w-interpolation";
            }else{
                cout << setw(25) << "interpolation";
            }

            cout << "   failed function: " << f->getDescription();
            cout << setw(10) << "observed: " << R.error << "  expected: " << tols[i] << endl;
        }
    }
    delete[] x;
    return bPass;
}

bool ExternalTester::testSurplusRefinement(const BaseFunction *f, TasmanianSparseGrid *grid, double tol, TypeRefinement rtype, const int np[], const double errs[], int max_iter ) const{
    for(int itr=0; itr<max_iter; itr++){
        TestResults R = getError(f, grid, type_internal_interpolation);
        if ( (R.num_points != np[itr]) || (R.error > errs[itr]) ){
            cout << setw(18) << "ERROR: FAILED refinement test at iteration: " << itr << endl;
            cout << " expected: " << np[itr] << "  " << errs[itr] << "   computed: " << R.num_points << "  " << R.error << endl;
            return false;
        }
        if (grid->isGlobal()){
            grid->setSurplusRefinement(tol, 0);
        }else if (grid->isSequence()){
            grid->setSurplusRefinement(tol, -1);
        }else{
            grid->setSurplusRefinement(tol, rtype);
        }
    }
    return true;
}
bool ExternalTester::testAnisotropicRefinement(const BaseFunction *f, TasmanianSparseGrid *grid, TypeDepth type, int min_growth, const int np[], const double errs[], int max_iter ) const{
    for(int itr=0; itr<max_iter; itr++){
        TestResults R = getError(f, grid, type_internal_interpolation);
        if ( (R.num_points != np[itr]) || (R.error > errs[itr]) ){
            cout << setw(18) << "ERROR: FAILED refinement test at iteration: " << itr << endl;
            cout << " expected: " << np[itr] << "  " << errs[itr] << "   computed: " << R.num_points << "  " << R.error << endl;
            return false;
        }
        if (grid->isGlobal()){
            grid->setAnisotropicRefinement(type, min_growth, 0);
        }else{
            grid->setAnisotropicRefinement(type, min_growth, -1);
        }
    }
    return true;
}

bool ExternalTester::testAllPWLocal() const{
    bool pass = true;
    int wfirst = 10, wsecond = 35, wthird = 15;

    { TasGrid::TypeOneDRule oned = TasGrid::rule_semilocalp;
    const int depths1[18] = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
    const double tols1[18] = { 1.E-03, 5.E-01, 5.E-01, 1.E-03, 1.E-03, 1.E-03, 1.E-07, 1.E-04, 1.E-04, 1.E-07, 1.E-05, 1.E-05, 1.E-07, 4.E-06, 4.E-06, 1.E-07, 4.E-06, 4.E-06 };
    if (testLocalPolynomialRule(&f21nx2, oned, depths1, tols1)){
        if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl; pass = false;
    }}
    { TasGrid::TypeOneDRule oned = TasGrid::rule_localp;
    const int depths1[18] = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
    const double tols1[18] = { 1.E-03, 5.E-01, 5.E-01, 1.E-03, 1.E-03, 1.E-03, 1.E-07, 1.E-04, 1.E-04, 1.E-07, 1.E-05, 1.E-05, 1.E-07, 4.E-06, 4.E-06, 1.E-07, 4.E-06, 4.E-06 };
    if (testLocalPolynomialRule(&f21nx2, oned, depths1, tols1)){
        if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl; pass = false;
    }}
    { TasGrid::TypeOneDRule oned = TasGrid::rule_localpb;
    const int depths1[18] = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
    const double tols1[18] = { 1.E-03, 5.E-01, 5.E-01, 1.E-03, 1.E-03, 1.E-03, 1.E-07, 1.E-04, 1.E-04, 1.E-07, 1.E-05, 1.E-05, 1.E-07, 9.E-06, 9.E-06, 1.E-06, 2.E-05, 2.E-05 };
    if (testLocalPolynomialRule(&f21sincosaxis, oned, depths1, tols1)){
        if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl; pass = false;
    }}
    { TasGrid::TypeOneDRule oned = TasGrid::rule_localp0;
    const int depths1[18] = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
    const double tols1[18] = { 1.E-03, 5.E-01, 5.E-01, 1.E-03, 2.E-04, 2.E-04, 1.E-09, 1.E-06, 1.E-06, 1.E-09, 3.E-08, 3.E-08, 1.E-09, 4.E-09, 4.E-09, 1.E-09, 4.E-09, 4.E-09 };
    if (testLocalPolynomialRule(&f21coscos, oned, depths1, tols1)){
        if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl; pass = false;
    }}
    { TasGrid::TasmanianSparseGrid grid; grid.makeLocalPolynomialGrid(2, 1, 4, 1);
        int *indx = 0, *pntr = 0;
        double *vals = 0;
        double *pnts = new double[20]; setRandomX(20, pnts);
        grid.evaluateSparseHierarchicalFunctions(pnts, 10, pntr, indx, vals);
        getError(&f21nx2, &grid, type_internal_interpolation); // this is done to load the values
        const double *coeff = grid.getHierarchicalCoefficients();
        double *y = new double[10];
        grid.evaluateBatch(pnts, 10, y);
        for(int i=0; i<10; i++){
            for(int j=pntr[i]; j<pntr[i+1]; j++){
                y[i] -= coeff[indx[j]] * vals[j];
            }
        }
        for(int i=0; i<10; i++){
            if (fabs(y[i]) > TSG_NUM_TOL){
                cout << "Error in evaluateSparseHierarchicalFunctions() (localp)" << endl;
                pass = false;
            }
        }
        delete[] indx; delete[] pntr; delete[] vals; delete[] pnts; delete[] y;
    }
    wfirst = 11; wsecond = 34;
    if (pass){
        cout << setw(wfirst) << "Rules" << setw(wsecond) << "local polynomial" << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Rules" << setw(wsecond) << "local polynomial" << setw(wthird) << "FAIL" << endl;
    }
    return pass;
}

bool ExternalTester::testLocalWaveletRule(const BaseFunction *f, const int depths[], const double tols[]) const{
    TasGrid::TasmanianSparseGrid grid;
    TestResults R;
    TestType tests[3] = { type_integration, type_nodal_interpolation, type_internal_interpolation };
    int orders[2] = { 1, 3 };
    double *x = new double[f->getNumInputs()]; setRandomX(f->getNumInputs(),x);
    bool bPass = true;
    for(int i=0; i<6; i++){
        grid.makeWaveletGrid(f->getNumInputs(), f->getNumOutputs(), depths[i], orders[i/3]);
        R = getError(f, &grid, tests[i%3], x);
        if (R.error > tols[i]){
            bPass = false;
            cout << setw(18) << "ERROR: FAILED";
            cout << setw(6) << TasGrid::OneDimensionalMeta::getIORuleString(rule_wavelet);
            cout << " order: " << orders[i%3];

            if (tests[i%3] == type_integration){
                cout << setw(25) << "integration test";
            }else if (tests[i%3] == type_nodal_interpolation){
                cout << setw(25) << "w-interpolation";
            }else{
                cout << setw(25) << "interpolation";
            }

            cout << "   failed function: " << f->getDescription();
            cout << setw(10) << "observed: " << R.error << "  expected: " << tols[i] << endl;
        }
    }
    delete[] x;
    return bPass;
}
bool ExternalTester::testAllWavelet() const{
    bool pass = true;
    const int depths1[6] = { 7, 7, 7, 5, 5, 5 };
    const double tols1[6] = { 5.E-05, 1.E-04, 1.E-04, 1.E-08, 1.E-07, 1.E-07 };
    int wfirst = 11, wsecond = 34, wthird = 15;
    if (testLocalWaveletRule(&f21nx2, depths1, tols1)){
        //cout << setw(wfirst) << "Rules" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(rule_wavelet) << setw(wthird) << "Pass" << endl;
        cout << setw(wfirst) << "Rules" << setw(wsecond) << "wavelet" << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(rule_wavelet) << setw(wthird) << "FAIL" << endl; pass = false;
    }{ TasGrid::TasmanianSparseGrid grid;
        grid.makeWaveletGrid(2, 1, 2, 1);
        int *indx = 0, *pntr = 0;
        double *vals = 0;
        double *pnts = new double[20]; setRandomX(20, pnts);
        grid.evaluateSparseHierarchicalFunctions(pnts, 10, pntr, indx, vals);
        getError(&f21nx2, &grid, type_internal_interpolation); // this is done to load the values
        const double *coeff = grid.getHierarchicalCoefficients();
        double *y = new double[10];
        grid.evaluateBatch(pnts, 10, y);
        for(int i=0; i<10; i++){
            for(int j=pntr[i]; j<pntr[i+1]; j++){
                y[i] -= coeff[indx[j]] * vals[j];
            }
        }
        for(int i=0; i<10; i++){
            if (fabs(y[i]) > TSG_NUM_TOL){
                cout << "Error in evaluateSparseHierarchicalFunctions() (wavelet)" << endl;
                cout << y[i] << endl;
                pass = false;
            }
        }
        double *v = new double[10 * grid.getNumPoints()];
        getError(&f21nx2, &grid, type_internal_interpolation);
        grid.evaluateHierarchicalFunctions(pnts, 10, v);
        coeff = grid.getHierarchicalCoefficients();
        grid.evaluateBatch(pnts, 10, y);
        for(int i=0; i<10; i++){
            for(int j=0; j<grid.getNumPoints(); j++){
                y[i] -= coeff[j] * v[i*grid.getNumPoints() + j];
            }
        }
        for(int i=0; i<10; i++){
            if (fabs(y[i]) > TSG_NUM_TOL){
                cout << "Error in getHierarchicalCoefficients() (wavelet)" << endl;
                cout << y[i] << endl;
                pass = false;
            }
        }
        delete[] indx;
        delete[] pntr;
        delete[] vals;
        delete[] pnts;
        delete[] y;
        delete[] v;
    }
    return pass;
}

bool ExternalTester::testAllFourier() const{
    bool pass = true;
    const int depths1[3] = { 5, 5, 5 };
    const int depths2[3] = { 4, 4, 4 };
    const double tols1[3] = { 1.E-11, 1.E-06, 1.E-06 };
    const double tols2[3] = { 1.E-11, 1.E-03, 1.E-03 };
    int wfirst = 11, wsecond = 34, wthird = 15;
    if (testGlobalRule(&f21expsincos, TasGrid::rule_fourier, 0, 0, 0, true, depths1, tols1) && testGlobalRule(&f21expsincos, TasGrid::rule_fourier, 0, 0, 0, true, depths2, tols2)){
        cout << setw(wfirst) << "Rules" << setw(wsecond) << "fourier" << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Rules" << setw(wsecond) << "fourier" << setw(wthird) << "FAIL" << endl; pass = false;
    }{ TasGrid::TasmanianSparseGrid grid;
        grid.makeFourierGrid(2, 1, 4, TasGrid::type_level);
        int num_eval = 10;
        int *indx = 0, *pntr = 0;
        double *vals = 0;
        double *pnts = new double[2*num_eval]; setRandomX(2*num_eval, pnts);
        for(int i=0; i<2*num_eval; i++) pnts[i] = 0.5*(pnts[i]+1.0);    // map to [0,1]^d canonical Fourier domain
        grid.evaluateSparseHierarchicalFunctions(pnts, num_eval, pntr, indx, vals);
        getError(&f21expsincos, &grid, type_internal_interpolation); // this is done to load the values
        const double *coeff = grid.getHierarchicalCoefficients();    // coeff = [fourier_coeff_1.real(), fourier_coeff_1.imag(), fourier_coeff_2.real(), ...]
        double *y = new double[num_eval];
        grid.evaluateBatch(pnts, num_eval, y);

        int num_points = grid.getNumPoints();
        int num_outputs = grid.getNumOutputs();
        for(int i=0; i<num_eval; i++){
            for(int j=pntr[i]; j<pntr[i+1]; j++){
                y[i] -= (coeff[indx[j]] * vals[2*j] - coeff[indx[j] + num_points*num_outputs]*vals[2*j+1]);
            }
        }
        for(int i=0; i<num_eval; i++){
            if (fabs(y[i]) > TSG_NUM_TOL){
                cout << "Error in evaluateSparseHierarchicalFunctions() (fourier)" << endl;
                cout << "y["<<i<<"] = "<<y[i] << endl;
                pass = false;
            }
        }

        double *v = new double[2 * num_eval * num_points];
        getError(&f21expsincos, &grid, type_internal_interpolation);
        grid.evaluateHierarchicalFunctions(pnts, num_eval, v);
        grid.evaluateBatch(pnts, num_eval, y);
        for(int i=0; i<num_eval; i++){
            for(int j=0; j<grid.getNumPoints(); j++){
                y[i] -= (coeff[j] * v[2*(i*num_points+j)] - coeff[j+num_points] * v[2*(i*num_points+j)+1]);
            }
        }
        for(int i=0; i<num_eval; i++){
            if (fabs(y[i]) > TSG_NUM_TOL){
                cout << "Error in getHierarchicalCoefficients() (fourier)" << endl;
                cout << "y["<<i<<"] = "<<y[i] << endl;
                pass = false;
            }
        }
        delete[] indx;
        delete[] pntr;
        delete[] vals;
        delete[] pnts;
        delete[] y;
        delete[] v;
    }
    return pass;
}

bool ExternalTester::testAllRefinement() const{
    TasmanianSparseGrid grid;
    bool pass = true;
    {
        const BaseFunction *f = &f21nx2;
        grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_leja);
        int np[13] = { 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 118, 130 };
        double err[13] = { 2, 2.E-1, 5.E-1, 2.E-2, 4.E-2, 2.E-3, 4.E-3, 2.E-4, 2.E-4, 2.E-5, 2.E-5, 8.E-7, 8.E-7 };
        if (!testSurplusRefinement(f, &grid, 1.E-6, refine_classic, np, err, 13)){
            cout << "ERROR: failed leja surplus refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21coscos;
        int np[9] = {    21,    24,    30,    39,    49,    60,    72,    79,    85 };
        double err[9] = { 2.E-1, 7.E-3, 2.E-2, 3.E-4, 6.E-4, 4.E-6, 9.E-6, 5.E-7, 5.E-7 };
        grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), 5, type_iptotal, rule_rleja);
        if (!testSurplusRefinement(f, &grid, 1.E-6, refine_classic, np, err, 9)){
            cout << "ERROR: failed rleja global surplus refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21coscos;
        int np[9] = {    21,    24,    30,    39,    49,    60,    72,    79,    85 };
        double err[9] = { 2.E-1, 7.E-3, 2.E-2, 3.E-4, 6.E-4, 4.E-6, 9.E-6, 5.E-7, 5.E-7 };
        grid.makeSequenceGrid(f->getNumInputs(), f->getNumOutputs(), 5, type_iptotal, rule_rleja);
        if (!testSurplusRefinement(f, &grid, 1.E-6, refine_classic, np, err, 9)){
            cout << "ERROR: failed rleja sequence surplus refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_semilocalp);
        int np[8] = { 29, 65, 145, 321, 705, 1521, 2753, 3569 };
        double err[8] = { 4.E-2, 1.E-2, 1.E-3, 2.E-4, 4.E-5, 5.E-6, 1.E-6, 5.E-7 };
        if (!testSurplusRefinement(f, &grid, 1.E-6, refine_classic, np, err, 8)){
            cout << "ERROR: failed semi-local classic refinement for " << f->getDescription() << endl;  pass = false;
        }
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_semilocalp);
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_semilocalp);
        int np[8] = { 29, 65, 145, 321, 705, 1521, 2753, 3569 };
        double err[8] = { 4.E-2, 1.E-2, 1.E-3, 2.E-4, 4.E-5, 5.E-6, 1.E-6, 5.E-7 };
        if (!testSurplusRefinement(f, &grid, 1.E-6, refine_parents_first, np, err, 8)){
            cout << "ERROR: failed semi-local parents refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 2, 2, rule_semilocalp);
        int np[6] = { 13, 29, 65, 145, 321, 545 };
        double err[6] = { 8.E-02, 5.E-02, 7.E-03, 2.E-03, 3.E-04, 6.E-05 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_direction_selective, np, err, 6)){
            cout << "ERROR: failed semi-local direction refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 2, 2, rule_semilocalp);
        int np[6] = { 13, 29, 65, 145, 321, 545 };
        double err[6] = { 8.E-02, 5.E-02, 7.E-03, 2.E-03, 3.E-04, 6.E-05 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_fds, np, err, 6)){
            cout << "ERROR: failed semi-local fds refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 3, 1, rule_localp);
        int np[10] = { 29, 65, 145, 321, 705, 1537, 3321, 6981, 13517, 19113 };
        double err[10] = { 4.E-2, 2.E-2, 6.E-3, 2.E-3, 6.E-4, 2.E-4, 6.E-5, 2.E-5, 6.E-6, 2.E-6 };
        if (!testSurplusRefinement(f, &grid, 1.E-6, refine_classic, np, err, 10)){
            cout << "ERROR: failed localp classic refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 3, 1, rule_localp);
        int np[10] = { 29, 65, 145, 321, 705, 1537, 3321, 6981, 13517, 19113 };
        double err[10] = { 4.E-2, 2.E-2, 6.E-3, 2.E-3, 6.E-4, 2.E-4, 6.E-5, 2.E-5, 6.E-6, 2.E-6 };
        if (!testSurplusRefinement(f, &grid, 1.E-6, refine_parents_first, np, err, 10)){
            cout << "ERROR: failed localp parents refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 2, 1, rule_localp);
        int np[8] = { 13, 29, 65, 145, 321, 673, 1233, 1433 };
        double err[8] = { 1.E-01, 5.E-02, 2.E-02, 5.E-03, 2.E-03, 6.E-04, 2.E-04, 1.E-04 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_direction_selective, np, err, 8)){
            cout << "ERROR: failed localp direction refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 2, 1, rule_localp);
        int np[8] = { 13, 29, 65, 145, 321, 673, 1233, 1433 };
        double err[8] = { 1.E-01, 5.E-02, 2.E-02, 5.E-03, 2.E-03, 6.E-04, 2.E-04, 1.E-04 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_fds, np, err, 8)){
            cout << "ERROR: failed localp fds refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21sincosaxis;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 3, 1, rule_localpb);
        int np[7] = { 37, 77, 157, 317, 637, 1277, 2317 };
        double err[7] = { 3.E-01, 5.E-02, 2.E-02, 4.E-03, 7.E-04, 3.E-04, 1.E-04 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_classic, np, err, 7)){
            cout << "ERROR: failed localp-boundary fds refinement for " << f->getDescription() << endl;
        }
    }{
        const BaseFunction *f = &f21coscos;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_localp0);
        int np[6] = { 49, 129, 321, 769, 1761, 2209 };
        double err[6] = { 2.E-3, 3.E-4, 5.E-5, 7.E-6, 8.E-7, 5.E-7 };
        if (!testSurplusRefinement(f, &grid, 1.E-6, refine_classic, np, err, 6)){
            cout << "ERROR: failed localp-zero classic refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21coscos;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_localp0);
        int np[6] = { 49, 129, 321, 769, 1761, 2209 };
        double err[6] = { 2.E-3, 3.E-4, 5.E-5, 7.E-6, 8.E-7, 5.E-7 };
        if (!testSurplusRefinement(f, &grid, 1.E-6, refine_parents_first, np, err, 6)){
            cout << "ERROR: failed localp-zero parents refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21coscos;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 2, 2, rule_localp0);
        int np[4] = { 17, 49, 129, 305 };
        double err[4] = { 7.E-03, 2.E-03, 4.E-04, 4.E-05 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_direction_selective, np, err, 4)){
            cout << "ERROR: failed localp-zero direction refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21coscos;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 2, 2, rule_localp0);
        int np[4] = { 17, 49, 129, 305 };
        double err[4] = { 7.E-03, 2.E-03, 4.E-04, 4.E-05 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_fds, np, err, 4)){
            cout << "ERROR: failed localp-zero fds refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 2, 0, rule_localp);
        int np[5] =     {    21,    81,   297, 1053,  3637 };
        double err[5] = { 3.E-1, 2.E-1, 6.E-2, 3E-2, 8.5E-3 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_fds, np, err, 5)){
            cout << "ERROR: failed pwc fds refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21coscos;
        grid.makeWaveletGrid(f->getNumInputs(), f->getNumOutputs(), 2, 1);
        int np[7] = { 49, 81, 193, 449, 993, 1921, 1937 };
        double err[7] = { 6.E-02, 3.E-02, 6.E-03, 3.E-03, 6.E-04, 3.E-04, 2.E-04 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_classic, np, err, 7)){
            cout << "ERROR: failed wavelet classic refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21coscos;
        grid.makeWaveletGrid(f->getNumInputs(), f->getNumOutputs(), 2, 1);
        int np[7] = { 49, 81, 193, 449, 993, 1921, 1937 };
        double err[7] = { 6.E-02, 3.E-02, 6.E-03, 3.E-03, 6.E-04, 3.E-04, 2.E-04 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_parents_first, np, err, 7)){
            cout << "ERROR: failed wavelet parents refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21nx2;
        grid.makeWaveletGrid(f->getNumInputs(), f->getNumOutputs(), 2, 1);
        int np[6] = { 49, 113, 257, 561, 1113, 1481 };
        double err[6] = { 6.E-02, 1.E-02, 5.E-03, 1.E-03, 5.E-04, 1.E-04 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_direction_selective, np, err, 6)){
            cout << "ERROR: failed wavelet direction refinement for " << f->getDescription() << endl;  pass = false;
        }
    }{
        const BaseFunction *f = &f21coscos;
        grid.makeWaveletGrid(f->getNumInputs(), f->getNumOutputs(), 2, 1);
        int np[7] = { 49, 81, 161, 385, 889, 1737, 1769 };
        double err[7] = { 6.E-02, 3.E-02, 6.E-03, 3.E-03, 6.E-04, 3.E-04, 2.E-04 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_fds, np, err, 7)){
            cout << "ERROR: failed wavelet fds refinement for " << f->getDescription() << endl;  pass = false;
        }
    }

    cout << "      Refinement                      surplus" << setw(15) << ((pass) ? "Pass" : "FAIL") << endl;

    bool pass2 = true;
    {
        const BaseFunction *f = &f21aniso;
        grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_leja);
        int np[40] = { 10, 15, 21, 28, 29, 30, 31, 32, 34, 35, 37, 40, 41, 45, 49, 54, 59, 64, 70, 77, 84, 92, 100, 108, 117, 126, 135, 145, 155, 165, 176, 187, 198, 210, 212, 224, 237, 248, 251, 263 };
        double errs[40] = { 9.04e-01, 4.24e-01, 5.73e-01, 2.78e-01, 3.15e-01, 2.49e-01, 3.00e-01, 8.85e-02, 9.30e-02, 9.67e-02, 2.06e-01, 3.03e-01, 5.24e-02, 4.63e-02, 5.85e-02, 5.11e-02, 9.80e-03, 2.71e-02, 5.42e-03, 7.85e-03, 6.21e-03, 5.41e-03, 2.56e-03, 3.32e-03, 5.18e-04, 6.14e-04, 3.66e-04, 4.87e-04, 8.19e-05, 2.58e-04, 5.76e-05, 5.54e-05, 5.22e-05, 4.89e-05, 4.68e-05, 8.92e-06, 2.20e-05, 5.56e-06, 5.14e-06, 5.79e-06 };
        if (!testAnisotropicRefinement(f, &grid, type_iptotal, 1, np, errs, 40)){
            cout << "ERROR: failed anisotropic refinement using leja iptotal nodes for " << f->getDescription() << endl;  pass2 = false;
        }
    }{
        const BaseFunction *f = &f21curved;
        grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_leja);
        int np[10] = { 10, 12, 17, 24, 32, 34, 41, 42, 57, 59 };
        double errs[10] = { 9.48e-03, 9.50e-03, 6.85e-03, 5.11e-04, 6.26e-05, 7.11e-06, 5.07e-06, 5.19e-06, 1.17e-08, 1.86e-08 };
        if (!testAnisotropicRefinement(f, &grid, type_ipcurved, 1, np, errs, 7)){
            cout << "ERROR: failed anisotropic refinement using leja ipcurved nodes for " << f->getDescription() << endl;  pass2 = false;
        }
    }{
        const BaseFunction *f = &f21curved;
        grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_clenshawcurtis);
        int np[3] = { 13, 21, 29 };
        double errs[3] = { 6.12e-04, 6.05e-04, 1.33e-08 };
        if (!testAnisotropicRefinement(f, &grid, type_ipcurved, 1, np, errs, 3)){
            cout << "ERROR: failed anisotropic refinement using clenshaw-curtis ipcurved nodes for " << f->getDescription() << endl;  pass2 = false;
        }
    }{
        const BaseFunction *f = &f21curved;
        grid.makeSequenceGrid(f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_leja);
        int np[10] = { 10, 12, 17, 24, 32, 34, 41, 42, 57, 59 };
        double errs[10] = { 9.48e-03, 9.50e-03, 6.85e-03, 5.11e-04, 6.26e-05, 7.11e-06, 5.07e-06, 5.19e-06, 1.17e-08, 1.86e-08 };
        if (!testAnisotropicRefinement(f, &grid, type_ipcurved, 1, np, errs, 7)){
            cout << "ERROR: failed anisotropic refinement using leja ipcurved nodes for " << f->getDescription() << endl;  pass2 = false;
        }
    }

    cout << "      Refinement                  anisotropic" << setw(15) << ((pass2) ? "Pass" : "FAIL") << endl;

    return (pass && pass2);
}

bool ExternalTester::testAllDomain() const{
    TasmanianSparseGrid grid;
    bool pass1 = true;

    cout << std::scientific; cout.precision(16);
    {
        const BaseFunction *f = &f21nx2aniso;
        int np[5] = {1, 3, 7, 15, 29};
        double errs[5] = {9.E-1, 9.E-1, 2.E-1, 2.E-1, 2.E-2 };
        int aiso_weights[2] = { 2, 1 };
        for(int i=0; i<5; i++){
            grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), i, TasGrid::type_level, TasGrid::rule_clenshawcurtis, aiso_weights);
            TestResults R = getError(f, &grid, type_internal_interpolation);
            if ((R.num_points != np[i]) || (R.error>errs[i])){
                cout << "Using clenshaw-curtis rule" << endl;
                cout << "Failed anisotropic grid test for " << f->getDescription() << "  number of points = " << R.num_points << "  expacted: " << np[i]
                     << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass1 = false;
            }
        }
    }{
        const BaseFunction *f = &f21nx2aniso;
        int np[5] = {5, 10, 15, 28, 37};
        double errs[5] = {3.E-1, 3.E-1, 5.E-2, 2.E-2, 2.E-3};
        int aiso_weights[2] = {2, 1};
        for(int i=0; i<5; i++){
            grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), i+2, TasGrid::type_level, TasGrid::rule_gausslegendre, aiso_weights);
            TestResults R = getError(f, &grid, type_integration);
            if ((R.num_points != np[i]) || (R.error>errs[i])){
                cout << "Using gauss-legendre rule" << endl;
                cout << "Failed anisotropic grid test for " << f->getDescription() << "  number of points = " << R.num_points << "  expacted: " << np[i]
                     << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass1 = false;
            }
        }
    }{
        const BaseFunction *f = &f21nx2aniso;
        int np[5] = {36, 42, 49, 56, 64};
        double errs[5] = {5.E-2, 5.E-2, 5.E-3, 5.E-3, 3.E-3 };
        int aiso_weights[2] = { 2, 1 };
        for(int i=0; i<5; i++){
            grid.makeSequenceGrid(f->getNumInputs(), f->getNumOutputs(), i+10, TasGrid::type_level, TasGrid::rule_leja, aiso_weights);
            TestResults R = getError(f, &grid, type_internal_interpolation);
            if ((R.num_points != np[i]) || (R.error>errs[i])){
                cout << "Using leja rule" << endl;
                cout << "Failed anisotropic grid test for " << f->getDescription() << "  number of points = " << R.num_points << "  expacted: " << np[i]
                     << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass1 = false;
            }
        }
    }

    cout << "      Domain                      anisotropic" << setw(15) << ((pass1) ? "Pass" : "FAIL") << endl;

    bool pass2 = true;

    {
        const BaseFunction *f = &f21expDomain;
        double errs[5] = {3.E-5, 2.E-7, 2.E-8, 6.E-10, 4.E-11 };
        double errs2[5] = {1.E-5, 2.E-6, 6.E-8, 6.E-9, 8.E-11 };
        double errs3[5] = {6.E-2, 1.E-2, 5.E-3, 6.E-4, 6.E-5 };
        double transform_a[2] = { 3.0, -3.0 };
        double transform_b[2] = { 4.0,  2.0 };
        for(int i=0; i<5; i++){
            grid.makeSequenceGrid(f->getNumInputs(), f->getNumOutputs(), i+5, TasGrid::type_level, TasGrid::rule_leja);
            grid.setDomainTransform(transform_a, transform_b);
            double *needed_points = grid.getNeededPoints();
            int num_needed = grid.getNumNeeded();
            TestResults R = getError(f, &grid, type_integration);
            if (R.error > errs[i]){
                cout << "Using leja rule" << endl;
                cout << "Failed domain transform test for " << f->getDescription() << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass2 = false;
            }
            double test_x[2] = {3.314, -1.71732};
            R = getError(f, &grid, type_nodal_interpolation, test_x);
            if (R.error > errs2[i]){
                cout << "Using leja rule" << endl;
                cout << "Failed domain transform test interpolation for " << f->getDescription() << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass2 = false;
            }
            R = getError(f, &grid, type_internal_interpolation, test_x);
            if (R.error > errs3[i]){
                cout << "Using leja rule" << endl;
                cout << "Failed domain transform test interpolation for " << f->getDescription() << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass2 = false;
            }
            double *loaded_points = grid.getLoadedPoints();
            for(int j=0; j<num_needed * f->getNumInputs(); j++){
                if (fabs(needed_points[j] - loaded_points[j]) > TSG_NUM_TOL){
                    cout << "Mismatch between needed and loaded points" << endl;
                    pass2 = false;
                }
            }
            delete[] loaded_points;
            delete[] needed_points;
        }
    }{
        const BaseFunction *f = &f21expDomain;
        double errs[5] = {7.E-1, 3.E-3, 6.E-3, 4.E-5, 4.E-6 };
        double transform_a[2] = { 3.0, -3.0 };
        double transform_b[2] = { 4.0,  2.0 };
        for(int i=0; i<5; i++){
            grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), i, TasGrid::type_iptotal, TasGrid::rule_clenshawcurtis);
            grid.setDomainTransform(transform_a, transform_b);
            TestResults R = getError(f, &grid, type_integration);
            if (R.error>errs[i]){
                cout << "Using clenshaw-curtis rule" << endl;
                cout << "Failed domain transform test for " << f->getDescription() << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass2 = false;
            }
        }
    }{
        const BaseFunction *f = &f21expDomain;
        double errs[5] = {7.E-1, 7.E-3, 3.E-5, 6.E-8, 7.E-11 };
        double transform_a[2] = { 3.0, -3.0 };
        double transform_b[2] = { 4.0,  2.0 };
        for(int i=0; i<5; i++){
            grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), i, TasGrid::type_iptotal, TasGrid::rule_gausslegendre);
            grid.setDomainTransform(transform_a, transform_b);
            TestResults R = getError(f, &grid, type_integration);
            if (R.error>errs[i]){
                cout << "Using gauss-legendre rule" << endl;
                cout << "Failed domain transform test for " << f->getDescription() << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass2 = false;
            }
        }
    }{
        const BaseFunction *f = &f21expDomain;
        double errs[5] = {7.E-1, 4.E-1, 4.E-3, 8.E-4, 7.E-7 };
        double transform_a[2] = { 3.0, -3.0 };
        double transform_b[2] = { 4.0,  2.0 };
        for(int i=0; i<5; i++){
            grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), i, 2, TasGrid::rule_localp);
            grid.setDomainTransform(transform_a, transform_b);
            TestResults R = getError(f, &grid, type_integration);
            if (R.error>errs[i]){
                cout << "Using localp rule" << endl;
                cout << "Failed domain transform test for " << f->getDescription() << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass2 = false;
            }
        }
    }{
        const BaseFunction *f = &f21coscos;
        double errs[5] = { 5.E-2, 7.E-3, 9.E-4, 1.E-13 };
        double transform_a[2] = { -1.0, -1.0 };     // canonical domain of Fourier grid is [0,1]
        double transform_b[2] = { 1.0,  1.0 };
        for(int i=0; i<3; i++){     // keeping depth low until we implement an FFT algorithm
            grid.makeFourierGrid(f->getNumInputs(), f->getNumOutputs(), i+3, TasGrid::type_level);
            grid.setDomainTransform(transform_a, transform_b);
            TestResults R = getError(f, &grid, type_integration);
            if (R.error>errs[i]){
                cout << "Using fourier rule" << endl;
                cout << "Failed domain transform test for " << f->getDescription() << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass2 = false;
            }
        }
    }

    cout << "      Domain                      transformed" << setw(15) << ((pass2) ? "Pass" : "FAIL") << endl;

    bool pass3 = true;

    {
        TasmanianSparseGrid gridc;
        const BaseFunction *f = &f21conformal;
        int asin_conformal[2] = {4, 4};
        for(int l=0; l<6; l++){
            grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), l+2, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
            gridc.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), l+2, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
            gridc.setConformalTransformASIN(asin_conformal);
            double *needed_points = grid.getNeededPoints();
            int num_needed = grid.getNumNeeded();
            TestResults R1 = getError(f, &grid, type_internal_interpolation);
            TestResults R2 = getError(f, &gridc, type_internal_interpolation);
            if (R1.num_points != R2.num_points){
                cout << "Failed in number of points for conformal mapping and clenshaw-curtis rule" << endl;
                pass3 = false;
            }
            if (R1.error < R2.error){
                cout << "Failed in error for conformal mapping and clenshaw-curtis rule" << endl;
                cout << "  standard error = " << R1.error << endl;
                cout << " conformal error = " << R2.error << endl;
                pass3 = false;
            }
            double y1, y2, y_true;
            grid.integrate(&y1);
            gridc.integrate(&y2);
            f->getIntegral(&y_true);
            if (fabs(y1 - y_true) < fabs(y2 - y_true)){
                cout << "Failed in error for conformal mapping and clenshaw-curtis rule" << endl;
                cout << "  standard error = " << fabs(y1 - y_true) << endl;
                cout << " conformal error = " << fabs(y2 - y_true) << endl;
                pass3 = false;
            }
            double *loaded_points = grid.getLoadedPoints();
            for(int j=0; j<num_needed * f->getNumInputs(); j++){
                if (fabs(needed_points[j] - loaded_points[j]) > TSG_NUM_TOL){
                    cout << "Mismatch between needed and loaded points" << endl;
                    pass2 = false;
                }
            }
            delete[] loaded_points;
            delete[] needed_points;
        }
    }{
        TasmanianSparseGrid gridc;
        const BaseFunction *f = &f21conformal;
        int asin_conformal[2] = { 4, 4 };
        for(int l=0; l<1; l++){
            grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), l+5, TasGrid::type_iptotal, TasGrid::rule_gausspatterson);
            gridc.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), l+5, TasGrid::type_iptotal, TasGrid::rule_gausspatterson);
            gridc.setConformalTransformASIN(asin_conformal);
            TestResults R1 = getError(f, &grid, type_integration);
            TestResults R2 = getError(f, &gridc, type_integration);
            if (R1.num_points != R2.num_points){
                cout << "Failed in number of points for conformal mapping and gauss-patterson rule" << endl;
                pass3 = false;
            }
            if (R1.error < R2.error){
                cout << "Failed in error for conformal mapping and gauss-patterson rule" << endl;
                cout << "  standard error = " << R1.error << endl;
                cout << " conformal error = " << R2.error << endl;
                pass3 = false;
            }
            R1 = getError(f, &grid, type_integration);
            R2 = getError(f, &gridc, type_integration);
            if (R1.error < R2.error){
                cout << "Failed in error for conformal mapping and gauss-patterson rule (integration)" << endl;
                cout << "  standard error = " << R1.error << endl;
                cout << " conformal error = " << R2.error << endl;
                pass3 = false;
            }
        }
    }{
        TasmanianSparseGrid gridc;
        const BaseFunction *f = &f21conformal;
        int asin_conformal[2] = { 4, 4 };
        for(int l=0; l<5; l++){
            grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), l+4, 3, TasGrid::rule_localp);
            gridc.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), l+4, 3, TasGrid::rule_localp);
            gridc.setConformalTransformASIN(asin_conformal);
            TestResults R1 = getError(f, &grid, type_internal_interpolation);
            TestResults R2 = getError(f, &gridc, type_internal_interpolation);
            if (R1.num_points != R2.num_points){
                cout << "Failed in number of points for conformal mapping and local polynomial rule" << endl;
                pass3 = false;
            }
            if (R1.error < R2.error){
                cout << "Failed in error for conformal mapping and local polynomial rule" << endl;
                cout << "  standard error = " << R1.error << endl;
                cout << " conformal error = " << R2.error << endl;
                pass3 = false;
            }
            R1 = getError(f, &grid, type_integration);
            R2 = getError(f, &gridc, type_integration);
            if (R1.error < R2.error){
                cout << "Failed in error for conformal mapping and local polynomial rule (integration)" << endl;
                cout << "  standard error = " << R1.error << endl;
                cout << " conformal error = " << R2.error << endl;
                pass3 = false;
            }
        }
    }

    cout << "      Domain                        conformal" << setw(15) << ((pass3) ? "Pass" : "FAIL") << endl;

    return (pass1 && pass2 && pass3);
}

bool ExternalTester::testAcceleration(const BaseFunction *f, TasmanianSparseGrid *grid) const{
    int dims = f->getNumInputs();
    int outs = f->getNumOutputs();
    if (grid->getNumNeeded() > 0){
        std::vector<double> points, vals(outs * grid->getNumPoints());
        grid->getNeededPoints(points);

        for(int i=0; i<grid->getNumPoints(); i++){
            f->eval(&(points[i*dims]), &(vals[i*outs]));
        }
        grid->loadNeededPoints(vals);
    }

    int num_x = 256; // for batched evaluations
    int num_fast = 16; // for fast evaluations, must be <= num_x
    std::vector<double> x(num_x * dims);
    setRandomX(num_x * dims, x.data());

    std::vector<double> test_y, baseline_y(outs * num_x);
    for(int i=0; i<num_x; i++) grid->evaluate(&(x[i*dims]), &(baseline_y[i*outs]));

    bool pass = true;
    TypeAcceleration acc[5] = {accel_none, accel_cpu_blas, accel_gpu_cublas, accel_gpu_cuda, accel_gpu_magma};
    int testGpuID = (gpuid == -1) ? 0 : gpuid;
    int c = 0;
    while(c < 5){
        grid->enableAcceleration(acc[c]);
        if (c > 1){ // gpu test
            grid->setGPUID(testGpuID);
        }
        //grid->printStats();
        //cout << "Testing Batch evaluations" << endl;

        test_y.resize(1); // makes sure that evaluate sets the right dimension

        //grid->favorSparseAcceleration(false); // switch between explicit test for the dense and sparse algorithms
        grid->evaluateBatch(x, test_y);

        double err = 0.0;
        for(int i=0; i<outs*num_x; i++) if (fabs(test_y[i] - baseline_y[i]) > err) err = fabs(test_y[i] - baseline_y[i]);

        if (err > 1.E-11){
            int tm = 0; err = 0.0;
            for(int i=0; i<outs*num_x; i++){ if (fabs(test_y[i] - baseline_y[i]) > err){ err = fabs(test_y[i] - baseline_y[i]); tm = i; }}
            cout << tm << "  " << test_y[tm] << "    " << baseline_y[tm] << endl;

            pass = false;
            cout << "Failed Batch evaluation for acceleration c = " << c << " gpuID = " << testGpuID << endl;
            cout << "Observed error: " << err << " for function: " << f->getDescription() << endl;
            grid->printStats();
            exit(1);
        }

        //cout << "Testing Fast evaluations" << endl;
        grid->evaluateFast(x, test_y);
        for(int i= 1; i<num_fast; i++){
            grid->evaluateFast(&(x[i*dims]), &(test_y[i*outs]));
        }

        err = 0.0;
        for(int i=0; i<outs*num_fast; i++) if (fabs(test_y[i] - baseline_y[i]) > err) err = fabs(test_y[i] - baseline_y[i]);

        if (err > 1.E-11){
            pass = false;
            cout << "Failed Fast evaluation for acceleration c = " << c << " gpuID = " << testGpuID << endl;
            cout << "Observed error: " << err << " for function: " << f->getDescription() << endl;
            grid->printStats();
        }

        if (c > 1){ // gpu test
            if (gpuid == -1){ // gpuid is not set, then cycle trough all GPUs
                testGpuID++;
                if (testGpuID >= grid->getNumGPUs()){
                    testGpuID = 0;
                    c++;
                }
            }else{
                c++;
            }
        }else{
            c++;
        }
    }

    //cout << "End of acceleration test." << endl;
    return pass;
}

bool ExternalTester::testGPU2GPUevaluations() const{
    #ifdef Tasmanian_ENABLE_CUDA
    // check back basis evaluations, x and result both sit on the GPU (using CUDA acceleration)
    TasGrid::TasmanianSparseGrid grid;
    int num_tests = 9;
    int dims = 3;
    TasGrid::TypeOneDRule pwp_rule[9] = {TasGrid::rule_localp, TasGrid::rule_localp0, TasGrid::rule_semilocalp, TasGrid::rule_localpb,
                                         TasGrid::rule_localp, TasGrid::rule_localp0, TasGrid::rule_semilocalp, TasGrid::rule_localpb,
                                         TasGrid::rule_localp};
    int order[9] = {1, 1, 1, 1, 2, 2, 2, 2, 0};
    std::vector<double> a = {3.0, 4.0, -10.0}, b = {5.0, 7.0, 2.0};

    bool pass = true;
    int gpu_index_first = (gpuid == -1) ? 0 : gpuid;
    int gpu_end_gpus = (gpuid == -1) ? grid.getNumGPUs() : gpuid+1;
    for(int t=0; t<num_tests; t++){
        grid.makeLocalPolynomialGrid(dims, 1, ((order[t] == 0) ? 4 : 7), order[t], pwp_rule[t]);

        grid.setDomainTransform(a, b);
        grid.enableAcceleration(TasGrid::accel_none);

        int nump = 2000;
        double *x = new double[dims*nump];
        std::vector<double> xt(dims * nump);
        setRandomX(dims*nump, x);
        for(int i=0; i<nump; i++){
            for(int j=0; j<dims; j++){
                xt[dims*i + j] = 0.5 * (b[j] - a[j]) * x[dims*i+j] + 0.5 * (b[j] + a[j]);
            }
        }
        delete[] x;

        // Dense version:
        std::vector<double> y_true_dense;
        grid.evaluateHierarchicalFunctions(xt, y_true_dense);
        // grid.printStats();
        // cout << "Memory requirements = " << (grid.getNumPoints() * nump * 8) / (1024 * 1024) << "MB" << endl;

        int *pntr = 0, *indx = 0;
        double *vals = 0;
        grid.evaluateSparseHierarchicalFunctions(xt.data(), nump, pntr, indx, vals);

        for(int gpuID=gpu_index_first; gpuID < gpu_end_gpus; gpuID++){
            bool dense_pass = true;

            cudaSetDevice(gpuID);
            double *gpux = TasGrid::TasCUDA::cudaSend<double>(xt);
            double *gpuy = TasGrid::TasCUDA::cudaNew<double>(grid.getNumPoints() * nump);

            grid.enableAcceleration(TasGrid::accel_gpu_cuda);
            grid.setGPUID(gpuID);
            grid.evaluateHierarchicalFunctionsGPU(gpux, nump, gpuy);

            double *y = new double[grid.getNumPoints() * nump];
            TasGrid::TasCUDA::cudaRecv<double>(grid.getNumPoints() * nump, gpuy, y);

            for(int i=0; i<grid.getNumPoints() * nump; i++){
                if (fabs(y[i] - y_true_dense[i]) > 1.E-11){
                    dense_pass = false;
                }
            }

            if (!dense_pass){
                cout << "Failed evaluateHierarchicalFunctionsGPU() when using grid: " << endl;
                grid.printStats();
            }
            pass = pass && dense_pass;

            TasGrid::TasCUDA::cudaDel<double>(gpuy);
            delete[] y;

            // Sparse version:
            bool sparse_pass = true;
            int *gpu_indx = 0, *gpu_pntr = 0, num_nz = 0;
            double *gpu_vals = 0;
            grid.enableAcceleration(TasGrid::accel_gpu_cuda);
            grid.setGPUID(gpuID);
            grid.evaluateSparseHierarchicalFunctionsGPU(gpux, nump, gpu_pntr, gpu_indx, gpu_vals, num_nz);

            std::vector<int> cpntr; TasGrid::TasCUDA::cudaRecv<int>(nump+1, gpu_pntr, cpntr);
            std::vector<int> cindx; TasGrid::TasCUDA::cudaRecv<int>(num_nz, gpu_indx, cindx);
            std::vector<double> cvals; TasGrid::TasCUDA::cudaRecv<double>(num_nz, gpu_vals, cvals);

            int num_points = grid.getNumPoints();
            std::vector<double> cpu_dense(nump * num_points, 0.0);
            std::vector<double> gpu_dense(nump * num_points, 0.0);
            for(int i=0; i<nump; i++){
                for(int j=pntr[i]; j < pntr[i+1]; j++){
                    cpu_dense[i * num_points + indx[j]] = vals[j];
                }
                for(int j=cpntr[i]; j < cpntr[i+1]; j++){
                    gpu_dense[i * num_points + cindx[j]] = cvals[j];
                }
            }
            auto igpu = gpu_dense.begin();
            for(auto c : cpu_dense){
                double diff = fabs(c - *igpu++);
                if (diff > 1.E-12){
                    cout << "ERROR: difference in sparse matrix: " << diff << endl;
                    sparse_pass = false;
                }
            }
            // simple test, works on Linux but Windows CUDA seems to round differently (keep in case of performance regression)
            //if (pntr[nump] != num_nz){
            //    cout << "ERROR: mismatch in the numnz from cuda: " << num_nz << " and cpu " << pntr[nump] << endl;
            //    grid.printStats();
            //    sparse_pass = false;
            //}
            //if (sparse_pass){
            //    for(int i=0; i<nump; i++){
            //        for(int j=pntr[i]; j<pntr[i+1]; j++){
            //            if (sparse_pass && (indx[j] != cindx[j])){
            //                cout << "ERROR: mismatch in index i = " << i << "   j = " << j << "  indx[j] = " << indx[j] << "  cindx[j] = " << cindx[j] << endl;
            //                sparse_pass = false;
            //            }
            //            if (sparse_pass && (fabs(vals[i] - cvals[i]) > 1.E-12)){
            //                cout << "ERROR: i = " << i << "  cpu value = " << vals[i] << "  gpu value = " << cvals[i] << endl;
            //                sparse_pass = false;
            //            }
            //        }
            //    }
            //}
            if (!sparse_pass){
                cout << "Failed evaluateSparseHierarchicalFunctionsGPU() when using grid: " << endl;
                grid.printStats();
            }

            pass = pass && sparse_pass;

            TasGrid::TasCUDA::cudaDel<double>(gpux);
            TasGrid::TasCUDA::cudaDel<int>(gpu_pntr);
            TasGrid::TasCUDA::cudaDel<int>(gpu_indx);
            TasGrid::TasCUDA::cudaDel<double>(gpu_vals);
        }

        delete[] pntr;
        delete[] indx;
        delete[] vals;
    }

    // Sequence Grid evaluations of the basis functions
    {int numx = 2020;

    std::vector<double> cpux(numx * dims);
    setRandomX(numx * dims, cpux.data());

    grid.makeSequenceGrid(dims, 0, 20, type_level, rule_rleja);
    //cout << "Memory requirements = " << (grid.getNumPoints() * numx * 8) / (1024 * 1024) << "MB" << endl;
    std::vector<double> truey;
    grid.evaluateHierarchicalFunctions(cpux, truey);

    for(int gpuID=gpu_index_first; gpuID < gpu_end_gpus; gpuID++){
        grid.makeSequenceGrid(dims, 0, 20, type_level, rule_rleja);
        cudaSetDevice(gpuID);
        grid.enableAcceleration(TasGrid::accel_gpu_cuda);
        grid.setGPUID(gpuID);

        double *gpux = TasCUDA::cudaSend(cpux);
        double *gpuy = TasCUDA::cudaNew<double>(numx * grid.getNumPoints());
        grid.evaluateHierarchicalFunctionsGPU(gpux, numx, gpuy);
        double *cpuy = new double[numx * grid.getNumPoints()];
        TasCUDA::cudaRecv<double>(numx * grid.getNumPoints(), gpuy, cpuy);

        double err = 0.0;
        for(int i=0; i<numx * grid.getNumPoints(); i++){
            //cout << cpuy[i] << "  " << truey[i] << "  " << fabs(cpuy[i] - truey[i]) << endl;
            if (err < fabs(cpuy[i] - truey[i])) err = fabs(cpuy[i] - truey[i]);
        }
        if (err > TSG_NUM_TOL){
            cout << "ERROR: failed Sequence grid GPU basis evaluations with error: " << err << endl;
            grid.printStats();
            pass = false;
        }
        delete[] cpuy;
        TasCUDA::cudaDel<double>(gpuy);
        TasCUDA::cudaDel<double>(gpux);
    }}

    // Fourier Grid evaluations of the basis functions
    {int numx = 2020;

    std::vector<double> cpux(numx * dims);
    setRandomX(numx * dims, cpux.data());
    for(int i=0; i<numx; i++){
        for(int j=0; j<dims; j++){
            cpux[dims*i + j] = cpux[dims*i + j] * (b[j] - a[j]) + a[j];
        }
    }

    grid.makeFourierGrid(dims, 0, 5, type_level); //cout << grid.getNumPoints() << endl;
    grid.setDomainTransform(a, b);
    //cout << "Memory requirements = " << (grid.getNumPoints() * numx * 16) / (1024 * 1024) << "MB" << endl;
    std::vector<double> truey;
    grid.evaluateHierarchicalFunctions(cpux, truey);

    for(int gpuID=gpu_index_first; gpuID < gpu_end_gpus; gpuID++){
        grid.makeFourierGrid(dims, 0, 5, type_level);
        grid.setDomainTransform(a, b);
        cudaSetDevice(gpuID);
        grid.enableAcceleration(TasGrid::accel_gpu_cuda);
        grid.setGPUID(gpuID);

        double *gpux = TasCUDA::cudaSend(cpux);
        double *gpuy = TasCUDA::cudaNew<double>(2 * numx * grid.getNumPoints());
        grid.evaluateHierarchicalFunctionsGPU(gpux, numx, gpuy);
        double *cpuy = new double[2 * numx * grid.getNumPoints()];
        TasCUDA::cudaRecv<double>(2 * numx * grid.getNumPoints(), gpuy, cpuy);

        double err = 0.0;
        for(int i=0; i<2 * numx * grid.getNumPoints(); i++){
            //cout << cpuy[i] << "  " << truey[i] << "  " << fabs(cpuy[i] - truey[i]) << endl;
            if (err < fabs(cpuy[i] - truey[i])) err = fabs(cpuy[i] - truey[i]);
        }
        if (err > TSG_NUM_TOL){
            cout << "ERROR: failed Sequence grid GPU basis evaluations with error: " << err << endl;
            grid.printStats();
            pass = false;
        }
        delete[] cpuy;
        TasCUDA::cudaDel<double>(gpuy);
        TasCUDA::cudaDel<double>(gpux);
    }}

    return pass;

    #else
    return true;
    #endif // Tasmanian_ENABLE_CUDA
}

bool ExternalTester::testAllAcceleration() const{
    const BaseFunction *f = &f23Kexpsincos;
    TasmanianSparseGrid grid;
    bool pass = true;

    int wsecond = 28, wthird = 15;

    grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), 5, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
    pass = pass && testAcceleration(f, &grid);
    if (pass){
        if (verbose) cout << "      Accelerated" << setw(wsecond) << "global" << setw(wthird) << "Pass" << endl;
    }else{
        cout << "      Accelerated" << setw(wsecond) << "global" << setw(wthird) << "FAIL" << endl;
    }

    grid.makeSequenceGrid(f->getNumInputs(), f->getNumOutputs(), 5, TasGrid::type_level, TasGrid::rule_leja);
    pass = pass && testAcceleration(f, &grid);
    if (pass){
        if (verbose) cout << "      Accelerated" << setw(wsecond) << "sequence" << setw(wthird) << "Pass" << endl;
    }else{
        cout << "      Accelerated" << setw(wsecond) << "sequence" << setw(wthird) << "FAIL" << endl;
    }

    // for the purpose of testing CUDA evaluations, test all three localp rules vs orders 0, 1, and 2
    // for order 0, regardless of the selected rule, thegrid should switch to localp
    TasGrid::TypeOneDRule pwp_rule[4] = {TasGrid::rule_localp, TasGrid::rule_localp0, TasGrid::rule_semilocalp, TasGrid::rule_localpb};
    for(int t=0; t<12; t++){
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), ((t / 4 == 0) ? 5 : 6), (t / 4), pwp_rule[t % 4]);
        pass = pass && testAcceleration(f, &grid);
    }
    // test cusparse sparse mat times dense vec used in accel_type cuda
    grid.makeLocalPolynomialGrid(f21nx2.getNumInputs(), f21nx2.getNumOutputs(), 5, 1, TasGrid::rule_localp);
    pass = pass && testAcceleration(&f21nx2, &grid);
    if (pass){
        if (verbose) cout << "      Accelerated" << setw(wsecond) << "local polynomial" << setw(wthird) << "Pass" << endl;
    }else{
        cout << "      Accelerated" << setw(wsecond) << "local polynomial" << setw(wthird) << "FAIL" << endl;
    }

    grid.makeFourierGrid(f->getNumInputs(), f->getNumOutputs(), 4, TasGrid::type_level);
    pass = pass && testAcceleration(f, &grid);
    if (pass){
        if (verbose) cout << "      Accelerated" << setw(wsecond) << "fourier" << setw(wthird) << "Pass" << endl;
    }else{
        cout << "      Accelerated" << setw(wsecond) << "fourier" << setw(wthird) << "FAIL" << endl;
    }

    #ifdef Tasmanian_ENABLE_CUDA
    pass = pass && testGPU2GPUevaluations();
    if (pass){
        if (verbose) cout << "      Accelerated" << setw(wsecond) << "gpu-to-gpu" << setw(wthird) << "Pass" << endl;
    }else{
        cout << "      Accelerated" << setw(wsecond) << "gpu-to-gpu" << setw(wthird) << "FAIL" << endl;
    }
    #else
    if (verbose) cout << "      Accelerated" << setw(wsecond) << "gpu-to-gpu" << setw(wthird) << "Skipped (needs Tasmanian_ENABLE_CUDA=ON)" << endl;
    #endif // Tasmanian_ENABLE_CUDA

    cout << "      Acceleration                        all" << setw(15) << ((pass) ? "Pass" : "FAIL") << endl;

    return pass;
}

#ifndef _TASMANIAN_WINDOWS_
#include <sys/time.h>
extern "C" double gettime(){
    struct timeval t;
    gettimeofday(&t, NULL);
    return ((double) t.tv_sec) + ((double) t.tv_usec) * 1.0E-6;
}
#else
    double gettime(){ return ((double) time(0)); }
#endif // _TASMANIAN_WINDOWS_
void loadGridValues(TasmanianSparseGrid *grid){ // for benchmark
    int dims = grid->getNumDimensions(), num_outputs = grid->getNumOutputs();
    int num_points = grid->getNumNeeded();
    if (num_points == 0) return;
    double *x = grid->getNeededPoints();
    double *v = new double[num_outputs * num_points];
    double s = 1.0 / ((double) (dims*dims));
    for(int i=0; i<num_points; i++){
        double nx2 = 0.0;
        for(int k=0; k<dims; k++) nx2 += x[i*dims +k] * x[i*dims +k];
        nx2 *= s;
        for(int j=0; j<num_outputs-1; j+=2){
            v[i*num_outputs + j]     = exp(-nx2);
            v[i*num_outputs + j + 1] = sin(nx2);
        }
        if (num_outputs % 2 == 1){
            v[i*num_outputs + num_outputs - 1] = exp(-nx2);
        }
    }
    //cout << "Computing surplusses = " << num_outputs << endl;
    //double start = gettime();
    grid->loadNeededPoints(v);
    //cout << gettime() - start << endl;
}

void ExternalTester::benchmark(int argc, const char **argv){
    if ((argc < 3) || (strcmp(argv[2],"help") == 0)){
        cout << "Accepted benchmarks:" << endl;
        cout << "./tasgrid -bench alpha <dims> <outs> <depth> <type> <rule> <batch size> <iterations> <gpu> <use fast>" << endl;
        return;
    }
    if (strcmp(argv[2],"alpha") == 0){
        if (argc > 11){
            cout << "./tasgrid -bench alpha <dims> <outs> <depth> <type> <rule> <batch size> <iterations> <gpu> fast-eval" << endl;
        }else{
            cout << "./tasgrid -bench alpha <dims> <outs> <depth> <type> <rule> <batch size> <iterations> <gpu> batch-eval" << endl;
        }
        cout << "Compare acceleration: if batch size is 0 use surpluses, otherwise use evaluateBatch()";
        if (argc < 11){ cout << endl; return; }
        int dims = atoi(argv[3]), outs = atoi(argv[4]), depth = atoi(argv[5]);
        TypeDepth d = OneDimensionalMeta::getIOTypeString(argv[6]);
        TypeOneDRule r = OneDimensionalMeta::getIORuleString(argv[7]);
        int num_x = atoi(argv[8]), num_runs = atoi(argv[9]), gpu = atoi(argv[10]);
        TasmanianSparseGrid *grid = new TasmanianSparseGrid();
        if (OneDimensionalMeta::isLocalPolynomial(r)){
            grid->makeLocalPolynomialGrid(dims, outs, depth, 2, r);
        }else if (OneDimensionalMeta::isSequence(r)){
            grid->makeSequenceGrid(dims, outs, depth, d, r);
        }else if (OneDimensionalMeta::isFourier(r)){
            grid->makeFourierGrid(dims, outs, depth, d);
        }else{
            grid->makeGlobalGrid(dims, outs, depth, d, r);
        }

        int np = grid->getNumPoints();
        cout << " with " << np << " points." << endl;
        int width = 15;
        cout << setw(24) << "CPU";
        if (gpu > -1){
            char *name = grid->getGPUName(gpu);
            cout << setw(2*width + width/2) << name;
            delete[] name;
        }else{
            cout << setw(width + width/2) << "CPU";
        }
        cout << endl;
        cout << setw(width) << "num outputs";
        TypeAcceleration cpu_tests[2] = {accel_none, accel_cpu_blas};
        TypeAcceleration gpu_tests[4] = {accel_cpu_blas, accel_gpu_cublas, accel_gpu_cuda, accel_gpu_magma};
        TypeAcceleration *tests;
        int num_tests;
        if (gpu > -1){
            num_tests = 4;
            tests = gpu_tests;
            cout << setw(width) << "cpu_blas" << setw(width) << "gpu_cublas" << setw(width) << "gpu_cuda" << setw(width) << "gpu_magma" << endl;
        }else{
            num_tests = 2;
            tests = cpu_tests;
            cout << setw(width) << "none" << setw(width) << "cpu_blas" << endl;
        }
        if (num_x > 0){
            double *x = new double[dims * num_x];
            cout << std::scientific;
            cout.precision(5);
            setRandomX(dims * num_x, x);
            for(int s=0; s<7; s++){
                if (OneDimensionalMeta::isLocalPolynomial(r)){
                    grid->makeLocalPolynomialGrid(dims, outs, depth, 2, r);
                    //grid->favorSparseAlgorithmForLocalPolynomials(true); // use sparse or dense algorithm, defaults to auto
                    //grid->favorSparseAlgorithmForLocalPolynomials(false);
                }else if (OneDimensionalMeta::isSequence(r)){
                    grid->makeSequenceGrid(dims, outs, depth, d, r);
                }else if (OneDimensionalMeta::isFourier(r)){
                    grid->makeFourierGrid(dims, outs, depth, d);
                }else{
                    grid->makeGlobalGrid(dims, outs, depth, d, r);
                }
                loadGridValues(grid);
                double *y = new double[outs * num_x];
                cout << setw(width) << outs;
                double start;

                for(int t=0; t<num_tests; t++){
                    grid->enableAcceleration(tests[t]);
                    if ((gpu > -1) && (t > 0)) grid->setGPUID(gpu);
                    grid->evaluateFast(x, y);
                    start = gettime();
                    //if ((num_tests == 2) || (t > -1))
                    if (argc == 11){
                        for(int i=0; i<num_runs; i++){
                            grid->evaluateBatch(x, num_x, y);
                        }
                    }else{
                        for(int i=0; i<num_runs; i++){
                            for(int j=0; j<num_x; j++){
                                grid->evaluateFast(&(x[j*dims]), &(y[j*outs]));
                            }
                        }
                    }
                    cout << setw(width) << ((int) ((gettime() - start) * 1000.0));
                }
                outs *= 2;
                cout << setw(width+5) << "milliseconds" << endl;
                delete[] y;
            }
            delete[] x;
        }
    }
}

void ExternalTester::debugTest(){
    cout << "Debug Test" << endl;
    cout << "Put here testing code and call this with ./tasgrid -test debug" << endl;
}

void ExternalTester::debugTestII(){
    cout << "Debug Test II" << endl;
    cout << "Put here testing code and call this with ./tasgrid -test db" << endl;
}

#endif
