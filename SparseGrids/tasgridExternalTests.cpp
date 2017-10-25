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

ExternalTester::ExternalTester(int in_num_mc) : num_mc(in_num_mc), verbose(false){ /*srand(time(0));*/ srand(10); }
ExternalTester::~ExternalTester(){}
void ExternalTester::resetRandomSeed(){ srand(time(0)); }

void ExternalTester::setVerbose(bool new_verbose){ verbose = new_verbose; }

void ExternalTester::setRandomX(int n, double x[]) const{
    for(int i=0; i<n; i++){
        x[i] = 2.0 * ((double) rand()) / ((double) RAND_MAX) -1.0;
    }
}

bool ExternalTester::Test() const{
    cout << endl << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "       Tasmanian Sparse Grids Module: Functionality Test" << endl;
    cout << "---------------------------------------------------------------------" << endl << endl;

    bool passAccel   = testAllAcceleration();
    bool passDomain  = testAllDomain();
    bool passRefine  = testAllRefinement();
    bool passGlobal  = testAllGlobal();
    bool passLocal   = testAllPWLocal();
    bool passWavelet = testAllWavelet();

    bool pass = passGlobal && passLocal && passWavelet && passRefine && passDomain && passAccel;
    //bool pass = passAccel;

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
        double *points = grid->getPoints(), *weights = 0;
        if (type == type_integration){
            weights = grid->getQuadratureWeights();
        }else{
            weights = grid->getInterpolationWeights(x);
        }

        double *y = new double[num_outputs];
        double *r = new double[num_outputs];  std::fill(r, r + num_outputs, 0.0);
//      Sequential version
//        for(int i=0; i<num_points; i++){
//            f->eval(&(points[i*num_dimensions]), y);
//            for(int k=0; k<num_outputs; k++){
//                r[k] += weights[i] * y[k];
//            }
//        }

        #pragma omp parallel
        {
            double *y_local = new double[num_outputs];
            double *r_local = new double[num_outputs];  std::fill(r_local, r_local + num_outputs, 0.0);
            #pragma omp for
            for(int i=0; i<num_points; i++){
                f->eval(&(points[i*num_dimensions]), y_local);
                for(int j=0; j<num_outputs; j++){
                    r_local[j] += weights[i] * y_local[j];
                };
            };

            #pragma omp critical
            {
                for(int j=0; j<num_outputs; j++){
                    r[j] += r_local[j];
                }
            }
            delete[] y_local;
            delete[] r_local;
        }

        double err = 0.0;
        if (type == type_integration){
            f->getIntegral(y);
        }else{
            f->eval(x, y);
        }
        for(int j=0; j<num_outputs; j++){
            err += fabs(y[j] - r[j]);
        };
        R.error = err;

        delete[] r;
        delete[] y;
        delete[] points;
        delete[] weights;
    }else if (type == type_internal_interpolation){
        // load needed points
        int num_needed_points = grid->getNumNeeded();
        if (num_needed_points > 0){
            double *needed_points = grid->getNeededPoints();
            double *values = new double[num_outputs * num_needed_points];

            for(int i=0; i<num_needed_points; i++){
                f->eval(&(needed_points[i*num_dimensions]), &(values[i*num_outputs]));
            }

            grid->loadNeededPoints(values);

            delete[] values;
            delete[] needed_points;
        }

        double *r = new double[num_outputs];  std::fill(r, r + num_outputs, 0.0);
        double *n = new double[num_outputs];  std::fill(n, n + num_outputs, 0.0);

        #pragma omp parallel
        {
        	double *x_local = new double[num_dimensions];
        	double *n_local = new double[num_dimensions];
        	double *r_local = new double[num_dimensions];
        	double *y       = new double[num_dimensions];
        	double *s       = new double[num_dimensions];

        	for(int i=0; i<num_outputs; i++){
        		y[i] = s[i] = r_local[i] = n_local[i] = 0.0;
        	}
            #pragma omp for
        	for(int k=0; k<num_mc; k++){
        		setRandomX(num_dimensions, x_local);
        		f->eval(x_local, y);
        		//std::fill(s, s + num_dimensions, 0.0);
        		grid->evaluate(x_local, s);
        		for(int j=0; j<num_outputs; j++){
                    double e = fabs(y[j] - s[j]);
                    if (r_local[j] < e) r_local[j] = e;
                    e = fabs(y[j]);
                    if (n_local[j] < e) n_local[j] = e;
        		}
        	}

            #pragma omp critical
        	{
        		for(int i=0; i<num_outputs; i++){
        			if (r[i] < r_local[i]) r[i] = r_local[i];
        			if (n[i] < n_local[i]) n[i] = n_local[i];
        		}
        	}

        	delete[] x_local;
        	delete[] n_local;
        	delete[] r_local;
        	delete[] y;
        	delete[] s;
        }

        double err = r[0] / n[0];
        for(int j=1; j<num_outputs; j++){
            err = (err > r[j] / n[j]) ? err : r[j] / n[j];
        }

        delete[] r;
        delete[] n;
        R.error = err;
    }
    R.num_points = grid->getNumPoints();
    return R;
}

bool ExternalTester::testGlobalRule(const BaseFunction *f, TasGrid::TypeOneDRule rule, const int *anisotropic, double alpha, double beta, const bool interpolation, const int depths[], const double tols[]) const{
    TasGrid::TasmanianSparseGrid grid;
    TestResults R;
    int num_global_tests = (interpolation) ? 3 : 1;
    TestType tests[3] = { type_integration, type_nodal_interpolation, type_internal_interpolation };
    TasGrid::TypeDepth type = TasGrid::type_iptotal;
    double *x = new double[f->getNumInputs()]; setRandomX(f->getNumInputs(),x);
    bool bPass = true;
    const char *custom_filename = (rule == rule_customtabulated) ? "GaussPattersonRule.table" : 0;
    for(int i=0; i<num_global_tests; i++){
        if (interpolation){
            grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), depths[i], type, rule, anisotropic, alpha, beta, custom_filename);
            R = getError(f, &grid, tests[i], x);
        }else{
            grid.makeGlobalGrid(f->getNumInputs(), 0, depths[i], type, rule, anisotropic, alpha, beta, custom_filename);
            R = getError(f, &grid, type_integration);
        }
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
    }
    //cout << "HERE 1" << TasGrid::OneDimensionalMeta::getIORuleString(rule) << endl;
    if (TasGrid::OneDimensionalMeta::isSequence(rule)){
        //cout << "HERE 2" << TasGrid::OneDimensionalMeta::getIORuleString(rule) << endl;
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
    }else if (rule == TasGrid::rule_chebyshev){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_chebyshev;
        const int depths1[3] = { 22, 22, 22 };
        const double tols1[3] = { 1.E-12, 1.E-10, 1.E-10 };
        const int depths2[3] = { 22, 22, 22 };
        const double tols2[3] = { 1.E-12, 1.E-09, 1.E-09 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_chebyshevodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_chebyshevodd;
        const int depths1[3] = { 22, 22, 22 };
        const double tols1[3] = { 1.E-12, 1.E-10, 1.E-10 };
        const int depths2[3] = { 22, 22, 22 };
        const double tols2[3] = { 1.E-12, 1.E-09, 1.E-09 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_leja){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_leja;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_lejaodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_lejaodd;
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
    }else if (rule == TasGrid::rule_rlejadouble2){
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
    }else if (rule == TasGrid::rule_rlejadouble4){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejadouble4;
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
    }else if (rule == TasGrid::rule_maxlebesgue){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_maxlebesgue;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_maxlebesgueodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_maxlebesgueodd;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_minlebesgue){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_minlebesgue;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_minlebesgueodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_minlebesgueodd;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_mindelta){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_mindelta;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_mindeltaodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_mindeltaodd;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
        const int depths2[3] = { 20, 20, 20 };
        const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausslegendre){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausslegendre;
        const int depths1[3] = { 20, 36, 38 };
        const double tols1[3] = { 1.E-10, 1.E-07, 1.E-07 };
        const int depths2[3] = { 24, 36, 36 };
        const double tols2[3] = { 1.E-10, 1.E-07, 1.E-07 };
        if (testGlobalRule(&f21nx2, oned, 0, alpha, beta, true, depths1, tols1) && testGlobalRule(&f21cos, oned, 0, alpha, beta, true, depths2, tols2)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausslegendreodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausslegendreodd;
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
        std::ifstream ftest("GaussPattersonRule.table");
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
    }else if (rule == TasGrid::rule_gausschebyshev1){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausschebyshev1;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 5.E-14, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGC1, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausschebyshev1odd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausschebyshev1odd;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 2.E-14, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGC1, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausschebyshev2){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausschebyshev2;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 1.1E-14, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGC2, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausschebyshev2odd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausschebyshev2odd;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 1.1E-14, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGC2, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gaussgegenbauer){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gaussgegenbauer;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 1.E-11, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGG, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gaussgegenbauerodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gaussgegenbauerodd;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 1.E-11, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGG, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gaussjacobi){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gaussjacobi;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 1.E-08, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGJ, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gaussjacobiodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gaussjacobiodd;
        const int depths1[3] = { 20, 20, 20 };
        const double tols1[3] = { 1.E-08, 1.E-05, 1.E-05 };
        if (testGlobalRule(&f21constGJ, oned, 0, alpha, beta, true, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausslaguerre){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausslaguerre;
        const int depths1[1] = { 20 };
        const double tols1[1] = { 1.E-08 };
        if (testGlobalRule(&f21constGGL, oned, 0, alpha, beta, false, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausslaguerreodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausslaguerreodd;
        const int depths1[1] = { 20 };
        const double tols1[1] = { 1.E-08 };
        if (testGlobalRule(&f21constGGL, oned, 0, alpha, beta, false, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausshermite){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausshermite;
        const int depths1[1] = { 20 };
        const double tols1[1] = { 1.E-09 };
        if (testGlobalRule(&f21constGH, oned, 0, alpha, beta, false, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
    }else if (rule == TasGrid::rule_gausshermiteodd){
        { TasGrid::TypeOneDRule oned = TasGrid::rule_gausshermiteodd;
        const int depths1[1] = { 20 };
        const double tols1[1] = { 1.E-09 };
        if (testGlobalRule(&f21constGH, oned, 0, alpha, beta, false, depths1, tols1)){
            if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
        }else{
            cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl;  pass = false;
        }}
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
        //grid.enableAcceleration(accel_gpu_cuda);
        //grid.setGPUID(0);
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
    // TODO: double-check the piece-wise constant rule, why it fails on Windows
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
    { TasGrid::TypeOneDRule oned = TasGrid::rule_localp0;
    const int depths1[18] = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
    const double tols1[18] = { 1.E-03, 5.E-01, 5.E-01, 1.E-03, 2.E-04, 2.E-04, 1.E-09, 1.E-06, 1.E-06, 1.E-09, 3.E-08, 3.E-08, 1.E-09, 4.E-09, 4.E-09, 1.E-09, 4.E-09, 4.E-09 };
    if (testLocalPolynomialRule(&f21coscos, oned, depths1, tols1)){
        if (verbose) cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "Pass" << endl;
    }else{
        cout << setw(wfirst) << "Rule" << setw(wsecond) << TasGrid::OneDimensionalMeta::getIORuleString(oned) << setw(wthird) << "FAIL" << endl; pass = false;
    }}
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
    }/*{
        const BaseFunction *f = &f21nx2;
        grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 2, 0, rule_semilocalp);
        int np[5] =     {    21,    81,   297, 1053,  3637 };
        double err[5] = { 3.E-1, 2.E-1, 6.E-2, 3E-2, 8.E-3 };
        if (!testSurplusRefinement(f, &grid, 1.E-4, refine_classic, np, err, 5)){
            cout << "ERROR: failed semi-localp classic refinement for " << f->getDescription() << endl;  pass = false;
        }
    }*/{
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
            //cout << "num points = " << R.num_points << "  error = " << R.error << endl;
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
            //cout << "num points = " << R.num_points << "  error = " << R.error << endl;
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
            //cout << "num points = " << R.num_points << "  error = " << R.error << endl;
        }
    }

    cout << "      Domain                      anisotropic" << setw(15) << ((pass1) ? "Pass" : "FAIL") << endl;

    bool pass2 = true;

    {
        const BaseFunction *f = &f21expDomain;
        double errs[5] = {3.E-5, 2.E-7, 2.E-8, 6.E-10, 4.E-11 };
        double transform_a[2] = { 3.0, -3.0 };
        double transform_b[2] = { 4.0,  2.0 };
        for(int i=0; i<5; i++){
            grid.makeSequenceGrid(f->getNumInputs(), f->getNumOutputs(), i+5, TasGrid::type_level, TasGrid::rule_leja);
            grid.setDomainTransform(transform_a, transform_b);
            TestResults R = getError(f, &grid, type_integration);
            if (R.error>errs[i]){
                cout << "Using leja rule" << endl;
                cout << "Failed domain transform test for " << f->getDescription() << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass2 = false;
            }
            //cout << "  error = " << R.error << endl;
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
            //cout << "  error = " << R.error << endl;
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
            //cout << "  error = " << R.error << endl;
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
                cout << "Using gauss-legendre rule" << endl;
                cout << "Failed domain transform test for " << f->getDescription() << "   error = " << R.error << "  expected: " << errs[i] << endl;
                     pass2 = false;
            }
            //cout << "  error = " << R.error << endl;
        }
    }

    cout << "      Domain                      transformed" << setw(15) << ((pass2) ? "Pass" : "FAIL") << endl;

    bool pass3 = true;

    {
        TasmanianSparseGrid gridc;
        const BaseFunction *f = &f21conformal;
        int asin_conformal[2] = { 4, 4 };
        for(int l=0; l<6; l++){
            grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), l+2, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
            gridc.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), l+2, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
            gridc.setConformalTransformASIN(asin_conformal);
            TestResults R1 = getError(f, &grid, type_internal_interpolation);
            TestResults R2 = getError(f, &gridc, type_internal_interpolation);
            //cout << R1.num_points << "  " << R2.num_points << endl;
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
        }
    }{
        TasmanianSparseGrid gridc;
        const BaseFunction *f = &f21conformal;
        int asin_conformal[2] = { 4, 4 };
        for(int l=0; l<15; l++){
            grid.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), l+5, TasGrid::type_iptotal, TasGrid::rule_gausspatterson);
            gridc.makeGlobalGrid(f->getNumInputs(), f->getNumOutputs(), l+5, TasGrid::type_iptotal, TasGrid::rule_gausspatterson);
            gridc.setConformalTransformASIN(asin_conformal);
            TestResults R1 = getError(f, &grid, type_integration);
            TestResults R2 = getError(f, &gridc, type_integration);
            //cout << R1.num_points << "  " << R2.num_points << endl;
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
            //cout << R1.num_points << "  " << R2.num_points << endl;
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
        }
    }

    cout << "      Domain                        conformal" << setw(15) << ((pass3) ? "Pass" : "FAIL") << endl;

    return (pass1 && pass2 && pass3);
}

bool ExternalTester::testAcceleration(const BaseFunction *f, TasmanianSparseGrid *grid) const{
    int dims = f->getNumInputs();
    int outs = f->getNumOutputs();
    //grid->printStats();
    if (grid->getNumNeeded() > 0){
        double *points = grid->getNeededPoints();
        double *vals = new double[outs * grid->getNumPoints()];

        for(int i=0; i<grid->getNumPoints(); i++){
            f->eval(&(points[i*dims]), &(vals[i*outs]));
        }
        grid->loadNeededPoints(vals);
        delete[] vals;
        delete[] points;
    }

    int num_x = 256; // for batched evaluations
    int num_fast = 16; // for fast evaluations, must be <= num_x
    double *x = new double[num_x * dims];
    setRandomX(num_x * dims, x);

    double *baseline_y = new double[outs * num_x];
    double *test_y = new double[outs * num_x];
    for(int i=0; i<num_x; i++) grid->evaluate(&(x[i*dims]), &(baseline_y[i*outs]));

    bool pass = true;
    TypeAcceleration acc[4] = {accel_none, accel_cpu_blas, accel_gpu_cublas, accel_gpu_cuda};
    int gpuID = 0;
    int c = 0;
    while(c < 4){
        //cout << "Settin c = " << c << endl;
        grid->enableAcceleration(acc[c]);
        if (c > 1){ // gpu test
            grid->setGPUID(gpuID);
        }
        //grid->printStats();

        for(int i=0; i<outs * num_x; i++ ) test_y[i] = 0.0;

        grid->evaluateBatch(x, num_x, test_y);

        double err = 0.0;
        for(int i=0; i<outs*num_x; i++) if (fabs(test_y[i] - baseline_y[i]) > err) err = fabs(test_y[i] - baseline_y[i]);

        if (err > 1.E-11){
            pass = false;
            cout << "Failed Batch evaluation for acceleration c = " << c << " gpuID = " << gpuID << endl;
            cout << "Observed error: " << err << " for function: " << f->getDescription() << endl;
            grid->printStats();
        }

        for(int i=0; i<num_fast; i++){
            grid->evaluateFast(&(x[i*dims]), &(test_y[i*outs]));
        }

        err = 0.0;
        for(int i=0; i<outs*num_fast; i++) if (fabs(test_y[i] - baseline_y[i]) > err) err = fabs(test_y[i] - baseline_y[i]);

        if (err > 1.E-11){
            pass = false;
            cout << "Failed Fast evaluation for acceleration c = " << c << " gpuID = " << gpuID << endl;
            cout << "Observed error: " << err << " for function: " << f->getDescription() << endl;
            grid->printStats();
        }

        if (c > 1){ // gpu test
            gpuID++;
            if (gpuID >= grid->getNumGPUs()){
                gpuID = 0;
                c++;
            }
        }else{
            c++;
        }
    }

    delete[] test_y;
    delete[] baseline_y;
    delete[] x;

    return pass;
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

    grid.makeLocalPolynomialGrid(f->getNumInputs(), f->getNumOutputs(), 7, 1, TasGrid::rule_localp);
    pass = pass && testAcceleration(f, &grid);
    if (pass){
        if (verbose) cout << "      Accelerated" << setw(wsecond) << "local polynomial" << setw(wthird) << "Pass" << endl;
    }else{
        cout << "      Accelerated" << setw(wsecond) << "local polynomial" << setw(wthird) << "FAIL" << endl;
    }

    //cout << "      Domain                      anisotropic" << setw(15) << ((pass) ? "Pass" : "FAIL") << endl;
    cout << "      Acceleration                        all" << setw(15) << ((pass) ? "Pass" : "FAIL") << endl;

    return pass;
}

#include <sys/time.h>
extern "C" double gettime(){
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1.0E-6;
}
void loadGridValues(TasmanianSparseGrid *grid){
    int dims = grid->getNumDimensions(), num_outputs = grid->getNumOutputs();
    int num_points = grid->getNumNeeded();
    if (num_points == 0) return;
    double *x = grid->getNeededPoints();
    double *v = new double[num_outputs * num_points];
    double s = 1.0 / ((double) (dims*dims));
    //#pragma omp parallel for
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
        //cout << exp(-nx2) << "   " << sin(nx2) << endl;
        //cout << x[2*i] << "   " << x[2*i+1] << endl;
    }
    //cout << "Computing surplusses = " << num_outputs << endl;
    //double start = gettime();
    grid->loadNeededPoints(v);
    //cout << gettime() - start << endl;
}

void ExternalTester::benchmark(int argc, const char **argv){
    if (strcmp(argv[2],"alpha") == 0){
        cout << "./tasgrid -bench alpha <dims> <outs> <depth> <type> <rule> <batch size> <iterations> <gpu>" << endl;
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
        TypeAcceleration gpu_tests[3] = {accel_cpu_blas, accel_gpu_cublas, accel_gpu_cuda};
        TypeAcceleration *tests;
        int num_tests;
        if (gpu > -1){
            num_tests = 3;
            tests = gpu_tests;
            //cout << setw(width) << "cpu_blas" << setw(width) << "gpu_cublas" << setw(width) << "gpu_cuda" << setw(width) << "gpu_magma" << endl;
            cout << setw(width) << "cpu_blas" << setw(width) << "gpu_cublas" << setw(width) << "gpu_cuda" << endl;
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
            for(int s=0; s<5; s++){
                if (OneDimensionalMeta::isLocalPolynomial(r)){
                    grid->makeLocalPolynomialGrid(dims, outs, depth, 2, r);
                }else if (OneDimensionalMeta::isSequence(r)){
                    grid->makeSequenceGrid(dims, outs, depth, d, r);
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
                    if ((num_tests == 2) || (t > 0))
                    for(int i=0; i<num_runs; i++){
                        grid->evaluateBatch(x, num_x, y);
                    }
                    cout << setw(width) << ((int) ((gettime() - start) * 100.0));
                }
                outs *= 2;
                cout << setw(width) << "seconds^{-2}" << endl;
                delete[] y;
            }
            delete[] x;
        }else{ // compare compute surpluses
            double *v = new double[np * outs * 16];
            cout << std::scientific;
            cout.precision(5);
            setRandomX(np * outs * 16, v);
            for(int s=0; s<5; s++){
                if (OneDimensionalMeta::isLocalPolynomial(r)){
                    grid->makeLocalPolynomialGrid(dims, outs, depth, 2, r);
                }else if (OneDimensionalMeta::isSequence(r)){
                    grid->makeSequenceGrid(dims, outs, depth, d, r);
                }else{
                    grid->makeGlobalGrid(dims, outs, depth, d, r);
                }
                cout << setw(width) << outs;
                double start;

                for(int t=0; t<num_tests; t++){
                    grid->enableAcceleration(tests[t]);
                    if ((gpu > -1) && (t > 0)) grid->setGPUID(gpu);
                    start = gettime();
                    if ((num_tests == 2) || (t > 0) || (num_x == 0))
                    for(int i=0; i<num_runs; i++){
                        grid->loadNeededPoints(v);
                    }
                    cout << setw(width) << ((int) ((gettime() - start) * 100.0));
                }
                outs *= 2;
                cout << setw(width) << "seconds^{-2}" << endl;
            }
            delete[] v;
        }
    }
}

void ExternalTester::debugTest(){
    cout << "Debug Test" << endl;
    cout << "Put here testing code and call this with ./tasgrid -test debug" << endl;

    double x[2] = { 0.33, -0.1133 };
    double y1[1], y2[1];

    TasmanianSparseGrid *grid = new TasmanianSparseGrid();
    grid->makeLocalPolynomialGrid(2, 1, 3, 1, rule_localp);
    grid->enableAcceleration(accel_gpu_cublas);
    grid->setGPUID(0);
    loadGridValues(grid);
    grid->evaluate(x, y1);

    cout << "-------------" << endl;
    grid->makeLocalPolynomialGrid(2, 1, 3, 1, rule_localp);
    loadGridValues(grid);
    grid->evaluate(x, y2);

    cout << std::scientific;
    cout.precision(16);
    cout << y1[0] - y2[0] << endl;

}

void ExternalTester::debugTestII(){
    cout << "Debug Test II" << endl;
    cout << "Put here testing code and call this with ./tasgrid -test db" << endl;

}

#endif
