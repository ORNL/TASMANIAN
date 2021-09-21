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

#include "testConstructSurrogate.hpp"
#include "testLoadUnstructured.hpp"
#include "testExoticQuadrature.hpp"

void debugTest() {
    // cout << "Debug Test (callable from the CMake build folder)" << endl;
    // cout << "Put testing code here and call with ./Addons/addontester debug" << endl;

    // Set weight function.
    double freq = 1.0;
    double phase_shift = 1.0;
    int dimension = 1;
    int level = 40;
    int shift = 1.0;
    int depth = 20;
    bool symmetric_weights = false;
    // Get exotic nodes and weights for the sinc function.
    auto sinc = [freq, phase_shift](double x)->double {return (fabs(freq * (x - phase_shift)) <= 1e-20) ?
                1.0 :
                sin(freq * (x - phase_shift)) / (freq * (x - phase_shift));};
    TasGrid::CustomTabulated ct1, ct2;
    TasGrid::gauss_quadrature_version = 1;
    ct1 = TasGrid::getExoticQuadrature(level, shift, sinc, 50 * depth + 1, "Sinc-Weighted Quadrature 1", symmetric_weights);
    TasGrid::gauss_quadrature_version = 2;
    ct2 = TasGrid::getExoticQuadrature(level, shift, sinc, 50 * depth + 1, "Sinc-Weighted Quadrature 2", symmetric_weights);
    // Form the grid.
    TasGrid::TasmanianSparseGrid sg1, sg2;
    sg1.makeGlobalGrid(dimension, 1, depth, TasGrid::type_qptotal, std::move(ct1));
    std::vector<double> quad_points1 = sg1.getPoints();
    std::vector<double> quad_weights1 = sg1.getQuadratureWeights();
    sg2.makeGlobalGrid(dimension, 1, depth, TasGrid::type_qptotal, std::move(ct2));
    std::vector<double> quad_points2 = sg2.getPoints();
    std::vector<double> quad_weights2 = sg2.getQuadratureWeights();
    // Print
    std::cout << "n1, n2" << std::endl;
    std::cout << quad_points1.size() << ", " << quad_points2.size() << "\n" << std::endl;
    if (!(quad_points1.size() == quad_points2.size())) {
        std::cout << "Node/Weight lengths differ!" << std::endl;
    }
    size_t n = std::min(quad_points1.size(), quad_points2.size()); 
    std::cout << "      p1[],       p2[],        err[]" << std::endl;
    for (size_t i=0; i<n; i++)
        std::cout << std::fixed << std::setw(10) << quad_points1[i] << ", " << std::setw(10) << quad_points2[i] << ", "
                  << std::setw(10) << std::scientific << std::fabs(quad_points1[i] - quad_points2[i]) << std::fixed << std::endl;
    std::cout << endl;
    std::cout << "      w1[],       w2[],        err[]" << std::endl;
    for (size_t i=0; i<n; i++)
        std::cout << std::fixed << std::setw(10) << quad_weights1[i] << ", " << std::setw(10) << quad_weights2[i] << ", "
                  << std::setw(10) << std::scientific << std::fabs(quad_weights1[i] - quad_weights2[i]) << std::fixed << std::endl;
    std::cout << endl;
}

int main(int argc, char const **argv){

    std::deque<std::string> args = stringArgs(argc, argv);

    bool debug = false;
    bool pass_all = true;
    bool verbose = false;
    int gpuid = -1;

    while(not args.empty()){
        if (args.front() == "debug") debug = true;
        if (hasInfo(args.front())) verbose = true;
        if (hasGpuID(args.front())){
            args.pop_front();
            gpuid = getGpuID(args);
        }
        args.pop_front();
    }

    if (debug){
        debugTest();
        return 0;
    }

    cout << "\n\n";
    cout << "---------------------------------------------------------------------" << endl;
    cout << "          Tasmanian Addons Module: Functionality Test" << endl;
    cout << "---------------------------------------------------------------------" << endl << endl;

    bool pass = testConstructSurrogate(verbose);
    cout << std::setw(40) << "Automated construction" << std::setw(10) << ((pass) ? "Pass" : "FAIL") << endl;
    pass_all = pass_all && pass;

    #if defined(Tasmanian_ENABLE_BLAS) || defined(Tasmanian_ENABLE_GPU)
    pass = true;
    pass = testLoadUnstructuredL2(verbose, gpuid);
    cout << std::setw(40) << "Unstructured construction" << std::setw(10) << ((pass) ? "Pass" : "FAIL") << endl;
    pass_all = pass_all && pass;
    #else
    cout << std::setw(40) << "Unstructured construction" << std::setw(10) << "skipping" << endl;
    gpuid *= 2; // no op to register the use of gpuid
    #endif

    // Tests for Exotic Quadrature.
    #ifdef Tasmanian_ENABLE_BLAS
    pass = true;
    pass = testExoticQuadrature();
    cout << std::setw(40) << "Exotic quadrature" << std::setw(10) << ((pass) ? "Pass" : "FAIL") << endl;
    pass_all = pass_all && pass;
    #else
    cout << std::setw(40) << "Exotic quadrature" << std::setw(10) << "skipping" << endl;
    #endif

    cout << "\n";
    if (pass){
        cout << "---------------------------------------------------------------------" << endl;
        cout << "               All Tests Completed Successfully" << endl;
        cout << "---------------------------------------------------------------------" << endl << endl;
    }else{
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl;
        cout << "         Some Tests Have Failed" << endl;
        cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl << endl;
    }

    return ((pass_all) ? 0 : 1);

}
