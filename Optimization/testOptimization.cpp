/*
 * Copyright (c) 2022, Miroslav Stoyanov & Weiwei Kong
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
 * conditions are met:
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
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND
 * IMPLIED. THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
 * THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL
 * ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE. THE USER ASSUMES
 * RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING
 * FROM OR ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */


#include "TasmanianOptimization.hpp"
#include "gridtestCLICommon.hpp"

void debugTest() {
    cout << "Debug Test (callable from the CMake build folder)" << endl;
    cout << "Put testing code here and call with ./Optimization/optimizationtester debug" << endl;

    // Objective function and stopping criterion.
    TasOptimization::ObjectiveFunction shc =
            [&](const std::vector<double> x) {
                return (4 - 2.1 * x[0]*x[0] + x[0]*x[0]*x[0]*x[0] / 3) * x[0]*x[0] +
                        x[0] * x[1] +
                        (-4 + 4 * x[1]*x[1]) * x[1]*x[1];};
    TasOptimization::StoppingCondition le_100_iter =
            [&](TasOptimization::OptimizationState os) {return (os.num_iterations > 100);};

    // State object.
    std::vector<double> x0 = {0.0, 0.0};
    TasOptimization::OptimizationState os(shc, x0);

    std::cout << "f(x0) = " << shc(x0) << std::endl;

    // Algorithm parameters.
    std::vector<double> lower = {-3, -2};
    std::vector<double> upper = {3, 2};
    double inertia_weight = 0.5;
    double cognitive_coeff = 2.0;
    double social_coeff = 2.0;
    int num_particles = 100;

    // Main run.
    TasOptimization::optimizeParticleSwarm(os, le_100_iter, lower, upper, inertia_weight, cognitive_coeff, social_coeff, num_particles);

    std::cout << "f(x) = " << shc(os.x) << std::endl;

}

int main(int argc, char const **argv){

    std::deque<std::string> args = stringArgs(argc, argv);

    bool debug = false;
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

}

