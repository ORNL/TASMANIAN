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

#include <iostream>

#include "tasdreamExternalTests.hpp"

using namespace std;

//enum TypeHelp{
//    help_generic,
//    help_benchmark
//};

//void printHelp(TypeHelp ht = help_generic);

int main(int argc, const char ** argv){
    //cout << " Phruuuuphrrr " << endl; // this is the sound that the Tasmanian devil makes

    bool debug = false;
    TypeDREAMTest test = test_all;

    DreamExternalTester tester;
    int k = 1;
    while(k < argc){
        if ((strcmp(argv[k],"verbose") == 0) || (strcmp(argv[k],"-verbose") == 0)){
            tester.showVerbose();
            tester.showStatsValues();
        }else if ((strcmp(argv[k],"v") == 0) || (strcmp(argv[k],"-v") == 0)){
            tester.showVerbose();
        }else if ((strcmp(argv[k],"random") == 0) || (strcmp(argv[k],"-random") == 0)){
            tester.useRandomRandomSeed();
        }else if ((strcmp(argv[k],"debug") == 0)) debug = true;
        else if ((strcmp(argv[k],"analytic") == 0)) test = test_analytic;
        else if ((strcmp(argv[k],"posterior") == 0)) test = test_posterior;
        else{
            cerr << "ERROR: Unknown option '" << argv[k] << "'" << endl;
            cerr << "   to see list of available options use: ./gridtest --help" << endl;
            return 1;
        }
        k++;
    }

    if (debug){
        testDebug();
        return 0;
    }

    return (tester.performTests(test)) ? 0 : 1;
}
