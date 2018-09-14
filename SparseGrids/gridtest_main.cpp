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

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <string.h>
#include <math.h>

#include "TasmanianSparseGrid.hpp"
#include "tasgridExternalTests.hpp"
#include "tasgridUnitTests.hpp"

using namespace std;
using namespace TasGrid;

int main(int argc, const char ** argv){

    //cout << " Phruuuuphrrr " << endl; // this is the sound that the Tasmanian devil makes

    // tests for speed
    if ((argc > 1) && ((strcmp(argv[1],"-benchmark") == 0) || (strcmp(argv[1],"-bench") == 0))){
        ExternalTester tester(1);
        tester.resetRandomSeed();
        tester.benchmark(argc, argv);
        return 0;
    }

    // testing
    bool debug = false;
    bool debugII = false;
    bool verbose = false;
    bool seed_reset = false;

    TestList test = test_all;
    UnitTests utest = unit_none;

    int gpuid = -1;
    int k = 1;
    while (k < argc){
        if ((strcmp(argv[k],"debug") == 0)) debug = true;
        else if ((strcmp(argv[k],"db") == 0)) debugII = true;
        else if ((strcmp(argv[k],"verbose") == 0)) verbose = true;
        else if ((strcmp(argv[k],"v") == 0)) verbose = true;
        else if ((strcmp(argv[k],"-v") == 0)) verbose = true;
        else if ((strcmp(argv[k],"-random") == 0)) seed_reset = true;
        else if ((strcmp(argv[k],"random") == 0)) seed_reset = true;
        else if ((strcmp(argv[k],"acceleration") == 0)) test = test_acceleration;
        else if ((strcmp(argv[k],"domain") == 0)) test = test_domain;
        else if ((strcmp(argv[k],"refinement") == 0)) test = test_refinement;
        else if ((strcmp(argv[k],"global") == 0)) test = test_global;
        else if ((strcmp(argv[k],"local") == 0)) test = test_local;
        else if ((strcmp(argv[k],"wavelet") == 0)) test = test_wavelet;
        else if ((strcmp(argv[k],"fourier") == 0)) test = test_fourier;
        else if ((strcmp(argv[k],"errors") == 0)) utest = unit_except;
        else if ((strcmp(argv[k],"api") == 0)) utest = unit_api;
        else if ((strcmp(argv[k],"c") == 0)) utest = unit_c;
        else if ((strcmp(argv[k],"cover") == 0)) utest = unit_cover;
        else if ((strcmp(argv[k],"-gpuid") == 0)){
            if (k+1 >= argc){
                cerr << "ERROR: -gpuid requires a valid number!" << endl;
                return 1;
            }
            gpuid = atoi(argv[k+1]);
            k++;
            if ((gpuid < -1) || (gpuid >= TasmanianSparseGrid::getNumGPUs())){
                cerr << "ERROR: -gpuid " << gpuid << " is not a valid gpuid!" << endl;
                cerr << "      see ./tasgrid -v for a list of detected GPUs." << endl;
                return 1;
            }
        }else{
            cerr << "ERROR: unknown option " << argv[k] << endl;
            return 1;
        }
        k++;
    }

    ExternalTester tester(1000);
    GridUnitTester utester;
    tester.setGPUID(gpuid);
    bool pass = true;
    if (debug){
        tester.debugTest();
    }else if (debugII){
        tester.debugTestII();
    }else{
        if (verbose) tester.setVerbose(true);
        if (verbose) utester.setVerbose(true);

        if (seed_reset) tester.resetRandomSeed();

        if (utest == unit_none){
            if (test == test_all) pass = pass && utester.Test(unit_all);
            pass = pass && tester.Test(test);
        }else{
            pass = pass && utester.Test(utest);
        }
    }
    return (pass) ? 0 : 1;
}
