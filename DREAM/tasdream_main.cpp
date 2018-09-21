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
#include <time.h>
#include <string.h>

#include "TasmanianDREAM.hpp"

#include "TasmanianSparseGrid.hpp"

#include "tasdreamExternalTests.hpp"
#include "tasdreamBenchmark.hpp"

using namespace std;
using namespace TasGrid;
using namespace TasDREAM;


enum TypeHelp{
    help_generic,
    help_benchmark
};

void printHelp(TypeHelp ht = help_generic);

int main(int argc, const char ** argv){
    //cout << " Phruuuuphrrr " << endl; // this is the sound that the Tasmanian devil makes

    if (argc < 2){
        cerr << "ERROR: no command specified" << endl << endl;
        printHelp();
        return 1;
    }

    // basic help
    if ((strcmp(argv[1],"--help") == 0)||(strcmp(argv[1],"-help") == 0)||(strcmp(argv[1],"help") == 0)||(strcmp(argv[1],"-h") == 0)){
        printHelp();
        return 0;
    }

    if (strcmp(argv[1],"-test") == 0){
        bool debug = false;
        bool debugII = false;
        bool verbose = false;
        bool seed_reset = false;

        TestList test = test_all;

        int k = 2;
        while (k < argc){
            if ((strcmp(argv[k],"debug") == 0)) debug = true;
            if ((strcmp(argv[k],"db") == 0)) debugII = true;
            if ((strcmp(argv[k],"verbose") == 0)) verbose = true;
            if ((strcmp(argv[k],"v") == 0)) verbose = true;
            if ((strcmp(argv[k],"-v") == 0)) verbose = true;
            if ((strcmp(argv[k],"-random") == 0)) seed_reset = true;
            if ((strcmp(argv[k],"random") == 0)) seed_reset = true;
            if ((strcmp(argv[k],"analytic") == 0)) test = test_analytic;
            if ((strcmp(argv[k],"model") == 0)) test = test_model;
            k++;
        }

        ExternalTester tester(1000);
        bool pass = true;
        if (debug){
            tester.Test(test);
        }else if (debugII){
            tester.Test(test);
        }else{
            if (verbose) tester.setVerbose(true);
            if (seed_reset) tester.resetRandomSeed();
            pass = tester.Test(test);
        }
        return (pass) ? 0 : 1;
    }

    if ((strcmp(argv[1],"-benchmark") == 0) || (strcmp(argv[1],"-bench") == 0)){
        if (argc < 3){
            cerr << "ERROR: missing benchmark name" << endl << endl;
            printHelp(help_benchmark);
            return 1;
        }
        if ((strcmp(argv[2],"--help") == 0) || (strcmp(argv[2],"-help") == 0) || (strcmp(argv[2],"help") == 0) || (strcmp(argv[2],"-h") == 0)){
            printHelp(help_benchmark);
            return 0;
        }
        if (strcmp(argv[2],"basic-alpha") == 0){
            if ((argc > 3) && ((strcmp(argv[3],"--help") == 0) || (strcmp(argv[3],"-help") == 0) || (strcmp(argv[3],"help") == 0) || (strcmp(argv[3],"-h") == 0))){
                printHelp(help_benchmark);
                return 0;
            }
            if ((argc != 9) && (argc != 10)){
                cerr << "ERROR: wrong parameter count" << endl << endl;
                printHelp(help_benchmark);
                return 1;
            }
            int num_outputs = atoi(argv[3]);
            int depth = atoi(argv[4]);
            int num_chains = atoi(argv[5]);
            int num_burnup = atoi(argv[6]);
            int num_mcmc = atoi(argv[7]);
            int gpuID = atoi(argv[8]);
            const char *outfilename = (argc == 10) ? argv[9] : 0;
            sharedBenchmarkBasicAlpha(num_outputs, depth, num_chains, num_burnup, num_mcmc, gpuID, outfilename);
            return 0;
        }else if (strcmp(argv[2],"basic-alpha-mpi") == 0){
            #ifdef MPI_VERSION
            if ((argc > 3) && ((strcmp(argv[3],"--help") == 0)||(strcmp(argv[3],"-help") == 0)||(strcmp(argv[3],"help") == 0)||(strcmp(argv[3],"-h") == 0))){
                printHelp(help_benchmark);
                return 0;
            }
            if ((argc != 8) && (argc != 9)){
                cerr << "ERROR: wrong parameter count" << endl << endl;
                printHelp(help_benchmark);
                return 1;
            }
            int num_outputs = atoi(argv[3]);
            int depth = atoi(argv[4]);
            int num_chains = atoi(argv[5]);
            int num_burnup = atoi(argv[6]);
            int num_mcmc = atoi(argv[7]);
            const char *outfilename = (argc == 9) ? argv[8] : 0;
            mpiBenchmarkBasicAlpha(num_outputs, depth, num_chains, num_burnup, num_mcmc, outfilename);
            return 0;
            #else
            cerr << "ERROR: benchmark 'basic-alpha-mpi' requires MPI support!" << endl;
            return 1;
            #endif
        }
    }else{
        cerr << "Unknown command!" << argv[1] << endl;
        return 1;
    }

    return 0;
}

void printHelp(TypeHelp ht){
    cout << endl;

    if (ht == help_generic){

        cout << " Usage: tasdream <command> <option1> <value1> <option2> <value2> ... " << endl << endl;

        cout << " -help\t"         << "\t\t-h,--help" << "\tshow verbose help options" << endl;
        cout << " -test\t"         << "\t\t\t"    << "\tperform a number of internal tests" << endl;

    }else if (ht == help_benchmark){

        cout << " Usage: tasdream -bench <benchmark_name> <parameter1> <parameter2> <parameter3> ..." << endl << endl;
        cout << " Available benchmarks:" << endl;
        cout << "  basic-alpha\t" << "\trequired parameters:" << endl;
        cout << "\t\t int:num_outputs, int:depth, int:num_chains, int:num_burnup, int:num_mcmc, int:gpuID" << endl;
        cout << "\t\t\tgpuID=-1 implies  using BLAS or no acceleration" << endl;
        cout << "  \t\t\toptional parameters:" << endl;
        cout << "\t\t string:outfilename" << endl << endl;
        #ifdef MPI_VERSION
        cout << "  basic-alpha-mpi" << "\trequired parameters:" << endl;
        cout << "\t\t int:num_outputs, int:depth, int:num_chains, int:num_burnup, int:num_mcmc" << endl;
        cout << "  \t\t\toptional parameters:" << endl;
        cout << "\t\t string:outfilename" << endl << endl;
        #else
        cout << "  basic-alpha-mpi" << "\tRequires build with mpi support" << endl;
        #endif

    }
}
