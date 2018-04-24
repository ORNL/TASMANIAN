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
#include "tasgridWrapper.hpp"

using namespace std;
using namespace TasGrid;

enum TypeHelp{
    help_generic,
    help_command,
    help_listtypes
};

void printHelp(TypeHelp ht = help_generic, TypeCommand com = command_none);

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

    // basic info, i.e., version, license, parallel support
    if ((strcmp(argv[1],"-version") == 0) || (strcmp(argv[1],"-v") == 0) || (strcmp(argv[1],"-info") == 0)){
        cout << "Tasmanian Sparse Grids  version: " << TasmanianSparseGrid::getVersion() << endl;
        cout << "                        license: " << TasmanianSparseGrid::getLicense() << endl;
        if (TasmanianSparseGrid::isOpenMPEnabled()){
            cout << "          OpenMP multithreading: Enabled" << endl;
        }else{
            cout << "          OpenMP multithreading: Disabled" << endl;
        }
        cout << "         Available acceleration: ";
        bool anyAcc = false, anyGPU = false;
        if (TasmanianSparseGrid::isAccelerationAvailable(accel_cpu_blas)){
            cout << AccelerationMeta::getIOAccelerationString(accel_cpu_blas);
            anyAcc = true;
        }
        if (TasmanianSparseGrid::isAccelerationAvailable(accel_gpu_cublas)){
            cout << " " << AccelerationMeta::getIOAccelerationString(accel_gpu_cublas);
            anyAcc = true;
            anyGPU = true;
        }
        if (TasmanianSparseGrid::isAccelerationAvailable(accel_gpu_cuda)){
            cout << " " << AccelerationMeta::getIOAccelerationString(accel_gpu_cuda);
            anyAcc = true;
            anyGPU = true;
        }
        if (!anyAcc){
            cout << " none";
        }
        cout << endl;
        if (anyGPU){
            int numGPUs = TasmanianSparseGrid::getNumGPUs();
            if (numGPUs > 0){
                cout << "                 Available GPUs:" << endl;
                for(int i=0; i<numGPUs; i++){
                    char *name = TasmanianSparseGrid::getGPUName(i);
                    int memory = TasmanianSparseGrid::getGPUMemory(i);
                    cout << setw(11) << i << ":" << setw(20) << name << " with" << setw(7) << memory << "MB of RAM" << endl;
                    delete[] name;
                }
            }else{
                cout << "        Available GPUs: none" << endl;
            }
        }

        cout << endl;
        return 0;
    }

    // tests for speed
    if ((strcmp(argv[1],"-benchmark") == 0) || (strcmp(argv[1],"-bench") == 0)){
        ExternalTester tester(1);
        tester.resetRandomSeed();
        tester.benchmark(argc, argv);
        return 0;
    }


    // help with interface commands
    if (strcmp(argv[1],"-listtypes") == 0){
        printHelp(help_listtypes);
        return 0;
    }

    // testing
    if (strcmp(argv[1],"-test") == 0){
        bool debug = false;
        bool debugII = false;
        bool verbose = false;
        bool seed_reset = false;

        TestList test = test_all;

        int gpuid = -1;
        int k = 2;
        while (k < argc){
            if ((strcmp(argv[k],"debug") == 0)) debug = true;
            if ((strcmp(argv[k],"db") == 0)) debugII = true;
            if ((strcmp(argv[k],"verbose") == 0)) verbose = true;
            if ((strcmp(argv[k],"v") == 0)) verbose = true;
            if ((strcmp(argv[k],"-v") == 0)) verbose = true;
            if ((strcmp(argv[k],"-random") == 0)) seed_reset = true;
            if ((strcmp(argv[k],"random") == 0)) seed_reset = true;
            if ((strcmp(argv[k],"acceleration") == 0)) test = test_acceleration;
            if ((strcmp(argv[k],"domain") == 0)) test = test_domain;
            if ((strcmp(argv[k],"refinement") == 0)) test = test_refinement;
            if ((strcmp(argv[k],"global") == 0)) test = test_global;
            if ((strcmp(argv[k],"local") == 0)) test = test_local;
            if ((strcmp(argv[k],"wavelet") == 0)) test = test_wavelet;
            if ((strcmp(argv[k],"-gpuid") == 0)){
                if (k+1 >= argc){
                    cerr << "ERROR: -gpuid required a valid number!" << endl;
                    return 1;
                }
                gpuid = atoi(argv[k+1]);
                k++;
                if ((gpuid < -1) || (gpuid >= TasmanianSparseGrid::getNumGPUs())){
                    cerr << "ERROR: -gpuid " << gpuid << " is not a valid gpuid!" << endl;
                    cerr << "      see ./tasgrid -v for a list of detected GPUs." << endl;
                    return 1;
                }
            }
            k++;
        }

        ExternalTester tester(1000);
        tester.setGPUID(gpuid);
        bool pass = true;
        if (debug){
            tester.debugTest();
        }else if (debugII){
            tester.debugTestII();
        }else{
            if (verbose) tester.setVerbose(true);
            if (seed_reset) tester.resetRandomSeed();
            pass = tester.Test(test);
        }
        return (pass) ? 0 : 1;
    }

    // doing actual work with a grid
    // get the command
    TasgridWrapper wrap;
    bool deprecatedMakeGrid = false;
    if (strcmp(argv[1],"-makegrid") == 0){
        cerr << "WARNING: -makegrid is a deprecated command, in the future use -makeglobal, -makesequence, -makelocalpoly or -makewavelet" << endl;
        deprecatedMakeGrid = true; // use the correct command by guessing the grid type form the 1-D rule, wait for the 1-D rule
    }else if ((strcmp(argv[1],"-makeglobal") == 0) || (strcmp(argv[1],"-mg") == 0)){
        wrap.setCommand(command_makeglobal);
    }else if ((strcmp(argv[1],"-makesequence") == 0) || (strcmp(argv[1],"-ms") == 0)){
        wrap.setCommand(command_makesequence);
    }else if ((strcmp(argv[1],"-makelocalpoly") == 0) || (strcmp(argv[1],"-mp") == 0)){
        wrap.setCommand(command_makelocalp);
    }else if ((strcmp(argv[1],"-makewavelet") == 0) || (strcmp(argv[1],"-mw") == 0)){
        wrap.setCommand(command_makewavelet);
    }else if ((strcmp(argv[1],"-makequadrature") == 0) || (strcmp(argv[1],"-mq") == 0)){
        wrap.setCommand(command_makequadrature);
    }else if ((strcmp(argv[1],"-makeupdate") == 0) || (strcmp(argv[1],"-mu") == 0)){
        wrap.setCommand(command_update);
    }else if ((strcmp(argv[1],"-setconformal") == 0) || (strcmp(argv[1],"-sc") == 0)){
        wrap.setCommand(command_setconformal);
    }else if ((strcmp(argv[1],"-getquadrature") == 0) || (strcmp(argv[1],"-gq") == 0)){
        wrap.setCommand(command_getquadrature);
    }else if ((strcmp(argv[1],"-getinterweights") == 0) || (strcmp(argv[1],"-gi") == 0)){
        wrap.setCommand(command_getinterweights);
    }else if ((strcmp(argv[1],"-getpoints") == 0) || (strcmp(argv[1],"-gp") == 0)){
        wrap.setCommand(command_getpoints);
    }else if ((strcmp(argv[1],"-getneededpoints") == 0) || (strcmp(argv[1],"-gn") == 0) || (strcmp(argv[1],"-getneeded") == 0)){
        if (strcmp(argv[1],"-getneededpoints") == 0){
            cerr << "WARNING: -getneededpoints is a deprecated command, use -getneeded instead" << endl;
        }
        wrap.setCommand(command_getneeded);
    }else if ((strcmp(argv[1],"-loadvalues") == 0) || (strcmp(argv[1],"-l") == 0)){
        wrap.setCommand(command_loadvalues);
    }else if ((strcmp(argv[1],"-evaluate") == 0) || (strcmp(argv[1],"-e") == 0)){
        wrap.setCommand(command_evaluate);
    }else if ((strcmp(argv[1],"-evalhierarchyd") == 0) || (strcmp(argv[1],"-ehd") == 0)){
        wrap.setCommand(command_evalhierarchical_dense);
    }else if ((strcmp(argv[1],"-evalhierarchys") == 0) || (strcmp(argv[1],"-ehs") == 0)){
        wrap.setCommand(command_evalhierarchical_sparse);
    }else if ((strcmp(argv[1],"-integrate") == 0) || (strcmp(argv[1],"-i") == 0)){
        wrap.setCommand(command_integrate);
    }else if ((strcmp(argv[1],"-getanisotropy") == 0) || (strcmp(argv[1],"-ga") == 0)){
        wrap.setCommand(command_getanisocoeff);
    }else if ((strcmp(argv[1],"-refinesurp") == 0) || (strcmp(argv[1],"-rs") == 0)){
        wrap.setCommand(command_refine_surp);
    }else if ((strcmp(argv[1],"-refineaniso") == 0) || (strcmp(argv[1],"-ra") == 0)){
        wrap.setCommand(command_refine_aniso);
    }else if ((strcmp(argv[1],"-refine") == 0) || (strcmp(argv[1],"-r") == 0)){
        wrap.setCommand(command_refine);
    }else if ((strcmp(argv[1],"-cancelrefine") == 0) || (strcmp(argv[1],"-cr") == 0)){
        wrap.setCommand(command_refine_clear);
    }else if ((strcmp(argv[1],"-mergerefine") == 0) || (strcmp(argv[1],"-mr") == 0)){
        wrap.setCommand(command_refine_merge);
    }else if (strcmp(argv[1],"-getpoly") == 0){
        wrap.setCommand(command_getpoly);
    }else if ((strcmp(argv[1],"-summary") == 0) || (strcmp(argv[1],"-s") == 0)){
        wrap.setCommand(command_summary);
    }else if ((strcmp(argv[1],"-getcoefficients") == 0) || (strcmp(argv[1],"-gc") == 0)) {
        wrap.setCommand(command_getcoefficients);
    }else if ((strcmp(argv[1],"-setcoefficients") == 0) || (strcmp(argv[1],"-sc") == 0)) {
        wrap.setCommand(command_setcoefficients);
    }else if (strcmp(argv[1],"-getpointsindexes") == 0){
        wrap.setCommand(command_getpointsindex);
    }else if (strcmp(argv[1],"-getneededindexes") == 0){
        wrap.setCommand(command_getneededindex);
    }else{
        cout << "ERROR: unknown command " << argv[1] << endl;
        printHelp();
        return 1;
    }

    // parse the parameters
    int k = 2;
    bool commandHelp = false;
    while(k < argc){
        if ((strcmp(argv[k],"--help") == 0)||(strcmp(argv[k],"-help") == 0)||(strcmp(argv[k],"-h") == 0)||(strcmp(argv[k],"help") == 0)){
            commandHelp = true;
            break;
        }else if ((strcmp(argv[k],"-if") == 0)||(strcmp(argv[k],"-inputfile") == 0)){
            if (k+1 < argc){
                if ((wrap.getCommand() == command_evaluate) || (wrap.getCommand() == command_getinterweights)){
                    wrap.setXFilename(argv[++k]);
                    cerr << "WARNING: -inputfile is a deprecated parameter, in future use -xfile" << endl;
                }else if (wrap.getCommand() == command_loadvalues){
                    wrap.setValsFilename(argv[++k]);
                    cerr << "WARNING: -inputfile is a deprecated parameter, in future use -valsfile" << endl;
                }else{
                    cerr << "WARNING: ignoring deprecated parameter -inputfile" << endl;
                }
            }else{
                cerr << "WARNING: -inputfile is a deprecated, see ./tasgrid -help for correct use in the future" << endl;
                cerr << "ERROR: must provide input filename!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-xf") == 0)||(strcmp(argv[k],"-xfile") == 0)){
            if (k+1 < argc){
                wrap.setXFilename(argv[++k]);
            }else{
                cerr << "ERROR: must provide file name with x values!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-vf") == 0)||(strcmp(argv[k],"-valsfile") == 0)){
            if (k+1 < argc){
                wrap.setValsFilename(argv[++k]);
            }else{
                cerr << "ERROR: must provide values file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-of") == 0)||(strcmp(argv[k],"-outfile") == 0)){
            if (k+1 < argc){
                wrap.setOutFilename(argv[++k]);
            }else{
                cerr << "ERROR: must provide output file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-gf") == 0)||(strcmp(argv[k],"-gridfile") == 0)){
            if (k+1 < argc){
                wrap.setGridFilename(argv[++k]);
            }else{
                cerr << "ERROR: must provide grid file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if (strcmp(argv[k],"-ascii") == 0){
            wrap.setUseASCII(true);
        }else if ((strcmp(argv[k],"-af") == 0)||(strcmp(argv[k],"-anisotropyfile") == 0)){
            if (k+1 < argc){
                wrap.setAnisoFilename(argv[++k]);
            }else{
                cerr << "ERROR: must provide anisotropy file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-tf") == 0)||(strcmp(argv[k],"-transformfile") == 0)){
            if (k+1 < argc){
                wrap.setTransformFilename(argv[++k]);
            }else{
                cerr << "ERROR: must provide transform file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if (strcmp(argv[k],"-conformalfile") == 0){
            if (k+1 < argc){
                wrap.setConformalFilename(argv[++k]);
            }else{
                cerr << "ERROR: must provide conformal transform file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-levellimitsfile") == 0) || (strcmp(argv[k],"-lf") == 0)){
            if (k+1 < argc){
                wrap.setLevelLimitsFilename(argv[++k]);
            }else{
                cerr << "ERROR: must provide conformal transform file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-cf") == 0)||(strcmp(argv[k],"-customfile") == 0)){
            if (k+1 < argc){
                wrap.setCustomFilename(argv[++k]);
            }else{
                cerr << "ERROR: must provide custom file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-p") == 0)||(strcmp(argv[k],"-print") == 0)){
            wrap.setPrintPoints(true);
        }else if ((strcmp(argv[k],"-dim") == 0)||(strcmp(argv[k],"-dimensions") == 0)){
            if (k+1 < argc){
                int n = atoi(argv[++k]);
                if (n < 1){
                    cerr << "ERROR: -dimensions takes a positive integer" << endl << endl;
                    return 1;
                }
                wrap.setNumDimensions(n);
            }else{
                cerr << "ERROR: must provide number of dimensions!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-out") == 0)||(strcmp(argv[k],"-outputs") == 0)){
            if (k+1 < argc){
                int n = atoi(argv[++k]);
                if (n < 0){
                    cerr << "ERROR: -outputs takes a non-negative integer" << endl << endl;
                    return 1;
                }
                wrap.setNumOutputs(n);
            }else{
                cerr << "ERROR: must provide number of outputs!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-dt") == 0)||(strcmp(argv[k],"-depth") == 0)){
            if (k+1 < argc){
                int n = atoi(argv[++k]);
                if (n < 0){
                    cerr << "ERROR: -depth takes a non-negative integer" << endl << endl;
                    return 1;
                }
                wrap.setNumDepth(n);
            }else{
                cerr << "ERROR: must provide valid depth!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-tt") == 0)||(strcmp(argv[k],"-type") == 0)){
            if (k+1 < argc){
                k++;
                TypeDepth d = OneDimensionalMeta::getIOTypeString(argv[k]);
                if (d == type_none){
                    cerr << "ERROR: " << argv[k] << " is not a valid type!!!  For help see: ./tasgrid -help" << endl << endl;
                    return 1;
                }
                wrap.setDepthType(d);
            }else{
                cerr << "ERROR: must provide valid -type!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-ct") == 0)||(strcmp(argv[k],"-conformaltype") == 0)){
            if (k+1 < argc){
                k++;
                TypeConformalMap c = TasgridWrapper::getConfromalType(argv[k]);
                if (c == conformal_none){
                    cerr << "ERROR: " << argv[k] << " is not a valid type!!!  For help see: ./tasgrid -help" << endl << endl;
                    return 1;
                }
                wrap.setConformalType(c);
            }else{
                cerr << "ERROR: must provide valid -conformaltype!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-onedim") == 0)||(strcmp(argv[k],"-1d") == 0)){
            if (k+1 < argc){
                k++;
                TypeOneDRule r = OneDimensionalMeta::getIORuleString(argv[k]);
                if (r == rule_none){
                    cerr << "ERROR: unrecognized rule: " << argv[k] << endl << endl;
                    return 1;
                }
                wrap.setRule(r);
            }else{
                cerr << "ERROR: must provide valid -onedim!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-order") == 0)||(strcmp(argv[k],"-or") == 0)){
            if (k+1 < argc){
                k++;
                int n = atoi(argv[k]);
                if (n < -1){
                    cerr << "ERROR: invalid order: " << n << "  order should be at least -1" << endl << endl;
                    return 1;
                }
                wrap.setOrder(n);
            }else{
                cerr << "ERROR: must provide valid -order!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if (strcmp(argv[k],"-alpha") == 0){
            if (k+1 < argc){
                k++;
                double alpha = atof(argv[k]);
                wrap.setAlpha(alpha);
            }else{
                cerr << "ERROR: must provide valid -alpha!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if (strcmp(argv[k],"-beta") == 0){
            if (k+1 < argc){
                k++;
                double beta = atof(argv[k]);
                wrap.setBeta(beta);
            }else{
                cerr << "ERROR: must provide valid -beta!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-tolerance") == 0)||(strcmp(argv[k],"-tol") == 0)){
            if (k+1 < argc){
                k++;
                double tol = atof(argv[k]);
                if (tol < 0.0){
                    cerr << "WARNING: tolerance set to: " << tol << "  which is a negative number, are you sure about this?" << endl;
                }
                wrap.setTolerance(tol);
            }else{
                cerr << "ERROR: must provide valid -tolerance!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-refout") == 0)||(strcmp(argv[k],"-rout") == 0)){
            if (k+1 < argc){
                k++;
                int rout = atoi(argv[k]);
                if (rout < 0){
                    cerr << "ERROR: the refinement output -refout/-rout must be non-negative!!!  For help see: ./tasgrid -help" << endl << endl;
                    return 1;
                }
                wrap.setRefOutput(rout);
            }else{
                cerr << "ERROR: must provide valid -refout!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-mingrowth") == 0)||(strcmp(argv[k],"-ming") == 0)){
            if (k+1 < argc){
                k++;
                int mg = atoi(argv[k]);
                if (mg < 1){
                    cerr << "ERROR: the minimum growth must be positive!!!  For help see: ./tasgrid -help" << endl << endl;
                    return 1;
                }
                wrap.setMinGrowth(mg);
            }else{
                cerr << "ERROR: must provide valid -mingrowth!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if ((strcmp(argv[k],"-reftype") == 0)||(strcmp(argv[k],"-rt") == 0)){
            if (k+1 < argc){
                k++;
                TypeRefinement r = OneDimensionalMeta::getIOTypeRefinementString(argv[k]);
                if (r == refine_none){
                    cerr << "ERROR: " << argv[k] << " is not a valid refinement type!!!  For help see: ./tasgrid -help" << endl << endl;
                    return 1;
                }
                wrap.setTypeRefinement(r);
            }else{
                cerr << "ERROR: must provide valid -reftype!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }else if (strcmp(argv[k],"-gpuid") == 0){
            if (k+1 < argc){
                k++;
                int g = atoi(argv[k]);
                wrap.setGPID(g);
            }else{
                cerr << "ERROR: must provide valid -gpuid  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
        }
        k++;
    }

    if (commandHelp){
        printHelp(help_command, wrap.getCommand());
        return 0;
    }

    if (deprecatedMakeGrid){
        if (OneDimensionalMeta::isGlobal(wrap.getRule())){
            wrap.setCommand(command_makeglobal);
        }else if (OneDimensionalMeta::isLocalPolynomial(wrap.getRule())){
            wrap.setCommand(command_makelocalp);
        }else if (OneDimensionalMeta::isWavelet(wrap.getRule())){
            wrap.setCommand(command_makewavelet);
        }else{
            cerr << "ERROR: deprecated -makegrid used but no -onedim specified" << endl;
            return 1;
        }
    }

    if (!wrap.executeCommand()){
        return 1;
    }

    return 0;
}


void printHelp(TypeHelp ht, TypeCommand com){

    cout << endl;

    if (ht == help_generic){

        cout << " Usage: tasgrid <command> <option1> <value1> <option2> <value2> ... " << endl << endl;

        cout << "Commands\t"       << "\tShorthand"   << "\tAction" << endl;
        cout << " -help\t"         << "\t\t-h,--help" << "\tshow verbose help options" << endl;
        cout << " -listtypes\t"    << "\t-lt\t"       << "\tlist accepted grid types and 1-D rules" << endl;
        cout << " -test\t"         << "\t\t\t"    << "\tperform a number of internal tests" << endl;
        cout << " -makeglobal\t"       << "\t-mg"     << "\t\tmake a grid from a global rule" << endl;
        cout << " -makesequence\t"     << "\t-ms"     << "\t\tmake a grid from a sequence rule" << endl;
        cout << " -makelocalpoly\t"    << "\t-mp"     << "\t\tmake a grid from a local polynomial rule" << endl;
        cout << " -makewavelet\t"      << "\t-mw"     << "\t\tmake a grid from a wavelet rule" << endl;
        cout << " -makequadrature"     << "\t-mq"     << "\t\tmake a quadrature" << endl;
        cout << " -makeupdate\t"       << "\t-mu"     << "\t\tupdates a new global or sequence grid" << endl;
        cout << " -setconformal\t"     << "\t-sc"     << "\t\tset conformal domain transform" << endl;
        cout << " -getquadrature\t"    << "\t-gq"     << "\t\toutput quadrature weights and points" << endl;
        cout << " -getinterweights"    << "\t-gi"     << "\t\toutput the interpolation weights" << endl;
        cout << " -getpoints\t"    << "\t-gp"     << "\t\toutput the points" << endl;
        cout << " -getneededpoints"    << "\t-gn"     << "\t\toutputs the points needing values to build an interpolant" << endl;
        cout << " -loadvalues\t"       << "\t-l"      << "\t\tload the values of the interpolated function" << endl;
        cout << " -evaluate\t"     << "\t-e"      << "\t\tevaluates the interpolant" << endl;
        cout << " -evalhierarchyd"     << "\t-ehd"      << "\t\tevaluates the hierarchical basis (dense output)" << endl;
        cout << " -evalhierarchys"     << "\t-ehs"      << "\t\tevaluates the hierarchical basis (sparse output)" << endl;
        cout << " -integrate\t"    << "\t-i"      << "\t\toutput the integral" << endl;
        cout << " -getanisotropy\t"    << "\t-ga"     << "\t\testimates the anisotropic coefficients" << endl;
        cout << " -refineaniso\t"      << "\t-ra"     << "\t\trefines the grid" << endl;
        cout << " -refinesurp\t"       << "\t-rs"     << "\t\trefines the grid" << endl;
        cout << " -refine\t"       << "\t-r"      << "\t\trefines the grid" << endl;
        cout << " -cancelrefine\t"     << "\t-cr"     << "\t\tdiscards the last refinement (unless values are already loaded)" << endl;
        cout << " -mergerefine\t"     << "\t-mr"     << "\t\tcombines the loaded and needed points and discards the loaded values" << endl;
        cout << " -getcoefficients"    << "\t-gc"     << "\t\tget the hierarchical coefficients of the grid" << endl;
        cout << " -setcoefficients"    << "\t-sc"     << "\t\tset the hierarchical coefficients of the grid" << endl;
        cout << " -getpoly\t"      << "\t"      << "\t\tget polynomial space" << endl;
        cout << " -summary\t"      << "\t-s"      << "\t\twrites short description" << endl << endl;

        cout << "Options\t\t"    << "\tShorthand"  << "\tValue"    << "\t\tAction" << endl;
        cout << " -help\t\t"     << "\thelp\t"     << "\t"         << "\t\tdisplay verbose information about this command" << endl;
        cout << " -dimensions\t"     << "\t-dim\t"     << "\t<int>"    << "\t\tset the number of dimensions" << endl;
        cout << " -outputs\t"    << "\t-out\t"     << "\t<int>"    << "\t\tset the number of outputs" << endl;
        cout << " -depth\t\t"    << "\t-dt\t"      << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)" << endl;
        cout << " -type\t\t"     << "\t-tt\t"      << "\t<type>"       << "\t\tset the type of the grid" << endl;
        cout << " -conformaltype\t"     << "\t-tt\t"      << "\t<type>"       << "\t\tset the type of the transformation" << endl;
        cout << " -onedim\t"     << "\t-1d\t"      << "\t<rule>"       << "\t\tset the one dimensional rule" << endl;
        cout << " -order\t\t"    << "\t-or\t"      << "\t<int>"    << "\t\tset the order for local polynomial and wavelet basis" << endl;
        cout << " -alpha\t\t"    << "\t\t"     << "\t<float>"      << "\t\tthe alpha parameter for Gegenbauer/Jacobi/Laguerre/Hermite quadrature" << endl;
        cout << " -beta\t\t"     << "\t\t"     << "\t<float>"      << "\t\tthe beta parameter for Jacobi quadrature" << endl;
        cout << " -tolerance\t"      << "\t-tol\t"     << "\t<float>"      << "\t\tset the tolerance for the refinement" << endl;
        cout << " -refout\t"     << "\t-rout\t"    << "\t<int>"    << "\t\tselect the output to use for the refinement" << endl;
        cout << " -mingrowth\t"      << "\t-ming\t"   << "\t<int>"    << "\t\tminimum number of new points" << endl;
        cout << " -reftype\t"    << "\t-rt\t"      << "\t<int>"    << "\t\tset the type of refinement" << endl;

        cout << " -gridfile\t"       << "\t-gf\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
        cout << " -xfile\t\t"    << "\t-xf\t"      << "\t<filename>"   << "\tset the name for the file with points" << endl;
        cout << " -valsfile\t"       << "\t-vf\t"      << "\t<filename>"   << "\tset the name for the file with values" << endl;
        cout << " -outputfile\t"     << "\t-of\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
        cout << " -anisotropyfile"   << "\t-af\t"      << "\t<filename>"   << "\tset the anisotropic weights" << endl;
        cout << " -transformfile\t"  << "\t-tf\t"      << "\t<filename>"   << "\tset the transformation of the domain" << endl;
        cout << " -conformalfile\t"  << "\t-tf\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain" << endl;
        cout << " -levellimitsfile"  << "\t-lf\t"      << "\t<filename>"   << "\tset the limits for the levels" << endl;
        cout << " -customfile\t"     << "\t-cf\t"      << "\t<filename>"   << "\tset the file with the custom-tabulated rule" << endl;
        cout << " -print\t\t"    << "\t-p\t"       << "\t<none>"       << "\t\tprint to standard output" << endl;
        cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format" << endl;

        cout << endl;
    }else if (ht == help_command){
        if (com == command_makeglobal){
            cout << "Commands\t"     << "\tShorthand" << "\tAction" << endl;
            cout << " -makeglobal\t"     << "\t-mg"       << "\t\tmake a grid from a global rule" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions" << endl;
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs" << endl;
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)" << endl;
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of the grid" << endl;
            cout << " -onedim\t"     << "\tyes\t"     << "\t<rule>"       << "\t\tset the one dimensional rule" << endl;
            cout << " \t\t\t\t\t"    << "\t\tmust use a global rule (see manual)" << endl;
            cout << " -alpha\t\t"    << "\tsometimes" << "\t<float>"      << "\t\tthe alpha parameter for Gegenbauer/Jacobi/Laguerre/Hermite quadrature" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired for those rules" << endl;
            cout << " -beta\t\t"     << "\tsometimes" << "\t<float>"      << "\t\tthe beta parameter for Jacobi quadrature" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired for Jacobi rule" << endl;
            cout << " -gridfile\t"       << "\tno\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -anisotropyfile"   << "\tno\t"      << "\t<filename>"   << "\tset the anisotropic weights" << endl;
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain" << endl;
            cout << " -customfile\t"     << "\tsometimes" << "\t<filename>"   << "\tset the file with the custom-tabulated rule" << endl;
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map" << endl;
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired for custom-tabulated rule" << endl;
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels" << endl;
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output" << endl;
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format" << endl << endl;
            cout << "Note: -outputfile or -print output the points of the grid" << endl;
            cout << "Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_makesequence){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -makesequence\t"   << "\t-ms"     << "\t\tmake a grid from a sequence rule" << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions" << endl << endl;
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs" << endl;
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)" << endl;
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of the grid" << endl;
            cout << " -onedim\t"     << "\tyes\t"     << "\t<rule>"       << "\t\tset the one dimensional rule" << endl;
            cout << " \t\t\t\t\t"    << "\t\tmust use a sequence rule (see manual)" << endl;
            cout << " -gridfile\t"       << "\tno\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -anisotropyfile"   << "\tno\t"      << "\t<filename>"   << "\tset the anisotropic weights" << endl;
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain" << endl;
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map" << endl;
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain" << endl;
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels" << endl;
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output" << endl;
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format" << endl << endl;
            cout << "Note: -outputfile or -print output the points of the grid" << endl;
            cout << "Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_makelocalp){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -makelocalpoly\t"  << "\t-mp"     << "\t\tmake a grid from a local polynomial rule" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions" << endl;
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs" << endl;
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)" << endl;
            cout << " -order\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the order for local polynomial basis" << endl;
            cout << " -onedim\t"     << "\tyes\t"     << "\t<rule>"       << "\t\tset the one dimensional rule" << endl;
            cout << " \t\t\t\t\t"    << "\t\tmust use a local polynomial rule (see manual)" << endl;
            cout << " -gridfile\t"       << "\tno\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain" << endl;
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map" << endl;
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain" << endl;
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels" << endl;
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output" << endl;
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format" << endl << endl;
            cout << "Note: -outputfile or -print output the points of the grid" << endl;
            cout << "Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_makewavelet){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -makewavelet\t"    << "\t-mw"     << "\t\tmake a grid from a wavelet rule" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions" << endl;
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs" << endl;
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)" << endl;
            cout << " -order\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the order for the wavelet basis" << endl;
            cout << " -gridfile\t"       << "\tno\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain" << endl;
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map" << endl;
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain" << endl;
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels" << endl;
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output" << endl;
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format" << endl << endl;
            cout << "Note: -outputfile or -print output the points of the grid" << endl;
            cout << "Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_makequadrature){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -makequadrature"   << "\t-mq"     << "\t\tmake a quadrature" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions" << endl;
            cout << " -depth\t\t"    << "\tyes\t"      << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)" << endl;
            cout << " -type\t\t"     << "\tsometimes"      << "\t<type>"       << "\t\tset the type of the grid" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired only for global rules (see manual)" << endl;
            cout << " -order\t\t"    << "\tsometimes"      << "\t<int>"    << "\t\tset the order for local polynomial basis" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired only for local polynomial and wavelet (see manual)" << endl;
            cout << " -onedim\t"     << "\tyes\t"      << "\t<rule>"       << "\t\tset the one dimensional rule" << endl;
            cout << " \t\t\t\t\t"    << "\t\tcan use any rule (see manual)" << endl;
            cout << " -alpha\t\t"    << "\tsometimes"  << "\t<float>"      << "\t\tthe alpha parameter for Gegenbauer/Jacobi/Laguerre/Hermite quadrature" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired for those rules" << endl;
            cout << " -beta\t\t"     << "\tsometimes"  << "\t<float>"      << "\t\tthe beta parameter for Jacobi quadrature" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired for Jacobi rule" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -anisotropyfile"   << "\tno\t"      << "\t<filename>"   << "\tset the anisotropic weights" << endl;
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain" << endl;
            cout << " -customfile\t"     << "\tsometimes"  << "\t<filename>"   << "\tset the file with the custom-tabulated rule" << endl;
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map" << endl;
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired for custom-tabulated rule" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the points and weights of the grid" << endl;
            cout << "Note: at least one of -outputfile, or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_update){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -makeupdate\t"     << "\t-mu"     << "\t\tupdates a new global or sequence grid" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -depth\t\t"    << "\tyes\t"      << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)" << endl;
            cout << " -type\t\t"     << "\tyes\t"      << "\t<type>"       << "\t\tset the type of the grid" << endl;
            cout << " -anisotropyfile"   << "\tno\t"      << "\t<filename>"   << "\tset the anisotropic weights" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the new points of the grid" << endl;
        }else if (com == command_setconformal){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -setconformal\t"     << "\t-sc"     << "\t\tset conformal transformation" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -conformaltype\t"  << "\tyes\t"      << "\t<type>"       << "\t\tset the type of the map" << endl;
            cout << " -conformalfile\t"  << "\tyes\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain" << endl;
        }else if (com == command_getquadrature){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -getquadrature"    << "\t-gq"     << "\t\tmake a quadrature" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the points and weights of the grid" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_getcoefficients){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -getcoefficients"    << "\t-gc"     << "\t\tget the hierarchical coefficients" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the hierarchical coefficients of the sparse grid" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_setcoefficients){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -setcoefficients"    << "\t-sc"     << "\t\tset the hierarchical coefficients" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the hierarchical coefficients of the sparse grid" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_getinterweights){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -getinterweights"  << "\t-gi"     << "\t\toutput the interpolation weights" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the interpolation weight for each point in the xfile, see equation (1.2) in the manual" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_getpoints){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -getpoints\t"      << "\t-gp"     << "\t\toutput the points" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the points of the grid" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_getneeded){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -getneededpoints"  << "\t-gn"     << "\t\toutputs the points needing values to build an interpolant" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the new points of the grid" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_loadvalues){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -loadvalues\t"     << "\t-l"      << "\t\tload the values of the interpolated function" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -valsfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the file with values" << endl << endl;
        }else if (com == command_evaluate){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -evaluate\t"       << "\t-e"      << "\t\tevaluates the interpolant" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points" << endl;
            cout << " -gpuid\t\t"    << "\tno\t"     << "\t<int>\t"   << "\tset the gpu to use for evaluations" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output values of the interpolant at the points specified in the xfile" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_evalhierarchical_dense){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -evalhierarchyd"       << "\t-ehd"      << "\t\tevaluates the hierarchical basis (dense output)" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output values of the hierarchical basis functions in the xfile" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_evalhierarchical_sparse){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -evalhierarchys"       << "\t-ehs"      << "\t\tevaluates the hierarchical basis (sparse output)" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output values of the hierarchical basis functions in the xfile" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_integrate){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -integrate\t"      << "\t-i"      << "\t\toutput the integral" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the integral if the loaded function, see equation (1.3) in the manual" << endl;
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output" << endl << endl;
        }else if (com == command_getanisocoeff){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -getanisotropy\t"    << "\t-ga"     << "\t\testimates the anisotropic coefficients" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of anisotropic coefficients" << endl;
            cout << " -refout\t"     << "\tsometimes" << "\t<int>"    << "\t\tselect the output to use" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired by global grids only" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the estimated anisotropic coefficients" << endl << endl;
        }else if (com == command_refine_aniso){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -refineaniso\t"    << "\t-ra"     << "\t\trefines the grid" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of anisotropic refinement" << endl;
            cout << " -mingrowth\t"      << "\tno\t"      << "\t<int>"    << "\t\tminimum number of new points (defaults to 1)" << endl;
            cout << " -refout\t"     << "\tsometimes" << "\t<int>"    << "\t\tselect the output to use for the refinement" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired by global grids, for sequence grids defaults to -1 (use all outputs)" << endl;
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl;
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format" << endl << endl;
            cout << "Note: -outputfile or -print output the new points of the grid" << endl << endl;
        }else if (com == command_refine_surp){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -refinesurp\t"     << "\t-rs"     << "\t\trefines the grid" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -tolerance\t"      << "\tyes\t"     << "\t<float>"      << "\t\tset the tolerance for the refinement" << endl;
            cout << " -reftype\t"    << "\tsometimes" << "\t<int>"    << "\t\tset the type of refinement" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired by local polynomial and wavelet grids" << endl;
            cout << " -refout\t"     << "\tsometimes" << "\t<int>"    << "\t\tselect the output to use for the refinement" << endl;
            cout << " \t\t\t\t\t"    << "\t\trequired by global grids, for sequence grids defaults to -1 (use all outputs)" << endl;
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels" << endl;
            cout << " -valsfile\t"       << "\t-vf\t"      << "\t<filename>"   << "\tset the correction weights for the surpluses" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl;
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format" << endl << endl;
            cout << "Note: -outputfile or -print output the new points of the grid" << endl << endl;
        }else if (com == command_refine){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -refine\t"     << "\t-r"    << "\t\trefines the grid" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << " -refine calls -refineaniso for Global and Sequence grids and -refinesurp otherwise" << endl;
            cout << " see \"-refineaniso help\" or \"-refinesurp help\" for the corresponding accepted options" << endl << endl;
            cout << "Note: -outputfile or -print output the new points of the grid" << endl << endl;
        }else if (com == command_refine_clear){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -cancelrefine\t"   << "\t-cr"     << "\t\tdiscards the last refinement (unless values are already loaded)" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format" << endl << endl;
        }else if (com == command_refine_merge){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -mergerefine\t"   << "\t-mr"     << "\t\tmerges the loaded and needed points and discards any loaded values" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format" << endl << endl;
        }else if (com == command_getpoly){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -getpoly\t"      << "\t\t"      << "\t\tget polynomial space" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tspecifies whether to use quadrature or interpolation" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
            cout << "Note: -outputfile or -print output the polynomial indexes" << endl << endl;
        }else if (com == command_summary){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -summary\t"    << "\t-s"      << "\t\twrites short description" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file" << endl << endl;
        }
    }else if (ht == help_listtypes){
        cout << "This only lists the strings associated with each option for spelling purposes." << endl;
        cout << "Refer to the manual for details about each type." << endl << endl;
        cout << "List of global grids types:" << endl;
        cout << "  level     curved     hyperbolic     tensor" << endl;
        cout << "  iptotal   ipcurved   iphyperbolic   iptensor" << endl;
        cout << "  qptotal   qpcurved   qphyperbolic   qptensor" << endl << endl;
        cout << "List of global grids rules:" << endl;
        cout << "         chebyshev          chebyshev-odd    clenshaw-curtis   clenshaw-curtis-zero" << endl;
        cout << "              leja               leja-odd              rleja              rleja-odd" << endl;
        cout << "     rleja-double2          rleja-double4      rleja-shifted     rleja-shifted-even" << endl;
        cout << "      max-lebesgue       max-lebesgue-odd       min-lebesgue       min-lebesgue-odd" << endl;
        cout << "         min-delta          min-delta-odd             fejer2" << endl;
        cout << "    gauss-legendre     gauss-legendre-odd    gauss-patterson       custom-tabulated" << endl;
        cout << "  gauss-gegenbauer   gauss-gegenbauer-odd       gauss-jacobi       gauss-jacobi-odd" << endl;
        cout << "    gauss-laguerre     gauss-laguerre-odd      gauss-hermite      gauss-hermite-odd" << endl;
        cout << "  gauss-chebyshev1   gauss-chebyshev1-odd   gauss-chebyshev2   gauss-chebyshev2-odd" << endl << endl;
        cout << "List of sequence grids rules:" << endl;
        cout << "              leja              rleja     rleja-shifted" << endl;
        cout << "      max-lebesgue       min-lebesgue         min-delta " << endl << endl;
        cout << "List of local polynomial grids rules:" << endl;
        cout << "      localp    localp-zero    semi-localp" << endl << endl;
        cout << "List of local wavelet grids rules:" << endl;
        cout << "      wavelet" << endl << endl;
        cout << "List of local polynomial and wavelet refinement types:" << endl;
        cout << "      classic     parents    direction    fds" << endl << endl;
        cout << "List of conformal maps:" << endl;
        cout << "      asin" << endl << endl;
    }
}
