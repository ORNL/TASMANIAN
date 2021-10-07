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

    std::deque<std::string> args = stringArgs(argc, argv);

    if (args.empty()){
        cerr << "ERROR: no command specified\n\n";
        printHelp();
        return 1;
    }

    // basic help
    if (hasHelp(args.front())){
        printHelp();
        return 0;
    }

    // print log files
    if(args.front() == "-log"){
        show_log();
        return 0;
    }
    if(args.front() == "-cmakelog"){
        show_cmake_log();
        return 0;
    }

    // basic info, i.e., version, license, parallel support
    if (hasInfo(args.front())){
        cout << "Tasmanian Sparse Grids  version: " << TasmanianSparseGrid::getVersion() << endl;
        if ((std::string(TasmanianSparseGrid::getGitCommitHash()).compare("Tasmanian git hash is not available here") != 0)
            && (std::string(TasmanianSparseGrid::getGitCommitHash()).find("Release") != 0)){
            cout << "                git commit hash: " << TasmanianSparseGrid::getGitCommitHash() << endl;
            cout << "                cmake cxx flags: " << TasmanianSparseGrid::getCmakeCxxFlags() << endl;
        }
        cout << "                        license: " << TasmanianSparseGrid::getLicense() << endl;
        if (TasmanianSparseGrid::isOpenMPEnabled()){
            cout << "          OpenMP multithreading: Enabled" << endl;
        }else{
            cout << "          OpenMP multithreading: Disabled" << endl;
        }
        std::string gpu_backend = "none";
        if (TasmanianSparseGrid::isCudaEnabled()) gpu_backend = "CUDA";
        if (TasmanianSparseGrid::isHipEnabled()) gpu_backend = "ROCm/HIP";
        if (TasmanianSparseGrid::isDpcppEnabled()) gpu_backend = "oneAPI/DPC++";
        cout << "          GPU backend framework: " << gpu_backend << endl;
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
        if (TasmanianSparseGrid::isAccelerationAvailable(accel_gpu_magma)){
            cout << " " << AccelerationMeta::getIOAccelerationString(accel_gpu_magma);
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
                    std::string name = TasmanianSparseGrid::getGPUName(i);
                    int memory = TasmanianSparseGrid::getGPUMemory(i);
                    cout << setw(11) << i << ":" << setw(20) << name << " with" << setw(7) << memory << "MB of RAM" << endl;
                }
            }else{
                cout << "        Available GPUs: none" << endl;
            }
        }

        cout << endl;
        return 0;
    }

    // help with interface commands
    if (args.front() == "-listtypes"){
        printHelp(help_listtypes);
        return 0;
    }

    // testing
    if (args.front() == "-test"){
        bool debug = false;
        bool debugII = false;
        bool verbose = false;
        bool seed_reset = false;

        TestList test = test_all;

        int gpuid = -1;
        args.pop_front();
        while (!args.empty()){
            if (args.front() == "debug") debug = true;
            if (args.front() == "db") debugII = true;
            if (hasInfo(args.front())) verbose = true;
            if (hasRandom(args.front())) seed_reset = true;
            TestList test_maybe = ExternalTester::hasTest(args.front());
            if (test_maybe != test_none) test = test_maybe;
            if ((args.front() == "-gpuid") || (args.front() == "-gpu")){
                args.pop_front();
                if (args.empty()){
                    cerr << "ERROR: -gpuid required a valid number!" << endl;
                    return 1;
                }
                gpuid = std::stoi(args.front());
                if ((gpuid < -1) || (gpuid >= TasmanianSparseGrid::getNumGPUs())){
                    cerr << "ERROR: -gpuid " << gpuid << " is not a valid gpuid!" << endl;
                    cerr << "      see ./tasgrid -v for a list of detected GPUs." << endl;
                    return 1;
                }
            }
            args.pop_front();
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
    auto command = TasgridWrapper::hasCommand(args.front());
    if (command == command_none){
        cout << "ERROR: unknown command " << args.front() << endl;
        printHelp();
        return 1;
    }
    wrap.setCommand(command);

    // parse the parameters
    args.pop_front();
    while(!args.empty()){
        if (hasHelp(args.front())){
            printHelp(help_command, command);
            return 0;
        }else if (args.front() == "-xf" || args.front() == "-xfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide file name with x values!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setXFilename(args.front());
        }else if (args.front() == "-vf" || args.front() == "-valsfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide values file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setValsFilename(args.front());
        }else if (args.front() == "-of" || args.front() == "-outputfile" || args.front() == "-outfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide output file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setOutFilename(args.front());
        }else if (args.front() == "-gf" || args.front() == "-gridfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide grid file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setGridFilename(args.front());
        }else if (args.front() =="-ascii"){
            wrap.setUseASCII(true);
        }else if (args.front() == "-af" || args.front() == "-anisotropyfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide anisotropy file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setAnisoFilename(args.front());
        }else if (args.front() == "-tf" || args.front() == "-transformfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide transform file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setTransformFilename(args.front());
        }else if (args.front() == "-conformalfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide conformal transform file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setConformalFilename(args.front());
        }else if (args.front() == "-lf" || args.front() == "-levellimitsfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide level limits file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setLevelLimitsFilename(args.front());
        }else if (args.front() == "-cf" || args.front() == "-customfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide custom file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setCustomFilename(args.front());
        }else if (args.front() == "-print" || args.front() == "-p"){
            wrap.setPrintPoints(true);
        }else if (args.front() == "-dim" || args.front() == "-dimensions"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide number of dimensions!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setNumDimensions(std::stoi(args.front()));
        }else if (args.front() == "-out" || args.front() == "-outputs"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide number of outputs!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setNumOutputs(std::stoi(args.front()));
        }else if (args.front() == "-dt" || args.front() == "-depth"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid depth!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setNumDepth(std::stoi(args.front()));
        }else if (args.front() == "-or" || args.front() == "-order"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid order!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setOrder(std::stoi(args.front()));
        }else if (args.front() == "-tt" || args.front() == "-type"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid depth type!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            TypeDepth depth_type = IO::getDepthTypeString(args.front());
            if (depth_type == type_none){
                cerr << "ERROR: " << args.front() << " is not a valid type!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setDepthType(depth_type);
        }else if (args.front() == "-1d" || args.front() == "-onedim"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -onedim!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            TypeOneDRule rule = IO::getRuleString(args.front());
            if (rule == rule_none){
                cerr << "ERROR: " << args.front() << " is not a valid rule!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setRule(rule);
        }else if (args.front() == "-ct" || args.front() == "-conformaltype"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -conformaltype!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            TypeConformalMap conformal_type = TasgridWrapper::getConfromalType(args.front().c_str());
            if (conformal_type == conformal_none){
                cerr << "ERROR: " << args.front() << " is not a valid conformal type!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setConformalType(conformal_type);
        }else if (args.front() == "-alpha"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -alpha!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setAlpha(std::stof(args.front()));
        }else if (args.front() == "-beta"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -beta!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setBeta(std::stof(args.front()));
        }else if (args.front() == "-tol" || args.front() == "-tolerance"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -tolerance!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setTolerance(std::stof(args.front()));
        }else if (args.front() == "-rout" || args.front() == "-refout"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -refout!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setRefOutput(std::stoi(args.front()));
        }else if (args.front() == "-ming" || args.front() == "-mingrowth"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -mingrowth!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setMinGrowth(std::stoi(args.front()));
        }else if (args.front() == "-rt" || args.front() == "-reftype"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid depth type!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            TypeRefinement ref  = IO::getTypeRefinementString(args.front());
            if (ref == refine_none){
                cerr << "ERROR: " << args.front() << " is not a valid refinement type!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setTypeRefinement(ref);
        }else if (args.front() == "-gpu" || args.front() == "-gpuid"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -gpuid  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setGPID(std::stoi(args.front()));
        }else if (args.front() == "-shift"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -shift!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setShift(std::stof(args.front()));
        }else if (args.front() == "-wf" || args.front() == "-weightfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide weight file name!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setWeightFilename(args.front());
        }else if (args.front() == "-desc" || args.front() == "-description"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide description!!!  For help see: ./tasgrid -help" << endl << endl;
                return 1;
            }
            wrap.setDescription(args.front());
        }else if (args.front() == "-symm" || args.front() == "-symmetric"){
            wrap.setIsSymmetric(true);
        }else{
            cout << "WARNING: ignoring unknown option: " << args.front() << "\n";
        }
        args.pop_front();
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
        cout << " -makefourier\t"      << "\t-mf"     << "\t\tmake a grid from a Fourier rule" << endl;
        cout << " -makequadrature"     << "\t-mq"     << "\t\tmake a quadrature" << endl;
        cout << " -makeexoquad\t"     << "\t-meq"     << "\t\tmake an exotic quadrature" << endl;
        cout << " -makeupdate\t"       << "\t-mu"     << "\t\tupdates an existing global/sequence/fourier grid" << endl;
        cout << " -setconformal\t"     << "\t-sc"     << "\t\tset conformal domain transform" << endl;
        cout << " -getquadrature\t"    << "\t-gq"     << "\t\toutput quadrature weights and points" << endl;
        cout << " -getinterweights"    << "\t-gi"     << "\t\toutput the interpolation weights" << endl;
        cout << " -getpoints\t"    << "\t-gp"     << "\t\toutput the points" << endl;
        cout << " -getneededpoints"    << "\t-gn"     << "\t\toutputs the points needing values to build an interpolant" << endl;
        cout << " -loadvalues\t"       << "\t-l"      << "\t\tload the values of the interpolated function" << endl;
        cout << " -evaluate\t"     << "\t-e"      << "\t\tevaluates the interpolant" << endl;
        cout << " -evalhierarchyd"     << "\t-ehd"      << "\t\tevaluates the hierarchical basis (dense output)" << endl;
        cout << " -evalhierarchys"     << "\t-ehs"      << "\t\tevaluates the hierarchical basis (sparse output)" << endl;
        cout << " -gethsupport\t"       << "\t-ghsup"     << "\t\tget the hierarchical support" << endl;
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
        cout << " -gpuid\t\t"    << "\t\t"     << "\t<int>\t"   << "\tset the gpu to use for evaluations" << endl;
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
        }else if (com == command_makefourier){
            cout << "Commands\t"     << "\tShorthand" << "\tAction" << endl;
            cout << " -makefourier\t"     << "\t-mf"       << "\t\tmake a grid from a Fourier rule" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions" << endl;
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs" << endl;
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)" << endl;
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of the grid" << endl;
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
        }else if (com == command_makeexoquad){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -makeexoquad"   << "\t\t-meq"     << "\t\tmake an exotic quadrature" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -depth\t\t"    << "\tyes\t"      << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)" << endl;
            cout << " -shift\t\t"    << "\tyes\t"      << "\t<float>"    << "\t\tset the shift of the weight function" << endl;
            cout << " -weightfile\t"    << "\tyes\t"      << "\t<filename>"    << "\tset the name of the file containing a" << endl;
            cout << " \t\t\t\t\t"    << "\t\tsurrogate/interpolant of the weight function;" << endl;
            cout << " \t\t\t\t\t"    << "\t\tmust be a TasmanianSparseGrid in ASCII format" << endl;
            cout << " -description\t"    << "\tyes\t"      << "\t<string>"    << "\tshort description of the quadrature" << endl;
            cout << " -symmetric\t"    << "\tno\t"      << "\t<none>"    << "\t\tdeclare that the weight function is symmetric"  << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
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
        }else if (com == command_gethsupport){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction" << endl;
            cout << " -gethsupport\t"       << "\t-ghsup"      << "\t\tget the hierarchical support" << endl << endl;
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction" << endl;
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file" << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file" << endl;
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output" << endl << endl;
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
        cout << "      localp    localp-zero    semi-localp   localp-boundary\n\n";
        cout << "List of local wavelet grids rules:" << endl;
        cout << "      wavelet" << endl << endl;
        cout << "List of local polynomial and wavelet refinement types:" << endl;
        cout << "      classic     parents    direction    fds" << endl << endl;
        cout << "List of conformal maps:" << endl;
        cout << "      asin" << endl << endl;
    }
}
