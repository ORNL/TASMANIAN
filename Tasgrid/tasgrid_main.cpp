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

#include "gridtestExternalTests.hpp"
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

    //cout << " Phruuuuphrrr \n"; // this is the sound that the Tasmanian devil makes

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
            cout << "          OpenMP multithreading: Enabled\n";
        }else{
            cout << "          OpenMP multithreading: Disabled\n";
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
                cout << "        Available GPUs: none\n";
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
                    cerr << "ERROR: -gpuid required a valid number!\n";
                    return 1;
                }
                gpuid = std::stoi(args.front());
                if ((gpuid < -1) || (gpuid >= TasmanianSparseGrid::getNumGPUs())){
                    cerr << "ERROR: -gpuid " << gpuid << " is not a valid gpuid!\n";
                    cerr << "      see ./tasgrid -v for a list of detected GPUs.\n";
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
                cerr << "ERROR: must provide file name with x values!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setXFilename(args.front());
        }else if (args.front() == "-vf" || args.front() == "-valsfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide values file name!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setValsFilename(args.front());
        }else if (args.front() == "-of" || args.front() == "-outputfile" || args.front() == "-outfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide output file name!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setOutFilename(args.front());
        }else if (args.front() == "-gf" || args.front() == "-gridfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide grid file name!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setGridFilename(args.front());
        }else if (args.front() =="-ascii"){
            wrap.setUseASCII(true);
        }else if (args.front() == "-af" || args.front() == "-anisotropyfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide anisotropy file name!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setAnisoFilename(args.front());
        }else if (args.front() == "-tf" || args.front() == "-transformfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide transform file name!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setTransformFilename(args.front());
        }else if (args.front() == "-conformalfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide conformal transform file name!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setConformalFilename(args.front());
        }else if (args.front() == "-lf" || args.front() == "-levellimitsfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide level limits file name!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setLevelLimitsFilename(args.front());
        }else if (args.front() == "-cf" || args.front() == "-customfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide custom file name!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setCustomFilename(args.front());
        }else if (args.front() == "-print" || args.front() == "-p"){
            wrap.setPrintPoints(true);
        }else if (args.front() == "-dim" || args.front() == "-dimensions"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide number of dimensions!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setNumDimensions(std::stoi(args.front()));
        }else if (args.front() == "-out" || args.front() == "-outputs"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide number of outputs!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setNumOutputs(std::stoi(args.front()));
        }else if (args.front() == "-dt" || args.front() == "-depth"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid depth!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setNumDepth(std::stoi(args.front()));
        }else if (args.front() == "-or" || args.front() == "-order"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid order!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setOrder(std::stoi(args.front()));
        }else if (args.front() == "-tt" || args.front() == "-type"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid depth type!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            TypeDepth depth_type = IO::getDepthTypeString(args.front());
            if (depth_type == type_none){
                cerr << "ERROR: " << args.front() << " is not a valid type!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setDepthType(depth_type);
        }else if (args.front() == "-1d" || args.front() == "-onedim"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -onedim!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            TypeOneDRule rule = IO::getRuleString(args.front());
            if (rule == rule_none){
                cerr << "ERROR: " << args.front() << " is not a valid rule!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setRule(rule);
        }else if (args.front() == "-ct" || args.front() == "-conformaltype"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -conformaltype!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            TypeConformalMap conformal_type = TasgridWrapper::getConfromalType(args.front().c_str());
            if (conformal_type == conformal_none){
                cerr << "ERROR: " << args.front() << " is not a valid conformal type!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setConformalType(conformal_type);
        }else if (args.front() == "-alpha"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -alpha!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setAlpha(std::stof(args.front()));
        }else if (args.front() == "-beta"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -beta!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setBeta(std::stof(args.front()));
        }else if (args.front() == "-tol" || args.front() == "-tolerance"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -tolerance!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setTolerance(std::stof(args.front()));
        }else if (args.front() == "-rout" || args.front() == "-refout"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -refout!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setRefOutput(std::stoi(args.front()));
        }else if (args.front() == "-ming" || args.front() == "-mingrowth"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -mingrowth!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setMinGrowth(std::stoi(args.front()));
        }else if (args.front() == "-rt" || args.front() == "-reftype"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid depth type!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            TypeRefinement ref  = IO::getTypeRefinementString(args.front());
            if (ref == refine_none){
                cerr << "ERROR: " << args.front() << " is not a valid refinement type!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setTypeRefinement(ref);
        }else if (args.front() == "-gpu" || args.front() == "-gpuid"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -gpuid  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setGPID(std::stoi(args.front()));
        }else if (args.front() == "-shift"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide valid -shift!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setShift(std::stof(args.front()));
        }else if (args.front() == "-wf" || args.front() == "-weightfile"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide weight file name!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setWeightFilename(args.front());
        }else if (args.front() == "-desc" || args.front() == "-description"){
            args.pop_front();
            if (args.empty()){
                cerr << "ERROR: must provide description!!!  For help see: ./tasgrid -help\n\n";
                return 1;
            }
            wrap.setDescription(args.front());
        }else if (args.front() == "-symm" || args.front() == "-symmetric"){
            wrap.setIsSymmetric(true);
        }else if (command == command_summary){
            wrap.setGridFilename(args.front());
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

        cout << " Usage: tasgrid <command> <option1> <value1> <option2> <value2> ... \n\n";

        cout << "Commands\t"       << "\tShorthand"   << "\tAction\n";
        cout << " -help\t"         << "\t\t-h,--help" << "\tshow verbose help options\n";
        cout << " -listtypes\t"    << "\t-lt\t"       << "\tlist accepted grid types and 1-D rules\n";
        cout << " -test\t"         << "\t\t\t"    << "\tperform a number of internal tests\n";
        cout << " -makeglobal\t"       << "\t-mg"     << "\t\tmake a grid from a global rule\n";
        cout << " -makesequence\t"     << "\t-ms"     << "\t\tmake a grid from a sequence rule\n";
        cout << " -makelocalpoly\t"    << "\t-mp"     << "\t\tmake a grid from a local polynomial rule\n";
        cout << " -makewavelet\t"      << "\t-mw"     << "\t\tmake a grid from a wavelet rule\n";
        cout << " -makefourier\t"      << "\t-mf"     << "\t\tmake a grid from a Fourier rule\n";
        cout << " -makequadrature"     << "\t-mq"     << "\t\tmake a quadrature\n";
        cout << " -makeexoquad\t"     << "\t-meq"     << "\t\tmake an exotic quadrature\n";
        cout << " -makeupdate\t"       << "\t-mu"     << "\t\tupdates an existing global/sequence/fourier grid\n";
        cout << " -setconformal\t"     << "\t-sc"     << "\t\tset conformal domain transform\n";
        cout << " -getquadrature\t"    << "\t-gq"     << "\t\toutput quadrature weights and points\n";
        cout << " -getinterweights"    << "\t-gi"     << "\t\toutput the interpolation weights\n";
        cout << " -getdiffweights"    << "\t-gd"     << "\t\toutput the differentiation weights\n";
        cout << " -getpoints\t"    << "\t-gp"     << "\t\toutput the points\n";
        cout << " -getneededpoints"    << "\t-gn"     << "\t\toutputs the points needing values to build an interpolant\n";
        cout << " -loadvalues\t"       << "\t-l"      << "\t\tload the values of the interpolated function\n";
        cout << " -evaluate\t"     << "\t-e"      << "\t\tevaluates the interpolant\n";
        cout << " -evalhierarchyd"     << "\t-ehd"      << "\t\tevaluates the hierarchical basis (dense output)\n";
        cout << " -evalhierarchys"     << "\t-ehs"      << "\t\tevaluates the hierarchical basis (sparse output)\n";
        cout << " -gethsupport\t"       << "\t-ghsup"     << "\t\tget the hierarchical support\n";
        cout << " -integrate\t"    << "\t-i"      << "\t\toutput the integral\n";
        cout << " -differentiate\t"    << "\t-d"      << "\t\tdifferentiates the interpolant\n";
        cout << " -getanisotropy\t"    << "\t-ga"     << "\t\testimates the anisotropic coefficients\n";
        cout << " -refineaniso\t"      << "\t-ra"     << "\t\trefines the grid\n";
        cout << " -refinesurp\t"       << "\t-rs"     << "\t\trefines the grid\n";
        cout << " -refine\t"       << "\t-r"      << "\t\trefines the grid\n";
        cout << " -cancelrefine\t"     << "\t-cr"     << "\t\tdiscards the last refinement (unless values are already loaded)\n";
        cout << " -mergerefine\t"     << "\t-mr"     << "\t\tcombines the loaded and needed points and discards the loaded values\n";
        cout << " -getcoefficients"    << "\t-gc"     << "\t\tget the hierarchical coefficients of the grid\n";
        cout << " -setcoefficients"    << "\t-sc"     << "\t\tset the hierarchical coefficients of the grid\n";
        cout << " -getpoly\t"      << "\t"      << "\t\tget polynomial space\n";
        cout << " -summary\t"      << "\t-s"      << "\t\twrites short description\n\n";

        cout << "Options\t\t"    << "\tShorthand"  << "\tValue"    << "\t\tAction\n";
        cout << " -help\t\t"     << "\thelp\t"     << "\t"         << "\t\tdisplay verbose information about this command\n";
        cout << " -dimensions\t"     << "\t-dim\t"     << "\t<int>"    << "\t\tset the number of dimensions\n";
        cout << " -outputs\t"    << "\t-out\t"     << "\t<int>"    << "\t\tset the number of outputs\n";
        cout << " -depth\t\t"    << "\t-dt\t"      << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)\n";
        cout << " -type\t\t"     << "\t-tt\t"      << "\t<type>"       << "\t\tset the type of the grid\n";
        cout << " -conformaltype\t"     << "\t-tt\t"      << "\t<type>"       << "\t\tset the type of the transformation\n";
        cout << " -onedim\t"     << "\t-1d\t"      << "\t<rule>"       << "\t\tset the one dimensional rule\n";
        cout << " -order\t\t"    << "\t-or\t"      << "\t<int>"    << "\t\tset the order for local polynomial and wavelet basis\n";
        cout << " -alpha\t\t"    << "\t\t"     << "\t<float>"      << "\t\tthe alpha parameter for Gegenbauer/Jacobi/Laguerre/Hermite quadrature\n";
        cout << " -beta\t\t"     << "\t\t"     << "\t<float>"      << "\t\tthe beta parameter for Jacobi quadrature\n";
        cout << " -tolerance\t"      << "\t-tol\t"     << "\t<float>"      << "\t\tset the tolerance for the refinement\n";
        cout << " -refout\t"     << "\t-rout\t"    << "\t<int>"    << "\t\tselect the output to use for the refinement\n";
        cout << " -mingrowth\t"      << "\t-ming\t"   << "\t<int>"    << "\t\tminimum number of new points\n";
        cout << " -reftype\t"    << "\t-rt\t"      << "\t<int>"    << "\t\tset the type of refinement\n";

        cout << " -gridfile\t"       << "\t-gf\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
        cout << " -xfile\t\t"    << "\t-xf\t"      << "\t<filename>"   << "\tset the name for the file with points\n";
        cout << " -valsfile\t"       << "\t-vf\t"      << "\t<filename>"   << "\tset the name for the file with values\n";
        cout << " -outputfile\t"     << "\t-of\t"      << "\t<filename>"   << "\tset the name for the output file\n";
        cout << " -anisotropyfile"   << "\t-af\t"      << "\t<filename>"   << "\tset the anisotropic weights\n";
        cout << " -transformfile\t"  << "\t-tf\t"      << "\t<filename>"   << "\tset the transformation of the domain\n";
        cout << " -conformalfile\t"  << "\t-tf\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain\n";
        cout << " -levellimitsfile"  << "\t-lf\t"      << "\t<filename>"   << "\tset the limits for the levels\n";
        cout << " -customfile\t"     << "\t-cf\t"      << "\t<filename>"   << "\tset the file with the custom-tabulated rule\n";
        cout << " -gpuid\t\t"    << "\t\t"     << "\t<int>\t"   << "\tset the gpu to use for evaluations\n";
        cout << " -print\t\t"    << "\t-p\t"       << "\t<none>"       << "\t\tprint to standard output\n";
        cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n";

        cout << endl;
    }else if (ht == help_command){
        if (com == command_makeglobal){
            cout << "Commands\t"     << "\tShorthand" << "\tAction\n";
            cout << " -makeglobal\t"     << "\t-mg"       << "\t\tmake a grid from a global rule\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions\n";
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs\n";
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)\n";
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of the grid\n";
            cout << " -onedim\t"     << "\tyes\t"     << "\t<rule>"       << "\t\tset the one dimensional rule\n";
            cout << " \t\t\t\t\t"    << "\t\tmust use a global rule (see manual)\n";
            cout << " -alpha\t\t"    << "\tsometimes" << "\t<float>"      << "\t\tthe alpha parameter for Gegenbauer/Jacobi/Laguerre/Hermite quadrature\n";
            cout << " \t\t\t\t\t"    << "\t\trequired for those rules\n";
            cout << " -beta\t\t"     << "\tsometimes" << "\t<float>"      << "\t\tthe beta parameter for Jacobi quadrature\n";
            cout << " \t\t\t\t\t"    << "\t\trequired for Jacobi rule\n";
            cout << " -gridfile\t"       << "\tno\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -anisotropyfile"   << "\tno\t"      << "\t<filename>"   << "\tset the anisotropic weights\n";
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain\n";
            cout << " -customfile\t"     << "\tsometimes" << "\t<filename>"   << "\tset the file with the custom-tabulated rule\n";
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map\n";
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain\n";
            cout << " \t\t\t\t\t"    << "\t\trequired for custom-tabulated rule\n";
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels\n";
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output\n";
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n\n";
            cout << "Note: -outputfile or -print output the points of the grid\n";
            cout << "Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_makesequence){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -makesequence\t"   << "\t-ms"     << "\t\tmake a grid from a sequence rule\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions\n\n";
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs\n";
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)\n";
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of the grid\n";
            cout << " -onedim\t"     << "\tyes\t"     << "\t<rule>"       << "\t\tset the one dimensional rule\n";
            cout << " \t\t\t\t\t"    << "\t\tmust use a sequence rule (see manual)\n";
            cout << " -gridfile\t"       << "\tno\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -anisotropyfile"   << "\tno\t"      << "\t<filename>"   << "\tset the anisotropic weights\n";
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain\n";
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map\n";
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain\n";
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels\n";
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output\n";
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n\n";
            cout << "Note: -outputfile or -print output the points of the grid\n";
            cout << "Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_makelocalp){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -makelocalpoly\t"  << "\t-mp"     << "\t\tmake a grid from a local polynomial rule\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions\n";
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs\n";
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)\n";
            cout << " -order\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the order for local polynomial basis\n";
            cout << " -onedim\t"     << "\tyes\t"     << "\t<rule>"       << "\t\tset the one dimensional rule\n";
            cout << " \t\t\t\t\t"    << "\t\tmust use a local polynomial rule (see manual)\n";
            cout << " -gridfile\t"       << "\tno\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain\n";
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map\n";
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain\n";
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels\n";
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output\n";
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n\n";
            cout << "Note: -outputfile or -print output the points of the grid\n";
            cout << "Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_makewavelet){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -makewavelet\t"    << "\t-mw"     << "\t\tmake a grid from a wavelet rule\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions\n";
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs\n";
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)\n";
            cout << " -order\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the order for the wavelet basis\n";
            cout << " -gridfile\t"       << "\tno\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain\n";
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map\n";
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain\n";
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels\n";
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output\n";
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n\n";
            cout << "Note: -outputfile or -print output the points of the grid\n";
            cout << "Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_makefourier){
            cout << "Commands\t"     << "\tShorthand" << "\tAction\n";
            cout << " -makefourier\t"     << "\t-mf"       << "\t\tmake a grid from a Fourier rule\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions\n";
            cout << " -outputs\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the number of outputs\n";
            cout << " -depth\t\t"    << "\tyes\t"     << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)\n";
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of the grid\n";
            cout << " -gridfile\t"       << "\tno\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -anisotropyfile"   << "\tno\t"      << "\t<filename>"   << "\tset the anisotropic weights\n";
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain\n";
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map\n";
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain\n";
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels\n";
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output\n";
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n\n";
            cout << "Note: -outputfile or -print output the points of the grid\n";
            cout << "Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_makequadrature){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -makequadrature"   << "\t-mq"     << "\t\tmake a quadrature\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -dimensions\t"     << "\tyes\t"     << "\t<int>"    << "\t\tset the number of dimensions\n";
            cout << " -depth\t\t"    << "\tyes\t"      << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)\n";
            cout << " -type\t\t"     << "\tsometimes"      << "\t<type>"       << "\t\tset the type of the grid\n";
            cout << " \t\t\t\t\t"    << "\t\trequired only for global rules (see manual)\n";
            cout << " -order\t\t"    << "\tsometimes"      << "\t<int>"    << "\t\tset the order for local polynomial basis\n";
            cout << " \t\t\t\t\t"    << "\t\trequired only for local polynomial and wavelet (see manual)\n";
            cout << " -onedim\t"     << "\tyes\t"      << "\t<rule>"       << "\t\tset the one dimensional rule\n";
            cout << " \t\t\t\t\t"    << "\t\tcan use any rule (see manual)\n";
            cout << " -alpha\t\t"    << "\tsometimes"  << "\t<float>"      << "\t\tthe alpha parameter for Gegenbauer/Jacobi/Laguerre/Hermite quadrature\n";
            cout << " \t\t\t\t\t"    << "\t\trequired for those rules\n";
            cout << " -beta\t\t"     << "\tsometimes"  << "\t<float>"      << "\t\tthe beta parameter for Jacobi quadrature\n";
            cout << " \t\t\t\t\t"    << "\t\trequired for Jacobi rule\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -anisotropyfile"   << "\tno\t"      << "\t<filename>"   << "\tset the anisotropic weights\n";
            cout << " -transformfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the transformation of the domain\n";
            cout << " -customfile\t"     << "\tsometimes"  << "\t<filename>"   << "\tset the file with the custom-tabulated rule\n";
            cout << " -conformaltype\t"  << "\tno\t"      << "\t<type>"       << "\t\tset the type of the map\n";
            cout << " -conformalfile\t"  << "\tno\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain\n";
            cout << " \t\t\t\t\t"    << "\t\trequired for custom-tabulated rule\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the points and weights of the grid\n";
            cout << "Note: at least one of -outputfile, or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_makeexoquad){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -makeexoquad"   << "\t\t-meq"     << "\t\tmake an exotic quadrature\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -depth\t\t"    << "\tyes\t"      << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)\n";
            cout << " -shift\t\t"    << "\tyes\t"      << "\t<float>"    << "\t\tset the shift of the weight function\n";
            cout << " -weightfile\t"    << "\tyes\t"      << "\t<filename>"    << "\tset the name of the file containing a\n";
            cout << " \t\t\t\t\t"    << "\t\tsurrogate/interpolant of the weight function;\n";
            cout << " \t\t\t\t\t"    << "\t\tmust be a TasmanianSparseGrid in ASCII format\n";
            cout << " -description\t"    << "\tyes\t"      << "\t<string>"    << "\tshort description of the quadrature\n";
            cout << " -symmetric\t"    << "\tno\t"      << "\t<none>"    << "\t\tdeclare that the weight function is symmetric"  << endl;
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the points and weights of the grid\n";
            cout << "Note: at least one of -outputfile, or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_update){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -makeupdate\t"     << "\t-mu"     << "\t\tupdates a new global or sequence grid\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -depth\t\t"    << "\tyes\t"      << "\t<int>"    << "\t\tset the depth of the grid (e.g. levels)\n";
            cout << " -type\t\t"     << "\tyes\t"      << "\t<type>"       << "\t\tset the type of the grid\n";
            cout << " -anisotropyfile"   << "\tno\t"      << "\t<filename>"   << "\tset the anisotropic weights\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the new points of the grid\n";
        }else if (com == command_setconformal){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -setconformal\t"     << "\t-sc"     << "\t\tset conformal transformation\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -conformaltype\t"  << "\tyes\t"      << "\t<type>"       << "\t\tset the type of the map\n";
            cout << " -conformalfile\t"  << "\tyes\t"      << "\t<filename>"   << "\tset the conformal transformation of the domain\n";
        }else if (com == command_getquadrature){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -getquadrature"    << "\t-gq"     << "\t\tmake a quadrature\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the points and weights of the grid\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_getcoefficients){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -getcoefficients"    << "\t-gc"     << "\t\tget the hierarchical coefficients\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the hierarchical coefficients of the sparse grid\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_setcoefficients){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -setcoefficients"    << "\t-sc"     << "\t\tset the hierarchical coefficients\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the hierarchical coefficients of the sparse grid\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_getinterweights){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -getinterweights"  << "\t-gi"     << "\t\toutput the interpolation weights\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the interpolation weight for each point in the xfile, see equation (1.2) in the manual\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_getdiffweights){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -getdiffweights"  << "\t-gd"     << "\t\toutput the differentiation weights\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"      << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the differentiation weight for each point in the xfile\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_getpoints){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -getpoints\t"      << "\t-gp"     << "\t\toutput the points\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the points of the grid\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_getneeded){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -getneededpoints"  << "\t-gn"     << "\t\toutputs the points needing values to build an interpolant\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the new points of the grid\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_loadvalues){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -loadvalues\t"     << "\t-l"      << "\t\tload the values of the interpolated function\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -valsfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the file with values\n\n";
        }else if (com == command_evaluate){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -evaluate\t"       << "\t-e"      << "\t\tevaluates the interpolant\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points\n";
            cout << " -gpuid\t\t"    << "\tno\t"     << "\t<int>\t"   << "\tset the gpu to use for evaluations\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output values of the interpolant at the points specified in the xfile\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_evalhierarchical_dense){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -evalhierarchyd"       << "\t-ehd"      << "\t\tevaluates the hierarchical basis (dense output)\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output values of the hierarchical basis functions in the xfile\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_evalhierarchical_sparse){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -evalhierarchys"       << "\t-ehs"      << "\t\tevaluates the hierarchical basis (sparse output)\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output values of the hierarchical basis functions in the xfile\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_gethsupport){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -gethsupport\t"       << "\t-ghsup"      << "\t\tget the hierarchical support\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_integrate){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -integrate\t"      << "\t-i"      << "\t\toutput the integral\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the integral if the loaded function, see equation (1.3) in the manual\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_differentiate){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -differentiate\t"       << "\t-d"      << "\t\tdifferentiates the interpolant\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -xfile\t\t"    << "\tyes\t"     << "\t<filename>"   << "\tset the name for the file with points\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print derivative (Jacobian matrix) of the interpolant at the points specified in the xfile\n";
            cout << "Note: at least one of -outputfile or -print must be specified, otherwise the command has no output\n\n";
        }else if (com == command_getanisocoeff){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -getanisotropy\t"    << "\t-ga"     << "\t\testimates the anisotropic coefficients\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of anisotropic coefficients\n";
            cout << " -refout\t"     << "\tsometimes" << "\t<int>"    << "\t\tselect the output to use\n";
            cout << " \t\t\t\t\t"    << "\t\trequired by global grids only\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the estimated anisotropic coefficients\n\n";
        }else if (com == command_refine_aniso){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -refineaniso\t"    << "\t-ra"     << "\t\trefines the grid\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"      << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tset the type of anisotropic refinement\n";
            cout << " -mingrowth\t"      << "\tno\t"      << "\t<int>"    << "\t\tminimum number of new points (defaults to 1)\n";
            cout << " -refout\t"     << "\tsometimes" << "\t<int>"    << "\t\tselect the output to use for the refinement\n";
            cout << " \t\t\t\t\t"    << "\t\trequired by global grids, for sequence grids defaults to -1 (use all outputs)\n";
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n";
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n\n";
            cout << "Note: -outputfile or -print output the new points of the grid\n\n";
        }else if (com == command_refine_surp){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -refinesurp\t"     << "\t-rs"     << "\t\trefines the grid\n\n";
            cout << "Accepted options:"  << endl;
            cout << "Options\t\t"    << "\tRequired"  << "\tValue"    << "\t\tAction\n";
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -tolerance\t"      << "\tyes\t"     << "\t<float>"      << "\t\tset the tolerance for the refinement\n";
            cout << " -reftype\t"    << "\tsometimes" << "\t<int>"    << "\t\tset the type of refinement\n";
            cout << " \t\t\t\t\t"    << "\t\trequired by local polynomial and wavelet grids\n";
            cout << " -refout\t"     << "\tsometimes" << "\t<int>"    << "\t\tselect the output to use for the refinement\n";
            cout << " \t\t\t\t\t"    << "\t\trequired by global grids, for sequence grids defaults to -1 (use all outputs)\n";
            cout << " -levellimitsfile"  << "\tno\t"      << "\t<filename>"   << "\tset the limits for the levels\n";
            cout << " -valsfile\t"       << "\t-vf\t"      << "\t<filename>"   << "\tset the correction weights for the surpluses\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n";
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n\n";
            cout << "Note: -outputfile or -print output the new points of the grid\n\n";
        }else if (com == command_refine){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -refine\t"     << "\t-r"    << "\t\trefines the grid\n\n";
            cout << "Accepted options:\n";
            cout << " -refine calls -refineaniso for Global and Sequence grids and -refinesurp otherwise\n";
            cout << " see \"-refineaniso help\" or \"-refinesurp help\" for the corresponding accepted options\n\n";
            cout << "Note: -outputfile or -print output the new points of the grid\n\n";
        }else if (com == command_refine_clear){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -cancelrefine\t"   << "\t-cr"     << "\t\tdiscards the last refinement (unless values are already loaded)\n\n";
            cout << "Accepted options:\n";
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n\n";
        }else if (com == command_refine_merge){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -mergerefine\t"   << "\t-mr"     << "\t\tmerges the loaded and needed points and discards any loaded values\n\n";
            cout << "Accepted options:\n";
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -ascii\t\t"    << "\t\t"       << "\t<none>"       << "\t\tuse ASCII grid file format\n\n";
        }else if (com == command_getpoly){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -getpoly\t"      << "\t\t"      << "\t\tget polynomial space\n\n";
            cout << "Accepted options:\n";
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file\n";
            cout << " -type\t\t"     << "\tyes\t"     << "\t<type>"       << "\t\tspecifies whether to use quadrature or interpolation\n";
            cout << " -outputfile\t"     << "\tno\t"      << "\t<filename>"   << "\tset the name for the output file\n";
            cout << " -print\t\t"    << "\tno\t"       << "\t<none>"       << "\t\tprint to standard output\n\n";
            cout << "Note: -outputfile or -print output the polynomial indexes\n\n";
        }else if (com == command_summary){
            cout << "Commands\t"     << "\tShorthand"   << "\tAction\n";
            cout << " -summary\t"    << "\t-s"      << "\t\twrites short description\n\n";
            cout << "Accepted options:\n";
            cout << " -gridfile\t"       << "\tyes\t"     << "\t<filename>"   << "\tset the name for the grid file\n\n";
            cout << "Note that 'tasgrid -s <filename>' is also accepted\n\n";
        }
    }else if (ht == help_listtypes){
        cout << "This only lists the strings associated with each option for spelling purposes.\n";
        cout << "Refer to the manual for details about each type.\n\n";
        cout << "List of global grids types:\n";
        cout << "  level     curved     hyperbolic     tensor\n";
        cout << "  iptotal   ipcurved   iphyperbolic   iptensor\n";
        cout << "  qptotal   qpcurved   qphyperbolic   qptensor\n\n";
        cout << "List of global grids rules:\n";
        cout << "         chebyshev          chebyshev-odd    clenshaw-curtis   clenshaw-curtis-zero\n";
        cout << "              leja               leja-odd              rleja              rleja-odd\n";
        cout << "     rleja-double2          rleja-double4      rleja-shifted     rleja-shifted-even\n";
        cout << "      max-lebesgue       max-lebesgue-odd       min-lebesgue       min-lebesgue-odd\n";
        cout << "         min-delta          min-delta-odd             fejer2\n";
        cout << "    gauss-legendre     gauss-legendre-odd    gauss-patterson       custom-tabulated\n";
        cout << "  gauss-gegenbauer   gauss-gegenbauer-odd       gauss-jacobi       gauss-jacobi-odd\n";
        cout << "    gauss-laguerre     gauss-laguerre-odd      gauss-hermite      gauss-hermite-odd\n";
        cout << "  gauss-chebyshev1   gauss-chebyshev1-odd   gauss-chebyshev2   gauss-chebyshev2-odd\n\n";
        cout << "List of sequence grids rules:\n";
        cout << "              leja              rleja     rleja-shifted\n";
        cout << "      max-lebesgue       min-lebesgue         min-delta \n\n";
        cout << "List of local polynomial grids rules:\n";
        cout << "      localp    localp-zero    semi-localp   localp-boundary\n\n";
        cout << "List of local wavelet grids rules:\n";
        cout << "      wavelet\n\n";
        cout << "List of local polynomial and wavelet refinement types:\n";
        cout << "      classic     parents    direction    fds\n\n";
        cout << "List of conformal maps:\n";
        cout << "      asin\n\n";
    }
}
