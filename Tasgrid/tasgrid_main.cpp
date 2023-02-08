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
        cout << "Tasmanian Sparse Grids  version: " << TasmanianSparseGrid::getVersion() << "\n";
        if ((std::string(TasmanianSparseGrid::getGitCommitHash()).compare("Tasmanian git hash is not available here") != 0)
            && (std::string(TasmanianSparseGrid::getGitCommitHash()).find("Release") != 0)){
            cout << "                git commit hash: " << TasmanianSparseGrid::getGitCommitHash() << "\n";
            cout << "                cmake cxx flags: " << TasmanianSparseGrid::getCmakeCxxFlags() << "\n";
        }
        cout << "                        license: " << TasmanianSparseGrid::getLicense() << "\n";
        if (TasmanianSparseGrid::isOpenMPEnabled()){
            cout << "          OpenMP multithreading: Enabled\n";
        }else{
            cout << "          OpenMP multithreading: Disabled\n";
        }
        std::string gpu_backend = "none";
        if (TasmanianSparseGrid::isCudaEnabled()) gpu_backend = "CUDA";
        if (TasmanianSparseGrid::isHipEnabled()) gpu_backend = "ROCm/HIP";
        if (TasmanianSparseGrid::isDpcppEnabled()) gpu_backend = "oneAPI/DPC++";
        cout << "          GPU backend framework: " << gpu_backend << "\n";
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
        cout << "\n";
        if (anyGPU){
            int numGPUs = TasmanianSparseGrid::getNumGPUs();
            if (numGPUs > 0){
                cout << "                 Available GPUs:" << "\n";
                for(int i=0; i<numGPUs; i++){
                    std::string name = TasmanianSparseGrid::getGPUName(i);
                    int memory = TasmanianSparseGrid::getGPUMemory(i);
                    cout << setw(11) << i << ":" << setw(20) << name << " with" << setw(7) << memory << "MB of RAM\n";
                }
            }else{
                cout << "        Available GPUs: none\n";
            }
        }

        cout << "\n";
        return 0;
    }

    // help with interface commands
    if (args.front() == "-listtypes" || args.front() == "-lt"){
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
        cout << "ERROR: unknown command " << args.front() << "\n";
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
        }else if (command == command_summary || command == command_using_construct){
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

    if (ht == help_generic){
        cout << R"help(
    Usage: tasgrid <command> <option1> <value1> <option2> <value2> ...

Commands                Shorthand    Action
 -help                  -h,--help    show verbose help options
 -listtypes             -lt          list accepted grid types and 1-D rules
 -test                               perform a number of internal tests
 -makeglobal            -mg          make a grid from a global rule
 -makesequence          -ms          make a grid from a sequence rule
 -makelocalpoly         -mp          make a grid from a local polynomial rule
 -makewavelet           -mw          make a grid from a wavelet rule
 -makefourier           -mf          make a grid from a Fourier rule
 -makequadrature        -mq          make a quadrature
 -makeexoquad           -meq         make an exotic quadrature
 -makeupdate            -mu          updates an existing global/sequence/fourier grid
 -setconformal          -sc          set conformal domain transform
 -getquadrature         -gq          output quadrature weights and points
 -getinterweights       -gi          output the interpolation weights
 -getdiffweights        -gd          output the differentiation weights
 -getpoints             -gp          output the points
 -getneededpoints       -gn          outputs the points needing values to build an interpolant
 -loadvalues            -l           load the values of the interpolated function
 -evaluate              -e           evaluates the interpolant
 -evalhierarchyd        -ehd         evaluates the hierarchical basis (dense output)
 -evalhierarchys        -ehs         evaluates the hierarchical basis (sparse output)
 -gethsupport           -ghsup       get the hierarchical support
 -integrate             -i           output the integral
 -differentiate         -d           differentiates the interpolant
 -getanisotropy         -ga          estimates the anisotropic coefficients
 -refineaniso           -ra          refines the grid
 -refinesurp            -rs          refines the grid
 -refine                -r           refines the grid
 -cancelrefine          -cr          discards the last refinement (unless values are already loaded)
 -mergerefine           -mr          combines the loaded and needed points and discards the loaded values
 -using-construct                    prints simple string indicating whether dynamic construction is in use
 -getconstructpnts      -gcp         get points for dynamic construction
 -loadconstructed       -lcp         load points for dynamic construction
 -getcoefficients       -gc          get the hierarchical coefficients of the grid
 -setcoefficients       -sc          set the hierarchical coefficients of the grid
 -getpoly                            get polynomial space
 -summary               -s           writes short description

Options                 Shorthand  Value         Action
 -help                  help                     display verbose information about this command
 -dimensions            -dim       <int>         set the number of dimensions
 -outputs               -out       <int>         set the number of outputs
 -depth                 -dt        <int>         set the depth of the grid (e.g. levels)
 -type                  -tt        <type>        set the type of the grid
 -conformaltype         -tt        <type>        set the type of the transformation
 -onedim                -1d        <rule>        set the one dimensional rule
 -order                 -or        <int>         set the order for local polynomial and wavelet basis
 -alpha                            <float>       the alpha parameter for Gegenbauer/Jacobi/Laguerre/Hermite quadrature
 -beta                             <float>       the beta parameter for Jacobi quadrature
 -tolerance             -tol       <float>       set the tolerance for the refinement
 -refout                -rout      <int>         select the output to use for the refinement
 -mingrowth             -ming      <int>         minimum number of new points
 -reftype               -rt        <int>         set the type of refinement

 -gridfile              -gf        <filename>    set the name for the grid file
 -xfile                 -xf        <filename>    set the name for the file with points
 -valsfile              -vf        <filename>    set the name for the file with values
 -outputfile            -of        <filename>    set the name for the output file
 -anisotropyfile        -af        <filename>    set the anisotropic weights
 -transformfile         -tf        <filename>    set the transformation of the domain
 -conformalfile         -tf        <filename>    set the conformal transformation of the domain
 -levellimitsfile       -lf        <filename>    set the limits for the levels
 -customfile            -cf        <filename>    set the file with the custom-tabulated rule
 -gpuid                            <int>         set the gpu to use for evaluations
 -print                 -p         <none>        print to standard output
 -ascii                            <none>        use ASCII grid file format

)help";
    }else if (ht == help_command){
        switch(com){
            case command_makeglobal:
                cout << R"help(
Commands             Shorthand    Action
 -makeglobal         -mg          make a grid from a global rule

Accepted options:
Options              Required     Value         Action
 -dimensions         yes          <int>         set the number of dimensions
 -outputs            yes          <int>         set the number of outputs
 -depth              yes          <int>         set the depth of the grid (e.g. levels)
 -type               yes          <type>        set the type of the grid
 -onedim             yes          <rule>        set the one dimensional rule
                                                must use a global rule (see manual)
 -alpha              sometimes    <float>       the alpha parameter for Gegenbauer/Jacobi/Laguerre/Hermite quadrature
                                                required for those rules
 -beta               sometimes    <float>       the beta parameter for Jacobi quadrature
                                                required for Jacobi rule
 -gridfile           no           <filename>    set the name for the grid file
 -outputfile         no           <filename>    set the name for the output file
 -anisotropyfile     no           <filename>    set the anisotropic weights
 -transformfile      no           <filename>    set the transformation of the domain
 -customfile         sometimes    <filename>    set the file with the custom-tabulated rule
 -conformaltype      no           <type>        set the type of the map
 -conformalfile      no           <filename>    set the conformal transformation of the domain
                                                required for custom-tabulated rule
 -levellimitsfile    no           <filename>    set the limits for the levels
 -print              no           <none>        print to standard output
 -ascii                           <none>        use ASCII grid file format
Note: -outputfile or -print output the points of the grid
Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output
)help";
                break;
            case command_makesequence:
                cout << R"help(
Commands             Shorthand    Action
 -makesequence       -ms          make a grid from a sequence rule

Accepted options:
Options              Required     Value         Action
 -dimensions         yes          <int>         set the number of dimensions
 -outputs            yes          <int>         set the number of outputs
 -depth              yes          <int>         set the depth of the grid (e.g. levels)
 -type               yes          <type>        set the type of the grid
 -onedim             yes          <rule>        set the one dimensional rule
                                                must use a sequence rule (see manual)
 -gridfile           no           <filename>    set the name for the grid file
 -outputfile         no           <filename>    set the name for the output file
 -anisotropyfile     no           <filename>    set the anisotropic weights
 -transformfile      no           <filename>    set the transformation of the domain
 -conformaltype      no           <type>        set the type of the map
 -conformalfile      no           <filename>    set the conformal transformation of the domain
 -levellimitsfile    no           <filename>    set the limits for the levels
 -print              no           <none>        print to standard output
 -ascii                           <none>        use ASCII grid file format
Note: -outputfile or -print output the points of the grid
Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output
)help";
                break;
            case command_makelocalp:
                cout << R"help(
Commands            Shorthand    Action
 -makelocalpoly     -mp          make a grid from a local polynomial rule

Accepted options:
Options             Required     Value         Action
 -dimensions        yes          <int>         set the number of dimensions
 -outputs           yes          <int>         set the number of outputs
 -depth             yes          <int>         set the depth of the grid (e.g. levels)
 -order             yes          <int>         set the order for local polynomial basis
 -onedim            yes          <rule>        set the one dimensional rule
                                               must use a local polynomial rule (see manual)
 -gridfile          no           <filename>    set the name for the grid file
 -outputfile        no           <filename>    set the name for the output file
 -transformfile     no           <filename>    set the transformation of the domain
 -conformaltype     no           <type>        set the type of the map
 -conformalfile     no           <filename>    set the conformal transformation of the domain
 -levellimitsfile   no           <filename>    set the limits for the levels
 -print             no           <none>        print to standard output
 -ascii                          <none>        use ASCII grid file format
Note: -outputfile or -print output the points of the grid
Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_makewavelet:
                cout << R"help(
Commands             Shorthand    Action
 -makewavelet        -mw          make a grid from a wavelet rule

Accepted options:
Options              Required     Value         Action
 -dimensions         yes          <int>         set the number of dimensions
 -outputs            yes          <int>         set the number of outputs
 -depth              yes          <int>         set the depth of the grid (e.g. levels)
 -order              yes          <int>         set the order for the wavelet basis
 -gridfile           no           <filename>    set the name for the grid file
 -outputfile         no           <filename>    set the name for the output file
 -transformfile      no           <filename>    set the transformation of the domain
 -conformaltype      no           <type>        set the type of the map
 -conformalfile      no           <filename>    set the conformal transformation of the domain
 -levellimitsfile    no           <filename>    set the limits for the levels
 -print              no           <none>        print to standard output
 -ascii                           <none>        use ASCII grid file format
Note: -outputfile or -print output the points of the grid
Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_makefourier:
                cout << R"help(
Commands             Shorthand    Action
 -makefourier        -mf          make a grid from a Fourier rule

Accepted options:
Options              Required     Value         Action
 -dimensions         yes          <int>         set the number of dimensions
 -outputs            yes          <int>         set the number of outputs
 -depth              yes          <int>         set the depth of the grid (e.g. levels)
 -type               yes          <type>        set the type of the grid
 -gridfile           no           <filename>    set the name for the grid file
 -outputfile         no           <filename>    set the name for the output file
 -anisotropyfile     no           <filename>    set the anisotropic weights
 -transformfile      no           <filename>    set the transformation of the domain
 -conformaltype      no           <type>        set the type of the map
 -conformalfile      no           <filename>    set the conformal transformation of the domain
 -levellimitsfile    no           <filename>    set the limits for the levels
 -print              no           <none>        print to standard output
 -ascii                           <none>        use ASCII grid file format
Note: -outputfile or -print output the points of the grid
Note: at least one of -gridfile, -outputfile, or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_makequadrature:
                cout << R"help(
Commands            Shorthand    Action
 -makequadrature    -mq          make a quadrature

Accepted options:
Options             Required     Value         Action
 -dimensions        yes          <int>         set the number of dimensions
 -depth             yes          <int>         set the depth of the grid (e.g. levels)
 -type              sometimes    <type>        set the type of the grid
                                               required only for global rules (see manual)
 -order             sometimes    <int>         set the order for local polynomial basis
                                               required only for local polynomial and wavelet (see manual)
 -onedim            yes          <rule>        set the one dimensional rule
                                               can use any rule (see manual)
 -alpha             sometimes    <float>       the alpha parameter for Gegenbauer/Jacobi/Laguerre/Hermite quadrature
                                               required for those rules
 -beta              sometimes    <float>       the beta parameter for Jacobi quadrature
                                               required for Jacobi rule
 -outputfile        no           <filename>    set the name for the output file
 -anisotropyfile    no           <filename>    set the anisotropic weights
 -transformfile     no           <filename>    set the transformation of the domain
 -customfile        sometimes    <filename>    set the file with the custom-tabulated rule
 -conformaltype     no           <type>        set the type of the map
 -conformalfile     no           <filename>    set the conformal transformation of the domain
                                               required for custom-tabulated rule
 -print             no           <none>        print to standard output
Note: -outputfile or -print output the points and weights of the grid
Note: at least one of -outputfile, or -print must be specified, otherwise the command has no output

)help";
                break;

            case command_makeexoquad:
                cout << R"help(
Commands         Shorthand    Action
 -makeexoquad    -meq         make an exotic quadrature

Accepted options:
Options          Required     Value         Action
 -depth          yes          <int>         set the depth of the grid (e.g. levels)
 -shift          yes          <float>       set the shift of the weight function
 -weightfile     yes          <filename>    set the name of the file containing a
                                            surrogate/interpolant of the weight function
                                            must be a TasmanianSparseGrid in ASCII format
 -description    yes          <string>      short description of the quadrature
 -symmetric      no           <none>        declare that the weight function is symmetric   endl
 -outputfile     no           <filename>    set the name for the output file
 -print          no           <none>        print to standard output
Note: -outputfile or -print output the points and weights of the grid
Note: at least one of -outputfile, or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_update:
                cout << R"help(
Commands            Shorthand    Action
 -makeupdate        -mu          updates a new global or sequence grid

Accepted options:
Options             Required     Value         Action
 -depth             yes          <int>         set the depth of the grid (e.g. levels)
 -type              yes          <type>        set the type of the grid
 -anisotropyfile    no           <filename>    set the anisotropic weights
 -outputfile        no           <filename>    set the name for the output file
 -print             no           <none>        print to standard output
Note: -outputfile or -print output the new points of the grid

)help";
                break;
            case command_setconformal:
                cout << R"help(
Commands           Shorthand    Action
 -setconformal     -sc          set conformal transformation

Accepted options:
Options            Required     Value         Action
 -conformaltype    yes          <type>        set the type of the map
 -conformalfile    yes          <filename>    set the conformal transformation of the domain

)help";
                break;
            case command_getquadrature:
                cout << R"help(
Commands           Shorthand    Action
 -getquadrature    -gq          make a quadrature

Accepted options:
Options            Required     Value         Action
 -gridfile         yes          <filename>    set the name for the grid file
 -outputfile       no           <filename>    set the name for the output file
 -print            no           <none>        print to standard output
Note: -outputfile or -print output the points and weights of the grid
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_getcoefficients:
                cout << R"help(
Commands             Shorthand    Action
 -getcoefficients    -gc          get the hierarchical coefficients

Accepted options:
Options              Required     Value         Action
 -gridfile           yes          <filename>    set the name for the grid file
 -outputfile         no           <filename>    set the name for the output file
 -print              no           <none>        print to standard output
Note: -outputfile or -print output the hierarchical coefficients of the sparse grid
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_setcoefficients:
                cout << R"help(
Commands             Shorthand    Action
 -setcoefficients    -sc          set the hierarchical coefficients

Accepted options:
Options              Required     Value         Action
 -gridfile           yes          <filename>    set the name for the grid file
 -outputfile         no           <filename>    set the name for the output file
 -print              no           <none>        print to standard output
Note: -outputfile or -print output the hierarchical coefficients of the sparse grid
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_getinterweights:
                cout << R"help(
Commands             Shorthand    Action
 -getinterweights    -gi          output the interpolation weights

Accepted options:
Options              Required     Value         Action
 -gridfile           yes          <filename>    set the name for the grid file
 -xfile              yes          <filename>    set the name for the file with points
 -outputfile         no           <filename>    set the name for the output file
 -print              no           <none>        print to standard output
Note: -outputfile or -print output the interpolation weight for each point in the xfile, see equation (1.2) in the manual
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_getdiffweights:
                cout << R"help(
Commands            Shorthand    Action
 -getdiffweights    -gd          output the differentiation weights

Accepted options:
Options             Required     Value         Action
 -gridfile          yes          <filename>    set the name for the grid file
 -xfile             yes          <filename>    set the name for the file with points
 -outputfile        no           <filename>    set the name for the output file
 -print             no           <none>        print to standard output
Note: -outputfile or -print output the differentiation weight for each point in the xfile
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_getpoints:
                cout << R"help(
Commands        Shorthand    Action
 -getpoints     -gp          output the points

Accepted options:
Options         Required     Value         Action
 -gridfile      yes          <filename>    set the name for the grid file
 -outputfile    no           <filename>    set the name for the output file
 -print         no           <none>        print to standard output
Note: -outputfile or -print output the points of the grid
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_getneeded:
                cout << R"help(
Commands             Shorthand    Action
 -getneededpoints    -gn          outputs the points needing values to build an interpolant

Accepted options:
Options              Required     Value         Action
 -gridfile           yes          <filename>    set the name for the grid file
 -outputfile         no           <filename>    set the name for the output file
 -print              no           <none>        print to standard output
Note: -outputfile or -print output the new points of the grid
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_loadvalues:
                cout << R"help(
Commands        Shorthand    Action
 -loadvalues    -l           provides values of the model outputs at the needed grid points

Accepted options:
Options         Required     Value         Action
 -gridfile      yes          <filename>    set the name for the grid file
 -valsfile      -vf          <filename>    set the name for the file with values
Note: the -valsfile must contains rows equal to the number of needed points or
      (if there are no needed points) the number of loaded points
      the number of columns must match the number of outputs set for the grid

)help";
                break;
            case command_evaluate:
                cout << R"help(
Commands        Shorthand    Action
 -evaluate      -e           evaluates the interpolant

Accepted options:
Options         Required     Value         Action
 -gridfile      yes          <filename>    set the name for the grid file
 -xfile         yes          <filename>    set the name for the file with points
 -gpuid         no           <int>         set the gpu to use for evaluations
 -outputfile    no           <filename>    set the name for the output file
 -print         no           <none>        print to standard output
Note: -outputfile or -print output values of the interpolant at the points specified in the xfile
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_evalhierarchical_dense:
                cout << R"help(
Commands            Shorthand    Action
 -evalhierarchyd    -ehd         evaluates the hierarchical basis (dense output)

Accepted options:
Options             Required     Value         Action
 -gridfile          yes          <filename>    set the name for the grid file
 -xfile             yes          <filename>    set the name for the file with points
 -outputfile        no           <filename>    set the name for the output file
 -print             no           <none>        print to standard output
Note: -outputfile or -print output values of the hierarchical basis functions in the xfile
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_evalhierarchical_sparse:
                cout << R"help(
Commands            Shorthand    Action
 -evalhierarchys    -ehs         evaluates the hierarchical basis (sparse output)

Accepted options:
Options             Required     Value         Action
 -gridfile          yes          <filename>    set the name for the grid file
 -xfile             yes          <filename>    set the name for the file with points
 -outputfile        no           <filename>    set the name for the output file
 -print             no           <none>        print to standard output
Note: -outputfile or -print output values of the hierarchical basis functions in the xfile
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_gethsupport:
                cout << R"help(
Commands         Shorthand    Action
 -gethsupport    -ghsup       get the hierarchical support

Accepted options:
Options          Required     Value         Action
 -gridfile       yes          <filename>    set the name for the grid file
 -outputfile     no           <filename>    set the name for the output file
 -print          no           <none>        print to standard output
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_integrate:
                cout << R"help(
Commands        Shorthand    Action
 -integrate     -i           output the integral

Accepted options:
Options         Required     Value         Action
 -gridfile      yes          <filename>    set the name for the grid file
 -outputfile    no           <filename>    set the name for the output file
 -print         no           <none>        print to standard output
Note: -outputfile or -print output the integral if the loaded function, see equation (1.3) in the manual
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_differentiate:
                cout << R"help(
Commands           Shorthand    Action
 -differentiate    -d           differentiates the interpolant

Accepted options:
Options            Required     Value         Action
 -gridfile         yes          <filename>    set the name for the grid file
 -xfile            yes          <filename>    set the name for the file with points
 -outputfile       no           <filename>    set the name for the output file
 -print            no           <none>        print to standard output
Note: -outputfile or -print derivative (Jacobian matrix) of the interpolant at the points specified in the xfile
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_getanisocoeff:
                cout << R"help(
Commands           Shorthand    Action
 -getanisotropy    -ga          estimates the anisotropic coefficients

Accepted options:
Options            Required     Value         Action
 -gridfile         yes          <filename>    set the name for the grid file
 -type             yes          <type>        set the type of anisotropic coefficients
 -refout           sometimes    <int>         select the output to use
                                              required by global grids only
 -outputfile       no           <filename>    set the name for the output file
 -print            no           <none>        print to standard output
Note: -outputfile or -print output the estimated anisotropic coefficients

)help";
                break;
            case command_refine_aniso:
                cout << R"help(
Commands             Shorthand    Action
 -refineaniso        -ra          refines the grid

Accepted options:
Options              Required     Value         Action
 -gridfile           yes          <filename>    set the name for the grid file
 -type               yes          <type>        set the type of anisotropic refinement
 -mingrowth          no           <int>         minimum number of new points (defaults to 1)
 -refout             sometimes    <int>         select the output to use for the refinement
                                                required by global grids, for sequence grids
                                                defaults to -1 (use all outputs)
 -levellimitsfile    no           <filename>    set the limits for the levels
 -outputfile         no           <filename>    set the name for the output file
 -print              no           <none>        print to standard output
 -ascii              <none>       use ASCII grid file format
Note: -outputfile or -print output the new points of the grid

)help";
                break;
            case command_refine_surp:
                cout << R"help(
Commands             Shorthand    Action
 -refinesurp         -rs          refines the grid

Accepted options:
Options              Required     Value     Action
 -gridfile           yes          <filename>    set the name for the grid file
 -tolerance          yes          <float>       set the tolerance for the refinement
 -reftype            sometimes    <int>         set the type of refinement
                                                required by local polynomial and wavelet grids
 -refout             sometimes    <int>         select the output to use for the refinement
                                                required by global grids, for sequence grids
                                                defaults to -1 (use all outputs)
 -levellimitsfile    no           <filename>    set the limits for the levels
 -valsfile           no           <filename>    set the correction weights for the surpluses
 -outputfile         no           <filename>    set the name for the output file
 -print              no           <none>        print to standard output
 -ascii                           <none>        use ASCII grid file format
Note: -outputfile or -print output the new points of the grid

)help";
                break;
            case command_refine:
                cout << R"help(
Commands      Shorthand    Action
 -refine      -r     refines the grid

 -refine calls -refineaniso for Global and Sequence grids and -refinesurp otherwise
 see "-refineaniso help" or "-refinesurp help" for the corresponding accepted options
Note: -outputfile or -print output the new points of the grid

)help";
                break;
            case command_refine_clear:
                cout << R"help(
Commands          Shorthand    Action
 -cancelrefine    -cr          discards the last refinement (unless values are already loaded)
Accepted options:
Options           Required     Value         Action
 -gridfile        yes          <filename>    set the name for the grid file
 -ascii                        <none>        use ASCII grid file format

)help";
                break;
            case command_refine_merge:
                cout << R"help(
Commands         Shorthand    Action
 -mergerefine    -mr          merges the loaded and needed points and discards any loaded values

Accepted options:
Options          Required     Value         Action
    -gridfile    yes          <filename>    set the name for the grid file
    -ascii       <none>       use ASCII grid file format

)help";
                break;
            case command_using_construct:
                cout << R"help(
Commands             Action
 -using-construct    writes short string indicating whether dynamic construction is on

Accepted options:
Options       Required     Value         Action
 -gridfile    yes          <filename>    set the name for the grid file
Note that 'tasgrid -using-construct <filename>' is also accepted

)help";
                break;
            case command_get_candidate_construction:
                cout << R"help(
Commands              Shorthand    Action
 -getconstructpnts    -gcp         get points for dynamic construction

Accepted options:
Options               Required     Value         Action
 -gridfile            yes          <filename>    set the name for the grid file
 -refout              sometimes    <int>         select the output to use for the refinement
 -anisotropyfile      sometimes    <filename>    set the anisotropic weights
 -levellimitsfile     no           <filename>    set the limits for the levels
 -outputfile          no           <filename>    set the name for the output file
 -print               no           <none>        print to standard output
 -ascii                            <none>        use ASCII grid file format

)help";
                break;
            case command_load_construction:
                cout << R"help(
Commands             Shorthand    Action
 -loadconstructed    -lcp         load points for dynamic construction

Accepted options:
Options              Required     Value         Action
 -gridfile           yes          <filename>    set the name for the grid file
 -xfile              yes          <filename>    set the name for the file with points
 -valsfile           yes          <filename>    set the name for the file with values

)help";
                break;
            case command_getpoly:
                cout << R"help(
Commands        Shorthand    Action
 -getpoly                    get polynomial space

Accepted options:
Options         Required     Value         Action
 -gridfile      yes          <filename>    set the name for the grid file
 -type          yes          <type>        specifies whether to use quadrature or interpolation
 -outputfile    no           <filename>    set the name for the output file
 -print         no           <none>        print to standard output
Note: -outputfile or -print output the polynomial indexes

)help";
                break;
            case command_summary:
                cout << R"help(
Commands      Shorthand    Action
 -summary     -s           writes short description

Accepted options:
Options       Required     Value         Action
 -gridfile    yes          <filename>    set the name for the grid file
Note that 'tasgrid -s <filename>' is also accepted

)help";
                break;
            case command_getpointsindex:
            case command_getneededindex:
                cout << R"help(
Commands              Action
 -getpointsindexes    returns the integer multi-index set of the points
 -getneededindexes    returns the integer multi-index set of the needed points
Note: these commands exist primarily for debugging purposes

Accepted options:
Options         Required     Value         Action
 -gridfile      yes          <filename>    set the name for the grid file
 -outputfile    no           <filename>    set the name for the output file
 -print         no           <none>        print to standard output
Note: -outputfile or -print derivative (Jacobian matrix) of the interpolant at the points specified in the xfile
Note: at least one of -outputfile or -print must be specified, otherwise the command has no output

)help";
                break;
            case command_none:
                throw std::runtime_error("ERROR: incorrect command for help");
        }
    }else if (ht == help_listtypes){
        cout << R"help(
This only lists the strings associated with each option for spelling purposes.
Refer to the manual for details about each type.

List of global grids types:
  level     curved     hyperbolic     tensor
  iptotal   ipcurved   iphyperbolic   iptensor
  qptotal   qpcurved   qphyperbolic   qptensor

List of global grids rules:
         chebyshev          chebyshev-odd    clenshaw-curtis   clenshaw-curtis-zero
              leja               leja-odd              rleja              rleja-odd
     rleja-double2          rleja-double4      rleja-shifted     rleja-shifted-even
      max-lebesgue       max-lebesgue-odd       min-lebesgue       min-lebesgue-odd
         min-delta          min-delta-odd             fejer2
    gauss-legendre     gauss-legendre-odd    gauss-patterson       custom-tabulated
  gauss-gegenbauer   gauss-gegenbauer-odd       gauss-jacobi       gauss-jacobi-odd
    gauss-laguerre     gauss-laguerre-odd      gauss-hermite      gauss-hermite-odd
  gauss-chebyshev1   gauss-chebyshev1-odd   gauss-chebyshev2   gauss-chebyshev2-odd

List of sequence grids rules:
              leja              rleja     rleja-shifted
      max-lebesgue       min-lebesgue         min-delta

List of local polynomial grids rules:
      localp    localp-zero    semi-localp   localp-boundary

List of local wavelet grids rules:
      wavelet

List of local polynomial and wavelet refinement types:
      classic     parents    direction    fds

List of conformal maps:
      asin

)help";
    }
}
