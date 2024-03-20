#include "benchMakeGrid.hpp"
#include "benchLoadNeeded.hpp"
#include "benchEvaluate.hpp"
#include "benchDifferentiate.hpp"
#include "benchInterpolationWeights.hpp"
#include "benchRefine.hpp"

void printHelp(BenchFuction test);

int main(int argc, const char** argv){

    //cout << " Phruuuuphrrr " << endl; // this is the sound that the Tasmanian devil makes

    std::deque<std::string> args = stringArgs(argc, argv);

    if (args.empty() || hasHelp(args.front())){
        printHelp(bench_none);
        return 0;
    }

    auto test = getTest(args.front());
    args.pop_front();
    if ((test == bench_none) || args.empty() || hasHelp(args.front())){
        printHelp(test);
        return (test == bench_none) ? 1 : 0;
    }

    bool pass = true; // check if the rest of the inputs are OK
    switch(test){
        case bench_make:
            pass = benchmark_makegrid(args);
            break;
        case bench_loadneeded:
            pass = benchmark_loadneeded(args);
            break;
        case bench_evaluate:
        case bench_evaluate_mixed:
            pass = benchmark_evaluate(args, (test == bench_evaluate_mixed));
            break;
        case bench_differentiate:
            pass = benchmark_differentiate(args);
            break;
        case bench_iweights:
            pass = benchmark_iweights(args);
            break;
        case bench_refine:
            pass = benchmark_refine(args);
            break;
        default:
            throw std::runtime_error("bench_main.cpp: invalid test type in switch statement!");
    }

    if (!pass) // if problem with inputs
        printHelp(test);

    return (pass) ? 0 : 1;
}

void printHelp(BenchFuction test){
    if (test == bench_none){
        cout << "\nusage: ./benchmark <function> <parameters>\n\n";
        cout << "functions: makegrid, loadneeded, evaluate(-mixed), differentiate, iweights, refine\n";
        cout << "\n see: ./benchmark <function> help\n";
    }else if (test == bench_make){
        cout << "\nusage: ./benchmark makegrid <grid> <dims> <depth> <type> <rule> <iters> <jumps> <aniso>\n\n";
        cout << "grid  : global, sequence, localp, wavelet, fourier\n";
        cout << "dims  : number of dimensions\n";
        cout << "depth : grid density\n";
        cout << "type  : level, iptotal, etc.; ignored if not used by the grid\n";
        cout << "rule  : rleja, clenshaw-curtis, etc.; ignored for wavelet and fourier grids\n";
        cout << "iters : number of times to repeat the function call\n";
        cout << "jumps : how many times to increment <depth> by 1\n";
        cout << "aniso : (optional) list of anisotropic weights and level limits\n";
        cout << "      : anisotropic weights come first (if used by the grid), then level limits\n";
    }else if (test == bench_loadneeded){
        cout << "\nusage: ./benchmark loadneeded <grid> <dims> <outs> <depth> <type> <rule> <order> <iters> <jumps> <acc> <gpu> <extra>\n\n";
        cout << "grid  : global, sequence, localp, wavelet, fourier\n";
        cout << "dims  : number of dimensions\n";
        cout << "outs  : number of outputs\n";
        cout << "depth : grid density\n";
        cout << "type  : level, iptotal, etc.; ignored if not used by the grid\n";
        cout << "rule  : rleja, clenshaw-curtis, etc.; ignored for wavelet and fourier grids\n";
        cout << "order : -1, 0, 1, 2; ignored if not used by the grid\n";
        cout << "iters : number of times to repeat the function call\n";
        cout << "jumps : how many times to double <outs>\n";
        cout << "acc   : acceleration type, e.g., gpu-cuda, cpu-blas, none, etc.\n";
        cout << "gpu   : cuda device ID; ignored for cpu acceleration\n";
        cout << "extra : (optional) sparse/dense flavor and/or list of anisotropic weights and level limits\n";
        cout << "      : anisotropic weights come first (if used by the grid), then level limits\n";
    }else if (test == bench_evaluate || test == bench_evaluate_mixed){
        cout << "\nusage: ./benchmark evaluate <grid> <dims> <outs> <depth> <type> <rule> <order> <batch> <iters> <jumps> <acc> <gpu> <extra>\n\n";
        cout << "grid  : global, sequence, localp, wavelet, fourier\n";
        cout << "dims  : number of dimensions\n";
        cout << "outs  : number of outputs\n";
        cout << "depth : grid density\n";
        cout << "type  : level, iptotal, etc.; ignored if not used by the grid\n";
        cout << "rule  : rleja, clenshaw-curtis, etc.; ignored for wavelet and fourier grids\n";
        cout << "order : -1, 0, 1, 2; ignored if not used by the grid\n";
        cout << "batch : number of points to use for the evaluate command\n";
        cout << "iters : number of times to repeat the function call\n";
        cout << "jumps : how many times to double <outs>\n";
        cout << "acc   : acceleration type, e.g., gpu-cuda, cpu-blas, none, etc.\n";
        cout << "gpu   : cuda device ID; ignored for cpu acceleration\n";
        cout << "extra : (optional) sparse/dense flavor and/or list of anisotropic weights and level limits\n";
        cout << "      : anisotropic weights come first (if used by the grid), then level limits\n";
    }else if (test == bench_differentiate){
        cout << "\nusage: ./benchmark differentiate <grid> <dims> <outs> <depth> <type> <rule> <order> <iters> <jumps>\n\n";
        cout << "grid  : global, sequence, localp, wavelet, fourier\n";
        cout << "dims  : number of dimensions\n";
        cout << "outs  : number of outputs\n";
        cout << "depth : grid density\n";
        cout << "type  : level, iptotal, etc.; ignored if not used by the grid\n";
        cout << "rule  : rleja, clenshaw-curtis, etc.; ignored for wavelet and fourier grids\n";
        cout << "order : -1, 0, 1, 2; ignored if not used by the grid\n";
        cout << "iters : number of times to repeat the function call\n";
        cout << "jumps : how many times to double <outs>\n";
        cout << "extra : (optional) sparse/dense flavor and/or list of anisotropic weights and level limits\n";
        cout << "      : anisotropic weights come first (if used by the grid), then level limits\n";
    }else if (test == bench_iweights){
        cout << "\nusage: ./benchmark iweights <grid> <dims> <depth> <type> <rule> <order> <iters> <jumps> <aniso>\n\n";
        cout << "grid  : global, sequence, localp, wavelet, fourier\n";
        cout << "dims  : number of dimensions\n";
        cout << "depth : grid density\n";
        cout << "type  : level, iptotal, etc.; ignored if not used by the grid\n";
        cout << "rule  : rleja, clenshaw-curtis, etc.; ignored for wavelet and Fourier grids\n";
        cout << "order : -1, 0, 1, 2; ignored if not used by the grid\n";
        cout << "iters : number of times to repeat the function call\n";
        cout << "jumps : how many times to increase <depth> by 1\n";
        cout << "aniso : (optional) list of anisotropic weights and level limits\n";
        cout << "      : anisotropic weights come first (if used by the grid), then level limits\n";
    }else if (test == bench_refine){
        cout << "\nusage: ./benchmark refine <grid> <dims> <outs> <depth> <type> <rule> <order> <ref-type-depth> <min-growth> <surp-tolerance> <surp-criteria> <output> <iters> <acc> <gpu> <extra>\n\n";
        cout << "grid  : global, sequence, localp, wavelet, fourier\n";
        cout << "dims  : number of dimensions\n";
        cout << "outs  : number of outputs\n";
        cout << "depth : grid density\n";
        cout << "type  : level, iptotal, etc.; ignored if not used by the grid\n";
        cout << "rule  : rleja, clenshaw-curtis, etc.; ignored for wavelet and fourier grids\n";
        cout << "order : -1, 0, 1, 2; ignored if not used by the grid\n\n";

        cout << "ref-type-depth : (anisotropic refinement) refinement type, e.g., iptotal, ipcurved\n";
        cout << "min-growth     : (anisotropic refinement) minumum number of refinement points, use 0 to switch to surplus refinement\n";
        cout << "surp-tolerance : (surplus refinement) tolerance\n";
        cout << "surp-criteria  : (surplus refinement) selection criteria, e.g., stable, fds\n";
        cout << "output         : (all refinement) output to use in the refinement\n\n";

        cout << "iters : number of times to repeat the function call\n";
        cout << "acc   : acceleration type, e.g., gpu-cuda, cpu-blas, none, etc.\n";
        cout << "gpu   : cuda device ID; ignored for cpu acceleration\n";
        cout << "extra : (optional) list of anisotropic weights and level limits\n";
        cout << "      : anisotropic weights come first (if used by the grid), then level limits\n";
    }
    cout << endl;
}
