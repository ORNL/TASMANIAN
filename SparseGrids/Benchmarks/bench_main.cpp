#include "benchMakeGrid.hpp"
#include "benchEvaluate.hpp"

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
    if (test == bench_make)
        pass = benchmark_makegrid(args);
    if (test == bench_evaluate)
        pass = benchmark_evaluate(args);

    if (!pass) // if problem with inputs
        printHelp(test);

    return (pass) ? 0 : 1;
}

void printHelp(BenchFuction test){
    if (test == bench_none){
        cout << "\nusage: ./benchmark <function> <parameters>\n\n";
        cout << "functions: evaluate\n";
        cout << "\n see: ./benchmark <function> help\n";
    }else if (test == bench_make){
        cout << "\nusage: ./benchmark makegrid <grid> <dims> <depth> <type> <rule> <iters> <jumps>\n\n";
        cout << "grid  : global, sequence, localp, wavelet, fourier\n";
        cout << "dims  : number of dimensions\n";
        cout << "depth : grid density\n";
        cout << "type  : level, iptotal, etc.; ignored if not used by the grid\n";
        cout << "rule  : rleja, clenshaw-curtis, etc.; ignored for wavelet and fourier grids\n";
        cout << "iters : number of times to repeat the function call\n";
        cout << "jumps : how many times to increment <depth> by 1\n";
    }else if (test == bench_evaluate){
        cout << "\nusage: ./benchmark evaluate <grid> <dims> <outs> <depth> <type> <rule> <order> <batch> <iters> <jumps> <acc> <gpu>\n\n";
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
    }
    cout << endl;
}
