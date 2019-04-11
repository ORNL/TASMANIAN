#include "benchEvaluate.hpp"

void printHelp(TestFuction test);

int main(int argc, char** argv){

    //cout << " Phruuuuphrrr " << endl; // this is the sound that the Tasmanian devil makes

    std::vector<std::string> args = stringArgs(argc, argv);

    if (args.empty() || hasHelp(args.back())){
        printHelp(test_none);
        return 0;
    }

    auto test = getTest(args.back());
    args.pop_back();
    if ((test == test_none) || args.empty() || hasHelp(args.back())){
        printHelp(test);
        return (test == test_none) ? 1 : 0;
    }

    if (test == test_evaluate){
        bench_evaluate(args);
    }

    return 0;
}

void printHelp(TestFuction test){
    if (test == test_none){
        cout << "\nusage: ./benchmark <function> <parameters>\n\n";
        cout << "functions: evaluate\n";
        cout << "\n see: ./benchmark <function> help\n";
    }else if (test == test_evaluate){
        cout << "\nusage: ./benchmark evaluate <grid> <dims> <outs> <type> <depth> <rule> <order> <batch> <iters> <jumps> <gpu>\n\n";
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
        cout << "gpu   : cuda device ID, set to -1 to use accel_cpu_blas or -2 for accel_none\n";
    }
    cout << endl;
}
