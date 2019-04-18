#ifndef _TASMANIAN_BENCHMARK_MAKEGRID_HPP
#define _TASMANIAN_BENCHMARK_MAKEGRID_HPP

#include "benchCommon.hpp"

bool benchmark_makegrid(std::deque<std::string> &args){
    if (args.size() != 7) return false;

    // report the test parameters to reference later
    cout << "makegrid";
    for(auto &s : args) cout << " " << s;
    cout << endl;

    TypeOneDRule grid_family = getGridFamily(args.front());
    if (grid_family == rule_none) return false;
    args.pop_front();

    auto riter = args.begin();
    int num_dimensions   = std::stoi(*riter++);
    int num_depth        = std::stoi(*riter++);
    TypeDepth dtype      = OneDimensionalMeta::getIOTypeString((*riter++).c_str());
    TypeOneDRule rule    = OneDimensionalMeta::getIORuleString((*riter++).c_str());
    int iteratons        = std::stoi(*riter++);
    int num_jumps        = std::stoi(*riter++);

    auto make_grid = [&]()->TasmanianSparseGrid{
        TasmanianSparseGrid grid;
        if (grid_family == rule_clenshawcurtis){
            grid.makeGlobalGrid(num_dimensions, 1, num_depth, dtype, rule);
        }else if (grid_family == rule_rleja){
            grid.makeSequenceGrid(num_dimensions, 1, num_depth, dtype, rule);
        }else if (grid_family == rule_localp){
            grid.makeLocalPolynomialGrid(num_dimensions, 1, num_depth, 1, rule);
        }else if (grid_family == rule_fourier){
            grid.makeFourierGrid(num_dimensions, 1, num_depth, dtype);
        }else if (grid_family == rule_wavelet){
            grid.makeWaveletGrid(num_dimensions, 1, num_depth, 1);
        }
        return grid;
    };

    num_jumps = std::max(num_jumps, 1); // make at least one test

    cout << setw(12) << "num points" << setw(12) << "miliseconds" << endl;
    for(int jump = 0; jump < num_jumps; jump++){
        auto grid = make_grid();
        //grid.printStats(); // uncomment to make sure the right grid is constructed

        auto time_start = std::chrono::system_clock::now();
        for(int i=0; i < iteratons; i++)
            auto grid_temp = make_grid();
        auto time_end = std::chrono::system_clock::now();

        long long elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
        int normalized = int( double(elapsed) / double(iteratons) );

        cout << setw(12) << grid.getNumPoints() << setw(12) << normalized << endl;
        num_depth += 1;
    }

    return true;
}

#endif
