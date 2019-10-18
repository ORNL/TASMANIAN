#ifndef _TASMANIAN_BENCHMARK_MAKEGRID_HPP
#define _TASMANIAN_BENCHMARK_MAKEGRID_HPP

#include "benchCommon.hpp"

bool benchmark_makegrid(std::deque<std::string> &args){
    if (args.size() < 7) return false;

    // report the test parameters to reference later
    cout << "makegrid";
    for(auto &s : args) cout << " " << s;
    cout << endl;

    auto grid_family = getGridFamily(args);
    if (grid_family == GridFamily::none) return false;

    int num_dimensions, num_depth, iteratons, num_jumps;
    TypeDepth dtype;
    TypeOneDRule rule;

    auto riter = readEntries(args.begin(), num_dimensions, num_depth, dtype, rule, iteratons, num_jumps);

    auto extra = extractWeightsLimits(grid_family, num_dimensions, dtype, riter, args.end());

    int num_outputs = 1, order = 1;
    auto make_grid = getLambdaMakeGrid(grid_family, num_dimensions, num_outputs, num_depth, dtype, rule, order, extra);

    num_jumps = std::max(num_jumps, 1); // make at least one test

    cout << setw(12) << "num points" << setw(12) << "miliseconds" << endl;
    for(int jump = 0; jump < num_jumps; jump++){
        auto grid = make_grid();
        //grid.printStats(); // uncomment to make sure the right grid is constructed

        cout << setw(12) << grid.getNumPoints() << setw(12)
             << testMethod(iteratons, [&](int)->void{ auto grid_temp = make_grid(); })
             << endl;
        num_depth += 1;
    }

    return true;
}

#endif
