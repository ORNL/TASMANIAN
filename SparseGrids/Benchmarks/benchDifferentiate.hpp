#ifndef _TASMANIAN_BENCHMARK_DIFFERENTIATE_HPP
#define _TASMANIAN_BENCHMARK_DIFFERENTIATE_HPP

#include "benchCommon.hpp"

bool benchmark_differentiate(std::deque<std::string> &args){
    if (args.size() < 9) return false;

    // report the test parameters to reference later
    cout << "differentiate";
    for(auto &s : args) cout << " " << s;

    auto grid_family = getGridFamily(args);
    if (grid_family == GridFamily::none) return false;

    int num_dimensions, num_outputs, num_depth, order, iteratons, num_jumps;
    TypeDepth dtype;
    TypeOneDRule rule;

    auto riter = readEntries(args.begin(), num_dimensions, num_outputs, num_depth, dtype, rule, order, iteratons, num_jumps);

    auto extra = extractWeightsLimits(grid_family, num_dimensions, dtype, riter, args.end());

    auto make_grid = getLambdaMakeGrid(grid_family, num_dimensions, num_outputs, num_depth, dtype, rule, order, extra);

    num_jumps = std::max(num_jumps, 1); // make at least one test

    std::vector<std::vector<double>> inputs = getRandomVectors<double>(iteratons, num_dimensions, 1);

    for(int jump = 0; jump < num_jumps; jump++){
        auto grid = make_grid();
        if (jump == 0) cout << " (points: " << grid.getNumPoints() << ")" << endl;
        loadGenericModel(grid);
        if (jump == 0) cout << "    note: reporting total time (ms), does not normalize (divide) by iterations\n";

        // make one dry-run, using RawTime since each call to differentiate is very fast
        std::vector<double> result;
        cout << setw(7) << testMethod<DryRun, RawTime>(iteratons, [&](int i)->void{ grid.differentiate(inputs[i], result); }) << endl;;
        num_outputs *= 2;
    }

    return true;
}

#endif
