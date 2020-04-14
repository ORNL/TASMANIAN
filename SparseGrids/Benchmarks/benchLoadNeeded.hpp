#ifndef _TASMANIAN_BENCHMARK_LOADNEEDED_HPP
#define _TASMANIAN_BENCHMARK_LOADNEEDED_HPP

#include "benchCommon.hpp"

bool benchmark_loadneeded(std::deque<std::string> &args){
    if (args.size() < 11) return false;

    // report the test parameters to reference later
    cout << "loadneeded";
    for(auto &s : args) cout << " " << s;

    auto grid_family = getGridFamily(args);
    if (grid_family == GridFamily::none) return false;

    int num_dimensions, num_outputs, num_depth, order, iteratons, num_jumps, device;
    TypeDepth dtype;
    TypeOneDRule rule;
    TypeAcceleration acc;

    auto riter = readEntries(args.begin(), num_dimensions, num_outputs, num_depth, dtype, rule, order, iteratons, num_jumps, acc, device);

    std::string flavor = checkFlavor(riter, args.end());

    auto extra = extractWeightsLimits(grid_family, num_dimensions, dtype, riter, args.end());

    auto make_grid = getLambdaMakeGrid(grid_family, num_dimensions, num_outputs, num_depth, dtype, rule, order, extra);

    for(int jump = 0; jump < num_jumps; jump++){
        auto grid = make_grid();
        //grid.printStats(); // uncomment to make sure the right grid is constructed

        if (jump == 0) cout << " (points: " << grid.getNumPoints() << ")" << endl;
        grid.enableAcceleration(acc, device);
        if (flavor != "auto")
            grid.favorSparseAcceleration((flavor == "sparse"));

        auto values = getGenericModel((size_t) grid.getNumDimensions(),
                                      (size_t) grid.getNumOutputs(),
                                      grid.getNeededPoints());

        cout << setw(7) << testMethod<DryRun>(iteratons, [&](int)->void{ grid.loadNeededPoints(values); }) << endl;
        num_outputs *= 2;
    }

    return true;
}

#endif
