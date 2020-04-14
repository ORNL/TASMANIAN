#ifndef _TASMANIAN_BENCHMARK_EVALUATE_HPP
#define _TASMANIAN_BENCHMARK_EVALUATE_HPP

#include "benchCommon.hpp"

bool benchmark_evaluate(std::deque<std::string> &args, bool use_mixed){
    if (args.size() < 12) return false;

    // report the test parameters to reference later
    cout << "evaluate";
    for(auto &s : args) cout << " " << s;

    auto grid_family = getGridFamily(args);
    if (grid_family == GridFamily::none) return false;

    int num_dimensions, num_outputs, num_depth, order, num_batch, iteratons, num_jumps, device;
    TypeDepth dtype;
    TypeOneDRule rule;
    TypeAcceleration acc;

    auto riter = readEntries(args.begin(), num_dimensions, num_outputs, num_depth, dtype, rule, order, num_batch, iteratons, num_jumps, acc, device);

    std::string flavor = checkFlavor(riter, args.end());

    auto extra = extractWeightsLimits(grid_family, num_dimensions, dtype, riter, args.end());

    auto make_grid = getLambdaMakeGrid(grid_family, num_dimensions, num_outputs, num_depth, dtype, rule, order, extra);

    num_jumps = std::max(num_jumps, 1); // make at least one test

    std::vector<std::vector<double>> inputs = getRandomVectors<double>(iteratons, num_dimensions, num_batch);
    std::vector<std::vector<float>> inputsf = getRandomVectors<float>((use_mixed) ? iteratons : 0, num_dimensions, num_batch);

    for(int jump = 0; jump < num_jumps; jump++){
        auto grid = make_grid();
        if (jump == 0) cout << " (points: " << grid.getNumPoints() << ")" << endl;
        loadGenericModel(grid);

        grid.enableAcceleration(acc, device);
        if (flavor != "auto")
            grid.favorSparseAcceleration((flavor == "sparse"));

        // make one dry-run
        std::vector<double> result;
        cout << setw(7) << testMethod<DryRun>(iteratons, [&](int i)->void{ grid.evaluateBatch(inputs[i], result); });

        if (use_mixed){
            std::vector<float> resultf;
            cout << setw(7) << testMethod<DryRun>(iteratons, [&](int i)->void{ grid.evaluateBatch(inputsf[i], resultf); });
        }
        cout << endl;
        num_outputs *= 2;
    }

    return true;
}

#endif
