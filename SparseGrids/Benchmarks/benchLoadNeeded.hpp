#ifndef _TASMANIAN_BENCHMARK_LOADNEEDED_HPP
#define _TASMANIAN_BENCHMARK_LOADNEEDED_HPP

#include "benchCommon.hpp"

bool benchmark_loadneeded(std::deque<std::string> &args){
    if (args.size() < 11) return false;

    // report the test parameters to reference later
    cout << "load";
    for(auto &s : args) cout << " " << s;

    TypeOneDRule grid_family = getGridFamily(args.front());
    if (grid_family == rule_none) return false;
    args.pop_front();

    auto riter = args.begin();
    int num_dimensions   = std::stoi(*riter++);
    int num_outputs      = std::stoi(*riter++);
    int num_depth        = std::stoi(*riter++);
    TypeDepth dtype      = OneDimensionalMeta::getIOTypeString((*riter++).c_str());
    TypeOneDRule rule    = OneDimensionalMeta::getIORuleString((*riter++).c_str());
    int order            = std::stoi(*riter++);
    int iteratons        = std::stoi(*riter++);
    int num_jumps        = std::stoi(*riter++);
    TypeAcceleration acc = AccelerationMeta::getIOAccelerationString((*riter++).c_str());
    int device           = std::stoi(*riter++);

    auto extra = extractWeightsLimits(grid_family, num_dimensions, dtype, riter, args.end());

    auto make_grid = getLambdaMakeGrid(grid_family, num_dimensions, num_outputs, num_depth, dtype, rule, order, extra);

    for(int jump = 0; jump < num_jumps; jump++){
        auto grid = make_grid();
        if (jump == 0) cout << " (points: " << grid.getNumPoints() << ")" << endl;
        grid.enableAcceleration(acc);
        if (AccelerationMeta::isAccTypeGPU(acc)) grid.setGPUID(device);

        auto values = getGenericModel((size_t) grid.getNumDimensions(),
                                      (size_t) grid.getNumOutputs(),
                                      grid.getNeededPoints());
        grid.loadNeededPoints(values); // dry-run

        auto time_start = std::chrono::system_clock::now();
        for(size_t i=0; i < (size_t) iteratons; i++)
            grid.loadNeededPoints(values);
        auto time_end = std::chrono::system_clock::now();

        long long elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
        int normalized = int( double(elapsed) / double(iteratons) );

        cout << setw(7) << normalized << endl;
        num_outputs *= 2;
    }

    return true;
}

#endif
