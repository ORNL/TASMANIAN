#ifndef _TASMANIAN_BENCHMARK_EVALUATE_HPP
#define _TASMANIAN_BENCHMARK_EVALUATE_HPP

#include "benchCommon.hpp"

bool benchmark_evaluate(std::deque<std::string> &args){
    if (args.size() < 12) return false;

    // report the test parameters to reference later
    cout << "evaluate";
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
    int num_batch        = std::stoi(*riter++);
    int iteratons        = std::stoi(*riter++);
    int num_jumps        = std::stoi(*riter++);
    TypeAcceleration acc = AccelerationMeta::getIOAccelerationString((*riter++).c_str());
    int device           = std::stoi(*riter++);

    auto extra = extractWeightsLimits(grid_family, num_dimensions, dtype, riter, args.end());

    auto make_grid = [&]()->TasmanianSparseGrid{
        TasmanianSparseGrid grid;
        if (grid_family == rule_clenshawcurtis){
            grid.makeGlobalGrid(num_dimensions, num_outputs, num_depth, dtype, rule, extra.first, 0.0, 0.0, nullptr, extra.second);
        }else if (grid_family == rule_rleja){
            grid.makeSequenceGrid(num_dimensions, num_outputs, num_depth, dtype, rule, extra.first, extra.second);
        }else if (grid_family == rule_localp){
            grid.makeLocalPolynomialGrid(num_dimensions, num_outputs, num_depth, order, rule, extra.second);
        }else if (grid_family == rule_fourier){
            grid.makeFourierGrid(num_dimensions, num_outputs, num_depth, dtype, extra.first, extra.second);
        }else if (grid_family == rule_wavelet){
            grid.makeWaveletGrid(num_dimensions, num_outputs, num_depth, order, extra.second);
        }
        return grid;
    };

    num_jumps = std::max(num_jumps, 1); // make at least one test
    std::vector<std::vector<double>> inputs((size_t) iteratons);
    int seed = 44;
    for(auto &inp : inputs) inp = getRandomVector(num_dimensions, num_batch, seed++);

    for(int jump = 0; jump < num_jumps; jump++){
        auto grid = make_grid();
        if (jump == 0) cout << " (points: " << grid.getNumPoints() << ")" << endl;
        loadGenericModel(grid);

        grid.enableAcceleration(acc);
        if (AccelerationMeta::isAccTypeGPU(acc)) grid.setGPUID(device);

        // make one dry-run
        std::vector<double> result;
        grid.evaluateBatch(inputs.back(), result);
        //grid.printStats(); // uncomment to make sure the right grid is constructed

        auto time_start = std::chrono::system_clock::now();
        for(size_t i=0; i < (size_t) iteratons; i++)
            grid.evaluateBatch(inputs[i], result);
        auto time_end = std::chrono::system_clock::now();

        long long elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
        int normalized = int( double(elapsed) / double(iteratons) );

        cout << setw(7) << normalized << endl;
        num_outputs *= 2;
    }

    return true;
}

#endif
