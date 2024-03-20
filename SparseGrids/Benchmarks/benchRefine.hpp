#ifndef _TASMANIAN_BENCHMARK_REFINE_HPP
#define _TASMANIAN_BENCHMARK_REFINE_HPP

#include "benchCommon.hpp"

bool benchmark_refine(std::deque<std::string> &args){
    if (args.size() < 15) return false;

    // report the test parameters to reference later
    cout << "refine";
    for(auto &s : args) cout << " " << s;

    auto grid_family = getGridFamily(args);
    if (grid_family == GridFamily::none) return false;

    int num_dimensions, num_outputs, num_depth, order, iteratons;
    TypeDepth dtype;
    TypeOneDRule rule;
    TypeAcceleration acc;
    int device;

    TypeDepth ref_aniso_type;
    int aniso_min_growth;
    double surp_tolerance;
    TypeRefinement surp_criteria;
    int output;

    auto riter = readEntries(args.begin(), num_dimensions, num_outputs, num_depth, dtype, rule, order,
                             ref_aniso_type, aniso_min_growth, surp_tolerance, surp_criteria, output, iteratons, acc, device);

    std::string flavor = checkFlavor(riter, args.end());

    auto extra = extractWeightsLimits(grid_family, num_dimensions, dtype, riter, args.end());

    auto make_grid = getLambdaMakeGrid(grid_family, num_dimensions, num_outputs, num_depth, dtype, rule, order, extra);

    auto grid = make_grid();
    cout << " (points: " << grid.getNumPoints() << ")\n";
    loadGenericModel(grid);

    grid.enableAcceleration(acc, device);
    if (flavor != "auto")
        grid.favorSparseAcceleration((flavor == "sparse"));

    auto test_lambda = [&](int)->void{
            if (aniso_min_growth < 1) {
                grid.setSurplusRefinement(surp_tolerance, surp_criteria, output);
            } else {
                grid.setAnisotropicRefinement(ref_aniso_type, aniso_min_growth, output);
            }
        };

    // make one dry-run
    cout << setw(7) << testMethod<DryRun>(iteratons, test_lambda) << endl;

    return true;
}

#endif
