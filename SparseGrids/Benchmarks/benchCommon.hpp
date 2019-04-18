#ifndef _TASMANIAN_BENCHMARK_COMMON_HPP
#define _TASMANIAN_BENCHMARK_COMMON_HPP

#include <chrono>

#include "tasgridCLICommon.hpp"

enum BenchFuction{
    bench_none,
    bench_make,
    bench_evaluate
};

BenchFuction getTest(std::string const &s){
    std::map<std::string, BenchFuction> str_to_test = {
        {"evaluate", bench_evaluate},
        {"makegrid", bench_make},
    };

    try{
        return str_to_test.at(s);
    }catch(std::out_of_range &){
        cout << "ERROR: Unknown test: " << s << endl;
        return bench_none;
    }
}

TypeOneDRule getGridFamily(std::string const &s){
    std::map<std::string, TypeOneDRule> str_to_rule = {
        {"global",   rule_clenshawcurtis},
        {"sequence", rule_rleja},
        {"localp",   rule_localp},
        {"fourier",  rule_fourier},
        {"wavelet",  rule_wavelet},
    };

    TypeOneDRule rule = rule_none;
    try{
        rule = str_to_rule.at(s);
    }catch(std::out_of_range &){
        cout << "ERROR: Unknown grid type: " << s << endl;

    }
    return rule;
}

template<typename IteratorToList>
std::pair<std::vector<int>, std::vector<int>>
extractWeightsLimits(TypeOneDRule grid_family, int num_dimensions, TypeDepth dtype,
                     IteratorToList &arg, IteratorToList const &argend){
    std::vector<int> anisotropic_weights;
    if (grid_family != rule_localp && grid_family != rule_wavelet){
        int num_weights = (OneDimensionalMeta::getControurType(dtype) == type_curved) ? 2 * num_dimensions : num_dimensions;
        for(int i=0; i<num_weights && arg != argend; i++)
            anisotropic_weights.push_back(std::stoi(*arg++));
    }
    std::vector<int> level_limits;
    for(int i=0; i<num_dimensions && arg != argend; i++)
        level_limits.push_back(std::stoi(*arg++));
    return std::make_pair(anisotropic_weights, level_limits);
}

std::vector<double> getRandomVector(int dim1, int dim2, long int seed){
    std::vector<double> x(Utils::size_mult(dim1, dim2));
    std::minstd_rand park_miller(seed);
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    for(auto &v : x) v = unif(park_miller);
    return x;
}

std::vector<double> getGenericModel(size_t num_dimensions, size_t num_outputs,
                                    std::vector<double> const &points){
    // model output k = k * exp(sum(x_1 ... x_dims))
    size_t num_points = points.size() / num_dimensions;
    std::vector<double> values(num_points * num_outputs);

    auto ip = points.begin();
    auto iv = values.begin();
    while(ip != points.end()){
        double exponent = std::exp( std::accumulate(ip, ip + num_dimensions, 0.0) );
        std::advance(ip, num_dimensions);
        for(size_t i = 0; i < num_outputs; i++)
            *iv++ = double(i) * exponent;
    }

    return values;
}

void loadGenericModel(TasmanianSparseGrid &grid){
    if (grid.getNumNeeded() == 0) return;

    auto values = getGenericModel((size_t) grid.getNumDimensions(),
                                  (size_t) grid.getNumOutputs(),
                                  grid.getNeededPoints());

    grid.loadNeededPoints(values);
}

#endif
