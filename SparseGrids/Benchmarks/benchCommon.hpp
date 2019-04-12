#ifndef _TASMANIAN_BENCHMARK_COMMON_HPP
#define _TASMANIAN_BENCHMARK_COMMON_HPP

#include <chrono>

#include "tasgridCLICommon.hpp"

enum TestFuction{
    test_none,
    test_evaluate
};

// skips the first argument, then goes over the others
std::vector<std::string> stringArgs(int argc, char** argv){
    std::vector<std::string> args;
    for(int i = argc-1; i > 0; i--)
        args.push_back(std::string(argv[i]));
    return args;
}

bool hasHelp(std::string const &arg){
    std::string lower(arg.size(), ' ');
    std::transform(arg.begin(), arg.end(), lower.begin(),
        [](char c)->char{
            return static_cast<char>(std::tolower(static_cast<int>(c)));
        });

    auto pos = lower.find("help");
    if ((pos < lower.size()) && (pos  + 4 <= lower.size()))
        return (lower.substr(pos, pos + 4).compare("help") == 0);

    return false;
}

TestFuction getTest(std::string const &s){
    std::map<std::string, TestFuction> str_to_test = {
        {"evaluate", test_evaluate},
    };

    TestFuction test = test_none;
    try{
        test = str_to_test.at(s);
    }catch(std::out_of_range &){
        cout << "ERROR: Unknown test: " << s << endl;

    }
    return test;
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
