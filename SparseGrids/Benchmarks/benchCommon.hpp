#ifndef _TASMANIAN_BENCHMARK_COMMON_HPP
#define _TASMANIAN_BENCHMARK_COMMON_HPP

#include <iomanip>
#include <string>
#include <chrono>
#include <cmath>

#include "TasmanianSparseGrid.hpp"

using std::cout;
using std::endl;

// skips the first argument, then goes over the others
std::vector<std::string> stringArgs(int argc, char** argv){
    std::vector<std::string> args;
    for(int i = 1; i < argc; i++)
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

#endif
