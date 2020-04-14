#ifndef _TASMANIAN_BENCHMARK_COMMON_HPP
#define _TASMANIAN_BENCHMARK_COMMON_HPP

#include <chrono>

#include "tasgridCLICommon.hpp"

enum BenchFuction{
    bench_none,
    bench_make,
    bench_loadneeded,
    bench_evaluate,
    bench_evaluate_mixed,
    bench_iweights
};

BenchFuction getTest(std::string const &s){
    std::map<std::string, BenchFuction> str_to_test = {
        {"evaluate", bench_evaluate},
        {"evaluate-mixed", bench_evaluate_mixed},
        {"loadneeded", bench_loadneeded},
        {"makegrid", bench_make},
        {"iweights", bench_iweights}
    };

    try{
        return str_to_test.at(s);
    }catch(std::out_of_range &){
        cout << "ERROR: Unknown test: " << s << endl;
        return bench_none;
    }
}

enum class GridFamily{
    none, global, sequence, localp, fourier, wavelet
};

GridFamily getGridFamily(std::string const &s){
    std::map<std::string, GridFamily> str_to_rule = {
        {"global",   GridFamily::global},
        {"sequence", GridFamily::sequence},
        {"localp",   GridFamily::localp},
        {"fourier",  GridFamily::fourier},
        {"wavelet",  GridFamily::wavelet},
    };

    GridFamily grid_family = GridFamily::none;
    try{
        grid_family = str_to_rule.at(s);
    }catch(std::out_of_range &){
        cout << "ERROR: Unknown grid type: " << s << endl;

    }
    return grid_family;
}

//! \brief Convert \b s.front() to a grid type indicated by the canonical rules, then \b s.pop_front() the processed string.
GridFamily getGridFamily(std::deque<std::string> &s){
    GridFamily grid_family = getGridFamily(s.front());
    s.pop_front();
    return grid_family;
}

//! \brief Convert a string to int and advance the iterator.
template<typename StringListIterator>
void readEntry(StringListIterator &iter, int &val){
    val = std::stoi(*iter++);
}
//! \brief Convert a string to TypeDepth and advance the iterator.
template<typename StringListIterator>
void readEntry(StringListIterator &iter, TypeDepth &val){
    val = IO::getDepthTypeString(*iter++);
}
//! \brief Convert a string to TypeOneDRule and advance the iterator.
template<typename StringListIterator>
void readEntry(StringListIterator &iter, TypeOneDRule &val){
    val = IO::getRuleString(*iter++);
}
//! \brief Convert a string to TypeAcceleration and advance the iterator.
template<typename StringListIterator>
void readEntry(StringListIterator &iter, TypeAcceleration &val){
    val = AccelerationMeta::getIOAccelerationString((*iter++).c_str());
}

//! \brief Template to terminate recursion, read one entry of type \b ValType and return the iterator.
template<typename StringListIterator, typename ValType>
StringListIterator readEntries(StringListIterator iter, ValType &val){
    readEntry(iter, val);
    return iter;
}
//! \brief Template to read multiple entries from a string list, returns an iterator past the last processed string.
template<typename StringListIterator, typename ValType1, typename ValType2, typename...Other>
StringListIterator readEntries(StringListIterator iter, ValType1 &val1, ValType2 &val2, Other & ...others){
    readEntry(iter, val1);
    return readEntries(iter, val2, others...);
}

//! \brief If the current iterator is pointing to a string dense/sparse return the string and advance, else return "auto" and do nothing.
template<typename StringListIterator>
std::string checkFlavor(StringListIterator &iter, StringListIterator argend){
    if (iter == argend) return "auto";
    if (*iter == "sparse"){
        iter++;
        return "sparse";
    }else if (*iter == "dense"){
        iter++;
        return "dense";
    }else{
        return "auto";
    }
}

template<typename IteratorToList>
std::pair<std::vector<int>, std::vector<int>>
extractWeightsLimits(GridFamily grid_family, int num_dimensions, TypeDepth dtype,
                     IteratorToList &arg, IteratorToList const &argend){
    std::vector<int> anisotropic_weights;
    if (grid_family != GridFamily::localp && grid_family != GridFamily::wavelet){
        int num_weights = (OneDimensionalMeta::getControurType(dtype) == type_curved) ? 2 * num_dimensions : num_dimensions;
        for(int i=0; i<num_weights && arg != argend; i++)
            anisotropic_weights.push_back(std::stoi(*arg++));
    }
    std::vector<int> level_limits;
    for(int i=0; i<num_dimensions && arg != argend; i++)
        level_limits.push_back(std::stoi(*arg++));
    return std::make_pair(anisotropic_weights, level_limits);
}

inline std::function<TasmanianSparseGrid()>
getLambdaMakeGrid(GridFamily grid_family, int const &num_dimensions, const int &num_outputs,
                  int const &num_depth, TypeDepth const &dtype, TypeOneDRule const &rule, int const &order,
                  std::pair<std::vector<int>, std::vector<int>> const &extra){
    if (grid_family == GridFamily::global){
        return [&]()->TasmanianSparseGrid{
            return makeGlobalGrid(num_dimensions, num_outputs, num_depth, dtype, rule, extra.first, 0.0, 0.0, nullptr, extra.second);
        };
    }else if (grid_family == GridFamily::sequence){
        return [&]()->TasmanianSparseGrid{
            return makeSequenceGrid(num_dimensions, num_outputs, num_depth, dtype, rule, extra.first, extra.second);
        };
    }else if (grid_family == GridFamily::localp){
        return [&]()->TasmanianSparseGrid{
            return makeLocalPolynomialGrid(num_dimensions, num_outputs, num_depth, order, rule, extra.second);
        };
    }else if (grid_family == GridFamily::fourier){
        return [&]()->TasmanianSparseGrid{
            return makeFourierGrid(num_dimensions, num_outputs, num_depth, dtype, extra.first, extra.second);
        };
    }else{ // default - wavelet
        return [&]()->TasmanianSparseGrid{
            return makeWaveletGrid(num_dimensions, num_outputs, num_depth, order, extra.second);
        };
    }
}

template<typename T = double>
std::vector<T> getRandomVector(int dim1, int dim2, long int seed){
    std::vector<T> x(Utils::size_mult(dim1, dim2));
    std::minstd_rand park_miller(seed);
    std::uniform_real_distribution<T> unif(-1.0, 1.0);
    for(auto &v : x) v = unif(park_miller);
    return x;
}

template<typename T>
std::vector<std::vector<T>> getRandomVectors(int num_vectors, int dim1, int dim2){
    std::vector<std::vector<T>> result((size_t) num_vectors);
    long int seed = 44;
    for(auto &r : result) r = getRandomVector<T>(dim1, dim2, seed++);
    return result;
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

struct DryRun{};
struct NoDryRun{};
template<typename use_dry_run = NoDryRun>
int testMethod(int iteratons, std::function<void(int)> test){
    if (std::is_same<use_dry_run, DryRun>::value)
        test(iteratons-1);

    auto time_start = std::chrono::system_clock::now();
    for(int i=0; i<iteratons; i++) test(i);
    auto time_end = std::chrono::system_clock::now();

    long long elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    return int( 0.5 + double(elapsed) / double(iteratons) );
}

#endif
