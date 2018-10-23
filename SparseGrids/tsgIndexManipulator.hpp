/*
 * Copyright (c) 2017, Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
 *    and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
 * THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
 * THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
 * IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TSG_INDEX_MANIPULATOR_HPP
#define __TSG_INDEX_MANIPULATOR_HPP

#include <functional>
#include <numeric>

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgOneDimensionalWrapper.hpp"
#include "tsgRuleLocalPolynomial.hpp"

namespace TasGrid{

namespace MultiIndexManipulations{

template<typename I>
void generateFullTensorSet(const std::vector<I> &num_entries, MultiIndexSet &set){
//! \internal
//! \brief generate a full tensor multi-index set with **num_entries**  in each direction
//! \ingroup TasmanianMultiIndexManipulations
    I num_total = 1;
    for(auto &l : num_entries) num_total *= l;
    size_t num_dimensions = num_entries.size();
    std::vector<I> indexes(((size_t) num_total) * num_dimensions);
    auto iter = indexes.rbegin();
    for(I i=num_total-1; i>=0; i--){
        I t = i;
        auto l = num_entries.rbegin();
        // in order to generate indexes in the correct order, the for loop must go backwards
        for(size_t j = 0; j<num_dimensions; j++){
            *iter++ = (I) (t % *l);
            t /= *l++;
        }
    }
    set.setIndexes(indexes);
}

template<typename I>
void generateLowerMultiIndexSet(std::function<bool(const std::vector<I> &index)> criteria, MultiIndexSet &set){
//! \internal
//! \brief generate the multi-index with entries satisfying the **criteria()**, assumes that the **criteria()** defines a lower set
//! \ingroup TasmanianMultiIndexManipulations
    size_t num_dimensions = set.getNumDimensions();
    size_t c = num_dimensions -1;
    bool outside = false;
    std::vector<I> root(num_dimensions, 0);
    std::vector<I> indexes;
    while( !(outside && (c == 0)) ){
        if (outside){
            for(size_t k=c; k<num_dimensions; k++) root[k] = 0;
            c--;
            root[c]++;
        }else{
            indexes.insert(indexes.end(), root.begin(), root.end());
            c = num_dimensions-1;
            root[c]++;
        }
        outside = criteria(root);
    }
    set.setIndexes(indexes);
}

template<typename I, bool completion>
void recursiveLoadPoints(std::function<bool(const std::vector<I> &index)> criteria, const MultiIndexSet &set, std::vector<MultiIndexSet> &level_sets){
//! \internal
//! \brief take the last set of **level_sets** append a new set of the parents or children that satisfy **criteria()**, repeat recursively until there are no more indexes to add
//! \ingroup TasmanianMultiIndexManipulations
//!
//! **level_sets** must contain at least one set, then this function considers all the children/parents of the entries of the last set (given by `level_sets.back()`)
//! and then creates a new set with only the children/parents that satisfy the **criteria()**. The new set is appended to the **level_sets**
//! (similar to `push_back()`, but using `resize()`). The recursion terminates at the set where all children/parents fail the **criteria()**.
//! * **I** defines the Int type for the set (should be `int` right now)
//! * **completion** equal `true` indicates the use of *parents* (i.e., do the lower completion), set to `false` indicates the use of *children*
    size_t num_dimensions = level_sets.back().getNumDimensions();
    bool adding = true;
    while(adding){
        Data2D<I> level;
        level.resize((int) num_dimensions, 0); // might be a problem if we use "getStrip()"
        int num_indexes = level_sets.back().getNumIndexes();
        for(int i=0; i<num_indexes; i++){
            std::vector<I> point(num_dimensions);
            std::copy_n(level_sets.back().getIndex(i), num_dimensions, point.data());
            if (completion){
                for(auto &p : point){
                    if (p != 0){
                        p--;
                        if (set.getSlot(point) == -1) level.appendStrip(point);
                        p++;
                    }
                }
            }else{
                for(auto &p : point){
                    p++;
                    if (criteria(point)) level.appendStrip(point);
                    p--;
                }
            }
        }

        adding = (level.getNumStrips() > 0);
        if (adding){
            level_sets.resize(level_sets.size() + 1);
            level_sets.back().setNumDimensions((int) num_dimensions);
            level_sets.back().addUnsortedInsexes(*level.getVector());
        }
    }
}

template<bool overwrite>
void unionSets(std::vector<MultiIndexSet> &level_sets, MultiIndexSet &set){
//! \internal
//! \brief **set** is combines with the union of the **level_sets**, if **overwrite** is **true** the content of **set** is discarded
//! \ingroup TasmanianMultiIndexManipulations
    int num_levels = (int) level_sets.size();
    while(num_levels > 2){
        int stride = num_levels / 2 + (((num_levels % 2) > 0) ? 1 : 0);
        for(int i=0; i<stride; i++){
            if (i + stride < num_levels) level_sets[i].addSortedInsexes(*level_sets[i + stride].getVector());
        }
        num_levels /= 2;
    }
    if (overwrite){
        set.move(level_sets[0]);
    }else{
        set.addSortedInsexes(*level_sets[0].getVector());
    }
    if (num_levels == 2) set.addSortedInsexes(*level_sets[1].getVector());
}

template<typename I>
void generateGeneralMultiIndexSet(std::function<bool(const std::vector<I> &index)> criteria, MultiIndexSet &set){
//! \internal
//! \brief generate the multi-index with entries satisfying the **criteria()**, the result is a lower set regardless of the **criteria()**
//! \ingroup TasmanianMultiIndexManipulations
    size_t num_dimensions = set.getNumDimensions();
    std::vector<MultiIndexSet> level_sets(1);
    level_sets[0].setNumDimensions((int) num_dimensions);
    std::vector<I> root(num_dimensions, 0);
    level_sets[0].setIndexes(root);

    recursiveLoadPoints<I, false>(criteria, set, level_sets);

    unionSets<true>(level_sets, set);

    int num = set.getNumIndexes();
    Data2D<I> completion;
    for(int i=0; i<num; i++){
        std::vector<I> point(num_dimensions);
        std::copy_n(set.getIndex(i), num_dimensions, point.data());
        for(auto &p : point){
            if (p != 0){
                p--;
                if (set.getSlot(point) != -1) completion.appendStrip(point);
                p++;
            }
        }
    }

    if (completion.getNumStrips() > 0){
        level_sets.resize(1);
        level_sets[0].reset();
        level_sets[0].setNumDimensions((int) num_dimensions);
        level_sets[0].addUnsortedInsexes(*completion.getVector());

        recursiveLoadPoints<I, true>(criteria, set, level_sets);

        unionSets<false>(level_sets, set);
    }
}

template<typename I>
void getProperWeights(size_t num_dimensions, I offset, TypeDepth type, const std::vector<I> &anisotropic_weights,
                      std::vector<I> &weights, I &normalized_offset, bool &known_lower){
//! \internal
//! \brief set an empty vector **weights** to the canonical weights of the type (if **anisotropic_weights** is empty) or the actual weights with scaling adjustment
//! \ingroup TasmanianMultiIndexManipulations
//!
//! Based on the **anisotropic_weights** and selection parameters, create canonical/normalized weights and run basic analysis on the resulting multi-index set.
//! * **num_dimensions** specifies the number of dimensions of the multi-index
//! * **type** is the selection type, e.g., type_iptotal
//! * **anisotropic_weights** is either an empty vector or a vector containing the user specified (non-normalized) weights
//! * *weights* contains either the normalized weights or the canonical normalized weights
//! * **normalized_offset** returns the normalized selection threshold
//! * **known_lower** returns true if the combination of weights and type will guarantee a lower multi-index set
    if ((type == type_tensor) || (type == type_iptensor) || (type == type_qptensor)){ // special case, full tensor
        weights.resize(num_dimensions, offset);
        if (!anisotropic_weights.empty()) for(size_t j=0; j<num_dimensions; j++) weights[j] *= anisotropic_weights[j];
        normalized_offset = 0; // dummy value, not used
        known_lower = true; // full-tensors are always lower
    }else{
        if (anisotropic_weights.empty()){
            if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
                weights.resize(2*num_dimensions);
                std::fill_n(weights.begin(), num_dimensions, 1);
                std::fill_n(&(weights[num_dimensions]), num_dimensions, 0);
            }else{
                weights.resize(num_dimensions, 1);
            }
            normalized_offset = offset;
            known_lower = true;
        }else{
            weights = anisotropic_weights; // copy assign
            normalized_offset = offset * *std::max_element(weights.begin(), weights.begin() + num_dimensions);
            if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
                known_lower = false;
                for(size_t i=0; i<num_dimensions; i++) if (weights[i] + weights[i+num_dimensions] < 0) known_lower = false;
            }else{
                known_lower = true;
            }
        }
    }
}

template<typename I, typename CacheType, TypeDepth TotalCurvedHyper>
void generateWeightedTensorsCached(const std::vector<I> &weights, CacheType normalized_offset, std::function<long long(I i)> rule_exactness, MultiIndexSet &mset){
//! \internal
//! \brief generates a multi-index set based on the combination of weights, rule_exactness, offset and type
//! \ingroup TasmanianMultiIndexManipulations
//!
//! Simplifies the calls to caching and generating multi-indexes
//! * **I** is the Int type of multi-index (use int)
//! * **CacheType** should wither match **I** (for level, iptotal and qptotal types) or **double** for curved and hyperbolic rules
//! * **TotalCurvedHyper** should be either `type_level`, `type_curved` or `type_hyperbolic`, indicating the selection strategy
//! * **rule_exactness()** handles the cases of quadrature, interpolation or level based selection
//! * **mset** is the resulting multi-index set
    size_t num_dimension = (size_t) mset.getNumDimensions();
    std::vector<std::vector<CacheType>> cache(num_dimension);
    for(size_t j=0; j<num_dimension; j++){
        CacheType xi = (CacheType) weights[j]; // anisotropic weight for this direction
        CacheType eta;
        if (TotalCurvedHyper == type_curved){
            eta = (CacheType) weights[j + num_dimension]; //curved correction
        }else if (TotalCurvedHyper == type_hyperbolic){
            eta = 1;
            if (std::any_of(weights.begin(), weights.end(), [](I w)->bool{ return (w != 1); })) // if not using canonical weights
                for(auto w : weights) eta += (CacheType) w;
        }
        I i = 0;
        CacheType weight1d = (TotalCurvedHyper == type_hyperbolic) ? 1 : 0;
        cache[j].push_back(weight1d);
        do{
            i++;
            long long exactness = 1 + rule_exactness(i - 1);
            if (TotalCurvedHyper == type_level){
                weight1d = ((CacheType) exactness) * xi;
            }else if (TotalCurvedHyper == type_curved){
                weight1d = ((CacheType) exactness) * xi + eta * log1p((CacheType) exactness);
            }else{
                weight1d = pow((CacheType) (1 + exactness), xi / eta);
                std::cout << weight1d << std::endl;
            }
            cache[j].push_back(weight1d);
        }while(ceil(weight1d) <= normalized_offset);
    }
    generateLowerMultiIndexSet<I>([&](const std::vector<I> &index) -> bool{
                                        CacheType w = 0;
                                        auto i = index.begin();
                                        if (TotalCurvedHyper == type_hyperbolic){
                                            w = 1;
                                            for(const auto &v : cache) w *= v[*i++];
                                        }else{
                                            for(const auto &v : cache) w += v[*i++];
                                        }
                                        return (ceil(w) > normalized_offset);
                                    }, mset);
}

template<typename I, typename CacheType>
void generateWeightedTensorsDynamicCached(const std::vector<I> &weights, CacheType normalized_offset, std::function<long long(I i)> rule_exactness, MultiIndexSet &mset){
//! \internal
//! \brief generates a multi-index set based on the combination of weights, rule_exactness, offset and type, assume cache has to be build dynamically (non-lower set case)
//! \ingroup TasmanianMultiIndexManipulations
//!
//! Simplifies the calls to caching and generating multi-indexes
//! * **I** is the Int type of multi-index (use int)
//! * **CacheType** should wither match **I** (for level, iptotal and qptotal types) or **double** for curved and hyperbolic rules
//! * **rule_exactness()** handles the cases of quadrature, interpolation or level based selection
//! * **mset** is the resulting multi-index set
    size_t num_dimension = (size_t) mset.getNumDimensions();
    std::vector<std::vector<CacheType>> cache(num_dimension);
    for(size_t j=0; j<num_dimension; j++) cache[j].push_back(0);
    generateGeneralMultiIndexSet<I>([&](const std::vector<I> &index) -> bool{
                                            CacheType w = 0;
                                            for(size_t j=0; j<num_dimension; j++){
                                                while(index[j] >= (I) cache[j].size()){
                                                    CacheType xi = (CacheType) weights[j]; // anisotropic weight for this direction
                                                    CacheType eta = (CacheType) weights[j + num_dimension];
                                                    I i = (I) cache[j].size();
                                                    long long exactness = 1 + rule_exactness(i - 1);
                                                    cache[j].push_back( ((CacheType) exactness) * xi + eta * log1p((CacheType) exactness) );
                                                }
                                                w += cache[j][index[j]];
                                            }
                                            return (ceil(w) <= normalized_offset);
                                        }, mset);
}

//! \internal
//! \brief generate the multi-index set defined by the parameters, **rule** cannot be custom
//! \ingroup TasmanianMultiIndexManipulations
void selectTensors(int offset, TypeDepth type, std::function<long long(int i)> rule_exactness, const std::vector<int> &anisotropic_weights, MultiIndexSet &set);

//! \internal
//! \brief returns a vector that is the sum of entries of each multi-index in the set
//! \ingroup TasmanianMultiIndexManipulations
void computeLevels(const MultiIndexSet &mset, std::vector<int> &level);

//! \internal
//! \brief returns a vector that is the maximum index in each direction and a total maximum index, **max_levels** must be empty
//! \ingroup TasmanianMultiIndexManipulations
void getMaxIndex(const MultiIndexSet &mset, std::vector<int> &max_levels, int &total_max);

//! \internal
//! \brief returns a Data2D structure where each strip holds the indexes of all parents (for all directions) associates with the index
//! \ingroup TasmanianMultiIndexManipulations
void computeDAGup(const MultiIndexSet &mset, Data2D<int> &parents);
}


class IndexManipulator{
public:
    IndexManipulator(int cnum_dimensions, const CustomTabulated* custom = 0);
    ~IndexManipulator();

    IndexSet* selectTensors(int offset, TypeDepth type, const std::vector<int> &anisotropic_weights, TypeOneDRule rule) const;
    IndexSet* selectTensors(const IndexSet *target_space, bool integration, TypeOneDRule rule) const;
    IndexSet* getLowerCompletion(const IndexSet *iset) const;

    void computeLevels(const IndexSet* iset, std::vector<int> &level) const; // returns a vector of its corresponding to the sum of entries of set

    void makeTensorWeights(const IndexSet* iset, std::vector<int> &weights) const;
    IndexSet* nonzeroSubset(const IndexSet* iset, const std::vector<int> &weights) const; // returns a subset corresponding to the non-zeros weights

    UnsortedIndexSet* tensorGenericPoints(const int levels[], const OneDimensionalWrapper *rule) const;
    IndexSet* generateGenericPoints(const IndexSet *tensors, const OneDimensionalWrapper *rule) const;

    IndexSet* removeIndexesByLimit(IndexSet *iset, const std::vector<int> &limits) const;

    void computeDAGup(const IndexSet *iset, Data2D<int> &parents) const;

    IndexSet* tensorNestedPoints(const int levels[], const OneDimensionalWrapper *rule) const;
    IndexSet* generateNestedPoints(const IndexSet *tensors, const OneDimensionalWrapper *rule) const;

    template<bool nested> void referencePoints(const int levels[], const OneDimensionalWrapper *rule, const IndexSet *points, std::vector<int> &refs) const{
        std::vector<int> num_points(num_dimensions);
        int num_total = 1; // this will be a subset of all points, no danger of overflow
        for(int j=0; j<num_dimensions; j++) num_points[j] = rule->getNumPoints(levels[j]);
        for(auto n : num_points) num_total *= n;

        refs.resize(num_total);
        std::vector<int> p(num_dimensions);

        for(int i=0; i<num_total; i++){
            int t = i;
            auto n = num_points.rbegin();
            for(int j=num_dimensions-1; j>=0; j--){
                p[j] = (nested) ? t % *n : rule->getPointIndex(levels[j], t % *n);
                t /= *n++;
            }
            refs[i] = points->getSlot(p);
        }
    }

    IndexSet* getPolynomialSpace(const IndexSet *tensors, TypeOneDRule rule, bool iexact) const;

    int getMinChildLevel(const IndexSet *iset, TypeDepth type, const std::vector<int> &weights, TypeOneDRule rule);
    // find the minimum level of a child of iset

    IndexSet* selectFlaggedChildren(const IndexSet *iset, const std::vector<bool> &flagged, const std::vector<int> &level_limits) const;

    void getMaxLevels(const IndexSet *iset, std::vector<int> &max_levels, int &total_max) const;
    int getMaxLevel(const IndexSet *iset) const;

    // use by Local Grids
    IndexSet* generatePointsFromDeltas(const IndexSet* deltas, std::function<int(int)> getNumPoints) const;

    void computeDAGupLocal(const IndexSet *iset, const BaseRuleLocalPolynomial *rule, Data2D<int> &parents) const;

protected:
    void getProperWeights(TypeDepth type, const std::vector<int> &anisotropic_weights, std::vector<int> &weights) const;

    template<TypeDepth type>
    long long getIndexWeight(const std::vector<int> &index, const std::vector<int> &weights, TypeOneDRule rule) const{
        long long l = ((type == type_hyperbolic) || (type == type_iphyperbolic) || (type == type_qphyperbolic)) ? 1 : 0;
        if (type == type_level){
            auto witer = weights.begin();
            for(auto i : index) l += i * *witer++;
        }else if (type == type_curved){
            auto walpha = weights.begin();
            auto wbeta = weights.begin(); std::advance(wbeta, num_dimensions);
            double c = 0.0;
            for(auto i : index){
                l += i * *walpha++;
                c += log1p((double) i) * *wbeta++;
            }
            l += (long long) ceil(c);
        }else if ((type == type_iptotal) || (type == type_qptotal)){
            auto witer = weights.begin();
            for(auto i : index) l += ((i > 0) ? (1 + ((type == type_iptotal) ? meta.getIExact(i-1, rule) : meta.getQExact(i-1, rule))) : 0) * *witer++;
        }else if ((type == type_ipcurved) || (type == type_qpcurved)){
            auto walpha = weights.begin();
            auto wbeta = weights.begin(); std::advance(wbeta, num_dimensions);
            double c = 0.0;
            for(auto i : index){
                long long pex = (i > 0) ? (1 + ((type == type_ipcurved) ? meta.getIExact(i-1, rule) : meta.getQExact(i-1, rule))) : 0;
                l += pex * *walpha++;
                c += log1p((double) pex) * *wbeta++;
            }
            l += (long long) ceil(c);
        }else if (type == type_iphyperbolic){
            double nweight = (double) weights[num_dimensions];
            double s = 1.0;
            auto witer = weights.begin();
            for(auto i : index) s *= pow((double) (i+1), ((double) *witer++) / nweight);
            l = (long long) ceil(s);
        }else{
            double nweight = (double) weights[num_dimensions];
            double s = 1.0;
            auto witer = weights.begin();
            for(auto i : index) s *= pow((i>0) ? (2.0 + (double) ((type == type_iphyperbolic) ? meta.getIExact(i-1, rule) : meta.getQExact(i-1, rule))) : 1.0,
                                         ((double) *witer++) / nweight);
            l = (long long) ceil(s);
        }
        return l;
    }

    long long getIndexWeight(const std::vector<int> &index, TypeDepth type, const std::vector<int> &weights, TypeOneDRule rule) const;

private:
    int num_dimensions;

    OneDimensionalMeta meta;
};

}

#endif
