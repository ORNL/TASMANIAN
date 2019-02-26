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

//! \internal
//! \file tsgIndexManipulator.hpp
//! \brief Algorithms for manipulating sets of multi-indexes.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianMultiIndexManipulations
//!
//! A series of templates, lambda, and regular functions that allow the manipulation of
//! multi-indexes and sets of multi-indexes. The algorithms include the selection of
//! multi-indexes according to a criteria, union of a vector of multi-index sets,
//! normalization of anisotropic weights, map parents and children within a set (which
//! generates a DAG structure), complete a set to satisfy the lower-property, etc.

/*!
 * \internal
 * \ingroup TasmanianSG
 * \addtogroup TasmanianMultiIndexManipulations Multi-Index manipulation algorithms
 *
 * \par Multi-Index manipulation algorithms
 * A series of templates, lambda, and regular functions that allow the manipulation of
 * multi-indexes and sets of multi-indexes. The algorithms include the selection of
 * multi-indexes according to a criteria, union of a vector of multi-index sets,
 * normalization of anisotropic weights, map parents and children within a set (which
 * generates a DAG structure), complete a set to satisfy the lower-property, etc.
 * \endinternal
 */

namespace TasGrid{

namespace MultiIndexManipulations{

//! \internal
//! \brief Generate a full tensor multi-index set with **num_entries**  in each direction
//! \ingroup TasmanianMultiIndexManipulations
template<typename I>
void generateFullTensorSet(const std::vector<I> &num_entries, MultiIndexSet &set){
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
    set.setNumDimensions((int) num_dimensions);
    set.setIndexes(indexes);
}

//! \internal
//! \brief Generate the multi-index with entries satisfying the **!outside()**, assumes that the **!outside()** defines a lower set
//! \ingroup TasmanianMultiIndexManipulations
template<typename I>
void generateLowerMultiIndexSet(std::function<bool(const std::vector<I> &index)> outside, MultiIndexSet &set){
    size_t num_dimensions = set.getNumDimensions();
    size_t c = num_dimensions -1;
    bool out = false;
    std::vector<I> root(num_dimensions, 0);
    std::vector<I> indexes;
    while( !(out && (c == 0)) ){
        if (out){
            for(size_t k=c; k<num_dimensions; k++) root[k] = 0;
            c--;
            root[c]++;
        }else{
            indexes.insert(indexes.end(), root.begin(), root.end());
            c = num_dimensions-1;
            root[c]++;
        }
        out = outside(root);
    }
    set.setIndexes(indexes);
}

//! \internal
//! \brief Take the last set of **level_sets** append a new set of the parents or children that satisfy **inside()**, repeat recursively until there are no more indexes to add
//! \ingroup TasmanianMultiIndexManipulations
template<typename I, bool completion>
void recursiveLoadPoints(std::function<bool(const std::vector<I> &index)> inside, const MultiIndexSet &set, std::vector<MultiIndexSet> &level_sets){
//! \b level_sets must contain at least one set, then this function considers all the children/parents of the entries of the last set (given by `level_sets.back()`)
//! and then creates a new set with only the children that satisfy the \b inside() or the parents (regardless of \b inside()). The new set is appended to the \b level_sets
//! (similar to `push_back()`, but using `resize()`). The recursion terminates at the set where all children fail the \b inside() or all parents are already included.
//! - \b I defines the Int type for the set (should be `int` right now)
//! - \b completion equal \b True indicates the use of *parents* (i.e., do the lower completion), set to \b False indicates the use of *children*
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
                        if (set.missing(point)) level.appendStrip(point);
                        p++;
                    }
                }
            }else{
                for(auto &p : point){
                    p++;
                    if (inside(point)) level.appendStrip(point);
                    p--;
                }
            }
        }

        adding = (level.getNumStrips() > 0);
        if (adding){
            level_sets.resize(level_sets.size() + 1);
            level_sets.back().setNumDimensions((int) num_dimensions);
            level_sets.back().addData2D(level);
        }
    }
}

//! \internal
//! \brief Take the union of all \b level_sets and either \b overwrite the \b set or if \b overwrite is \b False then add to the \b set
//! \ingroup TasmanianMultiIndexManipulations
template<bool overwrite>
void unionSets(std::vector<MultiIndexSet> &level_sets, MultiIndexSet &set){
    int num_levels = (int) level_sets.size();
    while(num_levels > 1){
        int stride = num_levels / 2 + (((num_levels % 2) > 0) ? 1 : 0);
        for(int i=0; i<stride; i++)
            if (i + stride < num_levels)
                level_sets[i].addMultiIndexSet(level_sets[i + stride]);
        num_levels = stride;
    }
    if (overwrite){
        set = std::move(level_sets[0]);
    }else{
        set.addMultiIndexSet(level_sets[0]);
    }
}

//! \internal
//! \brief If **set** is a lower set, then do nothing, otherwise add the the multi-indexes needed to make **set** a lower-set
//! \ingroup TasmanianMultiIndexManipulations
template<typename I>
void completeSetToLower(MultiIndexSet &set){
    size_t num_dimensions = set.getNumDimensions();
    int num = set.getNumIndexes();
    Data2D<I> completion;
    completion.resize((size_t) set.getNumDimensions(), 0);
    for(int i=0; i<num; i++){
        std::vector<I> point(num_dimensions);
        std::copy_n(set.getIndex(i), num_dimensions, point.data());
        for(auto &p : point){
            if (p != 0){
                p--;
                if (set.missing(point)) completion.appendStrip(point);
                p++;
            }
        }
    }

    if (completion.getNumStrips() > 0){
        std::vector<MultiIndexSet> level_sets(1);
        level_sets[0].setNumDimensions((int) num_dimensions);
        level_sets[0].addData2D(completion);

        recursiveLoadPoints<I, true>([](const std::vector<I> &) -> bool{ return true; }, set, level_sets);

        unionSets<false>(level_sets, set);
    }
}

//! \internal
//! \brief Generate the multi-index with entries satisfying the **criteria()**, the result is a lower set regardless of the **criteria()**
//! \ingroup TasmanianMultiIndexManipulations
template<typename I>
void generateGeneralMultiIndexSet(std::function<bool(const std::vector<I> &index)> criteria, MultiIndexSet &set){
    size_t num_dimensions = set.getNumDimensions();
    std::vector<MultiIndexSet> level_sets(1);
    level_sets[0].setNumDimensions((int) num_dimensions);
    std::vector<I> root(num_dimensions, 0);
    level_sets[0].setIndexes(root);

    recursiveLoadPoints<I, false>(criteria, set, level_sets);

    unionSets<true>(level_sets, set);

    completeSetToLower<I>(set);
}

//! \internal
//! \brief Set an empty vector **weights** to the canonical weights of the type (if **anisotropic_weights** is empty) or the actual weights with scaling adjustment
//! \ingroup TasmanianMultiIndexManipulations
template<typename I>
void getProperWeights(size_t num_dimensions, I offset, TypeDepth type, const std::vector<I> &anisotropic_weights,
                      std::vector<I> &weights, I &normalized_offset, bool &known_lower){
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
            normalized_offset = offset * *std::min_element(weights.begin(), weights.begin() + num_dimensions);
            if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
                known_lower = false;
                for(size_t i=0; i<num_dimensions; i++) if (weights[i] + weights[i+num_dimensions] < 0) known_lower = false;
            }else{
                known_lower = true;
            }
        }
    }
}

//! \internal
//! \brief Split a single \b weights vector into two vectors and computes some meta data.
//! \ingroup TasmanianMultiIndexManipulations
template<typename I>
void splitWeights(size_t num_dimensions, TypeDepth type, const std::vector<I> &weights,
                  std::vector<I> &proper_weights, std::vector<double> &curved_weights, double &hyper_denom, TypeDepth &contour_type){
//! Interpretation of the weight goes in two stages, first we decide between total-degree, curved, or hyperbolic types,
//! then the weights are split into linear and logarithmic (curved), and finally scaling factors are computed (hyperbolic only).
//! - \b num_dimensions is the number of dimensions of the problem
//! - \b type defines the contour described by the weights
//! - \b weights has size \b num_dimensions for linear and hyperbolic contours, and 2 times \b num_dimensions for curved contours;
//!   \b weights can be empty, in which case isotropic weights will be used in the outputs
//! - \b proper_weights returns the linear profile (linear and curved) or the unscaled hyperbolic profile
//! - \b curved_weights returns the logarithmic or curved profile, or the hyperbolic weights converged to double-precision numbers
//! - \b hyper_denom is the scaling factor for the hyperbolic weights
//! - \b contour_type reduces \b type to either \b type_level, \b type_curved, or \b type_hyperbolic
//! This function is used when constructing Global and Sequence grids with user specified weights.
    proper_weights = weights;
    contour_type = OneDimensionalMeta::getControurType(type);
    if (contour_type == type_hyperbolic){
        if (proper_weights.empty()){
            curved_weights = std::vector<double>(num_dimensions, 1.0);
            hyper_denom = 1.0;
        }else{
            curved_weights.resize(num_dimensions);
            std::transform(weights.begin(), weights.end(), curved_weights.begin(), [&](int i)->double{ return (double) i; });
            hyper_denom = (double) std::accumulate(weights.begin(), weights.end(), 1);
        }
    }else if (contour_type == type_curved){
        if (proper_weights.empty()){
            proper_weights = std::vector<int>(num_dimensions, 1);
            contour_type = type_level; // canonical curved weights have no curve, switch to linear profile
        }else{
            proper_weights.resize(num_dimensions); // saves the first num_dimensions entries
            curved_weights.resize(num_dimensions);
            std::transform(weights.begin() + num_dimensions, weights.end(), curved_weights.begin(), [&](int i)->double{ return (double) i; });
        }
    }else{
        if (proper_weights.empty()) proper_weights = std::vector<int>(num_dimensions, 1);
    }
}

//! \internal
//! \brief Compute the weight of a multi-index based on the split weights and the contour.
//! \ingroup TasmanianMultiIndexManipulations

//! Here \b exactness is the vector with the exactness values based on the rule and the space type, i.e., level vs curved.
inline double computeMultiIndexWeight(const std::vector<int> &exactness, const std::vector<int> &proper_weights,
                                      const std::vector<double> &curved_weights, double hyper_denom, TypeDepth contour_type){
    if (contour_type == type_level){
        return (double) std::inner_product(exactness.begin(), exactness.end(), proper_weights.data(), 0);
    }else if (contour_type == type_hyperbolic){
        double result = 1.0;
        auto itr = curved_weights.begin();
        for(auto e : exactness)
            result *= pow((double) (1.0 + e), *itr++ / hyper_denom);
        return result;
    }else{
        double result = (double) std::inner_product(exactness.begin(), exactness.end(), proper_weights.data(), 0);
        auto itr = curved_weights.begin();
        for(auto e : exactness)
            result += *itr++ * log1p((double) e);
        return result;
    }
}

//! \internal
//! \brief Generates a multi-index set based on the combination of weights, rule_exactness, offset and type
//! \ingroup TasmanianMultiIndexManipulations
template<typename I, typename CacheType, TypeDepth TotalCurvedHyper>
void generateWeightedTensorsCached(const std::vector<I> &weights, CacheType normalized_offset, std::function<long long(I i)> rule_exactness, MultiIndexSet &mset){
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

//! \internal
//! \brief Generates a multi-index set based on the combination of weights, rule_exactness, offset and type, assume cache has to be build dynamically (non-lower set case)
//! \ingroup TasmanianMultiIndexManipulations
template<typename I, typename CacheType>
void generateWeightedTensorsDynamicCached(const std::vector<I> &weights, CacheType normalized_offset, std::function<long long(I i)> rule_exactness, MultiIndexSet &mset){
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
//! \brief Generate the multi-index set defined by the parameters, **rule** cannot be custom
//! \ingroup TasmanianMultiIndexManipulations
void selectTensors(int offset, TypeDepth type, std::function<long long(int i)> rule_exactness, const std::vector<int> &anisotropic_weights, MultiIndexSet &set);

//! \internal
//! \brief Returns a vector that is the sum of entries of each multi-index in the set
//! \ingroup TasmanianMultiIndexManipulations
void computeLevels(const MultiIndexSet &mset, std::vector<int> &level);

//! \internal
//! \brief Returns a vector that is the maximum index in each direction and a total maximum index, **max_levels** must be empty
//! \ingroup TasmanianMultiIndexManipulations
void getMaxIndex(const MultiIndexSet &mset, std::vector<int> &max_levels, int &total_max);

//! \internal
//! \brief Returns a Data2D structure where each strip holds the indexes of the parents of indexes of mset (for each direction), using one-point-growth hierarchy
//! \ingroup TasmanianMultiIndexManipulations
void computeDAGup(const MultiIndexSet &mset, Data2D<int> &parents);

//! \internal
//! \brief Cache the single-indexes of the parents of the multi-indexes (points) in \b mset.
//! \ingroup TasmanianMultiIndexManipulations

//! Each node defined by a multi-index in \b mset can have one or more parents in each direction,
//! where the parent-offspring relation is defined by the \b rule. For each index in \b mset, the
//! Data2D structure \b parents will hold a strip with the location of each parent in \b mset
//! (or -1 if the parent is missing from \b mset).
void computeDAGup(const MultiIndexSet &mset, const BaseRuleLocalPolynomial *rule, Data2D<int> &parents);

//! \internal
//! \brief Using the **flagged** map, create **new_set** with the flagged children of **mset** but only if they obey the **level_limits**
//! \ingroup TasmanianMultiIndexManipulations
void selectFlaggedChildren(const MultiIndexSet &mset, const std::vector<bool> &flagged, const std::vector<int> &level_limits, MultiIndexSet &new_set);

//! \internal
//! \brief Replace \b mset with the subset of \b mset that obeys the \b level_limits
//! \ingroup TasmanianMultiIndexManipulations
void removeIndexesByLimit(const std::vector<int> &level_limits, MultiIndexSet &mset);

//! \internal
//! \brief Assuming that **tensors** describe a set of nested tensor operators, and each 1-D operators has **getNumPoints()**, then generate the actual points
//! \ingroup TasmanianMultiIndexManipulations

//! Assuming that we are working with a nested rule, then instead of generating the points for all tensor rules and taking the union,
//! it is much faster to generate just the surplus points and union those, i.e., the points associated with the surplus tensor operator.
//! Working with surplus points ensures that there is no repetition of points when merging the sets.
//! * \b tensors is a lower set of multi-indexes describing tensor rules
//! * \b getNumPoints() described the number of points for each 1-D rule
//! * On exit, \b points is the union of the points of all tensors
void generateNestedPoints(const MultiIndexSet &tensors, std::function<int(int)> getNumPoints, MultiIndexSet &points);

//! \internal
//! \brief Assuming that **tensors** describe a set of non-nested tensor operators described by the **wrapper**, then generate the actual **points**
//! \ingroup TasmanianMultiIndexManipulations

//! Assuming that we are working with a non-nested rule, then for each tensor we must generate the points and map them to the global indexing,
//! then take the union of all the tensors
//! * **tensors** is a set of tensor rules, non-necessarily lower
//! * **wrapper** described the one dimensional rules (most notably the level-order-to-global-index mapping)
//! * **points** is the union of the points of all tensors
void generateNonNestedPoints(const MultiIndexSet &tensors, const OneDimensionalWrapper &wrapper, MultiIndexSet &points);

//! \internal
//! \brief Given a tensor defined by **levels** find the references to all tensor points in the **points** set (assuming the standard order of the tensor entries)
//! \ingroup TasmanianMultiIndexManipulations
template<bool nested>
void referencePoints(const int levels[], const OneDimensionalWrapper &wrapper, const MultiIndexSet &points, std::vector<int> &refs){
    size_t num_dimensions = (size_t) points.getNumDimensions();
    std::vector<int> num_points(num_dimensions);
    int num_total = 1; // this will be a subset of all points, no danger of overflow
    for(size_t j=0; j<num_dimensions; j++) num_points[j] = wrapper.getNumPoints(levels[j]);
    for(auto n : num_points) num_total *= n;

    refs.resize(num_total);
    std::vector<int> p(num_dimensions);

    for(int i=0; i<num_total; i++){
        int t = i;
        auto n = num_points.rbegin();
        for(int j=(int) num_dimensions-1; j>=0; j--){
            p[j] = (nested) ? t % *n : wrapper.getPointIndex(levels[j], t % *n);
            t /= *n++;
        }
        refs[i] = points.getSlot(p);
    }
}

//! \internal
//! \brief Computes the weights for the tensor linear combination, **mset** is a lower multi-index set and **weight** is resized
//! \ingroup TasmanianMultiIndexManipulations
void computeTensorWeights(const MultiIndexSet &mset, std::vector<int> &weights);

//! \internal
//! \brief On exit, **active** is the set containing the multi-indexes of **mset** corresponding to non-zero **weights**
//! \ingroup TasmanianMultiIndexManipulations
void createActiveTensors(const MultiIndexSet &mset, const std::vector<int> &weights, MultiIndexSet &active);

//! \internal
//! \brief For a set of **tensors** compute the corresponding polynomial **space** assuming the 1D rules have given **exactness**
//! \ingroup TasmanianMultiIndexManipulations
void createPolynomialSpace(const MultiIndexSet &tensors, std::function<int(int)> exactness, MultiIndexSet &space);

//! \internal
//! \brief Assuming that \b mset is lower complete, return \b true if adding the \b point will preserve completeness.
//! \ingroup TasmanianMultiIndexManipulations
template<typename I> bool isLowerComplete(const std::vector<I> &point, const MultiIndexSet &mset){
    std::vector<I> dad = point;
    for(auto &d : dad){
        if (d > 0){
            d--;
            if (mset.missing(dad)) return false;
            d++;
        }
    }
    return true;
}

//! \internal
//! \brief For a set of \b tensors create an \b mset that contain the children of indexes in \b tensors that are missing from \b exclude and obey the \b level_limits.
//! \ingroup TasmanianMultiIndexManipulations
template<bool limited>
void addExclusiveChildren(const MultiIndexSet &tensors, const MultiIndexSet &exclude, const std::vector<int> level_limits, MultiIndexSet &mset){
    int num_dimensions = tensors.getNumDimensions();
    Data2D<int> tens;
    tens.resize(num_dimensions, 0);
    for(int i=0; i<tensors.getNumIndexes(); i++){ // add the new tensors (so long as they are not included in the initial grid)
        const int *t = tensors.getIndex(i);
        std::vector<int> kid(t, t + num_dimensions);
        auto ilimit = level_limits.begin();
        for(auto &k : kid){
            k++;
            if (exclude.missing(kid) && tensors.missing(kid)){ // if the kid is not to be excluded and if not included in the current set
                if (isLowerComplete(kid, tensors)){
                    if (limited){
                        if ((*ilimit == -1) || (k <= *ilimit))
                            tens.appendStrip(kid);
                        ilimit++;
                    }else{
                        tens.appendStrip(kid);
                    }
                }
            }
            k--;
        }
    }
    mset.setNumDimensions(num_dimensions); // set of new tensors
    if (tens.getVector().size() > 0) mset.addData2D(tens);
}

}

}

#endif
