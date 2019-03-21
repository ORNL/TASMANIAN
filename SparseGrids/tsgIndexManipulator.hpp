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

#include <numeric>
#include <iostream>

#include "tsgIndexSets.hpp"
#include "tsgOneDimensionalWrapper.hpp"
#include "tsgRuleLocalPolynomial.hpp"

/*!
 * \internal
 * \file tsgIndexManipulator.hpp
 * \brief Algorithms for manipulating sets of multi-indexes.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianMultiIndexManipulations
 *
 * A series of templates, lambda, and regular functions that allow the manipulation of
 * multi-indexes and sets of multi-indexes.
 * \endinternal
 */

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

/*!
 * \internal
 * \ingroup TasmanianMultiIndexManipulations
 * \brief Generate the multi-index with entries satisfying the \b inside(), assumes that \b inside() defines a lower-complete set.
 * \endinternal
 */
inline MultiIndexSet generateLowerMultiIndexSet(size_t num_dimensions, std::function<bool(const std::vector<int> &index)> inside){
    size_t c = num_dimensions -1;
    bool is_in = true;
    std::vector<int> root(num_dimensions, 0);
    std::vector<int> indexes;
    while(is_in || (c > 0)){
        if (is_in){
            indexes.insert(indexes.end(), root.begin(), root.end());
            c = num_dimensions-1;
            root[c]++;
        }else{
            std::fill(root.begin() + c, root.end(), 0);
            root[--c]++;
        }
        is_in = inside(root);
    }
    return MultiIndexSet(num_dimensions, indexes);
}

/*!
 * \internal
 * \ingroup TasmanianMultiIndexManipulations
 * \brief Expand the \b set with the minimum number of multi-indexes that will result in a lower complete set.
 * \endinternal
 */
void completeSetToLower(MultiIndexSet &set);

/*!
 * \internal
 * \ingroup TasmanianMultiIndexManipulations
 * \brief Holds anisotropic weights in a usable format.
 *
 * The external API accepts weigths as a single vector of integers or even an empty vector,
 * the \b ProperWeights holds weights that are always defined in a form suitable for computing.
 * \endinternal
 */
struct ProperWeights{
    //! \brief Interpret the single vector and form the proper weights, e.g., use defaults, scale hyperbolic weights, split \b type_ipcurved into linear and curved, etc.
    ProperWeights(size_t num_dimensions, TypeDepth type, std::vector<int> const &weights){
        contour = OneDimensionalMeta::getControurType(type);
        if (weights.empty()) linear = std::vector<int>(num_dimensions, 1);
        if (contour == type_level){ // linear weights only
            if (!weights.empty()){
                linear = weights; // simplest type
            }
        }else if (contour == type_curved){
            if (weights.empty()){
                curved = std::vector<double>(num_dimensions, 0.0); // cannot define default log-correction
            }else{
                linear = std::vector<int>(weights.begin(), weights.begin() + num_dimensions);
                curved = std::vector<double>(num_dimensions);
                std::transform(weights.begin() + num_dimensions, weights.end(), curved.begin(), [&](int i)->double{ return (double) i; });
            }
        }else{ // hyperbolic
            if (weights.empty()){
                curved = std::vector<double>(num_dimensions, 1.0);
            }else{
                linear = std::vector<int>(num_dimensions, 1); // used only for offset normalization
                double exponent_normalization = (double) std::accumulate(weights.begin(), weights.end(), 0);
                curved = std::vector<double>(num_dimensions);
                std::transform(weights.begin(), weights.end(), curved.begin(),
                               [&](int i)->double{ return ((double) i) / exponent_normalization; });
            }
        }
    }
    //! \brief Return \b true if the combination of weights and contour is guaranteed to give a lower-complete set.
    bool provenLower() const{
        if (contour == type_curved)
            for(size_t i=0; i<linear.size(); i++)
                if ((double)linear[i] + curved[i] < 0) return false;
        return true;
    }
    //! \brief Return the smallest linear weight, used for normalization of the offset.
    int minLinear() const{ return *std::min_element(linear.begin(), linear.end()); }
    //! \brief Return the number of dimensions.
    size_t getNumDimensions() const{ return linear.size(); }

    //! \brief Always equal to \b type_level, \b type_curved or \b type_hyperbolic, indicates the general contour and how the weights should be used.
    TypeDepth contour;
    //! \brief The linear components of the weights for \b type_level and \b type_curved contours; a vector of ones for \b type_hyperbolic.
    std::vector<int> linear;
    //! \brief The curved components of the weights \b type_curved contours; a scaled vector of exponents for \b type_hyperbolic; empty for \b type_level.
    std::vector<double> curved;
};

/*!
 * \internal
 * \brief Generate a lower complete multi-index set that satisfies the given properties.
 *
 * Using the \b num_dimensions select all multi-indexes with weight less than the \b offset.
 * If not empty, \b anisotropic_weights are used to form the weight and the offset is scaled by the minimum linear weight.
 * Rule exactness is used in place of the raw-index to cover the correct polynomial space, see references in \b TypeDepth.
 * \endinternal
 */
MultiIndexSet selectTensors(size_t num_dimensions, int offset, TypeDepth type, std::function<int(int i)> rule_exactness,
                            std::vector<int> const &anisotropic_weights, std::vector<int> const &level_limits);

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
Data2D<int> computeDAGup(MultiIndexSet const &mset);

//! \internal
//! \brief Cache the single-indexes of the parents of the multi-indexes (points) in \b mset.
//! \ingroup TasmanianMultiIndexManipulations

//! Each node defined by a multi-index in \b mset can have one or more parents in each direction,
//! where the parent-offspring relation is defined by the \b rule. For each index in \b mset, the
//! Data2D structure \b parents will hold a strip with the location of each parent in \b mset
//! (or -1 if the parent is missing from \b mset).
Data2D<int> computeDAGup(MultiIndexSet const &mset, const BaseRuleLocalPolynomial *rule);

/*!
 * \internal
 * \brief Returns a vector that is the sum of the one dimensional levels of each multi-index in the set.
 * \ingroup TasmanianMultiIndexManipulations
 * \endinternal
 */
std::vector<int> computeLevels(MultiIndexSet const &mset, BaseRuleLocalPolynomial const *rule);

/*!
 * \internal
 * \brief Will call \b apply() with the slot index in \b mset of each parent/child of \b point.
 * \ingroup TasmanianMultiIndexManipulations
 * \endinternal
 */
inline void touchAllImmediateRelatives(std::vector<int> &point, MultiIndexSet const &mset, BaseRuleLocalPolynomial const *rule, std::function<void(int i)> apply){
    int max_kids = rule->getMaxNumKids();
    for(auto &v : point){
        int save = v; // replace one by one each index of p with either parent or kid

        // check the parents
        v = rule->getParent(save);
        if (v > -1){
            int parent_index = mset.getSlot(point);
            if (parent_index > -1)
                apply(parent_index);
        }

        v = rule->getStepParent(save);
        if (v > -1){
            int parent_index = mset.getSlot(point);
            if (parent_index > -1)
                apply(parent_index);
        }

        for(int k=0; k<max_kids; k++){
            v = rule->getKid(save, k);
            if (v > -1){
                int kid_index = mset.getSlot(point);
                if (kid_index > -1)
                    apply(kid_index);
            }
        }

        v = save; // restore the original index for the next iteration
    }
}

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

/*!
 * \internal
 * \ingroup TasmanianMultiIndexManipulations
 * \brief For a set of \b tensors compute the corresponding polynomial space assuming the 1D rules have given \b exactness().
 * \endinternal
 */
MultiIndexSet createPolynomialSpace(const MultiIndexSet &tensors, std::function<int(int)> exactness);

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
    int num_dimensions = (int) tensors.getNumDimensions();
    Data2D<int> tens(num_dimensions, 0);
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
