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

#ifndef __TSG_HIERARCHY_MANIPULATOR_HPP
#define __TSG_HIERARCHY_MANIPULATOR_HPP

#include "tsgIndexSets.hpp"
#include "tsgRuleLocalPolynomial.hpp"

/*!
 * \internal
 * \file tsgHierarchyManipulator.hpp
 * \brief Algorithms for manipulating multi-indexes defined by hierarchy rules.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianHierarchyManipulations
 *
 * A series of templates, lambda, and regular functions that allow the manipulation of
 * multi-indexes and sets of multi-indexes defined by hierarchy rules.
 * \endinternal
 */

/*!
 * \internal
 * \ingroup TasmanianSG
 * \addtogroup TasmanianHierarchyManipulations Hierarchical multi-Index manipulation algorithms
 *
 * \par One Dimensional Hierarchy
 * The construction of Global, Sequence and Fourier grids the multi-index hierarchy is associated
 * with tensors and has very simple structure, in one dimension there is only one parent and
 * one child indicated by the previous and next index.
 * The Local-Polynomial and Wavelet girds use more complex hierarchies that have multiple
 * children and (in some cases) parents.
 * Many of the associated manipulation algorithms require a \b rule object that
 * describes the hierarchy, and additional algorithms are needed to handle the more complex
 * parent-child relations.
 *
 * \endinternal
 */

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianHierarchyManipulations
 * \brief Collection of algorithm to manipulate multi-indexes.
 */
namespace HierarchyManipulations{

/*!
 * \internal
 * \ingroup TasmanianHierarchyManipulations
 * \brief Cache the indexes slot numbers of the parents of the multi-indexes in \b mset.
 *
 * Each node defined by a multi-index in \b mset can have one or more parents in each direction,
 * where the parent-offspring relation is defined by the \b rule. For each index in \b mset, the
 * Data2D structure \b parents will hold a strip with the location of each parent in \b mset
 * (or -1 if the parent is missing from \b mset).
 * \endinternal
 */
Data2D<int> computeDAGup(MultiIndexSet const &mset, const BaseRuleLocalPolynomial *rule);

/*!
 * \internal
 * \ingroup TasmanianHierarchyManipulations
 * \brief Cache the indexes slot numbers of the children of the multi-indexes in \b mset.
 *
 * Each node defined by a multi-index in \b mset can have one or more children in each direction,
 * where the parent-offspring relation is defined by the \b rule. For each index in \b mset, the
 * returned Data2D structure will hold a strip with the location of each child in \b mset
 * (or -1 if the kid is missing from \b mset).
 * \endinternal
 */
Data2D<int> computeDAGDown(MultiIndexSet const &mset, const BaseRuleLocalPolynomial *rule);

/*!
 * \internal
 * \brief Returns a vector that is the sum of the one dimensional levels of each multi-index in the set.
 * \ingroup TasmanianHierarchyManipulations
 * \endinternal
 */
std::vector<int> computeLevels(MultiIndexSet const &mset, BaseRuleLocalPolynomial const *rule);

/*!
 * \internal
 * \brief Will call \b apply() with the slot index in \b mset of each parent/child of \b point.
 * \ingroup TasmanianHierarchyManipulations
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

/*!
 * \internal
 * \ingroup TasmanianHierarchyManipulations
 * \brief Return the largest subset of \b candidates such that adding it to \b current will result in a connected graph.
 *
 * \endinternal
 */
inline MultiIndexSet getLargestConnected(MultiIndexSet const &current, MultiIndexSet const &candidates, BaseRuleLocalPolynomial const *rule){
    if (candidates.empty()) return MultiIndexSet();
    auto num_dimensions = candidates.getNumDimensions();
    MultiIndexSet result; // start with an empty set
    if (current.empty()){
        if (candidates.missing(std::vector<int>(num_dimensions, 0))){
            return MultiIndexSet(); // current is empty, 0-th index is the only thing that can be added
        }else{
            result = MultiIndexSet(num_dimensions, std::vector<int>(num_dimensions, 0));
        }
    }
    int max_kids      = rule->getMaxNumKids();
    int max_relatives = rule->getMaxNumParents() + max_kids;
    bool loopon = true;
    while(loopon){
        Data2D<int> update(num_dimensions, 0);

        MultiIndexSet total = current;
        if (!result.empty()) total.addMultiIndexSet(result);
        for(int i=0; i<total.getNumIndexes(); i++){
            std::vector<int> relative(total.getIndex(i), total.getIndex(i) + num_dimensions);
            for(auto &r : relative){
                int k = r; // save the value
                for(int j=0; j<max_relatives; j++){
                    r = (j < max_kids) ? rule->getKid(k, j) : ((j - max_kids == 0) ? rule->getParent(k) : rule->getStepParent(k));
                    if ((r != -1) && !candidates.missing(relative) && total.missing(relative))
                        update.appendStrip(relative);
                }
                r = k;
            }
        }

        loopon = (update.getNumStrips() > 0);
        if (loopon) result.addMultiIndexSet(MultiIndexSet(update));
    }
    return result;
}

/*!
 * \internal
 * \ingroup TasmanianHierarchyManipulations
 * \brief Split the \b data into strips with given \b stride and returns \b Data2D structures grouped by \b levels, preserves the order.
 *
 * \endinternal
 */
template<typename T>
std::vector<Data2D<T>> splitByLevels(size_t stride, std::vector<T> const &data, std::vector<int> const &levels){
    size_t top_level = (size_t) *std::max_element(levels.begin(), levels.end());

    std::vector<Data2D<T>> split(top_level + 1, Data2D<T>(stride, 0));

    for( struct {int i; typename std::vector<T>::const_iterator idata;} v = {0, data.begin()};
         v.idata != data.end();
         v.i++, std::advance(v.idata, stride))
        split[levels[v.i]].appendStrip(v.idata);

    return split;
}

/*!
 * \internal
 * \ingroup TasmanianHierarchyManipulations
 * \brief Reorganize the \b points into sets of nodes that align in one-dimension, used for directional localp refinement.
 *
 * \endinternal
 */
class SplitDirections{
public:
    //! \brief Constructor, deinfe the set to split into directions.
    SplitDirections(const MultiIndexSet &points);
    //! \brief Destroy all data.
    ~SplitDirections(){}

    //! \brief Returns the number of one dimensional jobs.
    int getNumJobs() const{ return (int) job_pnts.size(); }
    //! \brief Return the direction for the \b job.
    int getJobDirection(int job) const{ return job_directions[job]; }
    //! \brief Return the number of points associated with the \b job.
    int getJobNumPoints(int job) const{ return (int) job_pnts[job].size(); }
    //! \brief Return the indexes of the points associated with the job.
    const int* getJobPoints(int job) const{ return job_pnts[job].data(); }

private:
    std::vector<int> job_directions;
    std::vector<std::vector<int>> job_pnts;
};

}

}

#endif
