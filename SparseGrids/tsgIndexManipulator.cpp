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

#ifndef __TSG_INDEX_MANIPULATOR_CPP
#define __TSG_INDEX_MANIPULATOR_CPP

#include "tsgIndexManipulator.hpp"

namespace TasGrid{

void MultiIndexManipulations::selectTensors(int offset, TypeDepth type, std::function<long long(int i)> rule_exactness, const std::vector<int> &anisotropic_weights, MultiIndexSet &mset){
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    std::vector<int> weights;
    int normalized_offset;
    bool known_lower;

    MultiIndexManipulations::getProperWeights<int>(num_dimensions, offset, type, anisotropic_weights, weights, normalized_offset, known_lower);

    if ((type == type_tensor) || (type == type_iptensor) || (type == type_qptensor)){ // special case, full tensor
        std::vector<int> levels(num_dimensions, 0); // how many levels to keep in each direction
        std::transform(weights.begin(), weights.end(), levels.begin(),
                        [&](int w)-> int{
                            int l = 0;
                            while(rule_exactness(l) < w) l++; // get the first level that covers the weight
                            return l+1;
                        });
        generateFullTensorSet<int>(levels, mset);
    }else if (known_lower){ // using lower-set logic (cache can be pre-computed)
        if ((type == type_level) || (type == type_iptotal) || (type == type_qptotal)){
            generateWeightedTensorsCached<int, int, type_level>(weights, normalized_offset, rule_exactness, mset);
        }else if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
            generateWeightedTensorsCached<int, double, type_curved>(weights, (double) normalized_offset, rule_exactness, mset);
        }else{ // some type_hyperbolic type
            generateWeightedTensorsCached<int, double, type_hyperbolic>(weights, (double) normalized_offset, rule_exactness, mset);
        }
    }else{
        generateWeightedTensorsDynamicCached<int, double>(weights, (double) normalized_offset, rule_exactness, mset);
    }
}

void MultiIndexManipulations::computeLevels(const MultiIndexSet &mset, std::vector<int> &level){
    int num_indexes = mset.getNumIndexes();
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    level.resize((size_t) num_indexes);
    #pragma omp parallel for
    for(int i=0; i<num_indexes; i++){
        const int* p = mset.getIndex(i);
        level[i] = std::accumulate(p, p + num_dimensions, 0);
    }
}

void MultiIndexManipulations::getMaxIndex(const MultiIndexSet &mset, std::vector<int> &max_levels, int &total_max){
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    max_levels.resize(num_dimensions, 0);
    int n = mset.getNumIndexes();
    for(int i=0; i<n; i++){
        const int* t = mset.getIndex(i);
        for(size_t j=0; j<num_dimensions; j++) if (max_levels[j] < t[j]) max_levels[j] = t[j];
    }
    total_max = *std::max_element(max_levels.begin(), max_levels.end());
}

void MultiIndexManipulations::computeDAGup(const MultiIndexSet &mset, Data2D<int> &parents){
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    int n = mset.getNumIndexes();
    parents.resize(mset.getNumDimensions(), n);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<n; i++){
        std::vector<int> dad(num_dimensions);
        std::copy_n(mset.getIndex(i), num_dimensions, dad.data());
        int *v = parents.getStrip(i);
        for(auto &d : dad){
            d--;
            *v = (d < 0) ? -1 : mset.getSlot(dad);
            d++;
            v++;
        }
    }
}

void MultiIndexManipulations::computeDAGup(const MultiIndexSet &mset, const BaseRuleLocalPolynomial *rule, Data2D<int> &parents){
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    int num_points = mset.getNumIndexes();
    if (rule->getMaxNumParents() > 1){ // allow for multiple parents and level 0 may have more than one node
        int max_parents = rule->getMaxNumParents() * (int) num_dimensions;
        parents.resize(max_parents, num_points, -1);
        int level0_offset = rule->getNumPoints(0);
        #pragma omp parallel for schedule(static)
        for(int i=0; i<num_points; i++){
            const int *p = mset.getIndex(i);
            std::vector<int> dad(num_dimensions);
            std::copy_n(p, num_dimensions, dad.data());
            int *pp = parents.getStrip(i);
            for(size_t j=0; j<num_dimensions; j++){
                if (dad[j] >= level0_offset){
                    int current = p[j];
                    dad[j] = rule->getParent(current);
                    pp[2*j] = mset.getSlot(dad);
                    while ((dad[j] >= level0_offset) && (pp[2*j] == -1)){
                        current = dad[j];
                        dad[j] = rule->getParent(current);
                        pp[2*j] = mset.getSlot(dad);
                    }
                    dad[j] = rule->getStepParent(current);
                    if (dad[j] != -1){
                        pp[2*j + 1] = mset.getSlot(dad);
                    }
                    dad[j] = p[j];
                }
            }
        }
    }else{ // this assumes that level zero has only one node
        parents.resize((int) num_dimensions, num_points);
        #pragma omp parallel for schedule(static)
        for(int i=0; i<num_points; i++){
            const int *p = mset.getIndex(i);
            std::vector<int> dad(num_dimensions);
            std::copy_n(p, num_dimensions, dad.data());
            int *pp = parents.getStrip(i);
            for(size_t j=0; j<num_dimensions; j++){
                if (dad[j] == 0){
                    pp[j] = -1;
                }else{
                    dad[j] = rule->getParent(dad[j]);
                    pp[j] = mset.getSlot(dad.data());
                    while((dad[j] != 0) && (pp[j] == -1)){
                        dad[j] = rule->getParent(dad[j]);
                        pp[j] = mset.getSlot(dad);
                    }
                    dad[j] = p[j];
                }
            }
        }
    }
}

std::vector<int> MultiIndexManipulations::computeLevels(MultiIndexSet const &mset, BaseRuleLocalPolynomial const *rule){
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    int num_points = mset.getNumIndexes();
    std::vector<int> level((size_t) num_points);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = mset.getIndex(i);
        int current_level = rule->getLevel(p[0]);
        for(size_t j=1; j<num_dimensions; j++){
            current_level += rule->getLevel(p[j]);
        }
        level[i] = current_level;
    }
    return level;
}

void MultiIndexManipulations::selectFlaggedChildren(const MultiIndexSet &mset, const std::vector<bool> &flagged, const std::vector<int> &level_limits, MultiIndexSet &new_set){
    new_set = MultiIndexSet(mset.getNumDimensions());
    size_t num_dimensions = (size_t) mset.getNumDimensions();

    Data2D<int> children_unsorted;
    children_unsorted.resize(mset.getNumDimensions(), 0);

    std::vector<int> kid(num_dimensions);

    int n = mset.getNumIndexes();
    if (level_limits.empty()){
        for(int i=0; i<n; i++){
            if (flagged[i]){
                std::copy_n(mset.getIndex(i), num_dimensions, kid.data());
                for(auto &k : kid){
                    k++;
                    if (mset.missing(kid)) children_unsorted.appendStrip(kid);
                    k--;
                }
            }
        }
    }else{
        for(int i=0; i<n; i++){
            if (flagged[i]){
                std::copy_n(mset.getIndex(i), num_dimensions, kid.data());
                auto ill = level_limits.begin();
                for(auto &k : kid){
                    k++;
                    if (((*ill == -1) || (k <= *ill)) && mset.missing(kid))
                        children_unsorted.appendStrip(kid);
                    k--;
                    ill++;
                }
            }
        }
    }
    if (children_unsorted.getNumStrips() > 0)
        new_set.addData2D(children_unsorted);
}

void MultiIndexManipulations::removeIndexesByLimit(const std::vector<int> &level_limits, MultiIndexSet &mset){
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    int num_indexes = mset.getNumIndexes();
    std::vector<int> keep;
    keep.reserve((size_t) num_indexes);
    auto imset = mset.getVector().begin();
    for(int i=0; i<num_indexes; i++){
        bool obey = true;
        for(auto l : level_limits){
            if ((l > -1) && (*imset > l)) obey = false;
            imset++;
        }
        if (obey) keep.push_back(i);
    }
    if (keep.size() == 0){
        mset = MultiIndexSet((int) num_dimensions);
    }else if (keep.size() < (size_t) num_indexes){
        std::vector<int> new_indexes(num_dimensions * keep.size());
        auto inew = new_indexes.begin();
        for(auto i : keep){
            std::copy_n(mset.getIndex(i), num_dimensions, inew);
            std::advance(inew, num_dimensions);
        }
        mset.setIndexes(new_indexes);
    }
}

void MultiIndexManipulations::generateNestedPoints(const MultiIndexSet &tensors, std::function<int(int)> getNumPoints, MultiIndexSet &points){
    size_t num_dimensions = (size_t) tensors.getNumDimensions();
    Data2D<int> raw_points;
    raw_points.resize((int) num_dimensions, 0);

    std::vector<int> num_points_delta(num_dimensions);
    std::vector<int> offsets(num_dimensions);
    std::vector<int> index(num_dimensions);

    for(int i=0; i<tensors.getNumIndexes(); i++){
        const int *p = tensors.getIndex(i);
        size_t num_total = 1;
        for(size_t j=0; j<num_dimensions; j++){
            num_points_delta[j] = getNumPoints(p[j]);
            if (p[j] > 0){
                offsets[j] = getNumPoints(p[j]-1);
                num_points_delta[j] -= offsets[j];
            }else{
                offsets[j] = 0;
            }

            num_total *= (size_t) num_points_delta[j];
        }

        for(size_t k=0; k<num_total; k++){
            size_t t = k;
            for(int j = (int) num_dimensions-1; j>=0; j--){
                index[j] = offsets[j] + (int) (t % num_points_delta[j]);
                t /= (size_t) num_points_delta[j];
            }
            raw_points.appendStrip(index);
        }
    }

    points = MultiIndexSet();
    points.setNumDimensions((int) num_dimensions);
    points.addData2D(raw_points);
}

void MultiIndexManipulations::generateNonNestedPoints(const MultiIndexSet &tensors, const OneDimensionalWrapper &wrapper, MultiIndexSet &points){
    size_t num_dimensions = (size_t) tensors.getNumDimensions();
    int num_tensors = tensors.getNumIndexes();
    std::vector<MultiIndexSet> point_tensors((size_t) num_tensors);

    #pragma omp parallel for
    for(int i=0; i<num_tensors; i++){
        std::vector<int> npoints(num_dimensions);
        const int *p = tensors.getIndex(i);
        for(size_t j=0; j<num_dimensions; j++)
            npoints[j] = wrapper.getNumPoints(p[j]);
        MultiIndexSet raw_points;
        MultiIndexManipulations::generateFullTensorSet<int>(npoints, raw_points);

        size_t j = 0;
        for(auto &g : raw_points.getVector()){
            g = wrapper.getPointIndex(p[j++], g); // remap local-order-to-global-index
            if (j == num_dimensions) j = 0;
        }
        point_tensors[i].setNumDimensions((int) num_dimensions);
        point_tensors[i].addUnsortedInsexes(raw_points.getVector());
    }

    MultiIndexManipulations::unionSets<true>(point_tensors, points);
}

void MultiIndexManipulations::computeTensorWeights(const MultiIndexSet &mset, std::vector<int> &weights){
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    int num_tensors = mset.getNumIndexes();

    std::vector<int> level;
    computeLevels(mset, level);
    int max_level = *std::max_element(level.begin(), level.end());

    Data2D<int> dag_down;
    dag_down.resize((int) num_dimensions, num_tensors);

    weights.resize(num_tensors);

    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_tensors; i++){
        std::vector<int> kid(num_dimensions);
        std::copy_n(mset.getIndex(i), num_dimensions, kid.data());

        int *ref_kids = dag_down.getStrip(i);
        for(size_t j=0; j<num_dimensions; j++){
            kid[j]++;
            ref_kids[j] = mset.getSlot(kid);
            kid[j]--;
        }

        if (level[i] == max_level) weights[i] = 1;
    }

    for(int l=max_level-1; l>=0; l--){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<num_tensors; i++){
            if (level[i] == l){
                std::vector<int> monkey_tail(max_level-l+1);
                std::vector<int> monkey_count(max_level-l+1);
                std::vector<bool> used(num_tensors, false);

                int current = 0;
                monkey_count[0] = 0;
                monkey_tail[0] = i;

                int sum = 0;

                while(monkey_count[0] < (int) num_dimensions){
                    if (monkey_count[current] < (int) num_dimensions){
                        int branch = dag_down.getStrip(monkey_tail[current])[monkey_count[current]];
                        if ((branch == -1) || (used[branch])){
                            monkey_count[current]++;
                        }else{
                            used[branch] = true;
                            sum += weights[branch];
                            monkey_count[++current] = 0;
                            monkey_tail[current] = branch;
                        }
                    }else{
                        monkey_count[--current]++;
                    }
                }

                weights[i] = 1 - sum;
            }
        }
    }
}

void MultiIndexManipulations::createActiveTensors(const MultiIndexSet &mset, const std::vector<int> &weights, MultiIndexSet &active){
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    size_t nz_weights = 0;
    for(auto w: weights) if (w != 0) nz_weights++;

    std::vector<int> indexes(nz_weights * num_dimensions);
    nz_weights = 0;
    auto iter = indexes.begin();
    auto iset = mset.getVector().begin();
    for(auto w: weights){
        if (w != 0){
            std::copy_n(iset, num_dimensions, iter);
            std::advance(iter, num_dimensions);
        }
        std::advance(iset, num_dimensions);
    }

    active.setNumDimensions((int) num_dimensions);
    active.setIndexes(indexes);
}

void MultiIndexManipulations::createPolynomialSpace(const MultiIndexSet &tensors, std::function<int(int)> exactness, MultiIndexSet &space){
    size_t num_dimensions = (size_t) tensors.getNumDimensions();
    int num_tensors = tensors.getNumIndexes();
    std::vector<MultiIndexSet> polynomial_tensors((size_t) num_tensors);

    #pragma omp parallel for
    for(int i=0; i<num_tensors; i++){
        std::vector<int> npoints(num_dimensions);
        const int *p = tensors.getIndex(i);
        for(size_t j=0; j<num_dimensions; j++)
            npoints[j] = exactness(p[j]) + 1;
        MultiIndexManipulations::generateFullTensorSet<int>(npoints, polynomial_tensors[i]);
    }

    MultiIndexManipulations::unionSets<true>(polynomial_tensors, space);
}

}

#endif
