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

namespace MultiIndexManipulations{

/*!
 * \internal
 * \ingroup TasmanianMultiIndexManipulations
 * \brief Generate a series of \b level_sets where each set has the parents/children of the previous one that satisfy the \b inside() criteria.
 *
 * On entry, \b level_sets must constain at least one set.
 * The function takes the last set in \b level_sets and adds a new set of
 * eithe the parents or children that also satisfy the \b inside() condition.
 * The process is repeated until the new set is empty.
 * \endinternal
 */
template<bool use_parents_direction>
void repeatAddIndexes(std::function<bool(const std::vector<int> &index)> inside, std::vector<MultiIndexSet> &level_sets){
    size_t num_dimensions = level_sets.back().getNumDimensions();
    bool adding = true;
    while(adding){
        Data2D<int> level((int) num_dimensions, 0);
        int num_indexes = level_sets.back().getNumIndexes();
        #pragma omp parallel for
        for(int i=0; i<num_indexes; i++){
            std::vector<int> point(num_dimensions);
            std::copy_n(level_sets.back().getIndex(i), num_dimensions, point.data());
            for(auto &p : point){
                p += (use_parents_direction) ? -1 : 1; // parents have lower index, children have higher indexes
                if ( (not use_parents_direction or p >= 0) and inside(point) ){
                    #pragma omp critical
                    {
                        level.appendStrip(point);
                    }
                }
                p -= (use_parents_direction) ? -1 : 1; // restore p
            }
        }

        adding = (level.getNumStrips() > 0);
        if (adding)
            level_sets.push_back(MultiIndexSet(level));
    }
}

/*!
 * \internal
 * \ingroup TasmanianMultiIndexManipulations
 * \brief Retuns the union of all \b level_sets, all sets are destroyed in the process.
 *
 * \endinternal
 */
inline MultiIndexSet unionSets(std::vector<MultiIndexSet> &level_sets){
    long long num_levels = level_sets.size();
    while(num_levels > 1){
        long long stride = num_levels / 2 + (((num_levels % 2) > 0) ? 1 : 0);
        #pragma omp parallel for
        for(long long i=0; i<stride; i++)
            if (i + stride < num_levels)
                level_sets[i] += level_sets[i + stride];
        num_levels = stride;
    }
    return std::move(level_sets[0]);
}

void completeSetToLower(MultiIndexSet &set){
    size_t num_dimensions = set.getNumDimensions();
    int num = set.getNumIndexes();
    Data2D<int> completion((int) num_dimensions, 0);
    #pragma omp parallel for
    for(int i=0; i<num; i++){
        std::vector<int> point(num_dimensions);
        std::copy_n(set.getIndex(i), num_dimensions, point.data());
        for(auto &p : point){
            if (p != 0){
                p--;
                if (set.missing(point)){
                    #pragma omp critical
                    {
                        completion.appendStrip(point);
                    }
                }
                p++;
            }
        }
    }

    if (completion.getNumStrips() > 0){
        std::vector<MultiIndexSet> level_sets = { MultiIndexSet(completion) };

        constexpr bool use_parents = true;
        repeatAddIndexes<use_parents>([&](std::vector<int> const &p) -> bool{ return set.missing(p); }, level_sets);

        set += unionSets(level_sets);
    }
}

/*!
 * \internal
 * \ingroup TasmanianMultiIndexManipulations
 * \brief Generate the minimum lower complete multi-index set that includes the indexes satisfying \b criteria(), assumes \b criteria() defines a connected set.
 *
 * \endinternal
 */
inline MultiIndexSet generateGeneralMultiIndexSet(size_t num_dimensions, std::function<bool(const std::vector<int> &index)> criteria){
    std::vector<MultiIndexSet> level_sets = { MultiIndexSet(num_dimensions, std::vector<int>(num_dimensions, 0)) };

    constexpr bool use_parents = false;
    repeatAddIndexes<use_parents>(criteria, level_sets);

    MultiIndexSet set = unionSets(level_sets);

    completeSetToLower(set);
    return set;
}

/*!
 * \internal
 * \brief Generate the multi-index of indexes with weighs less than the \b normalized_offset.
 *
 * The weight of an index uses the \b weights combined with the \b rule_exactness().
 * Called only when the set is guaranteed to be lower complete,
 * then the one dimensional weights can be cached prior to running the selection algorithm.
 *
 * If \b check_limits is \b false, then \b level_limits are ignored for speedup.
 * \endinternal
 */
template<bool check_limits>
MultiIndexSet selectLowerSet(ProperWeights const &weights, std::function<int(int i)> rule_exactness,
                             int normalized_offset, std::vector<int> const &level_limits){
    size_t num_dimensions = weights.getNumDimensions();
    if (weights.contour == type_level){
        auto cache = generateLevelWeightsCache<int, type_level, false>(weights, rule_exactness, normalized_offset);
        return generateLowerMultiIndexSet(num_dimensions,
                [&](std::vector<int> const &index)->bool{
                    if (check_limits) for(size_t j=0; j<num_dimensions; j++) if ((level_limits[j] > -1) && (index[j] > level_limits[j])) return false;
                    return (getIndexWeight<int, type_level>(index.data(), cache) <= normalized_offset);
                });
    }else if (weights.contour == type_curved){
        auto cache = generateLevelWeightsCache<double, type_curved, false>(weights, rule_exactness, normalized_offset);
        double noff = (double) normalized_offset;
        return generateLowerMultiIndexSet(num_dimensions,
                [&](std::vector<int> const &index)->bool{
                    if (check_limits) for(size_t j=0; j<num_dimensions; j++) if ((level_limits[j] > -1) && (index[j] > level_limits[j])) return false;
                    return (std::ceil(getIndexWeight<double, type_curved>(index.data(), cache)) <= noff);
                });
    }else{ // type_hyperbolic
        auto cache = generateLevelWeightsCache<double, type_hyperbolic, false>(weights, rule_exactness, normalized_offset);
        double noff = (double) normalized_offset;
        return generateLowerMultiIndexSet(num_dimensions,
                [&](std::vector<int> const &index)->bool{
                    if (check_limits) for(size_t j=0; j<num_dimensions; j++) if ((level_limits[j] > -1) && (index[j] > level_limits[j])) return false;
                    return (std::ceil(getIndexWeight<double, type_hyperbolic>(index.data(), cache)) <= noff);
                });
    }
}

/*!
 * \internal
 * \ingroup TasmanianMultiIndexManipulations
 * \brief Generates the minimum lower complete set that contains all indexes with weights less than \b normalized_offset.
 *
 * The weight of an index uses the \b weights combined with the \b rule_exactness().
 * Called only for contour \b type_curved and caches values on-the-fly.
 *
 * If \b check_limits is \b false, then \b level_limits are ignored for speedup.
 * \endinternal
 */
template<bool check_limits>
MultiIndexSet selectGeneralSet(ProperWeights const &weights, std::function<int(int i)> rule_exactness,
                                      int normalized_offset, std::vector<int> const &level_limits){
    size_t num_dimensions = weights.getNumDimensions();
    std::vector<std::vector<double>> cache(num_dimensions);
    for(size_t j=0; j<num_dimensions; j++) cache[j].push_back(0.0);
    double noff = (double) normalized_offset;
    return generateGeneralMultiIndexSet(num_dimensions,
                                        [&](std::vector<int> const &index) -> bool{
                                            if (check_limits) for(size_t j=0; j<num_dimensions; j++) if (index[j] > level_limits[j]) return false;
                                            double w = 0;
                                            for(size_t j=0; j<num_dimensions; j++){
                                                while(index[j] >= (int) cache[j].size()){
                                                    int exactness = 1 + rule_exactness((int)(cache[j].size() - 1));
                                                    cache[j].push_back( (double) weights.linear[j] * exactness +
                                                                         weights.curved[j] * std::log1p((double) exactness) );
                                                }
                                                w += cache[j][index[j]];
                                            }
                                            return (std::ceil(w) <= noff);
                                        });
}


MultiIndexSet selectTensors(size_t num_dimensions, int offset, TypeDepth type,
                                                     std::function<int(int i)> rule_exactness, std::vector<int> const &anisotropic_weights,
                                                     std::vector<int> const &level_limits){
    // special case of a full tensor selection
    if ((type == type_tensor) || (type == type_iptensor) || (type == type_qptensor)){ // special case, full tensor
        std::vector<int> max_exactness = (anisotropic_weights.empty()) ? std::vector<int>(num_dimensions, 1) : anisotropic_weights;
        for(auto &e : max_exactness) e *= offset;
        std::vector<int> num_points(num_dimensions, 0); // how many points to have in each direction
        std::transform(max_exactness.begin(), max_exactness.end(), num_points.begin(),
                        [&](int e)-> int{
                            int l = 0; // level
                            while(rule_exactness(l) < e) l++; // get the first level that covers the weight
                            return l+1; // adding extra one to change interpretation from 0-index "level" to 1-index "number of points"
                        });
        if (!level_limits.empty()){
            for(size_t j=0; j<num_dimensions; j++)
                if (level_limits[j] >= 0)
                    num_points[j] = std::min(num_points[j], level_limits[j]+1); // the +1 indicates switch from max-level to number-of-points
        }
        return generateFullTensorSet(num_points);
    }

    ProperWeights weights(num_dimensions, type, anisotropic_weights);

    int normalized_offset = offset * weights.minLinear();
    if (weights.provenLower()){ // if the set is guaranteed to be lower
        if (level_limits.empty()){
            return selectLowerSet<false>(weights, rule_exactness, normalized_offset, level_limits);
        }else{
            return selectLowerSet<true>(weights, rule_exactness, normalized_offset, level_limits);
        }
    }else{
        if (level_limits.empty()){
            return selectGeneralSet<false>(weights, rule_exactness, normalized_offset, level_limits);
        }else{
            return selectGeneralSet<true>(weights, rule_exactness, normalized_offset, level_limits);
        }
    }
}

std::vector<int> computeLevels(MultiIndexSet const &mset){
    // cannot add as inline to the public header due to the pragma and possible "unknown pragma" message
    int num_indexes = mset.getNumIndexes();
    size_t num_dimensions = mset.getNumDimensions();
    std::vector<int> levels((size_t) num_indexes);
    #pragma omp parallel for
    for(int i=0; i<num_indexes; i++){
        const int* p = mset.getIndex(i);
        levels[i] = std::accumulate(p, p + num_dimensions, 0);
    }
    return levels;
}

std::vector<int> getMaxIndexes(const MultiIndexSet &mset){
    size_t num_dimensions = mset.getNumDimensions();
    std::vector<int> max_index(num_dimensions, 0);
    int n = mset.getNumIndexes();

    #pragma omp parallel
    {
        std::vector<int> local_max_index(num_dimensions, 0);
        #pragma omp for
        for(int i=0; i<n; i++){
            const int* p = mset.getIndex(i);
            for(size_t j=0; j<num_dimensions; j++) if (local_max_index[j] < p[j]) local_max_index[j] = p[j];
        }
        #pragma omp critical
        {
            for(size_t j=0; j<num_dimensions; j++){
                if (max_index[j] < local_max_index[j]) max_index[j] = local_max_index[j];
            }
        }
    }

    return max_index;
}

Data2D<int> computeDAGup(MultiIndexSet const &mset){
    size_t num_dimensions = (size_t) mset.getNumDimensions();
    int n = mset.getNumIndexes();
    Data2D<int> parents(mset.getNumDimensions(), n);
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
    return parents;
}

MultiIndexSet selectFlaggedChildren(const MultiIndexSet &mset, const std::vector<bool> &flagged, const std::vector<int> &level_limits){
    size_t num_dimensions = mset.getNumDimensions();
    int n = mset.getNumIndexes();

    Data2D<int> children_unsorted(num_dimensions, 0);

    #ifdef _OPENMP
    #pragma omp parallel
    {
        Data2D<int> lrefined(num_dimensions, 0);
        std::vector<int> kid(num_dimensions);

        if (level_limits.empty()){
            #pragma omp for
            for(int i=0; i<n; i++){
                if (flagged[i]){
                    std::copy_n(mset.getIndex(i), num_dimensions, kid.data());
                    for(auto &k : kid){
                        k++;
                        if (mset.missing(kid)) lrefined.appendStrip(kid);
                        k--;
                    }
                }
            }
        }else{
            #pragma omp for
            for(int i=0; i<n; i++){
                if (flagged[i]){
                    std::copy_n(mset.getIndex(i), num_dimensions, kid.data());
                    auto ill = level_limits.begin();
                    for(auto &k : kid){
                        k++;
                        if (((*ill == -1) || (k <= *ill)) && mset.missing(kid))
                            lrefined.appendStrip(kid);
                        k--;
                        ill++;
                    }
                }
            }
        }

        #pragma omp critical
        {
            children_unsorted.append(lrefined);
        }
    }
    #else
    std::vector<int> kid(num_dimensions);

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
    #endif

    return MultiIndexSet(children_unsorted);
}

MultiIndexSet generateNestedPoints(const MultiIndexSet &tensors, std::function<int(int)> getNumPoints){
    size_t num_dimensions = (size_t) tensors.getNumDimensions();
    std::vector<MultiIndexSet> delta_sets((size_t) tensors.getNumIndexes());

    #pragma omp parallel for
    for(int i=0; i<tensors.getNumIndexes(); i++){
        Data2D<int> raw_points(num_dimensions, 0);

        std::vector<int> num_points_delta(num_dimensions);
        std::vector<int> offsets(num_dimensions);
        std::vector<int> index(num_dimensions);

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

        delta_sets[i] = MultiIndexSet(raw_points);
    }

    return unionSets(delta_sets);
}

MultiIndexSet generateNonNestedPoints(const MultiIndexSet &tensors, const OneDimensionalWrapper &wrapper){
    size_t num_dimensions = tensors.getNumDimensions();
    int num_tensors = tensors.getNumIndexes();
    std::vector<MultiIndexSet> point_tensors((size_t) num_tensors);

    #pragma omp parallel for
    for(int t=0; t<num_tensors; t++){
        std::vector<int> num_entries(num_dimensions);
        const int *p = tensors.getIndex(t);
        std::transform(p, p + num_dimensions, num_entries.begin(), [&](int l)->int{ return wrapper.getNumPoints(l); });

        int num_total = 1;
        for(auto &l : num_entries) num_total *= l;

        Data2D<int> raw_points(num_dimensions, num_total);
        auto iter = raw_points.rbegin();
        for(int i=num_total-1; i>=0; i--){
            int d = i;
            auto l = num_entries.rbegin();
            // in order to generate indexes in the correct order, the for loop must go backwards
            for(size_t j = 0; j<num_dimensions; j++){
                *iter++ = wrapper.getPointIndex(p[num_dimensions - j - 1], (d % *l));
                d /= *l++;
            }
        }

        point_tensors[t] = MultiIndexSet(raw_points);
    }

    return unionSets(point_tensors);
}

void resortIndexes(const MultiIndexSet &iset, std::vector<std::vector<int>> &map, std::vector<std::vector<int>> &lines1d) {
    int num_dimensions = iset.getNumDimensions();
    int num_tensors = iset.getNumIndexes();

    int num_levels = 1 + *std::max_element(iset.begin(), iset.end());

    auto match_outside_dim = [&](int d, int const*a, int const *b) -> bool {
        for(int j=0; j<d; j++)
            if (a[j] != b[j]) return false;
        for(int j=d+1; j<num_dimensions; j++)
            if (a[j] != b[j]) return false;
        return true;
    };

    map = std::vector<std::vector<int>>(num_dimensions, std::vector<int>(num_tensors));
    lines1d = std::vector<std::vector<int>>(num_dimensions);
    for(auto &ji : lines1d) ji.reserve(num_levels);

    #pragma omp parallel for
    for(int d=0; d<num_dimensions; d++) {
        // for each dimension, use the map to group indexes together
        std::iota(map[d].begin(), map[d].end(), 0);
        if (d != num_dimensions - 1) {
            std::sort(map[d].begin(), map[d].end(), [&](int a, int b)->bool{
                const int * idxa = iset.getIndex(a);
                const int * idxb = iset.getIndex(b);
                for(int j=0; j<num_dimensions; j++) {
                    if (j != d){
                        if (idxa[j] < idxb[j]) return true;
                        if (idxa[j] > idxb[j]) return false;
                    }
                }
                // lexigographical order, dimension d is the fastest moving one
                if (idxa[d] < idxb[d]) return true;
                if (idxa[d] > idxb[d]) return false;
                return false;
            });
        }

        if (num_dimensions == 1) {
            lines1d[d].push_back(0);
            lines1d[d].push_back(num_tensors);
        } else {
            int const *c_index = iset.getIndex(map[d][0]);
            lines1d[d].push_back(0);
            for(int i=1; i<num_tensors; i++) {
                if (not match_outside_dim(d, c_index, iset.getIndex(map[d][i]))) {
                    lines1d[d].push_back(i);
                    c_index = iset.getIndex(map[d][i]);
                }
            }
            lines1d[d].push_back(num_tensors);
        }
    }
}

std::vector<int> computeTensorWeights(MultiIndexSet const &mset){
    int num_dimensions = mset.getNumDimensions();
    int num_tensors = mset.getNumIndexes();

    if (num_dimensions == 1) {
        std::vector<int> weights(num_tensors, 0);
        weights.back() = 1;
        return weights;
    }

    std::vector<std::vector<int>> map;
    std::vector<std::vector<int>> lines1d;

    resortIndexes(mset, map, lines1d);

    Data2D<int> dag_down(num_dimensions, num_tensors);

    std::vector<int> weights(num_tensors, 0);

    // the row with contiguous indexes has a trivial solution
    auto const& last_jobs = lines1d[num_dimensions-1];
    for(int i=0; i<static_cast<int>(last_jobs.size() - 1); i++)
        weights[last_jobs[i+1] - 1] = 1;

    for(int d=num_dimensions-2; d>=0; d--) {
        #pragma omp parallel for
        for(int job = 0; job < static_cast<int>(lines1d[d].size() - 1); job++) {
            for(int i=lines1d[d][job+1]-2; i>=lines1d[d][job]; i--) {
                int &val = weights[map[d][i]];
                for(int j=i+1; j<lines1d[d][job+1]; j++) {
                    val -= weights[map[d][j]];
                }
            }
        }
    }

    return weights;
}

MultiIndexSet createPolynomialSpace(const MultiIndexSet &tensors, std::function<int(int)> exactness){
    size_t num_dimensions = (size_t) tensors.getNumDimensions();
    int num_tensors = tensors.getNumIndexes();
    std::vector<MultiIndexSet> polynomial_tensors((size_t) num_tensors);

    #pragma omp parallel for
    for(int i=0; i<num_tensors; i++){
        std::vector<int> npoints(num_dimensions);
        const int *p = tensors.getIndex(i);
        for(size_t j=0; j<num_dimensions; j++)
            npoints[j] = exactness(p[j]) + 1;
        polynomial_tensors[i] = generateFullTensorSet(npoints);
    }

    return unionSets(polynomial_tensors);
}

std::vector<int> inferAnisotropicWeights(AccelerationContext const *acceleration, TypeOneDRule rule, TypeDepth depth,
                                         MultiIndexSet const &points, std::vector<double> const &coefficients, double tol){

    int num_dimensions = static_cast<int>(points.getNumDimensions());
    int cols = (OneDimensionalMeta::getControurType(depth) == type_curved) ?
                2 * num_dimensions + 1 : num_dimensions + 1;

    int data_rows = static_cast<int>(std::count_if(coefficients.begin(), coefficients.end(), [=](double c)->bool{ return (std::abs(c) > tol); }));

    Data2D<double> A(data_rows + cols, cols, 0.0);
    std::vector<double> b(data_rows + cols, 0.0);

    int c = 0;
    for(int i=0; i<points.getNumIndexes(); i++){
        if (std::abs(coefficients[i]) > tol){
            int const *indx = points.getIndex(i);
            for(int j=0; j<num_dimensions; j++){
                A.getStrip(j)[c] = static_cast<double>(indx[j]);
            }
            A.getStrip(cols-1)[c] = 1.0;
            b[c++] = -std::log(std::abs(coefficients[i]));
        }
    }

    if (rule == rule_fourier){
        for(int j=0; j<num_dimensions; j++){
            double *cc = A.getStrip(j);
            for(int i=0; i<data_rows; i++)
                cc[i] = static_cast<double>((static_cast<int>(cc[i]) + 1) / 2);
        }
    }

    if (OneDimensionalMeta::getControurType(depth) == type_hyperbolic){
        for(int j=0; j<num_dimensions; j++){
            double *cc = A.getStrip(j);
            for(int i=0; i<data_rows; i++)
                cc[i] = std::log(cc[i] + 1.0);
        }
    }

    if (OneDimensionalMeta::getControurType(depth) == type_curved){
        for(int j=0; j<num_dimensions; j++){
            double *cc = A.getStrip(j);
            double *vv = A.getStrip(num_dimensions + j);
            for(int i=0; i<data_rows; i++)
                vv[i] = std::log(cc[i] + 1.0);
        }
    }

    // pad with a block of identity to provide regularization
    // adds 0.003^2 to the diagonal of the A^T * A matrix
    for(int j=0; j<cols; j++)
        A.getStrip(j)[data_rows + j] = 0.003;

    std::vector<double> x(cols);
    TasmanianDenseSolver::solveLeastSquares(acceleration, data_rows + cols, cols, A.getStrip(0), b.data(), x.data());

    std::vector<int> weights(cols - 1);
    for(size_t j=0; j<weights.size(); j++)
        weights[j] = static_cast<int>(x[j] * 1000.0 + 0.5);

    int max_weight = *std::max_element(weights.begin(), weights.begin() + num_dimensions);

    if (max_weight <= 0){ // all directions are diverging, default to isotropic total degree
        std::fill_n(weights.begin(), num_dimensions, 1);
        if (OneDimensionalMeta::getControurType(depth) == type_curved)
            std::fill(weights.begin() + num_dimensions, weights.end(), 0);
    }else{
        int min_weight = max_weight;
        for(int j=0; j<num_dimensions; j++)
            if (weights[j] > 0 and weights[j] < min_weight) min_weight = weights[j];  // find the smallest positive weight
        for(int j=0; j<num_dimensions; j++){
            if (weights[j] <= 0){
                weights[j] = min_weight;
                if (OneDimensionalMeta::getControurType(depth) == type_curved){
                    if (std::abs(weights[num_dimensions + j]) > weights[j])
                        weights[num_dimensions + j] = (weights[j] > 0.0) ? weights[num_dimensions + j] : -weights[num_dimensions + j];
                }
            }
        }
    }

    return weights;
}

} // MultiIndexManipulations

} // TasGrid

#endif
