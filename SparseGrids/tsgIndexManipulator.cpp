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

IndexManipulator::IndexManipulator(int cnum_dimensions, const CustomTabulated* custom) : num_dimensions(cnum_dimensions), meta(custom){}
IndexManipulator::~IndexManipulator(){}

long long IndexManipulator::getIndexWeight(const std::vector<int> &index, TypeDepth type, const std::vector<int> &weights, TypeOneDRule rule) const{
    switch (type){
        case type_level:         return getIndexWeight<type_level>(index, weights, rule);
        case type_curved:        return getIndexWeight<type_curved>(index, weights, rule);
        case type_iptotal:       return getIndexWeight<type_iptotal>(index, weights, rule);;
        case type_ipcurved:      return getIndexWeight<type_ipcurved>(index, weights, rule);
        case type_qptotal:       return getIndexWeight<type_qptotal>(index, weights, rule);
        case type_qpcurved:      return getIndexWeight<type_qpcurved>(index, weights, rule);
        case type_hyperbolic:    return getIndexWeight<type_hyperbolic>(index, weights, rule);
        case type_iphyperbolic:  return getIndexWeight<type_iphyperbolic>(index, weights, rule);
        case type_qphyperbolic:  return getIndexWeight<type_qphyperbolic>(index, weights, rule);
        default:
            return 0;
    }
}

IndexSet* IndexManipulator::selectTensors(int offset, TypeDepth type, const std::vector<int> &anisotropic_weights, TypeOneDRule rule) const{
    // construct the minimum tensor set that covers the target_space defined by offset, type, and anisotropic weights
    // consult the manual for the detailed definition of each case
    // used for constructing Global and Sequence grids
    // implements Corollary 1 of Theorem 1 from
    // Stoyanov, Webster: "A dynamically adaptive sparse grids method for quasi-optimal interpolation of multidimensional functions"
    // Computers & Mathematics with Applications, 71(11):2449–2465, 2016

    // This cheats a bit, but it handles the special case when we want a full tensor grid
    // instead of computing the grid in the standard gradual level by level way,
    // the tensor grid ca be computed directly with two for-loops
    if ((type == type_tensor) || (type == type_iptensor) || (type == type_qptensor)){ // special case, full tensor
        std::vector<int> levels; // how many levels to keep in each direction
        if (type == type_tensor){
            levels.resize(num_dimensions, offset);
            if (!anisotropic_weights.empty()) for(int j=0; j<num_dimensions; j++) levels[j] *= anisotropic_weights[j];
        }else if (type == type_iptensor){
            levels.resize(num_dimensions, 0);
            for(int j=0; j<num_dimensions; j++){
                long long target = offset;
                if (!anisotropic_weights.empty()) target *= anisotropic_weights[j];
                while(meta.getIExact(levels[j], rule) < target) levels[j]++;
            }
        }else{
            levels.resize(num_dimensions, 0);
            for(int j=0; j<num_dimensions; j++){
                long long target = offset;
                if (!anisotropic_weights.empty()) target *= anisotropic_weights[j];
                while(meta.getQExact(levels[j], rule) < target) levels[j]++;
            }
        }
        long long num_total = 1;
        for(auto &l : levels) num_total *= ++l;
        std::vector<int> index(num_total * num_dimensions);
        auto iter = index.rbegin();
        for(long long i=num_total-1; i>=0; i--){
            long long t = i;
            auto l = levels.rbegin();
            // in order to generate indexes in the correct order, the for loop must go backwards
            for(int j = num_dimensions-1; j>=0; j--){
                *iter++ = (int) (t % *l);
                t /= *l++;
            }
        }
        return new IndexSet(num_dimensions, index);
    }

    // non-tensor cases
    std::vector<int> weights;
    getProperWeights(type, anisotropic_weights, weights);

    // compute normalization and check if heavily curved
    long long normalized_offset = weights[0];
    for(int i=1; i<num_dimensions; i++){
        if (normalized_offset > weights[i]) normalized_offset = weights[i];
    }
    bool known_lower = true;
    if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
        for(int i=0; i<num_dimensions; i++) if (weights[i] + weights[i+num_dimensions] < 0) known_lower = false;

    }
    normalized_offset *= offset;
    IndexSet *total = 0;

    if (known_lower){ // use fast algorithm, but only works for sets guaranteed to be lower
        int c = num_dimensions -1;
        bool outside = false;
        std::vector<int> root(num_dimensions, 0);
        std::vector<int> index_dump;
        while( !(outside && (c == 0)) ){
            if (outside){
                for(int k=c; k<num_dimensions; k++) root[k] = 0;
                c--;
                root[c]++;
            }else{
                for(auto i : root) index_dump.push_back(i);
                c = num_dimensions-1;
                root[c]++;
            }
            outside = (getIndexWeight(root, type, weights, rule) > normalized_offset);
        }
        total = new IndexSet(num_dimensions, index_dump);
    }else{ // use slower algorithm, but more general
        GranulatedIndexSet **sets;
        std::vector<int> root(num_dimensions, 0);
        GranulatedIndexSet *set_level = new GranulatedIndexSet(num_dimensions, 1);  set_level->addIndex(root.data());
        total = new IndexSet(num_dimensions);
        bool adding = true;
        while(adding){
            int num_sets = (set_level->getNumIndexes() / 8) + 1;

            sets = new GranulatedIndexSet*[num_sets];

            for(int me=0; me<num_sets; me++){

                sets[me] = new GranulatedIndexSet(num_dimensions);
                std::vector<int> newp(num_dimensions);
                for(int i=me; i<set_level->getNumIndexes(); i+=num_sets){
                    const int *p = set_level->getIndex(i);
                    std::copy(p, p + num_dimensions, newp.data());
                    for(auto &j : newp){
                        j++;
                        if ((getIndexWeight(newp, type, weights, rule) <= normalized_offset) && (set_level->getSlot(newp.data()) == -1) && (total->getSlot(newp.data()) == -1)){
                            sets[me]->addIndex(newp.data());
                        }
                        j--;
                    }
                }
            }

            int warp = num_sets;
            while(warp > 1){
                #pragma omp parallel for schedule(dynamic)
                for(int i=0; i<warp-1; i+=2){
                    sets[i]->addGranulatedSet(sets[i+1]);
                    delete sets[i+1];
                    sets[i+1] = 0;
                }

                for(int i=1; i<warp/2; i++){
                    sets[i] = sets[2*i];
                    sets[2*i] = 0;
                }
                if (warp % 2 == 1){ sets[warp/2] = sets[warp-1]; sets[warp-1] = 0; }
                warp = warp / 2 + warp % 2;
            }

            adding = (sets[0]->getNumIndexes() > 0);
            total->addGranulatedSet(set_level);
            delete set_level;
            set_level = sets[0];
            sets[0] = 0;
            delete[] sets;
        }

        if (total == 0){
            total = new IndexSet(set_level);
        }else{
            total->addGranulatedSet(set_level);
        }
        delete set_level;

        // needed to preserved lower property
        IndexSet* completion = getLowerCompletion(total);
        if ((completion != 0) && (completion->getNumIndexes() > 0)){
            total->addIndexSet(completion);
            delete completion;
        }
    }

    return total;
}

IndexSet* IndexManipulator::selectTensors(const IndexSet *target_space, bool integration, TypeOneDRule rule) const{
    // construct the minimum tensor set that covers the target_space
    // used when projecting the Global grid onto Legendre polynomial basis
    // implements the clean version of Theorem 1 from
    // Stoyanov, Webster: "A dynamically adaptive sparse grids method for quasi-optimal interpolation of multidimensional functions"
    // Computers & Mathematics with Applications, 71(11):2449–2465, 2016
    // the other selectTensors() function cover the special case of Corollary 1
    GranulatedIndexSet **sets;
    std::vector<int> root(num_dimensions, 0);
    GranulatedIndexSet *set_level = new GranulatedIndexSet(num_dimensions, 1);  set_level->addIndex(root.data());
    IndexSet *total = new IndexSet(num_dimensions);
    bool adding = true;

    while(adding){
        int num_sets = (set_level->getNumIndexes() / 8) + 1;

        sets = new GranulatedIndexSet*[num_sets];

        #pragma omp parallel for schedule(dynamic)
        for(int me=0; me<num_sets; me++){

            sets[me] = new GranulatedIndexSet(num_dimensions);

            std::vector<int> newp(num_dimensions);
            std::vector<int> corner(num_dimensions);
            if (integration){
                for(int i=me; i<set_level->getNumIndexes(); i+=num_sets){
                    const int *p = set_level->getIndex(i);
                    std::copy(p, p + num_dimensions, newp.data());
                    auto iterp = newp.begin();
                    for(auto &c : corner){
                        c = (*iterp > 0) ? (meta.getQExact(*iterp -1, rule) + 1) : 0;
                        iterp++;
                    }
                    iterp = newp.begin();
                    int oldc;

                    for(auto &c : corner){
                        (*iterp)++; // increment the index of newp
                        oldc = c; // save the corner entry
                        c = meta.getQExact(*iterp -1, rule) + 1; // recompute the corner entry (the index here cannot be zero since we just incremented)

                        if ((target_space->getSlot(corner.data()) != -1) && (set_level->getSlot(newp.data()) == -1) && (total->getSlot(newp.data()) == -1)){
                            sets[me]->addIndex(newp); // add new index
                        }
                        c = oldc; // restore the corner entry
                        (*iterp)--; // restore the index entry
                        iterp++; // move to the next index
                    }
                }
            }else{
                for(int i=me; i<set_level->getNumIndexes(); i+=num_sets){
                    const int *p = set_level->getIndex(i);
                    std::copy(p, p + num_dimensions, newp.data());
                    auto iterp = newp.begin();
                    for(auto &c : corner){
                        c = (*iterp > 0) ? (meta.getIExact(*iterp -1, rule) + 1) : 0;
                        iterp++;
                    }
                    iterp = newp.begin();
                    int oldc;
                    for(auto &c : corner){
                        (*iterp)++;
                        oldc = c;
                        c = meta.getIExact(*iterp -1, rule) + 1;

                        if ((target_space->getSlot(corner.data()) != -1) && (set_level->getSlot(newp.data()) == -1) && (total->getSlot(newp.data()) == -1)){
                            sets[me]->addIndex(newp);
                        }
                        c = oldc;
                        (*iterp)--;
                        iterp++;
                    }
                }
            }
        }

        int warp = num_sets;
        while(warp > 1){
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<warp-1; i+=2){
                sets[i]->addGranulatedSet(sets[i+1]);
                delete sets[i+1];
            }

            for(int i=1; i<warp/2; i++){
                sets[i] = sets[2*i];
            }
            if (warp % 2 == 1){ sets[warp/2] = sets[warp-1];  };
            warp = warp / 2 + warp % 2;
        }

        adding = (sets[0]->getNumIndexes() > 0);
        total->addGranulatedSet(set_level);
        delete set_level;
        set_level = sets[0];
        delete[] sets;
    }

    if (total == 0){
        total = new IndexSet(set_level);
    }else{
        total->addGranulatedSet(set_level);
    }
    delete set_level;

    IndexSet* completion = getLowerCompletion(total);
    if ((completion != 0) && (completion->getNumIndexes() > 0)){
        total->addIndexSet(completion);
        delete completion;
    }

    return total;
}

IndexSet* IndexManipulator::getLowerCompletion(const IndexSet *iset) const{

    GranulatedIndexSet *set_level = new GranulatedIndexSet(num_dimensions);
    std::vector<int> dad(num_dimensions);
    for(int i=0; i<iset->getNumIndexes(); i++){
        const int* p = iset->getIndex(i);
        std::copy(p, p + num_dimensions, dad.data());
        for(auto &d : dad){
            d--;
            if ((d>-1) && (iset->getSlot(dad.data()) == -1)){
                set_level->addIndex(dad);
            }
            d++;
        }
    }

    if (set_level->getNumIndexes() > 0){
        bool adding = true;
        IndexSet *total = new IndexSet(num_dimensions);
        while(adding){
            total->addGranulatedSet(set_level);
            GranulatedIndexSet *set_future = new GranulatedIndexSet(num_dimensions);
            for(int i=0; i<set_level->getNumIndexes(); i++){
                const int* p = set_level->getIndex(i);
                std::copy(p, p + num_dimensions, dad.data());
                for(auto &d : dad){
                    d--;
                    if ((d>-1) && (iset->getSlot(dad.data()) == -1) && (total->getSlot(dad.data()) == -1)){
                        set_future->addIndex(dad);
                    }
                    d++;
                }
            }
            delete set_level;
            set_level = set_future;
            adding = (set_level->getNumIndexes() > 0);
        }
        total->addGranulatedSet(set_level);
        delete set_level;

        return total;
    }else{ // nothing to add, the set is lower
        delete set_level;
        return 0;
    }
}

void IndexManipulator::getProperWeights(TypeDepth type, const std::vector<int> &anisotropic_weights, std::vector<int> &weights) const{
    if (anisotropic_weights.empty()){
        if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
            weights.resize(2*num_dimensions, 1);
            std::fill(&(weights[num_dimensions]), &(weights[num_dimensions]) + num_dimensions, 0);
        }else if ((type == type_hyperbolic) || (type == type_iphyperbolic) || (type == type_qphyperbolic)){
            weights.resize(num_dimensions + 1, 1);
        }else{
            weights.resize(num_dimensions, 1);
        }
    }else{
        if ((type == type_hyperbolic) || (type == type_iphyperbolic) || (type == type_qphyperbolic)){
            weights.resize(num_dimensions + 1);
            std::copy(anisotropic_weights.begin(), anisotropic_weights.end(), weights.data());
            int wsum = 0;
            for(auto w : anisotropic_weights) wsum += w;
            weights[num_dimensions] = wsum;
        }else{
            weights = anisotropic_weights; // copy assign
        }
    }
}

void IndexManipulator::computeLevels(const IndexSet* iset, std::vector<int> &level) const{
    int num_indexes = iset->getNumIndexes();
    level.resize(num_indexes);
    #pragma omp parallel for
    for(int i=0; i<num_indexes; i++){
        const int* p = iset->getIndex(i);
        level[i] = p[0];
        for(int j=1; j<num_dimensions; j++){
            level[i] += p[j];
        }
    }
}

void IndexManipulator::makeTensorWeights(const IndexSet* iset, std::vector<int> &weights) const{
    int n = iset->getNumIndexes();
    std::vector<int> level;
    computeLevels(iset, level);
    weights.resize(n);

    int max_level = 0;
    for(auto l : level) if (l > max_level) max_level = l;

    std::vector<int> kids(((size_t) n) * ((size_t) num_dimensions));
    #pragma omp parallel for schedule(static)
    for(int i=0; i<n; i++){
        const int *p = iset->getIndex(i);
        std::vector<int> kid(num_dimensions);
        std::copy(p, p + num_dimensions, kid.data());
        int *ref_kids = &(kids[((size_t) i) * ((size_t) num_dimensions)]); // cannot use iterators with openmp

        for(int j=0; j<num_dimensions; j++){
            kid[j]++;
            ref_kids[j] = iset->getSlot(kid);
            kid[j]--;
        }
        if (level[i] == max_level){
            weights[i] = 1;
        }
    }

    // this is a hack to avoid the i*num_dimension + j indexing later
    // since kids are not sorted, only getIndex is a valid call here
    IndexSet ikids(num_dimensions, kids);

    for(int l=max_level-1; l>=0; l--){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<n; i++){
            if (level[i] == l){
                std::vector<int> monkey_tail(max_level-l+1);
                std::vector<int> monkey_count(max_level-l+1);
                std::vector<bool> used(n, false);

                int current = 0;
                monkey_count[0] = 0;
                monkey_tail[0] = i;

                int sum = 0;

                while(monkey_count[0] < num_dimensions){
                    if (monkey_count[current] < num_dimensions){
                        int branch = ikids.getIndex(monkey_tail[current])[monkey_count[current]];
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

IndexSet* IndexManipulator::nonzeroSubset(const IndexSet* iset, const std::vector<int> &weights) const{
    size_t nz_weights = 0;
    for(auto w: weights) if (w != 0) nz_weights++;

    std::vector<int> index(nz_weights * ((size_t) num_dimensions));
    nz_weights = 0;
    auto iter = index.begin();
    for(int i=0; i<iset->getNumIndexes(); i++){
        if (weights[i] != 0){
            const int *p = iset->getIndex(i);
            std::copy(p, p + num_dimensions, iter);
            std::advance(iter, num_dimensions);
        }
    }

    return new IndexSet(num_dimensions, index);
}

UnsortedIndexSet* IndexManipulator::tensorGenericPoints(const int levels[], const OneDimensionalWrapper *rule) const{
    std::vector<int> num_points(num_dimensions);
    int num_total = 1; // points here are a subset of the total, hence there is no chance for overflow in int
    for(int j=0; j<num_dimensions; j++){ num_points[j] = rule->getNumPoints(levels[j]); num_total *= num_points[j]; }

    UnsortedIndexSet* uset = new UnsortedIndexSet(num_dimensions, num_total);
    std::vector<int> p(num_dimensions);

    for(int i=0; i<num_total; i++){
        int t = i;
        for(int j=num_dimensions-1; j>=0; j--){
            p[j] = rule->getPointIndex(levels[j], t % num_points[j]);
            t /= num_points[j];
        }
        uset->addIndex(p.data());
    }

    return uset;
}

IndexSet* IndexManipulator::generateGenericPoints(const IndexSet *tensors, const OneDimensionalWrapper *rule) const{
    int num_tensors = tensors->getNumIndexes();
    UnsortedIndexSet **usets = new UnsortedIndexSet*[num_tensors];

    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<num_tensors; i++){
        usets[i] = tensorGenericPoints(tensors->getIndex(i), rule);
    }

    int warp = num_tensors / 2 + num_tensors % 2;
    IndexSet **psets = new IndexSet*[warp];
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<num_tensors/2; i++){
        psets[i] = new IndexSet(usets[2*i]);
        delete usets[2*i];
        psets[i]->addUnsortedSet(usets[2*i+1]);
        delete usets[2*i+1];
    }
    if (num_tensors % 2 == 1){
        psets[num_tensors/2] = new IndexSet(usets[num_tensors-1]);
        delete usets[num_tensors-1];
    }
    delete[] usets;

    while(warp > 1){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<warp-1; i+=2){
            psets[i]->addIndexSet(psets[i+1]);
            delete psets[i+1];
        }

        for(int i=1; i<warp/2; i++){
            psets[i] = psets[2*i];
        }
        if (warp % 2 == 1){ psets[warp/2] = psets[warp-1];  };
        warp = warp / 2 + warp % 2;
    }

    IndexSet *result = psets[0];
    delete[] psets;
    return result;
}

IndexSet* IndexManipulator::removeIndexesByLimit(IndexSet *iset, const std::vector<int> &limits) const{
    size_t c = 0, dims = iset->getNumDimensions();
    for(int i=0; i<iset->getNumIndexes(); i++){
        const int *idx = iset->getIndex(i);
        bool obeys = true;
        for(size_t j=0; j<dims; j++) if ((limits[j] > -1) && (idx[j] > limits[j])) obeys = false;
        if (obeys) c++;
    }
    if (c == (size_t) iset->getNumIndexes()) return 0;
    std::vector<int> new_idx(dims * c);
    c = 0;
    for(int i=0; i<iset->getNumIndexes(); i++){
        const int *idx = iset->getIndex(i);
        bool obeys = true;
        for(size_t j=0; j<dims; j++) if ((limits[j] > -1) && (idx[j] > limits[j])) obeys = false;
        if (obeys) std::copy(idx, idx + dims, &(new_idx[dims * (c++)]));
    }

    return new IndexSet((int) dims, new_idx);
}

void IndexManipulator::computeDAGup(const IndexSet *iset, Data2D<int> &parents) const{
    int n = iset->getNumIndexes();
    parents.resize(num_dimensions, n);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<n; i++){
        const int *p = iset->getIndex(i);
        std::vector<int> dad(num_dimensions);
        std::copy(p, p + num_dimensions, dad.data());
        int *v = parents.getStrip(i);
        for(auto &d : dad){
            d--;
            *v = (d < 0) ? -1 : iset->getSlot(dad.data());
            d++;
            v++;
        }
    }
}

IndexSet* IndexManipulator::tensorNestedPoints(const int levels[], const OneDimensionalWrapper *rule) const{
    std::vector<int> num_points(num_dimensions);
    std::vector<int> offsets(num_dimensions, 0);
    int num_total = 1; // this is only a subset of all points, no danger of overflow
    for(int j=0; j<num_dimensions; j++){
        if (levels[j] > 0) offsets[j] = rule->getNumPoints(levels[j]-1);
        num_points[j] = rule->getNumPoints(levels[j]) - offsets[j];
        num_total *= num_points[j];
    }
    std::vector<int> index(((size_t) num_dimensions) * ((size_t) num_total));

    auto iter = index.rbegin();
    for(int i=num_total-1; i>=0; i--){
        int t = i;
        auto n = num_points.rbegin();
        auto o = offsets.rbegin();
        for(int j=num_dimensions-1; j>=0; j--){
            *iter++ = *o++ + (t % *n);
            t /= *n++;
        }
    }

    return new IndexSet(num_dimensions, index);
}

IndexSet* IndexManipulator::generateNestedPoints(const IndexSet *tensors, const OneDimensionalWrapper *rule) const{
    int num_tensors = tensors->getNumIndexes();
    IndexSet **psets = new IndexSet*[num_tensors];
    for(int i=0; i<num_tensors; i++){
        psets[i] = tensorNestedPoints(tensors->getIndex(i), rule);
    }

    int warp = num_tensors;
    while(warp > 1){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<warp-1; i+=2){
            psets[i]->addIndexSet(psets[i+1]);
            delete psets[i+1];
        }

        for(int i=1; i<warp/2; i++){
            psets[i] = psets[2*i];
        }
        if (warp % 2 == 1){ psets[warp/2] = psets[warp-1];  };
        warp = warp / 2 + warp % 2;
    }

    IndexSet *result = psets[0];
    delete[] psets;
    return result;
}

IndexSet* IndexManipulator::getPolynomialSpace(const IndexSet *tensors, TypeOneDRule rule, bool iexact) const{
    // may optimize this to consider deltas only
    int num_tensors = tensors->getNumIndexes();
    IndexSet **sets = new IndexSet*[num_tensors];

    for(int t=0; t<num_tensors; t++){
        const int *ti = tensors->getIndex(t);
        std::vector<int> num_points(num_dimensions);
        if (iexact){
            for(int j=0; j<num_dimensions; j++) num_points[j] = meta.getIExact(ti[j], rule) + 1;
        }else{
            for(int j=0; j<num_dimensions; j++) num_points[j] = meta.getQExact(ti[j], rule) + 1;
        }

        int num_total = 1;
        for(auto n : num_points) num_total *= n;

        std::vector<int> poly(((size_t) num_total) * ((size_t) num_dimensions));

        auto iter = poly.rbegin();
        for(int i=num_total-1; i>=0; i--){
            int v = i;
            auto n = num_points.rbegin();
            for(int j=num_dimensions-1; j>=0; j--){
                *iter++ = v % *n;
                v /= *n++;
            }
        }

        sets[t] = new IndexSet(num_dimensions, poly);
    }

    int warp = num_tensors;
    while(warp > 1){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<warp-1; i+=2){
            sets[i]->addIndexSet(sets[i+1]);
            delete sets[i+1];
        }

        for(int i=1; i<warp/2; i++){
            sets[i] = sets[2*i];
        }
        if (warp % 2 == 1){ sets[warp/2] = sets[warp-1]; };
        warp = warp / 2 + warp % 2;
    }

    IndexSet *s = sets[0];
    delete[] sets;

    return s;
}

int IndexManipulator::getMinChildLevel(const IndexSet *iset, TypeDepth type, const std::vector<int> &weights, TypeOneDRule rule){
    int n = iset->getNumIndexes();

    std::vector<int> w;
    getProperWeights(type, weights, w);

    long long min_level = -1;

    #pragma omp parallel
    {
        long long local_level = -1;
        std::vector<int> kid(num_dimensions);

        #pragma omp for
        for(int i=0; i<n; i++){
            const int* p = iset->getIndex(i);
            std::copy(p, p + num_dimensions, kid.data());
            for(auto &k : kid){
                k++;
                if (iset->getSlot(kid) == -1){
                    long long l = getIndexWeight(kid, type, w, rule);
                    if ((local_level == -1) || (l < local_level)) local_level = l;
                }
                k--;
            }
        }

        #pragma omp critical
        {
            if ((min_level == -1) || ((local_level != -1) && (local_level<min_level))) min_level = local_level;
        }
    }

    long long min_weight = weights[0]; for(int j=1; j<num_dimensions; j++) if (min_weight > weights[j]) min_weight = weights[j];
    min_level /= min_weight;
    return (int) min_level;
}

IndexSet* IndexManipulator::selectFlaggedChildren(const IndexSet *iset, const std::vector<bool> &flagged, const std::vector<int> &level_limits) const{
    GranulatedIndexSet *next_level = new GranulatedIndexSet(num_dimensions);
    std::vector<int> kid(num_dimensions);

    int n = iset->getNumIndexes();
    if (level_limits.empty()){
        for(int i=0; i<n; i++){
            if (flagged[i]){
                const int* p = iset->getIndex(i);
                std::copy(p, p + num_dimensions, kid.data());
                for(auto &k : kid){
                    k++;
                    if (iset->getSlot(kid.data()) == -1){
                        next_level->addIndex(kid);
                    }
                    k--;
                }
            }
        }
    }else{
        for(int i=0; i<n; i++){
            if (flagged[i]){
                const int* p = iset->getIndex(i);
                std::copy(p, p + num_dimensions, kid.data());
                for(int j=0; j<num_dimensions; j++){
                    kid[j]++;
                    if (((level_limits[j] == -1) || (kid[j] <= level_limits[j])) && (iset->getSlot(kid.data()) == -1)){
                        next_level->addIndex(kid);
                    }
                    kid[j]--;
                }
            }
        }
    }
    if (next_level->getNumIndexes() == 0){
        delete next_level;
        return 0;
    }

    IndexSet *result = new IndexSet(next_level);
    delete next_level;

    return result;
}

void IndexManipulator::getMaxLevels(const IndexSet *iset, std::vector<int> &max_levels, int &total_max) const{
    max_levels.resize(num_dimensions, 0);
    int n = iset->getNumIndexes();
    for(int i=1; i<n; i++){
        const int* t = iset->getIndex(i);
        for(int j=0; j<num_dimensions; j++) if (max_levels[j] < t[j]) max_levels[j] = t[j];
    }
    total_max = 0;
    for(auto m : max_levels) if (total_max < m) total_max = m;
}
int IndexManipulator::getMaxLevel(const IndexSet *iset) const{
    int n = iset->getNumIndexes();
    int m = 0;
    for(int i=0; i<n; i++){
        const int* t = iset->getIndex(i);
        for(int j=0; j<num_dimensions; j++) if (m < t[j]) m = t[j];
    }
    return m;
}

IndexSet* IndexManipulator::generatePointsFromDeltas(const IndexSet* deltas, std::function<int(int)> getNumPoints) const{
    int num_points = 0;
    for(int i=0; i<deltas->getNumIndexes(); i++){
        const int *p = deltas->getIndex(i);
        int c = 1;
        for(int j=0; j<num_dimensions; j++){
            if (p[j] > 0){
                c *= (getNumPoints(p[j]) - getNumPoints(p[j]-1));
            }else{
                c *= getNumPoints(0);
            }
        }
        num_points += c;
    }

    UnsortedIndexSet *raw_points = new UnsortedIndexSet(num_dimensions, num_points);

    std::vector<int> num_points_delta(num_dimensions);
    std::vector<int> offsets(num_dimensions);
    std::vector<int> index(num_dimensions);

    for(int i=0; i<deltas->getNumIndexes(); i++){
        const int *p = deltas->getIndex(i);
        int num_total = 1;
        for(int j=0; j<num_dimensions; j++){
            num_points_delta[j] = getNumPoints(p[j]);
            if (p[j] > 0){
                offsets[j] = getNumPoints(p[j]-1);
                num_points_delta[j] -= offsets[j];
            }else{
                offsets[j] = 0;
            }

            num_total *= num_points_delta[j];
        }

        for(int k=0; k<num_total; k++){
            int t = k;
            for(int j=num_dimensions-1; j>=0; j--){
                index[j] = offsets[j] + t % num_points_delta[j];
                t /= num_points_delta[j];
            }
            raw_points->addIndex(index.data());
        }
    }

    IndexSet* iset = new IndexSet(raw_points);
    delete raw_points;

    return iset;
}

void IndexManipulator::computeDAGupLocal(const IndexSet *iset, const BaseRuleLocalPolynomial *rule, Data2D<int> &parents) const{
    int n = iset->getNumIndexes();
    if (rule->getMaxNumParents() > 1){ // allow for multiple parents and level 0 may have more than one node
        int max_parents = rule->getMaxNumParents()*num_dimensions;
        parents.resize(max_parents, n, -1);
        int level0_offset = rule->getNumPoints(0);
        #pragma omp parallel for schedule(static)
        for(int i=0; i<n; i++){
            const int *p = iset->getIndex(i);
            std::vector<int> dad(num_dimensions);
            std::copy(p, p + num_dimensions, dad.data());
            int *pp = parents.getStrip(i);
            for(int j=0; j<num_dimensions; j++){
                if (dad[j] >= level0_offset){
                    int current = p[j];
                    dad[j] = rule->getParent(current);
                    pp[2*j] = iset->getSlot(dad);
                    while ((dad[j] >= level0_offset) && (pp[2*j] == -1)){
                        current = dad[j];
                        dad[j] = rule->getParent(current);
                        pp[2*j] = iset->getSlot(dad);
                    }
                    dad[j] = rule->getStepParent(current);
                    if (dad[j] != -1){
                        pp[2*j + 1] = iset->getSlot(dad.data());
                    }
                    dad[j] = p[j];
                }
            }
        }
    }else{ // this assumes that level zero has only one node
        parents.resize(num_dimensions, n);
        #pragma omp parallel for schedule(static)
        for(int i=0; i<n; i++){
            const int *p = iset->getIndex(i);
            std::vector<int> dad(num_dimensions);
            std::copy(p, p + num_dimensions, dad.data());
            int *pp = parents.getStrip(i);
            for(int j=0; j<num_dimensions; j++){
                if (dad[j] == 0){
                    pp[j] = -1;
                }else{
                    dad[j] = rule->getParent(dad[j]);
                    pp[j] = iset->getSlot(dad.data());
                    while((dad[j] != 0) && (pp[j] == -1)){
                        dad[j] = rule->getParent(dad[j]);
                        pp[j] = iset->getSlot(dad);
                    }
                    dad[j] = p[j];
                }
            }
        }
    }
}

}

#endif
