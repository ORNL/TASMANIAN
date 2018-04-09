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

IndexManipulator::IndexManipulator(int cnum_dimensions, const CustomTabulated* custom) : num_dimensions(cnum_dimensions)
{
    meta = new OneDimensionalMeta(custom);
}
IndexManipulator::~IndexManipulator(){  delete meta;  }

int IndexManipulator::getLevel(const int index[], const int weights[]) const{
    int l = index[0] * weights[0];
    for(int j=1; j<num_dimensions; j++){
        l += index[j] * weights[j];
    }
    return l;
}
int IndexManipulator::getCurved(const int index[], const int weights[]) const{
    int l = index[0] * weights[0];
    double c = weights[0] * log1p(index[0]);
    for(int j=1; j<num_dimensions; j++){
        l += index[j] * weights[j];
        c += weights[j] * log1p(index[j]);
    }
    return ((int)ceil(((double) l) + c));
}
int IndexManipulator::getIPTotal(const int index[], const int weights[], TypeOneDRule rule) const{
    int l = ((index[0] > 0) ? (meta->getIExact(index[0]-1, rule) + 1) : 0) * weights[0];
    for(int j=1; j<num_dimensions; j++){
        l += ((index[j] > 0) ? (meta->getIExact(index[j]-1, rule) + 1) : 0) * weights[j];
    }
    return l;
}
int IndexManipulator::getIPCurved(const int index[], const int weights[], TypeOneDRule rule) const{
    int pex = (index[0] > 0) ? meta->getIExact(index[0]-1, rule) + 1 : 0;
    int l = pex * weights[0];
    double c = (weights[num_dimensions]) * log1p((double) (pex));
    //cout << "c = " << c << "   " << weights[num_dimensions] << "    " << pex << "   " << log1p((double) (pex)) << endl;
    for(int j=1; j<num_dimensions; j++){
        pex = (index[j] > 0) ? meta->getIExact(index[j]-1, rule) + 1 : 0;
        l += pex * weights[j];
        c += (weights[num_dimensions+j]) * log1p((double) (pex));
    }
    return ((int)ceil(((double) l) + c));
}
int IndexManipulator::getQPTotal(const int index[], const int weights[], TypeOneDRule rule) const{
    int l = ((index[0] > 0) ? (meta->getQExact(index[0]-1, rule) + 1) : 0) * weights[0];
    for(int j=1; j<num_dimensions; j++){
        l += ((index[j] > 0) ? (meta->getQExact(index[j]-1, rule) + 1) : 0) * weights[j];
    }
    return l;
}
int IndexManipulator::getQPCurved(const int index[], const int weights[], TypeOneDRule rule) const{
    int pex = (index[0] > 0) ? meta->getQExact(index[0]-1, rule) + 1 : 0;
    int l = pex * weights[0];
    double c = ((double) weights[num_dimensions]) * log((double) (pex+1));
    for(int j=1; j<num_dimensions; j++){
        pex = (index[j] > 0) ? meta->getQExact(index[j]-1, rule) + 1 : 0;
        l += pex * weights[j];
        c += (weights[num_dimensions+j]) * log((double) (pex+1));
    }
    return ((int)ceil(((double) l) + c));
}
int IndexManipulator::getHyperbolic(const int index[], const int weights[]) const{
    double l = pow((double) (index[0]+1), ((double) weights[0]) / ((double) weights[num_dimensions]));
    for(int j=1; j<num_dimensions; j++){
        l *= pow((double) (index[j]+1), ((double) weights[j]) / ((double) weights[num_dimensions]));
    }
    return ((int) ceil(l));
}
int IndexManipulator::getIPHyperbolic(const int index[], const int weights[], TypeOneDRule rule) const{
    double l = (index[0] > 0) ? pow((double) (meta->getIExact(index[0]-1, rule) + 2), ((double) weights[0]) / ((double) weights[num_dimensions])) : 1.0;
    for(int j=1; j<num_dimensions; j++){
        l *= (index[j] > 0) ? pow((double) (meta->getIExact(index[j]-1, rule) + 2), ((double) weights[j]) / ((double) weights[num_dimensions])) : 1.0;
    }
    return ((int) ceil(l));
}
int IndexManipulator::getQPHyperbolic(const int index[], const int weights[], TypeOneDRule rule) const{
    double l = (index[0] > 0) ? pow((double) (meta->getQExact(index[0]-1, rule) + 2), ((double) weights[0]) / ((double) weights[num_dimensions])) : 1.0;
    for(int j=1; j<num_dimensions; j++){
        l *= (index[j] > 0) ? pow((double) (meta->getQExact(index[j]-1, rule) + 2), ((double) weights[j]) / ((double) weights[num_dimensions])) : 1.0;
    }
    return ((int) ceil(l));
}

int IndexManipulator::getIndexWeight(const int index[], TypeDepth type, const int weights[], TypeOneDRule rule) const{
    switch (type){
        case type_level:         return getLevel(index, weights);
        case type_curved:        return getCurved(index, weights);
        case type_iptotal:       return getIPTotal(index, weights, rule);
        case type_ipcurved:      return getIPCurved(index, weights, rule);
        case type_qptotal:       return getQPTotal(index, weights, rule);
        case type_qpcurved:      return getQPCurved(index, weights, rule);
        case type_hyperbolic:    return getHyperbolic(index, weights);
        case type_iphyperbolic:  return getIPHyperbolic(index, weights, rule);
        case type_qphyperbolic:  return getQPHyperbolic(index, weights, rule);
        default:
            return 0;
    }
}

IndexSet* IndexManipulator::selectTensors(int offset, TypeDepth type, const int *anisotropic_weights, TypeOneDRule rule) const{
    // construct the minimum tensor set that covers the target_space defined by offset, type, and anisotropic weights
    // consult the manual for the detailed definition of each case
    // used for constructing Global and Sequence grids
    // implements Corrolary 1 of Theorem 1 from
    // Stoyanov, Webster: "A dynamically adaptive sparse grids method forquasi-optimal interpolation of multidimensional functions"
    // Computers & Mathematics with Applications, 71(11):2449–2465, 2016

    // This cheats a bit, but it handles the special case when we want a full tensor grid
    // instead of computing the grid in the standard gradual level by level way,
    // the tensor grid ca be computed directly with two for-loops
    if ((type == type_tensor) || (type == type_iptensor) || (type == type_qptensor)){ // special case, full tensor
        int *levels = new int[num_dimensions]; // how many levels to keep in each direction
        if (type == type_tensor){
            for(int j=0; j<num_dimensions; j++){
                levels[j] = offset;
                if (anisotropic_weights != 0) levels[j] *= anisotropic_weights[j];
            }
        }else if (type == type_iptensor){
            for(int j=0; j<num_dimensions; j++){
                int target = offset;
                if (anisotropic_weights != 0) target *= anisotropic_weights[j];
                levels[j] = 0;
                while(meta->getIExact(levels[j], rule) < target) levels[j]++;
            }
        }else{
            for(int j=0; j<num_dimensions; j++){
                int target = offset;
                if (anisotropic_weights != 0) target *= anisotropic_weights[j];
                levels[j] = 0;
                while(meta->getQExact(levels[j], rule) < target) levels[j]++;
            }
        }
        int num_total = ++levels[0];
        for(int j=1; j<num_dimensions; j++) num_total *= ++levels[j];
        int *index = new int[num_total * num_dimensions];
        for(int i=0; i<num_total; i++){
            int t = i;
            for(int j = num_dimensions-1; j>=0; j--){
                index[i*num_dimensions + j] = t % levels[j];
                t /= levels[j];
            }
        }
        IndexSet *tensors = new IndexSet(num_dimensions, num_total, index);
        delete[] levels;
        return tensors;
    }

    // non-tensor cases
    int *weights = getProperWeights(type, anisotropic_weights);

    // compute normalization and check if heavily curved
    int normalized_offset = weights[0];
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
        int *root = new int[num_dimensions];  std::fill(root, root + num_dimensions, 0);
        DumpIndexSet *index_dump = new DumpIndexSet(num_dimensions, 256);
        while( !(outside && (c == 0)) ){
            if (outside){
                for(int k=c; k<num_dimensions; k++) root[k] = 0;
                c--;
                root[c]++;
            }else{
                index_dump->addIndex(root);
                c = num_dimensions-1;
                root[c]++;
            }
            outside = (getIndexWeight(root, type, weights, rule) > normalized_offset);
        }
        delete[] root;
        int num_loaded = index_dump->getNumLoaded();
        int *res = index_dump->ejectIndexes();
        total = new IndexSet(num_dimensions, num_loaded, res);
        delete index_dump;

    }else{ // use slower algorithm, but more general
        GranulatedIndexSet **sets;
        int *root = new int[num_dimensions];  std::fill(root, root + num_dimensions, 0);
        GranulatedIndexSet *set_level = new GranulatedIndexSet(num_dimensions, 1);  set_level->addIndex(root);  delete[] root;
        total = new IndexSet(num_dimensions);
        bool adding = true;
        while(adding){
            int num_sets = (set_level->getNumIndexes() / 8) + 1;

            sets = new GranulatedIndexSet*[num_sets];

            //#pragma omp parallel for schedule(dynamic)
            for(int me=0; me<num_sets; me++){

                sets[me] = new GranulatedIndexSet(num_dimensions);
                int *newp = new int[num_dimensions];
                for(int i=me; i<set_level->getNumIndexes(); i+=num_sets){
                    const int *p = set_level->getIndex(i);  std::copy(p, p + num_dimensions, newp);
                    for(int j=0; j<num_dimensions; j++){
                        newp[j]++;
                        if ((getIndexWeight(newp, type, weights, rule) <= normalized_offset) && (set_level->getSlot(newp) == -1) && (total->getSlot(newp) == -1)){
                            sets[me]->addIndex(newp);
                        }
                        newp[j]--;
                    }
                }
                delete[] newp;

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

    delete[] weights;

    return total;
}

IndexSet* IndexManipulator::selectTensors(const IndexSet *target_space, bool integration, TypeOneDRule rule) const{
    // construct the minimum tensor set that covers the target_space
    // used when projecting the Global grid onto Lagendre polynomial basis
    // implements the clean version of Theorem 1 from
    // Stoyanov, Webster: "A dynamically adaptive sparse grids method forquasi-optimal interpolation of multidimensional functions"
    // Computers & Mathematics with Applications, 71(11):2449–2465, 2016
    // the other selectTensors() function cover the special caseof Corrolary 1
    GranulatedIndexSet **sets;
    int *root = new int[num_dimensions];  std::fill(root, root + num_dimensions, 0);
    GranulatedIndexSet *set_level = new GranulatedIndexSet(num_dimensions, 1);  set_level->addIndex(root);  delete[] root;
    IndexSet *total = new IndexSet(num_dimensions);
    bool adding = true;

    while(adding){
        int num_sets = (set_level->getNumIndexes() / 8) + 1;

        sets = new GranulatedIndexSet*[num_sets];

        #pragma omp parallel for schedule(dynamic)
        for(int me=0; me<num_sets; me++){

            int *corner = new int[num_dimensions];

            sets[me] = new GranulatedIndexSet(num_dimensions);
            int *newp = new int[num_dimensions];
            for(int i=me; i<set_level->getNumIndexes(); i+=num_sets){
                const int *p = set_level->getIndex(i);  std::copy(p, p + num_dimensions, newp);

                for(int j=0; j<num_dimensions; j++){
                    newp[j]++;

                    // compute the minimum polynomial for this delta
                    if (integration){
                        for(int k=0; k<num_dimensions; k++){
                            corner[k] = (newp[k] > 0) ? (meta->getQExact(newp[k]-1, rule) + 1) : 0;
                        }
                    }else{
                        for(int k=0; k<num_dimensions; k++){
                            corner[k] = (newp[k] > 0) ? (meta->getIExact(newp[k]-1, rule) + 1) : 0;
                        }
                    }

                    if ((target_space->getSlot(corner) != -1) && (set_level->getSlot(newp) == -1) && (total->getSlot(newp) == -1)){
                        sets[me]->addIndex(newp);
                    }
                    newp[j]--;
                }

            }
            delete[] newp;

            delete[] corner;
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

IndexSet* IndexManipulator::getLowerCompletion(const IndexSet *set) const{

    GranulatedIndexSet *set_level = new GranulatedIndexSet(num_dimensions);
    for(int i=0; i<set->getNumIndexes(); i++){
        const int* p = set->getIndex(i);
        int* dad = new int[num_dimensions];
        std::copy(p, p + num_dimensions, dad);
        for(int j=0; j<num_dimensions; j++){
            dad[j]--;
            if ((dad[j]>-1) && (set->getSlot(dad) == -1)){
                set_level->addIndex(dad);
            }
            dad[j]++;
        }
        delete[] dad;
    }

    if (set_level->getNumIndexes() > 0){
        bool adding = true;
        IndexSet *total = new IndexSet(num_dimensions);
        while(adding){
            total->addGranulatedSet(set_level);
            GranulatedIndexSet *set_future = new GranulatedIndexSet(num_dimensions);
            for(int i=0; i<set_level->getNumIndexes(); i++){
                const int* p = set_level->getIndex(i);
                int* dad = new int[num_dimensions];
                std::copy(p, p + num_dimensions, dad);
                for(int j=0; j<num_dimensions; j++){
                    dad[j]--;
                    if ((dad[j]>-1) && (set->getSlot(dad) == -1) && (total->getSlot(dad) == -1)){
                        set_future->addIndex(dad);
                    }
                    dad[j]++;
                }
                delete[] dad;
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

int* IndexManipulator::getProperWeights(TypeDepth type, const int *anisotropic_weights) const{
    int *weights;
    if (anisotropic_weights == 0){
        if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
            weights = new int[2*num_dimensions];
            std::fill(weights, weights + num_dimensions, 1);
            std::fill(&(weights[num_dimensions]), &(weights[num_dimensions]) + num_dimensions, 0);
        }else if ((type == type_hyperbolic) || (type == type_iphyperbolic) || (type == type_qphyperbolic)){
            weights = new int[num_dimensions+1];
            std::fill(weights, weights + num_dimensions+1, 1);
        }else{
            weights = new int[num_dimensions];
            std::fill(weights, weights + num_dimensions, 1);
        }
    }else{
        if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
            weights = new int[2*num_dimensions];
            std::copy(anisotropic_weights, anisotropic_weights + 2*num_dimensions, weights);
        }else if ((type == type_hyperbolic) || (type == type_iphyperbolic) || (type == type_qphyperbolic)){
            weights = new int[num_dimensions+1];
            std::copy(anisotropic_weights, anisotropic_weights + num_dimensions, weights);
            weights[num_dimensions] = weights[0];
            for(int j=1; j<num_dimensions; j++){ weights[num_dimensions] += weights[j]; }
        }else{
            weights = new int[num_dimensions];
            std::copy(anisotropic_weights, anisotropic_weights + num_dimensions, weights);
        }
    }
    return weights;
}

int* IndexManipulator::computeLevels(const IndexSet* set) const{
    int num_indexes = set->getNumIndexes();
    int *level = new int[num_indexes];
    #pragma omp parallel for
    for(int i=0; i<num_indexes; i++){
        const int* p = set->getIndex(i);
        level[i] = p[0];
        for(int j=1; j<num_dimensions; j++){
            level[i] += p[j];
        }
    }
    return level;
}

int* IndexManipulator::makeTensorWeights(const IndexSet* set) const{
    int n = set->getNumIndexes();
    int *level = computeLevels(set);
    int *weights = new int[n];

    int max_level = level[0];  for(int i=1; i<n; i++){  if (max_level < level[i]) max_level = level[i];  }

    int *kids = new int[n * num_dimensions];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<n; i++){
        const int *p = set->getIndex(i);
        int *kid = new int[num_dimensions];
        std::copy(p, p + num_dimensions, kid);

        for(int j=0; j<num_dimensions; j++){
            kid[j]++;
            kids[i*num_dimensions + j] = set->getSlot(kid);
            kid[j]--;
        }
        delete[] kid;
        if (level[i] == max_level){
            weights[i] = 1;
        }
    }

    for(int l=max_level-1; l>=0; l--){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<n; i++){
            if (level[i] == l){
                int *monkey_tail  = new int[max_level-l+1];
                int *monkey_count = new int[max_level-l+1];
                bool* used = new bool[n];
                std::fill(used, used + n, false);

                int current = 0;
                monkey_count[0] = 0;
                monkey_tail[0] = i;

                int sum = 0;

                while(monkey_count[0] < num_dimensions){
                    if (monkey_count[current] < num_dimensions){
                        int branch = kids[monkey_tail[current]*num_dimensions + monkey_count[current]];
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

                delete[] used;
                delete[] monkey_count;
                delete[] monkey_tail;
            }
        }
    }

    delete[] kids;
    delete[] level;
    return weights;
}

IndexSet* IndexManipulator::nonzeroSubset(const IndexSet* set, const int weights[]) const{
    int nz_weights = 0;
    for(int i=0; i<set->getNumIndexes(); i++){ if (weights[i] != 0) nz_weights++; }

    int *index = new int[nz_weights * num_dimensions];
    nz_weights = 0;
    for(int i=0; i<set->getNumIndexes(); i++){
        if (weights[i] != 0){
            const int *p = set->getIndex(i);
            std::copy(p, p + num_dimensions, &(index[num_dimensions * nz_weights++]));
        }
    }

    IndexSet *result = new IndexSet(num_dimensions, nz_weights, index);
    return result;
}

UnsortedIndexSet* IndexManipulator::tensorGenericPoints(const int levels[], const OneDimensionalWrapper *rule) const{
    int *num_points = new int[num_dimensions];
    int num_total = 1;
    for(int j=0; j<num_dimensions; j++){  num_points[j] = rule->getNumPoints(levels[j]); num_total *= num_points[j];  }

    UnsortedIndexSet* set = new UnsortedIndexSet(num_dimensions, num_total);
    int *p = new int[num_dimensions];

    for(int i=0; i<num_total; i++){
        int t = i;
        for(int j=num_dimensions-1; j>=0; j--){
            p[j] = rule->getPointIndex(levels[j], t % num_points[j]);
            t /= num_points[j];
        }
        set->addIndex(p);
    }

    delete[] p;
    delete[] num_points;

    return set;
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

int* IndexManipulator::referenceGenericPoints(const int levels[], const OneDimensionalWrapper *rule, const IndexSet *points) const{
    int *num_points = new int[num_dimensions];
    int num_total = 1;
    for(int j=0; j<num_dimensions; j++){  num_points[j] = rule->getNumPoints(levels[j]); num_total *= num_points[j];  }

    int* refs = new int[num_total];
    int *p = new int[num_dimensions];

    for(int i=0; i<num_total; i++){
        int t = i;
        for(int j=num_dimensions-1; j>=0; j--){
            p[j] = rule->getPointIndex(levels[j], t % num_points[j]);
            t /= num_points[j];
        }
        refs[i] = points->getSlot(p);
    }

    delete[] p;
    delete[] num_points;

    return refs;
}

IndexSet* IndexManipulator::removeIndexesByLimit(IndexSet *set, const int limits[]) const{
    int c = 0, dims = set->getNumDimensions();
    for(int i=0; i<set->getNumIndexes(); i++){
        const int *idx = set->getIndex(i);
        bool obeys = true;
        for(int j=0; j<dims; j++) if ((limits[j] > -1) && (idx[j] > limits[j])) obeys = false;
        if (obeys) c++;
    }
    if (c == set->getNumIndexes()) return 0;
    int *new_idx = new int[dims * c];
    c = 0;
    for(int i=0; i<set->getNumIndexes(); i++){
        const int *idx = set->getIndex(i);
        bool obeys = true;
        for(int j=0; j<dims; j++) if ((limits[j] > -1) && (idx[j] > limits[j])) obeys = false;
        if (obeys) std::copy(idx, idx + dims, &(new_idx[dims * (c++)]));
    }
    return (new IndexSet(dims, c, new_idx));
}
UnsortedIndexSet* IndexManipulator::removeIndexesByLimit(UnsortedIndexSet *set, const int limits[]) const{
    int c = 0, dims = set->getNumDimensions();
    for(int i=0; i<set->getNumIndexes(); i++){
        const int *idx = set->getIndex(i);
        bool obeys = true;
        for(int j=0; j<dims; j++) if ((limits[j] > -1) && (idx[j] > limits[j])) obeys = false;
        if (obeys) c++;
    }
    if (c == set->getNumIndexes()) return 0;
    UnsortedIndexSet* restricted = new UnsortedIndexSet(dims, c);
    for(int i=0; i<set->getNumIndexes(); i++){
        const int *idx = set->getIndex(i);
        bool obeys = true;
        for(int j=0; j<dims; j++) if ((limits[j] > -1) && (idx[j] > limits[j])) obeys = false;
        if (obeys) restricted->addIndex(idx);
    }
    return restricted;
}

int* IndexManipulator::computeDAGup(const IndexSet *set) const{
    int n = set->getNumIndexes();
    int *parents = new int[n * num_dimensions];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<n; i++){
        const int *p = set->getIndex(i);
        int *dad = new int[num_dimensions];
        std::copy(p, p + num_dimensions, dad);
        for(int j=0; j<num_dimensions; j++){
            dad[j]--;
            parents[i*num_dimensions + j] = (dad[j] < 0) ? -1 : set->getSlot(dad);
            dad[j]++;
        }
        delete[] dad;
    }
    return parents;
}

IndexSet* IndexManipulator::tensorNestedPoints(const int levels[], const OneDimensionalWrapper *rule) const{
    int *num_points = new int[num_dimensions];
    int *offsets = new int[num_dimensions];
    int num_total = 1;
    for(int j=0; j<num_dimensions; j++){
        num_points[j] = rule->getNumPoints(levels[j]);
        if (levels[j] > 0){
            offsets[j] = rule->getNumPoints(levels[j]-1);
            num_points[j] -= offsets[j];
        }else{
            offsets[j] = 0;
        }

        num_total *= num_points[j];
    }
    int *index = new int[num_dimensions * num_total];

    for(int i=0; i<num_total; i++){
        int t = i;
        for(int j=num_dimensions-1; j>=0; j--){
            index[i*num_dimensions + j] = offsets[j] + t % num_points[j];
            t /= num_points[j];
        }
    }

    delete[] offsets;
    delete[] num_points;

    IndexSet* set = new IndexSet(num_dimensions, num_total, index);
    return set;
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

int* IndexManipulator::referenceNestedPoints(const int levels[], const OneDimensionalWrapper *rule, const IndexSet *points) const{
    int *num_points = new int[num_dimensions];
    int num_total = 1;
    for(int j=0; j<num_dimensions; j++){  num_points[j] = rule->getNumPoints(levels[j]); num_total *= num_points[j];  }

    int* refs = new int[num_total];
    int *p = new int[num_dimensions];

    for(int i=0; i<num_total; i++){
        int t = i;
        for(int j=num_dimensions-1; j>=0; j--){
            p[j] = t % num_points[j];
            t /= num_points[j];
        }
        refs[i] = points->getSlot(p);
    }

    delete[] p;
    delete[] num_points;

    return refs;
}

IndexSet* IndexManipulator::getPolynomialSpace(const IndexSet *tensors, TypeOneDRule rule, bool iexact) const{
    // may optimize this to consider deltas only
    int num_tensors = tensors->getNumIndexes();
    IndexSet **sets = new IndexSet*[num_tensors];

    for(int t=0; t<num_tensors; t++){
        const int *ti = tensors->getIndex(t);
        int *num_points = new int[num_dimensions];
        if (iexact){
            for(int j=0; j<num_dimensions; j++)  num_points[j] = meta->getIExact(ti[j], rule) + 1;
        }else{
            for(int j=0; j<num_dimensions; j++)  num_points[j] = meta->getQExact(ti[j], rule) + 1;
        }

        int num_total = num_points[0];
        for(int j=1; j<num_dimensions; j++) num_total *= num_points[j];

        int *poly = new int[num_total * num_dimensions];

        for(int i=0; i<num_total; i++){
            int v = i;
            for(int j=num_dimensions-1; j>=0; j--){
                poly[i*num_dimensions + j] = v % num_points[j];
                v /= num_points[j];
            }
        }

        sets[t] = new IndexSet(num_dimensions, num_total, poly);
        delete[] num_points;
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
        if (warp % 2 == 1){ sets[warp/2] = sets[warp-1];  };
        warp = warp / 2 + warp % 2;
    }

    IndexSet *s = sets[0];
    delete[] sets;

    return s;
}

int IndexManipulator::getMinChildLevel(const IndexSet *set, TypeDepth type, const int weights[], TypeOneDRule rule){
    int n = set->getNumIndexes();

    int min_level = -1;

    #pragma omp parallel
    {
        int local_level = -1;
        int *kid = new int[num_dimensions];

        #pragma omp for
        for(int i=0; i<n; i++){
            const int* p = set->getIndex(i);
            std::copy(p, p + num_dimensions, kid);
            for(int j=0; j<num_dimensions; j++){
                kid[j]++;
                if (set->getSlot(kid) == -1){
                    int l = getIndexWeight(kid, type, weights, rule);
                    if ((local_level == -1) || (l<local_level)) local_level = l;
                }
                kid[j]--;
            }
        }

        #pragma omp critical
        {
            if ((min_level == -1) || ((local_level!=-1) && (local_level<min_level))) min_level = local_level;
        }

        delete[] kid;
    }

    int min_weight = weights[0];  for(int j=1; j<num_dimensions; j++) if (min_weight > weights[j]) min_weight = weights[j];
    min_level /= min_weight;
    return min_level;
}

IndexSet* IndexManipulator::selectFlaggedChildren(const IndexSet *set, const bool flagged[], const int *level_limits) const{
    GranulatedIndexSet *next_level = new GranulatedIndexSet(num_dimensions);
    int *kid = new int[num_dimensions];

    if (level_limits == 0){
        for(int i=0; i<set->getNumIndexes(); i++){
            if (flagged[i]){
                const int* p = set->getIndex(i);
                std::copy(p, p + num_dimensions, kid);
                for(int j=0; j<num_dimensions; j++){
                    kid[j]++;
                    if (set->getSlot(kid) == -1){
                        next_level->addIndex(kid);
                    }
                    kid[j]--;
                }
            }
        }
    }else{
        for(int i=0; i<set->getNumIndexes(); i++){
            if (flagged[i]){
                const int* p = set->getIndex(i);
                std::copy(p, p + num_dimensions, kid);
                for(int j=0; j<num_dimensions; j++){
                    kid[j]++;
                    if (((level_limits[j] == -1) || (kid[j] <= level_limits[j])) && (set->getSlot(kid) == -1)){
                        next_level->addIndex(kid);
                    }
                    kid[j]--;
                }
            }
        }
    }
    delete[] kid;
    if (next_level->getNumIndexes() == 0){
        delete next_level;
        return 0;
    }

    IndexSet *result = new IndexSet(next_level);
    delete next_level;

    return result;
}

void IndexManipulator::getMaxLevels(const IndexSet *set, int max_levels[], int &total_max) const{
    int * ml = new int[num_dimensions];
    std::fill(ml, ml + num_dimensions, 0);
    for(int i=1; i<set->getNumIndexes(); i++){
        const int* t = set->getIndex(i);
        for(int j=0; j<num_dimensions; j++){  if (t[j] > ml[j]) ml[j] = t[j];  }
    }
    total_max = ml[0];
    for(int j=1; j<num_dimensions; j++){  if (total_max < ml[j]) total_max = ml[j];  }
    if (max_levels != 0){
        std::copy(ml, ml+num_dimensions, max_levels);
    }
    delete[] ml;
}

UnsortedIndexSet* IndexManipulator::getToalDegreeDeltas(int level) const{
    int num_delta = 1;
    if (level > 0){
        int k = (level > num_dimensions) ? num_dimensions : level;
        int n = level + num_dimensions;
        for(int i=0; i<k; i++){
            num_delta *= (n - i);
            num_delta /= (i + 1);
        }
    }

    UnsortedIndexSet* deltas = new UnsortedIndexSet(num_dimensions, num_delta);
    int *index = new int[num_dimensions]; std::fill(index, index + num_dimensions, 0);
    deltas->addIndex(index);

    if (num_dimensions == 1){
        for(int l=1; l<=level; l++) deltas->addIndex(&l);
    }else{
        int *remainder = new int[num_dimensions];
        for(int l=1; l<=level; l++){
            remainder[0] = l;
            index[0] = 0;
            int current = 0;
            while(index[0] <= remainder[0]){
                if (index[current] > remainder[current]){
                    index[--current]++;
                }else{
                    if (current == num_dimensions-2){
                        index[current+1] = remainder[current] - index[current];
                        deltas->addIndex(index);
                        index[current]++;
                    }else{
                        remainder[current+1] = remainder[current] - index[current];
                        index[++current] = 0;
                    }
                }
            }
        }
        delete[] remainder;
    }
    delete[] index;

    return deltas;
}

IndexSet* IndexManipulator::generatePointsFromDeltas(const UnsortedIndexSet* deltas, const BaseRuleLocalPolynomial *rule) const{
    int num_points = 0;
    for(int i=0; i<deltas->getNumIndexes(); i++){
        const int *p = deltas->getIndex(i);
        int c = 1;
        for(int j=0; j<num_dimensions; j++){
            if (p[j] > 0){
                c *= (rule->getNumPoints(p[j]) - rule->getNumPoints(p[j]-1));
            }else{
                c *= rule->getNumPoints(0);
            }
        }
        num_points += c;
    }

    UnsortedIndexSet *raw_points = new UnsortedIndexSet(num_dimensions, num_points);

    int *num_points_delta = new int[num_dimensions];
    int *offsets = new int[num_dimensions];
    int *index = new int[num_dimensions];

    for(int i=0; i<deltas->getNumIndexes(); i++){
        const int *p = deltas->getIndex(i);
        int num_total = 1;
        for(int j=0; j<num_dimensions; j++){
            num_points_delta[j] = rule->getNumPoints(p[j]);
            if (p[j] > 0){
                offsets[j] = rule->getNumPoints(p[j]-1);
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
            raw_points->addIndex(index);
        }
    }

    delete[] index;
    delete[] offsets;
    delete[] num_points_delta;

    IndexSet* set = new IndexSet(raw_points);
    delete raw_points;

    return set;
}

int* IndexManipulator::computeDAGupLocal(const IndexSet *set, const BaseRuleLocalPolynomial *rule) const{
    int n = set->getNumIndexes();
    int *parents;
    if (rule->getMaxNumParents() > 1){
        int max_parents = rule->getMaxNumParents()*num_dimensions;
        parents = new int[n * max_parents];
        #pragma omp parallel for schedule(static)
        for(int i=0; i<n; i++){
            const int *p = set->getIndex(i);
            int *dad = new int[num_dimensions];
            std::copy(p, p + num_dimensions, dad);
            for(int j=0; j<num_dimensions; j++){
                if (dad[j] == 0){
                    parents[i*max_parents + 2*j    ] = -1;
                    parents[i*max_parents + 2*j + 1] = -1;
                }else{
                    int current = p[j];
                    dad[j] = rule->getParent(current);
                    parents[i*max_parents + 2*j] = set->getSlot(dad);
                    while ((dad[j] != 0) && (parents[i*max_parents + 2*j] == -1)){
                        current = dad[j];
                        dad[j] = rule->getParent(current);
                        parents[i*max_parents + 2*j] = set->getSlot(dad);
                    }
                    dad[j] = rule->getStepParent(current);
                    if (dad[j] != -1){
                        parents[i*max_parents + 2*j + 1] = set->getSlot(dad);
                    }else{
                        parents[i*max_parents + 2*j + 1] = -1;
                    }
                    dad[j] = p[j];
                }
            }
            delete[] dad;
        }
    }else{
        parents = new int[n * num_dimensions];
        #pragma omp parallel for schedule(static)
        for(int i=0; i<n; i++){
            const int *p = set->getIndex(i);
            int *dad = new int[num_dimensions];
            std::copy(p, p + num_dimensions, dad);
            for(int j=0; j<num_dimensions; j++){
                if (dad[j] == 0){
                    parents[i*num_dimensions + j] = -1;
                }else{
                    dad[j] = rule->getParent(dad[j]);
                    parents[i*num_dimensions + j] = set->getSlot(dad);
                    while((dad[j] != 0) && (parents[i*num_dimensions + j] == -1)){
                        dad[j] = rule->getParent(dad[j]);
                        parents[i*num_dimensions + j] = set->getSlot(dad);
                    }
                    dad[j] = p[j];
                }
            }
            delete[] dad;
        }
    }

    return parents;
}

}

#endif
