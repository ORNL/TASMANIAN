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

#ifndef __TASMANIAN_SPARSE_GRID_INDEX_SETS_CPP
#define __TASMANIAN_SPARSE_GRID_INDEX_SETS_CPP

#include "tsgIndexSets.hpp"

namespace TasGrid{

UnsortedIndexSet::UnsortedIndexSet(int cnum_dimensions, int cnum_slots) : num_indexes(0){
    num_dimensions = (size_t) cnum_dimensions;
    index.resize(num_dimensions * cnum_slots);
}
UnsortedIndexSet::~UnsortedIndexSet(){}

int UnsortedIndexSet::getNumDimensions() const{ return (int) num_dimensions; }

void UnsortedIndexSet::addIndex(const int p[]){
    std::copy(p, p + num_dimensions, &(index[num_dimensions * num_indexes++]));
}

void UnsortedIndexSet::getIndexesSorted(std::vector<int> &sorted) const{
    std::vector<int> list_source(num_indexes);
    std::vector<int> list_destination(num_indexes);

    for(int i=0; i<(int) (num_indexes - num_indexes%2); i+=2){
        if (compareIndexes(&(index[i*num_dimensions]), &(index[(i+1)*num_dimensions])) == type_abeforeb){
            list_source[i] = i;
            list_source[i+1] = i+1;
        }else{
            list_source[i] = i+1;
            list_source[i+1] = i;
        }
    }
    if (num_indexes%2 == 1) list_source[num_indexes-1] = (int) (num_indexes-1);

    size_t warp = 2;
    while(warp < num_indexes){
        size_t full_warps = 2*warp * ((num_indexes) / (2*warp));
        #pragma omp parallel for
        for(int i=0; i<(int) full_warps; i+=(int) (2*warp)){
            mergeLists(&(list_source[i]), warp, &(list_source[i+warp]), warp, &(list_destination[i]));
        }

        if (full_warps < num_indexes){
            size_t sizeA = (full_warps+warp < num_indexes) ? warp : num_indexes - full_warps; // try to use warp indexes, but do no exceed num_indexes
            size_t sizeB = (full_warps+2*warp < num_indexes) ? warp : num_indexes - full_warps - sizeA; // try to use warp indexes, but do no exceed num_indexes
            mergeLists(&(list_source[full_warps]), sizeA, &(list_source[full_warps+warp]), sizeB, &(list_destination[full_warps]));
        }
        std::swap(list_source, list_destination);
        warp *= 2;
    }

    sorted.resize(num_indexes * num_dimensions);
    for(size_t i=0; i<num_indexes; i++){
        std::copy(&(index[list_source[i] * num_dimensions]), &(index[list_source[i] * num_dimensions]) + num_dimensions, &(sorted[i*num_dimensions]));
    }
}

TypeIndexRelation UnsortedIndexSet::compareIndexes(const int a[], const int b[]) const{
    for(size_t i=0; i<num_dimensions; i++){
        if (a[i] < b[i]) return type_abeforeb;
        if (a[i] > b[i]) return type_bbeforea;
    }
    return type_asameb;
}

void UnsortedIndexSet::mergeLists(const int listA[], size_t sizeA, const int listB[], size_t sizeB, int destination[]) const{
    TypeIndexRelation relation;
    size_t offset_destination = 0, offsetA = 0, offsetB = 0;
    while( (offsetA < sizeA) || (offsetB < sizeB) ){
        if (offsetA >= sizeA){
            relation = type_bbeforea;
        }else if (offsetB >= sizeB){
            relation = type_abeforeb;
        }else{
            relation = compareIndexes(&(index[listA[offsetA]*num_dimensions]), &(index[listB[offsetB]*num_dimensions]));
        }
        if (relation == type_abeforeb){
            destination[offset_destination++] = listA[offsetA++];
        }else{
            destination[offset_destination++] = listB[offsetB++];
        }
    }
}

GranulatedIndexSet::GranulatedIndexSet(int cnum_dimensions, int cnum_slots) : num_dimensions((size_t) cnum_dimensions), num_indexes(0){
    index.reserve(num_dimensions * cnum_slots);
    imap.reserve(cnum_slots);
}

GranulatedIndexSet::~GranulatedIndexSet(){}

int GranulatedIndexSet::getNumDimensions() const{ return (int) num_dimensions; }
int GranulatedIndexSet::getNumIndexes() const{ return (int) num_indexes; }

int GranulatedIndexSet::getSlot(const int p[]) const{
    int sstart = 0, send = (int) (num_indexes - 1);
    int current = (sstart + send) / 2;
    TypeIndexRelation t;
    while (sstart <= send){
        t = compareIndexes(&(index[num_dimensions * imap[current]]), p);
        if (t == type_abeforeb){
            sstart = current+1;
        }else if (t == type_bbeforea){
            send = current-1;
        }else{
            return current;
        };
        current = (sstart + send) / 2;
    }
    return -1;
}

const int* GranulatedIndexSet::getIndex(int j) const{
    return &(index[imap[j] * num_dimensions]);
}

void GranulatedIndexSet::addIndex(const int p[]){
    int sstart = 0, send = (int) (num_indexes - 1); // start and end for the search
    int current = (sstart + send) / 2;
    TypeIndexRelation t;
    while (sstart <= send){
        t = compareIndexes(&(index[num_dimensions * imap[current]]), p);
        if (t == type_abeforeb){
            sstart = current+1;
        }else if (t == type_bbeforea){
            send = current-1;
        }else{
            return;
        };
        current = (sstart + send) / 2;
    }

    for(size_t i=0; i<num_dimensions; i++) index.push_back(p[i]);
    imap.push_back((int) num_indexes);

    size_t slot = sstart+1; slot = (slot > num_indexes) ? num_indexes : slot;
    while((slot>0) && (compareIndexes(p, &(index[imap[slot-1] * num_dimensions])) == type_abeforeb)){ slot--; }

    for(size_t i=num_indexes; i>slot; i--){ imap[i] = imap[i-1]; }
    imap[slot] = (int) num_indexes;
    num_indexes++;
}

void GranulatedIndexSet::addGranulatedSet(const GranulatedIndexSet *gset){
    mergeMapped(*(gset->getIndexes()), *(gset->getMap()), gset->getNumIndexes());
}
const std::vector<int>* GranulatedIndexSet::getIndexes() const{
    return &index;
}
const std::vector<int>* GranulatedIndexSet::getMap() const{
    return &imap;
}

TypeIndexRelation GranulatedIndexSet::compareIndexes(const int a[], const int b[]) const{
    for(size_t i=0; i<num_dimensions; i++){
        if (a[i] < b[i]) return type_abeforeb;
        if (a[i] > b[i]) return type_bbeforea;
    }
    return type_asameb;
}

void GranulatedIndexSet::mergeMapped(const std::vector<int> newIndex, const std::vector<int> newMap, int sizeNew){
    std::vector<int> oldIndex = std::move(index);
    size_t sizeOld = num_indexes;

    index.resize(0);
    index.reserve((sizeOld + sizeNew) * num_dimensions);

    TypeIndexRelation relation;
    num_indexes = 0;
    auto iterNew = newMap.begin();
    auto iterOld = imap.begin();
    while( (iterNew < newMap.end()) || (iterOld < imap.end()) ){
        if (iterNew >= newMap.end()){ // new is a
            relation = type_bbeforea;
        }else if (iterOld >= imap.end()){ // old is b
            relation = type_abeforeb;
        }else{
            relation = compareIndexes(&(newIndex[*iterNew * num_dimensions]), &(oldIndex[*iterOld * num_dimensions]));
        }
        const int *p;
        if (relation == type_abeforeb){
            p = &(newIndex[*iterNew * num_dimensions]);
            iterNew++;
        }else{
            p = &(oldIndex[*iterOld * num_dimensions]);
            iterOld++;
            if (relation == type_asameb) iterNew ++;
        }
        for(size_t i=0; i<num_dimensions; i++) index.push_back(p[i]);
        num_indexes++;
    }

    imap.resize(num_indexes);
    for(size_t i=0; i<num_indexes; i++){ imap[i] = (int) i; }
}

IndexSet::IndexSet(int cnum_dimensions){
    num_dimensions = (size_t) cnum_dimensions;
    cache_num_indexes = 0;
}
IndexSet::IndexSet(const UnsortedIndexSet *uset){
    num_dimensions = (size_t) uset->getNumDimensions();
    uset->getIndexesSorted(index);
    cacheNumIndexes();
}
IndexSet::IndexSet(const GranulatedIndexSet *gset){
    num_dimensions = (size_t) gset->getNumDimensions();
    size_t num_indexes = (size_t) gset->getNumIndexes();
    index.resize(num_dimensions * num_indexes);
    for(size_t i=0; i<num_indexes; i++){
        const int *p = gset->getIndex((int) i);
        std::copy(p, p + num_dimensions, &(index[i*num_dimensions]));
    }
    cacheNumIndexes();
}
IndexSet::IndexSet(const IndexSet *iset){
    num_dimensions = (size_t) iset->getNumDimensions();
    index = *iset->getIndexes();
    cacheNumIndexes();
}

IndexSet::IndexSet(int cnum_dimensions, std::vector<int> &cindex){
    num_dimensions = cnum_dimensions;
    index = std::move(cindex);
    cacheNumIndexes();
}
IndexSet::~IndexSet(){}

int IndexSet::getNumDimensions() const{ return (int) num_dimensions; }
int IndexSet::getNumIndexes() const{ return cache_num_indexes; }

void IndexSet::write(std::ofstream &ofs) const{
    size_t numi = index.size() / num_dimensions;
    ofs << num_dimensions << " " << numi;
    for(auto i : index) ofs << " " << i;
    ofs << std::endl;
}
void IndexSet::read(std::ifstream &ifs){
    size_t numi;
    ifs >> num_dimensions >> numi;
    index.resize(num_dimensions * numi);
    for(auto & i : index) ifs >> i;
    cacheNumIndexes();
}

void IndexSet::writeBinary(std::ofstream &ofs) const{
    int sizes[2];
    sizes[0] = (int) num_dimensions;
    sizes[1] = (int) (index.size() / num_dimensions);
    ofs.write((char*) sizes, 2*sizeof(int));
    ofs.write((char*) index.data(), index.size()*sizeof(int));
}
void IndexSet::readBinary(std::ifstream &ifs){
    int sizes[2];
    ifs.read((char*) sizes, 2*sizeof(int));
    num_dimensions = (size_t) sizes[0];
    size_t num_indexes = (size_t) sizes[1];
    index.resize(num_dimensions * num_indexes);
    ifs.read((char*) index.data(), num_dimensions*num_indexes*sizeof(int));
    cacheNumIndexes();
}

int IndexSet::getSlot(const int p[]) const{
    int sstart = 0, send = (int) (index.size() / num_dimensions - 1);
    int current = (sstart + send) / 2;
    TypeIndexRelation t;
    while (sstart <= send){
        t = compareIndexes(&(index[num_dimensions * current]), p);
        if (t == type_abeforeb){
            sstart = current+1;
        }else if (t == type_bbeforea){
            send = current-1;
        }else{
            return current;
        };
        current = (sstart + send) / 2;
    }
    return -1;
}
const int* IndexSet::getIndex(int i) const{
    return &(index[i*num_dimensions]);
}
const std::vector<int>* IndexSet::getIndexes() const{
    return &index;
}
void IndexSet::addUnsortedSet(const UnsortedIndexSet *uset){
    std::vector<int> uset_index;
    uset->getIndexesSorted(uset_index);
    mergeSet(uset_index);
}
void IndexSet::addGranulatedSet(const GranulatedIndexSet *gset){
    mergeMapped(gset->getIndexes(), gset->getMap());
}
void IndexSet::addIndexSet(const IndexSet *iset){
    mergeSet(iset->index);
}

IndexSet* IndexSet::diffSets(const IndexSet *iset) const{
    int numi = (int) (index.size() / num_dimensions);
    std::vector<int> slots(numi);

    #pragma omp parallel for schedule(static)
    for(int i=0; i<numi; i++){ // openmp needs a signed counter, go figure ... (this is safe conversion)
        slots[i] = iset->getSlot(&(index[i*num_dimensions]));
    }

    size_t n = 0;
    for(auto s : slots){
        if (s == -1) n++;
    }

    if (n == 0){
        return 0;
    }else{
        std::vector<int> diff_index(num_dimensions * n);
        n = 0;
        auto iterIndex = index.begin();
        auto iterDiff = diff_index.begin();
        for(auto s : slots){
            if (s == -1){
                std::copy(iterIndex, iterIndex + num_dimensions, iterDiff);
                std::advance(iterDiff, num_dimensions);
            }
            std::advance(iterIndex, num_dimensions);
        }

        IndexSet *diff_set = new IndexSet((int) num_dimensions, diff_index);
        return diff_set;
    }
}


TypeIndexRelation IndexSet::compareIndexes(const int a[], const int b[]) const{
    for(size_t i=0; i<num_dimensions; i++){
        if (a[i] < b[i]) return type_abeforeb;
        if (a[i] > b[i]) return type_bbeforea;
    }
    return type_asameb;
}
void IndexSet::mergeSet(const std::vector<int> newIndex){
    std::vector<int> oldIndex = std::move(index); // move assignment

    index.reserve(oldIndex.size() + newIndex.size());
    index.resize(0);

    TypeIndexRelation relation;
    auto iterNew = newIndex.begin();
    auto iterOld = oldIndex.begin();
    while( (iterNew < newIndex.end()) || (iterOld < oldIndex.end()) ){
        if (iterNew >= newIndex.end()){ // new is a
            relation = type_bbeforea;
        }else if (iterOld >= oldIndex.end()){ // old is b
            relation = type_abeforeb;
        }else{
            relation = compareIndexes(&*iterNew, &*iterOld);
        }
        if (relation == type_abeforeb){
            for(size_t i=0; i<num_dimensions; i++){
                index.push_back(*iterNew);
                iterNew++;
            }
        }else{
            for(size_t i=0; i<num_dimensions; i++){
                index.push_back(*iterOld);
                iterOld++;
            }
            if (relation == type_asameb) std::advance(iterNew, num_dimensions);
        }
    }
    cacheNumIndexes();
}

void IndexSet::mergeMapped(const std::vector<int> *newIndex, const std::vector<int> *imap){
    std::vector<int> oldIndex = std::move(index); // move assignment

    index.reserve(newIndex->size() + oldIndex.size());
    index.resize(0);

    TypeIndexRelation relation;
    auto iterNew = imap->begin();
    auto iterOld = oldIndex.begin();
    while( (iterNew < imap->end()) || (iterOld < oldIndex.end()) ){
        if (iterNew >= imap->end()){ // new is a
            relation = type_bbeforea;
        }else if (iterOld >= oldIndex.end()){ // old is b
            relation = type_abeforeb;
        }else{
            relation = compareIndexes(&((*newIndex)[(*iterNew) * num_dimensions]), &*iterOld);
        }
        const int *p;
        if (relation == type_abeforeb){
            p = &((*newIndex)[(*iterNew) * num_dimensions]);
            iterNew++;
        }else{
            p = &*iterOld;
            std::advance(iterOld, num_dimensions);
            if (relation == type_asameb) iterNew++;
        }
        for(size_t i=0; i<num_dimensions; i++) index.push_back(p[i]);
    }
    cacheNumIndexes();
}

void IndexSet::cacheNumIndexes(){ cache_num_indexes = (int) (index.size() / num_dimensions); }


MultiIndexSet::MultiIndexSet() : num_dimensions(0), cache_num_indexes(0){}
MultiIndexSet::MultiIndexSet(int cnum_dimensions)  : num_dimensions(cnum_dimensions), cache_num_indexes(0){}
MultiIndexSet::~MultiIndexSet(){}

void MultiIndexSet::write(std::ofstream &ofs) const{
    ofs << num_dimensions << " " << cache_num_indexes;
    for(auto i : indexes) ofs << " " << i;
    ofs << std::endl;
}
void MultiIndexSet::read(std::ifstream &ifs){
    reset();
    ifs >> num_dimensions >> cache_num_indexes;
    indexes.resize(num_dimensions * ((size_t) cache_num_indexes));
    for(auto &i : indexes) ifs >> i;
}

void MultiIndexSet::writeBinary(std::ofstream &ofs) const{
    int sizes[2];
    sizes[0] = (int) num_dimensions;
    sizes[1] = cache_num_indexes;
    ofs.write((char*) sizes, 2*sizeof(int));
    ofs.write((char*) indexes.data(), indexes.size() * sizeof(int));
}
void MultiIndexSet::readBinary(std::ifstream &ifs){
    int sizes[2];
    ifs.read((char*) sizes, 2*sizeof(int));
    num_dimensions = (size_t) sizes[0];
    cache_num_indexes = sizes[1];
    indexes.resize(num_dimensions * ((size_t) cache_num_indexes));
    ifs.read((char*) indexes.data(), indexes.size() * sizeof(int));
}

void MultiIndexSet::reset(){
    num_dimensions = 0;
    cache_num_indexes = 0;
    std::vector<int> tmp; // create empty vector
    std::swap(tmp, indexes); // swap with empty vector (destroy the empty vector)
}
void MultiIndexSet::copy(const MultiIndexSet &other){
    num_dimensions = other.num_dimensions;
    cache_num_indexes = other.cache_num_indexes;
    indexes = other.indexes;
}

bool MultiIndexSet::empty() const{ return indexes.empty(); }

void MultiIndexSet::setNumDimensions(int new_dimensions){
    reset();
    num_dimensions = (size_t) new_dimensions;
}
int MultiIndexSet::getNumDimensions() const{ return (int) num_dimensions; }
int MultiIndexSet::getNumIndexes() const{ return cache_num_indexes; }

void MultiIndexSet::setIndexes(std::vector<int> &new_indexes){
    indexes = std::move(new_indexes);
    cache_num_indexes = (int) (indexes.size() / num_dimensions);
}

void MultiIndexSet::addSortedInsexes(const std::vector<int> &addition){
    if (indexes.empty()){
        indexes.resize(addition.size());
        std::copy(addition.begin(), addition.end(), indexes.data());
    }else{
        std::vector<int> old_indexes;
        std::swap(old_indexes, indexes);
        std::vector<std::vector<int>::const_iterator> merge_map;
        merge_map.reserve(addition.size() + old_indexes.size());

        SetManipulations::push_merge_map<int, int>(old_indexes, addition,
                                                   [&](std::vector<int>::const_iterator &ia) -> void{ std::advance(ia, num_dimensions); },
                                                   [&](std::vector<int>::const_iterator &ib) -> void{ std::advance(ib, num_dimensions); },
                                                   [&](std::vector<int>::const_iterator ia, std::vector<int>::const_iterator ib) ->
                                                   TypeIndexRelation{
                                                       for(size_t j=0; j<num_dimensions; j++){
                                                           if (*ia   < *ib)   return type_abeforeb;
                                                           if (*ia++ > *ib++) return type_bbeforea;
                                                       }
                                                       return type_asameb;
                                                   },
                                                   [&](std::vector<int>::const_iterator ib, std::vector<std::vector<int>::const_iterator> &mmap) ->
                                                   void { mmap.push_back(ib); },
                                                   merge_map);

        indexes.resize(merge_map.size() * num_dimensions); // merge map will reference only indexes in both sets
        auto iindexes = indexes.begin();
        for(auto &i : merge_map){
            std::copy_n(i, num_dimensions, iindexes);
            std::advance(iindexes, num_dimensions);
        }
    }
    cache_num_indexes = (int) (indexes.size() / num_dimensions);
}
void MultiIndexSet::addUnsortedInsexes(const std::vector<int> &addition){
    size_t num = addition.size() / num_dimensions;
    std::vector<std::vector<int>::const_iterator> index_refs(num);
    auto iadd = addition.begin();
    for(auto &i : index_refs){
        i = iadd;
        std::advance(iadd, num_dimensions);
    }
    std::sort(index_refs.begin(), index_refs.end(),
              [&](std::vector<int>::const_iterator ia, std::vector<int>::const_iterator ib) ->
              bool{
                    for(size_t j=0; j<num_dimensions; j++){
                        if (*ia   < *ib)   return true;
                        if (*ia++ > *ib++) return false;
                    }
                    return false;
            });
    auto unique_end = std::unique(index_refs.begin(), index_refs.end(),
                                  [&](std::vector<int>::const_iterator ia, std::vector<int>::const_iterator ib) ->
                                  bool{
                                        for(size_t j=0; j<num_dimensions; j++) if (*ia++ != *ib++) return false;
                                        return true;
                                });
    index_refs.resize(std::distance(index_refs.begin(), unique_end));
    if (indexes.empty()){
        indexes.resize(index_refs.size() * num_dimensions);
        auto iindexes = indexes.begin();
        for(auto &i : index_refs){
            std::copy_n(i, num_dimensions, iindexes);
            std::advance(iindexes, num_dimensions);
        }
    }else{
        std::vector<int> old_indexes;
        std::swap(old_indexes, indexes);
        std::vector<std::vector<int>::const_iterator> merge_map;
        merge_map.reserve(addition.size() + old_indexes.size());

        SetManipulations::push_merge_map<int, std::vector<int>::const_iterator>(old_indexes, index_refs,
                [&](std::vector<int>::const_iterator &ia) -> void{ std::advance(ia, num_dimensions); },
                [&](std::vector<std::vector<int>::const_iterator>::const_iterator &ib) -> void{ ib++; },
                [&](std::vector<int>::const_iterator ia, std::vector<std::vector<int>::const_iterator>::const_iterator bin) ->
                TypeIndexRelation{
                    std::vector<int>::const_iterator ib = *bin;
                    for(size_t j=0; j<num_dimensions; j++){
                        if (*ia   < *ib)   return type_abeforeb;
                        if (*ia++ > *ib++) return type_bbeforea;
                    }
                    return type_asameb;
                },
                [&](std::vector<std::vector<int>::const_iterator>::const_iterator ib, std::vector<std::vector<int>::const_iterator> &mmap) ->
                void { mmap.push_back(*ib); },
                merge_map);

        indexes.resize(merge_map.size() * num_dimensions); // merge map will reference only indexes in both sets
        auto iindexes = indexes.begin();
        for(auto &i : merge_map){
            std::copy_n(i, num_dimensions, iindexes);
            std::advance(iindexes, num_dimensions);
        }
    }
    cache_num_indexes = (int) (indexes.size() / num_dimensions);
}

const std::vector<int>* MultiIndexSet::getVector() const{ return &indexes; }

int MultiIndexSet::getSlot(const int *p) const{
    size_t sstart = 0, send = (size_t)(cache_num_indexes - 1);
    size_t current = (sstart + send) / 2;
    while (sstart <= send){
        TypeIndexRelation t = [&](const int *a, const int *b) ->
            TypeIndexRelation{
                for(size_t j=0; j<num_dimensions; j++){
                    if (a[j] < b[j]) return type_abeforeb;
                    if (a[j] > b[j]) return type_bbeforea;
                }
                return type_asameb;
            }(&(indexes[current * num_dimensions]), p);
        if (t == type_abeforeb){
            sstart = current+1;
        }else if (t == type_bbeforea){
            send = current-1;
        }else{
            return (int) current;
        };
        current = (sstart + send) / 2;
    }
    return -1;
}

void MultiIndexSet::diffSets(const MultiIndexSet &substract, MultiIndexSet &result){
    result.reset();
    result.setNumDimensions((int) num_dimensions);

    std::vector<std::vector<int>::iterator> kept_indexes;

    auto ithis = indexes.begin();
    auto endthis = indexes.end();
    auto iother = substract.getVector()->begin();
    auto endother = substract.getVector()->end();

    while(ithis != endthis){
        if (iother == endother){
            kept_indexes.push_back(ithis);
            std::advance(ithis, num_dimensions);
        }else{
            TypeIndexRelation t = [&](std::vector<int>::iterator ia, std::vector<int>::const_iterator ib) ->
                                        TypeIndexRelation{
                                            for(size_t j=0; j<num_dimensions; j++){
                                                if (*ia   < *ib)   return type_abeforeb;
                                                if (*ia++ > *ib++) return type_bbeforea;
                                            }
                                            return type_asameb;
                                        }(ithis, iother);
            if (t == type_abeforeb){
                kept_indexes.push_back(ithis);
                std::advance(ithis, num_dimensions);
            }else{
                std::advance(iother, num_dimensions);
                if (t == type_asameb) std::advance(ithis, num_dimensions);
            }
        }
    }

    if (kept_indexes.size() > 0){
        std::vector<int> new_indexes(num_dimensions * kept_indexes.size());
        auto inew = new_indexes.begin();
        for(auto i : kept_indexes){
            std::copy_n(i, num_dimensions, inew);
            std::advance(inew, num_dimensions);
        }
        result.setIndexes(new_indexes);
    }
}

StorageSet::StorageSet() : num_outputs(0), num_values(0){}
StorageSet::StorageSet(int cnum_outputs, int cnum_values) : num_outputs((size_t) cnum_outputs), num_values((size_t) cnum_values){}
StorageSet::StorageSet(const StorageSet *storage) : values(0){
    num_outputs = storage->num_outputs;
    num_values = storage->num_values;
    values = storage->values; // copy assignment
}
StorageSet::~StorageSet(){}

void StorageSet::write(std::ofstream &ofs) const{
    ofs << num_outputs << " " << num_values;
    if (values.size() != 0){
        ofs << " 1";
        ofs << std::scientific; ofs.precision(17);
        for(auto v : values) ofs << " " << v;
    }else{
        ofs << " 0";
    }
    ofs << std::endl;
}
void StorageSet::read(std::ifstream &ifs){
    reset(); // empty values if the file doesn't contain vals
    int has_vals;
    ifs >> num_outputs >> num_values >> has_vals;
    if (has_vals == 1){
        values.resize(num_outputs * num_values);
        for(auto &v : values) ifs >> v;
    }
}
void StorageSet::writeBinary(std::ofstream &ofs) const{
    int num_out_vals[2];
    num_out_vals[0] = (int) num_outputs;
    num_out_vals[1] = (int) num_values;
    ofs.write((char*) num_out_vals, 2*sizeof(int));
    if (values.size() != 0){
        char flag = 'y'; ofs.write((char*) &flag, sizeof(char));
        ofs.write((char*) values.data(), values.size() * sizeof(double));
    }else{
        char flag = 'n'; ofs.write((char*) &flag, sizeof(char));
    }
}
void StorageSet::readBinary(std::ifstream &ifs){
    int num_out_vals[2];
    ifs.read((char*) num_out_vals, 2*sizeof(int));
    num_outputs = (size_t) num_out_vals[0];
    num_values = (size_t) num_out_vals[1];
    char flag; ifs.read((char*) &flag, sizeof(char));
    if (flag == 'y'){
        values.resize(num_outputs * num_values);
        ifs.read((char*) values.data(), values.size() * sizeof(double));
    }else{
        values.resize(0); // empty values if the file doesn't contain vals
    }
}

void StorageSet::reset(){
    num_outputs = 0;
    num_values = 0;
    std::vector<double> temp;
    std::swap(values, temp); // ensures that values is empty in size() and capacity()
}
void StorageSet::copy(const StorageSet &other){
    num_outputs = other.num_outputs;
    num_values = other.num_values;
    values = other.values;
}
void StorageSet::resize(int cnum_outputs, int cnum_values){
    reset();
    num_outputs = cnum_outputs;
    num_values = cnum_values;
}

int StorageSet::getNumOutputs() const{ return (int) num_outputs; }
const double* StorageSet::getValues(int i) const{ return &(values[i*num_outputs]); }
double* StorageSet::getValues(int i){ return &(values[i*num_outputs]); }
std::vector<double>* StorageSet::aliasValues(){ return &values; }

void StorageSet::setValues(const double vals[]){
    values.resize(num_outputs * num_values);
    std::copy_n(vals, num_values * num_outputs, values.data());
}
void StorageSet::setValues(std::vector<double> &vals){
    num_values = vals.size() / num_outputs;
    values = std::move(vals); // move assignment
}

void StorageSet::addValues(const IndexSet *old_set, const IndexSet *new_set, const double new_vals[]){
    std::vector<double> old_vals = std::move(values);
    size_t old_num_values = num_values;
    size_t new_num_values = new_set->getNumIndexes();
    int num_dimensions = old_set->getNumDimensions();

    TypeIndexRelation relation;
    values.resize(num_outputs * (old_num_values + new_num_values));
    size_t offsetOld = 0, offsetNew = 0;
    num_values = 0;
    while( (offsetNew < new_num_values) || (offsetOld < old_num_values) ){
        if (offsetNew >= new_num_values){ // new is a
            relation = type_bbeforea;
        }else if (offsetOld >= old_num_values){ // old is b
            relation = type_abeforeb;
        }else{
            relation = compareIndexes(num_dimensions, new_set->getIndex((int) offsetNew), old_set->getIndex((int) offsetOld));
        }
        if (relation == type_abeforeb){
            std::copy(&(new_vals[num_outputs * offsetNew]), &(new_vals[num_outputs * offsetNew]) + num_outputs, &(values[num_outputs * num_values++]));
            offsetNew++;
        }else{
            std::copy(&(old_vals[num_outputs * offsetOld]), &(old_vals[num_outputs * offsetOld]) + num_outputs, &(values[num_outputs * num_values++]));
            offsetOld++;
        }
    }
}
void StorageSet::addValues(const MultiIndexSet &old_set, const MultiIndexSet &new_set, const double new_vals[]){
    int num_old = old_set.getNumIndexes();
    int num_new = new_set.getNumIndexes();
    int num_dimensions = old_set.getNumDimensions();

    num_values += (size_t) num_new;
    std::vector<double> combined_values(num_values * num_outputs);

    int iold = 0, inew = 0;
    size_t off_vals = 0;
    auto ivals = values.begin();
    auto icombined = combined_values.begin();

    for(size_t i=0; i<num_values; i++){
        TypeIndexRelation relation;
        if (iold >= num_old){
            relation = type_abeforeb;
        }else if (inew >= num_new){
            relation = type_bbeforea;
        }else{
            relation = compareIndexes(num_dimensions, new_set.getIndex(inew), old_set.getIndex(iold));
        }
        if (relation == type_abeforeb){
            std::copy_n(&(new_vals[off_vals]), num_outputs, icombined);
            inew++;
            off_vals += num_outputs;
        }else{
            std::copy_n(ivals, num_outputs, icombined);
            iold++;
            std::advance(ivals, num_outputs);
        }
        std::advance(icombined, num_outputs);
    }
    std::swap(values, combined_values);
}

TypeIndexRelation StorageSet::compareIndexes(int num_dimensions, const int a[], const int b[]) const{
    for(int i=0; i<num_dimensions; i++){
        if (a[i] < b[i]) return type_abeforeb;
        if (a[i] > b[i]) return type_bbeforea;
    }
    return type_asameb;
}

}
#endif
