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

template<bool iomode>
void MultiIndexSet::write(std::ostream &os) const{
    if (cache_num_indexes > 0){
        IO::writeNumbers<iomode, IO::pad_rspace>(os, static_cast<int>(num_dimensions), cache_num_indexes);
        IO::writeVector<iomode, IO::pad_line>(indexes, os);
    }else{
        IO::writeNumbers<iomode, IO::pad_line>(os, static_cast<int>(num_dimensions), cache_num_indexes);
    }
}

template void MultiIndexSet::write<mode_ascii>(std::ostream &) const; // instantiate for faster build
template void MultiIndexSet::write<mode_binary>(std::ostream &) const;

void MultiIndexSet::addSortedIndexes(const std::vector<int> &addition){
    if (indexes.empty()){
        indexes = addition;
    }else{
        std::vector<int> old_indexes = std::move(indexes);
        indexes = std::vector<int>(); // reset the indexes

        std::vector<std::vector<int>::const_iterator> merge_map;
        merge_map.reserve(addition.size() + old_indexes.size());
        auto next = std::back_inserter(merge_map);

        auto compare = [&](std::vector<int>::const_iterator ia, std::vector<int>::const_iterator ib) ->
                            TypeIndexRelation{
                                for(size_t j=0; j<num_dimensions; j++){
                                    if (*ia   < *ib)   return type_abeforeb;
                                    if (*ia++ > *ib++) return type_bbeforea;
                                }
                                return type_asameb;
                            };

        // merge with three way compare
        auto ia = old_indexes.begin();
        auto ib = addition.begin();
        auto aend = old_indexes.end();
        auto bend = addition.end();
        while((ia != aend) || (ib != bend)){
            TypeIndexRelation relation;
            if (ib == bend){
                relation = type_abeforeb;
            }else if (ia == aend){
                relation = type_bbeforea;
            }else{
                relation = compare(ia, ib);
            }

            if (relation == type_bbeforea){
                next = ib;
                std::advance(ib, num_dimensions);
            }else{
                next = ia;
                std::advance(ia, num_dimensions);
                if (relation == type_asameb) std::advance(ib, num_dimensions);
            }
        }

        indexes.resize(merge_map.size() * num_dimensions); // merge map will reference only indexes in both sets
        auto iindexes = indexes.begin();
        for(auto &i : merge_map){
            std::copy_n(i, num_dimensions, iindexes);
            std::advance(iindexes, num_dimensions);
        }
    }
    cache_num_indexes = (int) (indexes.size() / num_dimensions);
}

MultiIndexSet::MultiIndexSet(Data2D<int> const &data) : num_dimensions((size_t) data.getStride()), cache_num_indexes(0){
    size_t num = (size_t) data.getNumStrips();
    if (num == 0) return; // nothing to do

    std::vector<std::vector<int>::const_iterator> index_refs(num);
    auto iadd = data.begin();
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

    indexes.resize(index_refs.size() * num_dimensions);
    auto iindexes = indexes.begin();
    for(auto &i : index_refs){
        std::copy_n(i, num_dimensions, iindexes);
        std::advance(iindexes, num_dimensions);
    }

    cache_num_indexes = (int) (indexes.size() / num_dimensions);
}

int MultiIndexSet::getSlot(const int *p) const{
    int sstart = 0, send = cache_num_indexes - 1;
    int current = (sstart + send) / 2;
    while (sstart <= send){
        TypeIndexRelation t = [&](const int *a, const int *b) ->
            TypeIndexRelation{
                for(size_t j=0; j<num_dimensions; j++){
                    if (a[j] < b[j]) return type_abeforeb;
                    if (a[j] > b[j]) return type_bbeforea;
                }
                return type_asameb;
            }(&(indexes[(size_t) current * num_dimensions]), p);
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

MultiIndexSet MultiIndexSet::operator -(const MultiIndexSet &substract) const{
    std::vector<std::vector<int>::const_iterator> kept_indexes;

    auto ithis = indexes.begin();
    auto endthis = indexes.end();
    auto iother = substract.begin();
    auto endother = substract.end();

    while(ithis != endthis){
        if (iother == endother){
            kept_indexes.push_back(ithis);
            std::advance(ithis, num_dimensions);
        }else{
            TypeIndexRelation t = [&](std::vector<int>::const_iterator ia, std::vector<int>::const_iterator ib) ->
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
        return MultiIndexSet(num_dimensions, std::move(new_indexes));
    }

    return MultiIndexSet();
}

void MultiIndexSet::removeIndex(const std::vector<int> &p){
    int slot = getSlot(p);
    if (slot > -1){
        indexes.erase(indexes.begin() + static_cast<size_t>(slot) * num_dimensions, indexes.begin() + static_cast<size_t>(slot) * num_dimensions + num_dimensions);
        cache_num_indexes--;
    }
}

template<bool iomode>
void StorageSet::write(std::ostream &os) const{
    IO::writeNumbers<iomode, IO::pad_rspace>(os, static_cast<int>(num_outputs), static_cast<int>(num_values));
    IO::writeFlag<iomode, IO::pad_auto>((values.size() != 0), os);
    if (values.size() != 0)
        IO::writeVector<iomode, IO::pad_line>(values, os);
}

template void StorageSet::write<mode_ascii>(std::ostream &) const;
template void StorageSet::write<mode_binary>(std::ostream &) const;

void StorageSet::addValues(const MultiIndexSet &old_set, const MultiIndexSet &new_set, const double new_vals[]){
    int num_old = old_set.getNumIndexes();
    int num_new = new_set.getNumIndexes();
    size_t num_dimensions = old_set.getNumDimensions();

    num_values += (size_t) num_new;
    std::vector<double> combined_values(num_values * num_outputs);

    int iold = 0, inew = 0;
    size_t off_vals = 0;
    auto ivals = values.begin();
    auto icombined = combined_values.begin();

    auto compareIndexes = [&](int const a[], int const b[])->
            TypeIndexRelation{
                for(size_t j=0; j<num_dimensions; j++){
                    if (a[j] < b[j]) return type_abeforeb;
                    if (a[j] > b[j]) return type_bbeforea;
                }
                return type_asameb;
            };

    for(size_t i=0; i<num_values; i++){
        TypeIndexRelation relation;
        if (iold >= num_old){
            relation = type_abeforeb;
        }else if (inew >= num_new){
            relation = type_bbeforea;
        }else{
            relation = compareIndexes(new_set.getIndex(inew), old_set.getIndex(iold));
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

}
#endif
