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

#ifndef __TASMANIAN_SPARSE_GRID_INDEX_SETS_HPP
#define __TASMANIAN_SPARSE_GRID_INDEX_SETS_HPP

#include "tsgEnumerates.hpp"
#include <vector>
#include <functional>
#include <algorithm>

namespace TasGrid{

class UnsortedIndexSet{ // assumes the elements would be added in a sorted order
public:
    UnsortedIndexSet(int cnum_dimensions, int cnum_slots);
    ~UnsortedIndexSet();

    int getNumDimensions() const;

    void addIndex(const int p[]);

    void getIndexesSorted(std::vector<int> &sorted) const;
    // returns a vector of num_dimensions X num_indexes of the sorted indexes

protected:
    TypeIndexRelation compareIndexes(const int a[], const int b[]) const;
    void mergeLists(const int listA[], size_t sizeA, const int listB[], size_t sizeB, int destination[]) const;

private:
    size_t num_dimensions;
    size_t num_indexes;
    std::vector<int> index;
};

class GranulatedIndexSet{ // optimized to add one index at a time
public:
    GranulatedIndexSet(int cnum_dimensions, int cnum_slots = 64);
    ~GranulatedIndexSet();

    int getNumDimensions() const;
    int getNumIndexes() const;

    void addIndex(const int p[]);
    inline void addIndex(const std::vector<int> p){ addIndex(p.data()); }
    int getSlot(const int p[]) const;

    const int* getIndex(int j) const;

    void addGranulatedSet(const GranulatedIndexSet *gset);

    const std::vector<int>* getIndexes() const;
    const std::vector<int>* getMap() const;

protected:
    TypeIndexRelation compareIndexes(const int a[], const int b[]) const;

    void mergeMapped(const std::vector<int> newIndex, const std::vector<int> newMap, int sizeNew);

private:
    size_t num_dimensions;
    size_t num_indexes;
    std::vector<int> index;
    std::vector<int> imap;
};

class IndexSet{ // rigid set but optimal in storage size
public:
    IndexSet(int cnum_dimensions);
    IndexSet(const UnsortedIndexSet *uset);
    IndexSet(const GranulatedIndexSet *gset);
    IndexSet(const IndexSet *iset);
    IndexSet(int cnum_dimensions, std::vector<int> &cindex); // move assignment
    ~IndexSet();

    int getNumDimensions() const;
    int getNumIndexes() const;

    void write(std::ofstream &ofs) const;
    void read(std::ifstream &ifs);

    void writeBinary(std::ofstream &ofs) const;
    void readBinary(std::ifstream &ifs);

    int getSlot(const int p[]) const;
    inline int getSlot(const std::vector<int> &p) const{ return getSlot(p.data()); }

    const int* getIndex(int i) const;
    const std::vector<int>* getIndexes() const;

    void addUnsortedSet(const UnsortedIndexSet *uset);
    void addGranulatedSet(const GranulatedIndexSet *gset);
    void addIndexSet(const IndexSet *iset);

    IndexSet* diffSets(const IndexSet *iset) const; // returns this set minus the points in set

protected:
    TypeIndexRelation compareIndexes(const int a[], const int b[]) const;

    void mergeSet(const std::vector<int> newIndex);
    void mergeMapped(const std::vector<int> *newIndex, const std::vector<int> *imap);

    void cacheNumIndexes();

private:
    size_t num_dimensions;
    int cache_num_indexes;
    std::vector<int> index;
};

namespace SetManipulations{
template<typename T, class B>
void push_merge_map(const std::vector<T> &a, const std::vector<B> &b,
                    std::function<void(typename std::vector<T>::const_iterator &ia)> iadvance,
                    std::function<void(typename std::vector<B>::const_iterator &ib)> ibdvance,
                    std::function<TypeIndexRelation(typename std::vector<T>::const_iterator ia, typename std::vector<B>::const_iterator ib)> compare,
                    std::function<void(typename std::vector<B>::const_iterator ib, std::vector<typename std::vector<T>::const_iterator> &merge_map)> pushB,
                    std::vector<typename std::vector<T>::const_iterator> &merge_map){
//! \internal
//! \brief merge two index-sets, creates a vector of iterators (pointers) to the multi-indexes in the merged (set union) index-set
//! \ingroup TasmanianMultiIndexSet
//!
//! Takes two multi-index sets described by **a** and **b**, generates a vector **merge_map** that describes the union set.
//! This operation is used when two multi-index sets are merged, which can happen which can happen with set **b** sorted or unsorted.
//! In the sorted case, the type for **a** and **b** is the same and **iadvance()** is the same as **ibadvance()**.
//! In the unsorted case, **b** will sorted references to the multi-indexes, i.e., B = std::vector<T>::const_iterator
//! * T is a type of multi-index entries, e.g., int, unsigned int, long long;
//! * **a** is a sorted set of multi-indexes of type T, each multi-index has size implicitly defined by **iadvance()**; entries associated to the same multi-index are adjacent
//! * **b** gives a description of the second set, in the simple case **b** has the same type as **a** (e.g., B == T); alternatively **b** can be a list of sorted references to the multi-indexes, e.g., when the **b** set was first sorted and then it is being merged;
//! * **iadvance(ia)** takes an iterator pointing to multi-index and moves it to the next multi-index, e.g., std::advance(ia, num_dimensions);
//! * **ibdvance(ib)** takes an iterator pointing to the multi-index described in **b** and moves it to the next multi-index, e.g., std::advance(ib, num_dimensions) or ib++;
//! * **compare(ia, ib)** returns the order relation between the multi-indexes, *type_abeforeb*, *type_asameb* and *type_bbeforea* (note that **ia** and **ib** always point to **a** and **b**);
//! * **pushB()** add the multi-index described by **b** to the **merge_map**, i.e., pushes **ib** or pushes deference of **ib** (note that **ib** always points to an element in **b**);
//! * on exit, **merge_map** holds the references to the multi-indexes of the union set
    auto ia = a.begin(), ib = b.begin();
    auto aend = a.end(), bend = b.end();
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
            pushB(ib, merge_map);
            ibdvance(ib);
        }else{
            merge_map.push_back(ia);
            iadvance(ia);
            if (relation == type_asameb) ibdvance(ib);
        }
    }
}
}

template<typename T>
class Data2D{
// this class is a work around of using indexing of type [i * stride + j] where i, j, and stride are int and can overflow
// the class internally uses size_t and can pass information back and forth between API calls
// the idea is to apply this to data stored as either std::vector<int/double> or double[]/int[]
public:
    Data2D() : stride(0), num_strips(0), data(0){}
    ~Data2D(){}

    void resize(int new_stride, int new_num_strips){
        stride = (size_t) new_stride;
        num_strips = (size_t) new_num_strips;
        vec.resize(stride * num_strips);
        data = vec.data();
        cdata = vec.data();
    }
    void resize(int new_stride, int new_num_strips, T val){
        stride = (size_t) new_stride;
        num_strips = (size_t) new_num_strips;
        vec.resize(stride * num_strips, val);
        data = vec.data();
        cdata = vec.data();
    }
    void load(int new_stride, int new_num_strips, T* new_data){
        stride = (size_t) new_stride;
        num_strips = (size_t) new_num_strips;
        data = new_data;
        cdata = data;
        vec.resize(0);
    }
    void cload(int new_stride, int new_num_strips, const T* new_data){
        stride = (size_t) new_stride;
        num_strips = (size_t) new_num_strips;
        cdata = new_data;
        data = 0;
        vec.resize(0);
    }

    T* getStrip(int i){ return &(data[i*stride]); }
    const T* getCStrip(int i) const{ return &(cdata[i*stride]); }
    int getStride() const{ return (int) stride; }
    int getNumStrips() const{ return (int) num_strips; }
    size_t getTotalEntries() const{ return stride * num_strips; }
    std::vector<T>* getVector(){ return &vec; }
    const std::vector<T>* getVector() const{ return &vec; }
    void clear(){
        stride = 0;
        num_strips = 0;
        vec.clear();
        vec.shrink_to_fit();
    }

    void appendStrip(const std::vector<T> &x){
        vec.insert(vec.end(), x.begin(), x.end());
        num_strips++;
    }

private:
    size_t stride, num_strips;
    T* data;
    const T* cdata;
    std::vector<T> vec;
};

class MultiIndexSet{
public:
    MultiIndexSet();
    MultiIndexSet(int cnum_dimensions);
    ~MultiIndexSet();

    void reset(); // empty and free memory
    void move(MultiIndexSet &other); // this becomes other, other becomes empty
    bool empty() const;

    void setNumDimensions(int new_dimensions);
    int getNumDimensions() const;
    int getNumIndexes() const;

    void setIndexes(std::vector<int> &new_indexes); // move assignment
    void addSortedInsexes(const std::vector<int> &addition);   // merge/copy assignment
    void addUnsortedInsexes(const std::vector<int> &addition); // sort/merge/copy assignment

    const std::vector<int>* getVector() const;
    int getSlot(const int *p) const;
    inline int getSlot(const std::vector<int> &p) const{ return getSlot(p.data()); }

private:
    size_t num_dimensions;
    int cache_num_indexes;
    std::vector<int> indexes;
};

class StorageSet{ // stores the values of the function
public:
    StorageSet(int cnum_outputs, int cnum_values);
    StorageSet(const StorageSet *storage);
    ~StorageSet();

    void write(std::ofstream &ofs) const;
    void read(std::ifstream &ifs);
    void writeBinary(std::ofstream &ofs) const;
    void readBinary(std::ifstream &ifs);

    const double* getValues(int i) const;
    double* getValues(int i);
    std::vector<double>* aliasValues(); // alternative to setValues()
    int getNumOutputs() const;

    void setValues(const double vals[]);
    void setValues(std::vector<double> &vals);
    void addValues(const IndexSet *old_set, const IndexSet *new_set, const double new_vals[]);

protected:
    TypeIndexRelation compareIndexes(int num_dimensions, const int a[], const int b[]) const;

private:
    size_t num_outputs, num_values; // kept as size_t to avoid conversions in products, but each one is small individually
    std::vector<double> values;
};

}

#endif
