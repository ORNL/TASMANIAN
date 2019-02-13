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

#include "tsgIOHelpers.hpp"
#include "tsgEnumerates.hpp"
#include <vector>
#include <functional>
#include <algorithm>

//! \internal
//! \file tsgIndexSets.hpp
//! \brief Data structures for storing multi-indexes and floating point values.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianSets
//!
//! Three classes for storage and easy access to multi-indexes and floating point values.
//! The majority of internal data-structures associated with the five types of grids
//! can be represented as two-dimensional data (or array/list of tuples), the three
//! classes here correspond to sorted or unsorted, integer or floating point data.

//! \internal
//! \defgroup TasmanianSets Set data structures.
//!
//! \par Fundamental data structures
//! Three classes for storage and easy access to multi-indexes and floating point values.
//! The majority of internal data-structures associated with the five types of grids
//! can be represented as two-dimensional data (or array/list of tuples), the three
//! classes here correspond to sorted or unsorted, integer or floating point data.

namespace TasGrid{

namespace SetManipulations{

//! \internal
//! \brief Merge two multi-index sets, creates a vector of iterators (pointers) to the multi-indexes in the merged index-set (set union)
//! \ingroup TasmanianSets

//! Takes two multi-index sets described by **a** and **b**, generates a vector **merge_map** that describes the union set.
//! This operation is used when two multi-index sets are merged, which can happen in two cases depending whether the set corresponding to **b** is sorted or not.
//! In the sorted case, the type for **a** and **b** is the same and **iadvance()** is the same as **ibadvance()**.
//! In the unsorted case, the variable **b** will hold sorted references to the multi-indexes, i.e., `B == std::vector<T>:: const_iterator`.
//! The second case saves one copy operation, i.e., the merge can be performed with the sorted references as opposed to the full sorted set.
//! * T is a type of multi-index entries, e.g., int, unsigned int, long long
//! * **a** is a sorted set of multi-indexes of type T, each multi-index has size implicitly defined by **iadvance()**; entries associated to the same multi-index are adjacent
//! * **b** gives a description of the second set, in the simple case **b** has the same type as **a** (e.g., B == T); alternatively **b** can be a list of sorted references to the multi-indexes, e.g., when the second set was first sorted externally to this function and now the sorted references are used for the merge
//! * **iadvance(ia)** takes an iterator pointing to multi-index and moves it to the next multi-index, e.g., std::advance(ia, num_dimensions)
//! * **ibdvance(ib)** takes an iterator pointing to the multi-index described in **b** and moves it to the next multi-index, e.g., std::advance(ib, num_dimensions) or ib++
//! * **compare(ia, ib)** returns the order relation between the multi-indexes, *type_abeforeb*, *type_asameb* and *type_bbeforea* (note that **ia** and **ib** always point to **a** and **b**)
//! * **pushB()** adds the multi-index described by **b** to the **merge_map**, i.e., pushes **ib** or pushes de-reference of **ib** (note that **ib** is an iterator to **b**)
//! * on exit, **merge_map** holds the references to the multi-indexes of the union set, note that the references are added with `std::vector::push_back` command
template<typename T, class B>
void push_merge_map(const std::vector<T> &a, const std::vector<B> &b,
                    std::function<void(typename std::vector<T>::const_iterator &ia)> iadvance,
                    std::function<void(typename std::vector<B>::const_iterator &ib)> ibdvance,
                    std::function<TypeIndexRelation(typename std::vector<T>::const_iterator ia, typename std::vector<B>::const_iterator ib)> compare,
                    std::function<void(typename std::vector<B>::const_iterator ib, std::vector<typename std::vector<T>::const_iterator> &merge_map)> pushB,
                    std::vector<typename std::vector<T>::const_iterator> &merge_map){
    auto ia = a.begin();
    auto ib = b.begin();
    auto aend = a.end();
    auto bend = b.end();
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

//! \brief Generic 2D data structure divided into contiguous strips of fixed length (similar to a matrix)
//! \ingroup TasmanianSets

//! Many of the internal data-structures relevant to sparse grids can be represented as two dimensional arrays.
//! The data is divided into **strips** of equal **stride**, e.g., number of multi-indexes and number of dimensions.
//! Internally the class uses either `std::vector` or const/non-const array pointers; the template parameter **T** could be anything
//! that can fit in a vector, but in practice the class is used primarily with `double` and `int`.
//! This panacea template handles multiple situations:
//! * allocate memory for the 2D data structure and access that memory with simple indexing avoiding the direct use of `[i * stride +j]`
//! * dynamically generate a data structure by appending strips to the back of the existing ones
//! * wrap around an existing data structure that can be either const or non-const array, this class preserves const correctness
template<typename T>
class Data2D{
public:
    //! \brief Default constructor makes an empty data-structure
    Data2D() : stride(0), num_strips(0), data(0){}
    //! \brief Default destructor
    ~Data2D(){}

    //! \brief Clear any existing data and allocate a new data-structure with given **stride** and number of **strips**
    void resize(int new_stride, int new_num_strips){
        stride = (size_t) new_stride;
        num_strips = (size_t) new_num_strips;
        vec.resize(stride * num_strips);
        data = vec.data();
        cdata = vec.data();
    }
    //! \brief Clear any existing data, allocate a new data-structure with given **stride** and number of **strips** and initializes the data to **val** (using `std::vector::resize()`)
    void resize(int new_stride, int new_num_strips, T val){
        stride = (size_t) new_stride;
        num_strips = (size_t) new_num_strips;
        vec.resize(stride * num_strips, val);
        data = vec.data();
        cdata = vec.data();
    }
    //! \brief Wrap around an existing array, the size of the array must be at least **new_stride** times **new_num_strips**
    void load(int new_stride, int new_num_strips, T* new_data){
        stride = (size_t) new_stride;
        num_strips = (size_t) new_num_strips;
        data = new_data;
        cdata = data;
        vec.resize(0);
    }
    //! \brief Wrap around an existing const array, the size of the array must be at least **new_stride** times **new_num_strips**
    void cload(int new_stride, int new_num_strips, const T* new_data){
        stride = (size_t) new_stride;
        num_strips = (size_t) new_num_strips;
        cdata = new_data;
        data = 0;
        vec.resize(0);
    }

    //! \brief Returns a reference to the **i**-th strip, cannot be called after **cload()** since the reference is non-const
    T* getStrip(int i){ return &(data[i*stride]); }
    //! \brief Returns a const reference to the **i**-th strip, can be called even after **cload()**
    const T* getCStrip(int i) const{ return &(cdata[i*stride]); }
    //! \brief Returns the stride
    int getStride() const{ return (int) stride; }
    //! \brief Returns the number of strips
    int getNumStrips() const{ return (int) num_strips; }
    //! \brief Returns **getStride()** times **getNumStrips()** but using `size_t` arithmetic avoiding a potential `int` overflow
    size_t getTotalEntries() const{ return stride * num_strips; }
    //! \brief Returns a reference to the `std::vector` that holds the internal data, can be used only after a call to **resize()**
    std::vector<T>* getVector(){ return &vec; }
    //! \brief Returns a const reference to the `std::vector` that holds the internal data, can be used only after a call to **resize()**
    const std::vector<T>* getVector() const{ return &vec; }
    //! \brief Clear all used data (redundant, remove)
    void clear(){
        stride = 0;
        num_strips = 0;
        vec.clear();
        vec.shrink_to_fit();
    }

    //! \brief Uses `std::vector::insert` to append a strip **x** to the existing data
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


//! \brief Class that stores multi-indexes in sorted (lexicographical) order
//! \ingroup TasmanianSets

//! At the core of each sparse grid, there are multiple multi-index sets.
//! The organization of the data is similar to the **Data2D<T>** class, but at any time the indexes
//! are stored in a lexicographical order. The main functionality provided here is:
//! * fast *O(log(n))* search utilizing the lexicographical order
//! * synchronization between multi-indexes and values (i.e., model outputs)
//! * adding or removing indexes while preserving the order
//! * basic file I/O
class MultiIndexSet{
public:
    //! \brief Default constructor, makes an empty set
    MultiIndexSet();
    //! \brief Constructor, makes an empty set and assigns the dimension of the future indexes
    MultiIndexSet(int cnum_dimensions);
    //! \brief Default destructor
    ~MultiIndexSet();

    //! \brief Write the set to ASCII or binary stream, use with std::ofstream and std::ifstream.

    //! The format consists of two `int` values corresponding to the number of dimensions and number of indexes,
    //! followed by all the entries of the array on a single line separated by a space, or dump of a single write command.
    template<bool useAscii> void write(std::ostream &os) const;

    //! \brief Read the from the stream, must know whether to use ASCII or binary format.

    //! Uses the same format as \b write<bool>
    template<bool useAscii> void read(std::istream &os);

    //! \brief Returns **true** if there are no multi-indexes in the set, **false** otherwise
    inline bool empty() const{ return indexes.empty(); }

    //! \brief Clears any currently loaded multi-indexes and set the dimensions to **new_dimensions**
    void setNumDimensions(int new_dimensions);
    //! \brief Returns the number of dimensions
    inline int getNumDimensions() const{ return (int) num_dimensions; }
    //! \brief Returns the number of indexes
    inline int getNumIndexes() const{ return cache_num_indexes; }

    //! \brief Overwrites the existing set of indexes using `std::move()`, note that **new_indexes** must be sorted before this call
    void setIndexes(std::vector<int> &new_indexes);
    //! \brief Add more indexes to the set, note that **addition** must be sorted before this call
    void addSortedInsexes(const std::vector<int> &addition);
    //! \brief Add more indexes to the set, **addition** is assumed unsorted (will be sorted first)
    void addUnsortedInsexes(const std::vector<int> &addition);
    //! \brief Add indexes from another **MultiIndexSet**, makes a call to **addSortedInsexes()**
    inline void addMultiIndexSet(const MultiIndexSet &addition){ addSortedInsexes(*addition.getVector()); }
    //! \brief Add indexes from a general **Data<>** structure, makes a call to **addUnsortedInsexes()**
    inline void addData2D(const Data2D<int> &addition){ addUnsortedInsexes(*addition.getVector()); }

    //! \brief Returns a const reference to the internal data
    inline const std::vector<int>* getVector() const{ return &indexes; }
    //! \brief Returns a reference to the internal data, must not modify the lexicographical order or the size of the vector
    inline std::vector<int>* getVector(){ return &indexes; } // used for remapping during tensor generic points

    //! \brief Returns the slot containing index **p**, returns `-1` if not found
    int getSlot(const int *p) const;
    //! \brief Returns the slot containing index **p**, returns `-1` if not found
    inline int getSlot(const std::vector<int> &p) const{ return getSlot(p.data()); }
    //! \brief Returns **true** if **p** is missing from the set, **false** otherwise
    inline bool missing(const std::vector<int> &p) const{ return (getSlot(p.data()) == -1); }

    //! \brief Returns the **i**-th index of the set, useful to loop over all indexes or to cross reference with values
    inline const int *getIndex(int i) const{ return &(indexes[((size_t) i) * num_dimensions]); }

    //! \brief A new ordered set is created in **result**, which holds the indexes from this set that are not present in **substract**
    //!
    //! The implementation uses an algorithm similar to merge with complexity linear in the number of multi-indexes of the two sets,
    //! i.e., does not use **missing()** which would add a logarithmic factor.
    void diffSets(const MultiIndexSet &substract, MultiIndexSet &result);

    //! \brief Removes \b p from the set (if exists).
    void removeIndex(const std::vector<int> &p);

    //! \brief Returns the maximum single index in the set.
    int getMaxIndex() const{ return (empty()) ? 0 : *std::max_element(indexes.begin(), indexes.end()); }

private:
    size_t num_dimensions;
    int cache_num_indexes;
    std::vector<int> indexes;
};

//! \brief Class that stores values, i.e., model outputs, the order of the values is in sync with the order of some **MultiIndexSet**
//! \ingroup TasmanianSets

//! The StorageSet stores the floating-point values (model outputs) in a contiguous order, suitable for interfacing with BLAS/cuBlas.
//! The values associated with one model simulation are adjacent to each other in a stride of size **num_outputs**.
//! The order of the values is kept in sync with a **MultiIndexSet** (usually the set **points** in the Grid classes),
//! Synchronization is achieved by merging values before merging multi-indexes, when
//! the **addValues()** is called with the old and new multi-index sets and the new values.
class StorageSet{
public:
    //! \brief Default constructor, makes an empty set
    StorageSet();
    //! \brief Default destructor
    ~StorageSet();

    //! \brief Write the set to ASCII or binary stream, use with std::ofstream and std::ifstream.

    //! The format consists of two `int` values corresponding to the number of dimensions and number of indexes,
    //! followed by all the entries of the array on a single line separated by a space, or dump of a single write command.
    template<bool useAscii> void write(std::ostream &os) const;

    //! \brief Read the from the stream, must know whether to use ASCII or binary format.

    //! Uses the same format as \b write<bool>
    template<bool useAscii> void read(std::istream &os);

    //! \brief Clear the existing values and assigns new dimensions, does not allocate memory for the new values
    void resize(int cnum_outputs, int cnum_values);

    //! \brief Returns const reference to the **i**-th value, the **i** index matches the corresponding **MultiIndexSet**
    const double* getValues(int i) const;
    //! \brief Returns reference to the **i**-th value, the **i** index matches the corresponding **MultiIndexSet**
    double* getValues(int i);
    //! \brief Returns reference to the internal data vector
    std::vector<double>* aliasValues(); // alternative to setValues()
    //! \brief Returns const reference to the internal data vector
    const std::vector<double>* aliasValues() const; // alternative to setValues()

    //! \brief Replace the existing values with a copy of **vals**, the size must be at least **num_outputs** times **num_values**
    void setValues(const double vals[]);
    //! \brief Replace the existing values with **vals** using `std::move()`, the size of **vals** must be **num_outputs** times **num_values**
    void setValues(std::vector<double> &vals);

    //! \brief Add more values to the set, the **old_set** and **new_set** define the order of the current values and the **new_vals**
    //!
    //! Add more values to an existing set of values, where **old_set** and **new_set** indicate the order.
    //! The **old_set** contains the ordered multi-indexes associated with the current values, the **new_set** corresponds to the order
    //! of the **new_vals**. After the merge, the order of the values will match that of the union of **old_set** and **new_set**.
    void addValues(const MultiIndexSet &old_set, const MultiIndexSet &new_set, const double new_vals[]);

protected:
    //! \brief Returns the relation between **a** and **b**, used for the merge
    TypeIndexRelation compareIndexes(int num_dimensions, const int a[], const int b[]) const;

private:
    size_t num_outputs, num_values; // kept as size_t to avoid conversions in products, but each one is small individually
    std::vector<double> values;
};

}

#endif
