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

namespace TasGrid{

class UnsortedIndexSet{ // assumes the elements would be added in a sorted order
public:
    UnsortedIndexSet(int cnum_dimensions, int cnum_slots);
    ~UnsortedIndexSet();

    int getNumDimensions() const;
    int getNumIndexes() const;

    void addIndex(const int p[]);
    const int* getIndex(int i) const;

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
    inline int getSlot(const std::vector<int> p) const{ return getSlot(p.data()); }

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

private:
    size_t num_dimensions;
    std::vector<int> index;
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
