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

namespace TasGrid{

class UnsortedIndexSet{ // assumes the elements would be added in a sorted order
public:
    UnsortedIndexSet(int cnum_dimensions, int cnum_slots);
    ~UnsortedIndexSet();

    int getNumDimensions() const;
    int getNumIndexes() const;

    void addIndex(const int p[]);
    const int* getIndex(int i) const;

    int* getIndexesSorted() const;
    // returns an array of num_dimensions X num_indexes of the sorted indexes

protected:
    TypeIndexRelation compareIndexes(const int a[], const int b[]) const;
    void merge(const int listA[], int sizeA, const int listB[], int sizeB, int destination[]) const;

private:
    int num_dimensions;
    int num_indexes;
    int *index;
};

class GranulatedIndexSet{ // optimized to add one index at a time
public:
    GranulatedIndexSet(int cnum_dimensions, int cnum_slots = 64);
    ~GranulatedIndexSet();

    int getNumDimensions() const;
    int getNumIndexes() const;

    void addIndex(const int p[]);
    int getSlot(const int p[]) const;

    const int* getIndex(int j) const;

    void addGranulatedSet(const GranulatedIndexSet *set);

    const int* getIndexes() const;
    const int* getMap() const;

protected:
    TypeIndexRelation compareIndexes(const int a[], const int b[]) const;

    void mergeMapped(const int newIndex[], const int newMap[], int sizeNew);

private:
    int num_dimensions;
    int num_indexes;
    int num_slots;
    int *index;
    int *map;
};

class DumpIndexSet{
public:
    DumpIndexSet(int cnum_dimensions, int initial_slots);
    ~DumpIndexSet();

    int getNumLoaded() const;

    void addIndex(const int p[]);
    int* ejectIndexes();

private:
    int num_dimensions, num_slots, num_loaded;
    int *index;
};

class IndexSet{ // rigid set but optimal in storage size
public:
    IndexSet(int cnum_dimensions);
    IndexSet(const UnsortedIndexSet *set);
    IndexSet(const GranulatedIndexSet *set);
    IndexSet(const IndexSet *set);
    IndexSet(int cnum_dimensions, int cnum_indexes, int* &cindex);
    ~IndexSet();

    int getNumDimensions() const;
    int getNumIndexes() const;

    void write(std::ofstream &ofs) const;
    void read(std::ifstream &ifs);

    void writeBinary(std::ofstream &ofs) const;
    void readBinary(std::ifstream &ifs);

    int getSlot(const int p[]) const;

    const int* getIndex(int i) const;

    void addUnsortedSet(const UnsortedIndexSet *set);
    void addGranulatedSet(const GranulatedIndexSet *set);
    void addIndexSet(const IndexSet *set);

    IndexSet* diffSets(const IndexSet *set) const; // returns this set minus the points in set

protected:
    TypeIndexRelation compareIndexes(const int a[], const int b[]) const;

    void merge(const int newIndex[], int sizeNew);
    void mergeMapped(const int newIndex[], const int map[], int sizeNew);

private:
    int num_dimensions;
    int num_indexes;
    int *index;
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
    double* aliasValues() const;
    int getNumOutputs() const;

    void setValues(const double vals[]);
    void setValuesPointer(double* &vals, int num_values);
    void addValues(const IndexSet *old_set, const IndexSet *new_set, const double new_vals[]);

protected:
    TypeIndexRelation compareIndexes(int num_dimensions, const int a[], const int b[]) const;

private:
    size_t num_outputs, num_values;
    double *values;
};

}

#endif
