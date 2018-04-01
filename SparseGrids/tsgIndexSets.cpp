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

UnsortedIndexSet::UnsortedIndexSet(int cnum_dimensions, int cnum_slots) : num_dimensions(cnum_dimensions), num_indexes(0){
    index = new int[num_dimensions * cnum_slots];
}
UnsortedIndexSet::~UnsortedIndexSet(){
    delete[] index;
}

int UnsortedIndexSet::getNumDimensions() const{ return num_dimensions; }
int UnsortedIndexSet::getNumIndexes() const{ return num_indexes; }

void UnsortedIndexSet::addIndex(const int p[]){
    std::copy(p, p + num_dimensions, &(index[num_dimensions * num_indexes++]));
}
const int* UnsortedIndexSet::getIndex(int i) const{
    return &(index[i * num_dimensions]);
}

int* UnsortedIndexSet::getIndexesSorted() const{
    int *list_source = new int[num_indexes];
    int *list_destination = new int[num_indexes];
    int *list_temp = 0;

    for(int i=0; i<num_indexes - num_indexes%2; i+=2){
        if (compareIndexes(&(index[i*num_dimensions]), &(index[(i+1)*num_dimensions])) == type_abeforeb){
            list_source[i] = i;
            list_source[i+1] = i+1;
        }else{
            list_source[i] = i+1;
            list_source[i+1] = i;
        }
    }
    if (num_indexes%2 == 1){
        list_source[num_indexes-1] = num_indexes-1;
    }

    int warp = 2;
    while(warp < num_indexes){
        int full_warps = 2*warp * ((num_indexes) / (2*warp));
        #pragma omp parallel for
        for(int i=0; i<full_warps; i+=2*warp){
            merge(&(list_source[i]), warp, &(list_source[i+warp]), warp, &(list_destination[i]));
        }

        if (full_warps < num_indexes){
            merge(&(list_source[full_warps]), ((full_warps+warp < num_indexes) ? warp : num_indexes - full_warps), &(list_source[full_warps+warp]), ((full_warps+2*warp < num_indexes) ? warp : num_indexes - full_warps - warp), &(list_destination[full_warps]));
        }
        list_temp = list_destination;
        list_destination = list_source;
        list_source = list_temp;
        list_temp = 0;
        warp *= 2;
    }

    int *result = new int[num_indexes * num_dimensions];
    for(int i=0; i<num_indexes; i++){
        std::copy(&(index[list_source[i] * num_dimensions]), &(index[list_source[i] * num_dimensions]) + num_dimensions, &(result[i*num_dimensions]));
    }

    delete[] list_source;
    delete[] list_destination;

    return result;
}

TypeIndexRelation UnsortedIndexSet::compareIndexes(const int a[], const int b[]) const{
    for(int i=0; i<num_dimensions; i++){
        if (a[i] < b[i]) return type_abeforeb;
        if (a[i] > b[i]) return type_bbeforea;
    }
    return type_asameb;
}

void UnsortedIndexSet::merge(const int listA[], int sizeA, const int listB[], int sizeB, int destination[]) const{
    TypeIndexRelation relation;
    int offset_destination = 0, offsetA = 0, offsetB = 0;
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

GranulatedIndexSet::GranulatedIndexSet(int cnum_dimensions, int cnum_slots) : num_dimensions(cnum_dimensions), num_indexes(0), num_slots(cnum_slots){
    index = new int[num_dimensions * num_slots];
    map = new int[num_slots];
}

GranulatedIndexSet::~GranulatedIndexSet(){
    delete[] index;
    delete[] map;
}

int GranulatedIndexSet::getNumDimensions() const{  return num_dimensions;  }
int GranulatedIndexSet::getNumIndexes() const{  return num_indexes;  }

int GranulatedIndexSet::getSlot(const int p[]) const{
    int start = 0, end = num_indexes - 1;
    int current = (start + end) / 2;
    TypeIndexRelation t;
    while (start <= end){
        t = compareIndexes(&(index[num_dimensions * map[current]]), p);
        if (t == type_abeforeb){
            start = current+1;
        }else if (t == type_bbeforea){
            end = current-1;
        }else{
            return current;
        };
        current = (start + end) / 2;
    }
    return -1;
}

const int* GranulatedIndexSet::getIndex(int j) const{
    return &(index[map[j] * num_dimensions]);
}

void GranulatedIndexSet::addIndex(const int p[]){
    int start = 0, end = num_indexes - 1;
    int current = (start + end) / 2;
    TypeIndexRelation t;
    while (start <= end){
        t = compareIndexes(&(index[num_dimensions * map[current]]), p);
        if (t == type_abeforeb){
            start = current+1;
        }else if (t == type_bbeforea){
            end = current-1;
        }else{
            return;
        };
        current = (start + end) / 2;
    }

    if (num_indexes == num_slots){
        int *old_index = index;
        int *old_map = map;
        num_slots *= 2;

        index = new int[num_dimensions * num_slots];
        map = new int[num_slots];

        std::copy(old_index, old_index + num_dimensions * num_indexes, index);
        std::copy(old_map, old_map + num_indexes, map);

        delete[] old_index;
        delete[] old_map;
    }

    std::copy(p, p + num_dimensions, &(index[num_indexes * num_dimensions]));

    int slot = start+1; slot = (slot > num_indexes) ? num_indexes : slot;
    while((slot>0) && (compareIndexes(p, &(index[map[slot-1] * num_dimensions])) == type_abeforeb)){ slot--; }

    for(int i=num_indexes; i>slot; i--){  map[i] = map[i-1];  }
    map[slot] = num_indexes;
    num_indexes++;
}

void GranulatedIndexSet::addGranulatedSet(const GranulatedIndexSet *set){
    mergeMapped(set->getIndexes(), set->getMap(), set->getNumIndexes());
}
const int* GranulatedIndexSet::getIndexes() const{
    return index;
}
const int* GranulatedIndexSet::getMap() const{
    return map;
}

TypeIndexRelation GranulatedIndexSet::compareIndexes(const int a[], const int b[]) const{
    for(int i=0; i<num_dimensions; i++){
        if (a[i] < b[i]) return type_abeforeb;
        if (a[i] > b[i]) return type_bbeforea;
    }
    return type_asameb;
}

void GranulatedIndexSet::mergeMapped(const int newIndex[], const int newMap[], int sizeNew){
    int *oldIndex = index;
    int sizeOld = num_indexes;
    num_slots = sizeOld + sizeNew;

    index = new int[num_dimensions * num_slots];

    TypeIndexRelation relation;
    num_indexes = 0;
    int offsetNew = 0, offsetOld = 0;
    while( (offsetNew < sizeNew) || (offsetOld < sizeOld) ){
        if (offsetNew >= sizeNew){ // new is a
            relation = type_bbeforea;
        }else if (offsetOld >= sizeOld){ // old is b
            relation = type_abeforeb;
        }else{
            relation = compareIndexes(&(newIndex[newMap[offsetNew] * num_dimensions]), &(oldIndex[map[offsetOld] * num_dimensions]));
        }
        if (relation == type_abeforeb){
            std::copy(&(newIndex[newMap[offsetNew] * num_dimensions]), &(newIndex[newMap[offsetNew] * num_dimensions]) + num_dimensions, &(index[num_indexes++ * num_dimensions]));
            offsetNew++;
        }
        if ((relation == type_bbeforea) || (relation == type_asameb)){
            std::copy(&(oldIndex[map[offsetOld] * num_dimensions]), &(oldIndex[map[offsetOld] * num_dimensions]) + num_dimensions, &(index[num_indexes++ * num_dimensions]));
            offsetOld++;
            offsetNew += (relation == type_asameb) ? 1 : 0;
        }
    }

    delete[] oldIndex;
    delete[] map;

    map = new int[num_slots];
    for(int i=0; i<num_indexes; i++){  map[i] = i;  }
}

DumpIndexSet::DumpIndexSet(int cnum_dimensions, int initial_slots) :
    num_dimensions(cnum_dimensions), num_slots(initial_slots), num_loaded(0), index(0){
    index = new int[num_dimensions * num_slots];
}
DumpIndexSet::~DumpIndexSet(){
    if (index != 0){ delete[] index; index = 0; }
}
int DumpIndexSet::getNumLoaded() const{ return num_loaded; }
void DumpIndexSet::addIndex(const int p[]){
    if (num_loaded == num_slots){
        int *old_index = index;
        index = new int[2 * num_dimensions * num_slots];
        std::copy(old_index, old_index + num_slots * num_dimensions, index);
        num_slots *= 2;
        delete[] old_index;
    }
    std::copy(p, p + num_dimensions, &(index[num_loaded * num_dimensions]));
    num_loaded++;
}
int* DumpIndexSet::ejectIndexes(){
    int *res;
    if (num_slots == num_loaded){
        res = index;
        index = 0;
        num_loaded = 0;
        return res;
    }else{
        res = new int[num_dimensions * num_loaded];
        std::copy(index, index + num_dimensions * num_loaded, res);
    }
    return res;
}

IndexSet::IndexSet(int cnum_dimensions) : num_dimensions(cnum_dimensions), num_indexes(0), index(0){}
IndexSet::IndexSet(const UnsortedIndexSet *set){
    num_dimensions = set->getNumDimensions();
    num_indexes = set->getNumIndexes();
    index = set->getIndexesSorted();
}
IndexSet::IndexSet(const GranulatedIndexSet *set){
    num_dimensions = set->getNumDimensions();
    num_indexes = set->getNumIndexes();
    index = new int[num_dimensions * num_indexes];
    for(int i=0; i<num_indexes; i++){
        const int *p = set->getIndex(i);
        std::copy(p, p + num_dimensions, &(index[i*num_dimensions]));
    }
}
IndexSet::IndexSet(const IndexSet *set){
    num_dimensions = set->getNumDimensions();
    num_indexes = set->getNumIndexes();
    index = new int[num_dimensions * num_indexes];
    const int *p = set->getIndex(0);
    std::copy(p, p + num_dimensions * num_indexes, index);
}

IndexSet::IndexSet(int cnum_dimensions, int cnum_indexes, int* &cindex) : num_dimensions(cnum_dimensions), num_indexes(cnum_indexes), index(cindex) {
    cindex = 0;
}
IndexSet::~IndexSet(){
    if (index != 0){  delete[] index;  }
}
int IndexSet::getNumDimensions() const{
    return num_dimensions;
}
int IndexSet::getNumIndexes() const{  return num_indexes;  }

void IndexSet::write(std::ofstream &ofs) const{
    ofs << num_dimensions << " " << num_indexes;
    for(int i=0; i<num_dimensions*num_indexes; i++){
        ofs << " " << index[i];
    }
    ofs << endl;
}
void IndexSet::read(std::ifstream &ifs){
    if (index != 0){ delete[] index; }
    ifs >> num_dimensions >> num_indexes;
    index = new int[num_dimensions * num_indexes];
    for(int i=0; i<num_dimensions*num_indexes; i++){
        ifs >> index[i];
    }
}

void IndexSet::writeBinary(std::ofstream &ofs) const{
    int sizes[2];
    sizes[0] = num_dimensions;
    sizes[1] = num_indexes;
    ofs.write((char*) sizes, 2*sizeof(int));
    ofs.write((char*) index, num_dimensions*num_indexes*sizeof(int));
}
void IndexSet::readBinary(std::ifstream &ifs){
    if (index != 0){ delete[] index; }
    int sizes[2];
    ifs.read((char*) sizes, 2*sizeof(int));
    num_dimensions = sizes[0];
    num_indexes = sizes[1];
    index = new int[num_dimensions * num_indexes];
    ifs.read((char*) index, num_dimensions*num_indexes*sizeof(int));
}

int IndexSet::getSlot(const int p[]) const{
    int start = 0, end = num_indexes - 1;
    int current = (start + end) / 2;
    TypeIndexRelation t;
    while (start <= end){
        t = compareIndexes(&(index[num_dimensions * current]), p);
        if (t == type_abeforeb){
            start = current+1;
        }else if (t == type_bbeforea){
            end = current-1;
        }else{
            return current;
        };
        current = (start + end) / 2;
    }
    return -1;
}
const int* IndexSet::getIndex(int i) const{
    return &(index[i*num_dimensions]);
}
void IndexSet::addUnsortedSet(const UnsortedIndexSet *set){
    int *set_index = set->getIndexesSorted();
    int set_num_indexes = set->getNumIndexes();
    merge(set_index, set_num_indexes);
    delete[] set_index;
}
void IndexSet::addGranulatedSet(const GranulatedIndexSet *set){
    const int *set_index = set->getIndexes();
    const int *set_map = set->getMap();
    int set_num_indexes = set->getNumIndexes();
    mergeMapped(set_index, set_map, set_num_indexes);
}
void IndexSet::addIndexSet(const IndexSet *set){
    const int *set_index = set->getIndex(0);
    int set_num_indexes = set->getNumIndexes();
    merge(set_index, set_num_indexes);
}

IndexSet* IndexSet::diffSets(const IndexSet *set) const{
    int *slots = new int[num_indexes];

    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_indexes; i++){
        slots[i] = set->getSlot(&(index[i*num_dimensions]));
    }

    int n = 0;
    for(int i=0; i<num_indexes; i++){
        if (slots[i] == -1) n++;
    }

    if (n == 0){
        delete[] slots;
        return 0;
    }else{
        int *diff_index = new int[num_dimensions * n];
        n = 0;
        for(int i=0; i<num_indexes; i++){
            if (slots[i] == -1) std::copy(&(index[i*num_dimensions]), &(index[i*num_dimensions]) + num_dimensions, &(diff_index[num_dimensions * n++]));
        }
        delete[] slots;

        IndexSet *diff_set = new IndexSet(num_dimensions, n, diff_index);
        return diff_set;
    }
}


TypeIndexRelation IndexSet::compareIndexes(const int a[], const int b[]) const{
    for(int i=0; i<num_dimensions; i++){
        if (a[i] < b[i]) return type_abeforeb;
        if (a[i] > b[i]) return type_bbeforea;
    }
    return type_asameb;
}
void IndexSet::merge(const int newIndex[], int sizeNew){
    int *oldIndex = index;
    int sizeOld = num_indexes;

    index = new int[num_dimensions * (sizeOld + sizeNew)];

    TypeIndexRelation relation;
    num_indexes = 0;
    int offsetNew = 0, offsetOld = 0;
    while( (offsetNew < sizeNew) || (offsetOld < sizeOld) ){
        if (offsetNew >= sizeNew){ // new is a
            relation = type_bbeforea;
        }else if (offsetOld >= sizeOld){ // old is b
            relation = type_abeforeb;
        }else{
            relation = compareIndexes(&(newIndex[offsetNew * num_dimensions]), &(oldIndex[offsetOld * num_dimensions]));
        }
        if (relation == type_abeforeb){
            std::copy(&(newIndex[offsetNew * num_dimensions]), &(newIndex[offsetNew * num_dimensions]) + num_dimensions, &(index[num_indexes++ * num_dimensions]));
            offsetNew++;
        }else{
            std::copy(&(oldIndex[offsetOld * num_dimensions]), &(oldIndex[offsetOld * num_dimensions]) + num_dimensions, &(index[num_indexes++ * num_dimensions]));
            offsetOld++;
            offsetNew += (relation == type_asameb) ? 1 : 0;
        }
    }

    delete[] oldIndex;
}

void IndexSet::mergeMapped(const int newIndex[], const int map[], int sizeNew){
    int *oldIndex = index;
    int sizeOld = num_indexes;

    index = new int[num_dimensions * (sizeOld + sizeNew)];

    TypeIndexRelation relation;
    num_indexes = 0;
    int offsetNew = 0, offsetOld = 0;
    while( (offsetNew < sizeNew) || (offsetOld < sizeOld) ){
        if (offsetNew >= sizeNew){ // new is a
            relation = type_bbeforea;
        }else if (offsetOld >= sizeOld){ // old is b
            relation = type_abeforeb;
        }else{
            relation = compareIndexes(&(newIndex[map[offsetNew] * num_dimensions]), &(oldIndex[offsetOld * num_dimensions]));
        }
        if (relation == type_abeforeb){
            std::copy(&(newIndex[map[offsetNew] * num_dimensions]), &(newIndex[map[offsetNew] * num_dimensions]) + num_dimensions, &(index[num_indexes++ * num_dimensions]));
            offsetNew++;
        }else{
            std::copy(&(oldIndex[offsetOld * num_dimensions]), &(oldIndex[offsetOld * num_dimensions]) + num_dimensions, &(index[num_indexes++ * num_dimensions]));
            offsetOld++;
            offsetNew += (relation == type_asameb) ? 1 : 0;
        }
    }

    delete[] oldIndex;
}

StorageSet::StorageSet(int cnum_outputs, int cnum_values) : num_outputs(cnum_outputs), num_values(cnum_values), values(0){}
StorageSet::StorageSet(const StorageSet *storage) : values(0){
    num_outputs = storage->num_outputs;
    num_values = storage->num_values;
    if (storage->values != 0){
        values = new double[num_outputs * num_values];
        std::copy(storage->values, storage->values + num_outputs*num_values, values);
    }
}
StorageSet::~StorageSet(){
    if (values != 0){ delete[] values; }
}

void StorageSet::write(std::ofstream &ofs) const{
    ofs << num_outputs << " " << num_values;
    if (values != 0){
        ofs << " 1";
        ofs << std::scientific; ofs.precision(17);
        for(size_t i=0; i<num_outputs*num_values; i++){
            ofs << " " << values[i];
        }
    }else{
        ofs << " 0";
    }
    ofs << endl;
}
void StorageSet::read(std::ifstream &ifs){
    if (values != 0){ delete[] values;  values = 0; }
    int has_vals;
    ifs >> num_outputs >> num_values >> has_vals;
    if (has_vals == 1){
        values = new double[num_outputs * num_values];
        for(size_t i=0; i<num_outputs*num_values; i++){
            ifs >> values[i];
        }
    }
}
void StorageSet::writeBinary(std::ofstream &ofs) const{
    int num_out_vals[2];
    num_out_vals[0] = (int) num_outputs;
    num_out_vals[1] = (int) num_values;
    ofs.write((char*) num_out_vals, 2*sizeof(int));
    if (values != 0){
        char flag = 'y'; ofs.write((char*) &flag, sizeof(char));
        ofs.write((char*) values, num_outputs * num_values * sizeof(double));
    }else{
        char flag = 'n'; ofs.write((char*) &flag, sizeof(char));
    }
}
void StorageSet::readBinary(std::ifstream &ifs){
    if (values != 0){ delete[] values;  values = 0; }
    int num_out_vals[2];
    ifs.read((char*) num_out_vals, 2*sizeof(int));
    num_outputs = (size_t) num_out_vals[0];
    num_values = (size_t) num_out_vals[1];
    char flag; ifs.read((char*) &flag, sizeof(char));
    if (flag == 'y'){
        values = new double[num_outputs * num_values];
        ifs.read((char*) values, num_outputs * num_values * sizeof(double));
    }
}

int StorageSet::getNumOutputs() const{ return num_outputs; }
const double* StorageSet::getValues(int i) const{ return &(values[i*num_outputs]); }
double* StorageSet::aliasValues() const{ return values; }

void StorageSet::setValues(const double vals[]){
    if (values == 0){ values = new double[num_outputs * num_values]; }
    std::copy(vals, vals + num_values * num_outputs, values);
}
void StorageSet::setValuesPointer(double* &vals, int cnum_values){
    num_values = cnum_values;
    if (values != 0) delete[] values;
    values = vals;
    vals = 0;
}

void StorageSet::addValues(const IndexSet *old_set, const IndexSet *new_set, const double new_vals[]){
    double *old_vals = values;
    size_t old_num_values = num_values;
    size_t new_num_values = new_set->getNumIndexes();
    int num_dimensions = old_set->getNumDimensions();

    TypeIndexRelation relation;
    values = new double[num_outputs * (old_num_values + new_num_values)];
    size_t offsetOld = 0, offsetNew = 0;
    num_values = 0;
    while( (offsetNew < new_num_values) || (offsetOld < old_num_values) ){
        if (offsetNew >= new_num_values){ // new is a
            relation = type_bbeforea;
        }else if (offsetOld >= old_num_values){ // old is b
            relation = type_abeforeb;
        }else{
            relation = compareIndexes(num_dimensions, new_set->getIndex(offsetNew), old_set->getIndex(offsetOld));
        }
        if (relation == type_abeforeb){
            std::copy(&(new_vals[num_outputs * offsetNew]), &(new_vals[num_outputs * offsetNew]) + num_outputs, &(values[num_outputs * num_values++]));
            offsetNew++;
        }else{
            std::copy(&(old_vals[num_outputs * offsetOld]), &(old_vals[num_outputs * offsetOld]) + num_outputs, &(values[num_outputs * num_values++]));
            offsetOld++;
        }
    }

    delete[] old_vals;
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
