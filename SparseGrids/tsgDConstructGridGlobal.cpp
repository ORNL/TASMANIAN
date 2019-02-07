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

#ifndef __TASMANIAN_SPARSE_GRID_DYNAMIC_CONST_GLOBAL_CPP
#define __TASMANIAN_SPARSE_GRID_DYNAMIC_CONST_GLOBAL_CPP

#include "tsgDConstructGridGlobal.hpp"

namespace TasGrid{

DynamicConstructorDataGlobal::DynamicConstructorDataGlobal(int cnum_dimensions, int cnum_outputs)
    : num_dimensions(cnum_dimensions), num_outputs(cnum_outputs){}
DynamicConstructorDataGlobal::~DynamicConstructorDataGlobal(){}

void DynamicConstructorDataGlobal::write(std::ofstream &ofs) const{
    std::vector<const TensorData*> tensor_refs;
    makeReverseReferenceVector<TensorData>(tensors, tensor_refs);

    ofs << tensor_refs.size() << std::endl;
    ofs << std::scientific;
    ofs.precision(16);
    for(auto d : tensor_refs){
        ofs << d->weight;
        for(auto i : d->tensor) ofs << " " << i;
        ofs << std::endl;
    }

    writeNodeDataList<std::ofstream, true>(data, ofs);
}
void DynamicConstructorDataGlobal::writeBinary(std::ofstream &ofs) const{
    std::vector<const TensorData*> tensor_refs;
    makeReverseReferenceVector<TensorData>(tensors, tensor_refs);

    int num = (int) tensor_refs.size();
    ofs.write((char*) &num, sizeof(int));
    for(auto d : tensor_refs){
        ofs.write((char*) &(d->weight), sizeof(double));
        ofs.write((char*) d->tensor.data(), num_dimensions * sizeof(int));
    }

    writeNodeDataList<std::ofstream, false>(data, ofs);
}
int DynamicConstructorDataGlobal::read(std::ifstream &ifs){
    int num, max_tensor = 0;
    ifs >> num; // get the number of tensors
    for(int i=0; i<num; i++){
        TensorData td;
        td.tensor.resize(num_dimensions);
        ifs >> td.weight;
        for(auto &t : td.tensor) ifs >> t;
        max_tensor = std::max(max_tensor, *std::max_element(td.tensor.begin(), td.tensor.end()));
        tensors.push_front(std::move(td));
    }
    ifs >> num; // get the number of data points
    for(int i=0; i<num; i++){
        NodeData nd;
        nd.point.resize(num_dimensions);
        for(auto &t : nd.point) ifs >> t;
        nd.value.resize(num_outputs);
        for(auto &t : nd.value) ifs >> t;
        data.push_front(std::move(nd));
    }
    return max_tensor;
}
int DynamicConstructorDataGlobal::readBinary(std::ifstream &ifs){
    int num, max_tensor = 0;
    ifs.read((char*) &num, sizeof(int)); // get the number of tensors
    for(int i=0; i<num; i++){
        TensorData td;
        ifs.read((char*) &td.weight, sizeof(double));
        td.tensor.resize(num_dimensions);
        ifs.read((char*) td.tensor.data(), num_dimensions * sizeof(int));
        max_tensor = std::max(max_tensor, *std::max_element(td.tensor.begin(), td.tensor.end()));
        tensors.push_front(std::move(td));
    }
    ifs.read((char*) &num, sizeof(int)); // get the number of data points
    for(int i=0; i<num; i++){
        NodeData nd;
        nd.point.resize(num_dimensions);
        ifs.read((char*) nd.point.data(), num_dimensions * sizeof(int));
        nd.value.resize(num_outputs);
        ifs.read((char*) nd.value.data(), num_outputs * sizeof(double));
        data.push_front(std::move(nd));
    }
    return max_tensor;
}

void DynamicConstructorDataGlobal::reloadPoints(std::function<int(int)> getNumPoints){
    auto t = tensors.begin();
    while(t != tensors.end()){
        std::vector<int> v = t->tensor;
        MultiIndexSet dummy_set(num_dimensions);
        dummy_set.setIndexes(v);
        t->points.setNumDimensions(num_dimensions);
        MultiIndexManipulations::generateNestedPoints(dummy_set, getNumPoints, t->points);
        t->loaded = std::vector<bool>(t->points.getNumIndexes(), false);
        t++;
    }

    auto p = data.begin(); // rebuild the loaded vectors
    while(p != data.end()){
        t = tensors.begin();
        while(t != tensors.end()){
            int i = t->points.getSlot(p->point);
            if (i != -1) t->loaded[i] = true;
            t++;
        }
        p++;
    }
}

void DynamicConstructorDataGlobal::clearTesnors(){
    auto p = tensors.before_begin();
    auto t = tensors.begin();
    while(t != tensors.end()){
        if (t->weight >= 0.0){
            tensors.erase_after(p);
            t = p;
        }else{
            p++;
        }
        t++;
    }
}

void DynamicConstructorDataGlobal::getInitialTensors(MultiIndexSet &set) const{
    Data2D<int> tens;
    tens.resize(num_dimensions, 0);
    auto t = tensors.begin();
    while(t != tensors.end()){
        if (t->weight < 0.0){
            tens.appendStrip((*t).tensor);
        }
        t++;
    }
    set = MultiIndexSet(num_dimensions);
    if (tens.getNumStrips() > 0)
        set.addData2D(tens);
}

void DynamicConstructorDataGlobal::addTensor(const int *tensor, std::function<int(int)> getNumPoints, double weight){
    TensorData t;
    std::vector<int> v(tensor, tensor + num_dimensions);
    t.tensor = v;
    MultiIndexSet dummy_set(num_dimensions);
    dummy_set.setIndexes(v);
    t.points.setNumDimensions(num_dimensions);
    MultiIndexManipulations::generateNestedPoints(dummy_set, getNumPoints, t.points);
    t.loaded = std::vector<bool>(t.points.getNumIndexes(), false);
    auto d = data.begin();
    while(d != data.end()){
        int slot = t.points.getSlot(d->point);
        if (slot != -1) t.loaded[slot] = true;
        d++;
    }
    t.weight = weight;
    tensors.push_front(std::move(t));
}

void DynamicConstructorDataGlobal::getNodesIndexes(std::vector<int> &inodes){
    inodes = std::vector<int>();
    tensors.sort([&](const TensorData &a, const TensorData &b)->bool{ return (a.weight < b.weight); });
    auto t = tensors.begin();
    while(t != tensors.end()){
        if (!t->loaded.empty()){
            for(int i=0; i<t->points.getNumIndexes(); i++){
                if (!t->loaded[i])
                    inodes.insert(inodes.end(), t->points.getIndex(i), t->points.getIndex(i) + num_dimensions);
            }
        }
        t++;
    }
}

bool DynamicConstructorDataGlobal::addNewNode(const std::vector<int> &point, const std::vector<double> &value){
    data.push_front({point, value});
    auto t = tensors.begin();
    while(t != tensors.end()){
        int slot = t->points.getSlot(point);
        if (slot != -1){
            t->loaded[slot] = true;
            if (std::all_of(t->loaded.begin(), t->loaded.end(), [](bool a)-> bool{ return a; })){ // if all entries are true
                t->loaded = std::vector<bool>();
                return true; // all points associated with this tensor have been loaded, signal to call ejectCompleteTensor()
            }else{
                return false; // point added to the tensor, but the tensor is not complete yet
            }
        }
        t++;
    }
    return false; // could not find a tensor that contains this node???
}

bool DynamicConstructorDataGlobal::ejectCompleteTensor(const MultiIndexSet &current_tensors, std::vector<int> &tensor, MultiIndexSet &points, std::vector<double> &vals){
    auto p = tensors.before_begin();
    auto t = tensors.begin();
    while(t != tensors.end()){
        if (t->loaded.empty()){ // empty loaded means all have been loaded
            if (MultiIndexManipulations::isLowerComplete(t->tensor, current_tensors)){
                tensor = t->tensor;
                points = t->points;
                vals.resize(points.getNumIndexes() *  num_outputs);
                auto d = data.before_begin();
                auto v = data.begin();
                int found = 0;
                while(found < points.getNumIndexes()){
                    int slot = points.getSlot(v->point);
                    if (slot == -1){
                        d++;
                        v++;
                    }else{
                        std::copy_n(v->value.begin(), num_outputs, &(vals[slot * num_outputs]));
                        data.erase_after(d);
                        v = d;
                        v++;
                        found++;
                    }
                }
                tensors.erase_after(p);
                return true;
            }
        }
        t++;
        p++;
    }
    return false; // could not find a complete tensor with parents present in tensors
}

}

#endif
