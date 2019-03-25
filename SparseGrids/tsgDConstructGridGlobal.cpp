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

DynamicConstructorDataGlobal::DynamicConstructorDataGlobal(size_t cnum_dimensions, size_t cnum_outputs)
    : num_dimensions(cnum_dimensions), num_outputs(cnum_outputs){}
DynamicConstructorDataGlobal::~DynamicConstructorDataGlobal(){}

template<bool useAscii> void DynamicConstructorDataGlobal::write(std::ostream &os) const{
    if (useAscii){ os << std::scientific; os.precision(17); }
    auto tensor_refs =  makeReverseReferenceVector(tensors);

    IO::writeNumbers<useAscii, IO::pad_line, int>(os, (int) tensor_refs.size());
    for(auto d : tensor_refs){
        IO::writeNumbers<useAscii, IO::pad_rspace, double>(os, d->weight);
        IO::writeVector<useAscii, IO::pad_line>(d->tensor, os);
    }
    writeNodeDataList<useAscii>(data, os);
}

template<bool useAscii> void DynamicConstructorDataGlobal::read(std::istream &is){
    int num_entries = IO::readNumber<useAscii, int>(is);

    for(int i=0; i<num_entries; i++){
        TensorData td;
        td.tensor.resize(num_dimensions);
        td.weight = IO::readNumber<useAscii, double>(is);
        IO::readVector<useAscii>(is, td.tensor);
        tensors.push_front(std::move(td));
    }

    data = readNodeDataList<useAscii>(num_dimensions, num_outputs, is);
}

template void DynamicConstructorDataGlobal::write<true>(std::ostream &) const; // instantiate for faster build
template void DynamicConstructorDataGlobal::write<false>(std::ostream &) const;
template void DynamicConstructorDataGlobal::read<true>(std::istream &);
template void DynamicConstructorDataGlobal::read<false>(std::istream &);

int DynamicConstructorDataGlobal::getMaxTensor() const{
    int max_tensor = 0;
    for(auto t = tensors.begin(); t != tensors.end(); t++)
        max_tensor = std::max(max_tensor, *std::max_element(t->tensor.begin(), t->tensor.end()));
    return max_tensor;
}

void DynamicConstructorDataGlobal::reloadPoints(std::function<int(int)> getNumPoints){
    auto t = tensors.begin();
    while(t != tensors.end()){
        std::vector<int> v = t->tensor;
        MultiIndexSet dummy_set(num_dimensions, v);
        t->points = MultiIndexManipulations::generateNestedPoints(dummy_set, getNumPoints);
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

MultiIndexSet DynamicConstructorDataGlobal::getInitialTensors() const{
    Data2D<int> tens(num_dimensions, 0);
    auto t = tensors.begin();
    while(t != tensors.end()){
        if (t->weight < 0.0){
            tens.appendStrip((*t).tensor);
        }
        t++;
    }
    return MultiIndexSet(tens);
}

void DynamicConstructorDataGlobal::addTensor(const int *tensor, std::function<int(int)> getNumPoints, double weight){
    TensorData t;
    std::vector<int> v(tensor, tensor + num_dimensions);
    t.tensor = v;
    MultiIndexSet dummy_set(num_dimensions, v);
    t.points = MultiIndexManipulations::generateNestedPoints(dummy_set, getNumPoints);
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
