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

#ifndef __TASMANIAN_SPARSE_GRID_DYNAMIC_CONST_GLOBAL_HPP
#define __TASMANIAN_SPARSE_GRID_DYNAMIC_CONST_GLOBAL_HPP

#include <forward_list>

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgIndexManipulator.hpp"

//! \file tsgDConstructGridGlobal.hpp
//! \brief Class and data types used for dynamic construction of Global grids.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianRefinement
//!
//! Defines the storage class and struct types used for
//! construction of Global grids. The class is more complex than other
//! examples because rules with growth of more than one point need
//! to collect several samples before the grid can be extended with another
//! tensor.

//! \internal
//! \defgroup TasmanianRefinement Refinement related classes, structures, and methods.
//!
//! \par Adaptive Refinement
//! The most efficient sparse grids approximation strategies are adaptive,
//! i.e., the grid is constructed to best capture the behavior of the specific
//! target model. Adapting a grid requires analysis of the existing loaded
//! model values followed by a decision on how to best expand the grid.
//! Sometimes the complexity of the refinement algorithms requires additional
//! data storage and methods implemented outside of the main grid classes.

namespace TasGrid{

//! \internal
//! \brief Holds the pair of point index and model value, the struct is used in a \b std::forward_list.
//! \ingroup TasmanianRefinement

//! Global grids can be expanded by a full tensor index at a time,
//! which requires model outputs data from multiple nodes.
//! However, the user provides the model data only a single node at a time.
//! Thus, node data has to be stored separately until it can be included in the grid.
//!
//! The same data-structure is used by the sequence grids
//! to store the initial set of samples until a lower set of point can be constructed.
struct NodeData{
    //! \brief The multi-index of the point.
    std::vector<int> point;
    //! \brief The values of the model outputs at the point.
    std::vector<double> value;
};

//! \internal
//! \brief Holds the description of a single tensor, points associated with the delta and whether model data has been loaded.
//! \ingroup TasmanianRefinement

//! Each candidate tensor is described by the tensor multi-index and the points associated with the corresponding delta.
//! Using the delta definition guarantees that each point appears in only one tensor.
//! A vector of booleans indicates which of the points have been \b loaded,
//! i.e., whether the user has already provided model data for the corresponding points.
//! If the \b loaded vector is empty, then all the needed data has been loaded,
//! i.e., equivalent to all entries being \b true but faster than traversing a vector.
//! The tensor \b weight is also stored here, where negative weights are used to indicate the initial sparse grid
//! which ensures that the initial set of nodes will always come first in order.
//! The tensor data is stored in a \b std::forward_list for easy random insertions and deletions.
struct TensorData{
    //! \brief Multi-index of the tensor.
    std::vector<int> tensor;
    //! \brief The points associated with the delta multi-index.
    MultiIndexSet points;
    //! \brief For each point \b false if the model data is not yet available for the point, true otherwise.
    std::vector<bool> loaded;
    //! \brief The weight indicates the relative importance of the tensor.
    double weight;
};

//! \internal
//! \brief Takes a list and creates a vector of references in reversed order, needed for I/O so that read can be followed by push_front.
//! \ingroup TasmanianRefinement
template <class T>
void makeReverseReferenceVector(const std::forward_list<T> &list, std::vector<const T*> &refs){
    size_t num_entries = (size_t) std::distance(list.begin(), list.end());
    refs.resize(num_entries);
    auto p = list.begin();
    auto r = refs.rbegin();
    while(p != list.end()) *r++ = &*p++;
}

//! \internal
//! \brief Writes a NodeData forward_list to a file using the lambdas to write ints and doubles.
//! \ingroup TasmanianRefinement
template<class STREAMCONCEPT>
void writeNodeDataList(const std::forward_list<NodeData> &data, STREAMCONCEPT &ofs, std::function<void(int, STREAMCONCEPT &)> write_an_int,
                       std::function<void(const NodeData *, STREAMCONCEPT &)> write_node){
    std::vector<const NodeData*> data_refs;
    makeReverseReferenceVector<NodeData>(data, data_refs);

    write_an_int((int) data_refs.size(), ofs);
    for(auto d : data_refs) write_node(d, ofs);
}

//! \internal
//! \brief Writes a NodeData forward_list to a file using either binary or ascii format.
//! \ingroup TasmanianRefinement
template<class STREAMCONCEPT, bool useAscii>
void writeNodeDataList(const std::forward_list<NodeData> &data, STREAMCONCEPT &ofs){
    if (useAscii){
        writeNodeDataList<STREAMCONCEPT>(data, ofs, [&](int i, STREAMCONCEPT &f)->void{ f << i << std::endl; },
                    [&](const NodeData *n, STREAMCONCEPT &f)->void{
                        for(auto i : n->point) f << i << " ";
                        f << n->value[0]; // do now allow extra blank spaces
                        for(size_t i = 1; i < n->value.size(); i++) f << " " << n->value[i];
                        f << std::endl;
                    });
    }else{
        writeNodeDataList<STREAMCONCEPT>(data, ofs, [&](int i, STREAMCONCEPT &f)->void{ f.write((char*) &i, sizeof(int)); },
                [&](const NodeData *n, STREAMCONCEPT &f)->void{
                    f.write((char*) n->point.data(), n->point.size() * sizeof(int));
                    f.write((char*) n->value.data(), n->value.size() * sizeof(double));
                });
    }
}

//! \internal
//! \brief Reads a NodeData forward_list from a file using either binary or ascii format.
//! \ingroup TasmanianRefinement
template<class STREAMCONCEPT, bool useAscii>
void readNodeDataList(int num_dimensions, int num_outputs, STREAMCONCEPT &ifs, std::forward_list<NodeData> &data){
    int num_nodes;
    if (useAscii){
        ifs >> num_nodes; // get the number of data points
        for(int i=0; i<num_nodes; i++){
            NodeData nd;
            nd.point.resize(num_dimensions);
            for(auto &t : nd.point) ifs >> t;
            nd.value.resize(num_outputs);
            for(auto &t : nd.value) ifs >> t;
            data.push_front(std::move(nd));
        }
    }else{
        ifs.read((char*) &num_nodes, sizeof(int)); // get the number of data points
        for(int i=0; i<num_nodes; i++){
            NodeData nd;
            nd.point.resize(num_dimensions);
            ifs.read((char*) nd.point.data(), num_dimensions * sizeof(int));
            nd.value.resize(num_outputs);
            ifs.read((char*) nd.value.data(), num_outputs * sizeof(double));
            data.push_front(std::move(nd));
        }
    }
}


//! \internal
//! \brief Helper class that stores data from dynamic construction of a Global grid.
//! \ingroup TasmanianRefinement

//! The class stores candidate tensors with corresponding weights
//! as well as the model data provided by the user.
//! When enough data has been computed to complete a tensor,
//! the tensor can be ejected with the points and model data returned in a format
//! that is easy to incorporate within the data structures of the \b GridGlobal class.
class DynamicConstructorDataGlobal{
public:
    //! \brief The only constructor, requires that the dimension and the number of model outputs is specified.
    DynamicConstructorDataGlobal(int cnum_dimensions, int cnum_outputs);
    //! \brief Default destructor, release all used memory and erase all stored data.
    ~DynamicConstructorDataGlobal();

    //! \brief Write the data to a file in ascii format.
    void write(std::ofstream &ofs) const;
    //! \brief Write the data to a file in binary format.
    void writeBinary(std::ofstream &ofs) const;
    //! \brief Read the data from a file in ascii format, returns the maximum tensor index so the 1D wrapper cache structure can be updated.
    int read(std::ifstream &ifs);
    //! \brief Read the data from a file in binary format, returns the maximum tensor index so the 1D wrapper cache structure can be updated.
    int readBinary(std::ifstream &ifs);

    //! \brief Called after read, reinitializes the points and loaded structures for the tensors.
    void reloadPoints(std::function<int(int)> getNumPoints);

    //! \brief Delete the tensors with non-negative weights, i.e., clear all but the tensors selected by the initial grid.
    void clearTesnors();

    //! \brief Get a set of all tensors with negative weight, i.e., the tensors selected for the initial run.
    void getInitialTensors(MultiIndexSet &set) const;

    //! \brief Add a new tensor to the candidates, with the given \b weight and using \b getNumPoints() to define the associated nodes.
    void addTensor(const int *tensor, std::function<int(int)> getNumPoints, double weight);

    //! \brief Get the node indexes of the points associated with the candidate tensors, the order is the same as the tensors sorted by ascending weight.
    void getNodesIndexes(std::vector<int> &inodes);

    //! \brief Add a new data point with the index and the value, returns \b true if there is enough data to complete a tensor.
    bool addNewNode(const std::vector<int> &point, const std::vector<double> &value); // returns whether a tensor is complete

    //! \brief Return a completed tensor with parent-tensors included in tensors, returns \b true if such tensor has been found.
    bool ejectCompleteTensor(const MultiIndexSet &current_tensors, std::vector<int> &tensor, MultiIndexSet &points, std::vector<double> &vals);

private:
    int num_dimensions, num_outputs;
    std::forward_list<NodeData> data;
    std::forward_list<TensorData> tensors;
};

}

#endif
