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

#include "tsgIndexManipulator.hpp"

/*!
 * \internal
 * \file tsgDConstructGridGlobal.hpp
 * \brief Class and data types used for dynamic construction of Global grids.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianRefinement
 *
 * Defines the storage class and struct types used for
 * construction of Global grids. The class is more complex than other
 * examples because rules with growth of more than one point need
 * to collect several samples before the grid can be extended with another
 * tensor.
 * \endinternal
 */

/*!
 * \internal
 * \ingroup TasmanianSG
 * \addtogroup TasmanianRefinement Refinement related classes, structures, and methods
 *
 * \par Adaptive Refinement
 * The most efficient sparse grids approximation strategies are adaptive,
 * i.e., the grid is constructed to best capture the behavior of the specific
 * target model. Adapting a grid requires analysis of the existing loaded
 * model values followed by a decision on how to best expand the grid.
 * Sometimes the complexity of the refinement algorithms requires additional
 * data storage and methods implemented outside of the main grid classes.
 * \endinternal
 */

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianRefinement
 * \brief Holds the pair of point index and model value, the struct is used in a \b std::forward_list.
 *
 * Dynamic construction allows the data for the points to be given as input in any arbitrary order.
 * However, points cannot be added to any grid unless some conditions are satisfied,
 * e.g., the resulting set must be connected or lower complete.
 * The nodes that cannot be included yet are stored into a \b std::forward_list of \b NodeData structs,
 * the point is kept in a multi-index integer format and the value into a single vector.
 * \endinternal
 */
struct NodeData{
    //! \brief The multi-index of the point.
    std::vector<int> point;
    //! \brief The values of the model outputs at the point.
    std::vector<double> value;
};

/*!
 * \internal
 * \ingroup TasmanianRefinement
 * \brief Holds the description of a single tensor candidate for inclusion into the grid.
 *
 * Some grid, e.g., Global Grids, cannot include a single point at a time, but only enough data to form a tensor.
 * The candidate tensors are stored in a \b \b std::forward_list of \b TensorData until all the points
 * associated with the tensor are provided.
 * \endinternal
 */
struct TensorData{
    //! \brief Multi-index of the tensor.
    std::vector<int> tensor;
    //! \brief The points associated with the surplus operator, the tensor can be included only when all points are available.
    MultiIndexSet points;
    //! \brief For each point \b false if the model data is not yet available for the point, true otherwise.
    std::vector<bool> loaded;
    //! \brief The weight indicates the relative importance of the tensor.
    double weight;
};

/*!
 * \internal
 * \ingroup TasmanianRefinement
 * \brief Takes a list and creates a vector of references in reversed order, needed for I/O so that read can be followed by push_front.
 * \endinternal
 */
template <class T>
std::vector<const T*> makeReverseReferenceVector(const std::forward_list<T> &list){
    size_t num_entries = (size_t) std::distance(list.begin(), list.end());
    std::vector<const T*> refs(num_entries);
    auto p = list.begin();
    auto r = refs.rbegin();
    while(p != list.end()) *r++ = &*p++;
    return refs;
}

/*!
 * \internal
 * \ingroup TasmanianRefinement
 * \brief Writes a NodeData std::forward_list to a file using either binary or ascii format.
 * \endinternal
 */
template<bool useAscii>
void writeNodeDataList(const std::forward_list<NodeData> &data, std::ostream &os){
    if (useAscii){ os << std::scientific; os.precision(17); }

    auto data_refs = makeReverseReferenceVector(data);

    IO::writeNumbers<useAscii, IO::pad_line>(os, (int) data_refs.size());
    for(auto d : data_refs){
        IO::writeVector<useAscii, IO::pad_rspace>(d->point, os);
        IO::writeVector<useAscii, IO::pad_line>(d->value, os);
    }
}

/*!
 * \internal
 * \ingroup TasmanianRefinement
 * \brief Reads a NodeData std::forward_list from a file using either binary or ascii format.
 * \endinternal
 */
template<bool useAscii>
std::forward_list<NodeData> readNodeDataList(size_t num_dimensions, size_t num_outputs, std::istream &is){
    std::forward_list<NodeData> data;
    int num_nodes = IO::readNumber<useAscii, int>(is);

    for(int i=0; i<num_nodes; i++){
        data.emplace_front(NodeData{
                            std::vector<int>(num_dimensions), // point
                            std::vector<double>(num_outputs)  // value
                           });
        IO::readVector<useAscii>(is, data.front().point);
        IO::readVector<useAscii>(is, data.front().value);
    }

    return data;
}

/*!
 * \internal
 * \ingroup TasmanianRefinement
 * \brief Helper class that stores data from dynamic construction of a Global grid.
 *
 * The class stores candidate tensors with corresponding weights
 * as well as the model data provided by the user.
 * When enough data has been computed to complete a tensor,
 * the tensor can be ejected with the points and model data returned in a format
 * that is easy to incorporate within the data structures of the \b GridGlobal class.
 */
class DynamicConstructorDataGlobal{
public:
    //! \brief The only constructor, requires that the dimension and the number of model outputs is specified.
    DynamicConstructorDataGlobal(size_t cnum_dimensions, size_t cnum_outputs);
    //! \brief Default destructor, release all used memory and erase all stored data.
    ~DynamicConstructorDataGlobal();

    //! \brief Write the data to a stream using ascii or binary format.
    template<bool useAscii> void write(std::ostream &os) const;
    //! \brief Read the data from a stream using ascii or binary format.
    template<bool useAscii> void read(std::istream &is);

    //! \brief Returns the maximum index of any of the stored tensors.
    int getMaxTensor() const;

    //! \brief Called after read, reinitializes the points and loaded structures for the tensors.
    void reloadPoints(std::function<int(int)> getNumPoints);

    //! \brief Delete the tensors with non-negative weights, i.e., clear all but the tensors selected by the initial grid.
    void clearTesnors();

    //! \brief Get a set of all tensors with negative weight, i.e., the tensors selected for the initial run.
    MultiIndexSet getInitialTensors() const;

    //! \brief Add a new tensor to the candidates, with the given \b weight and using \b getNumPoints() to define the associated nodes.
    void addTensor(const int *tensor, std::function<int(int)> getNumPoints, double weight);

    //! \brief Get the node indexes of the points associated with the candidate tensors, the order is the same as the tensors sorted by ascending weight.
    void getNodesIndexes(std::vector<int> &inodes);

    //! \brief Add a new data point with the index and the value, returns \b true if there is enough data to complete a tensor.
    bool addNewNode(const std::vector<int> &point, const std::vector<double> &value); // returns whether a tensor is complete

    //! \brief Returns a new set of tensors, points and values that can be added to the current tensors.
    void ejectCompleteTensor(MultiIndexSet const &current_tensors, MultiIndexSet &new_tensors, MultiIndexSet &new_points, StorageSet &vals);

private:
    size_t num_dimensions, num_outputs;
    std::forward_list<NodeData> data;
    std::forward_list<TensorData> tensors;
};

/*!
 * \internal
 * \brief Holds a std::forward_list of pairs of points indexes and values, and a MultiIndexSet of initial nodes.
 * \ingroup TasmanianRefinement
 *
 * Used for Sequence and Local Polynomial grids, where points can be added to the grid only when they form
 * a lower complete set (sequence) or a connected graph (local polynomial).
 * However, in order to facilitate parallelism, significantly large number of candidate points should be considered at any time.
 * A large initial grid will allow parallelism, but nodes may be computed in any arbitrary order
 * hence data has to be stored temporarily and wait until it can be incorporated into the grid.
 * \endinternal
 */
struct SimpleConstructData{
    //! \brief A list of pair indicating point and model output values.
    std::forward_list<NodeData> data;
    //! \brief Keeps track of the initial point set, so those can be computed first.
    MultiIndexSet initial_points;
    //! \brief Save to a file in either ascii or binary format.
    template<bool useAscii>
    void write(std::ostream &os) const{
        initial_points.write<useAscii>(os);
        writeNodeDataList<useAscii>(data, os);
    }
    //! \brief Remove \b points from the \b data and return a vector of the values in the \b points order.
    std::vector<double> extractValues(MultiIndexSet const &points){
        size_t num_outputs = data.front().value.size();
        int num_points = points.getNumIndexes();
        Data2D<double> result(num_outputs, num_points);
        auto p = data.before_begin();
        auto d = data.begin();
        while(d != data.end()){
            int slot = points.getSlot(d->point);
            if (slot != -1){ // found a point
                std::copy_n(d->value.begin(), num_outputs, result.getStrip(slot));
                data.erase_after(p);
                d = p;
                d++;
            }else{
                p++;
                d++;
            }
        }
        return result.getVector();
    }
};

/*!
 * \internal
 * \ingroup TasmanianRefinement
 * \brief Creates a unique_ptr holding a \b SimpleConstructData that has be read form the stream.
 *
 * \endinternal
 */
template<bool useAscii>
std::unique_ptr<SimpleConstructData> readSimpleConstructionData(size_t num_dimensions, size_t num_outputs, std::istream &is){
    std::unique_ptr<SimpleConstructData> dynamic_values = std::unique_ptr<SimpleConstructData>(new SimpleConstructData);
    dynamic_values->initial_points.read<useAscii>(is);
    dynamic_values->data = readNodeDataList<useAscii>(num_dimensions, num_outputs, is);
    return dynamic_values;
}

}

#endif
