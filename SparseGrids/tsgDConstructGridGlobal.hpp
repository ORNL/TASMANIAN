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

namespace TasGrid{

struct NodeData{
    std::vector<int> point;
    std::vector<double> value;
};

struct TensorData{
    std::vector<int> tensor;
    MultiIndexSet points;
    std::vector<bool> loaded;
    double weight;
};

class DynamicConstructorDataGlobal{
public:
    DynamicConstructorDataGlobal(int cnum_dimensions, int cnum_outputs);
    ~DynamicConstructorDataGlobal();

    //! \brief Delete the tensors with non-negative weights, i.e., clear all by the tensors selected by the initial grid.
    void clearTesnors();

    //! \brief Get a set of all tensors with negative weight, i.e., the tensors selected for the initial run.
    void getInitialTensors(MultiIndexSet &set) const;

    //! \brief Add a new tensor to the candidates, with the given \b weight and using \b getNumPoints() to define the associated nodes.
    void addTensor(const int *tensor, std::function<int(int)> getNumPoints, double weight);

    //! \brief Get the node indexes of the points associated with the candidate tensors, the order matches the weights of the tensors.
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
