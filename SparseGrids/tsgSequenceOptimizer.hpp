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

#ifndef __TASMANIAN_SPARSE_GRID_SEQUENCE_OPTIMIZER_HPP
#define __TASMANIAN_SPARSE_GRID_SEQUENCE_OPTIMIZER_HPP

#include "tsgEnumerates.hpp"

/*!
 * \internal
 * \file tsgSequenceOptimizer.hpp
 * \brief Algorithms for computing greedy sequences.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSets
 *
 * \endinternal
 */

/*!
 * \internal
 * \ingroup TasmanianSG
 * \addtogroup TasmanianSequenceOpt Greedy Sequences
 *
 * \par Sequence One Dimensional Rules
 * Nested rules that grow with one point per level allow great flexibility in
 * constructing sparse grids; however, controlling the Lebesgue constant for
 * such sequences is a challenge. The standard way to construct a sequence is
 * to follow a greedy optimization problem, start with a few initial nodes
 * (usually 0, 1, -1) then add new nodes at the maximum of a specific functional.
 * Thus, nodes are computed one at a time and one-dimensional optimization
 * algorithms are needed.
 * \endinternal
 */

namespace TasGrid{

namespace Optimizer{

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Simple pair of numbers for the \b node and the \b value of the functional at the node.
 */
struct OptimizerResult{
    //! \brief Node where the functional was evaluated.
    double node;
    //! \brief Value of the functional.
    double value;
};

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief For the given \b rule and set of \b nodes, compute the next node using the greedy procedure.
 */
template<TypeOneDRule rule> double getNextNode(std::vector<double> const &nodes);

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Get the hard-coded pre-computed nodes.
 *
 * Some nodes are very expensive to compute, thus we store a pre-computed set
 * that contains enough nodes for most applications.
 */
std::vector<double> getPrecomputedMinLebesgueNodes();
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Get the hard-coded pre-computed nodes.
 *
 * Some nodes are very expensive to compute, thus we store a pre-computed set
 * that contains enough nodes for most applications.
 */
std::vector<double> getPrecomputedMinDeltaNodes();

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Get \b n nodes for the given sequence \b rule, either compute or use pre-computed.
 */
template<TypeOneDRule rule> std::vector<double> getGreedyNodes(int n);

}

}

#endif
