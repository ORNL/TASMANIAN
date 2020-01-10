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

#ifndef __TASMANIAN_SPARSE_HARDCODED_RULES_HPP
#define __TASMANIAN_SPARSE_HARDCODED_RULES_HPP

/*!
 * \internal
 * \file tsgHardCodedTabulatedRules.hpp
 * \brief Hard-coded nodes and weights.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSets
 *
 * Some rules are hard to compute on the fly, thus points and weights are hard-coded.
 * \endinternal
 */

#include "tsgCoreOneDimensional.hpp"

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianCoreOneDimensional
 * \brief Rule with hard-corded tabulated points and weights.
 *
 * The Gauss-Patterson rule combines nested points with the optimality of the Gauss quadratures.
 * While not as powerful as Gauss-Legendre in one and tow dimensions, the rule wins in
 * polynomial space of exactness with respect to integration for 3 or more dimensions.
 * However, the points and weights are very hard to compute and most algorithms are very
 * ill-conditioned. Thus, it is beneficial to have the points for the first 9 levels hard-coded.
 *
 * Note that Gauss-Legendre rule combined with a full tensor grid gives highest exactness per
 * number of points with respect to integration in one and two dimensions.
 * \endinternal
 */
class TableGaussPatterson{
public:
    //! \brief Constructor, loads the nodes into the internal data structures.
    TableGaussPatterson();
    //! \brief Destrutor, cleans all memory.
    ~TableGaussPatterson() = default;

    //! \brief Return the number of hard-coded levels.
    static int getNumLevels(){ return 9; }
    //! \brief Returns the nodes for the \b level, note that the nodes are nested.
    std::vector<double> getNodes(int level) const;
    //! \brief Return the quadrature weight for \b level and given \b point.
    double getWeight(int level, int point) const;

protected:
    //! \brief Load the nodes into the local data-strutures.
    void loadNodes();
    //! \brief Load the weights into the local data-structures.
    void loadWeights();

private:
    std::vector<double> nodes; // contains the x-coordinate of each sample point
    std::vector<double> weights; // contains the weight associated with each level

    std::vector<int> weights_offsets;
};

}

#endif
