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

#ifndef __TSG_ONE_DIMENSIONAL_WRAPPER_HPP
#define __TSG_ONE_DIMENSIONAL_WRAPPER_HPP

#include "tsgSequenceOptimizer.hpp"
#include "tsgHardCodedTabulatedRules.hpp"

//! \internal
//! \file tsgOneDimensionalWrapper.hpp
//! \brief Cache for one dimensional rules.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianCoreOneDimensional
//!
//! A class to cache one dimensional rules, nodes, weight, etc.

namespace TasGrid{

//! \brief A class to cache one dimensional rules, nodes, weight, meta-data, etc.
//! \ingroup TasmanianCoreOneDimensional
class OneDimensionalWrapper{
public:
    //! \brief Default constructor, create an emptry cache.
    OneDimensionalWrapper();
    //! \brief Default destructor.
    ~OneDimensionalWrapper();

    //! \brief Load the cache.

    //! Load \b max_level of cache for a rule given by \b crule. If the rule is not
    //! custom tabulated, then \b custom is not used (could be empty).
    //! Similarly, \b alpha and \b beta are used only by rules that use the corresponding
    //! transformation parameters.
    void load(const CustomTabulated &custom, int max_level, TypeOneDRule crule, double alpha, double beta);

    //! \brief Overload that skips the custom rule altogether, more convenient in Fourier grids.
    void load(int max_level, TypeOneDRule crule, double alpha, double beta){
        load(CustomTabulated(), max_level, crule, alpha, beta);
    }

    //! \brief Get the number of points for the \b level, can only querry loaded levels.
    int getNumPoints(int level) const;
    //! \brief Get the global index of point \b j on \b level (used only by non-nested rules).
    int getPointIndex(int level, int j) const;

    //! \brief Get number of loaded nodes.
    int getNumNodes() const{ return (int) unique.size(); }
    //! \brief Get the canonical coordinate of the node with global index \b j.
    double getNode(int j) const{ return unique[j]; }
    //! \brief Get the quadrature weight of the \b j-th node on the \b level (for non-nested rules, using index local to the level)
    double getWeight(int level, int j) const;

    //! \brief Get an array of all nodes associated with the \b level (the array is ordered by local index).
    const double* getNodes(int level) const;
    //! \brief Get an array of the Lagrange coefficients associated with the \b level, the order matches \b getNodes().
    const double* getCoefficients(int level) const;

    //! \brief Get a reference to the offsets of all point counts used by the \b CacheLagrange class.
    const std::vector<int>& getPointsCount() const{ return pntr; }

    //! \brief Get the loaded rule.
    TypeOneDRule getType() const;
    //! \brief Get the number of cached levels.
    int getNumLevels() const;

    //! \brief Return reference to all unique nodes.
    const std::vector<double>& getUnique() const{ return unique; }
    //! \brief Return all nodes concatenated per level.
    std::vector<double> getAllNodes() const{ return Utils::mergeVectors(nodes); }
    //! \brief Return all coefficients concatenated per level.
    std::vector<double> getAllCoeff() const{ return Utils::mergeVectors(coeff); }
    //! \brief Return reference to the number of nodes per level.
    const std::vector<int>& getNumNodesPerLevel() const{ return num_points; }
    //! \brief Return cumulative points per level.
    std::vector<int> getOffsetNodesPerLevel() const{
        std::vector<int> offsets(num_points.size());
        offsets[0] = 0;
        for(size_t i=1; i<num_points.size(); i++)
            offsets[i] = offsets[i-1] + num_points[i-1];
        return offsets;
    }

private:
    bool isNonNested;
    int num_levels;
    TypeOneDRule rule;

    std::vector<int> num_points; // give the number of points per level
    std::vector<int> pntr; // gives the cumulative offset of each level
    std::vector<int> indx; // gives a list of the points associated with each level

    std::vector<std::vector<double>> weights; // contains the weight associated with each level
    std::vector<std::vector<double>> nodes; // contains all nodes for each level
    std::vector<double> unique; // contains the x-coordinate of each sample point

    std::vector<std::vector<double>> coeff; // the coefficients of the Lagrange
};

}

#endif
