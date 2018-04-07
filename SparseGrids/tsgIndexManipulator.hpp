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

#ifndef __TSG_INDEX_MANIPULATOR_HPP
#define __TSG_INDEX_MANIPULATOR_HPP

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgOneDimensionalWrapper.hpp"
#include "tsgRuleLocalPolynomial.hpp"

namespace TasGrid{

class IndexManipulator{
public:
    IndexManipulator(int cnum_dimensions, const CustomTabulated* custom = 0);
    ~IndexManipulator();

    int getIndexWeight(const int index[], TypeDepth type, const int weights[], TypeOneDRule rule) const;

    IndexSet* selectTensors(int offset, TypeDepth type, const int *anisotropic_weights, TypeOneDRule rule) const;
    IndexSet* selectTensors(const IndexSet *target_space, bool integration, TypeOneDRule rule) const;
    IndexSet* getLowerCompletion(const IndexSet *set) const;

    int* computeLevels(const IndexSet* set) const; // returns a vector of its corresponding to the sum of entries of set

    int* makeTensorWeights(const IndexSet* set) const;
    IndexSet* nonzeroSubset(const IndexSet* set, const int weights[]) const; // returns a subset corresponding to the non-zeros weights

    UnsortedIndexSet* tensorGenericPoints(const int levels[], const OneDimensionalWrapper *rule) const;
    IndexSet* generateGenericPoints(const IndexSet *tensors, const OneDimensionalWrapper *rule) const;
    int* referenceGenericPoints(const int levels[], const OneDimensionalWrapper *rule, const IndexSet *points) const;

    IndexSet* removeIndexesByLimit(IndexSet *set, const int limits[]) const;
    UnsortedIndexSet* removeIndexesByLimit(UnsortedIndexSet *set, const int limits[]) const;

    int* computeDAGup(const IndexSet *set) const;

    IndexSet* tensorNestedPoints(const int levels[], const OneDimensionalWrapper *rule) const;
    IndexSet* generateNestedPoints(const IndexSet *tensors, const OneDimensionalWrapper *rule) const;
    int* referenceNestedPoints(const int levels[], const OneDimensionalWrapper *rule, const IndexSet *points) const;

    IndexSet* getPolynomialSpace(const IndexSet *tensors, TypeOneDRule rule, bool iexact) const;

    int getMinChildLevel(const IndexSet *set, TypeDepth type, const int weights[], TypeOneDRule rule);

    IndexSet* selectFlaggedChildren(const IndexSet *set, const bool flagged[], const int *level_limits = 0) const;

    void getMaxLevels(const IndexSet *set, int max_levels[], int &total_max) const;

    // use by Local Grids
    UnsortedIndexSet* getToalDegreeDeltas(int level) const;
    IndexSet* generatePointsFromDeltas(const UnsortedIndexSet* deltas, const BaseRuleLocalPolynomial *rule) const;

    int* computeDAGupLocal(const IndexSet *set, const BaseRuleLocalPolynomial *rule) const;

protected:
    int getLevel(const int index[], const int weights[]) const;
    int getCurved(const int index[], const int weights[]) const;
    int getIPTotal(const int index[], const int weights[], TypeOneDRule rule) const;
    int getIPCurved(const int index[], const int weights[], TypeOneDRule rule) const;
    int getQPTotal(const int index[], const int weights[], TypeOneDRule rule) const;
    int getQPCurved(const int index[], const int weights[], TypeOneDRule rule) const;
    int getHyperbolic(const int index[], const int weights[]) const;
    int getIPHyperbolic(const int index[], const int weights[], TypeOneDRule rule) const;
    int getQPHyperbolic(const int index[], const int weights[], TypeOneDRule rule) const;

    int* getProperWeights(TypeDepth type, const int *anisotropic_weights) const;

private:
    int num_dimensions;

    OneDimensionalMeta *meta;
};

}

#endif
