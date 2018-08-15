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

#include <functional>

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

    IndexSet* selectTensors(int offset, TypeDepth type, const std::vector<int> &anisotropic_weights, TypeOneDRule rule) const;
    IndexSet* selectTensors(const IndexSet *target_space, bool integration, TypeOneDRule rule) const;
    IndexSet* getLowerCompletion(const IndexSet *iset) const;

    void computeLevels(const IndexSet* iset, std::vector<int> &level) const; // returns a vector of its corresponding to the sum of entries of set

    void makeTensorWeights(const IndexSet* iset, std::vector<int> &weights) const;
    IndexSet* nonzeroSubset(const IndexSet* iset, const std::vector<int> &weights) const; // returns a subset corresponding to the non-zeros weights

    UnsortedIndexSet* tensorGenericPoints(const int levels[], const OneDimensionalWrapper *rule) const;
    IndexSet* generateGenericPoints(const IndexSet *tensors, const OneDimensionalWrapper *rule) const;
    int* referenceGenericPoints(const int levels[], const OneDimensionalWrapper *rule, const IndexSet *points) const;

    IndexSet* removeIndexesByLimit(IndexSet *iset, const std::vector<int> &limits) const;

    void computeDAGup(const IndexSet *iset, Data2D<int> &parents) const;

    IndexSet* tensorNestedPoints(const int levels[], const OneDimensionalWrapper *rule) const;
    IndexSet* generateNestedPoints(const IndexSet *tensors, const OneDimensionalWrapper *rule) const;
    int* referenceNestedPoints(const int levels[], const OneDimensionalWrapper *rule, const IndexSet *points) const;

    IndexSet* getPolynomialSpace(const IndexSet *tensors, TypeOneDRule rule, bool iexact) const;

    int getMinChildLevel(const IndexSet *iset, TypeDepth type, const std::vector<int> &weights, TypeOneDRule rule);
    // find the minimum level of a child of iset

    IndexSet* selectFlaggedChildren(const IndexSet *iset, const std::vector<bool> &flagged, const std::vector<int> &level_limits) const;

    void getMaxLevels(const IndexSet *iset, std::vector<int> &max_levels, int &total_max) const;
    int getMaxLevel(const IndexSet *iset) const;

    // use by Local Grids
    IndexSet* generatePointsFromDeltas(const IndexSet* deltas, std::function<int(int)> getNumPoints) const;

    void computeDAGupLocal(const IndexSet *iset, const BaseRuleLocalPolynomial *rule, Data2D<int> &parents) const;

protected:
    void getProperWeights(TypeDepth type, const std::vector<int> &anisotropic_weights, std::vector<int> &weights) const;

    template<TypeDepth type>
    long long getIndexWeight(const std::vector<int> &index, const std::vector<int> &weights, TypeOneDRule rule) const{
        long long l = ((type == type_hyperbolic) || (type == type_iphyperbolic) || (type == type_qphyperbolic)) ? 1 : 0;
        if (type == type_level){
            auto witer = weights.begin();
            for(auto i : index) l += i * *witer++;
        }else if (type == type_curved){
            auto walpha = weights.begin();
            auto wbeta = weights.begin(); std::advance(wbeta, num_dimensions);
            double c = 0.0;
            for(auto i : index){
                l += i * *walpha++;
                c += log1p((double) i) * *wbeta++;
            }
            l += (long long) ceil(c);
        }else if ((type == type_iptotal) || (type == type_qptotal)){
            auto witer = weights.begin();
            for(auto i : index) l += ((i > 0) ? 1 + ((type == type_iptotal) ? meta.getIExact(i-1, rule) : meta.getQExact(i-1, rule)) : 0) * *witer++;
        }else if ((type == type_ipcurved) || (type == type_qpcurved)){
            auto walpha = weights.begin();
            auto wbeta = weights.begin(); std::advance(wbeta, num_dimensions);
            double c = 0.0;
            for(auto i : index){
                long long pex = (i > 0) ? 1 + ((type == type_ipcurved) ? meta.getIExact(i-1, rule) : meta.getQExact(i-1, rule)) : 0;
                l += pex * *walpha++;
                c += log1p((double) pex) * *wbeta++;
            }
            l += (long long) ceil(c);
        }else if (type == type_iphyperbolic){
            double nweight = (double) weights[num_dimensions];
            double s = 1.0;
            auto witer = weights.begin();
            for(auto i : index) s *= pow((double) (i+1), ((double) *witer++) / nweight);
            l = (long long) ceil(s);
        }else{
            double nweight = (double) weights[num_dimensions];
            double s = 1.0;
            auto witer = weights.begin();
            for(auto i : index) s *= pow((i>0) ? 2.0 + (double) ((type == type_iphyperbolic) ? meta.getIExact(i-1, rule) : meta.getQExact(i-1, rule)) : 1.0,
                                         ((double) *witer++) / nweight);
            l = (long long) ceil(s);
        }
        return l;
    }

    long long getIndexWeight(const std::vector<int> &index, TypeDepth type, const std::vector<int> &weights, TypeOneDRule rule) const;

private:
    int num_dimensions;

    OneDimensionalMeta meta;
};

}

#endif
