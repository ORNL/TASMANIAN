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

#ifndef __TSG_CACHE_LAGRANGE_HPP
#define __TSG_CACHE_LAGRANGE_HPP

#include "tsgEnumerates.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgOneDimensionalWrapper.hpp"

namespace TasGrid{

template <typename T>
class CacheLagrange{
public:
    CacheLagrange(int num_dimensions, const int max_levels[], const OneDimensionalWrapper *crule, const double x[]) : rule(crule){
        cache = new T*[num_dimensions];
        int full_cache = rule->getPointsCount(max_levels[0] + 1);
        for(int j=1; j<num_dimensions; j++){
            full_cache += rule->getPointsCount(max_levels[j]+1);
        }
        cache[0] = new T[full_cache];
        for(int j=0; j<num_dimensions-1; j++){
            cache[j+1] = &(cache[j][rule->getPointsCount(max_levels[j]+1)]);
        }

        for(int dim=0; dim<num_dimensions; dim++){
            for(int level=0; level <= max_levels[dim]; level++){
                const double *nodes = rule->getNodes(level);
                const double *coeff = rule->getCoefficients(level);
                int num_points = rule->getNumPoints(level);

                T *c = &(cache[dim][rule->getPointsCount(level)]);
                c[0] = 1.0;
                for(int j=0; j<num_points-1; j++){
                    c[j+1] = (x[dim] - nodes[j]) * c[j];
                }
                T w = (rule->getType() == rule_clenshawcurtis0) ? (x[dim] - 1.0) * (x[dim] + 1.0) : 1.0;
                c[num_points-1] *= w * coeff[num_points-1];
                for(int j=num_points-2; j>=0; j--){
                    w *= (x[dim] - nodes[j+1]);
                    c[j] *= w * coeff[j];
                }
            }
        }
    }
    ~CacheLagrange(){
        if (cache != 0){
            delete[] cache[0];
            delete[] cache;
        }
        rule = 0;
    }

    T getLagrange(int dimension, int level, int local) const{
        return cache[dimension][rule->getPointsCount(level) + local];
    }

private:
    T **cache;
    const OneDimensionalWrapper *rule;
};


}

#endif
