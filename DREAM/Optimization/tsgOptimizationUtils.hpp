/*
 * Copyright (c) 2022, Miroslav Stoyanov & Weiwei Kong
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
 * conditions are met:
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
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND
 * IMPLIED. THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
 * THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL
 * ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE. THE USER ASSUMES
 * RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING
 * FROM OR ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TASMANIAN_OPTIM_ENUMERATES_HPP
#define __TASMANIAN_OPTIM_ENUMERATES_HPP

#include "TasmanianDREAM.hpp"
#include <assert.h>

namespace TasOptimization {

inline void checkVarSize(std::string method_name, std::string var_name, int var_size, int exp_size) {
    if (var_size != exp_size) {
        throw std::runtime_error("Size of " + var_name + " (" + std::to_string(var_size) + ") in the function " + method_name +
                                 "() is not equal to its expected value of (" + std::to_string(exp_size) + ")");
    }
}

using ObjectiveFunctionSingle = std::function<double(const std::vector<double> &x)>;

using ObjectiveFunction = std::function<void(const std::vector<double> &x_batch, std::vector<double> &fval_batch)>;

inline ObjectiveFunction makeObjectiveFunction(int ndim, ObjectiveFunctionSingle f_single) {
    return [=](const std::vector<double> &x_values, std::vector<double> &fval_values)->void {
        int num_points = x_values.size() / ndim;
        std::vector<double> x(ndim);
        for (int i=0; i<num_points; i++) {
            std::copy(x_values.begin() + i * ndim, x_values.begin() + (i + 1) * ndim, x.begin());
            fval_values[i] = f_single(x);
        }
    };
}

// If inside_batch[i] is true, then fval_batch[i] is the function value of a point in the domain; otherwise, if inside_batch[i]
// is false, then fval_batch[i] is std::numeric_limits<double>::max().
using ObjectiveFunctionConstrained = std::function<void(const std::vector<double> &x_batch, std::vector<double> &fval_batch,
                                                        std::vector<bool> &inside_batch)>;

inline ObjectiveFunctionConstrained makeObjectiveFunctionConstrained(int num_dimensions, ObjectiveFunction f,
                                                                     TasDREAM::DreamDomain inside) {
    return [=](const std::vector<double> &x_batch, std::vector<double> &fval_batch, std::vector<bool> &inside_batch)->void {
        // Collect and apply the domain information given by inside() and x_batch.
        size_t num_particles(fval_batch.size()), num_inside(0);
        std::vector<double> candidate(num_dimensions), is_inside(num_particles), inside_points;
        for (size_t i=0; i<num_particles; i++) {
            fval_batch[i] = std::numeric_limits<double>::max();
            std::copy_n(x_batch.begin() + i * num_dimensions, num_dimensions, candidate.begin());
            inside_batch[i] = inside(candidate);
            if (inside_batch[i]) {
                std::copy_n(candidate.begin(), num_dimensions, std::back_inserter(inside_points));
                num_inside++;
            }
        }
        // Evaluate f on the inside points and copy the resulting values to fval_batch.
        std::vector<double> inside_vals(num_inside);
        f(inside_points, inside_vals);
        int j = 0;
        for (size_t i=0; i<num_particles; i++) {
            if (inside_batch[i])
                fval_batch[i] = inside_vals[j++];
        }
    };
}

} // End namespace

#endif
