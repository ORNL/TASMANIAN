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

#ifndef __TASMANIAN_GRADIENT_DESCENT_CPP
#define __TASMANIAN_GRADIENT_DESCENT_CPP

#include "tsgGradientDescent.hpp"

namespace TasOptimization {

GradientDescentState::GradientDescentState(const std::vector<double> &x, const double stepsize) :
        num_dimensions((int) x.size()), stepsize(stepsize), candidate(x) {};

void GradientDescent(const ObjectiveFunction &f, const GradientFunction &g, const ProjectionFunction &proj, const int num_iterations,
                     GradientDescentState &state, const std::vector<double> &line_search_coeffs) {

    if (line_search_coeffs.size() != 0 and line_search_coeffs.size() != 2)
        throw std::runtime_error("ERROR: in GradientDescent(), expects line_search_coeffs.size() == 2 if non-empty");

    int num_dimensions = state.num_dimensions;
    std::vector<double> &candidate = state.getCandidateRef();

    if (line_search_coeffs.empty()) {
        // Constant stepsize scheme (does not have convergence guarantees).
        std::vector<double> gradient(num_dimensions);
        for (int i=0; i<num_iterations; i++) {
            g(candidate, gradient);
            for (int j=0; j<num_dimensions; j++)
                candidate[j] -= gradient[j] * state.stepsize;
            proj(candidate, candidate);
        }
        state.setCandidate(candidate);
    } else {
        // Variable stepsize scheme (has convergence guarantees).
        std::vector<double> current_gradient(num_dimensions), current_fval(1), next_candidate(num_dimensions), next_fval(1);
        double stepsize(state.stepsize), current_iteration(0), lhs(0), rhs(0);
        while(current_iteration < num_iterations) {
            // Find the next stepsize/candidate point by making sure it satisfies the well-known descent inequality (see
            // γ_u in the referenced paper above).
            f(candidate, current_fval);
            g(candidate, current_gradient);
            do {
                if (current_iteration >= num_iterations) return;
                lhs = 0;
                rhs = 0;
                for (int j=0; j<num_dimensions; j++)
                    next_candidate[j] = candidate[j] - current_gradient[j] * stepsize;
                proj(next_candidate, next_candidate);
                f(next_candidate, next_fval);
                lhs += next_fval[0] - current_fval[0];
                for (int j=0; j<num_dimensions; j++) {
                    double delta = next_candidate[j] - candidate[j];
                    rhs += delta * delta / (2.0 * stepsize);
                    lhs -= current_gradient[j] * delta;
                }
                stepsize /= line_search_coeffs[0];
                state.setStepsize(stepsize);
                current_iteration++;
            } while (lhs > rhs + TasGrid::Maths::num_tol);
            // Update the key iterates manually for the next loop.
            std::swap(candidate, next_candidate);
            std::swap(current_fval, next_fval);
            state.setCandidate(candidate);
            // Optimistic stepsize update (see γ_d in the referenced paper above).
            stepsize *= line_search_coeffs[1] * line_search_coeffs[0];
            state.setStepsize(stepsize);
        }
    }

}

}

#endif
