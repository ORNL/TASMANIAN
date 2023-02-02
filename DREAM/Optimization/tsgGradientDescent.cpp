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

OptimizationStatus GradientDescent(const ObjectiveFunctionSingle &func, const GradientFunctionSingle &grad,
                                   const ProjectionFunctionSingle &proj, const double increase_coeff,
                                   const double decrease_coeff, const int max_iterations, const double tolerance,
                                   GradientDescentState &state) {

    OptimizationStatus status{0, tolerance + 1.0}; // {performed_iterations, residual}
    size_t num_dimensions = state.getNumDimensions();
    std::vector<double> x0 = state.x;
    std::vector<double> gx0(num_dimensions), gx(num_dimensions), z0(num_dimensions), xStep(num_dimensions);
    double fx(func(state.x)), fx0, fxStep;
    grad(state.x, gx);

    state.adaptive_stepsize /= increase_coeff; // Offset the first iteration.
    while(status.residual > tolerance and status.performed_iterations < max_iterations) {
        // Iteration swaps must be performed first (instead of last) to avoid reverting x to x0 after work has been completed.
        std::swap(x0, state.x);
        std::swap(fx0, fx);
        std::swap(gx0, gx);
        // Optimistic stepsize update (see γ_d in the referenced paper above).
        state.adaptive_stepsize *= increase_coeff;
        // Find the next stepsize/candidate point by making sure it satisfies the well-known descent inequality (see γ_u in the
        // reference paper associated with this function).
        double lhs(0), rhs(0);
        do {
            if (status.performed_iterations >= max_iterations) return status;
            lhs = 0;
            rhs = 0;
            for (size_t j=0; j<num_dimensions; j++)
                z0[j] = x0[j] - gx0[j] * state.adaptive_stepsize;
            proj(z0, xStep);
            fxStep = func(xStep);
            lhs += fxStep - fx0;
            for (size_t j=0; j<num_dimensions; j++) {
                double delta = xStep[j] - x0[j];
                rhs += delta * delta / (2.0 * state.adaptive_stepsize);
                lhs -= gx0[j] * delta;
            }
            state.adaptive_stepsize /= decrease_coeff;
            status.performed_iterations++;
        } while (lhs > rhs + TasGrid::Maths::num_tol);
        std::swap(xStep, state.x);
        std::swap(fxStep, fx);
        state.adaptive_stepsize *= decrease_coeff; // Offset the do-while loop.
        // Compute residual := ||(x0-x) / stepsize + gx - gx0||_2.
        grad(state.x, gx);
        status.residual = computeStationarityResidual(state.x, x0, gx, gx0, state.adaptive_stepsize);
    }
    return status;
}

OptimizationStatus GradientDescent(const ObjectiveFunctionSingle &func, const GradientFunctionSingle &grad,
                                   const double increase_coeff, const double decrease_coeff, const int max_iterations,
                                   const double tolerance, GradientDescentState &state) {
    // Wrapper to the proximal version with projection function == identity function.
    return GradientDescent(func, grad, identity, increase_coeff, decrease_coeff, max_iterations, tolerance, state);
}

OptimizationStatus GradientDescent(const GradientFunctionSingle &grad, const double stepsize, const int max_iterations,
                                   const double tolerance, std::vector<double> &state) {

    OptimizationStatus status{0, tolerance + 1.0}; // {performed_iterations, residual}
    size_t num_dimensions = state.size();
    std::vector<double> gx(num_dimensions);
    grad(state, gx);

    while (status.residual > tolerance and status.performed_iterations < max_iterations) {
        for (size_t j=0; j<num_dimensions; j++)
            state[j] -= gx[j] * stepsize;
        status.performed_iterations++;
        // Compute residual := ||grad(x)||_2.
        grad(state, gx);
        status.residual = 0.0;
        for (size_t j=0; j<num_dimensions; j++)
            status.residual += gx[j] * gx[j];
        status.residual = std::sqrt(status.residual);
    }
    return status;
}

}

#endif
