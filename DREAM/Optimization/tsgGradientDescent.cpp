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
                                   const double decrease_coeff, const int num_iterations, const double tolerance,
                                   GradientDescentState &state) {

    OptimizationStatus status;
    status.performed_iterations = 0;
    status.residual = std::numeric_limits<double>::max();

    int num_dimensions = state.num_dimensions;
    std::vector<double> &x = state.getXRef();
    std::vector<double> x0(x), gx0(num_dimensions), gx(num_dimensions), z0(num_dimensions), xStep(num_dimensions);
    double fx(func(x)), fx0, fxStep;
    grad(x, gx);

    while(status.performed_iterations < num_iterations) {
        // Iteration swaps must be performed first (instead of last) to avoid reverting x to x0 after work has been completed.
        std::swap(x0, x);
        std::swap(fx0, fx);
        std::swap(gx0, gx);
        // Find the next stepsize/candidate point by making sure it satisfies the well-known descent inequality (see γ_u in the
        // reference paper associated with this function).
        double lhs = 0;
        double rhs = 0;
        do {
            if (status.performed_iterations >= num_iterations) break;
            for (int j=0; j<num_dimensions; j++)
                z0[j] = x0[j] - gx0[j] * state.adaptive_stepsize;
            proj(z0, xStep);
            fxStep = func(xStep);
            lhs += fxStep - fx0;
            for (int j=0; j<num_dimensions; j++) {
                double delta = xStep[j] - x0[j];
                rhs += delta * delta / (2.0 * state.adaptive_stepsize);
                lhs -= gx0[j] * delta;
            }
            state.adaptive_stepsize /= decrease_coeff;
            status.performed_iterations++;
        } while (lhs > rhs + TasGrid::Maths::num_tol);
        std::swap(xStep, x);
        std::swap(fxStep, fx);
        state.adaptive_stepsize *= decrease_coeff; // offset the do-while loop.
        // Check for approximate stationarity.
        grad(x, gx);
        status.residual = compute_stationarity_residual(x, x0, gx, gx0, state.adaptive_stepsize);
        if (status.residual <= tolerance) break;
        // Optimistic stepsize update (see γ_d in the referenced paper above).
        state.adaptive_stepsize *= increase_coeff;
    }
    return status;
    
}

OptimizationStatus GradientDescent(const ObjectiveFunctionSingle &func, const GradientFunctionSingle &grad,
                                   const double increase_coeff, const double decrease_coeff, const int num_iterations,
                                   const double tolerance, GradientDescentState &state) {
    // Wrapper to the proximal version with projection function == identity function.
    return GradientDescent(func, grad, identity, increase_coeff, decrease_coeff, num_iterations, tolerance, state);
};

OptimizationStatus GradientDescent(const GradientFunctionSingle &grad, const double stepsize, const int num_iterations,
                                   const double tolerance, std::vector<double> &state) {
    OptimizationStatus status;
    status.performed_iterations = 0;
    status.residual = std::numeric_limits<double>::max();
    std::vector<double> x0(state), gx0(state.size()), gx(state.size());
    grad(state, gx);

    for (int i=0; i<num_iterations; i++) {
        std::swap(state, x0);
        std::swap(gx, gx0);
        for (size_t j=0; j<state.size(); j++)
            state[j] -= gx0[j] * stepsize;
        status.performed_iterations++;
        grad(state, gx);
        // Computes the stationarity residual as ||grad(x)||_2.
        status.residual = 0.0;
        for (size_t j=0; j<state.size(); j++)
            status.residual += gx[j] * gx[j];
        status.residual = std::sqrt(status.residual);
        if (status.residual <= tolerance) break;
    }
    return status;
}

}

#endif
