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

GradientDescentState::GradientDescentState(const std::vector<double> &x0, const double lambda0) :
        num_dimensions((int) x0.size()), adaptive_stepsize(lambda0), x(x0) {};

OptimizationStatus GradientDescent(const ObjectiveFunctionSingle &f, const GradientFunctionSingle &g, const ProjectionFunctionSingle &proj,
                                   GradientDescentState &state, const std::vector<double> &line_search_coeffs, const int num_iterations,
                                   const double tolerance) {

    if (line_search_coeffs.size() != 2)
        throw std::runtime_error("ERROR: in GradientDescent(), expects line_search_coeffs.size() == 2");

    OptimizationStatus status;
    status.num_iterations = 0;
    status.stationarity_residual = std::numeric_limits<double>::max();

    int num_dimensions = state.num_dimensions;
    std::vector<double> &x = state.getXRef();
    std::vector<double> x_prev(x), g_at_x_prev(g(x_prev)), g_at_x(num_dimensions), next_x(num_dimensions);
    double f_at_x_prev(f(x_prev)), f_at_x(0);

    while(status.num_iterations < num_iterations) {
        // Find the next stepsize/candidate point by making sure it satisfies the well-known descent inequality (see γ_u in the
        // reference paper associated with this function).
        double lhs = 0;
        double rhs = 0;
        do {
            if (status.num_iterations >= num_iterations)
                break;
            for (int j=0; j<num_dimensions; j++)
                x[j] = x_prev[j] - g_at_x_prev[j] * state.adaptive_stepsize;
            x = proj(x);
            f_at_x = f(x);
            lhs += f_at_x - f_at_x_prev;
            for (int j=0; j<num_dimensions; j++) {
                double delta = x[j] - x_prev[j];
                rhs += delta * delta / (2.0 * state.adaptive_stepsize);
                lhs -= g_at_x_prev[j] * delta;
            }
            state.adaptive_stepsize /= line_search_coeffs[0];
            status.num_iterations++;
        } while (lhs > rhs + TasGrid::Maths::num_tol);
        state.adaptive_stepsize *= line_search_coeffs[0]; // offset the do-while loop.
        // Check for approximate stationarity.
        g_at_x = g(x);
        status.stationarity_residual = compute_stationarity_residual(x, x_prev, g_at_x, g_at_x_prev, state.adaptive_stepsize);
        if (status.stationarity_residual <= tolerance)
            break;
        // Update the key iterates manually for the next loop.
        std::swap(x_prev, x);
        std::swap(f_at_x_prev, f_at_x);
        std::swap(g_at_x_prev, g_at_x);
        // Optimistic stepsize update (see γ_d in the referenced paper above).
        state.adaptive_stepsize *= line_search_coeffs[1];
    }
    return status;
    
}

OptimizationStatus GradientDescent(const ObjectiveFunctionSingle &f, const GradientFunctionSingle &g, GradientDescentState &state,
                                   const std::vector<double> &line_search_coeffs,  const int num_iterations, const double tolerance) {
    // Wrapper to the proximal version with projection function == identity function.
    return GradientDescent(f, g, identity, state, line_search_coeffs, num_iterations, tolerance);
};

OptimizationStatus GradientDescent(const GradientFunctionSingle &g, std::vector<double> &x, const double stepsize, const int num_iterations,
                                   const double tolerance) {
    OptimizationStatus status;
    status.num_iterations = 0;
    status.stationarity_residual = std::numeric_limits<double>::max();
    std::vector<double> x_prev(x), g_at_x_prev(g(x_prev)), g_at_x(x.size());
    for (int i=0; i<num_iterations; i++) {
        for (size_t j=0; j<x.size(); j++)
            x[j] -= g_at_x_prev[j] * stepsize;
        g_at_x = g(x);
        status.num_iterations++;
        status.stationarity_residual = compute_stationarity_residual(x, x_prev, g_at_x, g_at_x_prev, stepsize);
        if (status.stationarity_residual <= tolerance)
            break;
        std::swap(x, x_prev);
        std::swap(g_at_x, g_at_x_prev);
    }
    return status;
}

}

#endif
