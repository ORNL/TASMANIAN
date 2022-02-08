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

#ifndef __TASMANIAN_PARTICLE_SWARM_SOLVER_CPP
#define __TASMANIAN_PARTICLE_SWARM_SOLVER_CPP

#include "tsgParticleSwarm.hpp"

namespace TasOptimization {


ParticleSwarmState::ParticleSwarmState(int ndim, int npart):
        particle_positions(std::vector<double>(npart * ndim)), particle_velocities(std::vector<double>(npart * ndim)),
        best_particle_positions(std::vector<double>(npart * ndim)), best_position(std::vector<double>(ndim)),
        num_dimensions(ndim), num_particles(npart), num_iterations(0), initialized(false) {};


ParticleSwarmState::ParticleSwarmState(int ndim, std::vector<double> &&pp, std::vector<double> &&pv):
        num_dimensions(ndim), num_iterations(0), initialized(true) {
    assert(particle_positions.size() % ndim == 0);
    int npart = particle_positions.size() / ndim;
    particle_positions = std::move(pp);
    particle_velocities = std::move(pv);
    best_particle_positions = std::vector<double>(npart * ndim);
    best_position = std::vector<double>(ndim);
};


void ParticleSwarmState::addParticlesInsideBox(const std::vector<double> &box_lower, const std::vector<double> &box_upper,
                                               const std::function<double(void)> get_random01) {
    checkVarSize("box lower bounds", box_lower.size(), num_dimensions);
    checkVarSize("box upper bounds", box_upper.size(), num_dimensions);
    for (int i=0; i<num_particles; i++) {
        for (int j=0; j<num_dimensions; j++) {
            double range = std::fabs(box_upper[j] - box_lower[j]);
            int idx = i * num_dimensions + j;
            particle_positions[idx] = range * get_random01() + box_lower[j];
            particle_velocities[idx] = 2 * range * get_random01() - range;
        }
    }
    initialized = true;
}


void ParticleSwarm(ObjectiveFunction f, int max_iterations, TasDREAM::DreamDomain inside, ParticleSwarmState &state,
                   double inertia_weight, double cognitive_coeff, double social_coeff, std::function<double(void)> get_random01) {

    // Only run the algorithm on properly initialized states.
    if (!state.initialized) {
        throw std::runtime_error("Particles have not been initialized in the input state object");
    }

    // Initialize helper variables and functions.
    size_t num_dimensions = (size_t) state.getNumDimensions();
    size_t num_particles = (size_t) state.getNumParticles();
    std::vector<double> &particle_positions = state.getParticlePositionsRef();
    std::vector<double> &particle_velocities = state.getParticleVelocitiesRef();
    std::vector<double> &best_particle_positions = state.getBestParticlePositionsRef();
    std::vector<double> &best_position = state.getBestPositionRef();
    std::vector<double> prev_best_vals(num_particles), cur_best_vals(num_particles);
    std::vector<bool> prev_is_inside(num_particles), cur_is_inside(num_particles);
    double best_fval = std::numeric_limits<double>::max();
    ObjectiveFunctionConstrained f_constrained = makeObjectiveFunctionConstrained(num_dimensions, f, inside);

    // Create a lambda that updates the best particle and swarm positions based on function values and domain information.
    auto update_best_positions = [&](const std::vector<double> prev_vals, const std::vector<double> cur_vals,
                                     const std::vector<bool> prev_domain, const std::vector<bool> cur_domain)->void {
        for (size_t i=0; i<num_particles; i++) {
            bool is_smaller = (cur_vals[i] < prev_vals[i]) || !prev_domain[i];
            if (cur_domain[i] && is_smaller) {
                std::copy(particle_positions.begin() + i * num_dimensions, particle_positions.begin() + (i + 1) * num_dimensions,
                          best_particle_positions.begin() + i * num_dimensions);
                if (cur_vals[i] < best_fval) {
                    std::copy(particle_positions.begin() + i * num_dimensions, particle_positions.begin() + (i + 1) * num_dimensions,
                              best_position.begin());
                    best_fval = cur_vals[i];
                }
            }
        }
    };

    // Initialize the best particle and swarm function values, while updating the input state.
    f_constrained(best_particle_positions, prev_best_vals, prev_is_inside);
    f_constrained(particle_positions, cur_best_vals, cur_is_inside);
    update_best_positions(prev_best_vals, cur_best_vals, prev_is_inside, cur_is_inside);
    prev_best_vals = cur_best_vals;
    prev_is_inside = cur_is_inside;

    // Main algorithm starts here.
    while(state.getNumIterations() < max_iterations) {
        for (size_t i=0; i<num_particles; i++) {
            for (size_t j=0; j<num_dimensions; j++) {
                double rp = get_random01();
                double rg = get_random01();
                int idx = i * num_dimensions + j;
                particle_velocities[idx] = inertia_weight * particle_velocities[idx] +
                                           cognitive_coeff * rp * (best_particle_positions[idx] - particle_positions[idx]) +
                                           social_coeff * rg * (best_position[j] - particle_positions[idx]);
            }
            for (size_t j=0; j<num_dimensions; j++) {
                int idx = i * num_dimensions + j;
                particle_positions[idx] += particle_velocities[idx];
            }
        }
        f_constrained(particle_positions, cur_best_vals, cur_is_inside);
        update_best_positions(prev_best_vals, cur_best_vals, prev_is_inside, cur_is_inside);
        prev_best_vals = cur_best_vals;
        prev_is_inside = cur_is_inside;
        state.addIterations(1);
    }
}


} // End namespace

#endif
