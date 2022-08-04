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

#ifndef __TASMANIAN_PARTICLE_SWARM_CPP
#define __TASMANIAN_PARTICLE_SWARM_CPP

#include "tsgParticleSwarm.hpp"

namespace TasOptimization {


ParticleSwarmState::ParticleSwarmState(int cnum_dimensions, int cnum_particles):
        positions_initialized(false), velocities_initialized(false), best_positions_initialized(false), cache_initialized(false),
        num_dimensions(cnum_dimensions), num_particles(cnum_particles),
        particle_positions(std::vector<double>(num_particles * num_dimensions)),
        particle_velocities(std::vector<double>(num_particles * num_dimensions)),
        best_particle_positions(std::vector<double>((num_particles + 1) * num_dimensions)),
        cache_particle_fvals(std::vector<double>(num_particles, std::numeric_limits<double>::max())),
        cache_best_particle_fvals(std::vector<double>(num_particles + 1, std::numeric_limits<double>::max())),
        cache_particle_inside(std::vector<bool>(num_particles, false)),
        cache_best_particle_inside(std::vector<bool>(num_particles + 1, false)) {};


ParticleSwarmState::ParticleSwarmState(int cnum_dimensions, std::vector<double> &&pp, std::vector<double> &&pv):
        positions_initialized(true), velocities_initialized(true), best_positions_initialized(false), cache_initialized(false),
        num_dimensions(cnum_dimensions), num_particles(pp.size() / num_dimensions),
        particle_positions(std::move(pp)), particle_velocities(std::move(pv)),
        best_particle_positions(std::vector<double>((num_particles + 1) * num_dimensions)),
        cache_particle_fvals(std::vector<double>(num_particles, std::numeric_limits<double>::max())),
        cache_best_particle_fvals(std::vector<double>(num_particles + 1, std::numeric_limits<double>::max())),
        cache_particle_inside(std::vector<bool>(num_particles, false)),
        cache_best_particle_inside(std::vector<bool>(num_particles + 1, false)) {};


void ParticleSwarmState::initializeParticlesInsideBox(const double box_lower[], const double box_upper[],
                                                      const std::function<double(void)> get_random01) {
    for (int i=0; i<num_particles * num_dimensions; i++) {
        double range = std::fabs(box_upper[i % num_dimensions] - box_lower[i % num_dimensions]);
        particle_positions[i] = range * get_random01() + box_lower[i % num_dimensions];
        particle_velocities[i] = 2 * range * get_random01() - range;
    }
    positions_initialized = true;
    velocities_initialized = true;
}

void ParticleSwarmState::initializeParticlesInsideBox(const std::vector<double> &box_lower, const std::vector<double> &box_upper,
                                                      const std::function<double(void)> get_random01) {
    checkVarSize("ParticleSwarmState::initializeParticlesInsideBox", "box lower bounds", box_lower.size(), num_dimensions);
    checkVarSize("ParticleSwarmState::initializeParticlesInsideBox", "box upper bounds", box_upper.size(), num_dimensions);
    ParticleSwarmState::initializeParticlesInsideBox(box_lower.data(), box_upper.data(), get_random01);
}

void ParticleSwarm(const ObjectiveFunction f, const TasDREAM::DreamDomain inside, const double inertia_weight,
                   const double cognitive_coeff, const double social_coeff, const int num_iterations,
                   ParticleSwarmState &state, const std::function<double(void)> get_random01) {

    // Only run the algorithm on properly initialized states.
    if (!state.positions_initialized) {
        throw std::runtime_error("Particle positions have not been initialized in the input state object");
    }
    if (!state.velocities_initialized) {
        throw std::runtime_error("Particle velocities have not been initialized in the input state object");
    }

    // Initialize helper variables and functions.
    size_t num_dimensions = (size_t) state.getNumDimensions();
    size_t num_particles = (size_t) state.getNumParticles();
    std::vector<double> &particle_positions = state.getParticlePositionsRef();
    std::vector<double> &particle_velocities = state.getParticleVelocitiesRef();
    std::vector<double> &best_particle_positions = state.getBestParticlePositionsRef();
    std::vector<double> &cache_particle_fvals = state.getCacheParticleFValsRef();
    std::vector<double> &cache_best_particle_fvals = state.getCacheBestParticleFValsRef();
    std::vector<bool> &cache_particle_inside = state.getCacheParticleInsideRef();
    std::vector<bool> &cache_best_particle_inside = state.getCacheBestParticleInsideRef();

    // Create a lambda that converts f to a constrained version that only evaluates points inside the domain. This lambda also
    // writes to a bool vector whose i-th entry is true if particle i is in the domain.
    auto f_constrained = [=](const std::vector<double> &x_batch, std::vector<double> &fval_batch, std::vector<bool> &inside_batch)->void {
        // Collect and apply the domain information given by inside() and x_batch.
        size_t num_batch(fval_batch.size()), num_inside(0);
        std::vector<double> candidate(num_dimensions), is_inside(num_batch), inside_points;
        for (size_t i=0; i<num_batch; i++) {
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
        if (num_inside > 0)
            f(inside_points, inside_vals);
        int j = 0;
        for (size_t i=0; i<num_batch; i++) {
            if (inside_batch[i])
                fval_batch[i] = inside_vals[j++];
        }
    };

    // Create a lambda that updates the best particle positions and fvals (in the cache).
    auto update = [&]()->void {
        for (size_t i=0; i<num_particles; i++) {
            bool is_smaller = (cache_particle_fvals[i] < cache_best_particle_fvals[i]) or !cache_best_particle_inside[i];
            if (cache_particle_fvals[i] and is_smaller) {
                std::copy_n(particle_positions.begin() + i * num_dimensions, num_dimensions, best_particle_positions.begin() + i * num_dimensions);
                cache_best_particle_fvals[i] = cache_particle_fvals[i];
                if (cache_best_particle_fvals[i] < cache_best_particle_fvals[num_particles]) {
                    std::copy_n(particle_positions.begin() + i * num_dimensions, num_dimensions,
                                best_particle_positions.begin() + num_particles * num_dimensions);
                    cache_best_particle_fvals[num_particles] = cache_best_particle_fvals[i];
                    cache_best_particle_inside[num_particles] = true;
                }
            }
        }
    };

    // Set up the cache and best particle positions.
    if (!state.cache_initialized) {
        f_constrained(particle_positions, cache_particle_fvals, cache_particle_inside);
        if (state.best_positions_initialized) {
            f_constrained(best_particle_positions, cache_best_particle_fvals, cache_best_particle_inside);
        }
        state.cache_initialized = true;
    }
    update();
    state.best_positions_initialized = true;

    // Main algorithm starts here.
    for(int iter=0; iter<num_iterations; iter++) {
        for (size_t i=0; i<num_particles * num_dimensions; i++) {
            particle_velocities[i] = inertia_weight * particle_velocities[i] +
                                     cognitive_coeff * get_random01() * (best_particle_positions[i] - particle_positions[i]) +
                                     social_coeff * get_random01() * (
                                         // Best overall position.
                                         best_particle_positions[num_particles * num_dimensions + i % num_dimensions] -
                                         particle_positions[i]);
        }
        for (size_t i=0; i< num_particles * num_dimensions; i++) {
            particle_positions[i] += particle_velocities[i];
        }
        f_constrained(particle_positions, cache_particle_fvals, cache_particle_inside);
        update();
    }
}


} // End namespace

#endif
