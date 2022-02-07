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

void ParticleSwarmState::setParticlePositions(const std::vector<double> &pp) {
    if (pp.size() != (size_t) num_dimensions * num_particles) {
        throw std::runtime_error("Size of particle positions (" + std::to_string(pp.size()) + ") is not equal to the " +
                                 "number of dimensions (" + std::to_string(num_dimensions) + ") times the number of " +
                                 "particles (" + std::to_string(num_particles) + ")");
    }
    particle_positions = pp;
}


void ParticleSwarmState::setParticleVelocities(const std::vector<double> &pv) {
    if (pv.size() != (size_t) num_dimensions * num_particles) {
        throw std::runtime_error("Size of particle velocities (" + std::to_string(pv.size()) + ") is not equal to the " +
                                 "number of dimensions (" + std::to_string(num_dimensions) + ") times the number of " +
                                 "particles (" + std::to_string(num_particles) + ")");
    }
    particle_velocities = pv;
}


void ParticleSwarmState::setBestParticlePositions(const std::vector<double> &bpp) {
    if (bpp.size() != (size_t) num_dimensions * num_particles) {
        throw std::runtime_error("Size of best particle positions (" + std::to_string(bpp.size()) + ") is not equal to the " +
                                 "number of dimensions (" + std::to_string(num_dimensions) + ") times the number of " +
                                 "particles (" + std::to_string(num_particles) + ")");
    }
    best_particle_positions = bpp;
}

void ParticleSwarmState::setBestPosition(const std::vector<double> &bp) {
    if (bp.size() != (size_t) num_dimensions) {
        throw std::runtime_error("Size of best position (" + std::to_string(bp.size()) + ") is not equal to the number of " +
                                 "dimensions (" + std::to_string(num_dimensions) + ")");
    }
    best_position = bp;
}

void ParticleSwarmState::addParticlesInsideBox(const int box_num_particles, const std::vector<double> &box_lower,
                                               const std::vector<double> &box_upper, const std::function<double(void)> get_random01) {
    if (box_lower.size() != (size_t) num_dimensions) {
        throw std::runtime_error("The size of the box lower bounds (" + std::to_string(box_lower.size()) + ") does not equal " +
                                 "to the number of dimensions (" + std::to_string(num_dimensions) + ")");
    }
    if (box_upper.size() != (size_t) num_dimensions) {
        throw std::runtime_error("The size of the box upper bounds (" + std::to_string(box_upper.size()) + ") does not equal " +
                                 "to the number of dimensions (" + std::to_string(num_dimensions) + ")");
    }
    for (int i=0; i<box_num_particles; i++) {
        for (int j=0; j<num_dimensions; j++) {
            double range = std::fabs(box_upper[j] - box_lower[j]);
            particle_positions.push_back(range * get_random01() + box_lower[j]);
            particle_velocities.push_back(2 * range * get_random01() - range);
        }
    }
    num_particles += box_num_particles;
}

void ParticleSwarm(ObjectiveFunction f, int max_iterations, TasDREAM::DreamDomain inside, ParticleSwarmState &state,
                   double inertia_weight, double cognitive_coeff, double social_coeff,
                   std::function<double(void)> get_random01) {

    // Initialize.
    size_t num_dimensions = (size_t) state.getNumDimensions();
    size_t num_particles = (size_t) state.getNumParticles();
    std::vector<double>& particle_positions = state.getParticlePositionsRef();
    std::vector<double>& particle_velocities = state.getParticleVelocitiesRef();
    std::vector<double>& best_particle_positions = state.getBestParticlePositionsRef();
    if (best_particle_positions.size() == 0) {
        best_particle_positions = particle_positions;
    }
    std::vector<double>& best_position = state.getBestPositionRef();
    std::vector<double> prev_fvals(num_particles), fvals(num_particles), candidate(num_dimensions);
    double best_fval = std::numeric_limits<double>::max();
    if (best_position.size() == 0) {
        f(particle_positions, fvals);
        for (size_t i=0; i<num_particles; i++) {
            std::copy(particle_positions.begin() + i * num_dimensions, particle_positions.begin() + (i + 1) * num_dimensions,
                      candidate.begin());
            if (inside(candidate) && fvals[i] < best_fval) {
                best_position = candidate;
                best_fval = fvals[i];
            }
        }
    }
    f(best_particle_positions, prev_fvals);

    // Main optimization loop.
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
                particle_positions[i * num_dimensions + j] += particle_velocities[i * num_dimensions + j];
            }
        }

        // Update the best positions.
        f(particle_positions, fvals);
        for (size_t i=0; i<num_particles; i++) {
            std::copy(particle_positions.begin() + i * num_dimensions, particle_positions.begin() + (i + 1) * num_dimensions,
                      candidate.begin());
            if (inside(candidate) && fvals[i] < prev_fvals[i]) {
                std::copy(particle_positions.begin() + i * num_dimensions, particle_positions.begin() + (i+1) * num_dimensions,
                          best_particle_positions.begin() + i * num_dimensions);
                if (fvals[i] < best_fval) {
                    best_position = candidate;
                    best_fval = fvals[i];
                }
            }
        }
        state.addIterations(1);
    } // End while
}

} // End namespace

#endif
