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

#include "tsgParticleSwarmSolver.hpp"

namespace TasOptimization {

ParticleSwarmSolver::ParticleSwarmSolver(int num_particles, int seed) :
        inertia_weight(0.5), cognitive_coeff(2), social_coeff(2), num_particles(num_particles), seed(seed) {}

ParticleSwarmSolver::ParticleSwarmSolver(double inertia_weight, double cognitive_coeff, double social_coeff, int num_particles, int seed) :
        inertia_weight(inertia_weight), cognitive_coeff(cognitive_coeff), social_coeff(social_coeff), num_particles(num_particles), seed(seed) {}

void ParticleSwarmSolver::optimize() {

    // TODO: Add input checks (e.g., lower, upper).
    // TODO: Add runtime updates.

    // Initialize the particles, their positions & velocities, and the swarm position.
    std::vector<double> lower = getLowerBounds();
    std::vector<double> upper = getUpperBounds();
    size_t num_dimensions = getNumDimensions();
    ObjectiveFunction f = getObjectiveFunction();
    std::default_random_engine generator(seed);
    std::vector<std::vector<double>> particle(num_particles, std::vector<double>(num_dimensions));
    std::vector<std::vector<double>> particle_velocity(num_particles, std::vector<double>(num_dimensions));
    for (size_t j=0; j<num_dimensions; j++) {
        std::uniform_real_distribution<double> dist_ul1(lower[j], upper[j]);
        std::uniform_real_distribution<double> dist_ul2(-std::fabs(upper[j]-lower[j]), std::fabs(upper[j]-lower[j]));
        for (int i=0; i<num_particles; i++) {
            particle[i][j] = dist_ul1(generator);
            particle_velocity[i][j] = dist_ul2(generator);
        }
    }
    std::vector<std::vector<double>> best_particle_position(num_particles);
    std::vector<double> best_position(num_dimensions);
    for (int i=0; i<num_particles; i++) {
        best_particle_position[i] = particle[i];
        if (f(best_particle_position[i]) < f(best_position)) {
            best_position = best_particle_position[i];
        }
    }
    setX(best_position);

    // Main optimization loop.
    std::uniform_real_distribution<double> dist_01(0, 1);
    while(getStatus() != iteration_limit && getStatus() != time_limit) {
        for (int i=0; i<num_particles; i++) {
            for (size_t j=0; j<num_dimensions; j++) {
                double rp = dist_01(generator);
                double rg = dist_01(generator);
                particle_velocity[i][j] = inertia_weight * particle_velocity[i][j] +
                                          cognitive_coeff * rp * (best_particle_position[i][j] - particle[i][j]) +
                                          social_coeff * rg * (best_position[j] - particle[i][j]);
            }
            for (size_t j=0; j<num_dimensions; j++) {
                particle[i][j] = particle[i][j] + particle_velocity[i][j];
            }
            if (f(particle[i]) < f(best_particle_position[i])) {
                best_particle_position[i] = particle[i];
                if (f(particle[i]) < f(best_position)) {
                    best_position = particle[i];
                    setX(best_position);
                }
            }
        }
        addIterations(1);
    } // End while
}

}

#endif
