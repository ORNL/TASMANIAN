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

#ifndef __TASMANIAN_PARTICLE_SWARM_SOLVER_HPP
#define __TASMANIAN_PARTICLE_SWARM_SOLVER_HPP

#include "tsgOptimizationUtils.hpp"

namespace TasOptimization {

class ParticleSwarmState {
  public:
    ParticleSwarmState() = delete;
    ParticleSwarmState(int num_dimensions, int num_particles);
    // pp = particle positions, pv = particle velocities
    ParticleSwarmState(int num_dimensions, std::vector<double> &&pp, std::vector<double> &&pv);

    inline int getNumDimensions() const {return num_dimensions;}
    inline int getNumParticles() const {return num_particles;}
    inline std::vector<double> getParticlePositions() const {return particle_positions;}
    inline std::vector<double> getParticleVelocities() const {return particle_velocities;}
    inline std::vector<double> getBestParticlePositions() const {return best_particle_positions;}
    inline std::vector<double> getBestPosition() const {
        std::vector<double> best_position(num_dimensions);
        std::copy_n(best_particle_positions.begin() + num_particles * num_dimensions, num_dimensions, best_position.begin());
        return best_position;
    }
    inline std::vector<bool> getStateVector() const {
        return {positions_initialized, velocities_initialized, best_positions_initialized, cache_initialized};
    }

    void setNumDimensions(const int nd) {num_dimensions = nd;};
    void setNumParticles(const int np) {num_particles = np;};
    void setParticlePositions(const std::vector<double> &pp) {
        checkVarSize("ParticleSwarmState::setParticlePositions", "particle position", pp.size(), num_dimensions * num_particles);
        particle_positions = pp;
        positions_initialized = true;
    }
    void setParticlePositions(std::vector<double> &&pp) {
        checkVarSize("ParticleSwarmState::setParticlePositions", "particle positions", pp.size(), num_dimensions * num_particles);
        particle_positions = std::move(pp);
        positions_initialized = true;
    }
    void setParticleVelocities(const std::vector<double> &pv) {
        checkVarSize("ParticleSwarmState::setParticleVelocities", "particle velocities", pv.size(), num_dimensions * num_particles);
        particle_velocities = pv;
        velocities_initialized = true;
    }
    void setParticleVelocities(std::vector<double> &&pv) {
        checkVarSize("ParticleSwarmState::setParticleVelocities", "particle velocities", pv.size(), num_dimensions * num_particles);
        particle_velocities = std::move(pv);
        velocities_initialized = true;
    }
    void setBestParticlePositions(const std::vector<double> &bpp) {
        checkVarSize("ParticleSwarmState::setBestParticlePositions", "best particle positions", bpp.size(), num_dimensions * (num_particles + 1));
        best_particle_positions = bpp;
        best_positions_initialized = true;
    }
    void setBestParticlePositions(std::vector<double> &&bpp) {
        checkVarSize("ParticleSwarmState::setBestParticlePositions", "best particle positions", bpp.size(), num_dimensions * (num_particles + 1));
        best_particle_positions = std::move(bpp);
        best_positions_initialized = true;
    }

    void clearBestParticles() {
        best_positions_initialized = false;
        best_particle_positions = std::vector<double>((num_particles + 1) * num_dimensions);
    }
    void clearCache() {
        cache_initialized = false;
        cache_particle_fvals = std::vector<double>(num_particles);
        cache_particle_inside = std::vector<bool>(num_particles);
        cache_best_particle_fvals = std::vector<double>(num_particles + 1);
        cache_best_particle_inside = std::vector<bool>(num_particles + 1);
    }

    void initializeParticlesInsideBox(const std::vector<double> &box_lower, const std::vector<double> &box_upper,
                                      const std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);

    friend void ParticleSwarm(ObjectiveFunction f, int max_iterations, TasDREAM::DreamDomain inside, ParticleSwarmState &state,
                              double inertia_weight, double cognitive_coeff, double social_coeff, std::function<double(void)> get_random01);

  protected:
    inline std::vector<double> &getParticlePositionsRef() {return particle_positions;}
    inline std::vector<double> &getParticleVelocitiesRef() {return particle_velocities;}
    inline std::vector<double> &getBestParticlePositionsRef() {return best_particle_positions;}
    inline std::vector<double> &getCacheParticleFValsRef() {return cache_particle_fvals;}
    inline std::vector<double> &getCacheBestParticleFValsRef() {return cache_best_particle_fvals;}
    inline std::vector<bool> &getCacheParticleInsideRef() {return cache_particle_inside;}
    inline std::vector<bool> &getCacheBestParticleInsideRef() {return cache_best_particle_inside;}

    bool positions_initialized, velocities_initialized, best_positions_initialized, cache_initialized;
    int num_dimensions, num_particles;
    std::vector<double> particle_positions, particle_velocities, best_particle_positions, cache_particle_fvals, cache_best_particle_fvals;
    std::vector<bool> cache_particle_inside, cache_best_particle_inside;
};

// Forward declarations.
void ParticleSwarm(ObjectiveFunction f, int max_iterations, TasDREAM::DreamDomain inside, ParticleSwarmState &state,
                   double inertia_weight, double cognitive_coeff, double social_coeff,
                   std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);

}
#endif
