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
    ParticleSwarmState();
    ~ParticleSwarmState() = default;

    inline std::vector<double> getParticlePositions() const {return particle_positions;}
    inline std::vector<double> getParticleVelocities() const {return particle_velocities;}
    inline std::vector<double> getBestParticlePositions() const {return best_particle_positions;}
    inline std::vector<double> getBestPosition() const {return best_position;}
    inline int getNumIterations() {return num_iterations;}
    inline int getNumDimensions() {return num_dimensions;}
    inline int getNumParticles() {return num_particles;}

    void setParticlePositions(std::vector<double> &pp);
    void setParticleVelocities(std::vector<double> &pv);
    void setBestParticlePositions(std::vector<double> &bpp);
    void setBestPosition(std::vector<double> &bp);
    void setNumIterations(int it) {num_iterations = it;}
    void setNumDimensions(int nd) {num_dimensions = nd;};
    void setNumParticles(int np) {num_particles = np;};

    inline void addIterations(int k) {num_iterations += k;}
    void addParticlesInsideBox(int num_particles, std::vector<double> &lower, std::vector<double> &upper,
                               std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);

  private:
    std::vector<double> particle_positions, particle_velocities, best_particle_positions, best_position;
    int num_iterations, num_dimensions, num_particles;
};

void ParticleSwarm(BatchedObjectiveFunction f, int max_iterations, TasDREAM::DreamDomain inside, ParticleSwarmState &state,
                   double inertia_weight = 0.5, double cognitive_coeff = 2, double social_coeff = 0.5,
                   std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);

}
#endif
