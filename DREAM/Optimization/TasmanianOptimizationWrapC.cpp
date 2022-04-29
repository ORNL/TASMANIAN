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

#ifndef __TASMANIAN_OPTIMIZATION_WRAPC_CPP
#define __TASMANIAN_OPTIMIZATION_WRAPC_CPP

#include "tsgParticleSwarm.hpp"

// --------------------------- C Interface for use with Python ctypes and potentially other C codes --------------------------- //

// C Function Pointer Aliases
using tsg_dream_random = double (*)();
using tsg_optim_dom_fn = bool   (*)(const int, const double[]);
using tsg_optim_obj_fn = void   (*)(const int, const int, const double[], double[]);

namespace TasOptimization{

extern "C" {

    // Particle Swarm State.
    void* tsgParticleSwarmState_Construct(int num_dimensions, int num_particles) {
        return (void*) new ParticleSwarmState(num_dimensions, num_particles);
    }
    void tsgParticleSwarmState_Destruct(void* state) {
        delete reinterpret_cast<ParticleSwarmState*>(state);
    }
    int tsgParticleSwarmState_GetNumDimensions(void* state) {
        return reinterpret_cast<ParticleSwarmState*>(state)->getNumDimensions();
    }
    int tsgParticleSwarmState_GetNumParticles(void* state) {
        return reinterpret_cast<ParticleSwarmState*>(state)->getNumParticles();
    }
    void tsgParticleSwarmState_GetParticlePositions(void* state, double pp[]) {
        reinterpret_cast<ParticleSwarmState*>(state)->getParticlePositions(pp);
    }
    void tsgParticleSwarmState_GetParticleVelocities(void* state, double pv[]) {
        reinterpret_cast<ParticleSwarmState*>(state)->getParticleVelocities(pv);
    }
    void tsgParticleSwarmState_GetBestParticlePositions(void* state, double bpp[]) {
        reinterpret_cast<ParticleSwarmState*>(state)->getBestParticlePositions(bpp);
    }
    void tsgParticleSwarmState_GetBestPosition(void* state, double bp[]) {
        reinterpret_cast<ParticleSwarmState*>(state)->getBestPosition(bp);
    }
    void tsgParticleSwarmState_SetParticlePositions(void* state, const double pp[]) {
        reinterpret_cast<ParticleSwarmState*>(state)->setParticlePositions(pp);
    }
    void tsgParticleSwarmState_SetParticleVelocities(void* state, const double pv[]) {
        reinterpret_cast<ParticleSwarmState*>(state)->setParticleVelocities(pv);
    }
    void tsgParticleSwarmState_SetBestParticlePositions(void* state, const double bpp[]) {
        reinterpret_cast<ParticleSwarmState*>(state)->setBestParticlePositions(bpp);
    }
    void tsgParticleSwarmState_ClearBestParticles(void* state) {
        reinterpret_cast<ParticleSwarmState*>(state)->clearBestParticles();
    }
    void tsgParticleSwarmState_ClearCache(void* state) {
        reinterpret_cast<ParticleSwarmState*>(state)->clearCache();
    }
    void tsgParticleSwarmState_InitializeParticlesInsideBox(void* state, const double box_lower[], const double box_upper[],
                                                            const tsg_dream_random get_random01) {
        reinterpret_cast<ParticleSwarmState*>(state)->initializeParticlesInsideBox(box_lower, box_upper, get_random01);
    }

    // Particle Swarm Algorithm.
    void tsgParticleSwarm(const tsg_optim_obj_fn f_ptr, const int num_iterations, const tsg_optim_dom_fn inside_ptr, void *state,
                          const double inertia_weight, const double cognitive_coeff, const double social_coeff,
                          const tsg_dream_random get_random01_ptr) {
        auto f_cpp = [&](const std::vector<double> &x_batch, std::vector<double> &fval_batch)->void {
            int num_batch = fval_batch.size();
            int num_dims = x_batch.size() / num_batch;
            (*f_ptr)(num_dims, num_batch, x_batch.data(), fval_batch.data());
        };
        auto inside_cpp = [&](const std::vector<double> &x)->bool {
            int num_dims = x.size();
            return (*inside_ptr)(num_dims, x.data());
        };
        ParticleSwarm(f_cpp, num_iterations, inside_cpp, *(reinterpret_cast<ParticleSwarmState*>(state)), inertia_weight,
                      cognitive_coeff, social_coeff, *get_random01_ptr);
    }

}

}

#endif
