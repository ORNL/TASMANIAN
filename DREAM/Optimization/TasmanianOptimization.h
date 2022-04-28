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

#ifndef __TASMANIAN_OPTIMIZATION_H
#define __TASMANIAN_OPTIMIZATION_H

using tsg_dream_random = double (*)();
using tsg_dream_domain = int (*)(int, const double[]); // first argument is num_dimensions
using tsg_optim_obj_fn = void (*)(int, int, const double[], double[]); // first argument is num_dimensions, second argument is num_points

// ------------------------------------------ C Interface for TasmanianOptimization ------------------------------------------ //
// NOTE: you need to explicitly call the constructor and destructor and, in all cases, void *state is a pointer to a C++ class

// See tsgParticleSwarm.hpp for the corresponding C++ interface.
void* tsgParticleSwarmState_Construct(int num_dimensions, int num_particles);
void tsgParticleSwarmState_Destruct(void* state);
int tsgParticleSwarmState_GetNumDimensions(void* state);
int tsgParticleSwarmState_GetNumParticles(void* state);
double* tsgParticleSwarmState_GetParticlePositions(void* state);
double* tsgParticleSwarmState_GetParticleVelocities(void* state);
double* tsgParticleSwarmState_GetBestParticlePositions(void* state);
double* tsgParticleSwarmState_GetBestPosition(void* state);
void tsgParticleSwarmState_SetNumDimensions(void* state, const int nd);
void tsgParticleSwarmState_SetNumParticles(void* state, const int np);
void tsgParticleSwarmState_SetParticlePositions(void* state, const double* pp);
void tsgParticleSwarmState_SetParticleVelocities(void* state, const double* pv);
void tsgParticleSwarmState_SetBestParticlePositions(void* state, const double* bpp);
void tsgParticleSwarmState_ClearBestParticles(void* state);
void tsgParticleSwarmState_ClearCache(void* state);
void tsgParticleSwarmState_InitializeParticlesInsideBox(void* state, const double* box_lower, const double* box_upper, const tsg_dream_random get_random01);
void tsgParticleSwarm(const tsg_optim_obj_fn f, const int num_iterations, const tsg_dream_domain inside, void *state, const double inertia_weight,
                      const double cognitive_coeff, const double social_coeff, const tsg_dream_random get_random01);

#endif
