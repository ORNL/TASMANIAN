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
#include "tsgGradientDescent.hpp"

// --------------------------- C Interface for use with Python ctypes and potentially other C codes --------------------------- //

// C Function Pointer Aliases
using tsg_dream_random         = double (*)();
using tsg_optim_dom_fn         = int    (*)(const int, const double[], int[]);
using tsg_optim_obj_fn         = void   (*)(const int, const int, const double[], double[], int[]);
using tsg_optim_obj_fn_single  = double (*)(const int, const double[], int[]);
using tsg_optim_grad_fn_single = void   (*)(const int, const double[], double[], int[]);
using tsg_optim_proj_fn_single = void   (*)(const int, const double[], double[], int[]);

namespace TasOptimization{

// Helper methods to generate some C++ functions from C functions.
ObjectiveFunctionSingle convert_C_obj_fn_single(tsg_optim_obj_fn_single func_ptr, std::string err_msg) {
    return [=](const std::vector<double> &x_single)->double {
        int err_code = 0;
        int num_dims = x_single.size();
        double result = (*func_ptr)(num_dims, x_single.data(), &err_code);
        if (err_code != 0) throw std::runtime_error(err_msg);
        return result;
    };
}
GradientFunctionSingle convert_C_grad_fn_single(tsg_optim_grad_fn_single grad_ptr, std::string err_msg) {
    return [=](const std::vector<double> &x_single, std::vector<double> &grad)->void {
        int err_code = 0;
        int num_dims = x_single.size();
        (*grad_ptr)(num_dims, x_single.data(), grad.data(), &err_code);
        if (err_code != 0) throw std::runtime_error(err_msg);
    };
}
ProjectionFunctionSingle convert_C_proj_fn_single(tsg_optim_proj_fn_single proj_ptr, std::string err_msg) {
    return [=](const std::vector<double> &x_single, std::vector<double> &proj)->void {
        int err_code = 0;
        int num_dims = x_single.size();
        (*proj_ptr)(num_dims, x_single.data(), proj.data(), &err_code);
        if (err_code != 0) throw std::runtime_error(err_msg);
    };
}

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

    int tsgParticleSwarmState_IsPositionInitialized(void* state) {
        return reinterpret_cast<ParticleSwarmState*>(state)->isPositionInitialized();
    }
    int tsgParticleSwarmState_IsVelocityInitialized(void* state) {
        return reinterpret_cast<ParticleSwarmState*>(state)->isVelocityInitialized();
    }
    int tsgParticleSwarmState_IsBestPositionInitialized(void* state) {
        return reinterpret_cast<ParticleSwarmState*>(state)->isBestPositionInitialized();
    }
    int tsgParticleSwarmState_IsCacheInitialized(void* state) {
        return reinterpret_cast<ParticleSwarmState*>(state)->isCacheInitialized();
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
                                                            const char* random_type, const int random_seed, tsg_dream_random random_callback) {
        // Create the U[0,1] random number generator.
        std::minstd_rand park_miller((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed);
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        std::string rtype(random_type);
        auto randgen = [&]()->
                       std::function<double(void)>{
            if (rtype == "default") {
                srand((unsigned int) ((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed));
                return [&]()->double{ return TasDREAM::tsgCoreUniform01(); };
            } else if (rtype == "minstd_rand") {
                return [&]()->double{ return unif(park_miller); };
            } else {
                return [&]()->double{ return random_callback(); };
            }
        }();
        reinterpret_cast<ParticleSwarmState*>(state)->initializeParticlesInsideBox(box_lower, box_upper, randgen);
    }

    // Particle Swarm Algorithm.
    void tsgParticleSwarm(const tsg_optim_obj_fn f_ptr, const tsg_optim_dom_fn inside_ptr, const double inertia_weight,
                          const double cognitive_coeff, const double social_coeff, const int num_iterations,  void *state,
                          const char* random_type, const int random_seed, tsg_dream_random random_callback, int *err) {
        *err = 1;
        // Create the U[0,1] random number generator.
        std::minstd_rand park_miller((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed);
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        std::string rtype(random_type);
        auto randgen = [&]()->
                       std::function<double(void)>{
            if (rtype == "default") {
                srand((unsigned int) ((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed));
                return [&]()->double{ return TasDREAM::tsgCoreUniform01(); };
            } else if (rtype == "minstd_rand") {
                return [&]()->double{ return unif(park_miller); };
            } else {
                return [&]()->double{ return random_callback(); };
            }
        }();
        auto f_cpp = [=](const std::vector<double> &x_batch, std::vector<double> &fval_batch)->void {
            int err_code = 0;
            int num_batch = fval_batch.size();
            int num_dims = x_batch.size() / num_batch;
            (*f_ptr)(num_dims, num_batch, x_batch.data(), fval_batch.data(), &err_code);
            if (err_code != 0) throw std::runtime_error("The Python objective function callback returned an error in tsgParticleSwarm()");
        };
        auto inside_cpp = [=](const std::vector<double> &x)->bool {
            int err_code = 0;
            int num_dims = x.size();
            bool inside = (*inside_ptr)(num_dims, x.data(), &err_code);
            if (err_code != 0) throw std::runtime_error("The Python domain function callback returned an error in tsgParticleSwarm()");
            return inside;
        };
        try {
            ParticleSwarm(f_cpp, inside_cpp, inertia_weight, cognitive_coeff, social_coeff, num_iterations,
                          *(reinterpret_cast<ParticleSwarmState*>(state)), randgen);
            *err = 0; // Success
        } catch (std::runtime_error &) {}
    }

    // Gradient Descent State.
    void* tsgGradientDescentState_Construct(const int num_dimensions, const double x0[], const double initial_stepsize) {
        return (void*) new GradientDescentState(std::vector<double>(x0, x0 + num_dimensions), initial_stepsize);
    }
    void tsgGradientDescentState_Destruct(void* state) {
        delete reinterpret_cast<GradientDescentState*>(state);
    }
    int tsgGradientDescentState_GetNumDimensions(void* state) {
        return reinterpret_cast<GradientDescentState*>(state)->getNumDimensions();
    }
    double tsgGradientDescentState_GetAdaptiveStepsize(void* state) {
        return reinterpret_cast<GradientDescentState*>(state)->getAdaptiveStepsize();
    }
    void tsgGradientDescentState_GetX(void* state, double x_out[]) {
        reinterpret_cast<GradientDescentState*>(state)->getX(x_out);
    }
    void tsgGradientDescentState_SetAdaptiveStepsize(void* state, const double new_stepsize) {
        reinterpret_cast<GradientDescentState*>(state)->setAdaptiveStepsize(new_stepsize);
    }
    void tsgGradientDescentState_SetX(void* state, double x_new[]) {
        reinterpret_cast<GradientDescentState*>(state)->setX(x_new);
    }

    // Adaptive Stepsize Projected Gradient Descent Algorithm.
    OptimizationStatus tsgGradientDescent_AdaptProj(const tsg_optim_obj_fn_single func_ptr, const tsg_optim_grad_fn_single grad_ptr,
                                                    const tsg_optim_proj_fn_single proj_ptr, const double increase_coeff,
                                                    const double decrease_coeff, const int max_iterations, const double tolerance,
                                                    void* state, int* err) {

        *err = 1;
        // Convert C functions to safe C++ functions.
        ObjectiveFunctionSingle func_cpp = convert_C_obj_fn_single(func_ptr, "The Python objective function callback returned an error in tsgGradientDescent()");
        GradientFunctionSingle grad_cpp = convert_C_grad_fn_single(grad_ptr, "The Python gradient function callback returned an error in tsgGradientDescent()");
        ProjectionFunctionSingle proj_cpp = convert_C_proj_fn_single(proj_ptr, "The Python projection function callback returned an error in tsgGradientDescent()");

        // Main call and error handling.
        OptimizationStatus status;
        try {
            status = GradientDescent(func_cpp, grad_cpp, proj_cpp, increase_coeff, decrease_coeff, max_iterations, tolerance,
                                     *(reinterpret_cast<GradientDescentState*>(state)));
            *err = 0; // Success
        } catch (std::runtime_error &) {}
        return status;
    }

    // Adaptive Stepsize (Unconstrained) Gradient Descent Algorithm.
    OptimizationStatus tsgGradientDescent_Adapt(const tsg_optim_obj_fn_single func_ptr, const tsg_optim_grad_fn_single grad_ptr,
                                                const double increase_coeff, const double decrease_coeff, const int max_iterations,
                                                const double tolerance, void* state, int* err) {

        *err = 1;
        // Convert C functions to safe C++ functions.
        ObjectiveFunctionSingle func_cpp = convert_C_obj_fn_single(func_ptr, "The Python objective function callback returned an error in tsgGradientDescent()");
        GradientFunctionSingle grad_cpp = convert_C_grad_fn_single(grad_ptr, "The Python gradient function callback returned an error in tsgGradientDescent()");

        // Main call and error handling.
        OptimizationStatus status;
        try {
            status = GradientDescent(func_cpp, grad_cpp, increase_coeff, decrease_coeff, max_iterations, tolerance,
                                     *(reinterpret_cast<GradientDescentState*>(state)));
            *err = 0; // Success
        } catch (std::runtime_error &) {}
        return status;
    }

    // Constant Stepsize (Unconstrained) Gradient Descent Algorithm.
    OptimizationStatus tsgGradientDescent_Const(const tsg_optim_grad_fn_single grad_ptr, const double stepsize, const int max_iterations,
                                                const double tolerance, void* state, int* err) {

        *err = 1;
        // Convert C functions to safe C++ functions.
        GradientFunctionSingle grad_cpp = convert_C_grad_fn_single(grad_ptr, "The Python gradient function callback returned an error in tsgGradientDescent()");

        // Main call and error handling.
        OptimizationStatus status;
        try {
            status = GradientDescent(grad_cpp, stepsize, max_iterations, tolerance, *(reinterpret_cast<GradientDescentState*>(state)));
            *err = 0; // Success
        } catch (std::runtime_error &) {}
        return status;
    }

}

}

#endif
