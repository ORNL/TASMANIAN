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

#ifndef __TASMANIAN_PARTICLE_SWARM_HPP
#define __TASMANIAN_PARTICLE_SWARM_HPP

#include "tsgOptimizationUtils.hpp"

/*!
 * \internal
 * \file tsgParticleSwarm.hpp
 * \brief Particle swarm state and algorithm.
 * \author Weiwei Kong & Miroslav Stoyanov
 * \ingroup TasmanianOptimization
 *
 * Definition of the particle swarm state class and the particle swarm algorithm.
 * \endinternal
 */

namespace TasOptimization {

/*!
 * \brief Keeps the history of a particle swarm run.
 * \ingroup OptimizationState
 *
 * \par Particle Swarm State
 * The particle swarm state consists of \b num_particles vectors of size \b num_dimensions, where both of these scalars are either
 * specified or inferred from input data. Each vector is associated with a particular particle's position \b particle_position,
 * velocity \b particle_velocity, and best known position \b best_particle_position from a previous invocation of the particle swarm
 * algorithm ParticleSwarm(). The last \b num_dimensions entries of \b best_particle position contain the position of the
 * overall best position in the swarm.
 *
 * \par Particle Swarm Cache
 * When the particle swarm state is modified by invoking the particle swarm algorithm ParticleSwarm() with a particular objective
 * function and domain, the results of invocation are stored in several cache variables prefixed by \b cache_. This is done
 * to speed up subsequent calls of ParticleSwarm() on this modified swarm state, under the assumption that the objective function
 * and domain from the previous call remains unchanged.
 *
 * \par
 * To apply the ParticleSwarm() algorithm on a swarm state with an objective function or domain that is different from its
 * last invocation of ParticleSwarm(), the cache must be cleared using the class method clearCache(). The behavior of the
 * particle swarm algorithm without doing this step is otherwise \a undefined.
 */
class ParticleSwarmState {
public:
    //! \brief The particle swarm state must be initialized with data that defines number of particles and dimensions.
    ParticleSwarmState() = delete;
    //! \brief Constructor for a particle swarm state with the number of particles and dimensions.
    ParticleSwarmState(const int num_dimensions, const int num_particles);
    //! \brief Constructor for a particle swarm state with the number of dimensions and the number of particles inferred from
    //! a set of input particle positions \b pp and a set of input particle velocities \b pv.
    ParticleSwarmState(const int num_dimensions, std::vector<double> &&pp, std::vector<double> &&pv);

    //! \brief Return the number of dimensions.
    inline int getNumDimensions() const {return num_dimensions;}
    //! \brief Return the number of particles.
    inline int getNumParticles() const {return num_particles;}
    //! \brief Return the particle positions.
    inline void getParticlePositions(double pp[]) const {std::copy_n(particle_positions.begin(), num_particles * num_dimensions, pp);}
    inline std::vector<double> getParticlePositions() const {return particle_positions;}
    //! \brief Return the particle velocities.
    inline void getParticleVelocities(double pv[]) const {std::copy_n(particle_velocities.begin(), num_particles * num_dimensions, pv);}
    inline std::vector<double> getParticleVelocities() const {return particle_velocities;}
    //! \brief Return the previously best known particle positions.
    inline void getBestParticlePositions(double bpp[]) const {std::copy_n(best_particle_positions.begin(), (num_particles + 1) * num_dimensions, bpp);}
    inline std::vector<double> getBestParticlePositions() const {return best_particle_positions;}
    //! \brief Return the previously best known position in the swarm.
    inline void getBestPosition(double bp[]) const {std::copy_n(best_particle_positions.begin() + num_particles * num_dimensions, num_dimensions, bp);}
    inline std::vector<double> getBestPosition() const {
        std::vector<double> best_position(num_dimensions);
        getBestPosition(best_position.data());
        return best_position;
    }
    //! \brief Return the state vector.
    inline void getStateVector(bool sv[]) const {
        sv[0] = positions_initialized;
        sv[1] = velocities_initialized;
        sv[2] = best_positions_initialized;
        sv[3] = cache_initialized;
    }
    inline std::vector<bool> getStateVector() const {
        return {positions_initialized, velocities_initialized, best_positions_initialized, cache_initialized};
    }

    //! \brief Set the particle positions.
    void setParticlePositions(const double pp[]) {std::copy_n(pp, num_dimensions * num_particles, particle_positions.begin());}
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
    //! \brief Set the particle velocities.
    void setParticleVelocities(const double pv[]) {std::copy_n(pv, num_dimensions * num_particles, particle_velocities.begin());}
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
    //! \brief Set the previously best known particle velocities.
    void setBestParticlePositions(const double bpp[]) {std::copy_n(bpp, num_dimensions * (num_particles + 1), best_particle_positions.begin());}
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

    //! \brief Clear the previously best known particle velocities.
    void clearBestParticles() {
        best_positions_initialized = false;
        best_particle_positions = std::vector<double>((num_particles + 1) * num_dimensions);
    }
    //! \brief Clear the particle swarm cache.
    void clearCache() {
        cache_initialized = false;
        std::fill(cache_particle_fvals.begin(), cache_particle_fvals.end(), 0.0);
        std::fill(cache_particle_inside.begin(), cache_particle_inside.end(), false);
        std::fill(cache_best_particle_fvals.begin(), cache_best_particle_fvals.end(), 0.0);
        std::fill(cache_best_particle_inside.begin(), cache_best_particle_inside.end(), false);
    }

    /*! \brief Randomly initializes all of the particle positions and velocities.
     *
     * The i-th component of each particle's position is uniformly sampled from the interval [\b box_lower[i], \b box_upper[i]].
     * The i-th velocity of each particle's velocity is uniformly sampled from the interval [-R, R] where R =
     * abs(\b box_upper[i] - \b box_lower[i]). The uniform [0,1] random number generator used in the sampling is specified
     * by \b get_random01.
     */
    void initializeParticlesInsideBox(const double box_lower[], const double box_upper[],
                                      const std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);
    void initializeParticlesInsideBox(const std::vector<double> &box_lower, const std::vector<double> &box_upper,
                                      const std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);


    /*! \brief Applies the particle swarm algorithm to a particle swarm state.
     * \ingroup OptimizationAlgorithm Particle Swarm Algorithm
     *
     * Runs \b num_iterations of the particle swarm algorithm to a particle swarm \b state to minimize the function \b f over the
     * domain \b inside. The parameters of the algorithm are \b inertia_weight , \b cognitive_coeff , and \b social_coeff. The
     * uniform [0,1] random number generator used by the algorithm is \b get_random01.
     *
     * \param f Objective function to be minimized
     * \param num_iterations number of iterations to perform
     * \param inside indicates whether a given point is inside or outside of the domain of interest
     * \param state holds the state of the particles, e.g., positions and velocities, see TasOptimization::ParticleSwarmState
     * \param inertia_weight inertial weight for the particle swarm algorithm
     * \param cognitive_coeff cognitive coefficient for the particle swarm algorithm
     * \param social_coeff social coefficient for the particle swarm algorithm
     * \param get_random01 random number generator, defaults to rand()
     *
     * \throws std::runtime_error if either the positions or the velocities of the \b state have not been initialized
     */
    friend void ParticleSwarm(const ObjectiveFunction f, const int num_iterations, const TasDREAM::DreamDomain inside,
                              ParticleSwarmState &state, const double inertia_weight, const double cognitive_coeff,
                              const double social_coeff, const std::function<double(void)> get_random01);

protected:
    #ifndef __TASMANIAN_DOXYGEN_SKIP_INTERNAL
    //! \brief Returns a reference to the particle positions.
    inline std::vector<double> &getParticlePositionsRef() {return particle_positions;}
    //! \brief Returns a reference to the particle velocities.
    inline std::vector<double> &getParticleVelocitiesRef() {return particle_velocities;}
    //! \brief Returns a reference to the previously best known particle positions.
    inline std::vector<double> &getBestParticlePositionsRef() {return best_particle_positions;}
    //! \brief Returns a reference to the cached function values.
    inline std::vector<double> &getCacheParticleFValsRef() {return cache_particle_fvals;}
    //! \brief Returns a reference to the previously best known function values.
    inline std::vector<double> &getCacheBestParticleFValsRef() {return cache_best_particle_fvals;}
    //! \brief Returns a reference to a boolean hash map of the domain of the loaded particles.
    inline std::vector<bool> &getCacheParticleInsideRef() {return cache_particle_inside;}
    //! \brief Returns a reference to a boolean hash map of the domain of the previously best known loaded particles.
    inline std::vector<bool> &getCacheBestParticleInsideRef() {return cache_best_particle_inside;}
    #endif

private:
    bool positions_initialized, velocities_initialized, best_positions_initialized, cache_initialized;
    int num_dimensions, num_particles;
    std::vector<double> particle_positions, particle_velocities, best_particle_positions, cache_particle_fvals, cache_best_particle_fvals;
    std::vector<bool> cache_particle_inside, cache_best_particle_inside;
};

// Forward declarations.
void ParticleSwarm(const ObjectiveFunction f, const int num_iterations, const TasDREAM::DreamDomain inside, ParticleSwarmState &state,
                   const double inertia_weight, const double cognitive_coeff, const double social_coeff,
                   const std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);

}
#endif
