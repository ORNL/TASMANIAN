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
 * \ingroup OptimizationState
 * \brief Stores the information about a particle swarm.
 *
 * \par class ParticleSwarmState
 * A state class associated with the Particle Swarm State algorithm.
 *
 * \par Constructors and Copy/Move assignment
 * All constructors requires the number of dimensions and particles,
 * those \b cannot be modified but can be retrieved by get methods.
 * The class is movable and copyable by constructor or operator=.
 * - ParticleSwarmState()
 * - getNumDimensions()
 * - getNumParticles()
 *
 * \par Data layout
 * The logical organization of the particle data follows conventions similar to the Sparse Grid module.
 * Namely, the data is organized into strips of size \b num_dimensions and there is one strip
 * for each particle. Each particle has several vectors: position, best-known position and velocity;
 * as well as cached values for the objective functional at the position and best-position.
 *
 * \par Set and Get Particle Positions and Velocities
 * Each particle of the swarm is associated with a (possibly empty) position and velocity vector. Methods are also available to
 * determine if the positions and velocities are empty (uninitialized) or non-empty (initialized).
 * - getParticlePositions(), setParticlePositions()
 * - getParticleVelocities(), setParticleVelocities()
 * - isPositionInitialized(), isVelocityInitialized()
 *
 * \par Set and Get Best Particle Positions
 * A \b best known position will be assigned to each particle during the execution of the algorithm.
 * The \b best know position for the entire swarm is also available and it is the current "solution",
 * i.e., the point with smallest observed functional.
 * - getBestPosition()
 * - getBestParticlePositions(), setBestParticlePositions()
 * - clearBestParticles()
 * - isBestPositionInitialized()
 *
 * \par Initializing a Random Swarm
 * Helper methods are available for initializing random particle positions and velocities.
 * - initializeParticlesInsideBox()
 * - isBestPositionInitialized()
 *
 * \par Particle Swarm Algorithm
 * If all of particles of the swarm are initialized with velocities and positions, then the particle swarm
 * optimization algorithm can be used on the swarm with the goal of minimizing a particular objective function. The best
 * point, with respect to the given objective function, can be obtained from the swarm using the getBestPosition() method.
 * - ParticleSwarm()
 *
 * \par
 * More information about the particle swarm algorithm can be found in the following paper:
 *
 * \par
 * > M. R. Bonyadi and Z. Michalewicz, "Particle Swarm Optimization for Single Objective Continuous Space Problems: A Review," in
 * > <em>Evolutionary Computation</em>, vol. 25, no. 1, pp. 1-54, March 2017, doi: 10.1162/EVCO_r_00180.
 *
 * \par Clearing the Cache
 * After an optimization method is applied to the swarm, certain data related to the objective function are stored in cache variables
 * inside this class. If an optimization method needs to be applied to the swarm with a different objective function than the one used
 * to generate the cache, the cache \b must be cleared.
 * - clearCache()
 * - isCacheInitialized()
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
    //! \brief Copy constructor.
    ParticleSwarmState(const ParticleSwarmState &source) = default;
    //! \brief Move constructor.
    ParticleSwarmState(ParticleSwarmState &&source) = default;

    //! \brief Move assignment.
    ParticleSwarmState& operator=(ParticleSwarmState &&source) = default;

    //! \brief Return the number of dimensions.
    inline int getNumDimensions() const {return num_dimensions;}
    //! \brief Return the number of particles.
    inline int getNumParticles() const {return num_particles;}
    //! \brief Return the particle positions.
    inline void getParticlePositions(double pp[]) const {std::copy_n(particle_positions.begin(), num_particles * num_dimensions, pp);}
    //! \brief Return the particle positions, vector overload.
    inline std::vector<double> getParticlePositions() const {return particle_positions;}
    //! \brief Return the particle velocities.
    inline void getParticleVelocities(double pv[]) const {std::copy_n(particle_velocities.begin(), num_particles * num_dimensions, pv);}
    //! \brief Return the particle velocities, vector overload.
    inline std::vector<double> getParticleVelocities() const {return particle_velocities;}
    /*!
     * \brief Return the best known particle positions.
     *
     * The method actually returns positions with one more than the total number of particles,
     * the last \b num_dimensions entries of this vector contain the best particle position of the entire swarm.
     */
    inline void getBestParticlePositions(double bpp[]) const {std::copy_n(best_particle_positions.begin(), (num_particles + 1) * num_dimensions, bpp);}
    //! \brief Return the best known particle positions, vector overload.
    inline std::vector<double> getBestParticlePositions() const {return best_particle_positions;}
    //! \brief Loads the best known position in the swarm.
    inline void getBestPosition(double bp[]) const {std::copy_n(best_particle_positions.begin() + num_particles * num_dimensions, num_dimensions, bp);}
    //! \brief Returns the best knows position in the swarm.
    inline std::vector<double> getBestPosition() const {
        return std::vector<double>(best_particle_positions.begin() + num_particles * num_dimensions,
                                   best_particle_positions.begin() + num_particles * num_dimensions + num_dimensions);
    }

    //! \brief Returns true if the particle positions have been initialized.
    inline bool isPositionInitialized() const {return positions_initialized;}
    //! \brief Returns true if the particle velocities have been initialized.
    inline bool isVelocityInitialized() const {return velocities_initialized;}
    //! \brief Returns true if the best particle positions have been initialized.
    inline bool isBestPositionInitialized() const {return best_positions_initialized;}
    //! \brief Returns true if the cache has been initialized.
    inline bool isCacheInitialized() const {return cache_initialized;}
    //! \brief Return the complete state vector.
    inline std::vector<bool> getStateVector() const {
        return {positions_initialized, velocities_initialized, best_positions_initialized, cache_initialized};
    }

    //! \brief Set the particle positions, raw-array variant.
    void setParticlePositions(const double pp[]) {
        std::copy_n(pp, num_dimensions * num_particles, particle_positions.begin());
        positions_initialized = true;
    }
    //! \brief Set the particle positions, vector variant.
    void setParticlePositions(const std::vector<double> &pp) {
        checkVarSize("ParticleSwarmState::setParticlePositions", "particle position", pp.size(), num_dimensions * num_particles);
        particle_positions = pp;
        positions_initialized = true;
    }
    //! \brief Set the particle positions, with a move.
    void setParticlePositions(std::vector<double> &&pp) {
        checkVarSize("ParticleSwarmState::setParticlePositions", "particle positions", pp.size(), num_dimensions * num_particles);
        particle_positions = std::move(pp);
        positions_initialized = true;
    }
    //! \brief Set the particle velocities.
    void setParticleVelocities(const double pv[]) {
        std::copy_n(pv, num_dimensions * num_particles, particle_velocities.begin());
        velocities_initialized = true;
    }
    //! \brief Sets the best position.
    void setParticleVelocities(const std::vector<double> &pv) {
        checkVarSize("ParticleSwarmState::setParticleVelocities", "particle velocities", pv.size(), num_dimensions * num_particles);
        particle_velocities = pv;
        velocities_initialized = true;
    }
    //! \brief Sets the best position, with a move.
    void setParticleVelocities(std::vector<double> &&pv) {
        checkVarSize("ParticleSwarmState::setParticleVelocities", "particle velocities", pv.size(), num_dimensions * num_particles);
        particle_velocities = std::move(pv);
        velocities_initialized = true;
    }
    //! \brief Set the previously best known particle velocities.
    void setBestParticlePositions(const double bpp[]) {
        std::copy_n(bpp, num_dimensions * (num_particles + 1), best_particle_positions.begin());
        best_positions_initialized = true;
    }
    //! \brief Sets the best position per particle.
    void setBestParticlePositions(const std::vector<double> &bpp) {
        checkVarSize("ParticleSwarmState::setBestParticlePositions", "best particle positions", bpp.size(), num_dimensions * (num_particles + 1));
        best_particle_positions = bpp;
        best_positions_initialized = true;
    }
    //! \brief Sets the best position per particle, allows for a move.
    void setBestParticlePositions(std::vector<double> &&bpp) {
        checkVarSize("ParticleSwarmState::setBestParticlePositions", "best particle positions", bpp.size(), num_dimensions * (num_particles + 1));
        best_particle_positions = std::move(bpp);
        best_positions_initialized = true;
    }

    //! \brief Clear the previously best known particle velocities.
    void clearBestParticles() {
        best_positions_initialized = false;
        std::fill(best_particle_positions.begin(), best_particle_positions.end(), 0.0);
    }
    //! \brief Clear the particle swarm cache.
    void clearCache() {
        cache_initialized = false;
        std::fill(cache_particle_fvals.begin(), cache_particle_fvals.end(), 0.0);
        std::fill(cache_particle_inside.begin(), cache_particle_inside.end(), false);
        std::fill(cache_best_particle_fvals.begin(), cache_best_particle_fvals.end(), 0.0);
        std::fill(cache_best_particle_inside.begin(), cache_best_particle_inside.end(), false);
    }

    /*! \brief Randomly initializes all of the particle positions and velocities inside of a box.
     *
     * The i-th component of each particle's position is uniformly sampled from the interval [\b box_lower[i], \b box_upper[i]].
     * The i-th velocity of each particle's velocity is uniformly sampled from the interval [-R, R] where R =
     * abs(\b box_upper[i] - \b box_lower[i]). The uniform [0,1] random number generator used in the sampling is specified
     * by \b get_random01.
     */
    void initializeParticlesInsideBox(const double box_lower[], const double box_upper[],
                                      const std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);
    //! \brief Randomly initializes all of the particles, vector API overload.
    void initializeParticlesInsideBox(const std::vector<double> &box_lower, const std::vector<double> &box_upper,
                                      const std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);

    friend void ParticleSwarm(const ObjectiveFunction f, const TasDREAM::DreamDomain inside, const double inertia_weight,
                              const double cognitive_coeff, const double social_coeff, const int num_iterations,
                              ParticleSwarmState &state, const std::function<double(void)> get_random01);

private:
    bool positions_initialized, velocities_initialized, best_positions_initialized, cache_initialized;
    int num_dimensions, num_particles;
    std::vector<double> particle_positions, particle_velocities, best_particle_positions, cache_particle_fvals, cache_best_particle_fvals;
    std::vector<bool> cache_particle_inside, cache_best_particle_inside;
};

// Forward declarations.

/*!
 * \brief Applies the classic particle swarm algorithm to a particle swarm state.
 * \ingroup OptimizationAlgorithm
 *
 * Runs \b num_iterations of the particle swarm algorithm to a particle swarm \b state to minimize the function \b f over the
 * domain \b inside. The parameters of the algorithm are \b inertia_weight , \b cognitive_coeff , and \b social_coeff. The
 * uniform [0,1] random number generator used by the algorithm is \b get_random01.
 *
 * \param f Objective function to be minimized
 * \param inside indicates whether a given point is inside or outside of the domain of interest
 * \param inertia_weight inertial weight for the particle swarm algorithm
 * \param cognitive_coeff cognitive coefficient for the particle swarm algorithm
 * \param social_coeff social coefficient for the particle swarm algorithm
 * \param num_iterations number of iterations to perform
 * \param state holds the state of the particles, e.g., positions and velocities, see TasOptimization::ParticleSwarmState
 * \param get_random01 random number generator, defaults to rand()
 *
 * \throws std::runtime_error if either the positions or the velocities of the \b state have not been initialized
 */
void ParticleSwarm(const ObjectiveFunction f, const TasDREAM::DreamDomain inside, const double inertia_weight,
                   const double cognitive_coeff, const double social_coeff, const int num_iterations,
                   ParticleSwarmState &state, const std::function<double(void)> get_random01 = TasDREAM::tsgCoreUniform01);

}
#endif
