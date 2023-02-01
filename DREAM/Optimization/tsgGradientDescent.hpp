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

#ifndef __TASMANIAN_GRADIENT_DESCENT_HPP
#define __TASMANIAN_GRADIENT_DESCENT_HPP

#include "tsgOptimizationUtils.hpp"

/*!
 * \internal
 * \file tsgGradientDescent.hpp
 * \brief Gradient descent state and algorithm.
 * \author Weiwei Kong & Miroslav Stoyanov
 * \ingroup TasmanianOptimization
 *
 * Definition of the gradient descent state class and the particle swarm algorithm.
 * \endinternal
 */

namespace TasOptimization {

/*!
 * \ingroup OptimizationState
 * \brief Stores the information about a gradient descent run.
 *
 * \par class GradientDescentState
 * A state class associated with the Gradient Descent algorithm.
 *
 * \par Constructors and Copy/Move assignment
 * Constructors require information about an initial candidate point and initial step-size.
 * The class is movable and copyable by constructor or operator=.
 * Once set, the number of dimensions \b cannot be modified.
 * - GradientDescentState()
 * - getNumDimensions()
 *
 * \par Set and Get Candidate and Stepsize
 * The Gradient Descent algorithm has one parameter, the step-size,
 * and the current (best) point.
 * The two have set/get methods with overloads for both std::vector and raw-arrays.
 * Using common math conventions, the current best point is designated by \b X.
 * - getX(), setX()
 * - getAdaptiveStepsize(), setAdaptiveStepsize()
 *
 * \par Gradient Descent Algorithm
 * See TasOptimization::GradientDescent()
 *
 * \par
 * More information about the gradient descent algorithm can be found in the following paper:
 *
 * \par
 * > Nesterov, Y. (2013). Gradient methods for minimizing composite functions. <em>Mathematical programming</em>, 140(1), 125-161.
 *
 */
class GradientDescentState {
  public:
    //! \brief The default constructor is NOT allowed.
    GradientDescentState() = delete;
    //! \brief Constructor for a gradient descent state with the initial candidate \b x and stepsize \b lambda0.
    GradientDescentState(const std::vector<double> &x0, const double initial_stepsize) :
        adaptive_stepsize(initial_stepsize), x(x0) {};

    //! \brief Copy constructor.
    GradientDescentState(const GradientDescentState &source) = default;
    //! \brief Move constructor.
    GradientDescentState(GradientDescentState &&source) = default;

    //! \brief Move assignment.
    GradientDescentState& operator=(GradientDescentState &&source) = default;
    //! \brief Copy assignment.
    GradientDescentState& operator=(GradientDescentState &source) = default;

    //! \brief Implicit conversion to the current candidate \b x by reference.
    inline operator std::vector<double>&() {return x;};

    //! \brief Return the number of dimensions.
    inline size_t getNumDimensions() const {return x.size();}
    //! \brief Return the stepsize.
    inline double getAdaptiveStepsize() const {return adaptive_stepsize;}
    //! \brief Return the current candidate point.
    inline void getX(double x_out[]) const {std::copy_n(x.begin(), x.size(), x_out);}
    //! \brief Overload for when the output is a vector.
    inline std::vector<double> getX() const {return x;}

    //! \brief Set the stepsize.
    inline void setAdaptiveStepsize(const double new_stepsize) {adaptive_stepsize = new_stepsize;}
    //! \brief Set the current candidate point.
    inline void setX(const double x_new[]) {std::copy_n(x_new, x.size(), x.begin());}
    //! \brief Overload for when the input is a vector.
    inline void setX(const std::vector<double> &x_new) {
        checkVarSize("GradientDescentState::setCandidate", "candidate point", x_new.size(), x.size());
        x = x_new;
    }

    friend OptimizationStatus GradientDescent(const ObjectiveFunctionSingle &func, const GradientFunctionSingle &grad,
                                              const ProjectionFunctionSingle &proj, const double increase_coeff,
                                              const double decrease_coeff, const int max_iterations, const double tolerance,
                                              GradientDescentState &state);
    friend OptimizationStatus GradientDescent(const ObjectiveFunctionSingle &func, const GradientFunctionSingle &grad,
                                              const double increase_coeff, const double decrease_coeff, const int max_iterations,
                                              const double tolerance, GradientDescentState &state);
    friend OptimizationStatus GradientDescent(const GradientFunctionSingle &grad, const double stepsize, const int max_iterations,
                                              const double tolerance, std::vector<double> &state);

  private:
    double adaptive_stepsize;
    std::vector<double> x;
};

/*!
 * \brief Applies the constant step-size gradient descent algorithm for functions with unbounded domains.
 * \ingroup OptimizationAlgorithm
 *
 * Minimize a function with gradient \b g over an unconstrained domain.
 * Perform work until reaching the desired tolerance (measured in the stationarity residual),
 * or until \b max_iterations is reached.
 * See also TasOptimization::computeStationarityResidual()
 *
 * \param grad Gradient of the objective functional
 * \param stepsize is the step-size of the algorithm
 * \param max_iterations is the maximum number of iterations to perform
 * \param tolerance Stationarity tolerance; the algorithm terminates when the stationarity residual computed by
 *        TasOptimization::computeStationarityResidual() is less than or equal to \b tolerance
 * \param state contains the current iterate and returns the best iterate.
 *        This algorithm does not use the adaptive step-size, so the state can be just a vector,
 *        but the signature accepts a GradientDescentState with an automatic conversion.
 *
 * \returns TasOptimization::OptimizationStatus struct that contains information about the last iterate.
 */
OptimizationStatus GradientDescent(const GradientFunctionSingle &grad, const double stepsize, const int max_iterations,
                                   const double tolerance, std::vector<double> &state);

/*!
 * \brief Applies the adaptive gradient descent algorithm on unrestricted domain.
 * \ingroup OptimizationAlgorithm Adaptive Non-Proximal Gradient Descent Algorithm
 *
 * Similar to the constant step-size algorithm GradientDescent() but applying an adaptive stepping.
 * This method is guaranteed to converge to a stationary point
 * if the gradient of \b f is Lipschitz continuous on its domain.
 * The algorithm is known as Non-Proximal, i.e., no restriction is applied to the domain
 * which implies either work on an unbounded domain or the starting point and the minimum
 * are sufficiently far from the boundary and the restriction is not needed.
 *
 * This variant requires the value of the functional that is to be minimized, in addition to the gradient.
 * There are two control parameters \b increase_coeff and \b decrease_coeff that guide the rate
 * at which the step-size is adjusted.
 * The parameters can affect the convergence rate, but not the final result.
 *
 * \param func is the objective function to be minimized
 * \param grad is the gradient of the objective function
 * \param increase_coeff Controls how quickly the step-size is increased; should be greater than 1
 * \param decrease_coeff Controls how quickly the step-size is decreased; should be greater than 1
 * \param max_iterations Maximum number of iterations to perform
 * \param tolerance same as in GradientDescent()
 * \param state Holds the state of the gradient descent algorithm, including the current iterate and the current adaptive step-size.
 *
 * \returns TasOptimization::OptimizationStatus struct that contains information about the last iterate.
 */
OptimizationStatus GradientDescent(const ObjectiveFunctionSingle &func, const GradientFunctionSingle &grad,
                                   const double increase_coeff,  const double decrease_coeff, const int max_iterations,
                                   const double tolerance, GradientDescentState &state);
/*!
 * \brief Applies the adaptive gradient descent algorithm on a restricted domain.
 * \ingroup OptimizationAlgorithm
 *
 * Similar to the adaptive step-size algorithm on the unrestricted domain,
 * but it uses a projection function to constrain each iterate to a user-defined domain.
 *
 * The \b proj function computes the orthogonal projection of a point inside the domain,
 * e.g., restricts the point to a hypercube.
 */
OptimizationStatus GradientDescent(const ObjectiveFunctionSingle &func, const GradientFunctionSingle &grad,
                                   const ProjectionFunctionSingle &proj, const double increase_coeff,
                                   const double decrease_coeff, const int max_iterations, const double tolerance,
                                   GradientDescentState &state);

}

#endif
