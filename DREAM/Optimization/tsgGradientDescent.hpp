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
 * An object of this class represents the components of a gradient descent run and all key parts of the swarm can be read through
 * member methods. Methods are also provided to modify a subset of its parts.
 *
 * \par Constructors and Copy/Move assignment
 * All constructors require information about an initial candidate point and initial stepsize. A method is available to retrieve
 * the number of dimensions. No default constructor is available for this class, but copy/move constructors and assignment operator=
 * overloads are available. The number of dimensions \b cannot be modified.
 * - GradientDescentState()
 * - getNumDimensions()
 *
 * \par Set and Get Candidate and Stepsize
 * Throughout the run of gradient descent, the candidate point and stepsize may change. Methods are available to retrieve these
 * data and set them manually.
 * - getCandidate(), setCandidate()
 * - getStepsize(), setStepsize()
 *
 * \par Gradient Descent Algorithm
 * The gradient descent optimization algorithm can be used on a gradient descent algorithm with the goal of minimizing a particular
 * objective function. The current candidate point, with respect to the given objective function, can be obtained from the state
 * using the getCandidate() method.
 * - GradientDescent()
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
    //! \brief The default constructor is allowed.
    GradientDescentState() = delete;
    //! \brief Constructor for a gradient descent state with the initial candidate and stepsize.
    GradientDescentState(const std::vector<double> &x, const double stepsize);
    //! \brief Copy constructor.
    GradientDescentState(const GradientDescentState &source) = default;
    //! \brief Move constructor.
    GradientDescentState(GradientDescentState &&source) = default;

    //! \brief Move assignment.
    GradientDescentState& operator=(GradientDescentState &&source) = default;

    //! \brief Return the number of dimensions.
    inline int getNumDimensions() const {return num_dimensions;}
    //! \brief Return the stepsize.
    inline double getStepsize() const {return stepsize;}
    //! \brief Return the current candidate point.
    inline void getCandidate(double x[]) const {std::copy_n(candidate.begin(), num_dimensions, x);}
    inline std::vector<double> getCandidate() const {return candidate;}

    //! \brief Set the stepsize.
    inline void setStepsize(const double ss) {stepsize = ss;}
    //! \brief Set the current candidate point.
    inline void setCandidate(const double x[]) {std::copy_n(x, num_dimensions, candidate.begin());}
    inline void setCandidate(const std::vector<double> &x) {
        checkVarSize("GradientDescentState::setCandidate", "candidate point", x.size(), num_dimensions);
        candidate = x;
    }

    /*! \brief Applies the classic (proximal) gradient descent algorithm to a gradient descent state.
     * \ingroup OptimizationAlgorithm Gradient Descent Algorithm
     *
     * Runs \b num_iterations of the gradient algorithm to a gradient descent \b state to minimize the function \b f over the
     * a domain implied by a given projection function. The parameters of the algorithm are (optional) line search coefficients.
     *
     * \param f Objective function to be minimized
     * \param g Gradient of the objective function
     * \param proj Projection function that orthogonally projects points into the domain of the objective function
     * \param num_iterations number of iterations to perform
     * \param state holds the state of the gradient descent algorithm, see TasOptimization::GradientDescentState
     * \param line_search_coeffs vector of floats which may be empty or of size=2; if it is empty, then a constant stepsize scheme
     *        is used; if it is non-empty, then an adaptive stepsize scheme is used that adapts to the local curvature of the
     *        function; in this case, the first (resp. second) entry controls quickly the stepsize is decreased (resp. increased)
     *        and both entries must be greater than 1
     */
    friend void GradientDescent(const ObjectiveFunction &f, const GradientFunction &g, const ProjectionFunction &proj,
                                const int num_iterations, GradientDescentState &state, const std::vector<double> &line_search_coeffs);

  protected:
    #ifndef __TASMANIAN_DOXYGEN_SKIP_INTERNAL
    //! \brief Returns a reference to the current candidate point.
    inline std::vector<double> &getCandidateRef() {return candidate;}
    #endif

  private:
    int num_dimensions;
    double stepsize;
    std::vector<double> candidate;
};

// Forward declarations.
void GradientDescent(const ObjectiveFunction &f, const GradientFunction &g, const ProjectionFunction &proj,
                     const int num_iterations, GradientDescentState &state,  const std::vector<double> &line_search_coeffs = {});

}

#endif
