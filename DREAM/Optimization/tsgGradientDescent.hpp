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

namespace TasOptimization {

class GradientDescentState {
  public:
    GradientDescentState() = delete;
    GradientDescentState(const std::vector<double> x, const double ss, std::vector<double> lsc = {});
    GradientDescentState(const GradientDescentState &source) = default;
    GradientDescentState(GradientDescentState &&source) = default;
    GradientDescentState& operator=(GradientDescentState &&source) = default;

    inline int getNumDimensions() const {return num_dimensions;}
    inline double getStepsize() const {return stepsize;}

    inline void getCandidate(double x[]) const {std::copy_n(candidate.begin(), num_dimensions, x);}
    inline std::vector<double> getCandidate() const {return candidate;}

    inline void getLineSearchCoeffs(double c[]) const {std::copy_n(line_search_coeffs.begin(), 2, c);}
    inline std::vector<double> getLineSearchCoeffs() const {return line_search_coeffs;}

    inline void setCandidate(const double x[]) {std::copy_n(x, num_dimensions, candidate.begin());}
    inline void setCandidate(const std::vector<double> &x) {
        checkVarSize("GradientDescentState::setCandidate", "candidate point", x.size(), num_dimensions);
        candidate = x;
    }

    inline void setLineSearchCoeffs(const double c[]) {std::copy_n(c, 2, line_search_coeffs.begin());}
    inline void setLineSearchCoeffs(const std::vector<double> &c) {
        checkVarSize("GradientDescentState::setLineSearchCoeffs", "line search coefficients", c.size(), 2);
        line_search_coeffs = c;
    }

    friend void GradientDescent(const ObjectiveFunction f, const GradientFunction grad, const ProjectionFunction proj,
                                const int num_iterations, GradientDescentState &state);

  protected:
    inline std::vector<double> &getCandidateRef() {return candidate;}

  private:
    int num_dimensions;
    double stepsize;
    std::vector<double> candidate, line_search_coeffs;
};

// Forward declarations.
void GradientDescent(const ObjectiveFunction f, const GradientFunction grad, const ProjectionFunction proj,
                     const int num_iterations, GradientDescentState &state);

}

#endif
