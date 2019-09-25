/*
 * Copyright (c) 2017, Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
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
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
 * THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
 * THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
 * IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TASMANIAN_DREAM_CORE_PDF_HPP
#define __TASMANIAN_DREAM_CORE_PDF_HPP

/*!
 * \internal
 * \file tsgDreamCorePDF.hpp
 * \brief Gives the unscaled formulas for several probability distributions.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAM
 * \endinternal
 */

#include "tsgDreamEnumerates.hpp"

namespace TasDREAM{

//! \brief Returns the unscaled probability density of \b distribution (defined by \b params) at the point \b x.
//! \ingroup DREAMPDF

//! Variadric template that defines different one dimensional probability density function.
//! Each function is defined by the \b distribution type and a set of parameters \b params,
//! but the densities are not normalized (this is not needed by the DREAM algorithm).
//! \b Returns the value of the density function at \b x, but \b assumes \b x \b is \b in \b the \b domain,
//! i.e., if the probability distribution is restricted to an interval no check is performed.
//! For example, for uniform distribution over (a, b), the function will always return 1.
//!
//! Both the regular and log-form of the density can be computed based on \b form.
//!
//! \par Uniform distribution, dist_uniform,
//! No additional parameters are needed, always returns 1, e.g.,
//! \code double foo = getDensity<dist_uniform>(0.2); // foo will be 1.0 \endcode
//!
//! \par Gaussian distribution, dist_gaussian
//! Uses two parameters, mean and variance, defined in that order, e.g.,
//! \code double y = getDensity<dist_gaussian>(x, M, V); \endcode
//! then \f$ y = \exp(- 0.5 (x - M)^2 / V) \f$.
//! Note that this also works for truncated Gaussian distribution, since the range and scale are not considered.
//!
//! \par Exponential distribution, dist_exponential
//! Uses two parameters, lower bound and rate (or inverse scale), e.g.,
//! \code double y = getDensity<dist_exponential>(x, x0, R); \endcode
//! then \f$ y = \exp(- R * (x - x_0)) \f$.
//!
//! \par Beta distribution, dist-beta
//! Uses four parameters, lower and upper bounds and two shapes, e.g.,
//! \code double y = getDensity<dist_beta>(x, x0, x1, alpha, beta); \endcode
//! then \f$ y = (x - x_0)^{\alpha - 1} (x_1 - x_0)^{\beta - 1} \f$.
//!
//! \par Gamma distribution, dist-beta
//! Uses three parameters, lower bounds, shape and rate, e.g.,
//! \code double y = getDensity<dist_gamma>(x, x0, alpha, beta); \endcode
//! then \f$ y = (x - x_0)^{\alpha - 1} \exp(-\beta (x - x_0)) \f$.
//!
template<TypeDistribution distribution, TypeSamplingForm form = regform, typename... Params>
double getDensity(double x, Params... params){
    std::vector<typename std::tuple_element<0, std::tuple<Params...>>::type> ParameterArray = {params...};
    if (form == regform){
        if (distribution == dist_gaussian){
            return std::exp(-0.5 * (x - ParameterArray[0]) * (x - ParameterArray[0]) / ParameterArray[1]);
        }else if (distribution == dist_exponential){
            return std::exp(-ParameterArray[1] * (x - ParameterArray[0]));
        }else if (distribution == dist_beta){
            return std::pow(x - ParameterArray[0], ParameterArray[2] - 1.0) * std::pow(ParameterArray[1] - x, ParameterArray[3] - 1.0);
        }else if (distribution == dist_gamma){
            return std::pow(x - ParameterArray[0], ParameterArray[1] - 1.0) * std::exp(- ParameterArray[2] * (x - ParameterArray[0]));
        }else{ // uniform
            return 1.0;
        }
    }else{
        if (distribution == dist_gaussian){
            return -0.5 * (x - ParameterArray[0]) * (x - ParameterArray[0]) / ParameterArray[1];
        }else if (distribution == dist_exponential){
            return -ParameterArray[1] * (x - ParameterArray[0]);
        }else if (distribution == dist_beta){
            return std::log(x - ParameterArray[0]) * (ParameterArray[2] - 1.0) + std::log(ParameterArray[1] - x) * (ParameterArray[3] - 1.0);
        }else if (distribution == dist_gamma){
            return std::log(x - ParameterArray[0]) * (ParameterArray[1] - 1.0) - ParameterArray[2] * (x - ParameterArray[0]);
        }else{ // uniform
            return 0.0;
        }
    }
}

}


#endif
