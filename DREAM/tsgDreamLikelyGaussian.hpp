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

#ifndef __TASMANIAN_DREAM_LIKELY_GAUSS_HPP
#define __TASMANIAN_DREAM_LIKELY_GAUSS_HPP

#include <numeric>
#include <stdexcept>

#include "tsgDreamLikelihoodCore.hpp"

//! \file tsgDreamLikelyGaussian.hpp
//! \brief Several likelihood implementations based on Gaussian noise.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAM

namespace TasDREAM{

//! \brief Implements likelihood under the assumption of isotropic white noise.
//! \ingroup DREAMLikelihood

//! \par Gaussian Likelihood
//! The general formula for Gaussian likelihood is \f$ L(y | d_1 \cdots d_n) = \exp\left( 0.5 sum_{i=1}^n (y - d_i)^T \Sigma^{-1} (y - d_i) \right)\f$
//! where the \f$ d_i \f$ are the data observations, \b y is the model output, and \f$ \Sigma \f$ is the noise covariance matrix.
//!
//! \par Isotropic Gaussian
//! The simplest isotopic case of Gaussian likelihood assumes that the covariance matrix is scaled identity,
//! i.e., the noise has no correlation and has the same magnitude for each model output,
//! then inverting the covariance matrix reduces to multiplying by the inverse of the noise variance (magnitude).
//! Also, the sum corresponding to the multiple data samples can be replaced by scaled operation on the data mean (average).
class LikelihoodGaussIsotropic : public TasmanianLikelihood{
public:
    //! \brief Default constructor for convenience, an object constructed with the default cannot be used until \b setData() is called.
    LikelihoodGaussIsotropic(){}
    //! \brief Constructs the class and calls \b setData().
    LikelihoodGaussIsotropic(double variance, const std::vector<double> &data_mean, double num_observe = 1.0){ setData(variance, data_mean, num_observe); }
    //! \brief Default destructor.
    ~LikelihoodGaussIsotropic(){}

    //! \brief Set the noise magnitude (\b varaince) the observed data (\b data_mean) and number of observations (\b num_observe).

    //! The \b variance must be a positive number, the \b data_mean must have the same size as the number of model outputs.
    //! The \b num_observe defaults to 1.0 and while it is not restircted to integer values (since it just multiplies by the inverse of the variance),
    //! \b num_observe must be positive.
    void setData(double variance, const std::vector<double> &data_mean, double num_observe = 1.0);

    //! \brief Compute the likelihood of a set of model outputs.
    void getLikelihood(TypeSamplingForm form, const std::vector<double> &model, std::vector<double> &likely) const;

    //! \brief Returns the size of the \b data_mean vector (for error checking purposes).
    int getNumOuputs() const{ return (int) data.size(); }

private:
    std::vector<double> data;
    double scale;
};


}

#endif
