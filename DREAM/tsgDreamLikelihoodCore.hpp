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

#ifndef __TASMANIAN_DREAM_LIKELY_CORE_HPP
#define __TASMANIAN_DREAM_LIKELY_CORE_HPP

#include "tsgDreamEnumerates.hpp"

/*!
 * \internal
 * \file tsgDreamLikelihoodCore.hpp
 * \brief The interface mother-class for the likelihood classes.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAM
 *
 * Defines the TasmanianLikelihood class which implements the relation between model outputs and measurement data.
 * \endinternal
 */

/*!
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMLikelihood Likelihood definitions
 *
 * Classes that define the likelihood relation between an observed data and a model output.
 */

namespace TasDREAM{

//! \brief Interface for the likelihood classes.
//! \ingroup DREAMLikelihood

//! \par Likelihood
//! In the framework of Bayesian inference the likelihood is a measure of how likely is a specific outcome (i.e., model output), given some observed data.
//! The likelihood class (usually) contains the data and implements a specific formula that measure the discrepancy.
//! The \b TasmanianLikelihood class is virtual, but an inherited class is required for the \b SampleDREAMPosterior() methods.
//!
//! \par Included Likelihood Formulas
//! The chose of likelihood is very problem dependent and Tasmanian implements only a few commonly used cases,
//! mostly relying on the assumption that the data is contaminated with white noise.
//! However, the user can implement any other likelihood by simply inheriting from this class.
class TasmanianLikelihood{
public:
    //! \brief Empty default constructor.
    TasmanianLikelihood(){}
    //! \brief Empty virtual destructor.
    virtual ~TasmanianLikelihood(){}

    //! \brief Purely virtual method used by \b SampleDREAMPosterior(), computes the likelihood of multiple model values.

    //! The \b model vector is the same as the output of the model lambda in \b SampleDREAMPosterior() or the \b TasGrid::TasmanianSparseGrid::evaluateBatch(),
    //! the model realizations are stored contiguously in strides of length equal to the number of model outputs.
    //! The \b likely vector is pre-allocated with size matching the number of model realizations under considerations (no resize is needed),
    //! this function must populate the \b likely entries with the corresponding values of the likelihood.
    //! Note that model.size() / likely.size() will divide evenly and will equal the number of model realizations.
    virtual void getLikelihood(TypeSamplingForm form, const std::vector<double> &model, std::vector<double> &likely) const = 0;

    //! \brief Overload for raw-arrays, for interface purposes mostly, never called from C++ directly.
    virtual void getLikelihood(TypeSamplingForm form, double const model[], int num_samples, double likely[]) const = 0;

    //! \brief Return the number of expected model outputs.
    virtual int getNumOutputs() const = 0;

    //! \brief Automatically convert the likelihood into input for TasDREAM::posterior().
    virtual operator std::function<void(TypeSamplingForm, const std::vector<double> &, std::vector<double> &)>() const{
        return [&](TypeSamplingForm form, const std::vector<double> &model, std::vector<double> &likely)->void{
            getLikelihood(form, model, likely);
        };
    };

};

}

#endif
