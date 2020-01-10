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

#include "tsgDreamLikelihoodCore.hpp"

/*!
 * \internal
 * \file tsgDreamLikelyGaussian.hpp
 * \brief Several likelihood implementations based on Gaussian noise.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAM
 * \endinternal
 */

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
    LikelihoodGaussIsotropic() : scale(0.0){}
    //! \brief Constructs the class and calls \b setData().
    LikelihoodGaussIsotropic(double variance, const std::vector<double> &data_mean, size_t num_observe = 1){ setData(variance, data_mean, num_observe); }
    //! \brief Default destructor.
    ~LikelihoodGaussIsotropic() = default;

    /*!
     * \brief Set the noise magnitude (\b varaince) the observed data (\b data_mean) and number of observations (\b num_observe).
     *
     * Set the parameters of the likelihood.
     * \param variance must be a positive number indicating the magnitude of the noise.
     * \param data_mean must have the same size as the number of model outputs and hold
     *      the average of all measurements.
     * \param num_observe must be a positive integer indicating the number of samples.
     */
    void setData(double variance, const std::vector<double> &data_mean, size_t num_observe = 1);

    //! \brief Compute the likelihood of a set of model outputs.
    void getLikelihood(TypeSamplingForm form, const std::vector<double> &model, std::vector<double> &likely) const override final;

    //! \brief Overload for raw-arrays, for interface purposes mostly, e.g., python.
    void getLikelihood(TypeSamplingForm form, double const model[], int num_samples, double likely[]) const override final;

    //! \brief Returns the size of the \b data_mean vector (for error checking purposes).
    int getNumOutputs() const override{ return (int) data.size(); }

    /*!
     * \brief Writes the data for a portion of the outputs into a stream.
     *
     * The likelihood object does not store the raw inputs to setData(), instead optimized data-structures are used.
     * This method writes either entire likelihood or the optimized data for a portion of the outputs.
     *
     * \param os is the stream where the data will be written.
     * \param outputs_begin is the first output to include in the write process.
     * \param outputs_end is one more than the last output to write,
     *      use -1 to indicate all outputs after \b output_begin.
     *
     * This method is used by the MPI scatter likelihood template.
     */
    void write(std::ostream &os, int outputs_begin = 0, int outputs_end = -1) const{
        if (outputs_end < 0) outputs_end = getNumOutputs();
        outputs_end = std::min(std::max(outputs_begin + 1, outputs_end), getNumOutputs());
        int num_entries = outputs_end - outputs_begin;
        TasGrid::IO::writeNumbers<TasGrid::mode_binary, TasGrid::IO::pad_none>(os, num_entries);
        TasGrid::IO::writeNumbers<TasGrid::mode_binary, TasGrid::IO::pad_none>(os, scale);
        os.write((char*) &data[outputs_begin], num_entries * sizeof(double));
    }

    //! \brief Reads the data from a stream, assumes write() has been used first.
    void read(std::istream &is){
        int num_entries = TasGrid::IO::readNumber<TasGrid::IO::mode_binary_type, int>(is);
        scale = TasGrid::IO::readNumber<TasGrid::IO::mode_binary_type, double>(is);
        data = std::vector<double>((size_t) num_entries);
        TasGrid::IO::readVector<TasGrid::IO::mode_binary_type>(is, data);
    }

private:
    std::vector<double> data;
    double scale;
};

/*!
 * \brief Implements likelihood under the assumption of anisotropic white noise.
 * \ingroup DREAMLikelihood
 *
 * \par Gaussian Likelihood
 * The general formula for Gaussian likelihood is \f$ L(y | d_1 \cdots d_n) = \exp\left( 0.5 sum_{i=1}^n (y - d_i)^T \Sigma^{-1} (y - d_i) \right)\f$
 * where the \f$ d_i \f$ are the data observations, \b y is the model output, and \f$ \Sigma \f$ is the noise covariance matrix.
 *
 * \par Anisotropic Gaussian
 * A more advanced version of the Gaussian likelihood associated each model output
 * with noise of different magnitude. The inversion of the covariance reduces to
 * division of each input by the corresponding magnitude
 * (i.e., there is no matrix inversion).
 * Similarly to the simple case, the sum corresponding to the multiple data samples
 * can be replaced by scaled operation on the data mean (average).
 */
class LikelihoodGaussAnisotropic : public TasmanianLikelihood{
public:
    //! \brief Default constructor for convenience, an object constructed with the default cannot be used until \b setData() is called.
    LikelihoodGaussAnisotropic() = default;
    //! \brief Constructs the class and calls \b setData().
    LikelihoodGaussAnisotropic(std::vector<double> const &variance, std::vector<double> const &data_mean, size_t num_observe = 1){ setData(variance, data_mean, num_observe); }
    //! \brief Default destructor.
    ~LikelihoodGaussAnisotropic() = default;

    /*!
     * \brief Set the noise magnitude (\b variance) the observed data (\b data_mean) and number of observations (\b num_observe).
     *
     * \param variance is a vector with size equal to the number of model outputs.
     *      Each entry represents the noise magnitude (variance) associated with that output.
     * \param data_mean is the average of all available observations of the data.
     * \param num_observe is the number of observations used to compute the \b data_mean.
     */
    void setData(std::vector<double> const &variance, std::vector<double> const &data_mean, size_t num_observe = 1);

    //! \brief Compute the likelihood of a set of model outputs.
    void getLikelihood(TypeSamplingForm form, std::vector<double> const &model, std::vector<double> &likely) const override final;

    //! \brief Overload for raw-arrays, for interface purposes mostly, e.g., python.
    void getLikelihood(TypeSamplingForm form, double const model[], int num_samples, double likely[]) const override final;

    //! \brief Returns the size of the \b data_mean vector (for error checking purposes).
    int getNumOutputs() const override{ return (int) noise_variance.size(); }

    /*!
     * \brief Writes the data for a portion of the outputs into a stream.
     *
     * See LikelihoodGaussIsotropic::write().
     */
    void write(std::ostream &os, int outputs_begin = 0, int outputs_end = -1) const{
        if (outputs_end < 0) outputs_end = getNumOutputs();
        outputs_end = std::min(std::max(outputs_begin + 1, outputs_end), getNumOutputs());
        int num_entries = outputs_end - outputs_begin;
        TasGrid::IO::writeNumbers<TasGrid::mode_binary, TasGrid::IO::pad_none>(os, num_entries);
        os.write((char*) &data_by_variance[outputs_begin], num_entries * sizeof(double));
        os.write((char*) &noise_variance[outputs_begin], num_entries * sizeof(double));
    }

    //! \brief Reads the data from a stream, assumes write() has been used first.
    void read(std::istream &is){
        int num_entries = TasGrid::IO::readNumber<TasGrid::IO::mode_binary_type, int>(is);
        data_by_variance = std::vector<double>((size_t) num_entries);
        noise_variance   = std::vector<double>((size_t) num_entries);
        TasGrid::IO::readVector<TasGrid::IO::mode_binary_type>(is, data_by_variance);
        TasGrid::IO::readVector<TasGrid::IO::mode_binary_type>(is, noise_variance);
    }

private:
    std::vector<double> data_by_variance;
    std::vector<double> noise_variance;
};


}

#endif
