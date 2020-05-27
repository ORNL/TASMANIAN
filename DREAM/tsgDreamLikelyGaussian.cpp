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

#ifndef __TASMANIAN_DREAM_LIKELY_GAUSS_CPP
#define __TASMANIAN_DREAM_LIKELY_GAUSS_CPP

#include "tsgDreamLikelyGaussian.hpp"
#ifdef Tasmanian_ENABLE_BLAS
#include "tsgBlasWrappers.hpp"
#endif

namespace TasDREAM{

void LikelihoodGaussIsotropic::setData(double variance, const std::vector<double> &data_mean, size_t num_observe){
    if (variance <= 0.0) throw std::runtime_error("ERROR: LikelihoodGaussIsotropic, should have positive varience.");
    if (data_mean.empty()) throw std::runtime_error("ERROR: LikelihoodGaussIsotropic, emptry data vector.");

    data = data_mean;
    scale = -0.5 * double(num_observe) / variance;
}

void LikelihoodGaussIsotropic::getLikelihood(TypeSamplingForm form, const std::vector<double> &model, std::vector<double> &likely) const{
    getLikelihood(form, model.data(), (int) (model.size() / data.size()), likely.data());
}
void LikelihoodGaussIsotropic::getLikelihood(TypeSamplingForm form, double const model[], int num_samples, double likely[]) const{
    int num_outputs = (int) data.size();
    Utils::Wrapper2D<const double> wrapped_model(num_outputs, model);
    #ifdef Tasmanian_ENABLE_BLAS
    for(int i=0; i<num_samples; i++)
        likely[i] = scale * TasBLAS::norm2_2(num_outputs, wrapped_model.getStrip(i));
    TasBLAS::gemv('T', num_outputs, num_samples, -2.0 * scale, model, num_outputs, data.data(), 1, 1.0, likely, 1);
    #else
    for(int i=0; i<num_samples; i++)
        likely[i] = scale * (std::inner_product(wrapped_model.getStrip(i), wrapped_model.getStrip(i) + num_outputs, wrapped_model.getStrip(i), 0.0)
                             - 2.0 * std::inner_product(wrapped_model.getStrip(i), wrapped_model.getStrip(i) + num_outputs, data.data(), 0.0));
    #endif
    if (form == regform) for(int i=0; i<num_samples; i++) likely[i] = std::exp(likely[i]);
}

void LikelihoodGaussAnisotropic::setData(std::vector<double> const &variance, std::vector<double> const &data_mean, size_t num_observe){
    if (variance.size() != data_mean.size()) throw std::invalid_argument("ERROR: LikelihoodGaussAnisotropic, should have variance and data with same size.");

    double scale = -0.5 * double(num_observe);
    noise_variance = std::vector<double>(variance.size());
    data_by_variance = std::vector<double>(variance.size());
    for(size_t i=0; i<variance.size(); i++){
        noise_variance[i] = scale / variance[i];
        data_by_variance[i] = scale * data_mean[i] / variance[i];
    }
}

void LikelihoodGaussAnisotropic::getLikelihood(TypeSamplingForm form, std::vector<double> const &model, std::vector<double> &likely) const{
    getLikelihood(form, model.data(), (int) (model.size() / data_by_variance.size()), likely.data());
}

void LikelihoodGaussAnisotropic::getLikelihood(TypeSamplingForm form, double const model[], int num_samples, double likely[]) const{
    int num_outputs = (int) data_by_variance.size();
    Utils::Wrapper2D<const double> wrapped_model(num_outputs, model);
    #ifdef Tasmanian_ENABLE_BLAS
    for(int i=0; i<num_samples; i++){
        const double *sample = wrapped_model.getStrip(i);
        likely[i] = 0.0;
        for(int k=0; k<num_outputs; k++){
            likely[i] += sample[k] * sample[k] * noise_variance[k];
        }
    }
    TasBLAS::gemv('T', num_outputs, num_samples, -2.0, model, num_outputs, data_by_variance.data(), 1, 1.0, likely, 1);
    #else
    for(int i=0; i<num_samples; i++){
        const double *sample = wrapped_model.getStrip(i);
        likely[i] = 0.0;
        for(int k=0; k<num_outputs; k++){
            likely[i] += sample[k] * sample[k] * noise_variance[k] - 2.0 * sample[k] * data_by_variance[k];
        }
    }
    #endif
    if (form == regform) for(int i=0; i<num_samples; i++) likely[i] = std::exp(likely[i]);
}

extern "C"{ // for python purposes
void *tsgMakeLikelihoodGaussIsotropic(int num_outputs, double variance, double const data[], int num_samples){
    return (void*) new LikelihoodGaussIsotropic(variance, std::vector<double>(data, data + num_outputs), (size_t) num_samples);
}
void *tsgMakeLikelihoodGaussAnisotropic(int num_outputs, double const variance[], double const data[], int num_samples){
    return (void*) new LikelihoodGaussAnisotropic(std::vector<double>(variance, variance + num_outputs),
                                                  std::vector<double>(data, data + num_outputs), (size_t) num_samples);
}
void tsgGetLikelihood(void *likelihood, int form, double const model[], int num_samples, double likely[]){
    ((TasmanianLikelihood*) likelihood)->getLikelihood(IO::intToForm(form), model, num_samples, likely);
}
int tsgGetNumOutputsLikelihood(void *likelihood){
    return ((TasmanianLikelihood*) likelihood)->getNumOutputs();
}
void tsgDeleteLikelihood(void *likelihood){
    delete ((TasmanianLikelihood*) likelihood);
}

}

}

#endif
