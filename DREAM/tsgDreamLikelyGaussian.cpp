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
#include "tsgDreamInternalBlas.hpp"

namespace TasDREAM{

void LikelihoodGaussIsotropic::setData(double variance, const std::vector<double> &data_mean, size_t num_observe){
    if (variance <= 0.0) throw std::runtime_error("ERROR: LikelihoodGaussIsotropic, should have positive varience.");
    if (data_mean.empty()) throw std::runtime_error("ERROR: LikelihoodGaussIsotropic, emptry data vector.");

    data = data_mean;
    scale = -0.5 * double(num_observe) / variance;
}

void LikelihoodGaussIsotropic::getLikelihood(TypeSamplingForm form, const std::vector<double> &model, std::vector<double> &likely) const{
    auto im = model.begin();
    #ifdef Tasmanian_ENABLE_BLAS
    int num_outputs = (int) data.size();
    int num_points = (int) (model.size() / data.size());

    for(auto &l : likely){
        l = scale * TasBLAS::dnrm2squared(num_outputs, &*im);
        std::advance(im, num_outputs);
    }
    TasBLAS::dgemtv(num_outputs, num_points, model.data(), data.data(), likely.data(), -2.0 * scale, 1.0);
    #else
    size_t num_outputs = data.size();
    for(auto &l : likely){
        l = scale * (std::inner_product(im, im + num_outputs, im, 0.0) - 2.0 * std::inner_product(im, im + num_outputs, data.data(), 0.0));
        std::advance(im, num_outputs);
    }
    #endif
    if (form == regform) for(auto &l : likely) l = std::exp(l);
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
    auto im = model.begin();
    #ifdef Tasmanian_ENABLE_BLAS
    int num_outputs = (int) data_by_variance.size();
    int num_points = (int) (model.size() / data_by_variance.size());

    for(auto &l : likely){
        l = 0.0;
        for(auto const &v : noise_variance){
            l += *im * *im * v;
            im++;
        }
    }
    TasBLAS::dgemtv(num_outputs, num_points, model.data(), data_by_variance.data(), likely.data(), -2.0, 1.0);
    #else
    for(auto &l : likely){
        l = 0.0;
        auto id = data_by_variance.begin();
        for(auto const &v : noise_variance){
            l += *im * *im * v - 2.0 * *im * *id;
            im++; id++;
        }
    }
    #endif
    if (form == regform) for(auto &l : likely) l = std::exp(l);
}

}

#endif
