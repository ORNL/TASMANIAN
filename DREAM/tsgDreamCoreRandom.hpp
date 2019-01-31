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

#ifndef __TASMANIAN_DREAM_CORE_RANDOM_HPP
#define __TASMANIAN_DREAM_CORE_RANDOM_HPP

//! \file tsgDreamCoreRandom.hpp
//! \brief Core random sampling methods
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAM
//! \internal
//!
//! Implements several methods for random sampling and defines specific probability density functions.

#include <math.h>

namespace TasDREAM{

//! \internal
//! \brief Default random sampler, using \b rand() divided by \b RAND_MAX
//! \ingroup TasmanianDREAM
    
//! Generates random numbers uniformly distributed in (0, 1), uses the \b rand() command.
inline double tsgCoreUniform01(){ return ((double) rand()) / ((double) RAND_MAX); }

//! \brief Add a correction to every entry in \b x, use uniform samples over (-\b magnitude, \b magnitude).
//! \ingroup TasmanianDREAM

//! The function \b get_random01() returns random numbers distributed over (0, 1).
inline void applyUniformUpdate(std::vector<double> &x, double magnitude, std::function<double(void)> get_random01 = tsgCoreUniform01){
    if (magnitude == 0.0) return;
    for(auto &v : x) v += magnitude * (2.0 * get_random01() -1.0);
}

//! \brief  Add a correction to every entry in \b x, sue Gaussian distribution with zero mean and standard deviation equal to \b magnitude.
//! \ingroup TasmanianDREAM

//! The function \b get_random01() returns random numbers distributed over (0, 1).
//! Gaussian numbers are generated using the Box-Muller algorithm.
inline void applyGaussianUpdate(std::vector<double> &x, double magnitude, std::function<double(void)> get_random01 = tsgCoreUniform01){
    if (magnitude == 0.0) return;
    bool tictoc = false;
    double g = 0.0;
    for(auto &v : x){
        tictoc = !tictoc;
        if (tictoc){
            double r = magnitude * sqrt(-2.0 * log(get_random01())), t = 2.0 * M_PI * get_random01(); // radius and angle
            v += r * cos(t);
            g = r * sin(t);
        }else{
            v += g;
        }
    }
}

//! \brief Generate uniform random samples in the hypercube defined by \b lower and \b upper limits.
//! \ingroup TasmanianDREAM

//! The size of the \b lower and \b upper must match.
//! The output vector \b x will be resized to match \b num_samples times \b upper.size(), and the values will be overwritten.
//! The function \b get_random01() returns random numbers distributed over (0, 1).
inline void genUniformSamples(const std::vector<double> &lower, const std::vector<double> &upper, int num_samples, std::vector<double> &x, std::function<double(void)> get_random01 = tsgCoreUniform01){
    if (lower.size() != upper.size()) throw std::runtime_error("ERROR: genUniformSamples() requires lower and upper vectors with matching size.");
    if (x.size() != lower.size() * num_samples) x.resize(lower.size() * num_samples);
    for(auto &v : x) v = get_random01();
    
    std::vector<double> length(lower.size());
    std::transform(lower.begin(), lower.end(), upper.begin(), length.begin(), [&](double l, double u)->double{ return (u - l); });
    
    auto ix = x.begin();
    while(ix != x.end()){
        auto ilow = lower.begin();
        for(auto l : length){
            *ix *= l;
            *ix++ += *ilow++;
        }
    }
}

}

#endif
