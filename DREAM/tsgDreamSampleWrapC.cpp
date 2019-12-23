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

#ifndef __TASMANIAN_DREAM_SAMPLE_WRAPC_CPP
#define __TASMANIAN_DREAM_SAMPLE_WRAPC_CPP

#include "tsgDreamSample.hpp"

namespace TasDREAM{

using tsg_dream_pdf     = void (*)(int, int, const double[], double[], int*);
using tsg_dream_domain  =  int (*)(int, const double[]);
using tsg_dream_iupdate = void (*)(int, double[], int*);
using tsg_dream_dupdate = double (*)();
using tsg_dream_random  = double (*)();

std::function<bool(const std::vector<double> &x)>
getSpecifiedDomain(int num_dimensions, void *domain_grid, double *domain_lower, double *domain_upper, tsg_dream_domain domain_callback){
    if (domain_grid != nullptr){
        return reinterpret_cast<TasGrid::TasmanianSparseGrid*>(domain_grid)->getDomainInside();
    }else if (domain_upper != nullptr){
        return hypercube(std::vector<double>(domain_lower, domain_lower + num_dimensions),
                         std::vector<double>(domain_upper, domain_upper + num_dimensions));
    }else{
        return [=](std::vector<double> const &x)->
        bool{
            return (domain_callback((int) x.size(), x.data()) != 0);
        };
    }
}

std::function<double(void)>
getSpecifiedDifferentialUpdate(int dupdate_percent, tsg_dream_dupdate dupdate_callback){
    if (dupdate_percent >= 0){
        return [=]()->double{ return double(dupdate_percent) / 100.0; };
    }else{
        return [=]()->double{ return dupdate_callback(); };
    }
}


extern "C"{

void tsgGenUniformSamples(int num_dimensions, int num_samples, double const lower[], double const upper[],
                          const char* random_type, int random_seed, tsg_dream_random random_callback, double *samples){

    std::minstd_rand park_miller((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    srand((unsigned int) ((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed));
    std::string rtype(random_type);

    auto randgen = [&]()->
    std::function<double(void)>{
        if (rtype == "default"){
            return [&]()->double{ return tsgCoreUniform01(); };
        }else if (rtype == "minstd_rand"){
            return [&]()->double{ return unif(park_miller); };
        }else{
            return [&]()->double{ return random_callback(); };
        }
    }();

    std::vector<double> result = TasDREAM::genUniformSamples(Utils::copyArray(lower, num_dimensions),
                                                             Utils::copyArray(upper, num_dimensions),
                                                             num_samples, randgen);
    std::copy(result.begin(), result.end(), samples);
}

void tsgGenGaussianSamples(int num_dimensions, int num_samples, double const mean[], double const deviation[],
                           const char* random_type, int random_seed, tsg_dream_random random_callback, double *samples){

    std::minstd_rand park_miller((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    srand((unsigned int) ((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed));
    std::string rtype(random_type);

    auto randgen = [&]()->
    std::function<double(void)>{
        if (rtype == "default"){
            return [&]()->double{ return tsgCoreUniform01(); };
        }else if (rtype == "minstd_rand"){
            return [&]()->double{ return unif(park_miller); };
        }else{
            return [&]()->double{ return random_callback(); };
        }
    }();

    std::vector<double> result = TasDREAM::genGaussianSamples(Utils::copyArray(mean, num_dimensions),
                                                              Utils::copyArray(deviation, num_dimensions),
                                                              num_samples, randgen);
    std::copy(result.begin(), result.end(), samples);
}

void tsgDreamSample(int form,
                    int num_burnup, int num_collect,
                    tsg_dream_pdf distribution,
                    void* state_pntr,
                    void *domain_grid, double domain_lower[], double dommain_upper[], tsg_dream_domain domain_callback,
                    const char* iupdate_type, double iupdate_magnitude, tsg_dream_iupdate iupdate_callback,
                    int dupdate_percent, tsg_dream_dupdate dupdate_callback,
                    const char* random_type, int random_seed, tsg_dream_random random_callback, int *err){
    *err = 1;
    TasmanianDREAM& state = *reinterpret_cast<TasmanianDREAM*>(state_pntr);

    int num_dimensions = (int) state.getNumDimensions();

    auto domain = getSpecifiedDomain(num_dimensions, domain_grid, domain_lower, dommain_upper, domain_callback);

    TypeDistribution dist = IO::getDistributionString(iupdate_type);

    auto diff_update = getSpecifiedDifferentialUpdate(dupdate_percent, dupdate_callback);

    std::minstd_rand park_miller((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    srand((unsigned int) ((random_seed == -1) ? static_cast<long unsigned>(std::time(nullptr)) : random_seed));
    std::string rtype(random_type);

    auto randgen = [&]()->
    std::function<double(void)>{
        if (rtype == "default"){
            return [&]()->double{ return tsgCoreUniform01(); };
        }else if (rtype == "minstd_rand"){
            return [&]()->double{ return unif(park_miller); };
        }else{
            return [&]()->double{ return random_callback(); };
        }
    }();

    try{
        if (dist == dist_null){
            if (IO::intToForm(form) == regform){
                SampleDREAM<regform>(num_burnup, num_collect, [&](const std::vector<double> &candidates, std::vector<double> &values)->
                void{
                    int num_samples = (int) candidates.size() / num_dimensions;
                    int error_code = 0;
                    distribution(num_samples, num_dimensions, candidates.data(), values.data(), &error_code);
                    if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgDreamSample()");
                }, domain, state, [&](std::vector<double> &x)->
                void{
                    int error_code = 0;
                    iupdate_callback((int) x.size(), x.data(), &error_code);
                    if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgDreamSample()");
                }, diff_update, randgen);
            }else{
                SampleDREAM<logform>(num_burnup, num_collect, [&](const std::vector<double> &candidates, std::vector<double> &values)->
                void{
                    int num_samples = (int) candidates.size() / num_dimensions;
                    int error_code = 0;
                    distribution(num_samples, num_dimensions, candidates.data(), values.data(), &error_code);
                    if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgDreamSample()");
                }, domain, state, [&](std::vector<double> &x)->
                void{
                    int error_code = 0;
                    iupdate_callback((int) x.size(), x.data(), &error_code);
                    if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgDreamSample()");
                }, diff_update, randgen);
            }
        }else{
            if (IO::intToForm(form) == regform){
                SampleDREAM<regform>(num_burnup, num_collect, [&](const std::vector<double> &candidates, std::vector<double> &values)->
                void{
                    int num_samples = (int) candidates.size() / num_dimensions;
                    int error_code = 0;
                    distribution(num_samples, num_dimensions, candidates.data(), values.data(), &error_code);
                    if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgDreamSample()");
                }, domain, state, dist, iupdate_magnitude, diff_update, randgen);
            }else{
                SampleDREAM<logform>(num_burnup, num_collect, [&](const std::vector<double> &candidates, std::vector<double> &values)->
                void{
                    int num_samples = (int) candidates.size() / num_dimensions;
                    int error_code = 0;
                    distribution(num_samples, num_dimensions, candidates.data(), values.data(), &error_code);
                    if (error_code != 0) throw std::runtime_error("The Python callback returned an error in analysis()");
                }, domain, state, dist, iupdate_magnitude, diff_update, randgen);
            }
        }
        *err = 0; // success
    }catch(std::runtime_error &){} // *err will remain 1
}

}

}

#endif
