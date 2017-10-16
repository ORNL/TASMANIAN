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

#include "TasmanianSparseGrid.hpp"

#include "tdrEnumerates.hpp"

namespace TasDREAM{

class BaseUniform{
public:
    BaseUniform();
    virtual ~BaseUniform();

    virtual double getSample01() const = 0;
};

class CppUniformSampler : public BaseUniform{
public:
    CppUniformSampler();
    ~CppUniformSampler();

    double getSample01() const; // returns a smaple uniformly distributed on 01
};

class BasePDF{
public:
    BasePDF();
    virtual ~BasePDF();
    virtual void overwriteBaseUnifrom(const BaseUniform *new_uniform) = 0;

    virtual bool isBoundedBelow() const = 0;
    virtual bool isBoundedAbove() const = 0;

    virtual double getBoundBelow() const = 0;
    virtual double getBoundAbove() const = 0;

    virtual double getSample() const = 0;
    virtual double getDensity(double x) const = 0;
    virtual double getDensityLog(double x) const = 0;

    virtual TypeDistribution getType() const = 0;
};

class UniformPDF : public BasePDF{
public:
    UniformPDF(double lower_bound, double upper_bound);
    ~UniformPDF();
    void overwriteBaseUnifrom(const BaseUniform *new_uniform);

    bool isBoundedBelow() const;
    bool isBoundedAbove() const;

    double getBoundBelow() const;
    double getBoundAbove() const;

    double getSample() const;
    double getDensity(double x) const;
    double getDensityLog(double x) const;

    TypeDistribution getType() const;

private:
    double lower, upper, slope;

    const BaseUniform *core;
    CppUniformSampler unifrom_cpp;
};

class GaussianPDF : public BasePDF{
public:
    GaussianPDF(double gaussian_mean, double gaussian_variance);
    ~GaussianPDF();
    void overwriteBaseUnifrom(const BaseUniform *new_uniform);

    bool isBoundedBelow() const;
    bool isBoundedAbove() const;

    double getBoundBelow() const;
    double getBoundAbove() const;

    double getSample() const;
    double getDensity(double x) const;
    double getDensityLog(double x) const;

    TypeDistribution getType() const;

    void setThreadSave(bool thread_safe);

private:
    double mean, variance, deviation, scale, scaled_precision;
    mutable double saved_value; // reduce the calls to core->getSample01()
    mutable bool saved;
    bool threaded;

    const BaseUniform *core;
    CppUniformSampler unifrom_cpp;
};

class TruncatedGaussianPDF : public BasePDF{ // the PDF is not scaled properly to 1.0, still good to be used as a prior
public:
    TruncatedGaussianPDF(double gaussian_mean, double gaussian_variance, double lower_bound, double upper_bound);
    ~TruncatedGaussianPDF();
    void overwriteBaseUnifrom(const BaseUniform *new_uniform);

    bool isBoundedBelow() const;
    bool isBoundedAbove() const;

    double getBoundBelow() const;
    double getBoundAbove() const;

    double getSample() const;
    double getDensity(double x) const;
    double getDensityLog(double x) const;

    TypeDistribution getType() const;

protected:
    double getStandardTail(double x) const;

private:
    double mean, variance, lower, upper;
    double deviation, scale, scaled_precision;

    const BaseUniform *core;
    CppUniformSampler unifrom_cpp;
};

//class WeibullPDF : public BasePDF{
//public:
//    WeibullPDF(double in_scale, double in_shape, double in_lower);
//    ~WeibullPDF();
//    void overwriteBaseUnifrom(const BaseUniform *new_uniform);
//
//    bool isBoundedBelow() const;
//    bool isBoundedAbove() const;
//
//    double getBoundBelow() const;
//    double getBoundAbove() const;
//
//    double getSample() const;
//    double getDensity(double x) const;
//    double getDensityLog(double x) const;
//
//    TypeDistribution getType() const;
//
//private:
//    double scale, shape, lower; // parameters
//    double ssratio, shapem1, shapeinv; // pre-computed constants
//
//    const BaseUniform *core;
//    CppUniformSampler unifrom_cpp;
//};

class ExponentialPDF : public BasePDF{
public:
    ExponentialPDF(double exponential_rate, double exponential_lower);
    ~ExponentialPDF();
    void overwriteBaseUnifrom(const BaseUniform *new_uniform);

    bool isBoundedBelow() const;
    bool isBoundedAbove() const;

    double getBoundBelow() const;
    double getBoundAbove() const;

    double getSample() const;
    double getDensity(double x) const;
    double getDensityLog(double x) const;

    TypeDistribution getType() const;

private:
    double rate, lower; // parameters

    const BaseUniform *core;
    CppUniformSampler unifrom_cpp;
};

// Gamma must come before Beta because of the sampling algorithm
class GammaPDF : public BasePDF{
public:
    GammaPDF(double lower_bound, double gamma_shape, double gamma_rate); // convection shape = alpha, scale = beta
    ~GammaPDF();
    void overwriteBaseUnifrom(const BaseUniform *new_uniform);

    bool isBoundedBelow() const;
    bool isBoundedAbove() const;

    double getBoundBelow() const;
    double getBoundAbove() const;

    double getSample() const;
    double getDensity(double x) const;
    double getDensityLog(double x) const;

    TypeDistribution getType() const;

private:
    double lower, shape, rate, shapem1, constant, constantLog;

    const BaseUniform *core;
    CppUniformSampler unifrom_cpp;
};

class BetaPDF : public BasePDF{
public:
    BetaPDF(double lower_bound, double upper_bound, double beta_shape_alpha, double beta_shape_beta);
    ~BetaPDF();
    void overwriteBaseUnifrom(const BaseUniform *new_uniform);

    bool isBoundedBelow() const;
    bool isBoundedAbove() const;

    double getBoundBelow() const;
    double getBoundAbove() const;

    double getSample() const;
    double getDensity(double x) const;
    double getDensityLog(double x) const;

    TypeDistribution getType() const;

private:
    double lower, upper, shape_alpha, shape_beta, constant, constantLog;
    GammaPDF *gamma_alpha, *gamma_beta;

    const BaseUniform *core;
    CppUniformSampler unifrom_cpp;
};

class SparseGridDomainToPDF{
public:
    SparseGridDomainToPDF();
    ~SparseGridDomainToPDF();

    static void assumeDefaultPDF(const TasGrid::TasmanianSparseGrid *grid, BasePDF **priors);
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Likelihood Section
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BaseLikelihood{
public:
    BaseLikelihood();
    virtual ~BaseLikelihood();

    virtual double* getLikelihood(int num_model, const double *model, int num_data, const double *data, double *likelihood = 0, bool useLogForm = true) = 0;
};

class GaussianLikelihood : public BaseLikelihood{
public:
    GaussianLikelihood(int outputs, TypeLikelihood likelihood, const double covariance[], int data_entries, const double data[]);
    ~GaussianLikelihood();

    double* getLikelihood(int num_model, const double *model, int num_data = 0, const double *data = 0, double *likelihood = 0, bool useLogForm = true);

private:
    int num_outputs;
    double model_scale;
    TypeLikelihood likely_type;
    double *covariance_cache;
    double *data_cache;
};


//double sample_gamma(double alpha){
//    // George S. Fishman, Sampling from the gamma distribution on a computer, Communications of the ACM, vol. 19, num. 7, pg. 407-409, 1976
//    double alphap = alpha - 1.0;
//    double w, v;
//    do{
//        v = -log(uniform01());
//        w = -log(uniform01());
//    }while(w <= alphap * (v - log(v) - 1.0));
//    return alpha *v;
//}

}


#endif
