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

#ifndef __TASMANIAN_DREAM_CORE_PDF_CPP
#define __TASMANIAN_DREAM_CORE_PDF_CPP

#include "tdrCorePDF.hpp"

#include "tsgHiddenExternals.hpp"

namespace TasDREAM{

BaseUniform::BaseUniform(){}
BaseUniform::~BaseUniform(){}

CppUniformSampler::CppUniformSampler(){}
CppUniformSampler::~CppUniformSampler(){}
double CppUniformSampler::getSample01() const{ return ((double) rand()) / ((double) RAND_MAX); }


//////////////////////////////////////////////////////////////////////////////////////////////////
//
// Build in Probability Distributions
//
//////////////////////////////////////////////////////////////////////////////////////////////////

BasePDF::BasePDF(){}
BasePDF::~BasePDF(){}

UniformPDF::UniformPDF(double lower_bound, double upper_bound) : lower(lower_bound), upper(upper_bound) {  core = &unifrom_cpp;  slope = upper - lower; }
UniformPDF::~UniformPDF(){ core = 0; }
void UniformPDF::overwriteBaseUnifrom(const BaseUniform *new_uniform){ core = new_uniform; }
bool UniformPDF::isBoundedBelow() const{ return true; }
bool UniformPDF::isBoundedAbove() const{ return true; }
double UniformPDF::getBoundBelow() const{ return lower; }
double UniformPDF::getBoundAbove() const{ return upper; }
double UniformPDF::getSample() const{ return lower + core->getSample01() * slope; }
double UniformPDF::getDensity(double) const{ return 1.0 / slope; }
double UniformPDF::getDensityLog(double) const{ return -log(slope); }
TypeDistribution UniformPDF::getType() const{ return dist_uniform; }


GaussianPDF::GaussianPDF(double gaussian_mean, double gaussian_variance) : mean(gaussian_mean), variance(gaussian_variance), threaded(false)
{
    core = &unifrom_cpp;
    deviation = sqrt(variance);
    scale = 1.0 / sqrt(2.0 * M_PI * variance);
    scaled_precision = -1.0 / (2.0 * variance);
    saved_value = 0.0; saved = false;
}
GaussianPDF::~GaussianPDF(){ core = 0; }
void GaussianPDF::overwriteBaseUnifrom(const BaseUniform *new_uniform){ core = new_uniform; }
bool GaussianPDF::isBoundedBelow() const{ return false; }
bool GaussianPDF::isBoundedAbove() const{ return false; }
double GaussianPDF::getBoundBelow() const{ return 0.0; }
double GaussianPDF::getBoundAbove() const{ return 0.0; }
double GaussianPDF::getSample() const{
    if (threaded) return mean + deviation * sqrt(-2.0 * log(core->getSample01())) * cos(2.0 * M_PI * core->getSample01());
    saved = !saved;
    if (saved){ // actually not-saved
        double uniform1 = core->getSample01();
        double uniform2 = core->getSample01();
        saved_value = mean + deviation * sqrt(-2.0 * log(uniform1)) * sin(2.0 * M_PI * uniform2);
        return mean + deviation * sqrt(-2.0 * log(uniform1)) * cos(2.0 * M_PI * uniform2);
    }else{
        return saved_value;
    }
}
double GaussianPDF::getDensity(double x) const{ return scale * exp(scaled_precision*(x-mean)*(x-mean)); }
double GaussianPDF::getDensityLog(double x) const{ return log(scale) + scaled_precision*(x-mean)*(x-mean); }
TypeDistribution GaussianPDF::getType() const{ return dist_gaussian; }
void GaussianPDF::setThreadSave(bool thread_safe){ threaded = thread_safe; }

TruncatedGaussianPDF::TruncatedGaussianPDF(double gaussian_mean, double gaussian_variance, double lower_bound, double upper_bound)
    : mean(gaussian_mean), variance(gaussian_variance), lower(lower_bound), upper(upper_bound)
{
    core = &unifrom_cpp;
    deviation = sqrt(variance);
    scale = 1.0 / sqrt(2.0 * M_PI * variance); // correct the scale by 1 / ( 1 - tail_left - tail_right )
    scaled_precision = -1.0 / (2.0 * variance);
    double tail_left, tail_right;
    if (mean > lower){
        tail_left = deviation * getStandardTail((lower - mean) / deviation);
    }else{
        tail_left = 1.0 - deviation * getStandardTail((mean - lower) / deviation);
    }
    if (mean < upper){
        tail_right = deviation * getStandardTail((mean - upper) / deviation);
    }else{
        tail_right = 1.0 - deviation * getStandardTail((upper - mean) / deviation);
    }
    scale /= (1.0 - tail_left - tail_right);
}
TruncatedGaussianPDF::~TruncatedGaussianPDF(){  core = 0; }
void TruncatedGaussianPDF::overwriteBaseUnifrom(const BaseUniform *new_uniform){ core = new_uniform; }
bool TruncatedGaussianPDF::isBoundedBelow() const{ return true; }
bool TruncatedGaussianPDF::isBoundedAbove() const{ return true; }
double TruncatedGaussianPDF::getBoundBelow() const{  return lower; }
double TruncatedGaussianPDF::getBoundAbove() const{  return upper; }
double TruncatedGaussianPDF::getSample() const{
    double candidate = mean + deviation * sqrt(-2.0 * log(core->getSample01())) * cos(2.0 * M_PI * core->getSample01());
    while ((candidate < lower) || (candidate > upper)){
        candidate = mean + deviation * sqrt(-2.0 * log(core->getSample01())) * cos(2.0 * M_PI * core->getSample01());
    }
    return candidate;
}
double TruncatedGaussianPDF::getDensity(double x) const{
    return scale * exp(scaled_precision*(x-mean)*(x-mean));
}
double TruncatedGaussianPDF::getDensityLog(double x) const{
    return log(scale) + scaled_precision*(x-mean)*(x-mean);
}
TypeDistribution TruncatedGaussianPDF::getType() const{ return dist_truncated_gaussian; }
double TruncatedGaussianPDF::getStandardTail(double x) const{
    double t = 1.0 / (1.0 + 0.2316419 * x);
    double cdf = 1.330274439 * t - 1.821255978;
    cdf *= t;
    cdf += 1.781477937;
    cdf *= t;
    cdf -= 0.356563782;
    cdf *= t;
    cdf += 0.319381530;
    cdf *= t;
    return t * exp(-0.5 * x * x) / sqrt(2.0 * M_PI);
}


//WeibullPDF::WeibullPDF(double in_scale, double in_shape, double in_lower) : scale(in_scale), shape(in_shape), lower(in_lower)
//{ core = &unifrom_cpp; ssratio = shape / scale; shapem1 = shape-1.0; shapeinv = 1.0 / shape; }
//WeibullPDF::~WeibullPDF(){  core = 0; }
//void WeibullPDF::overwriteBaseUnifrom(const BaseUniform *new_uniform){ core = new_uniform; }
//bool WeibullPDF::isBoundedBelow() const{ return true; }
//bool WeibullPDF::isBoundedAbove() const{ return false; }
//double WeibullPDF::getBoundBelow() const{  return lower;  }
//double WeibullPDF::getBoundAbove() const{  return 0.0;  }
//double WeibullPDF::getSample() const{  return lower + scale * pow(-log(1 - core->getSample01()), shapeinv);  }
//double WeibullPDF::getDensity(double x) const{
//    double norm_x = (x - lower) / scale;
//    //return (norm_x <= 0) ? 0.0 : ssratio * pow(norm_x, shapem1) * exp(-pow(norm_x, shape));
//    return ssratio * pow(norm_x, shapem1) * exp(-pow(norm_x, shape));
//}
//double WeibullPDF::getDensityLog(double x) const{
//    double norm_x = (x - lower) / scale;
//    return log(ssratio) + shapem1 * log(norm_x) - pow(norm_x, shape);
//}
//TypeDistribution WeibullPDF::getType() const{ return dist_weibull; }


ExponentialPDF::ExponentialPDF(double exponential_rate, double lower_bound) : rate(exponential_rate), lower(lower_bound)
{ core = &unifrom_cpp; }
ExponentialPDF::~ExponentialPDF(){  core = 0; }
void ExponentialPDF::overwriteBaseUnifrom(const BaseUniform *new_uniform){ core = new_uniform; }
bool ExponentialPDF::isBoundedBelow() const{return true; }
bool ExponentialPDF::isBoundedAbove() const{return false; }
double ExponentialPDF::getBoundBelow() const{ return lower; }
double ExponentialPDF::getBoundAbove() const{ return 0.0; }
double ExponentialPDF::getSample() const{ return lower + (log(core->getSample01()) - log(rate)) / rate; }
double ExponentialPDF::getDensity(double x) const{
    double norm_x = rate * (lower - x);
    return rate * exp(norm_x);
}
double ExponentialPDF::getDensityLog(double x) const{
    double norm_x = rate * (lower - x)  ;
    return log(rate) + norm_x;
}
TypeDistribution ExponentialPDF::getType() const{ return dist_exponential; }


GammaPDF::GammaPDF(double lower_bound, double gamma_shape, double gamma_rate) : lower(lower_bound), shape(gamma_shape), rate(gamma_rate){
    core = &unifrom_cpp;
    shapem1 = shape - 1.0;
    constant = pow(rate, shape) / tgamma(shape);
    constantLog = shape * log(rate) - lgamma(shape);
}
GammaPDF::~GammaPDF(){}
void GammaPDF::overwriteBaseUnifrom(const BaseUniform *new_uniform){ core = new_uniform; }
bool GammaPDF::isBoundedBelow() const{ return true; }
bool GammaPDF::isBoundedAbove() const{ return false; }
double GammaPDF::getBoundBelow() const{ return lower; }
double GammaPDF::getBoundAbove() const{ return 0.0; }
double GammaPDF::getSample() const{
    double v = -log(core->getSample01());
    double w = -log(core->getSample01());

    while(w <= shapem1 * (v - log(v) - 1.0)){
        v = -log(core->getSample01());
        w = -log(core->getSample01());
    }

    return lower + shape * v / rate;
}
double GammaPDF::getDensity(double x) const{ return constant * pow(x - lower, shapem1) * exp(-rate * (x - lower)); }
double GammaPDF::getDensityLog(double x) const{ return constantLog + shapem1 * log(x-lower) - rate * (x - lower); }
TypeDistribution GammaPDF::getType() const{ return dist_gamma; }


BetaPDF::BetaPDF(double lower_bound, double upper_bound, double beta_shape_alpha, double beta_shape_beta) :
                lower(lower_bound), upper(upper_bound), shape_alpha(beta_shape_alpha), shape_beta(beta_shape_beta){
    constant = tgamma(shape_alpha + shape_beta) / (tgamma(shape_alpha) * tgamma(shape_beta));
    constantLog = lgamma(shape_alpha + shape_beta) - lgamma(shape_alpha) - lgamma(shape_beta);
    gamma_alpha = new GammaPDF(0, shape_alpha, 1.0);
    gamma_beta = new GammaPDF(0, shape_beta, 1.0);
}
BetaPDF::~BetaPDF(){
    if (gamma_alpha != 0){ delete gamma_alpha; gamma_alpha = 0; }
    if (gamma_beta != 0){ delete gamma_beta; gamma_beta = 0; }
}
void BetaPDF::overwriteBaseUnifrom(const BaseUniform *new_uniform){
    core = new_uniform;
    gamma_alpha->overwriteBaseUnifrom(core);
    gamma_beta->overwriteBaseUnifrom(core);
}
bool BetaPDF::isBoundedBelow() const{ return true; }
bool BetaPDF::isBoundedAbove() const{ return true; }
double BetaPDF::getBoundBelow() const{ return lower; }
double BetaPDF::getBoundAbove() const{ return upper; }
double BetaPDF::getSample() const{
    double X = gamma_alpha->getSample();
    double Y = gamma_beta->getSample();
    return lower + (upper - lower) * X / (X + Y);
}
double BetaPDF::getDensity(double x) const{
    double x_canonical = (x - lower) / (upper - lower);
    return constant * pow(x_canonical, shape_alpha - 1.0) * pow(1.0 - x_canonical, shape_beta - 1.0);
}
double BetaPDF::getDensityLog(double x) const{
    double x_canonical = (x - lower) / (upper - lower);
    return constantLog + (shape_alpha - 1.0) * log(x_canonical) + (shape_beta - 1.0) * log(1.0 - x_canonical);
}
TypeDistribution BetaPDF::getType() const{ return dist_beta; }

void SparseGridDomainToPDF::assumeDefaultPDF(const TasGrid::TasmanianSparseGrid *grid, std::vector<BasePDF*> &priors){
    int num_dimensions = grid->getNumDimensions();
    if (priors.size() < (size_t) num_dimensions) priors.resize(num_dimensions);

    TasGrid::TypeOneDRule rule = grid->getRule();

    std::vector<double> a, b;
    if (grid->isSetDomainTransfrom()){
        grid->getDomainTransform(a, b);
    }

    if ((rule == TasGrid::rule_gausslaguerre) || (rule == TasGrid::rule_gausslaguerreodd)){ // Exponential distribution
        if (grid->isSetDomainTransfrom()){
            for(int i=0; i<num_dimensions; i++) priors[i] = new ExponentialPDF(b[i], a[i]);
        }else{
            for(int i=0; i<num_dimensions; i++) priors[i] = new ExponentialPDF(1.0, 0.0);
        }
    }else if ((rule == TasGrid::rule_gausshermite) || (rule == TasGrid::rule_gausshermiteodd)){
        if (grid->isSetDomainTransfrom()){
            for(int i=0; i<num_dimensions; i++) priors[i] = new GaussianPDF(a[i], b[i]);
        }else{
            for(int i=0; i<num_dimensions; i++) priors[i] = new GaussianPDF(0.0, 1.0);
        }
    }else if ((rule == TasGrid::rule_localp0) || (rule == TasGrid::rule_clenshawcurtis0)){
        if (grid->isSetDomainTransfrom()){
            for(int i=0; i<num_dimensions; i++) priors[i] = new BetaPDF(a[i], b[i], 2.0, 2.0);
        }else{
            for(int i=0; i<num_dimensions; i++) priors[i] = new BetaPDF(-1.0, 1.0, 2.0, 2.0);
        }
    }else{ // no clue what to do with custom rule, must be overwritten later
        if (grid->isSetDomainTransfrom()){
            for(int i=0; i<num_dimensions; i++) priors[i] = new UniformPDF(a[i], b[i]);
        }else{
            for(int i=0; i<num_dimensions; i++) priors[i] = new UniformPDF(-1.0, 1.0);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//
// Build in Likelihood function
//
//////////////////////////////////////////////////////////////////////////////////////////////////

BaseLikelihood::BaseLikelihood(){}
BaseLikelihood::~BaseLikelihood(){}

GaussianLikelihood::GaussianLikelihood(int outputs, TypeLikelihood likelihood, const double covariance[], int data_entries, const double data[])
    : num_outputs(outputs), model_scale(-(double) data_entries), likely_type(likelihood), covariance_cache(0), data_cache(0)
{
    if (likely_type == likely_gauss_scale){
        covariance_cache = new double[1];
        covariance_cache[0] = 1.0 / covariance[0];
        model_scale *= covariance_cache[0];
    }else if (likely_type == likely_gauss_scale){
        covariance_cache = new double[num_outputs];
        for(int i=0; i<num_outputs; i++) covariance_cache[i] = 1.0 / covariance[i];
    }else if (likely_type == likely_gauss_dense){
        covariance_cache = new double[num_outputs * num_outputs];
        std::copy(covariance, covariance + num_outputs * num_outputs, covariance_cache);
        TasGrid::TasBLAS::cholUTU(num_outputs, covariance_cache); // lazy cholesky, add some error checking here
    } // when more likelihood types are added, add else error check. For now we get invalid type error so no need to check

    data_cache = new double[num_outputs];
    if (data_entries == 1){
        for(int i=0; i<num_outputs; i++) data_cache[i] = 2.0 * data[i];
    }else{
        double *scale_vector = new double[data_entries];
        std::fill(scale_vector, scale_vector + data_entries, 2.0);
        TasGrid::TasBLAS::dgemv(num_outputs, data_entries, data, scale_vector, data_cache);
        delete[] scale_vector;
    }
}
GaussianLikelihood::~GaussianLikelihood(){
    if (covariance_cache != 0){ delete[] covariance_cache; covariance_cache = 0; }
    if (data_cache != 0){ delete[] data_cache; data_cache = 0; }
}

void GaussianLikelihood::getLikelihood(int num_model, const double *model, std::vector<double> &likelihood, int, const double*, bool useLogForm){
    likelihood.resize(num_model);
    switch (likely_type){
    case likely_gauss_scale:{
        for(int i=0; i<num_model; i++){
            likelihood[i] = model_scale * TasGrid::TasBLAS::ddot(num_outputs, &(model[i*num_outputs])) + covariance_cache[0] * TasGrid::TasBLAS::ddot(num_outputs, data_cache, &(model[i*num_outputs]));
        }
        break;}
    case likely_gauss_diagonal:{
        double *scaled_model = new double[num_outputs];
        for(int i=0; i<num_model; i++){
            for(int j=0; j<num_outputs; j++) scaled_model[j] = covariance_cache[j] * model[i*num_outputs + j];
            likelihood[i] = model_scale * TasGrid::TasBLAS::ddot(num_outputs, scaled_model, &(model[i*num_outputs])) + TasGrid::TasBLAS::ddot(num_outputs, data_cache, scaled_model);
        }
        delete[] scaled_model;
        break;}
    case likely_gauss_dense:{
        #ifdef TASMANIAN_CPU_BLAS // leverage level 3 operations (uses more RAM)
        double *cov_by_model = new double[num_model * num_outputs];
        std::copy(model, model + num_model * num_outputs, cov_by_model);
        TasGrid::TasBLAS::dtrsm_LUTN(num_outputs, num_model, covariance_cache, cov_by_model);
        TasGrid::TasBLAS::dtrsm_LUNN(num_outputs, num_model, covariance_cache, cov_by_model);
        TasGrid::TasBLAS::dgemtv(num_outputs, num_model, cov_by_model, data_cache, likelihood.data());
        for(int i=0; i<num_model; i++) likelihood[i] += model_scale * TasGrid::TasBLAS::ddot(num_outputs, &(cov_by_model[i*num_outputs]), &(model[i*num_outputs]));
        delete[] cov_by_model;
        #else // conserve RAM, cannot use level 3 BLAS anyway
        double *cov_by_model = new double[num_outputs];
        for(int i=0; i<num_model; i++){
            std::copy(&(model[i*num_outputs]), &(model[i*num_outputs]) + num_outputs, cov_by_model);
            TasGrid::TasBLAS::dtrsv_UTN(num_outputs, covariance_cache, cov_by_model);
            TasGrid::TasBLAS::dtrsv_UNN(num_outputs, covariance_cache, cov_by_model);
            likelihood[i] = model_scale * TasGrid::TasBLAS::ddot(num_outputs, cov_by_model, &(model[i*num_outputs])) + TasGrid::TasBLAS::ddot(num_outputs, data_cache, &(model[i*num_outputs]));
        }
        delete[] cov_by_model;
        #endif // TASMANIAN_CPU_BLAS
        break;}
    }
    if (!useLogForm) for(auto &l : likelihood) l = exp(l);
}

}

#endif

