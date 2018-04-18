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

#ifndef __TASMANIAN_SPARSE_GRID_CPP
#define __TASMANIAN_SPARSE_GRID_CPP

#include "TasmanianSparseGrid.hpp"

#include "tsgCudaMacros.hpp"

#if defined(Tasmanian_ENABLE_CUBLAS) || defined(Tasmanian_ENABLE_CUDA)
#define _TASMANIAN_SETGPU cudaSetDevice(gpuID);
#endif // defined

namespace TasGrid{

const char* TasmanianSparseGrid::getVersion(){ return TASMANIAN_VERSION_STRING; }
const char* TasmanianSparseGrid::getLicense(){ return TASMANIAN_LICENSE; }
int TasmanianSparseGrid::getVersionMajor(){ return TASMANIAN_VERSION_MAJOR; }
int TasmanianSparseGrid::getVersionMinor(){ return TASMANIAN_VERSION_MINOR; }
bool TasmanianSparseGrid::isOpenMPEnabled(){
    #ifdef _OPENMP
    return true;
    #else
    return false;
    #endif // _OPENMP
}

TasmanianSparseGrid::TasmanianSparseGrid() : base(0), global(0), sequence(0), pwpoly(0), wavelet(0), domain_transform_a(0), domain_transform_b(0),
                                             conformal_asin_power(0), llimits(0), acceleration(accel_none), gpuID(0), acc_domain(0), logstream(0){
#ifndef TASMANIAN_XSDK
    logstream = &cerr;
#endif // TASMANIAN_XSDK
#ifdef Tasmanian_ENABLE_BLAS
    acceleration = accel_cpu_blas;
#endif // Tasmanian_ENABLE_BLAS
}
TasmanianSparseGrid::TasmanianSparseGrid(const TasmanianSparseGrid &source) : base(0), global(0), sequence(0), pwpoly(0), wavelet(0),
                                    domain_transform_a(0), domain_transform_b(0), conformal_asin_power(0), llimits(0),
                                    acceleration(accel_none), gpuID(0), acc_domain(0), logstream(0)
{
    copyGrid(&source);
#ifndef TASMANIAN_XSDK
    logstream = &cerr;
#endif // TASMANIAN_XSDK
#ifdef Tasmanian_ENABLE_BLAS
    acceleration = accel_cpu_blas;
#endif // Tasmanian_ENABLE_BLAS
}
TasmanianSparseGrid::~TasmanianSparseGrid(){
    clear();
}

void TasmanianSparseGrid::clear(){
    if (base != 0) base->clearAccelerationData();
    if (global != 0){ delete global; global = 0; }
    if (sequence != 0){ delete sequence; sequence = 0; }
    if (pwpoly != 0){ delete pwpoly; pwpoly = 0; }
    if (wavelet != 0){ delete wavelet; wavelet = 0; }
    if (domain_transform_a != 0){ delete[] domain_transform_a; domain_transform_a = 0; }
    if (domain_transform_b != 0){ delete[] domain_transform_b; domain_transform_b = 0; }
    if (conformal_asin_power != 0){ delete[] conformal_asin_power; conformal_asin_power = 0; }
    if (llimits != 0){ delete[] llimits; llimits = 0; }
    base = 0;
#ifndef TASMANIAN_XSDK
    logstream = &cerr;
#else
    logstream = 0;
#endif // TASMANIAN_XSDK
#ifdef Tasmanian_ENABLE_BLAS
    acceleration = accel_cpu_blas;
#else
    acceleration = accel_none;
#endif // Tasmanian_ENABLE_BLAS
#if defined (Tasmanian_ENABLE_CUBLAS) || defined (Tasmanian_ENABLE_CUDA)
    gpuID = 0;
#endif // Tasmanian_ENABLE_CUBLAS || Tasmanian_ENABLE_CUDA
    if (acc_domain != 0){ delete acc_domain; acc_domain = 0; }
}

void TasmanianSparseGrid::setErrorLog(std::ostream *os){ logstream = os; }
void TasmanianSparseGrid::disableLog(){ logstream = 0; }

void TasmanianSparseGrid::write(const char *filename, bool binary) const{
    std::ofstream ofs;
    if (binary){
        ofs.open(filename, std::ios::out | std::ios::binary);
    }else{
        ofs.open(filename);
    }
    write(ofs, binary);
    ofs.close();
}
bool TasmanianSparseGrid::read(const char *filename){
    std::ifstream ifs;
    char TSG[3];
    bool binary_format = false;
    ifs.open(filename, std::ios::in | std::ios::binary);
    ifs.read(TSG, 3 * sizeof(char));
    if ((TSG[0] == 'T') && (TSG[1] == 'S') && (TSG[2] == 'G')){
        binary_format = true;
    }
    ifs.close();
    if (binary_format){
        ifs.open(filename, std::ios::in | std::ios::binary);
    }else{
        ifs.open(filename);
    }
    bool isGood = read(ifs, binary_format);
    ifs.close();
    return isGood;
}

void TasmanianSparseGrid::write(std::ofstream &ofs, bool binary) const{
    if (binary){
        writeBinary(ofs);
    }else{
        writeAscii(ofs);
    }
}
bool TasmanianSparseGrid::read(std::ifstream &ifs, bool binary){
    if (binary){
        return readBinary(ifs);
    }else{
        return readAscii(ifs);
    }
}

void TasmanianSparseGrid::makeGlobalGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights, double alpha, double beta, const char* custom_filename, const int *level_limits){
    clear();
    global = new GridGlobal();
    global->makeGrid(dimensions, outputs, depth, type, rule, anisotropic_weights, alpha, beta, custom_filename, level_limits);
    base = global;
    if (level_limits != 0){
        llimits = new int[dimensions];
        std::copy(level_limits, level_limits + dimensions, llimits);
    }
}
void TasmanianSparseGrid::makeSequenceGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights, const int *level_limits){
    if (outputs < 1){
        if (logstream != 0){ (*logstream) << "ERROR: makeSequenceGrid is called with zero outputs, for zero outputs use makeGlobalGrid instead" << endl; }
        return;
    }
    if (OneDimensionalMeta::isSequence(rule)){
        clear();
        sequence = new GridSequence();
        sequence->makeGrid(dimensions, outputs, depth, type, rule, anisotropic_weights, level_limits);
        base = sequence;
        if (level_limits != 0){
            llimits = new int[dimensions];
            std::copy(level_limits, level_limits + dimensions, llimits);
        }
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: makeSequenceGrid is called with rule " << OneDimensionalMeta::getIORuleString(rule) << " which is not a sequence rule" << endl; }
    }
}
void TasmanianSparseGrid::makeLocalPolynomialGrid(int dimensions, int outputs, int depth, int order, TypeOneDRule rule, const int *level_limits){
    if ((rule != rule_localp) && (rule != rule_localp0) && (rule != rule_semilocalp)){
        if (logstream != 0){
            (*logstream) << "ERROR: makeLocalPolynomialGrid is called with rule " << OneDimensionalMeta::getIORuleString(rule) << " which is not a local polynomial rule" << endl;
            (*logstream) << "       use either " << OneDimensionalMeta::getIORuleString(rule_localp) << " or " << OneDimensionalMeta::getIORuleString(rule_semilocalp)
                         << " or " << OneDimensionalMeta::getIORuleString(rule_localp0)  << endl;
        }
        return;
    }
    if (order < -1){
        if (logstream != 0){ (*logstream) << "ERROR: makeLocalPolynomialGrid is called with order " << order << ", but the order cannot be less than -1." << endl; }
        return;
    }
    clear();
    pwpoly = new GridLocalPolynomial();
    pwpoly->makeGrid(dimensions, outputs, depth, order, rule, level_limits);
    base = pwpoly;
    if (level_limits != 0){
        llimits = new int[dimensions];
        std::copy(level_limits, level_limits + dimensions, llimits);
    }
}
void TasmanianSparseGrid::makeWaveletGrid(int dimensions, int outputs, int depth, int order, const int *level_limits){
    if ((order != 1) && (order != 3)){
        if (logstream != 0){ (*logstream) << "ERROR: makeWaveletGrid is called with order " << order << ", but wavelets are implemented only for orders 1 and 3." << endl; }
        return;
    }
    clear();
    wavelet = new GridWavelet();
    wavelet->makeGrid(dimensions, outputs, depth, order, level_limits);
    base = wavelet;
    if (level_limits != 0){
        llimits = new int[dimensions];
        std::copy(level_limits, level_limits + dimensions, llimits);
    }
}
void TasmanianSparseGrid::copyGrid(const TasmanianSparseGrid *source){
    clear();
    if (source->global != 0){
        global = new GridGlobal(*(source->global));
        base = global;
    }else if (source->sequence != 0){
        sequence = new GridSequence(*(source->sequence));
        base = sequence;
    }else if (source->pwpoly != 0){
        pwpoly = new GridLocalPolynomial(*(source->pwpoly));
        base = pwpoly;
    }else if (source->wavelet != 0){
        wavelet = new GridWavelet(*(source->wavelet));
        base = wavelet;
    }
    if (source->domain_transform_a != 0){
        setDomainTransform(source->domain_transform_a, source->domain_transform_b);
    }
    if (source->conformal_asin_power != 0){
        setConformalTransformASIN(source->conformal_asin_power);
    }
    if (source->llimits != 0){
        llimits = new int[base->getNumDimensions()];
        std::copy(source->llimits, source->llimits + base->getNumDimensions(), llimits);
    }
}

void TasmanianSparseGrid::updateGlobalGrid(int depth, TypeDepth type, const int *anisotropic_weights, const int *level_limits){
    if (global != 0){
        global->updateGrid(depth, type, anisotropic_weights, level_limits);
        if (level_limits != 0){
            if (llimits == 0) llimits = new int[global->getNumDimensions()];
            std::copy(level_limits, level_limits + global->getNumDimensions(), llimits);
        }
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: updateGlobalGrid called, but the grid is not global" << endl; }
    }
}
void TasmanianSparseGrid::updateSequenceGrid(int depth, TypeDepth type, const int *anisotropic_weights, const int *level_limits){
    if (sequence != 0){
        sequence->updateGrid(depth, type, anisotropic_weights, level_limits);
        if (level_limits != 0){
            if (llimits == 0) llimits = new int[sequence->getNumDimensions()];
            std::copy(level_limits, level_limits + sequence->getNumDimensions(), llimits);
        }
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: updateSequenceGrid called, but the grid is not sequence" << endl; }
    }
}

double TasmanianSparseGrid::getAlpha() const{
    return (global != 0) ? global->getAlpha() : 0.0;
}
double TasmanianSparseGrid::getBeta() const{
    return (global != 0) ? global->getBeta() : 0.0;
}
int TasmanianSparseGrid::getOrder() const{
    return (pwpoly != 0) ? pwpoly->getOrder() : ((wavelet != 0) ? wavelet->getOrder() : -1);
}

int TasmanianSparseGrid::getNumDimensions() const{ return (base == 0) ? 0 : base->getNumDimensions(); }
int TasmanianSparseGrid::getNumOutputs() const{ return (base == 0) ? 0 : base->getNumOutputs(); }
TypeOneDRule TasmanianSparseGrid::getRule() const{ return (base == 0) ? rule_none :base->getRule(); }
const char* TasmanianSparseGrid::getCustomRuleDescription() const{ return (global != 0) ? global->getCustomRuleDescription() : ""; }

int TasmanianSparseGrid::getNumLoaded() const{ return (base == 0) ? 0 : base->getNumLoaded(); }
int TasmanianSparseGrid::getNumNeeded() const{ return (base == 0) ? 0 : base->getNumNeeded(); }
int TasmanianSparseGrid::getNumPoints() const{ return (base == 0) ? 0 : base->getNumPoints(); }

double* TasmanianSparseGrid::getLoadedPoints() const{
    double *x = base->getLoadedPoints();
    formTransformedPoints(base->getNumLoaded(), x);
    return x;
}
void TasmanianSparseGrid::getLoadedPoints(double *x) const{
    base->getLoadedPoints(x);
    formTransformedPoints(base->getNumLoaded(), x);
}
double* TasmanianSparseGrid::getNeededPoints() const{
    double *x = base->getNeededPoints();
    formTransformedPoints(base->getNumNeeded(), x);
    return x;
}
void TasmanianSparseGrid::getNeededPoints(double *x) const{
    base->getNeededPoints(x);
    formTransformedPoints(base->getNumNeeded(), x);
}
double* TasmanianSparseGrid::getPoints() const{
    double *x = base->getPoints();
    formTransformedPoints(base->getNumPoints(), x);
    return x;
}
void TasmanianSparseGrid::getPoints(double *x) const{
    base->getPoints(x);
    formTransformedPoints(base->getNumPoints(), x);
}

double* TasmanianSparseGrid::getQuadratureWeights() const{
    double *w = base->getQuadratureWeights();
    mapConformalWeights(base->getNumDimensions(), base->getNumPoints(), w);
    if (domain_transform_a != 0){
        double scale = getQuadratureScale(base->getNumDimensions(), base->getRule());
        #pragma omp parallel for schedule(static)
        for(int i=0; i<getNumPoints(); i++) w[i] *= scale;
    }
    return w;
}
void TasmanianSparseGrid::getQuadratureWeights(double *weights) const{
    base->getQuadratureWeights(weights);
    mapConformalWeights(base->getNumDimensions(), base->getNumPoints(), weights);
    if (domain_transform_a != 0){
        double scale = getQuadratureScale(base->getNumDimensions(), base->getRule());
        #pragma omp parallel for schedule(static)
        for(int i=0; i<getNumPoints(); i++) weights[i] *= scale;
    }
}
double* TasmanianSparseGrid::getInterpolationWeights(const double x[]) const{
    double *x_tmp = 0;
    double *w = base->getInterpolationWeights(formCanonicalPoints(x, x_tmp, 1));
    clearCanonicalPoints(x_tmp);
    return w;
}
void TasmanianSparseGrid::getInterpolationWeights(const double x[], double *weights) const{
    double *x_tmp = 0;
    base->getInterpolationWeights(formCanonicalPoints(x, x_tmp, 1), weights);
    clearCanonicalPoints(x_tmp);
}

void TasmanianSparseGrid::loadNeededPoints(const double *vals){
    #if defined(Tasmanian_ENABLE_CUBLAS) || defined(Tasmanian_ENABLE_CUDA)
    if (AccelerationMeta::isAccTypeGPU(acceleration)){
        _TASMANIAN_SETGPU
    }
    #endif
    base->loadNeededPoints(vals, acceleration);
}

void TasmanianSparseGrid::evaluate(const double x[], double y[]) const{
    double *x_tmp = 0;
    base->evaluate(formCanonicalPoints(x, x_tmp, 1), y);
    clearCanonicalPoints(x_tmp);
}
void TasmanianSparseGrid::evaluateFast(const double x[], double y[]) const{
    double *x_tmp = 0;
    const double *x_canonical = formCanonicalPoints(x, x_tmp, 1);
    switch (acceleration){
        case accel_gpu_default:
        case accel_gpu_fullmemory:
        case accel_gpu_cublas:
            #ifdef Tasmanian_ENABLE_CUBLAS
            _TASMANIAN_SETGPU
            base->evaluateFastGPUcublas(x_canonical, y, logstream);
            break;
            #endif // Tasmanian_ENABLE_CUBLAS
        case accel_gpu_cuda:
            #ifdef Tasmanian_ENABLE_CUDA
            _TASMANIAN_SETGPU
            base->evaluateFastGPUcuda(x_canonical, y, logstream);
            break;
            #elif defined(Tasmanian_ENABLE_CUBLAS)
            _TASMANIAN_SETGPU
            base->evaluateFastGPUcublas(x_canonical, y, logstream);
            break;
            #endif // Tasmanian_ENABLE_CUDA
        case accel_cpu_blas:
            #ifdef Tasmanian_ENABLE_BLAS
            base->evaluateFastCPUblas(x_canonical, y);
            break;
            #endif // Tasmanian_ENABLE_BLAS
        default:
            base->evaluate(x_canonical, y);
            break;
    }
    clearCanonicalPoints(x_tmp);
}

void TasmanianSparseGrid::evaluateBatch(const double x[], int num_x, double y[]) const{
    double *x_tmp = 0;
    const double *x_canonical = formCanonicalPoints(x, x_tmp, num_x);
    //for(int i=0; i<4; i++){
    //    for(int j=0; j<base->getNumDimensions(); j++){
    //        cout << x[i*base->getNumDimensions() + j] << "  ";
    //    }
    //    cout << endl;
    //}
    switch (acceleration){
        case accel_gpu_default:
        case accel_gpu_fullmemory:
        case accel_gpu_cublas:
            #ifdef Tasmanian_ENABLE_CUBLAS
            _TASMANIAN_SETGPU
            base->evaluateBatchGPUcublas(x_canonical, num_x, y, logstream);
            break;
            #endif // Tasmanian_ENABLE_CUBLAS
        case accel_gpu_cuda:
            #ifdef Tasmanian_ENABLE_CUDA
            _TASMANIAN_SETGPU
            base->evaluateBatchGPUcuda(x_canonical, num_x, y, logstream);
            break;
            #elif defined(Tasmanian_ENABLE_CUBLAS)
            _TASMANIAN_SETGPU
            base->evaluateBatchGPUcublas(x_canonical, num_x, y, logstream);
            break;
            #endif // Tasmanian_ENABLE_CUDA
        case accel_cpu_blas:
            #ifdef Tasmanian_ENABLE_BLAS
            base->evaluateBatchCPUblas(x_canonical, num_x, y);
            break;
            #endif // Tasmanian_ENABLE_BLAS
        default:
            base->evaluateBatch(x_canonical, num_x, y);
            break;
    }
    clearCanonicalPoints(x_tmp);
}
void TasmanianSparseGrid::integrate(double q[]) const{
    if (conformal_asin_power != 0){
        int num_points = base->getNumPoints();
        double *correction = new double[num_points];  std::fill(correction, correction + num_points, 1.0);
        mapConformalWeights(base->getNumDimensions(), num_points, correction);
        base->integrate(q, correction);
        delete[] correction;
    }else{
        base->integrate(q, 0);
    }
    if (domain_transform_a != 0){
        double scale = getQuadratureScale(base->getNumDimensions(), base->getRule());
        for(int k=0; k<getNumOutputs(); k++) q[k] *= scale;
    }
}

bool TasmanianSparseGrid::isGlobal() const{
    return (global != 0);
}
bool TasmanianSparseGrid::isSequence() const{
    return (sequence != 0);
}
bool TasmanianSparseGrid::isLocalPolynomial() const{
    return (pwpoly != 0);
}
bool TasmanianSparseGrid::isWavelet() const{
    return (wavelet != 0);
}

void TasmanianSparseGrid::setDomainTransform(const double a[], const double b[]){
    if ((base == 0) || (base->getNumDimensions() == 0)){
        if (logstream != 0){ (*logstream) << "ERROR: cannot call setDomainTransform on uninitialized grid!" << endl; }
        return;
    }
    clearDomainTransform();
    if (acc_domain != 0){ delete acc_domain; acc_domain = 0; }
    int num_dimensions = base->getNumDimensions();
    domain_transform_a = new double[num_dimensions];  std::copy(a, a + num_dimensions, domain_transform_a);
    domain_transform_b = new double[num_dimensions];  std::copy(b, b + num_dimensions, domain_transform_b);
}
bool TasmanianSparseGrid::isSetDomainTransfrom() const{
    return (domain_transform_a != 0);
}
void TasmanianSparseGrid::clearDomainTransform(){
    if (domain_transform_a != 0){ delete[] domain_transform_a; domain_transform_a = 0; }
    if (domain_transform_b != 0){ delete[] domain_transform_b; domain_transform_b = 0; }
}
void TasmanianSparseGrid::getDomainTransform(double a[], double b[]) const{
    if ((base == 0) || (base->getNumDimensions() == 0) || (domain_transform_a == 0)){
        if (logstream != 0){ (*logstream) << "ERROR: cannot call getDomainTransform on uninitialized grid or if no transform has been set!" << endl; }
        return;
    }
    int num_dimensions = base->getNumDimensions();
    std::copy(domain_transform_a, domain_transform_a + num_dimensions, a);
    std::copy(domain_transform_b, domain_transform_b + num_dimensions, b);
}

void TasmanianSparseGrid::mapCanonicalToTransformed(int num_dimensions, int num_points, TypeOneDRule rule, double x[]) const{
    if ((rule == rule_gausslaguerre) || (rule == rule_gausslaguerreodd)){
        for(int i=0; i<num_points * num_dimensions; i++){
            int j = i % num_dimensions;
            x[i] /= domain_transform_b[j];
            x[i] += domain_transform_a[j];
        }
    }else if ((rule == rule_gausshermite) || (rule == rule_gausshermiteodd)){
        double *sqrt_b = new double[num_dimensions];
        for(int j=0; j<num_dimensions; j++) sqrt_b[j] = sqrt(domain_transform_b[j]);
        for(int i=0; i<num_points * num_dimensions; i++){
            int j = i % num_dimensions;
            x[i] /= sqrt_b[j];
            x[i] += domain_transform_a[j];
        }
        delete[] sqrt_b;
    }else{
        double *rate = new double[num_dimensions];
        double *shift = new double[num_dimensions];
        for(int j=0; j<num_dimensions; j++){
            rate[j]  = 0.5* (domain_transform_b[j] - domain_transform_a[j]);
            shift[j] = 0.5* (domain_transform_b[j] + domain_transform_a[j]);
        }
        for(int i=0; i<num_points * num_dimensions; i++){
            int j = i % num_dimensions;
            x[i] *= rate[j];
            x[i] += shift[j];
        }
        delete[] rate;
        delete[] shift;
    }
}
void TasmanianSparseGrid::mapTransformedToCanonical(int num_dimensions, int num_points, TypeOneDRule rule, double x[]) const{
    if ((rule == rule_gausslaguerre) || (rule == rule_gausslaguerreodd)){
        for(int i=0; i<num_points * num_dimensions; i++){
            int j = i % num_dimensions;
            x[i] -= domain_transform_a[j];
            x[i] *= domain_transform_b[j];
        }
    }else if ((rule == rule_gausshermite) || (rule == rule_gausshermiteodd)){
        double *sqrt_b = new double[num_dimensions];
        for(int j=0; j<num_dimensions; j++) sqrt_b[j] = sqrt(domain_transform_b[j]);
        for(int i=0; i<num_points * num_dimensions; i++){
            int j = i % num_dimensions;
            x[i] -= domain_transform_a[j];
            x[i] *= sqrt_b[j];
        }
        delete[] sqrt_b;
    }else{
        double *rate = new double[num_dimensions];
        double *shift = new double[num_dimensions];
        for(int j=0; j<num_dimensions; j++){
            rate[j]  = 2.0 / (domain_transform_b[j] - domain_transform_a[j]);
            shift[j] = (domain_transform_b[j] + domain_transform_a[j]) / (domain_transform_b[j] - domain_transform_a[j]);
        }
        for(int i=0; i<num_points * num_dimensions; i++){
            int j = i % num_dimensions;
            x[i] *= rate[j];
            x[i] -= shift[j];
        }
        delete[] rate;
        delete[] shift;
    }
}
double TasmanianSparseGrid::getQuadratureScale(int num_dimensions, TypeOneDRule rule) const{
    double scale = 1.0;
    if ((rule == rule_gausschebyshev1) || (rule == rule_gausschebyshev2) || (rule == rule_gaussgegenbauer) || (rule == rule_gaussjacobi)){
        double alpha = (rule == rule_gausschebyshev1) ? -0.5 : (rule == rule_gausschebyshev2) ? 0.5 : global->getAlpha();
        double beta = (rule == rule_gausschebyshev1) ? -0.5 : (rule == rule_gausschebyshev2) ? 0.5 : (rule == rule_gaussgegenbauer) ? global->getAlpha() : global->getBeta();
        for(int j=0; j<num_dimensions; j++) scale *= pow(0.5*(domain_transform_b[j] - domain_transform_a[j]), alpha + beta + 1.0);
    }else if ((rule == rule_gausslaguerre) || (rule == rule_gausslaguerreodd)){
        for(int j=0; j<num_dimensions; j++) scale *= pow(domain_transform_b[j], -(1.0 + global->getAlpha()));
    }else if ((rule == rule_gausshermite) || (rule == rule_gausshermiteodd)){
        double power = -0.5 * (1.0 + global->getAlpha());
        for(int j=0; j<num_dimensions; j++) scale *= pow(domain_transform_b[j], power);
    }else{
        for(int j=0; j<num_dimensions; j++) scale *= (domain_transform_b[j] - domain_transform_a[j]) / 2.0;
    }
    return scale;
}

void TasmanianSparseGrid::setConformalTransformASIN(const int truncation[]){
    if ((base == 0) || (base->getNumDimensions() == 0)){
        if (logstream != 0){ (*logstream) << "ERROR: cannot call setConformalTransformASIN on uninitialized grid!" << endl; }
        return;
    }
    clearConformalTransform();
    int num_dimensions = base->getNumDimensions();
    conformal_asin_power = new int[num_dimensions];  std::copy(truncation, truncation + num_dimensions, conformal_asin_power);
}
bool TasmanianSparseGrid::isSetConformalTransformASIN() const{ return (conformal_asin_power != 0); }
void TasmanianSparseGrid::clearConformalTransform(){
    if (conformal_asin_power != 0){ delete[] conformal_asin_power; conformal_asin_power = 0; }
}
void TasmanianSparseGrid::getConformalTransformASIN(int truncation[]) const{
    if ((base == 0) || (base->getNumDimensions() == 0) || (conformal_asin_power == 0)){
        if (logstream != 0){ (*logstream) << "ERROR: cannot call getDomainTransform on uninitialized grid or if no transform has been set!" << endl; }
        return;
    }
    int num_dimensions = base->getNumDimensions();
    std::copy(conformal_asin_power, conformal_asin_power + num_dimensions, truncation);
}

void TasmanianSparseGrid::mapConformalCanonicalToTransformed(int num_dimensions, int num_points, double x[]) const{
    if (conformal_asin_power != 0){
        // precompute constants, transform is sum exp(c_k + p_k * log(x))
        double **c = new double*[num_dimensions], **p = new double*[num_dimensions];
        int sum_powers = 0; for(int j=0; j<num_dimensions; j++) sum_powers += (conformal_asin_power[j] + 1);
        c[0] = new double[sum_powers]; p[0] = new double[sum_powers];
        sum_powers = 0;
        for(int j=1; j<num_dimensions; j++){
            sum_powers += (conformal_asin_power[j-1] + 1);
            c[j] = &(c[0][sum_powers]);
            p[j] = &(p[0][sum_powers]);
        }
        double lgamma_half = lgamma(0.5);
        double *cm = new double[num_dimensions];
        for(int j=0; j<num_dimensions; j++){
            double factorial = 0.0;
            cm[j] = 0.0;
            for(int k=0; k<=conformal_asin_power[j]; k++){
                p[j][k] = (double)(2*k+1);
                c[j][k] = lgamma(0.5 + ((double) k)) - lgamma_half - log(p[j][k]) - factorial;
                cm[j] += exp(c[j][k]);
                factorial += log((double)(k+1));
            }
        }
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++){
                if (x[i*num_dimensions+j] != 0.0){ // zero maps to zero and makes the log unstable
                    double sign = (x[i*num_dimensions+j] > 0.0) ? 1.0 : -1.0;
                    double logx = log(fabs(x[i*num_dimensions+j]));
                    x[i*num_dimensions+j] = 0.0;
                    for(int k=0; k<=conformal_asin_power[j]; k++){
                        x[i*num_dimensions+j] += exp(c[j][k] + p[j][k] * logx);
                    }
                    x[i*num_dimensions+j] *= sign / cm[j];
                }
            }
        }
        delete[] cm;
        delete[] c[0]; delete[] c;
        delete[] p[0]; delete[] p;
    }
}
void TasmanianSparseGrid::mapConformalTransformedToCanonical(int num_dimensions, int num_points, double x[]) const{
    if (conformal_asin_power != 0){
        // precompute constants, transform is sum exp(c_k + p_k * log(x))
        double **c = new double*[num_dimensions], **p = new double*[num_dimensions], **dc = new double*[num_dimensions], **dp = new double*[num_dimensions];
        int sum_powers = 0; for(int j=0; j<num_dimensions; j++) sum_powers += (conformal_asin_power[j] + 1);
        c[0] = new double[sum_powers]; p[0] = new double[sum_powers]; dc[0] = new double[sum_powers]; dp[0] = new double[sum_powers];
        sum_powers = 0;
        for(int j=1; j<num_dimensions; j++){
            sum_powers += (conformal_asin_power[j-1] + 1);
            c[j] = &(c[0][sum_powers]);
            p[j] = &(p[0][sum_powers]);
            dc[j] = &(dc[0][sum_powers]);
            dp[j] = &(dp[0][sum_powers]);
        }
        double lgamma_half = lgamma(0.5);
        double *cm = new double[num_dimensions];
        for(int j=0; j<num_dimensions; j++){
            double factorial = 0.0;
            cm[j] = 0.0;
            for(int k=0; k<=conformal_asin_power[j]; k++){
                p[j][k] = (double)(2*k+1);
                c[j][k] = lgamma(0.5 + ((double) k)) - lgamma_half - log(p[j][k]) - factorial;
                cm[j] += exp(c[j][k]);
                dp[j][k] = (double)(2*k);
                dc[j][k] = lgamma(0.5 + ((double) k)) - lgamma_half - factorial;
                factorial += log((double)(k+1));
            }
        }
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++){
                if (x[i*num_dimensions+j] != 0.0){ // zero maps to zero and makes the log unstable
                    double sign = (x[i*num_dimensions+j] > 0.0) ? 1.0 : -1.0;
                    x[i*num_dimensions+j] = fabs(x[i*num_dimensions+j]);
                    double b = x[i*num_dimensions+j];
                    double logx = log(x[i*num_dimensions+j]);
                    double r = x[i*num_dimensions+j];
                    double dr = 1.0;
                    for(int k=1; k<=conformal_asin_power[j]; k++){
                        r  += exp( c[j][k] +  p[j][k] * logx);
                        dr += exp(dc[j][k] + dp[j][k] * logx);
                   }
                    r /= cm[j];
                    r -= b; // transformed_x -b = 0
                    while(fabs(r) > TSG_NUM_TOL){
                        x[i*num_dimensions+j] -= r * cm[j] / dr;

                        logx = log(fabs(x[i*num_dimensions+j]));
                        r = x[i*num_dimensions+j]; dr = 1.0;
                        for(int k=1; k<=conformal_asin_power[j]; k++){
                            r  += exp( c[j][k] +  p[j][k] * logx);
                            dr += exp(dc[j][k] + dp[j][k] * logx);
                       }
                        r /= cm[j];
                        r -= b;
                   }
                    x[i*num_dimensions+j] *= sign;
                }
            }
        }
        delete[] cm;
        delete[] c[0]; delete[] c;
        delete[] p[0]; delete[] p;
        delete[] dc[0]; delete[] dc;
        delete[] dp[0]; delete[] dp;
    }
}
void TasmanianSparseGrid::mapConformalWeights(int num_dimensions, int num_points, double weights[]) const{
    if (conformal_asin_power != 0){
        // precompute constants, transform is sum exp(c_k + p_k * log(x))
        double *x = base->getPoints();
        double **c = new double*[num_dimensions], **p = new double*[num_dimensions];
        int sum_powers = 0; for(int j=0; j<num_dimensions; j++) sum_powers += (conformal_asin_power[j] + 1);
        c[0] = new double[sum_powers]; p[0] = new double[sum_powers];
        sum_powers = 0;
        for(int j=1; j<num_dimensions; j++){
            sum_powers += (conformal_asin_power[j-1] + 1);
            c[j] = &(c[0][sum_powers]);
            p[j] = &(p[0][sum_powers]);
        }
        double lgamma_half = lgamma(0.5);
        double *cm = new double[num_dimensions];
        for(int j=0; j<num_dimensions; j++){
            double factorial = 0.0;
            cm[j] = 0.0;
            for(int k=0; k<=conformal_asin_power[j]; k++){
                p[j][k] = (double)(2*k);
                c[j][k] = lgamma(0.5 + ((double) k)) - lgamma_half - factorial;
                factorial += log((double)(k+1));
                cm[j] += exp(c[j][k] - log((double)(2*k+1)));
            }
        }
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++){
                if (x[i*num_dimensions+j] != 0.0){ // derivative at zero is 1/cm[j] and zero makes the log unstable
                    double logx = log(fabs(x[i*num_dimensions+j]));
                    double trans = 1.0;
                    for(int k=1; k<=conformal_asin_power[j]; k++){
                        trans += exp(c[j][k] + p[j][k] * logx);
                   }
                    weights[i] *= trans / cm[j];
               }else{
                    weights[i] /= cm[j];
                }
            }
        }
        delete[] cm;
        delete[] c[0]; delete[] c;
        delete[] p[0]; delete[] p;
        delete[] x;
    }
}

const double* TasmanianSparseGrid::formCanonicalPoints(const double *x, double* &x_temp, int num_x) const{
    if ((domain_transform_a != 0) || (conformal_asin_power != 0)){
        int num_dimensions = base->getNumDimensions();
        x_temp = new double[num_dimensions*num_x]; std::copy(x, x + num_dimensions*num_x, x_temp);
        mapConformalTransformedToCanonical(num_dimensions, num_x, x_temp);
        if (domain_transform_a != 0) mapTransformedToCanonical(num_dimensions, num_x, base->getRule(), x_temp);
        return x_temp;
    }else{
        return x;
    }
}
void TasmanianSparseGrid::clearCanonicalPoints(double* &x_temp) const{
    if ((domain_transform_a != 0) || (conformal_asin_power != 0)){
        delete[] x_temp;
    }
}
void TasmanianSparseGrid::formTransformedPoints(int num_points, double x[]) const{
    mapConformalCanonicalToTransformed(base->getNumDimensions(), num_points, x); // internally switch based on the conformal transform
    if (domain_transform_a != 0){ // check the basic domain
        mapCanonicalToTransformed(base->getNumDimensions(), num_points, base->getRule(), x);
    }
}

#ifdef Tasmanian_ENABLE_CUDA
const double* TasmanianSparseGrid::formCanonicalPointsGPU(const double *gpu_x, double* &gpu_x_temp, int num_x) const{
    if (domain_transform_a != 0){
        if (acc_domain == 0) acc_domain = new AccelerationDomainTransform(base->getNumDimensions(), domain_transform_a, domain_transform_b, logstream);
        gpu_x_temp = acc_domain->getCanonicalPoints(base->getNumDimensions(), num_x, gpu_x);
        return gpu_x_temp;
    }else{
        return gpu_x;
    }
}
#else
const double* TasmanianSparseGrid::formCanonicalPointsGPU(const double *, double* &, int) const{ return 0; }
#endif // Tasmanian_ENABLE_CUDA

void TasmanianSparseGrid::clearLevelLimits(){
    if (llimits != 0){ delete[] llimits; llimits = 0; }
}
void TasmanianSparseGrid::getLevelLimits(int *limits) const{
    if (llimits == 0){
        if ((base != 0) && (base->getNumDimensions() > 0))
            std::fill(limits, limits + base->getNumDimensions(), -1);
    }else{
        std::copy(llimits, llimits + base->getNumDimensions(), limits);
    }
}

void TasmanianSparseGrid::setAnisotropicRefinement(TypeDepth type, int min_growth, int output, const int *level_limits){
    if (level_limits != 0){
        if (llimits == 0) llimits = new int[base->getNumDimensions()];
        std::copy(level_limits, level_limits + base->getNumDimensions(), llimits);
    }
    if (sequence != 0){
        sequence->setAnisotropicRefinement(type, min_growth, output, llimits);
    }else if (global != 0){
        if (OneDimensionalMeta::isNonNested(global->getRule())){
            if (logstream != 0){ (*logstream) << "ERROR: setAnisotropicRefinement called for a global grid with non-nested rule" << endl; }
        }else{
            global->setAnisotropicRefinement(type, min_growth, output, llimits);
        }
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: setAnisotropicRefinement called for grid that is neither sequence nor Global with sequence rule" << endl; }
    }
}
int* TasmanianSparseGrid::estimateAnisotropicCoefficients(TypeDepth type, int output){
    if (sequence != 0){
        return sequence->estimateAnisotropicCoefficients(type, output);
    }else if (global != 0){
        if (OneDimensionalMeta::isNonNested(global->getRule())){
            if (logstream != 0){ (*logstream) << "ERROR: estimateAnisotropicCoefficients called for a global grid with non-nested rule" << endl; }
        }else{
            return global->estimateAnisotropicCoefficients(type, output);
        }
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: estimateAnisotropicCoefficients called for grid that is neither sequence nor Global with sequence rule" << endl; }
    }
    return 0;
}
void TasmanianSparseGrid::setSurplusRefinement(double tolerance, int output, const int *level_limits){
    if (level_limits != 0){
        if (llimits == 0) llimits = new int[base->getNumDimensions()];
        std::copy(level_limits, level_limits + base->getNumDimensions(), llimits);
    }
    if (sequence != 0){
        sequence->setSurplusRefinement(tolerance, output, llimits);
    }else if (global != 0){
        if (OneDimensionalMeta::isSequence(global->getRule())){
            global->setSurplusRefinement(tolerance, output, llimits);
        }else{
            if (logstream != 0){ (*logstream) << "ERROR: setSurplusRefinement called for a global grid with non-sequence rule" << endl; }
        }
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: setSurplusRefinement(double, int) called for grid that is neither sequence nor Global with sequence rule" << endl; }
    }
}
void TasmanianSparseGrid::setSurplusRefinement(double tolerance, TypeRefinement criteria, int output, const int *level_limits, const double *scale_correction){
    if (level_limits != 0){
        if (llimits == 0) llimits = new int[base->getNumDimensions()];
        std::copy(level_limits, level_limits + base->getNumDimensions(), llimits);
    }
    if (pwpoly != 0){
        pwpoly->setSurplusRefinement(tolerance, criteria, output, llimits, scale_correction);
    }else if (wavelet != 0){
        wavelet->setSurplusRefinement(tolerance, criteria, output, llimits);
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: setSurplusRefinement(double, TypeRefinement) called for grid that is neither local polynomial nor wavelet" << endl; }
    }
}
void TasmanianSparseGrid::clearRefinement(){
    base->clearRefinement();
}
void TasmanianSparseGrid::mergeRefinement(){
    base->mergeRefinement();
}
void TasmanianSparseGrid::removePointsByHierarchicalCoefficient(double tolerance, int output, const double *scale_correction){
    if (pwpoly == 0){
        if (logstream != 0){ (*logstream) << "ERROR: removePointsBySurplus() called for a grid that is not local polynomial." << endl; }
        return;
    }else{
        if (pwpoly->removePointsByHierarchicalCoefficient(tolerance, output, scale_correction) == 0){
            clear();
        }
    }
}

void TasmanianSparseGrid::evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const{
    double *x_tmp = 0;
    base->evaluateHierarchicalFunctions(formCanonicalPoints(x, x_tmp, num_x), num_x, y);
    clearCanonicalPoints(x_tmp);
}
#ifdef Tasmanian_ENABLE_CUDA
void TasmanianSparseGrid::evaluateHierarchicalFunctionsGPU(const double gpu_x[], int cpu_num_x, double gpu_y[]) const{
    double *gpu_temp_x = 0;
    const double *gpu_canonical_x = formCanonicalPointsGPU(gpu_x, gpu_temp_x, cpu_num_x);
    pwpoly->buildDenseBasisMatrixGPU(gpu_canonical_x, cpu_num_x, gpu_y, logstream);
    if (gpu_temp_x != 0) TasCUDA::cudaDel<double>(gpu_temp_x, logstream);
}
void TasmanianSparseGrid::evaluateSparseHierarchicalFunctionsGPU(const double gpu_x[], int cpu_num_x, int* &gpu_pntr, int* &gpu_indx, double* &gpu_vals, int &num_nz) const{
    double *gpu_temp_x = 0;
    const double *gpu_canonical_x = formCanonicalPointsGPU(gpu_x, gpu_temp_x, cpu_num_x);
    pwpoly->buildSparseBasisMatrixGPU(gpu_canonical_x, cpu_num_x, gpu_pntr, gpu_indx, gpu_vals, num_nz, logstream);
    if (gpu_temp_x != 0) TasCUDA::cudaDel<double>(gpu_temp_x, logstream);
}
#else
void TasmanianSparseGrid::evaluateHierarchicalFunctionsGPU(const double*, int, double*) const{
    if (logstream != 0) (*logstream) << "ERROR: evaluateHierarchicalFunctionsGPU() called, but the library wasn't compiled with Tasmanian_ENABLE_CUDA=ON!" << endl;
}
void TasmanianSparseGrid::evaluateSparseHierarchicalFunctionsGPU(const double*, int, int*&, int*&, double*&, int&) const{
    if (logstream != 0) (*logstream) << "ERROR: evaluateSparseHierarchicalFunctionsGPU() called, but the library wasn't compiled with Tasmanian_ENABLE_CUDA=ON!" << endl;
}
#endif

void TasmanianSparseGrid::evaluateSparseHierarchicalFunctions(const double x[], int num_x, int* &pntr, int* &indx, double* &vals) const{
    double *x_tmp = 0;
    const double *x_canonical = formCanonicalPoints(x, x_tmp, num_x);
    if (pwpoly != 0){
        pwpoly->buildSpareBasisMatrix(x_canonical, num_x, 32, pntr, indx, vals);
    }else if (wavelet != 0){
        int num_points = base->getNumPoints();
        double *dense_vals = new double[num_points * num_x];
        wavelet->evaluateHierarchicalFunctions(x_canonical, num_x, dense_vals);
        int num_nz = 0;
        for(int i=0; i<num_points * num_x; i++) if (dense_vals[i] != 0.0) num_nz++;
        pntr = new int[num_x+1];
        indx = new int[num_nz];
        vals = new double[num_nz];
        num_nz = 0;
        for(int i=0; i<num_x; i++){
            pntr[i] = num_nz;
            for(int j=0; j<num_points; j++){
                if (dense_vals[i*num_points + j] != 0){
                    indx[num_nz] = j;
                    vals[num_nz++] = dense_vals[i*num_points + j];
                }
            }
        }
        pntr[num_x] = num_nz;
        delete[] dense_vals;
    }else{
        int num_points = base->getNumPoints();
        vals = new double[num_x * num_points];
        base->evaluateHierarchicalFunctions(x_canonical, num_x, vals);
        pntr = new int[num_x + 1];
        pntr[0] = 0;
        for(int i=0; i<num_x; i++) pntr[i+1] = pntr[i] + num_points;
        indx  = new int[num_x * num_points];
        for(int i=0; i<num_x; i++){
            for(int j=0; j<num_points; j++) indx[i*num_points + j] = j;
        }
    }
    clearCanonicalPoints(x_tmp);
}
int TasmanianSparseGrid::evaluateSparseHierarchicalFunctionsGetNZ(const double x[], int num_x) const{
    double *x_tmp = 0;
    const double *x_canonical = formCanonicalPoints(x, x_tmp, num_x);
    int num_nz = 0;
    if (pwpoly != 0){
        num_nz = pwpoly->getSpareBasisMatrixNZ(x_canonical, num_x, 32);
    }else if (wavelet != 0){
        int num_points = base->getNumPoints();
        double *dense_vals = new double[num_points * num_x];
        wavelet->evaluateHierarchicalFunctions(x_canonical, num_x, dense_vals);
        for(int i=0; i<num_points*num_x; i++) if (dense_vals[i] != 0.0) num_nz++;
        delete[] dense_vals;
    }else{
        return num_x * base->getNumPoints();
    }
    return num_nz;
    clearCanonicalPoints(x_tmp);
}
void TasmanianSparseGrid::evaluateSparseHierarchicalFunctionsStatic(const double x[], int num_x, int pntr[], int indx[], double vals[]) const{
    double *x_tmp = 0;
    const double *x_canonical = formCanonicalPoints(x, x_tmp, num_x);
    if (pwpoly != 0){
        pwpoly->buildSpareBasisMatrixStatic(x_canonical, num_x, 32, pntr, indx, vals);
    }else if (wavelet != 0){
        int num_points = base->getNumPoints();
        double *dense_vals = new double[num_points * num_x];
        wavelet->evaluateHierarchicalFunctions(x_canonical, num_x, dense_vals);
        int num_nz = 0;
        for(int i=0; i<num_points*num_x; i++) if (dense_vals[i] != 0.0) num_nz++;
        num_nz = 0;
        for(int i=0; i<num_x; i++){
            pntr[i] = num_nz;
            for(int j=0; j<num_points; j++){
                if (dense_vals[i*num_points + j] != 0){
                    indx[num_nz] = j;
                    vals[num_nz++] = dense_vals[i*num_points + j];
                }
            }
        }
        pntr[num_x] = num_nz;
        delete[] dense_vals;
    }else{
        int num_points = base->getNumPoints();
        base->evaluateHierarchicalFunctions(x_canonical, num_x, vals);
        pntr[0] = 0;
        for(int i=0; i<num_x; i++) pntr[i+1] = pntr[i] + num_points;
        for(int i=0; i<num_x; i++){
            for(int j=0; j<num_points; j++) indx[i*num_points + j] = j;
        }
    }
    clearCanonicalPoints(x_tmp);
}

void TasmanianSparseGrid::setHierarchicalCoefficients(const double c[]){
    base->setHierarchicalCoefficients(c, acceleration, logstream);
}

void TasmanianSparseGrid::getGlobalPolynomialSpace(bool interpolation, int &num_indexes, int* &poly) const{
    if (global != 0){
        global->getPolynomialSpace(interpolation, num_indexes, poly);
    }else if (sequence != 0){
        sequence->getPolynomialSpace(interpolation, num_indexes, poly);
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: getGlobalPolynomialSpace() called for a grid that is neither Global nor Sequence." << endl; }
        num_indexes = 0;
    }
}
const double* TasmanianSparseGrid::getHierarchicalCoefficients() const{
    if (pwpoly != 0){
        return pwpoly->getSurpluses();
    }else if (wavelet != 0){
        return wavelet->getSurpluses();
    }else if (sequence != 0){
        return sequence->getSurpluses();
    }else if (global != 0){
        return global->getLoadedValues();
    }else{
        return 0;
    }
}
const int* TasmanianSparseGrid::getPointsIndexes() const{
    if (pwpoly != 0){
        return pwpoly->getPointIndexes();
    }else if (wavelet != 0){
        return wavelet->getPointIndexes();
    }else if (global != 0){
        return global->getPointIndexes();
    }else if (sequence != 0){
        return sequence->getPointIndexes();
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: getPointIndexes() called for a grid that is neither local polynomial nor wavelet nor sequence." << endl; }
        return 0;
    }
}
const int* TasmanianSparseGrid::getNeededIndexes() const{
    if (pwpoly != 0){
        return pwpoly->getNeededIndexes();
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: getPointIndexes() called for a grid that is not local polynomial." << endl; }
        return 0;
    }
}

void TasmanianSparseGrid::printStats() const{
    printGridStats(&cout);
}
void TasmanianSparseGrid::printStatsLog() const{
    if (logstream == 0) return;
    printGridStats(logstream);
}

void TasmanianSparseGrid::printGridStats(std::ostream *os) const{
    if (os == 0) return;
    using std::setw;

    const int L1 = 20;
    (*os) << endl;
    (*os) << setw(L1) << "Grid Type:" << "  ";
    if (isGlobal()) (*os) << "Global";
    if (isSequence()) (*os) << "Sequence";
    if (isLocalPolynomial()) (*os) << "Local Polynomial";
    if (isWavelet()) (*os) << "Wavelets";
    if (!(isGlobal() || isSequence() || isLocalPolynomial() || isWavelet())) (*os) << "none";
    (*os) << endl;

    (*os) << setw(L1) << "Dimensions:" << "   " << getNumDimensions() << endl;
    (*os) << setw(L1) << "Outputs:" << "   " << getNumOutputs() << endl;
    if (getNumOutputs() == 0){
        (*os) << setw(L1) << "Nodes:" << "   " << getNumPoints() << endl;
    }else{
        (*os) << setw(L1) << "Loaded nodes:" << "   " << getNumLoaded() << endl;
        (*os) << setw(L1) << "Needed nodes:" << "   " << getNumNeeded() << endl;
    }
    (*os) << setw(L1) << "Rule:" << "  " << OneDimensionalMeta::getHumanString(getRule()) << endl;
    if (getRule() == rule_customtabulated){
        (*os) << setw(L1) << "Description:" << "  " << getCustomRuleDescription() << endl;
    }
    if (isSetDomainTransfrom()){
        (*os) << setw(L1) << "Domain:" << "  Custom" << endl;
    }else{
        (*os) << setw(L1) << "Domain:" << "  Canonical" << endl;
    }

    if (isGlobal()){
        TypeOneDRule rr = getRule();
        if ((rr == rule_gaussgegenbauer) || (rr == rule_gausslaguerre) || (rr == rule_gausshermite) || (rr == rule_gaussgegenbauerodd) || (rr == rule_gausshermiteodd) ){
            (*os) << setw(L1) << "Alpha:" << "   " << getAlpha() << endl;
        }
        if (rr == rule_gaussjacobi){
            (*os) << setw(L1) << "Alpha:" << "   " << getAlpha() << endl;
            (*os) << setw(L1) << "Beta:" << "   " << getBeta() << endl;
        }
    }else if (isSequence()){
        // sequence rules are simple, nothing to specify here
    }else if (isLocalPolynomial()){
        (*os) << setw(L1) << "Order:" << "   " << getOrder() << endl;
    }else if (isWavelet()){
        (*os) << setw(L1) << "Order:" << "   " << getOrder() << endl;
    }else{
        // empty grid, show nothing, just like the sequence grid
    }
    (*os) << setw(L1) << "Acceleration:" << "  " << AccelerationMeta::getIOAccelerationString(acceleration) << endl;
    if (AccelerationMeta::isAccTypeGPU(acceleration)){
        (*os) << setw(L1) << "GPU:" << "  " << getGPUID() << endl;
    }

    (*os) << endl;
}

void TasmanianSparseGrid::writeAscii(std::ofstream &ofs) const{
    ofs << "TASMANIAN SG " << getVersion() << endl;
    ofs << "WARNING: do not edit this manually" << endl;
    if (global != 0){
        ofs << "global" << endl;
        global->write(ofs);
    }else if (sequence != 0){
        ofs << "sequence" << endl;
        sequence->write(ofs);
    }else if (pwpoly != 0){
        ofs << "localpolynomial" << endl;
        pwpoly->write(ofs);
    }else if (wavelet != 0){
        ofs << "wavelet" << endl;
        wavelet->write(ofs);
    }else{
        ofs << "empty" << endl;
    }
    if (domain_transform_a != 0){
        ofs << "custom" << endl;
        ofs << std::scientific; ofs.precision(17);
        for(int j=0; j<base->getNumDimensions(); j++){
            ofs << domain_transform_a[j] << " " << domain_transform_b[j] << endl;
       }
    }else{
        ofs << "canonical" << endl;
    }
    if (conformal_asin_power != 0){
        ofs << "asinconformal" << endl;
        ofs << conformal_asin_power[0];
        for(int j=1; j<base->getNumDimensions(); j++){
            ofs << " " << conformal_asin_power[j];
       }
        ofs << endl;
    }else{
        ofs << "nonconformal" << endl;
    }
    if (llimits != 0){
        ofs << "limited" << endl;
        ofs << llimits[0];
        for(int j=1; j<base->getNumDimensions(); j++) ofs << " " << llimits[j];
        ofs << endl;
    }else{
        ofs << "unlimited" << endl;
    }
    ofs << "TASMANIAN SG end" << endl;
}
void TasmanianSparseGrid::writeBinary(std::ofstream &ofs) const{
    const char *TSG = "TSG5"; // last char indicates version
    ofs.write(TSG, 4 * sizeof(char)); // mark Tasmanian files
    char flag;
    // use Integers to indicate grid types, empty 'e', global 'g', sequence 's', pwpoly 'p', wavelet 'w'
    if (global != 0){
        flag = 'g'; ofs.write(&flag, sizeof(char));
        global->writeBinary(ofs);
    }else if (sequence != 0){
        flag = 's'; ofs.write(&flag, sizeof(char));
        sequence->writeBinary(ofs);
    }else if (pwpoly != 0){
        flag = 'p'; ofs.write(&flag, sizeof(char));
        pwpoly->writeBinary(ofs);
    }else if (wavelet != 0){
        flag = 'w'; ofs.write(&flag, sizeof(char));
        wavelet->writeBinary(ofs);
    }else{
        flag = 'e'; ofs.write(&flag, sizeof(char));
    }
    // domain transform: custom 'y', canonical: 'n'
    if (domain_transform_a != 0){
        flag = 'y'; ofs.write(&flag, sizeof(char));
        ofs.write((char*) domain_transform_a, base->getNumDimensions() * sizeof(double));
        ofs.write((char*) domain_transform_b, base->getNumDimensions() * sizeof(double));
    }else{
        flag = 'n'; ofs.write(&flag, sizeof(char));
    }
    // conformal transforms: none 'n', asin 'a'
    if (conformal_asin_power != 0){
        flag = 'a'; ofs.write(&flag, sizeof(char));
        ofs.write((char*) conformal_asin_power, base->getNumDimensions() * sizeof(int));
    }else{
        flag = 'n'; ofs.write(&flag, sizeof(char));
    }
    if (llimits != 0){
        flag = 'y'; ofs.write(&flag, sizeof(char));
        ofs.write((char*) llimits, base->getNumDimensions() * sizeof(int));
    }else{
        flag = 'n'; ofs.write(&flag, sizeof(char));
    }
    flag = 'e'; ofs.write(&flag, sizeof(char)); // E stands for END
}
bool TasmanianSparseGrid::readAscii(std::ifstream &ifs){
    std::string T;
    ifs >> T;  if (!(T.compare("TASMANIAN") == 0)){ if (logstream != 0){ (*logstream) << "ERROR: wrong file format, first word in not 'TASMANIAN'" << endl; } return false; }
    ifs >> T;  if (!(T.compare("SG") == 0)){ if (logstream != 0){ (*logstream) << "ERROR: wrong file format, second word in not 'SG'" << endl; } return false; }
    getline(ifs, T); T.erase(0,1); if (!(T.compare(getVersion()) == 0)){ if (logstream != 0){ (*logstream) << "WARNING: Version mismatch, possibly undefined behavior!" << endl; } }
    getline(ifs, T); if (!(T.compare("WARNING: do not edit this manually") == 0)){ if (logstream != 0){ (*logstream) << "ERROR: wrong file format, did not match 'WARNING: do not edit this manually'" << endl; } return false; }
    ifs >> T;
    if (T.compare("global") == 0){
        clear();
        global = new GridGlobal();
        global->read(ifs);
        base = global;
        getline(ifs, T);
    }else if (T.compare("sequence") == 0){
        clear();
        sequence = new GridSequence();
        sequence->read(ifs);
        base = sequence;
        getline(ifs, T);
    }else if (T.compare("localpolynomial") == 0){
        clear();
        pwpoly = new GridLocalPolynomial();
        pwpoly->read(ifs);
        base = pwpoly;
        getline(ifs, T);
    }else if (T.compare("wavelet") == 0){
        clear();
        wavelet = new GridWavelet();
        wavelet->read(ifs);
        base = wavelet;
        getline(ifs, T);
    }else if (T.compare("empty") == 0){
        clear();
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: wrong file format!" << endl; }
        return false;
    }
    getline(ifs, T);
    bool using_v3_format = false; // for compatibility with version 3.0
    bool using_v4_format = false; // for compatibility with version 3.1 and 4.0
    bool using_v5_format = false; // for compatibility with version 3.1 and 4.0
    if (T.compare("custom") == 0){
        domain_transform_a = new double[base->getNumDimensions()];
        domain_transform_b = new double[base->getNumDimensions()];
        for(int j=0; j<base->getNumDimensions(); j++){
            ifs >> domain_transform_a[j] >> domain_transform_b[j];
       }
        getline(ifs, T);
    }else if (T.compare("canonical") == 0){
        // do nothing, no domain transform
    }else if (T.compare("TASMANIAN SG end") == 0){
        // for compatibility with version 3.0 and the missing domain transform
        using_v3_format = true;
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: wrong file format! Domain is not specified!" << endl; }
        return false;
    }
    if (!using_v3_format){
        getline(ifs, T);
        if (T.compare("nonconformal") == 0){
            // do nothing, no conformal transform
        }else if (T.compare("asinconformal") == 0){
            conformal_asin_power = new int[base->getNumDimensions()];
            for(int j=0; j<base->getNumDimensions(); j++){
                ifs >> conformal_asin_power[j];
            }
            getline(ifs, T);
        }else if (T.compare("TASMANIAN SG end") == 0){
            // for compatibility with version 4.0/4.1 and the missing conformal maps
            using_v4_format = true;
        }else{
            if (logstream != 0){ (*logstream) << "ERROR: wrong file format! Conformal mapping is not specified!" << endl; }
            return false;
        }
        if (!using_v4_format){
            getline(ifs, T);
            if (T.compare("limited") == 0){
                llimits = new int[base->getNumDimensions()];
                for(int j=0; j<base->getNumDimensions(); j++) ifs >> llimits[j];
                getline(ifs, T);
            }else if (T.compare("unlimited") == 0){
                llimits = 0;
            }else if (T.compare("TASMANIAN SG end") == 0){
                using_v5_format = true;
            }else{
                if (logstream != 0){ (*logstream) << "WARNING: file did not end with 'TASMANIAN SG end', this may result in undefined behavior (v5)" << endl; }
            }
            if (!using_v5_format){
                getline(ifs, T);
                if (!(T.compare("TASMANIAN SG end") == 0)){
                    if (logstream != 0){ (*logstream) << "WARNING: file did not end with 'TASMANIAN SG end', this may result in undefined behavior (v51)" << endl; }
                    exit(1);
                }
            }
        }
    }

    return true;
}
bool TasmanianSparseGrid::readBinary(std::ifstream &ifs){
    char *TSG = new char[4];
    ifs.read(TSG, 4*sizeof(char));
    if ((TSG[0] != 'T') || (TSG[1] != 'S') || (TSG[2] != 'G')){
        if (logstream != 0){ (*logstream) << "ERROR: wrong file format, first 3 bytes are not 'TSG'" << endl; }
        return false;
    }
    if (TSG[3] != '5'){
        if (logstream != 0){ (*logstream) << "WARNING: grid seems to be saved in newer format, possibly undefined behavior" << endl; }
    }
    ifs.read(TSG, sizeof(char)); // what type of grid is it?
    if (TSG[0] == 'g'){
        clear();
        global = new GridGlobal();
        global->readBinary(ifs);
        base = global;
    }else if (TSG[0] == 's'){
        clear();
        sequence = new GridSequence();
        sequence->readBinary(ifs);
        base = sequence;
    }else if (TSG[0] == 'p'){
        clear();
        pwpoly = new GridLocalPolynomial();
        pwpoly->readBinary(ifs);
        base = pwpoly;
    }else if (TSG[0] == 'w'){
        clear();
        wavelet = new GridWavelet();
        wavelet->readBinary(ifs);
        base = wavelet;
    }else if (TSG[0] == 'e'){
        clear();
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: wrong grid type, wrong file format!" << endl; }
        return false;
    }
    ifs.read(TSG, sizeof(char)); // linear domain transform?
    if (TSG[0] == 'y'){
        domain_transform_a = new double[base->getNumDimensions()];
        domain_transform_b = new double[base->getNumDimensions()];
        ifs.read((char*) domain_transform_a, base->getNumDimensions() * sizeof(double));
        ifs.read((char*) domain_transform_b, base->getNumDimensions() * sizeof(double));
    }else if (TSG[0] != 'n'){
        if (logstream != 0){ (*logstream) << "ERROR: wrong domain type, wrong file format!" << endl; }
        return false;
    }
    ifs.read(TSG, sizeof(char)); // conformal domain transform?
    if (TSG[0] == 'a'){
        conformal_asin_power = new int[base->getNumDimensions()];
        ifs.read((char*) conformal_asin_power, base->getNumDimensions() * sizeof(int));
    }else if (TSG[0] != 'n'){
        if (logstream != 0){ (*logstream) << "ERROR: wrong conformal transform, wrong file format!" << endl; }
        return false;
    }
    ifs.read(TSG, sizeof(char)); // limits
    if (TSG[0] == 'y'){
        llimits = new int[base->getNumDimensions()];
        ifs.read((char*) llimits, base->getNumDimensions() * sizeof(int));
    }else if (TSG[0] != 'n'){
        if (logstream != 0){ (*logstream) << "ERROR: wrong level limits, wrong file format!" << endl; }
        return false;
    }
    ifs.read(TSG, sizeof(char)); // end character
    if (TSG[0] != 'e'){
        if (logstream != 0){ (*logstream) << "WARNING: the file did not terminate properly! Possibly undefined behavior." << endl; }
    }
    delete[] TSG;
    return true;
}

void TasmanianSparseGrid::enableAcceleration(TypeAcceleration acc){
    if (acc != acceleration){
        if (base != 0) base->clearAccelerationData();
        acceleration = acc;
        if (acc_domain != 0){ delete acc_domain; acc_domain = 0; }
    }
}
void TasmanianSparseGrid::forceSparseAlgorithmForLocalPolynomials(){
    if (pwpoly != 0) pwpoly->setForceSparse();
}
TypeAcceleration TasmanianSparseGrid::getAccelerationType() const{
    return acceleration;
}
bool TasmanianSparseGrid::isAccelerationAvailable(TypeAcceleration acc){
    switch (acc){
        case accel_none:   return true;
        #ifdef Tasmanian_ENABLE_BLAS
        case accel_cpu_blas:   return true;
        #else
        case accel_cpu_blas:   return false;
        #endif // Tasmanian_ENABLE_BLAS

        #ifdef Tasmanian_ENABLE_CUBLAS
        case accel_gpu_cublas:   return true;
        #else
        case accel_gpu_cublas:   return false;
        #endif // Tasmanian_ENABLE_CUBLAS

        #ifdef Tasmanian_ENABLE_CUDA
        case accel_gpu_cuda:   return true;
        #else
        case accel_gpu_cuda:   return false;
        #endif // Tasmanian_ENABLE_CUDA

        #ifdef TASMANIAN_ENABLE_MAGMA
        case accel_gpu_magma:   return true;
        #else
        case accel_gpu_magma:   return false;
        #endif // TASMANIAN_ENABLE_MAGMA

        #if defined(Tasmanian_ENABLE_CUDA) || defined(Tasmanian_ENABLE_CUBLAS) || defined(TASMANIAN_ENABLE_MAGMA)
        case accel_gpu_default:   return true;
        case accel_gpu_fullmemory:   return true;
        #else
        case accel_gpu_default:   return false;
        case accel_gpu_fullmemory:   return false;
        #endif // TASMANIAN_ENABLE_GPU
        default: return false;
    }
}

void TasmanianSparseGrid::setGPUID(int new_gpuID){
    if (new_gpuID != gpuID){
        if (AccelerationMeta::isAccTypeGPU(acceleration)){ // if using GPU acceleration
            if (base != 0) base->clearAccelerationData();
        }
        gpuID = new_gpuID;
    }
}
int TasmanianSparseGrid::getGPUID() const{ return gpuID; }

int TasmanianSparseGrid::getNumGPUs(){
    #if defined(Tasmanian_ENABLE_CUDA) || defined(Tasmanian_ENABLE_CUBLAS)
    int gpu_count = 0;
    cudaGetDeviceCount(&gpu_count);
    return gpu_count;
    #else
    return 0;
    #endif // Tasmanian_ENABLE_CUDA || Tasmanian_ENABLE_CUBLAS
}

#if defined(Tasmanian_ENABLE_CUDA) || defined(Tasmanian_ENABLE_CUBLAS)
int TasmanianSparseGrid::getGPUMemory(int gpu){
    if (gpu < 0) return 0;
    int gpu_count = 0;
    cudaGetDeviceCount(&gpu_count);
    if (gpu >= gpu_count) return 0;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, gpu);
    unsigned long int memB = prop.totalGlobalMem;
    return (int) (memB / 1048576);
}
char* TasmanianSparseGrid::getGPUName(int gpu){
    char *name = new char[1];
    name[0] = '\0';
    if (gpu < 0) return name;
    int gpu_count = 0;
    cudaGetDeviceCount(&gpu_count);
    if (gpu >= gpu_count) return name;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, gpu);

    int c = 0; while(prop.name[c] != '\0'){ c++; }
    delete[] name;
    name = new char[c+1];
    for(int i=0; i<c; i++){ name[i] = prop.name[i]; }
    name[c] = '\0';

    return name;
}
#else
int TasmanianSparseGrid::getGPUMemory(int){ return 0; }
char* TasmanianSparseGrid::getGPUName(int){
    char *name = new char[1];
    name[0] = '\0';
    return name;
}
#endif // Tasmanian_ENABLE_CUDA || Tasmanian_ENABLE_CUBLAS

// ------------ C Interface for use with Python ctypes and potentially other C codes -------------- //
extern "C" {
void* tsgConstructTasmanianSparseGrid(){ return (void*) new TasmanianSparseGrid(); }
void tsgDestructTasmanianSparseGrid(void *grid){ delete ((TasmanianSparseGrid*) grid); }

void tsgCopyGrid(void *destination, void *source){ ((TasmanianSparseGrid*) destination)->copyGrid(((TasmanianSparseGrid*) source)); }

const char* tsgGetVersion(){ return TasmanianSparseGrid::getVersion(); }
const char* tsgGetLicense(){ return TasmanianSparseGrid::getLicense(); }
int tsgGetVersionMajor(){ return TasmanianSparseGrid::getVersionMajor(); }
int tsgGetVersionMinor(){ return TasmanianSparseGrid::getVersionMinor(); }
int tsgIsOpenMPEnabled(){ return (TasmanianSparseGrid::isOpenMPEnabled()) ? 1 : 0; }

void tsgErrorLogCerr(void *grid){ ((TasmanianSparseGrid*) grid)->setErrorLog(&cerr); }
void tsgDisableErrorLog(void *grid){ ((TasmanianSparseGrid*) grid)->disableLog(); }

void tsgWrite(void *grid, const char* filename){ ((TasmanianSparseGrid*) grid)->write(filename); }
void tsgWriteBinary(void *grid, const char* filename){ ((TasmanianSparseGrid*) grid)->write(filename, true); }
int tsgRead(void *grid, const char* filename){
    bool result = ((TasmanianSparseGrid*) grid)->read(filename);
    return result ? 0 : 1;
}

void tsgMakeGlobalGrid(void *grid, int dimensions, int outputs, int depth, const char * sType, const char *sRule, const int *anisotropic_weights, double alpha, double beta, const char* custom_filename, const int *limit_levels){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    TypeOneDRule rule = OneDimensionalMeta::getIORuleString(sRule);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    if (rule == rule_none){ cerr << "WARNING: incorrect rule type: " << sType << ", defaulting to clenshaw-curtis." << endl; }
    #endif // _TASMANIAN_DEBUG_
    ((TasmanianSparseGrid*) grid)->makeGlobalGrid(dimensions, outputs, depth, depth_type, rule, anisotropic_weights, alpha, beta, custom_filename, limit_levels);
}
void tsgMakeSequenceGrid(void *grid, int dimensions, int outputs, int depth, const char *sType, const char *sRule, const int *anisotropic_weights, const int *limit_levels){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    TypeOneDRule rule = OneDimensionalMeta::getIORuleString(sRule);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    if (rule == rule_none){ cerr << "WARNING: incorrect rule type: " << sRule << ", defaulting to clenshaw-curtis." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    if (rule == rule_none){ rule = rule_clenshawcurtis; }
    ((TasmanianSparseGrid*) grid)->makeSequenceGrid(dimensions, outputs, depth, depth_type, rule, anisotropic_weights, limit_levels);
}
void tsgMakeLocalPolynomialGrid(void *grid, int dimensions, int outputs, int depth, int order, const char *sRule, const int *limit_levels){
    TypeOneDRule rule = OneDimensionalMeta::getIORuleString(sRule);
    #ifdef _TASMANIAN_DEBUG_
    if (rule == rule_none){ cerr << "WARNING: incorrect rule type: " << sRule << ", defaulting to localp." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (rule == rule_none){ rule = rule_localp; }
    ((TasmanianSparseGrid*) grid)->makeLocalPolynomialGrid(dimensions, outputs, depth, order, rule, limit_levels);
}
void tsgMakeWaveletGrid(void *grid, int dimensions, int outputs, int depth, int order, const int *limit_levels){
    ((TasmanianSparseGrid*) grid)->makeWaveletGrid(dimensions, outputs, depth, order, limit_levels);
}

void tsgUpdateGlobalGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->updateGlobalGrid(depth, depth_type, anisotropic_weights, limit_levels);
}
void tsgUpdateSequenceGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->updateSequenceGrid(depth, depth_type, anisotropic_weights, limit_levels);
}

double tsgGetAlpha(void *grid){ return ((TasmanianSparseGrid*) grid)->getAlpha(); }
double tsgGetBeta(void *grid){ return ((TasmanianSparseGrid*) grid)->getBeta(); }
int tsgGetOrder(void *grid){ return ((TasmanianSparseGrid*) grid)->getOrder(); }
int tsgGetNumDimensions(void *grid){ return ((TasmanianSparseGrid*) grid)->getNumDimensions(); }
int tsgGetNumOutputs(void *grid){ return ((TasmanianSparseGrid*) grid)->getNumOutputs(); }
const char* tsgGetRule(void *grid){ return  OneDimensionalMeta::getIORuleString( ((TasmanianSparseGrid*) grid)->getRule() ); }
const char* tsgGetCustomRuleDescription(void *grid){ return ((TasmanianSparseGrid*) grid)->getCustomRuleDescription(); }

int tsgGetNumLoaded(void *grid){ return ((TasmanianSparseGrid*) grid)->getNumLoaded(); }
int tsgGetNumNeeded(void *grid){ return ((TasmanianSparseGrid*) grid)->getNumNeeded(); }
int tsgGetNumPoints(void *grid){ return ((TasmanianSparseGrid*) grid)->getNumPoints(); }

void tsgGetLoadedPointsStatic(void *grid, double *x){ return ((TasmanianSparseGrid*) grid)->getLoadedPoints(x); }
double* tsgGetLoadedPoints(void *grid){
    if (((TasmanianSparseGrid*) grid)->getNumLoaded() == 0){
        return 0;
    }
    double *x = (double*) malloc(((TasmanianSparseGrid*) grid)->getNumLoaded() * ((TasmanianSparseGrid*) grid)->getNumDimensions() * sizeof(double));
    ((TasmanianSparseGrid*) grid)->getLoadedPoints(x);
    return x;
}
void tsgGetNeededPointsStatic(void *grid, double *x){ return ((TasmanianSparseGrid*) grid)->getNeededPoints(x); }
double* tsgGetNeededPoints(void *grid){
    if (((TasmanianSparseGrid*) grid)->getNumNeeded() == 0){
        return 0;
    }
    double *x = (double*) malloc(((TasmanianSparseGrid*) grid)->getNumNeeded() * ((TasmanianSparseGrid*) grid)->getNumDimensions() * sizeof(double));
    ((TasmanianSparseGrid*) grid)->getNeededPoints(x);
    return x;
}
void tsgGetPointsStatic(void *grid, double *x){ return ((TasmanianSparseGrid*) grid)->getPoints(x); }
double* tsgGetPoints(void *grid){
    if (((TasmanianSparseGrid*) grid)->getNumPoints() == 0){
        return 0;
    }
    double *x = (double*) malloc(((TasmanianSparseGrid*) grid)->getNumPoints() * ((TasmanianSparseGrid*) grid)->getNumDimensions() * sizeof(double));
    ((TasmanianSparseGrid*) grid)->getPoints(x);
    return x;
}

void tsgGetQuadratureWeightsStatic(void *grid, double *weights){ return ((TasmanianSparseGrid*) grid)->getQuadratureWeights(weights); }
double* tsgGetQuadratureWeights(void *grid){
    double *w = (double*) malloc(((TasmanianSparseGrid*) grid)->getNumPoints() * sizeof(double));
    ((TasmanianSparseGrid*) grid)->getQuadratureWeights(w);
    return w;
}
void tsgGetInterpolationWeightsStatic(void *grid, const double *x, double *weights){ ((TasmanianSparseGrid*) grid)->getInterpolationWeights(x, weights); }
double* tsgGetInterpolationWeights(void *grid, const double *x){
    double *w = (double*) malloc(((TasmanianSparseGrid*) grid)->getNumPoints() * sizeof(double));
    ((TasmanianSparseGrid*) grid)->getInterpolationWeights(x, w);
    return w;
}

void tsgLoadNeededPoints(void *grid, const double *vals){ ((TasmanianSparseGrid*) grid)->loadNeededPoints(vals); }

void tsgEvaluate(void *grid, const double *x, double *y){ ((TasmanianSparseGrid*) grid)->evaluate(x, y); }
void tsgEvaluateFast(void *grid, const double *x, double *y){ ((TasmanianSparseGrid*) grid)->evaluateFast(x, y); }
void tsgIntegrate(void *grid, double *q){ ((TasmanianSparseGrid*) grid)->integrate(q); }

void tsgEvaluateBatch(void *grid, const double *x, int num_x, double *y){ ((TasmanianSparseGrid*) grid)->evaluateBatch(x, num_x, y); }

void tsgBatchGetInterpolationWeightsStatic(void *grid, const double *x, int num_x, double *weights){
    TasmanianSparseGrid* tsg = (TasmanianSparseGrid*) grid;
    int iNumDim = tsg->getNumDimensions(), iNumPoints = tsg->getNumPoints();
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        tsg->getInterpolationWeights(&(x[i*iNumDim]), &(weights[i*iNumPoints]));
    }
}
double* tsgBatchGetInterpolationWeights(void *grid, const double *x, int num_x){
    double *weights = (double*) malloc(num_x * ((TasmanianSparseGrid*) grid)->getNumPoints() * sizeof(double));
    tsgBatchGetInterpolationWeightsStatic(grid, x, num_x, weights);
    return weights;
}

int tsgIsGlobal(void *grid){ return (((TasmanianSparseGrid*) grid)->isGlobal() ? 1 : 0); }
int tsgIsSequence(void *grid){ return (((TasmanianSparseGrid*) grid)->isSequence() ? 1 : 0); }
int tsgIsLocalPolynomial(void *grid){ return (((TasmanianSparseGrid*) grid)->isLocalPolynomial() ? 1 : 0); }
int tsgIsWavelet(void *grid){ return (((TasmanianSparseGrid*) grid)->isWavelet() ? 1 : 0); }

void tsgSetDomainTransform(void *grid, const double a[], const double b[]){ ((TasmanianSparseGrid*) grid)->setDomainTransform(a, b); }
int tsgIsSetDomainTransfrom(void *grid){ return (((TasmanianSparseGrid*) grid)->isSetDomainTransfrom() ? 1 : 0); }
void tsgClearDomainTransform(void *grid){ ((TasmanianSparseGrid*) grid)->clearDomainTransform(); }
void tsgGetDomainTransform(void *grid, double a[], double b[]){ ((TasmanianSparseGrid*) grid)->getDomainTransform(a, b); }

void tsgSetConformalTransformASIN(void *grid, const int truncation[]){ ((TasmanianSparseGrid*) grid)->setConformalTransformASIN(truncation); }
int tsgIsSetConformalTransformASIN(void *grid){ return (((TasmanianSparseGrid*) grid)->isSetConformalTransformASIN()) ? 1 : 0; }
void tsgClearConformalTransform(void *grid){ ((TasmanianSparseGrid*) grid)->clearConformalTransform(); }
void tsgGetConformalTransformASIN(void *grid, int truncation[]){ ((TasmanianSparseGrid*) grid)->getConformalTransformASIN(truncation); }

void tsgClearLevelLimits(void *grid){ ((TasmanianSparseGrid*) grid)->clearLevelLimits(); }
void tsgGetLevelLimits(void *grid, int *limits){ ((TasmanianSparseGrid*) grid)->getLevelLimits(limits); }

void tsgSetAnisotropicRefinement(void *grid, const char * sType, int min_growth, int output, const int *level_limits){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->setAnisotropicRefinement(depth_type, min_growth, output, level_limits);
}
int* tsgEstimateAnisotropicCoefficients(void *grid, const char * sType, int output, int *num_coefficients){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    *num_coefficients = ((TasmanianSparseGrid*) grid)->getNumDimensions();
    if ((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)){
        *num_coefficients *= 2;
    }
    int *coeff = ((TasmanianSparseGrid*) grid)->estimateAnisotropicCoefficients(depth_type, output);
    int *result = (int*) malloc((*num_coefficients) * sizeof(int));
    for(int i=0; i<*num_coefficients; i++) result[i] = coeff[i];
    delete[] coeff;
    return result;
}
void tsgEstimateAnisotropicCoefficientsStatic(void *grid, const char * sType, int output, int *coefficients){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    int num_coefficients = ((TasmanianSparseGrid*) grid)->getNumDimensions();
    if ((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)){
        num_coefficients *= 2;
    }
    int *coeff = ((TasmanianSparseGrid*) grid)->estimateAnisotropicCoefficients(depth_type, output);
    for(int i=0; i<num_coefficients; i++) coefficients[i] = coeff[i];
    delete[] coeff;
}
void tsgSetGlobalSurplusRefinement(void *grid, double tolerance, int output, const int *level_limits){
    ((TasmanianSparseGrid*) grid)->setSurplusRefinement(tolerance, output, level_limits);
}
void tsgSetLocalSurplusRefinement(void *grid, double tolerance, const char * sRefinementType, int output, const int *level_limits){
    TypeRefinement ref_type = OneDimensionalMeta::getIOTypeRefinementString(sRefinementType);
    #ifdef _TASMANIAN_DEBUG_
    if (ref_type == refine_none){ cerr << "WARNING: incorrect refinement type: " << sRefinementType << ", defaulting to type_classic." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (ref_type == refine_none){ ref_type = refine_classic; }
    ((TasmanianSparseGrid*) grid)->setSurplusRefinement(tolerance, ref_type, output, level_limits);
}
void tsgClearRefinement(void *grid){
    ((TasmanianSparseGrid*) grid)->clearRefinement();
}
void tsgMergeRefinement(void *grid){
    ((TasmanianSparseGrid*) grid)->mergeRefinement();
}
void tsgRemovePointsByHierarchicalCoefficient(void *grid, double tolerance, int output, const double *scale_correction){
    ((TasmanianSparseGrid*) grid)->removePointsByHierarchicalCoefficient(tolerance, output, scale_correction);
}

void tsgEvaluateHierarchicalFunctions(void *grid, const double *x, int num_x, double *y){
    ((TasmanianSparseGrid*) grid)->evaluateHierarchicalFunctions(x, num_x, y);
}
void tsgEvaluateSparseHierarchicalFunctions(void *grid, const double x[], int num_x, int **pntr, int **indx, double **vals){
    int num_nz = ((TasmanianSparseGrid*) grid)->evaluateSparseHierarchicalFunctionsGetNZ(x, num_x);
    *pntr = (int*) malloc((num_x+1) * sizeof(int));
    *indx = (int*) malloc(num_nz * sizeof(int));
    *vals = (double*) malloc(num_nz * sizeof(double));
    ((TasmanianSparseGrid*) grid)->evaluateSparseHierarchicalFunctionsStatic(x, num_x, *pntr, *indx, *vals);
}
int tsgEvaluateSparseHierarchicalFunctionsGetNZ(void *grid, const double x[], int num_x){
    return ((TasmanianSparseGrid*) grid)->evaluateSparseHierarchicalFunctionsGetNZ(x, num_x);
}
void tsgEvaluateSparseHierarchicalFunctionsStatic(void *grid, const double x[], int num_x, int *pntr, int *indx, double *vals){
    ((TasmanianSparseGrid*) grid)->evaluateSparseHierarchicalFunctionsStatic(x, num_x, pntr, indx, vals);
}
const double* tsgGetHierarchicalCoefficients(void *grid){
    return ((TasmanianSparseGrid*) grid)->getHierarchicalCoefficients();
}
void tsgGetHierarchicalCoefficientsStatic(void *grid, double *coeff){
    int num_points = ((TasmanianSparseGrid*) grid)->getNumPoints();
    int num_outputs = ((TasmanianSparseGrid*) grid)->getNumOutputs();
    if ((num_points == 0) || (num_outputs == 0)) return;
    const double *surp = ((TasmanianSparseGrid*) grid)->getHierarchicalCoefficients();
    std::copy(surp, surp + num_outputs * num_points, coeff);
}
void tsgSetHierarchicalCoefficients(void *grid, const double *c){
    ((TasmanianSparseGrid*) grid)->setHierarchicalCoefficients(c);
}

// to be called from Python only, must later call delete[] on the pointer
int* tsgPythonGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes){
    int *indx = 0;;
    ((TasmanianSparseGrid*) grid)->getGlobalPolynomialSpace((interpolation != 0), *num_indexes, indx);
    return indx;
}
// to be used in C, creates a C pointer (requires internal copy of data)
void tsgGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes, int **indexes){
    int *indx = 0, num_ind, num_dims = ((TasmanianSparseGrid*) grid)->getNumDimensions();
    ((TasmanianSparseGrid*) grid)->getGlobalPolynomialSpace((interpolation != 0), num_ind, indx);
    *num_indexes = num_ind;
    if (indx != 0){
        *indexes = (int*) malloc(((*num_indexes) * num_dims) * sizeof(int));
        for(int i=0; i<((*num_indexes) * num_dims); i++) *indexes[i] = indx[i];
        delete[] indx;
    }
}

void tsgPrintStats(void *grid){ ((TasmanianSparseGrid*) grid)->printStats(); }

void tsgEnableAcceleration(void *grid, const char *accel){ ((TasmanianSparseGrid*) grid)->enableAcceleration(AccelerationMeta::getIOAccelerationString(accel)); }
//int tsgGetAccelerationTypeInt(void *grid){ return AccelerationMeta::getIOAccelerationInt(((TasmanianSparseGrid*) grid)->getAccelerationType()); } // int to acceleration type
const char* tsgGetAccelerationType(void *grid){ return AccelerationMeta::getIOAccelerationString(((TasmanianSparseGrid*) grid)->getAccelerationType()); }

void tsgSetGPUID(void *grid, int gpuID){ ((TasmanianSparseGrid*) grid)->setGPUID(gpuID); }
int tsgGetGPUID(void *grid){ return ((TasmanianSparseGrid*) grid)->getGPUID(); }
int tsgGetNumGPUs(){ return TasmanianSparseGrid::getNumGPUs(); }
int tsgGetGPUMemory(int gpu){ return TasmanianSparseGrid::getGPUMemory(gpu); }
int tsgIsAccelerationAvailable(const char *accel){ return (TasmanianSparseGrid::isAccelerationAvailable(AccelerationMeta::getIOAccelerationString(accel))) ? 1 : 0; }
void tsgGetGPUName(int gpu, int num_buffer, char *buffer, int *num_actual){
    // gpu is the gpuID, num_buffer is the size of *buffer, num_actual returns the actual number of chars
    char* name = TasmanianSparseGrid::getGPUName(gpu);
    int c = 0;
    while ((name[c] != '\0') && (c < num_buffer-1)){
        buffer[c] = name[c];
        c++;
    }
    buffer[c] = '\0';
    *num_actual = c;
}

void tsgDeleteInts(int *p){ delete[] p; }

}
}

#endif
