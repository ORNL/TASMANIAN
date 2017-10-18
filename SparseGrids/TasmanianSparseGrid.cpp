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

#if defined(TASMANIAN_CUBLAS) || defined(TASMANIAN_CUDA)
#include <cuda_runtime_api.h>
#include <cuda.h>
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
bool TasmanianSparseGrid::isCudaEnabled(){
    #ifdef TASMANIAN_CUBLAS
    return true;
    #else
    return false;
    #endif // TASMANIAN_CUBLAS
}
bool TasmanianSparseGrid::isBLASEnabled(){
    #ifdef TASMANIAN_CPU_BLAS
    return true;
    #else
    return false;
    #endif // TASMANIAN_CPU_BLAS
}

TasmanianSparseGrid::TasmanianSparseGrid() : base(0), global(0), sequence(0), pwpoly(0), wavelet(0), domain_transform_a(0), domain_transform_b(0),
                                             conformal_asin_power(0), acceleration(accel_none), gpuID(0), logstream(0){
#ifndef TASMANIAN_XSDK
    logstream = &cerr;
#endif
#ifdef TASMANIAN_CPU_BLAS
    acceleration = accel_cpu_blas;
#endif // TASMANIAN_XSDK
}
TasmanianSparseGrid::TasmanianSparseGrid(const TasmanianSparseGrid &source) : base(0), global(0), sequence(0), pwpoly(0), wavelet(0),
                                    domain_transform_a(0), domain_transform_b(0), conformal_asin_power(0), acceleration(accel_none), gpuID(0), logstream(0)
{
    copyGrid(&source);
#ifndef TASMANIAN_XSDK
    logstream = &cerr;
#endif
#ifdef TASMANIAN_CPU_BLAS
    acceleration = accel_cpu_blas;
#endif // TASMANIAN_CPU_BLAS
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
    base = 0;
#ifndef TASMANIAN_XSDK
    logstream = &cerr;
#else
    logstream = 0;
#endif // TASMANIAN_XSDK
#ifdef TASMANIAN_CPU_BLAS
    acceleration = accel_cpu_blas;
#else
    acceleration = accel_none;
#endif // TASMANIAN_CPU_BLAS
#ifdef TASMANIAN_CUBLAS
    gpuID = 0;
#endif // TASMANIAN_CUBLAS
}

void TasmanianSparseGrid::setErrorLog(std::ostream *os){ logstream = os; }
void TasmanianSparseGrid::disableLog(){ logstream = 0; }

void TasmanianSparseGrid::write(const char *filename, bool binary) const{
    std::ofstream ofs;
    if (binary){
        ofs.open(filename, std::ios::out | std::ios::binary);
        writeBinary(ofs);
    }else{
        ofs.open(filename);
        writeAscii(ofs);
    }
    ofs.close();
}
bool TasmanianSparseGrid::read(const char *filename, bool binary){
    std::ifstream ifs;
    bool isGood;
    if (binary){
        ifs.open(filename, std::ios::in | std::ios::binary);
        isGood = readBinary(ifs);
    }else{
        ifs.open(filename);
        isGood = readAscii(ifs);
    }
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

void TasmanianSparseGrid::makeGlobalGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights, double alpha, double beta, const char* custom_filename){
    clear();
    global = new GridGlobal();
    global->makeGrid(dimensions, outputs, depth, type, rule, anisotropic_weights, alpha, beta, custom_filename);
    base = global;
}
void TasmanianSparseGrid::makeSequenceGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights){
    if (outputs < 1){
        if (logstream != 0){ (*logstream) << "ERROR: makeSequenceGrid is called with zero outputs, for zero outputs use makeGlobalGrid instead" << endl; }
        return;
    }
    if (OneDimensionalMeta::isSequence(rule)){
        clear();
        sequence = new GridSequence();
        sequence->makeGrid(dimensions, outputs, depth, type, rule, anisotropic_weights);
        base = sequence;
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: makeSequenceGrid is called with rule " << OneDimensionalMeta::getIORuleString(rule) << " which is not a sequence rule" << endl; }
    }
}
void TasmanianSparseGrid::makeLocalPolynomialGrid(int dimensions, int outputs, int depth, int order, TypeOneDRule rule){
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
    pwpoly->makeGrid(dimensions, outputs, depth, order, rule);
    base = pwpoly;
}
void TasmanianSparseGrid::makeWaveletGrid(int dimensions, int outputs, int depth, int order){
    if ((order != 1) && (order != 3)){
        if (logstream != 0){ (*logstream) << "ERROR: makeWaveletGrid is called with order " << order << ", but wavelets are implemented only for orders 1 and 3." << endl; }
        return;
    }
    clear();
    wavelet = new GridWavelet();
    wavelet->makeGrid(dimensions, outputs, depth, order);
    base = wavelet;
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
}

void TasmanianSparseGrid::updateGlobalGrid(int depth, TypeDepth type, const int *anisotropic_weights){
    if (global != 0){
        global->updateGrid(depth, type, anisotropic_weights);
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: updateGlobalGrid called, but the grid is not global" << endl; }
    }
}
void TasmanianSparseGrid::updateSequenceGrid(int depth, TypeDepth type, const int *anisotropic_weights){
    if (sequence != 0){
        sequence->updateGrid(depth, type, anisotropic_weights);
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
    mapConformalCanonicalToTransformed(base->getNumDimensions(), base->getNumLoaded(), x);
    if (domain_transform_a != 0){
        mapCanonicalToTransformed(base->getNumDimensions(), base->getNumLoaded(), base->getRule(), x);
    }
    return x;
}
void TasmanianSparseGrid::getLoadedPoints(double *x) const{
    base->getLoadedPoints(x);
    mapConformalCanonicalToTransformed(base->getNumDimensions(), base->getNumLoaded(), x);
    if (domain_transform_a != 0){
        mapCanonicalToTransformed(base->getNumDimensions(), base->getNumLoaded(), base->getRule(), x);
    }
}
double* TasmanianSparseGrid::getNeededPoints() const{
    double *x = base->getNeededPoints();
    mapConformalCanonicalToTransformed(base->getNumDimensions(), base->getNumNeeded(), x);
    if (domain_transform_a != 0){
        mapCanonicalToTransformed(base->getNumDimensions(), base->getNumNeeded(), base->getRule(), x);
    }
    return x;
}
void TasmanianSparseGrid::getNeededPoints(double *x) const{
    base->getNeededPoints(x);
    mapConformalCanonicalToTransformed(base->getNumDimensions(), base->getNumNeeded(), x);
    if (domain_transform_a != 0){
        mapCanonicalToTransformed(base->getNumDimensions(), base->getNumNeeded(), base->getRule(), x);
    }
}
double* TasmanianSparseGrid::getPoints() const{
    double *x = base->getPoints();
    mapConformalCanonicalToTransformed(base->getNumDimensions(), base->getNumPoints(), x);
    if (domain_transform_a != 0){
        mapCanonicalToTransformed(base->getNumDimensions(), base->getNumPoints(), base->getRule(), x);
    }
    return x;
}
void TasmanianSparseGrid::getPoints(double *x) const{
    base->getPoints(x);
    mapConformalCanonicalToTransformed(base->getNumDimensions(), base->getNumPoints(), x);
    if (domain_transform_a != 0){
        mapCanonicalToTransformed(base->getNumDimensions(), base->getNumPoints(), base->getRule(), x);
    }
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
    if ((domain_transform_a == 0) && (conformal_asin_power == 0)){
        return base->getInterpolationWeights(x);
    }else{
        int num_dimensions = base->getNumDimensions();
        double *x_canonical = new double[num_dimensions];  std::copy(x, x + num_dimensions, x_canonical);
        mapConformalTransformedToCanonical(num_dimensions, 1, x_canonical);
        if (domain_transform_a != 0) mapTransformedToCanonical(num_dimensions, base->getRule(), x_canonical);
        double *w = base->getInterpolationWeights(x_canonical);
        delete[] x_canonical;
        return w;
    }
}
void TasmanianSparseGrid::getInterpolationWeights(const double x[], double *weights) const{
    if ((domain_transform_a == 0) && (conformal_asin_power == 0)){
        base->getInterpolationWeights(x, weights);
    }else{
        int num_dimensions = base->getNumDimensions();
        double *x_canonical = new double[num_dimensions];  std::copy(x, x + num_dimensions, x_canonical);
        mapConformalTransformedToCanonical(num_dimensions, 1, x_canonical);
        if (domain_transform_a != 0) mapTransformedToCanonical(num_dimensions, base->getRule(), x_canonical);
        base->getInterpolationWeights(x_canonical, weights);
        delete[] x_canonical;
    }
}

void TasmanianSparseGrid::loadNeededPoints(const double *vals){ base->loadNeededPoints(vals); }

void TasmanianSparseGrid::evaluate(const double x[], double y[]) const{
    const double *x_effective = x;
    if ((domain_transform_a != 0) || (conformal_asin_power != 0)){
        int num_dimensions = base->getNumDimensions();
        double *x_canonical = new double[num_dimensions]; std::copy(x, x + num_dimensions, x_canonical);
        mapConformalTransformedToCanonical(num_dimensions, 1, x_canonical);
        if (domain_transform_a != 0) mapTransformedToCanonical(num_dimensions, base->getRule(), x_canonical);
        x_effective = x_canonical;
    }
    base->evaluate(x_effective, y);
    if (x_effective != x) delete[] x_effective;
}
void TasmanianSparseGrid::evaluateFast(const double x[], double y[]) const{
    const double *x_effective = x;
    if ((domain_transform_a != 0) || (conformal_asin_power != 0)){
        int num_dimensions = base->getNumDimensions();
        double *x_canonical = new double[num_dimensions]; std::copy(x, x + num_dimensions, x_canonical);
        mapConformalTransformedToCanonical(num_dimensions, 1, x_canonical);
        if (domain_transform_a != 0) mapTransformedToCanonical(num_dimensions, base->getRule(), x_canonical);
        x_effective = x_canonical;
    }
    switch (acceleration){
        case accel_gpu_default:
        case accel_gpu_fullmemory:
        case accel_gpu_cublas:
            #ifdef TASMANIAN_CUBLAS
            _TASMANIAN_SETGPU
            base->evaluateFastGPUcublas(x_effective, y, logstream);
            break;
            #endif // TASMANIAN_CUBLAS
        case accel_gpu_cuda:
            #ifdef TASMANIAN_CUDA
            _TASMANIAN_SETGPU
            base->evaluateFastGPUcuda(x_effective, y, logstream);
            break;
            #elif defined(TASMANIAN_CUBLAS)
            _TASMANIAN_SETGPU
            base->evaluateFastGPUcublas(x_effective, y, logstream);
            break;
            #endif // TASMANIAN_CUDA
        case accel_cpu_blas:
            #ifdef TASMANIAN_CPU_BLAS
            base->evaluateFastCPUblas(x_effective, y);
            break;
            #endif // TASMANIAN_CPU_BLAS
        default:
            base->evaluate(x_effective, y);
            break;
    }
    if (x_effective != x) delete[] x_effective;
}

void TasmanianSparseGrid::evaluateBatch(const double x[], int num_x, double y[]) const{
    const double *x_effective = x;
    if ((domain_transform_a != 0) || (conformal_asin_power != 0)){
        int num_dimensions = base->getNumDimensions();
        double *x_canonical = new double[num_dimensions*num_x]; std::copy(x, x + num_dimensions*num_x, x_canonical);
        mapConformalTransformedToCanonical(num_dimensions, num_x, x_canonical);
        if (domain_transform_a != 0) mapTransformedToCanonical(num_dimensions, num_x, base->getRule(), x_canonical);
        x_effective = x_canonical;
    }
    switch (acceleration){
        case accel_gpu_default:
        case accel_gpu_fullmemory:
        case accel_gpu_cublas:
            #ifdef TASMANIAN_CUBLAS
            _TASMANIAN_SETGPU
            base->evaluateBatchGPUcublas(x_effective, num_x, y, logstream);
            break;
            #endif // TASMANIAN_CUBLAS
        case accel_gpu_cuda:
            #ifdef TASMANIAN_CUDA
            _TASMANIAN_SETGPU
            base->evaluateBatchGPUcuda(x_effective, num_x, y, logstream);
            break;
            #elif defined(TASMANIAN_CUBLAS)
            _TASMANIAN_SETGPU
            base->evaluateBatchGPUcublas(x_effective, num_x, y, logstream);
            break;
            #endif // TASMANIAN_CUDA
        case accel_cpu_blas:
            #ifdef TASMANIAN_CPU_BLAS
            base->evaluateBatchCPUblas(x_effective, num_x, y);
            break;
            #endif // TASMANIAN_CPU_BLAS
        default:
            base->evaluateBatch(x_effective, num_x, y);
            break;
    }
    if (x_effective != x) delete[] x_effective;
}
void TasmanianSparseGrid::integrate(double q[]) const{
    if (conformal_asin_power != 0){
        int num_points = base->getNumPoints();
        double *correction = new double[num_points];  std::fill(correction, correction + num_points, 1.0);
        mapConformalWeights(base->getNumDimensions(), num_points, correction);
        base->integrate(q, correction);
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
    if (rule == rule_gausslaguerre){
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
void TasmanianSparseGrid::mapTransformedToCanonical(int num_dimensions, TypeOneDRule rule, double x[]) const{
    if (rule == rule_gausslaguerre){
        for(int j=0; j<num_dimensions; j++){
            x[j] -= domain_transform_a[j];
            x[j] *= domain_transform_b[j];
        }
    }else if ((rule == rule_gausshermite) || (rule == rule_gausshermiteodd)){
        for(int j=0; j<num_dimensions; j++){
            x[j] -= domain_transform_a[j];
            x[j] *= sqrt(domain_transform_b[j]);
        }
    }else{
        for(int j=0; j<num_dimensions; j++){
            x[j] *= 2.0;
            x[j] -= (domain_transform_b[j] + domain_transform_a[j]);
            x[j] /= (domain_transform_b[j] - domain_transform_a[j]);
        }
    }
}
void TasmanianSparseGrid::mapTransformedToCanonical(int num_dimensions, int num_points, TypeOneDRule rule, double x[]) const{
    if (rule == rule_gausslaguerre){
        for(int i=0; i<num_points; i++){
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
    }else if (rule == rule_gausslaguerre){
        for(int j=0; j<num_dimensions; j++) scale *= pow(domain_transform_b[j], global->getAlpha() + 1.0);
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

void TasmanianSparseGrid::setAnisotropicRefinement(TypeDepth type, int min_growth, int output){
    if (sequence != 0){
        sequence->setAnisotropicRefinement(type, min_growth);
    }else if (global != 0){
        if (OneDimensionalMeta::isNonNested(global->getRule())){
            if (logstream != 0){ (*logstream) << "ERROR: setAnisotropicRefinement called for a global grid with non-nested rule" << endl; }
        }else{
            global->setAnisotropicRefinement(type, min_growth, output);
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
void TasmanianSparseGrid::setSurplusRefinement(double tolerance, int output){
    if (sequence != 0){
        sequence->setSurplusRefinement(tolerance, output);
    }else if (global != 0){
        if (OneDimensionalMeta::isSequence(global->getRule())){
            global->setSurplusRefinement(tolerance, output);
        }else{
            if (logstream != 0){ (*logstream) << "ERROR: setSurplusRefinement called for a global grid with non-sequence rule" << endl; }
        }
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: setSurplusRefinement(double, int) called for grid that is neither sequence nor Global with sequence rule" << endl; }
    }
}
void TasmanianSparseGrid::setSurplusRefinement(double tolerance, TypeRefinement criteria, int output){
    if (pwpoly != 0){
        pwpoly->setSurplusRefinement(tolerance, criteria, output);
    }else if (wavelet != 0){
        wavelet->setSurplusRefinement(tolerance, criteria, output);
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: setSurplusRefinement(double, TypeRefinement) called for grid that is neither local polynomial nor wavelet" << endl; }
    }
}
void TasmanianSparseGrid::clearRefinement(){
    base->clearRefinement();
}
void TasmanianSparseGrid::removePointsBySurplus(double tolerance, int output){
    if (pwpoly == 0){
        if (logstream != 0){ (*logstream) << "ERROR: removePointsBySurplus() called for a grid that is not local polynomial." << endl; }
        return;
    }else{
        if (pwpoly->removePointsBySurplus(tolerance, output) == 0){
            clear();
        }
    }
}

double* TasmanianSparseGrid::evalHierarchicalFunctions(const double x[]) const{
    return base->evalHierarchicalFunctions(x);
}
void TasmanianSparseGrid::setHierarchicalCoefficients(const double c[]){
    base->setHierarchicalCoefficients(c);
}

int* TasmanianSparseGrid::getGlobalPolynomialSpace(bool interpolation, int &num_indexes) const{
    if (global != 0){
        return global->getPolynomialSpace(interpolation, num_indexes);
    }else if (sequence != 0){
        return sequence->getPolynomialSpace(interpolation, num_indexes);
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: getGlobalPolynomialSpace() called for a grid that is neither Global nor Sequence." << endl; }
        num_indexes = 0;
        return 0;
    }
}
const double* TasmanianSparseGrid::getSurpluses() const{
    if (pwpoly != 0){
        return pwpoly->getSurpluses();
    }else if (wavelet != 0){
        return wavelet->getSurpluses();
    }else if (sequence != 0){
        return sequence->getSurpluses();
    }else{
        if (logstream != 0){ (*logstream) << "ERROR: getSurplusses() called for a grid that is neither local polynomial nor wavelet nor sequence." << endl; }
        return 0;
    }
}
const int* TasmanianSparseGrid::getPointsIndexes() const{
    if (pwpoly != 0){
        //return ((pwpoly->getNumNeeded()>0) ? pwpoly->getNeededIndexes() : pwpoly->getPointIndexes());
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
        //cout << "Writing seq bin" << endl;
        flag = 's'; ofs.write(&flag, sizeof(char));
        sequence->writeBinary(ofs);
    }else if (pwpoly != 0){
        //cout << "Writing pwp bin" << endl;
        flag = 'p'; ofs.write(&flag, sizeof(char));
        pwpoly->writeBinary(ofs);
    }else if (wavelet != 0){
        //cout << "Writing wavelet bin" << endl;
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
        //cout << "Writing conformal" << endl;
        flag = 'a'; ofs.write(&flag, sizeof(char));
        ofs.write((char*) conformal_asin_power, base->getNumDimensions() * sizeof(int));
    }else{
        flag = 'n'; ofs.write(&flag, sizeof(char));
    }
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
                //cout << conformal_asin_power[j] << endl;
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
            if (!(T.compare("TASMANIAN SG end") == 0)){ if (logstream != 0){ (*logstream) << "WARNING: file did not end with 'TASMANIAN SG end', this may result in undefined behavior" << endl; } }
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
        //cout << "Reading seq bin" << endl;
        sequence = new GridSequence();
        sequence->readBinary(ifs);
        base = sequence;
    }else if (TSG[0] == 'p'){
        clear();
        //cout << "Reading pwp bin" << endl;
        pwpoly = new GridLocalPolynomial();
        pwpoly->readBinary(ifs);
        base = pwpoly;
    }else if (TSG[0] == 'w'){
        clear();
        //cout << "Reading wavelet bin" << endl;
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
        //cout << "Reading conformal" << endl;
        conformal_asin_power = new int[base->getNumDimensions()];
        ifs.read((char*) conformal_asin_power, base->getNumDimensions() * sizeof(int));
    }else if (TSG[0] != 'n'){
        if (logstream != 0){ (*logstream) << "ERROR: wrong conformal transform, wrong file format!" << endl; }
        return false;
    }
    delete[] TSG;
    return true;
}

void TasmanianSparseGrid::enableAcceleration(TypeAcceleration acc){
    if (acc != acceleration){
        if (base != 0) base->clearAccelerationData();
        acceleration = acc;
    }
}
TypeAcceleration TasmanianSparseGrid::getAccelerationType() const{
    return acceleration;
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
    #ifdef TASMANIAN_CUBLAS
    int gpu_count = 0;
    cudaGetDeviceCount(&gpu_count);
    return gpu_count;
    #else
    return 0;
    #endif // TASMANIAN_CUBLAS
}
int TasmanianSparseGrid::getGPUmemory(int gpu){
    #ifdef TASMANIAN_CUBLAS
    if (gpu < 0) return 0;
    int gpu_count = 0;
    cudaGetDeviceCount(&gpu_count);
    if (gpu >= gpu_count) return 0;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, gpu);
    unsigned long int memB = prop.totalGlobalMem;
    return (int) (memB / (1048576));
    #else
    return 0;
    #endif // TASMANIAN_CUBLAS
}
const char* TasmanianSparseGrid::getGPUname(int gpu){
    #ifdef TASMANIAN_CUBLAS
    if (gpu < 0) return "";
    int gpu_count = 0;
    cudaGetDeviceCount(&gpu_count);
    if (gpu >= gpu_count) return "";
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, gpu);
    char *name = new char[257];
    for(int i=0; i<256; i++) name[i] = prop.name[i];
    name[256] = '\0';
    return name;
    #else
    return "";
    #endif // TASMANIAN_CUBLAS
}

// ------------ C Interface for use with Python ctypes and potentially other C codes -------------- //
extern "C" {
void* tsgConstructTasmanianSparseGrid(){ return (void*) new TasmanianSparseGrid(); }
void tsgDestructTasmanianSparseGrid(void *grid){ delete ((TasmanianSparseGrid*) grid); }

void tsgCopyGrid(void *destination, void *source){ ((TasmanianSparseGrid*) destination)->copyGrid(((TasmanianSparseGrid*) source)); }

const char* tsgGetVersion(){ return TasmanianSparseGrid::getVersion(); }
const char* tsgGetLicense(){ return TasmanianSparseGrid::getLicense(); }
int tsgGetVersionMajor(){ return TasmanianSparseGrid::getVersionMajor(); }
int tsgGetVersionMinor(){ return TasmanianSparseGrid::getVersionMinor(); }
int tsgIsCudaEnabled(){ return (TasmanianSparseGrid::isCudaEnabled()) ? 0 : 1; }
int tsgIsBLASEnabled(){ return (TasmanianSparseGrid::isBLASEnabled()) ? 0 : 1; }
int tsgIsOpenMPEnabled(){ return (TasmanianSparseGrid::isOpenMPEnabled()) ? 0 : 1; }

void tsgErrorLogCerr(void *grid){ ((TasmanianSparseGrid*) grid)->setErrorLog(&cerr); }
void tsgDisableErrorLog(void *grid){ ((TasmanianSparseGrid*) grid)->disableLog(); }

void tsgWrite(void *grid, const char* filename){ ((TasmanianSparseGrid*) grid)->write(filename); }
int tsgRead(void *grid, const char* filename){
    bool result = ((TasmanianSparseGrid*) grid)->read(filename);
    return result ? 0 : 1;
}
void tsgWriteBinary(void *grid, const char* filename){ ((TasmanianSparseGrid*) grid)->write(filename, true); }
int tsgReadBinary(void *grid, const char* filename){
    bool result = ((TasmanianSparseGrid*) grid)->read(filename, true);
    return result ? 0 : 1;
}

void tsgMakeGlobalGrid(void *grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights, double alpha, double beta, const char* custom_filename){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    TypeOneDRule rule = OneDimensionalMeta::getIORuleString(sRule);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    if (rule == rule_none){ cerr << "WARNING: incorrect rule type: " << sType << ", defaulting to clenshaw-curtis." << endl; }
    #endif // _TASMANIAN_DEBUG_
    ((TasmanianSparseGrid*) grid)->makeGlobalGrid(dimensions, outputs, depth, depth_type, rule, anisotropic_weights, alpha, beta, custom_filename);
}
void tsgMakeSequenceGrid(void *grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    TypeOneDRule rule = OneDimensionalMeta::getIORuleString(sRule);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    if (rule == rule_none){ cerr << "WARNING: incorrect rule type: " << sRule << ", defaulting to clenshaw-curtis." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    if (rule == rule_none){ rule = rule_clenshawcurtis; }
    ((TasmanianSparseGrid*) grid)->makeSequenceGrid(dimensions, outputs, depth, depth_type, rule, anisotropic_weights);
}
void tsgMakeLocalPolynomialGrid(void *grid, int dimensions, int outputs, int depth, int order, const char * sRule){
    TypeOneDRule rule = OneDimensionalMeta::getIORuleString(sRule);
    #ifdef _TASMANIAN_DEBUG_
    if (rule == rule_none){ cerr << "WARNING: incorrect rule type: " << sRule << ", defaulting to localp." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (rule == rule_none){ rule = rule_localp; }
    ((TasmanianSparseGrid*) grid)->makeLocalPolynomialGrid(dimensions, outputs, depth, order, rule);
}
void tsgMakeWaveletGrid(void *grid, int dimensions, int outputs, int depth, int order){
    ((TasmanianSparseGrid*) grid)->makeWaveletGrid(dimensions, outputs, depth, order);
}

void tsgUpdateGlobalGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->updateGlobalGrid(depth, depth_type, anisotropic_weights);
}
void tsgUpdateSequenceGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->updateSequenceGrid(depth, depth_type, anisotropic_weights);
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

double* tsgBatchGetInterpolationWeights(void *grid, const double *x, int num_x){
    TasmanianSparseGrid* tsg = (TasmanianSparseGrid*) grid;
    int iNumDim = tsg->getNumDimensions(), iNumPoints = tsg->getNumPoints();
    double *weights = (double*) malloc(num_x * iNumPoints * sizeof(double));
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        double *w = tsg->getInterpolationWeights(&(x[i*iNumDim]));
        std::copy(w, w + iNumPoints, &(weights[i*iNumPoints]));
        delete[] w;
    }
    return weights;
}

int tsgIsGlobal(void *grid){ return (((TasmanianSparseGrid*) grid)->isGlobal() ? 0 : 1); }
int tsgIsSequence(void *grid){ return (((TasmanianSparseGrid*) grid)->isSequence() ? 0 : 1); }
int tsgIsLocalPolynomial(void *grid){ return (((TasmanianSparseGrid*) grid)->isLocalPolynomial() ? 0 : 1); }
int tsgIsWavelet(void *grid){ return (((TasmanianSparseGrid*) grid)->isWavelet() ? 0 : 1); }

void tsgSetDomainTransform(void *grid, const double a[], const double b[]){ ((TasmanianSparseGrid*) grid)->setDomainTransform(a, b); }
int tsgIsSetDomainTransfrom(void *grid){ return (((TasmanianSparseGrid*) grid)->isSetDomainTransfrom() ? 0 : 1); }
void tsgClearDomainTransform(void *grid){ ((TasmanianSparseGrid*) grid)->clearDomainTransform(); }
void tsgGetDomainTransform(void *grid, double a[], double b[]){ ((TasmanianSparseGrid*) grid)->getDomainTransform(a, b); }

void tsgSetConformalTransformASIN(void *grid, const int truncation[]){ ((TasmanianSparseGrid*) grid)->setConformalTransformASIN(truncation); }
int tsgIsSetConformalTransformASIN(void *grid){ return (((TasmanianSparseGrid*) grid)->isSetConformalTransformASIN()) ? 0 : 1; }
void tsgClearConformalTransform(void *grid){ ((TasmanianSparseGrid*) grid)->clearConformalTransform(); }
void tsgGetConformalTransformASIN(void *grid, int truncation[]){ ((TasmanianSparseGrid*) grid)->getConformalTransformASIN(truncation); }

void tsgSetAnisotropicRefinement(void *grid, const char * sType, int min_growth, int output){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->setAnisotropicRefinement(depth_type, min_growth, output);
}
int* tsgEstimateAnisotropicCoefficients(void *grid, const char * sType, int output, int *num_coefficients){
    TypeDepth depth_type = OneDimensionalMeta::getIOTypeString(sType);
    #ifdef _TASMANIAN_DEBUG_
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (depth_type == type_none){ depth_type = type_iptotal; }
    if ((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)){
        *num_coefficients = 2 * (((TasmanianSparseGrid*) grid)->getNumDimensions());
    }else{
        *num_coefficients = ((TasmanianSparseGrid*) grid)->getNumDimensions();
    }
    int *coeff = ((TasmanianSparseGrid*) grid)->estimateAnisotropicCoefficients(depth_type, output);
    int *result = (int*) malloc((*num_coefficients) * sizeof(int));
    for(int i=0; i<*num_coefficients; i++) result[i] = coeff[i];
    delete[] coeff;
    return result;
}
void tsgSetGlobalSurplusRefinement(void *grid, double tolerance, int output){
    ((TasmanianSparseGrid*) grid)->setSurplusRefinement(tolerance, output);
}
void tsgSetLocalSurplusRefinement(void *grid, double tolerance, const char * sRefinementType, int output){
    TypeRefinement ref_type = OneDimensionalMeta::getIOTypeRefinementString(sRefinementType);
    #ifdef _TASMANIAN_DEBUG_
    if (ref_type == refine_none){ cerr << "WARNING: incorrect refinement type: " << sRefinementType << ", defaulting to type_classic." << endl; }
    #endif // _TASMANIAN_DEBUG_
    if (ref_type == refine_none){ ref_type = refine_classic; }
    ((TasmanianSparseGrid*) grid)->setSurplusRefinement(tolerance, ref_type, output);
}
void tsgClearRefinement(void *grid){
    ((TasmanianSparseGrid*) grid)->clearRefinement();
}
void tsgRemovePointsBySurplus(void *grid, double tolerance, int output){
    ((TasmanianSparseGrid*) grid)->removePointsBySurplus(tolerance, output);
}

double* tsgEvalHierarchicalFunctions(void *grid, const double *x){
    return ((TasmanianSparseGrid*) grid)->evalHierarchicalFunctions(x);
}
double* tsgBatchEvalHierarchicalFunctions(void *grid, const double *x, int num_x){
    TasmanianSparseGrid* tsg = (TasmanianSparseGrid*) grid;
    int iNumDim = tsg->getNumDimensions(), iNumPoints = tsg->getNumPoints();
    double *vals = new double[num_x * iNumPoints];
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        double *v = tsg->evalHierarchicalFunctions(&(x[i*iNumDim]));
        std::copy(v, v + iNumPoints, &(vals[i*iNumPoints]));
        delete[] v;
   }
        return vals;
}
void tsgSetHierarchicalCoefficients(void *grid, const double *c){
    ((TasmanianSparseGrid*) grid)->setHierarchicalCoefficients(c);
}
const double* tsgGetSurpluses(void *grid){
    return ((TasmanianSparseGrid*) grid)->getSurpluses();
}

int* tsgGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes){
    int ni, *idx = ((TasmanianSparseGrid*) grid)->getGlobalPolynomialSpace((interpolation == 0), ni);
    *num_indexes = ni;
    return idx;
}

void tsgPrintStats(void *grid){ ((TasmanianSparseGrid*) grid)->printStats(); }

void tsgEnableAcceleration(void *grid, const char *accel){ ((TasmanianSparseGrid*) grid)->enableAcceleration(AccelerationMeta::getIOAccelerationString(accel)); }
int tsgGetAccelerationType(void *grid){ return AccelerationMeta::getIOAccelerationInt(((TasmanianSparseGrid*) grid)->getAccelerationType()); } // int to acceleration type

void tsgSetGPUID(void *grid, int gpuID){ ((TasmanianSparseGrid*) grid)->setGPUID(gpuID); }
int tsgGetGPUID(void *grid){ return ((TasmanianSparseGrid*) grid)->getGPUID(); }
int tsgGetNumGPUs(){ return TasmanianSparseGrid::getNumGPUs(); }
int tsgGetGPUmemory(int gpu){ return TasmanianSparseGrid::getGPUmemory(gpu); }
const char* tsgGetGPUname(int gpu){ return TasmanianSparseGrid::getGPUname(gpu); }

void tsgDeleteDoubles(double *p){ free(p); }
void tsgDeleteInts(int *p){ delete[] p; }
void tsgDeleteChars(char *p){ delete[] p; }
}
}

#endif
