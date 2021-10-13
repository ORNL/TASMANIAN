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

#ifndef __TASMANIAN_SPARSE_GRID_WRAPC_CPP
#define __TASMANIAN_SPARSE_GRID_WRAPC_CPP

#include "TasmanianSparseGrid.hpp"

// ------------ C Interface for use with Python ctypes and potentially other C codes -------------- //
using std::cerr;
using std::endl;

namespace TasGrid{

extern "C" {
void* tsgConstructTasmanianSparseGrid(){ return (void*) new TasmanianSparseGrid(); }
void tsgDestructTasmanianSparseGrid(void *grid){ delete ((TasmanianSparseGrid*) grid); }

void tsgCopyGrid(void *destination, void *source){ ((TasmanianSparseGrid*) destination)->copyGrid(((TasmanianSparseGrid*) source)); }
void tsgCopySubGrid(void *destination, void *source, int outputs_begin, int outputs_end){
    ((TasmanianSparseGrid*) destination)->copyGrid(((TasmanianSparseGrid*) source), outputs_begin, outputs_end);
}

const char* tsgGetVersion(){ return TasmanianSparseGrid::getVersion(); }
const char* tsgGetLicense(){ return TasmanianSparseGrid::getLicense(); }
int tsgGetVersionMajor(){ return TasmanianSparseGrid::getVersionMajor(); }
int tsgGetVersionMinor(){ return TasmanianSparseGrid::getVersionMinor(); }
int tsgIsOpenMPEnabled(){ return (TasmanianSparseGrid::isOpenMPEnabled()) ? 1 : 0; }
int tsgIsCudaEnabled(){ return (TasmanianSparseGrid::isCudaEnabled()) ? 1 : 0; }
int tsgIsHipEnabled(){ return (TasmanianSparseGrid::isHipEnabled()) ? 1 : 0; }
int tsgIsDpcppEnabled(){ return (TasmanianSparseGrid::isDpcppEnabled()) ? 1 : 0; }

void tsgWrite(void *grid, const char* filename){ ((TasmanianSparseGrid*) grid)->write(filename, mode_ascii); }
void tsgWriteBinary(void *grid, const char* filename){ ((TasmanianSparseGrid*) grid)->write(filename, mode_binary); }
int tsgRead(void *grid, const char* filename){
    try{
        ((TasmanianSparseGrid*) grid)->read(filename);
        return 1;
    }catch(std::runtime_error &e){
        cerr << e.what() << endl;
        return 0;
    }
}

void tsgMakeGlobalGrid(void *grid, int dimensions, int outputs, int depth, const char * sType, const char *sRule, const int *anisotropic_weights, double alpha, double beta, const char* custom_filename, const int *limit_levels){
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    TypeOneDRule rule = IO::getRuleString(sRule);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    if (rule == rule_none){ cerr << "WARNING: incorrect rule type: " << sType << ", defaulting to clenshaw-curtis." << endl; }
    #endif // NDEBUG
    ((TasmanianSparseGrid*) grid)->makeGlobalGrid(dimensions, outputs, depth, depth_type, rule, anisotropic_weights, alpha, beta, custom_filename, limit_levels);
}
void tsgMakeSequenceGrid(void *grid, int dimensions, int outputs, int depth, const char *sType, const char *sRule, const int *anisotropic_weights, const int *limit_levels){
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    TypeOneDRule rule = IO::getRuleString(sRule);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    if (rule == rule_none){ cerr << "WARNING: incorrect rule type: " << sRule << ", defaulting to clenshaw-curtis." << endl; }
    #endif // NDEBUG
    if (depth_type == type_none){ depth_type = type_iptotal; }
    if (rule == rule_none){ rule = rule_clenshawcurtis; }
    ((TasmanianSparseGrid*) grid)->makeSequenceGrid(dimensions, outputs, depth, depth_type, rule, anisotropic_weights, limit_levels);
}
void tsgMakeLocalPolynomialGrid(void *grid, int dimensions, int outputs, int depth, int order, const char *sRule, const int *limit_levels){
    TypeOneDRule rule = IO::getRuleString(sRule);
    #ifndef NDEBUG
    if (rule == rule_none){ cerr << "WARNING: incorrect rule type: " << sRule << ", defaulting to localp." << endl; }
    #endif // NDEBUG
    if (rule == rule_none){ rule = rule_localp; }
    ((TasmanianSparseGrid*) grid)->makeLocalPolynomialGrid(dimensions, outputs, depth, order, rule, limit_levels);
}
void tsgMakeWaveletGrid(void *grid, int dimensions, int outputs, int depth, int order, const int *limit_levels){
    ((TasmanianSparseGrid*) grid)->makeWaveletGrid(dimensions, outputs, depth, order, limit_levels);
}
void tsgMakeFourierGrid(void *grid, int dimensions, int outputs, int depth, const char *sType, const int *anisotropic_weights, const int *limit_levels){
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_level." << endl; }
    #endif // NDEBUG
    if (depth_type == type_none){ depth_type = type_level; }
    ((TasmanianSparseGrid*) grid)->makeFourierGrid(dimensions, outputs, depth, depth_type, anisotropic_weights, limit_levels);
}
void tsgMakeGridFromCustomTabulated(void *grid, int dimension, int outputs, int depth, const char *sType, void *custom_tabulated,
                                    const int *anisotropic_weights, const int *limit_levels) {
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // NDEBUG
    ((TasmanianSparseGrid*) grid)->makeGlobalGrid(dimension, outputs, depth, depth_type,
                                                  CustomTabulated(*reinterpret_cast<TasGrid::CustomTabulated*>(custom_tabulated)),
                                                  anisotropic_weights, limit_levels);
}

void tsgUpdateGlobalGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels){
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // NDEBUG
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->updateGlobalGrid(depth, depth_type, anisotropic_weights, limit_levels);
}
void tsgUpdateSequenceGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels){
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // NDEBUG
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->updateSequenceGrid(depth, depth_type, anisotropic_weights, limit_levels);
}
void tsgUpdateFourierGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels){
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // NDEBUG
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->updateFourierGrid(depth, depth_type, anisotropic_weights, limit_levels);
}

double tsgGetAlpha(void *grid){ return ((TasmanianSparseGrid*) grid)->getAlpha(); }
double tsgGetBeta(void *grid){ return ((TasmanianSparseGrid*) grid)->getBeta(); }
int tsgGetOrder(void *grid){ return ((TasmanianSparseGrid*) grid)->getOrder(); }
int tsgGetNumDimensions(void *grid){ return ((TasmanianSparseGrid*) grid)->getNumDimensions(); }
int tsgGetNumOutputs(void *grid){ return ((TasmanianSparseGrid*) grid)->getNumOutputs(); }
char* tsgGetRule(void *grid){
    std::string cppstring = IO::getRuleString( ((TasmanianSparseGrid*) grid)->getRule() );
    char *cstring = new char[cppstring.size() + 1];
    for(size_t i=0; i<cppstring.size(); i++) cstring[i] = cppstring[i];
    cstring[cppstring.size()] = '\0';
    return cstring;
}
void tsgCopyRuleChars(void *grid, int buffer_size, char *name, int *num_actual){
    std::string cppstring = IO::getRuleString( ((TasmanianSparseGrid*) grid)->getRule() );
    size_t max_num = std::min((size_t) buffer_size - 1, cppstring.size());
    std::copy_n(cppstring.begin(), max_num, name);
    name[max_num] = '\0';
    *num_actual = (int) max_num;
}
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

void tsgLoadNeededValues(void *grid, const double *vals){ ((TasmanianSparseGrid*) grid)->loadNeededValues(vals); }

const double* tsgGetLoadedValues(void *grid){
    return ((TasmanianSparseGrid*) grid)->getLoadedValues();
}
void tsgGetLoadedValuesStatic(void *grid, double *values){
    int num_points = ((TasmanianSparseGrid*) grid)->getNumPoints();
    int num_outputs = ((TasmanianSparseGrid*) grid)->getNumOutputs();
    if ((num_points == 0) || (num_outputs == 0)) return;
    const double *vals = ((TasmanianSparseGrid*) grid)->getLoadedValues();
    std::copy(vals, vals + Utils::size_mult(num_outputs, num_points), values);
}

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
int tsgIsFourier(void *grid){ return (((TasmanianSparseGrid*) grid)->isFourier() ? 1 : 0); }

void tsgSetDomainTransform(void *grid, const double a[], const double b[]){ ((TasmanianSparseGrid*) grid)->setDomainTransform(a, b); }
int tsgIsSetDomainTransfrom(void *grid){ return (((TasmanianSparseGrid*) grid)->isSetDomainTransfrom() ? 1 : 0); }
void tsgClearDomainTransform(void *grid){ ((TasmanianSparseGrid*) grid)->clearDomainTransform(); }
void tsgGetDomainTransform(void *grid, double a[], double b[]){ ((TasmanianSparseGrid*) grid)->getDomainTransform(a, b); }

void tsgSetConformalTransformASIN(void *grid, const int truncation[]){
    ((TasmanianSparseGrid*) grid)->setConformalTransformASIN(Utils::copyArray(truncation, ((TasmanianSparseGrid*) grid)->getNumDimensions()));
}
int tsgIsSetConformalTransformASIN(void *grid){ return (((TasmanianSparseGrid*) grid)->isSetConformalTransformASIN()) ? 1 : 0; }
void tsgClearConformalTransform(void *grid){ ((TasmanianSparseGrid*) grid)->clearConformalTransform(); }
void tsgGetConformalTransformASIN(void *grid, int truncation[]){
    auto truncation_vector = ((TasmanianSparseGrid*) grid)->getConformalTransformASIN();
    std::copy(truncation_vector.begin(), truncation_vector.end(), truncation);
}

void tsgClearLevelLimits(void *grid){ ((TasmanianSparseGrid*) grid)->clearLevelLimits(); }
void tsgGetLevelLimits(void *grid, int *limits){
    auto llimits = ((TasmanianSparseGrid*) grid)->getLevelLimits();
    if (llimits.empty()){
        std::fill_n(limits, ((TasmanianSparseGrid*) grid)->getNumDimensions(), -1);
    }else{
        std::copy(llimits.begin(), llimits.end(), limits);
    }
}

void tsgSetAnisotropicRefinement(void *grid, const char * sType, int min_growth, int output, const int *level_limits){
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // NDEBUG
    if (depth_type == type_none){ depth_type = type_iptotal; }
    ((TasmanianSparseGrid*) grid)->setAnisotropicRefinement(depth_type, min_growth, output, level_limits);
}
int* tsgEstimateAnisotropicCoefficients(void *grid, const char * sType, int output, int *num_coefficients){
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // NDEBUG
    if (depth_type == type_none){ depth_type = type_iptotal; }
    *num_coefficients = ((TasmanianSparseGrid*) grid)->getNumDimensions();
    if ((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)){
        *num_coefficients *= 2;
    }
    auto coeff = ((TasmanianSparseGrid*) grid)->estimateAnisotropicCoefficients(depth_type, output);
    int *result = (int*) malloc((*num_coefficients) * sizeof(int));
    for(int i=0; i<*num_coefficients; i++) result[i] = coeff[i];
    return result;
}
void tsgEstimateAnisotropicCoefficientsStatic(void *grid, const char * sType, int output, int *coefficients){
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // NDEBUG
    if (depth_type == type_none){ depth_type = type_iptotal; }
    int num_coefficients = ((TasmanianSparseGrid*) grid)->getNumDimensions();
    if ((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)){
        num_coefficients *= 2;
    }
    auto coeff = ((TasmanianSparseGrid*) grid)->estimateAnisotropicCoefficients(depth_type, output);
    for(int i=0; i<num_coefficients; i++) coefficients[i] = coeff[i];
}
void tsgSetGlobalSurplusRefinement(void *grid, double tolerance, int output, const int *level_limits){
    ((TasmanianSparseGrid*) grid)->setSurplusRefinement(tolerance, output, level_limits);
}
void tsgSetLocalSurplusRefinement(void *grid, double tolerance, const char * sRefinementType, int output, const int *level_limits, const double *scale_correction){
    TypeRefinement ref_type = IO::getTypeRefinementString(sRefinementType);
    #ifndef NDEBUG
    if (ref_type == refine_none){ cerr << "WARNING: incorrect refinement type: " << sRefinementType << ", defaulting to type_classic." << endl; }
    #endif // NDEBUG
    if (ref_type == refine_none){ ref_type = refine_classic; }
    ((TasmanianSparseGrid*) grid)->setSurplusRefinement(tolerance, ref_type, output, level_limits, scale_correction);
}
void tsgClearRefinement(void *grid){
    ((TasmanianSparseGrid*) grid)->clearRefinement();
}
void tsgMergeRefinement(void *grid){
    ((TasmanianSparseGrid*) grid)->mergeRefinement();
}
void tsgBeginConstruction(void *grid){
    ((TasmanianSparseGrid*) grid)->beginConstruction();
}
int tsgIsUsingConstruction(void *grid){
    return (((TasmanianSparseGrid*) grid)->isUsingConstruction()) ? 1 : 0;
}
void* tsgGetCandidateConstructionPointsVoidPntr(void *grid, const char *sType, int output, const int *anisotropic_weights, const int *limit_levels){ // internal use only
    TypeDepth depth_type = IO::getDepthTypeString(sType);
    #ifndef NDEBUG
    if (depth_type == type_none){ cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl; }
    #endif // NDEBUG
    if (depth_type == type_none){ depth_type = type_iptotal; }
    size_t dims = (size_t) ((TasmanianSparseGrid*) grid)->getNumDimensions();
    std::vector<double>* vecx = (std::vector<double>*) new std::vector<double>();
    std::vector<int> veclimits;
    if (limit_levels != nullptr) veclimits = std::vector<int>(limit_levels, limit_levels + dims);
    if (anisotropic_weights == nullptr){
        *vecx = ((TasmanianSparseGrid*) grid)->getCandidateConstructionPoints(depth_type, output, veclimits);
    }else{
        std::vector<int> vecweights(anisotropic_weights, anisotropic_weights +
                                    (((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)) ? 2*dims : dims));
        *vecx = ((TasmanianSparseGrid*) grid)->getCandidateConstructionPoints(depth_type, vecweights, veclimits);
    }
    return (void*) vecx;
}
void* tsgGetCandidateConstructionPointsSurplusVoidPntr(void *grid, double tolerance, const char *sRefType, int output, const int *limit_levels, const double *scale_correction){ // internal use only
    TypeRefinement ref_type = IO::getTypeRefinementString(sRefType);
    #ifndef NDEBUG
    if (ref_type == refine_none){ cerr << "WARNING: incorrect depth type: " << sRefType << ", defaulting to refine_classic." << endl; }
    #endif // NDEBUG
    if (ref_type == refine_none){ ref_type = refine_classic; }
    size_t dims = (size_t) ((TasmanianSparseGrid*) grid)->getNumDimensions();
    std::vector<double>* vecx = (std::vector<double>*) new std::vector<double>();
    std::vector<int> veclimits;
    if (limit_levels != nullptr) veclimits = std::vector<int>(limit_levels, limit_levels + dims);
    std::vector<double> vecscale;
    if (scale_correction != nullptr){
        size_t active_outputs = (size_t) (output == -1) ? ((TasmanianSparseGrid*) grid)->getNumOutputs() : 1;
        vecscale = std::vector<double>(scale_correction, scale_correction + ((size_t) ((TasmanianSparseGrid*) grid)->getNumLoaded() * active_outputs));
    }
    *vecx = ((TasmanianSparseGrid*) grid)->getCandidateConstructionPoints(tolerance, ref_type, output, veclimits, vecscale);
    return (void*) vecx;
}
void tsgGetCandidateConstructionPoints(void *grid, const char *sType, int output, const int *anisotropic_weights, const int *limit_levels, int *num_points, double **x){
    size_t dims = (size_t) ((TasmanianSparseGrid*) grid)->getNumDimensions();
    std::vector<double>* vecx = (std::vector<double>*) tsgGetCandidateConstructionPointsVoidPntr(grid, sType, output, anisotropic_weights, limit_levels);
    *num_points = (int)(vecx->size() / dims);
    *x = (double*) malloc(vecx->size() * sizeof(double));
    std::copy_n(vecx->data(), vecx->size(), *x);
    delete vecx;
}
void tsgGetCandidateConstructionSurplusPoints(void *grid, double tolerance, const char *sRefType, int output, const int *limit_levels, const double *scale_correction, int *num_points, double **x){
    size_t dims = (size_t) ((TasmanianSparseGrid*) grid)->getNumDimensions();
    std::vector<double>* vecx = (std::vector<double>*) tsgGetCandidateConstructionPointsSurplusVoidPntr(grid, tolerance, sRefType, output, limit_levels, scale_correction);
    *num_points = (int)(vecx->size() / dims);
    *x = (double*) malloc(vecx->size() * sizeof(double));
    std::copy_n(vecx->data(), vecx->size(), *x);
    delete vecx;
}
int tsgGetCandidateConstructionPointsPythonGetNP(void *grid, const void *vecx){
    return (int) (((std::vector<double>*) vecx)->size() / ((size_t) ((TasmanianSparseGrid*) grid)->getNumDimensions()));
}
void tsgGetCandidateConstructionPointsPythonStatic(const void *vecx, double *x){ std::copy_n(((std::vector<double>*) vecx)->data(), ((std::vector<double>*) vecx)->size(), x); }
void tsgGetCandidateConstructionPointsPythonDeleteVect(void *vecx){ delete ((std::vector<double>*) vecx); }
void tsgLoadConstructedPoint(void *grid, const double *x, int numx, const double *y){
    ((TasmanianSparseGrid*) grid)->loadConstructedPoints(x, numx, y);
}
void tsgFinishConstruction(void *grid){
    ((TasmanianSparseGrid*) grid)->finishConstruction();
}

void tsgRemovePointsByHierarchicalCoefficient(void *grid, double tolerance, int output, const double *scale_correction){
    ((TasmanianSparseGrid*) grid)->removePointsByHierarchicalCoefficient(tolerance, output, scale_correction);
}
void tsgRemovePointsByHierarchicalCoefficientHardCutoff(void *grid, int num_new, int output, const double *scale_correction){
    ((TasmanianSparseGrid*) grid)->removePointsByHierarchicalCoefficient(num_new, output, scale_correction);
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
void tsgGetHierarchicalSupportStatic(void *grid, double support[]){
    std::vector<double> sup = ((TasmanianSparseGrid*) grid)->getHierarchicalSupport();
    std::copy(sup.begin(), sup.end(), support);
}
const double* tsgGetHierarchicalCoefficients(void *grid){
    return ((TasmanianSparseGrid*) grid)->getHierarchicalCoefficients();
}
void tsgGetHierarchicalCoefficientsStatic(void *grid, double *coeff){
    ((TasmanianSparseGrid*) grid)->getHierarchicalCoefficientsStatic(coeff);
}
void tsgSetHierarchicalCoefficients(void *grid, const double *c){
    ((TasmanianSparseGrid*) grid)->setHierarchicalCoefficients(c);
}
double* tsgIntegrateHierarchicalFunctions(void *grid){
    double *x = (double*) malloc(((TasmanianSparseGrid*) grid)->getNumPoints() * sizeof(double));
    ((TasmanianSparseGrid*) grid)->integrateHierarchicalFunctions(x);
    return x;
}
void tsgIntegrateHierarchicalFunctionsStatic(void *grid, double *integrals){
    ((TasmanianSparseGrid*) grid)->integrateHierarchicalFunctions(integrals);
}

// to be called from Python only, must later call delete[] on the pointer
int* tsgPythonGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes){
    std::vector<int> space = ((TasmanianSparseGrid*) grid)->getGlobalPolynomialSpace((interpolation != 0));
    int *indx = new int[space.size()];
    std::copy(space.begin(), space.end(), indx);
    *num_indexes = (int) space.size() / ((TasmanianSparseGrid*) grid)->getNumDimensions();
    return indx;
}
// to be used in C, creates a C pointer (requires internal copy of data)
void tsgGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes, int **indexes){
    std::vector<int> space = ((TasmanianSparseGrid*) grid)->getGlobalPolynomialSpace((interpolation != 0));
    *num_indexes = (int) space.size() / ((TasmanianSparseGrid*) grid)->getNumDimensions();
    if (!space.empty()){
        *indexes = (int*) malloc(space.size() * sizeof(int));
        std::copy(space.begin(), space.end(), *indexes);
    }
}

void tsgPrintStats(void *grid){ ((TasmanianSparseGrid*) grid)->printStats(); }

void tsgEnableAcceleration(void *grid, const char *accel){ ((TasmanianSparseGrid*) grid)->enableAcceleration(AccelerationMeta::getIOAccelerationString(accel)); }
void tsgEnableAccelerationGPU(void *grid, const char *accel, int gpu){ ((TasmanianSparseGrid*) grid)->enableAcceleration(AccelerationMeta::getIOAccelerationString(accel), gpu); }
//int tsgGetAccelerationTypeInt(void *grid){ return AccelerationMeta::getIOAccelerationInt(((TasmanianSparseGrid*) grid)->getAccelerationType()); } // int to acceleration type
const char* tsgGetAccelerationType(void *grid){ return AccelerationMeta::getIOAccelerationString(((TasmanianSparseGrid*) grid)->getAccelerationType()); }

void tsgSetGPUID(void *grid, int gpuID){ ((TasmanianSparseGrid*) grid)->setGPUID(gpuID); }
int tsgGetGPUID(void *grid){ return ((TasmanianSparseGrid*) grid)->getGPUID(); }
int tsgGetNumGPUs(){ return TasmanianSparseGrid::getNumGPUs(); }
int tsgGetGPUMemory(int gpu){ return TasmanianSparseGrid::getGPUMemory(gpu); }
int tsgIsAccelerationAvailable(const char *accel){ return (TasmanianSparseGrid::isAccelerationAvailable(AccelerationMeta::getIOAccelerationString(accel))) ? 1 : 0; }
void tsgGetGPUName(int gpu, int num_buffer, char *buffer, int *num_actual){
    // gpu is the gpuID, num_buffer is the size of *buffer, num_actual returns the actual number of chars
    if (num_buffer == 0) return;
    std::string name = TasmanianSparseGrid::getGPUName(gpu);

    size_t chars = std::min((size_t) (num_buffer - 1), name.size());
    std::copy(name.begin(), name.begin() + chars, buffer);
    buffer[chars] = '\0';

    *num_actual = (int) chars;
}

void tsgDeleteInts(int *p){ delete[] p; }

void* tsgConstructCustomTabulated(){ return (void*) new CustomTabulated(); }
void tsgDestructCustomTabulated(void* ct){ delete ((CustomTabulated*) ct); }

void tsgWriteCustomTabulated(void *ct, const char* filename){
    std::ofstream ofs(filename, std::ios::out);
    if (!ofs.good()) std::cerr << "ERROR: must provide valid filename!" << std::endl;
    ((CustomTabulated*) ct)->write<false>(ofs); // false == mode_ascii
}
int tsgReadCustomTabulated(void *ct, const char* filename){
    try{
        ((CustomTabulated*) ct)->read(filename);
        return 1;
    }catch(std::runtime_error &e){
        cerr << e.what() << endl;
        return 0;
    }catch(std::invalid_argument &e){
        cerr << e.what() << endl;
        return 0;
    }
}

int tsgGetNumLevelsCustomTabulated(void* ct){ return ((CustomTabulated*) ct)->getNumLevels(); }
int tsgGetNumPointsCustomTabulated(void* ct, const int level){ return ((CustomTabulated*) ct)->getNumPoints(level); }
int tsgGetIExactCustomTabulated(void* ct, const int level){ return ((CustomTabulated*) ct)->getIExact(level); }
int tsgGetQExactCustomTabulated(void* ct, const int level){ return ((CustomTabulated*) ct)->getQExact(level); }
const char* tsgGetDescriptionCustomTabulated(void* ct) { return ((CustomTabulated*) ct)->getDescription(); }

void tsgGetWeightsNodesStaticCustomTabulated(void* ct, int level, double* w, double* x) {((CustomTabulated*) ct)->getWeightsNodes(level, w, x);
}

// Note: cnodes and cweights are passed as 1D arrays, but represent a list of vectors.
void* tsgMakeCustomTabulatedFromData(const int cnum_levels, const int* cnum_nodes, const int* cprecision, const double* cnodes,
                                     const double* cweights, char* cdescription) {
    std::vector<std::vector<double>> vec_cnodes(cnum_levels), vec_cweights(cnum_levels);
    int ptr_idx = 0;
    for (int l=0; l<cnum_levels; l++) {
        vec_cnodes[l] = std::vector<double>(&cnodes[ptr_idx], &cnodes[ptr_idx] + cnum_nodes[l]);
        vec_cweights[l] = std::vector<double>(&cweights[ptr_idx], &cweights[ptr_idx] + cnum_nodes[l]);
        ptr_idx += cnum_nodes[l];
    }
    return new TasGrid::CustomTabulated(cnum_levels, std::vector<int>(cnum_nodes, cnum_nodes + cnum_levels),
                                        std::vector<int>(cprecision, cprecision + cnum_levels), std::move(vec_cnodes),
                                        std::move(vec_cweights), cdescription);
}

}
}

#endif
