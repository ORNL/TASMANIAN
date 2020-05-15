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

#ifndef __TASMANIAN_SPARSE_GRID_H
#define __TASMANIAN_SPARSE_GRID_H

// ------------ C Interface for TasmanianSparseGrid -------------- //
// NOTE: you need to explicitly call the constructor and destructor
//       in all cases, void *grid is a pointer to a C++ class

void* tsgConstructTasmanianSparseGrid();
void tsgDestructTasmanianSparseGrid(void *grid);
void tsgCopyGrid(void *destination, void *source);
void tsgCopySubGrid(void *destination, void *source, int outputs_begin, int outputs_end);
const char* tsgGetVersion();
const char* tsgGetLicense();
int tsgGetVersionMajor();
int tsgGetVersionMinor();
int tsgIsOpenMPEnabled();
void tsgWrite(void *grid, const char* filename);
void tsgWriteBinary(void *grid, const char* filename);
int tsgRead(void *grid, const char* filename);
void tsgMakeGlobalGrid(void *grid, int dimensions, int outputs, int depth, const char * sType, const char *sRule, const int *anisotropic_weights, double alpha, double beta, const char* custom_filename, const int *limit_levels);
void tsgMakeSequenceGrid(void *grid, int dimensions, int outputs, int depth, const char *sType, const char *sRule, const int *anisotropic_weights, const int *limit_levels);
void tsgMakeLocalPolynomialGrid(void *grid, int dimensions, int outputs, int depth, int order, const char *sRule, const int *limit_levels);
void tsgMakeWaveletGrid(void *grid, int dimensions, int outputs, int depth, int order, const int *limit_levels);
void tsgMakeFourierGrid(void *grid, int dimensions, int outputs, int depth, const char *sType, const int *anisotropic_weights, const int *limit_levels);
void tsgUpdateGlobalGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels);
void tsgUpdateSequenceGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels);
void tsgUpdateFourierGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels);
double tsgGetAlpha(void *grid);
double tsgGetBeta(void *grid);
int tsgGetOrder(void *grid);
int tsgGetNumDimensions(void *grid);
int tsgGetNumOutputs(void *grid);
char* tsgGetRule(void *grid);
void tsgCopyRuleChars(void *grid, int buffer_size, char *name, int *num_actual);
const char* tsgGetCustomRuleDescription(void *grid);
int tsgGetNumLoaded(void *grid);
int tsgGetNumNeeded(void *grid);
int tsgGetNumPoints(void *grid);
void tsgGetLoadedPointsStatic(void *grid, double *x);
double* tsgGetLoadedPoints(void *grid);
void tsgGetNeededPointsStatic(void *grid, double *x);
double* tsgGetNeededPoints(void *grid);
void tsgGetPointsStatic(void *grid, double *x);
double* tsgGetPoints(void *grid);
void tsgGetQuadratureWeightsStatic(void *grid, double *weights);
double* tsgGetQuadratureWeights(void *grid);
void tsgGetInterpolationWeightsStatic(void *grid, const double *x, double *weights);
double* tsgGetInterpolationWeights(void *grid, const double *x);
void tsgLoadNeededPoints(void *grid, const double *vals);
const double* tsgGetLoadedValues(void *grid);
void tsgGetLoadedValuesStatic(void *grid, double *values);
void tsgEvaluate(void *grid, const double *x, double *y);
void tsgEvaluateFast(void *grid, const double *x, double *y);
void tsgIntegrate(void *grid, double *q);
void tsgEvaluateBatch(void *grid, const double *x, int num_x, double *y);
void tsgBatchGetInterpolationWeightsStatic(void *grid, const double *x, int num_x, double *weights);
double* tsgBatchGetInterpolationWeights(void *grid, const double *x, int num_x);
int tsgIsGlobal(void *grid);
int tsgIsSequence(void *grid);
int tsgIsLocalPolynomial(void *grid);
int tsgIsWavelet(void *grid);
int tsgIsFourier(void *grid);
void tsgSetDomainTransform(void *grid, const double a[], const double b[]);
int tsgIsSetDomainTransfrom(void *grid);
void tsgClearDomainTransform(void *grid);
void tsgGetDomainTransform(void *grid, double a[], double b[]);
void tsgSetConformalTransformASIN(void *grid, const int truncation[]);
int tsgIsSetConformalTransformASIN(void *grid);
void tsgClearConformalTransform(void *grid);
void tsgGetConformalTransformASIN(void *grid, int truncation[]);
void tsgClearLevelLimits(void *grid);
void tsgGetLevelLimits(void *grid, int *limits);
void tsgSetAnisotropicRefinement(void *grid, const char * sType, int min_growth, int output, const int *level_limits);
int* tsgEstimateAnisotropicCoefficients(void *grid, const char * sType, int output, int *num_coefficients);
void tsgEstimateAnisotropicCoefficientsStatic(void *grid, const char * sType, int output, int *coefficients);
void tsgSetGlobalSurplusRefinement(void *grid, double tolerance, int output, const int *level_limits);
void tsgSetLocalSurplusRefinement(void *grid, double tolerance, const char * sRefinementType, int output, const int *level_limits, const double *scale_correction);
void tsgClearRefinement(void *grid);
void tsgMergeRefinement(void *grid);
void tsgBeginConstruction(void *grid);
int tsgIsUsingConstruction(void *grid);
void* tsgGetCandidateConstructionPointsVoidPntr(void *grid, const char *sType, int output, const int *anisotropic_weights, const int *limit_levels);
void* tsgGetCandidateConstructionPointsSurplusVoidPntr(void *grid, double tolerance, const char *sRefType, int output, const int *limit_levels, const double *scale_correction);
void tsgGetCandidateConstructionPoints(void *grid, const char *sType, int output, const int *anisotropic_weights, const int *limit_levels, int *num_points, double **x);
void tsgGetCandidateConstructionSurplusPoints(void *grid, double tolerance, const char *sRefType, int output, const int *limit_levels, const double *scale_correction, int *num_points, double **x);
int tsgGetCandidateConstructionPointsPythonGetNP(void *grid, const void *vecx);
void tsgGetCandidateConstructionPointsPythonStatic(const void *vecx, double *x);
void tsgGetCandidateConstructionPointsPythonDeleteVect(void *vecx);
void tsgLoadConstructedPoint(void *grid, const double *x, int numx, const double *y);
void tsgFinishConstruction(void *grid);
void tsgRemovePointsByHierarchicalCoefficient(void *grid, double tolerance, int output, const double *scale_correction);
void tsgEvaluateHierarchicalFunctions(void *grid, const double *x, int num_x, double *y);
void tsgEvaluateSparseHierarchicalFunctions(void *grid, const double x[], int num_x, int **pntr, int **indx, double **vals);
int tsgEvaluateSparseHierarchicalFunctionsGetNZ(void *grid, const double x[], int num_x);
void tsgEvaluateSparseHierarchicalFunctionsStatic(void *grid, const double x[], int num_x, int *pntr, int *indx, double *vals);
void tsgGetHierarchicalSupportStatic(void *grid, double support[]);
const double* tsgGetHierarchicalCoefficients(void *grid);
void tsgGetHierarchicalCoefficientsStatic(void *grid, double *coeff);
void tsgSetHierarchicalCoefficients(void *grid, const double *c);
double* tsgIntegrateHierarchicalFunctions(void *grid);
void tsgIntegrateHierarchicalFunctionsStatic(void *grid, double *integrals);
// to be called from Python only, must later call delete[] on the pointer
int* tsgPythonGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes);
// to be used in C, creates a C pointer (requires internal copy of data)
void tsgGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes, int **indexes);
void tsgPrintStats(void *grid);
void tsgEnableAcceleration(void *grid, const char *accel);
void tsgEnableAccelerationGPU(void *grid, const char *accel, int gpu);
//int tsgGetAccelerationTypeInt(void *grid){ return AccelerationMeta::getIOAccelerationInt(((TasmanianSparseGrid*) grid)->getAccelerationType()); } // int to acceleration type
const char* tsgGetAccelerationType(void *grid);
void tsgSetGPUID(void *grid, int gpuID);
int tsgGetGPUID(void *grid);
int tsgGetNumGPUs();
int tsgGetGPUMemory(int gpu);
int tsgIsAccelerationAvailable(const char *accel);
void tsgGetGPUName(int gpu, int num_buffer, char *buffer, int *num_actual);
void tsgDeleteInts(int *p);

#endif
