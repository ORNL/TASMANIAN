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

#ifdef TSG_DLL
#define TSG_API __declspec(dllexport)
#else
#ifdef TSG_DYNAMIC
#define TSG_API __declspec(dllimport)
#else
#define TSG_API
#endif
#endif

TSG_API void* tsgConstructTasmanianSparseGrid();
TSG_API void tsgDestructTasmanianSparseGrid(void *grid);
TSG_API void tsgCopyGrid(void *destination, void *source);
TSG_API const char* tsgGetVersion();
TSG_API const char* tsgGetLicense();
TSG_API int tsgGetVersionMajor();
TSG_API int tsgGetVersionMinor();
TSG_API int tsgIsOpenMPEnabled();
TSG_API void tsgErrorLogCerr(void *grid);
TSG_API void tsgDisableErrorLog(void *grid);
TSG_API void tsgWrite(void *grid, const char* filename);
TSG_API void tsgWriteBinary(void *grid, const char* filename);
TSG_API int tsgRead(void *grid, const char* filename);
TSG_API void tsgMakeGlobalGrid(void *grid, int dimensions, int outputs, int depth, const char * sType, const char *sRule, const int *anisotropic_weights, double alpha, double beta, const char* custom_filename, const int *limit_levels);
TSG_API void tsgMakeSequenceGrid(void *grid, int dimensions, int outputs, int depth, const char *sType, const char *sRule, const int *anisotropic_weights, const int *limit_levels);
TSG_API void tsgMakeLocalPolynomialGrid(void *grid, int dimensions, int outputs, int depth, int order, const char *sRule, const int *limit_levels);
TSG_API void tsgMakeWaveletGrid(void *grid, int dimensions, int outputs, int depth, int order, const int *limit_levels);
TSG_API void tsgUpdateGlobalGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels);
TSG_API void tsgUpdateSequenceGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels);
TSG_API double tsgGetAlpha(void *grid);
TSG_API double tsgGetBeta(void *grid);
TSG_API int tsgGetOrder(void *grid);
TSG_API int tsgGetNumDimensions(void *grid);
TSG_API int tsgGetNumOutputs(void *grid);
TSG_API const char* tsgGetRule(void *grid);
TSG_API const char* tsgGetCustomRuleDescription(void *grid);
TSG_API int tsgGetNumLoaded(void *grid);
TSG_API int tsgGetNumNeeded(void *grid);
TSG_API int tsgGetNumPoints(void *grid);
TSG_API void tsgGetLoadedPointsStatic(void *grid, double *x);
TSG_API double* tsgGetLoadedPoints(void *grid);
TSG_API void tsgGetNeededPointsStatic(void *grid, double *x);
TSG_API double* tsgGetNeededPoints(void *grid);
TSG_API void tsgGetPointsStatic(void *grid, double *x);
TSG_API double* tsgGetPoints(void *grid);
TSG_API void tsgGetQuadratureWeightsStatic(void *grid, double *weights);
TSG_API double* tsgGetQuadratureWeights(void *grid);
TSG_API void tsgGetInterpolationWeightsStatic(void *grid, const double *x, double *weights);
TSG_API double* tsgGetInterpolationWeights(void *grid, const double *x);
TSG_API void tsgLoadNeededPoints(void *grid, const double *vals);
TSG_API void tsgEvaluate(void *grid, const double *x, double *y);
TSG_API void tsgEvaluateFast(void *grid, const double *x, double *y);
TSG_API void tsgIntegrate(void *grid, double *q);
TSG_API void tsgEvaluateBatch(void *grid, const double *x, int num_x, double *y);
TSG_API void tsgBatchGetInterpolationWeightsStatic(void *grid, const double *x, int num_x, double *weights);
TSG_API double* tsgBatchGetInterpolationWeights(void *grid, const double *x, int num_x);
TSG_API int tsgIsGlobal(void *grid);
TSG_API int tsgIsSequence(void *grid);
TSG_API int tsgIsLocalPolynomial(void *grid);
TSG_API int tsgIsWavelet(void *grid);
TSG_API void tsgSetDomainTransform(void *grid, const double a[], const double b[]);
TSG_API int tsgIsSetDomainTransfrom(void *grid);
TSG_API void tsgClearDomainTransform(void *grid);
TSG_API void tsgGetDomainTransform(void *grid, double a[], double b[]);
TSG_API void tsgSetConformalTransformASIN(void *grid, const int truncation[]);
TSG_API int tsgIsSetConformalTransformASIN(void *grid);
TSG_API void tsgClearConformalTransform(void *grid);
TSG_API void tsgGetConformalTransformASIN(void *grid, int truncation[]);
TSG_API void tsgClearLevelLimits(void *grid);
TSG_API void tsgGetLevelLimits(void *grid, int *limits);
TSG_API void tsgSetAnisotropicRefinement(void *grid, const char * sType, int min_growth, int output, const int *level_limits);
TSG_API int* tsgEstimateAnisotropicCoefficients(void *grid, const char * sType, int output, int *num_coefficients);
TSG_API void tsgEstimateAnisotropicCoefficientsStatic(void *grid, const char * sType, int output, int *coefficients);
TSG_API void tsgSetGlobalSurplusRefinement(void *grid, double tolerance, int output, const int *level_limits);
TSG_API void tsgSetLocalSurplusRefinement(void *grid, double tolerance, const char * sRefinementType, int output, const int *level_limits);
TSG_API void tsgClearRefinement(void *grid);
TSG_API void tsgMergeRefinement(void *grid);
TSG_API void tsgRemovePointsByHierarchicalCoefficient(void *grid, double tolerance, int output, const double *scale_correction);
TSG_API void tsgEvaluateHierarchicalFunctions(void *grid, const double *x, int num_x, double *y);
TSG_API void tsgEvaluateSparseHierarchicalFunctions(void *grid, const double x[], int num_x, int **pntr, int **indx, double **vals);
TSG_API int tsgEvaluateSparseHierarchicalFunctionsGetNZ(void *grid, const double x[], int num_x);
TSG_API void tsgEvaluateSparseHierarchicalFunctionsStatic(void *grid, const double x[], int num_x, int *pntr, int *indx, double *vals);
TSG_API const double* tsgGetHierarchicalCoefficients(void *grid);
TSG_API void tsgGetHierarchicalCoefficientsStatic(void *grid, double *coeff);
TSG_API void tsgSetHierarchicalCoefficients(void *grid, const double *c);
// to be called from Python only, must later call delete[] on the pointer
TSG_API int* tsgPythonGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes);
// to be used in C, creates a C pointer (requires internal copy of data)
TSG_API void tsgGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes, int **indexes);
TSG_API void tsgPrintStats(void *grid);
TSG_API void tsgEnableAcceleration(void *grid, const char *accel);
//int tsgGetAccelerationTypeInt(void *grid){ return AccelerationMeta::getIOAccelerationInt(((TasmanianSparseGrid*) grid)->getAccelerationType()); } // int to acceleration type
TSG_API const char* tsgGetAccelerationType(void *grid);
TSG_API void tsgSetGPUID(void *grid, int gpuID);
TSG_API int tsgGetGPUID(void *grid);
TSG_API int tsgGetNumGPUs();
TSG_API int tsgGetGPUMemory(int gpu);
TSG_API int tsgIsAccelerationAvailable(const char *accel);
TSG_API void tsgGetGPUName(int gpu, int num_buffer, char *buffer, int *num_actual);
TSG_API void tsgDeleteInts(int *p);

#endif
