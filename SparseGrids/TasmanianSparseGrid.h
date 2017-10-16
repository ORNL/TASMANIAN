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
//       in all cases, void * grid is a pointer to a C++ class
void* tsgConstructTasmanianSparseGrid();
void tsgDestructTasmanianSparseGrid(void * grid);

void tsgCopyGrid(void * destination, void * source);

const char* tsgGetVersion();
const char* tsgGetLicense();
int tsgGetVersionMajor();
int tsgGetVersionMinor();
int tsgIsCudaEnabled();
int tsgIsBLASEnabled();
int tsgIsOpenMPEnabled();

void tsgErrorLogCerr(void *grid);
void tsgDisableErrorLog(void *grid);

void tsgWrite(void * grid, const char* filename);
int tsgRead(void * grid, const char* filename);
void tsgWriteBinary(void * grid, const char* filename);
int tsgReadBinary(void * grid, const char* filename);

void tsgMakeGlobalGrid(void * grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights, double alpha, double beta, const char* custom_filename);
void tsgMakeSequenceGrid(void * grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights);
void tsgMakeLocalPolynomialGrid(void * grid, int dimensions, int outputs, int depth, int order, const char * sRule);
void tsgMakeWaveletGrid(void * grid, int dimensions, int outputs, int depth, int order);

void tsgUpdateGlobalGrid(void * grid, int depth, const char * sType, const int *anisotropic_weights);
void tsgUpdateSequenceGrid(void * grid, int depth, const char * sType, const int *anisotropic_weights);

double tsgGetAlpha(void * grid);
double tsgGetBeta(void * grid);
int tsgGetOrder(void * grid);
int tsgGetNumDimensions(void * grid);
int tsgGetNumOutputs(void * grid);
const char * tsgGetRule(void * grid);
const char * tsgGetCustomRuleDescription(void * grid);

int tsgGetNumLoaded(void * grid);
int tsgGetNumNeeded(void * grid);
int tsgGetNumPoints(void * grid);

double* tsgGetLoadedPoints(void * grid);
void tsgGetLoadedPointsStatic(void * grid, double *x);
double* tsgGetNeededPoints(void * grid);
void tsgGetNeededPointsStatic(void * grid, double *x);
double* tsgGetPoints(void * grid);
void tsgGetPointsStatic(void * grid, double *x);

double* tsgGetQuadratureWeights(void * grid);
void tsgGetQuadratureWeights(void * grid, double *weights);
double* tsgGetInterpolationWeights(void * grid, const double *x);
void tsgGetInterpolationWeights(void * grid, const double *x, double *weights);

void tsgLoadNeededPoints(void * grid, const double *vals);

void tsgEvaluate(void * grid, const double *x, double *y);
void tsgEvaluateFast(void *grid, const double *x, double *y);
void tsgIntegrate(void * grid, double *q);
void tsgEvaluateBatch(void * grid, const double *x, int num_x, double *y);
double* tsgBatchGetInterpolationWeights(void * grid, const double *x, int num_x);

int tsgIsGlobal(void * grid);
int tsgIsSequence(void * grid);
int tsgIsLocalPolynomial(void * grid);
int tsgIsWavelet(void * grid);

void tsgSetDomainTransform(void * grid, const double a[], const double b[]);
int tsgIsSetDomainTransfrom(void * grid);
void tsgClearDomainTransform(void * grid);
void tsgGetDomainTransform(void * grid, double a[], double b[]);

void tsgSetConformalTransformASIN(void * grid, const int truncation[]);
int tsgIsSetConformalTransformASIN(void * grid);
void tsgClearConformalTransform(void * grid);
void tsgGetConformalTransformASIN(void * grid, int truncation[]);

void tsgSetAnisotropicRefinement(void * grid, const char * sType, int min_growth, int output);
int* tsgEstimateAnisotropicCoefficients(void * grid, const char * sType, int output, int *num_coefficients);
void tsgSetGlobalSurplusRefinement(void * grid, double tolerance, int output);
void tsgSetLocalSurplusRefinement(void * grid, double tolerance, const char * sRefinementType, int output);
void tsgClearRefinement(void * grid);

double* tsgEvalHierarchicalFunctions(void * grid, const double *x);
double* tsgBatchEvalHierarchicalFunctions(void * grid, const double *x, int num_x);
void tsgSetHierarchicalCoefficients(void * grid, const double *c);
const double* tsgGetSurpluses(void *grid);

int* tsgGetGlobalPolynomialSpace(void * grid, int interpolation, int *num_indexes);

void tsgPrintStats(void *grid);

void tsgEnableAcceleration(void *grid, const char *accel);
int tsgGetAccelerationType(void *grid);

void tsgSetGPUID(void *grid, int gpuID);
int tsgGetGPUID(void *grid);
int tsgGetNumGPUs();
int tsgGetGPUmemory(int gpu);
const char* tsgGetGPUname(int gpu);

#endif
