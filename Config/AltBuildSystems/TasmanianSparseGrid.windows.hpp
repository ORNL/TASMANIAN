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

#ifndef __TASMANIAN_SPARSE_GRID_HPP
#define __TASMANIAN_SPARSE_GRID_HPP

#ifdef TSG_DLL
#define TSG_API __declspec(dllexport)
#else
#ifdef TSG_DYNAMIC
#define TSG_API __declspec(dllimport)
#else
#define TSG_API
#endif
#endif

#include "tasmanianConfig.hpp"

#include "tsgEnumerates.hpp"

#include "tsgGridGlobal.hpp"
#include "tsgGridSequence.hpp"
#include "tsgGridLocalPolynomial.hpp"
#include "tsgGridWavelet.hpp"

#include <iomanip> // only needed for printStats()

namespace TasGrid{

class TSG_API TasmanianSparseGrid{
public:
    TasmanianSparseGrid();
    TasmanianSparseGrid(const TasmanianSparseGrid &source);
    ~TasmanianSparseGrid();

    static const char* getVersion(); // human readable
    static int getVersionMajor();
    static int getVersionMinor();
    static const char* getLicense(); // human readable
    static bool isOpenMPEnabled();

    void setErrorLog(std::ostream *os);
    void disableLog();

    void write(const char *filename, bool binary = false) const;
    bool read(const char *filename);

    void write(std::ofstream &ofs, bool binary = false) const;
    bool read(std::ifstream &ifs, bool binary = false);

    void makeGlobalGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights = 0, double alpha = 0.0, double beta = 0.0, const char* custom_filename = 0, const int *level_limits = 0);
    void makeSequenceGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights = 0, const int *level_limits = 0);
    void makeLocalPolynomialGrid(int dimensions, int outputs, int depth, int order = 1, TypeOneDRule rule = rule_localp, const int *level_limits = 0);
    void makeWaveletGrid(int dimensions, int outputs, int depth, int order = 1, const int *level_limits = 0);
    void copyGrid(const TasmanianSparseGrid *source);

    void updateGlobalGrid(int depth, TypeDepth type, const int *anisotropic_weights = 0, const int *level_limits = 0);
    void updateSequenceGrid(int depth, TypeDepth type, const int *anisotropic_weights = 0, const int *level_limits = 0);

    double getAlpha() const;
    double getBeta() const;
    int getOrder() const;

    int getNumDimensions() const;
    int getNumOutputs() const;
    TypeOneDRule getRule() const;
    const char* getCustomRuleDescription() const; // used only for Global Grids with rule_customtabulated

    int getNumLoaded() const;
    int getNumNeeded() const;
    int getNumPoints() const; // returns the number of loaded points unless no points are loaded, then returns the number of needed points

    double* getLoadedPoints() const;
    double* getNeededPoints() const;
    double* getPoints() const; // returns the loaded points unless no points are loaded, then returns the needed points
    void getLoadedPoints(double *x) const;
    void getNeededPoints(double *x) const;
    void getPoints(double *x) const; // returns the loaded points unless no points are loaded, then returns the needed points

    double* getQuadratureWeights() const;
    double* getInterpolationWeights(const double x[]) const;
    void getQuadratureWeights(double weights[]) const;
    void getInterpolationWeights(const double x[], double weights[]) const;

    void loadNeededPoints(const double *vals);

    void evaluate(const double x[], double y[]) const;
    void evaluateFast(const double x[], double y[]) const; // evaluate that is potentially not thread safe!
    void evaluateBatch(const double x[], int num_x, double y[]) const; // uses acceleration, OpenMP, BLAS, GPU, etc.
    void integrate(double q[]) const;

    bool isGlobal() const;
    bool isSequence() const;
    bool isLocalPolynomial() const;
    bool isWavelet() const;

    void setDomainTransform(const double a[], const double b[]); // set the ranges of the box
    bool isSetDomainTransfrom() const;
    void clearDomainTransform();
    void getDomainTransform(double a[], double b[]) const;

    void setConformalTransformASIN(const int truncation[]);
    bool isSetConformalTransformASIN() const;
    void clearConformalTransform();
    void getConformalTransformASIN(int truncation[]) const;

    void clearLevelLimits(); // level limits will be set anew if non-null vector is given to refine command
    void getLevelLimits(int *limits) const; // static, assume limits is already allocated with length getNumDimensions()

    void setAnisotropicRefinement(TypeDepth type, int min_growth, int output, const int *level_limits = 0);
    int* estimateAnisotropicCoefficients(TypeDepth type, int output);
    void setSurplusRefinement(double tolerance, int output, const int *level_limits = 0);
    void setSurplusRefinement(double tolerance, TypeRefinement criteria, int output = -1, const int *level_limits = 0, const double *scale_correction = 0); // -1 indicates using all outputs
    void clearRefinement();
    void mergeRefinement();

    const double* getHierarchicalCoefficients() const; // formerly getSurpluses();
    void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const;
    void evaluateSparseHierarchicalFunctions(const double x[], int num_x, int* &pntr, int* &indx, double* &vals) const;
    void setHierarchicalCoefficients(const double c[]);

    void getGlobalPolynomialSpace(bool interpolation, int &num_indexes, int* &poly) const;

    void printStats() const;

    void printStatsLog() const;

    void enableAcceleration(TypeAcceleration acc);
    void forceSparseAlgorithmForLocalPolynomials();
    TypeAcceleration getAccelerationType() const;
    static bool isAccelerationAvailable(TypeAcceleration acc);

    // CUDA management functions
    void setGPUID(int new_gpuID);
    int getGPUID() const;
    static int getNumGPUs();
    static int getGPUMemory(int gpu); // returns the MB of a given GPU
    static char* getGPUName(int gpu); // returns a null-terminated char array

    // Do not use these functions within C++, not as efficient as evaluateSparseHierarchicalFunctions()
    // these functions are needed for interfaces with other languages
    int evaluateSparseHierarchicalFunctionsGetNZ(const double x[], int num_x) const;
    void evaluateSparseHierarchicalFunctionsStatic(const double x[], int num_x, int pntr[], int indx[], double vals[]) const;


    // EXPERIMENTAL: works only for LocalPolynomial grids (will be moved to sequences and others)
    void removePointsByHierarchicalCoefficient(double tolerance, int output = -1, const double *scale_correction = 0); // have python error tests, but needs math consistency test

    void evaluateHierarchicalFunctionsGPU(const double gpu_x[], int cpu_num_x, double gpu_y[]) const; // EXPERIMENTAL, works only for localPolynomial with order 1 or 2
    void evaluateSparseHierarchicalFunctionsGPU(const double gpu_x[], int cpu_num_x, int* &gpu_pntr, int* &gpu_indx, double* &gpu_vals, int &num_nz) const;

    // TODO
    // int* getIndexOfRefinedPoints(); // tells you the index of the refined points in the big vector after mergeRefinement (or loadNeededPoints())

    // WARNING: the functions below are mostly for debugging and research purposes
    // do not modify the returned pointers, unless you really know what you are doing
    const int* getPointsIndexes() const;
    const int* getNeededIndexes() const;

    // Expose the internal pointers, for testing and debugging only, must be suppressed for an official release
    //inline GridGlobal*          getGridGlobalPointer(){ return global; }
    //inline GridSequence*        getGridSequencePointer(){ return sequence; }
    //inline GridLocalPolynomial* getGridLocalPolyPointer(){ return pwpoly; }

protected:
    void clear();

    void printGridStats(std::ostream *os) const;

    void mapCanonicalToTransformed(int num_dimensions, int num_points, TypeOneDRule rule, double x[]) const;
    void mapTransformedToCanonical(int num_dimensions, int num_points, TypeOneDRule rule, double x[]) const;
    double getQuadratureScale(int num_dimensions, TypeOneDRule rule) const;

    void mapConformalCanonicalToTransformed(int num_dimensions, int num_points, double x[]) const;
    void mapConformalTransformedToCanonical(int num_dimensions, int num_points, double x[]) const;
    void mapConformalWeights(int num_dimensions, int num_points, double weights[]) const;

    const double* formCanonicalPoints(const double *x, double* &x_temp, int num_x) const;
    const double* formCanonicalPointsGPU(const double *gpu_x, double* &gpu_x_temp, int num_x) const;
    void clearCanonicalPoints(double* &x_temp) const;
    void clearCanonicalPointsGPU(double* &x_temp) const;
    void formTransformedPoints(int num_points, double x[]) const; // when calling get***Points()

    void writeAscii(std::ofstream &ofs) const;
    bool readAscii(std::ifstream &ifs);

    void writeBinary(std::ofstream &ofs) const;
    bool readBinary(std::ifstream &ifs);

private:
    BaseCanonicalGrid *base;

    GridGlobal *global;
    GridSequence *sequence;
    GridLocalPolynomial *pwpoly;
    GridWavelet *wavelet;

    double *domain_transform_a, *domain_transform_b;
    int *conformal_asin_power;
    int *llimits;

    TypeAcceleration acceleration;
    int gpuID;
    mutable AccelerationDomainTransform *acc_domain;

    std::ostream *logstream;
};

extern "C" TSG_API void* tsgConstructTasmanianSparseGrid();
extern "C" TSG_API void tsgDestructTasmanianSparseGrid(void *grid);
extern "C" TSG_API void tsgCopyGrid(void *destination, void *source);
extern "C" TSG_API const char* tsgGetVersion();
extern "C" TSG_API const char* tsgGetLicense();
extern "C" TSG_API int tsgGetVersionMajor();
extern "C" TSG_API int tsgGetVersionMinor();
extern "C" TSG_API int tsgIsOpenMPEnabled();
extern "C" TSG_API void tsgErrorLogCerr(void *grid);
extern "C" TSG_API void tsgDisableErrorLog(void *grid);
extern "C" TSG_API void tsgWrite(void *grid, const char* filename);
extern "C" TSG_API void tsgWriteBinary(void *grid, const char* filename);
extern "C" TSG_API int tsgRead(void *grid, const char* filename);
extern "C" TSG_API void tsgMakeGlobalGrid(void *grid, int dimensions, int outputs, int depth, const char * sType, const char *sRule, const int *anisotropic_weights, double alpha, double beta, const char* custom_filename, const int *limit_levels);
extern "C" TSG_API void tsgMakeSequenceGrid(void *grid, int dimensions, int outputs, int depth, const char *sType, const char *sRule, const int *anisotropic_weights, const int *limit_levels);
extern "C" TSG_API void tsgMakeLocalPolynomialGrid(void *grid, int dimensions, int outputs, int depth, int order, const char *sRule, const int *limit_levels);
extern "C" TSG_API void tsgMakeWaveletGrid(void *grid, int dimensions, int outputs, int depth, int order, const int *limit_levels);
extern "C" TSG_API void tsgUpdateGlobalGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels);
extern "C" TSG_API void tsgUpdateSequenceGrid(void *grid, int depth, const char * sType, const int *anisotropic_weights, const int *limit_levels);
extern "C" TSG_API double tsgGetAlpha(void *grid);
extern "C" TSG_API double tsgGetBeta(void *grid);
extern "C" TSG_API int tsgGetOrder(void *grid);
extern "C" TSG_API int tsgGetNumDimensions(void *grid);
extern "C" TSG_API int tsgGetNumOutputs(void *grid);
extern "C" TSG_API const char* tsgGetRule(void *grid);
extern "C" TSG_API const char* tsgGetCustomRuleDescription(void *grid);
extern "C" TSG_API int tsgGetNumLoaded(void *grid);
extern "C" TSG_API int tsgGetNumNeeded(void *grid);
extern "C" TSG_API int tsgGetNumPoints(void *grid);
extern "C" TSG_API void tsgGetLoadedPointsStatic(void *grid, double *x);
extern "C" TSG_API double* tsgGetLoadedPoints(void *grid);
extern "C" TSG_API void tsgGetNeededPointsStatic(void *grid, double *x);
extern "C" TSG_API double* tsgGetNeededPoints(void *grid);
extern "C" TSG_API void tsgGetPointsStatic(void *grid, double *x);
extern "C" TSG_API double* tsgGetPoints(void *grid);
extern "C" TSG_API void tsgGetQuadratureWeightsStatic(void *grid, double *weights);
extern "C" TSG_API double* tsgGetQuadratureWeights(void *grid);
extern "C" TSG_API void tsgGetInterpolationWeightsStatic(void *grid, const double *x, double *weights);
extern "C" TSG_API double* tsgGetInterpolationWeights(void *grid, const double *x);
extern "C" TSG_API void tsgLoadNeededPoints(void *grid, const double *vals);
extern "C" TSG_API void tsgEvaluate(void *grid, const double *x, double *y);
extern "C" TSG_API void tsgEvaluateFast(void *grid, const double *x, double *y);
extern "C" TSG_API void tsgIntegrate(void *grid, double *q);
extern "C" TSG_API void tsgEvaluateBatch(void *grid, const double *x, int num_x, double *y);
extern "C" TSG_API void tsgBatchGetInterpolationWeightsStatic(void *grid, const double *x, int num_x, double *weights);
extern "C" TSG_API double* tsgBatchGetInterpolationWeights(void *grid, const double *x, int num_x);
extern "C" TSG_API int tsgIsGlobal(void *grid);
extern "C" TSG_API int tsgIsSequence(void *grid);
extern "C" TSG_API int tsgIsLocalPolynomial(void *grid);
extern "C" TSG_API int tsgIsWavelet(void *grid);
extern "C" TSG_API void tsgSetDomainTransform(void *grid, const double a[], const double b[]);
extern "C" TSG_API int tsgIsSetDomainTransfrom(void *grid);
extern "C" TSG_API void tsgClearDomainTransform(void *grid);
extern "C" TSG_API void tsgGetDomainTransform(void *grid, double a[], double b[]);
extern "C" TSG_API void tsgSetConformalTransformASIN(void *grid, const int truncation[]);
extern "C" TSG_API int tsgIsSetConformalTransformASIN(void *grid);
extern "C" TSG_API void tsgClearConformalTransform(void *grid);
extern "C" TSG_API void tsgGetConformalTransformASIN(void *grid, int truncation[]);
extern "C" TSG_API void tsgClearLevelLimits(void *grid);
extern "C" TSG_API void tsgGetLevelLimits(void *grid, int *limits);
extern "C" TSG_API void tsgSetAnisotropicRefinement(void *grid, const char * sType, int min_growth, int output, const int *level_limits);
extern "C" TSG_API int* tsgEstimateAnisotropicCoefficients(void *grid, const char * sType, int output, int *num_coefficients);
extern "C" TSG_API void tsgEstimateAnisotropicCoefficientsStatic(void *grid, const char * sType, int output, int *coefficients);
extern "C" TSG_API void tsgSetGlobalSurplusRefinement(void *grid, double tolerance, int output, const int *level_limits);
extern "C" TSG_API void tsgSetLocalSurplusRefinement(void *grid, double tolerance, const char * sRefinementType, int output, const int *level_limits);
extern "C" TSG_API void tsgClearRefinement(void *grid);
extern "C" TSG_API void tsgMergeRefinement(void *grid);
extern "C" TSG_API void tsgRemovePointsByHierarchicalCoefficient(void *grid, double tolerance, int output, const double *scale_correction);
extern "C" TSG_API void tsgEvaluateHierarchicalFunctions(void *grid, const double *x, int num_x, double *y);
extern "C" TSG_API void tsgEvaluateSparseHierarchicalFunctions(void *grid, const double x[], int num_x, int **pntr, int **indx, double **vals);
extern "C" TSG_API int tsgEvaluateSparseHierarchicalFunctionsGetNZ(void *grid, const double x[], int num_x);
extern "C" TSG_API void tsgEvaluateSparseHierarchicalFunctionsStatic(void *grid, const double x[], int num_x, int *pntr, int *indx, double *vals);
extern "C" TSG_API const double* tsgGetHierarchicalCoefficients(void *grid);
extern "C" TSG_API void tsgGetHierarchicalCoefficientsStatic(void *grid, double *coeff);
extern "C" TSG_API void tsgSetHierarchicalCoefficients(void *grid, const double *c);
// to be called from Python only, must later call delete[] on the pointer
extern "C" TSG_API int* tsgPythonGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes);
// to be used in C, creates a C pointer (requires internal copy of data)
extern "C" TSG_API void tsgGetGlobalPolynomialSpace(void *grid, int interpolation, int *num_indexes, int **indexes);
extern "C" TSG_API void tsgPrintStats(void *grid);
extern "C" TSG_API void tsgEnableAcceleration(void *grid, const char *accel);
//int tsgGetAccelerationTypeInt(void *grid){ return AccelerationMeta::getIOAccelerationInt(((TasmanianSparseGrid*) grid)->getAccelerationType()); } // int to acceleration type
extern "C" TSG_API const char* tsgGetAccelerationType(void *grid);
extern "C" TSG_API void tsgSetGPUID(void *grid, int gpuID);
extern "C" TSG_API int tsgGetGPUID(void *grid);
extern "C" TSG_API int tsgGetNumGPUs();
extern "C" TSG_API int tsgGetGPUMemory(int gpu);
extern "C" TSG_API int tsgIsAccelerationAvailable(const char *accel);
extern "C" TSG_API void tsgGetGPUName(int gpu, int num_buffer, char *buffer, int *num_actual);
extern "C" TSG_API void tsgDeleteInts(int *p);

}

#endif
