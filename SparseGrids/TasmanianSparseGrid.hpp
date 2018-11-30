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

#include "TasmanianConfig.hpp"

#include "tsgEnumerates.hpp"

#include "tsgGridGlobal.hpp"
#include "tsgGridSequence.hpp"
#include "tsgGridLocalPolynomial.hpp"
#include "tsgGridWavelet.hpp"
#include "tsgGridFourier.hpp"

#include <iomanip> // only needed for printStats()

namespace TasGrid{

class TasmanianSparseGrid{
public:
    TasmanianSparseGrid();
    TasmanianSparseGrid(const TasmanianSparseGrid &source);
    ~TasmanianSparseGrid();

    static const char* getVersion(); // human readable
    static int getVersionMajor();
    static int getVersionMinor();
    static const char* getLicense(); // human readable
    static const char* getGitCommitHash();
    static const char* getCmakeCxxFlags();
    static bool isOpenMPEnabled();

    void write(const char *filename, bool binary = false) const;
    void read(const char *filename); // auto-check if format is binary or ascii

    void write(std::ofstream &ofs, bool binary = false) const;
    void read(std::ifstream &ifs, bool binary = false);

    void makeGlobalGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights = 0, double alpha = 0.0, double beta = 0.0, const char* custom_filename = 0, const int *level_limits = 0);
    void makeGlobalGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const std::vector<int> &anisotropic_weights, double alpha = 0.0, double beta = 0.0, const char* custom_filename = 0, const std::vector<int> &level_limits = std::vector<int>());

    void makeSequenceGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights = 0, const int *level_limits = 0);
    void makeSequenceGrid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const std::vector<int> &anisotropic_weights, const std::vector<int> &level_limits = std::vector<int>());

    void makeLocalPolynomialGrid(int dimensions, int outputs, int depth, int order = 1, TypeOneDRule rule = rule_localp, const int *level_limits = 0);
    void makeLocalPolynomialGrid(int dimensions, int outputs, int depth, int order, TypeOneDRule rule, const std::vector<int> &level_limits);

    void makeWaveletGrid(int dimensions, int outputs, int depth, int order = 1, const int *level_limits = 0);
    void makeWaveletGrid(int dimensions, int outputs, int depth, int order, const std::vector<int> &level_limits);

    void makeFourierGrid(int dimensions, int outputs, int depth, TypeDepth type, const int* anisotropic_weights = 0, const int* level_limits = 0);
    void makeFourierGrid(int dimensions, int outputs, int depth, TypeDepth type, const std::vector<int> &anisotropic_weights, const std::vector<int> &level_limits = std::vector<int>());

    void copyGrid(const TasmanianSparseGrid *source);

    void updateGlobalGrid(int depth, TypeDepth type, const int *anisotropic_weights = 0, const int *level_limits = 0);
    void updateGlobalGrid(int depth, TypeDepth type, const std::vector<int> &anisotropic_weights, const std::vector<int> &level_limits = std::vector<int>());

    void updateSequenceGrid(int depth, TypeDepth type, const int *anisotropic_weights = 0, const int *level_limits = 0);
    void updateSequenceGrid(int depth, TypeDepth type, const std::vector<int> &anisotropic_weights, const std::vector<int> &level_limits = std::vector<int>());

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

    void getLoadedPoints(double *x) const; // using static memory, assuming x has size num_dimensions X get***Points()
    void getNeededPoints(double *x) const;
    void getPoints(double *x) const; // returns the loaded points unless no points are loaded, then returns the needed points

    void getLoadedPoints(std::vector<double> &x) const; // dynamic memory, resizes x
    void getNeededPoints(std::vector<double> &x) const;
    void getPoints(std::vector<double> &x) const; // returns the loaded points unless no points are loaded, then returns the needed points

    double* getQuadratureWeights() const;
    double* getInterpolationWeights(const double x[]) const;

    void getQuadratureWeights(double weights[]) const; // static memory, assumes that weights has size getNumPoints()
    void getInterpolationWeights(const double x[], double weights[]) const; // static memory, assumes that weights has size getNumPoints()

    void getQuadratureWeights(std::vector<double> &weights) const; // dynamic memory, resizes weights
    void getInterpolationWeights(const std::vector<double> &x, std::vector<double> &weights) const; // dynamic memory, resizes weights

    void loadNeededPoints(const double *vals); // no error checking
    void loadNeededPoints(const std::vector<double> &vals); // checks if vals has size num_outputs X getNumNeeded()

    void evaluate(const double x[], double y[]) const; // has size num_dimensions, y has size num_outputs
    void evaluateFast(const double x[], double y[]) const; // evaluate that is potentially not thread safe!
    void evaluateBatch(const double x[], int num_x, double y[]) const; // uses acceleration, OpenMP, BLAS, GPU, etc., x is num_dimensions X num_x, y is num_outputs X num_x
    void integrate(double q[]) const; // y has size num_outputs

    // same as above, but num_x = x.size() / num_dimensions, and y is resized
    void evaluate(const std::vector<double> &x, std::vector<double> &y) const;
    void evaluateFast(const std::vector<double> &x, std::vector<double> &y) const;
    void evaluateBatch(const std::vector<double> &x, std::vector<double> &y) const;
    void integrate(std::vector<double> &q) const;

    bool isGlobal() const;
    bool isSequence() const;
    bool isLocalPolynomial() const;
    bool isWavelet() const;
    bool isFourier() const;
    inline bool empty() const{ return (base.get() == nullptr); }

    void setDomainTransform(const double a[], const double b[]); // set the ranges of the box, a[] and b[] must have size num_dimensions
    bool isSetDomainTransfrom() const;
    void clearDomainTransform();
    void getDomainTransform(double a[], double b[]) const;
    void setDomainTransform(const std::vector<double> &a, const std::vector<double> &b);
    void getDomainTransform(std::vector<double> &a, std::vector<double> &b) const;

    void setConformalTransformASIN(const int truncation[]);
    bool isSetConformalTransformASIN() const;
    void clearConformalTransform();
    void getConformalTransformASIN(int truncation[]) const;

    void clearLevelLimits(); // level limits will be set anew if non-null vector is given to refine command
    void getLevelLimits(int *limits) const; // static, assume limits is already allocated with length num_dimensions
    void getLevelLimits(std::vector<int> &limits) const; // allocates the vector

    void setAnisotropicRefinement(TypeDepth type, int min_growth, int output, const int *level_limits = 0);
    void setAnisotropicRefinement(TypeDepth type, int min_growth, int output, const std::vector<int> &level_limits);

    int* estimateAnisotropicCoefficients(TypeDepth type, int output);
    void estimateAnisotropicCoefficients(TypeDepth type, int output, std::vector<int> &weights);

    void setSurplusRefinement(double tolerance, int output, const int *level_limits = 0);
    void setSurplusRefinement(double tolerance, int output, const std::vector<int> &level_limits);

    void setSurplusRefinement(double tolerance, TypeRefinement criteria, int output = -1, const int *level_limits = 0, const double *scale_correction = 0); // -1 indicates using all outputs
    void setSurplusRefinement(double tolerance, TypeRefinement criteria, int output, const std::vector<int> &level_limits, const std::vector<double> &scale_correction = std::vector<double>()); // -1 indicates using all outputs
    // if using only one output, scale_correction has size getNumPoints(), if using all outputs, scale_correction has size num_outputs X getNumDimensions()
    // correction is multiplied by the relative coefficient (normalized across outputs) when comparing against the tolerance
    // empty scale_correction is equivalent to passing constant 1.0 vector

    void clearRefinement();
    void mergeRefinement(); // merges all points but resets all loaded values to 0.0

    //! \brief Begin a dynamic construction procedure, also calls \b clearRefinement() (cheap to call, only sets some flags)
    void beginConstruction();
    //! \brief Generate a sorted list of points weighted by descending importance (expensive call, equivalent to set-refinement)
    void getCandidateConstructionPoints(std::vector<double> &x);
    //! \brief Add the value of a single point (if the tensor of the point is not complete, the grid will not be updated but the value will be stored)
    void loadConstructedPoint(const std::vector<double> &x, const std::vector<double> &y);
    //! \brief End the procedure, clears flags and unused constructed points, can go back to using regular refinement
    void finishConstruction();

    const double* getHierarchicalCoefficients() const; // formerly getSurpluses();  returns an alias to internal data structure
    void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const;
    void evaluateSparseHierarchicalFunctions(const double x[], int num_x, int* &pntr, int* &indx, double* &vals) const;
    void evaluateSparseHierarchicalFunctions(const std::vector<double> &x, std::vector<int> &pntr, std::vector<int> &indx, std::vector<double> &vals) const;
    void setHierarchicalCoefficients(const double c[]);

    void evaluateHierarchicalFunctions(const std::vector<double> &x, std::vector<double> &y) const;
    void setHierarchicalCoefficients(const std::vector<double> &c);

    void getGlobalPolynomialSpace(bool interpolation, int &num_indexes, int* &poly) const;

    void printStats(std::ostream &os = std::cout) const;

    void enableAcceleration(TypeAcceleration acc);
    void favorSparseAcceleration(bool favor);
    TypeAcceleration getAccelerationType() const;
    static bool isAccelerationAvailable(TypeAcceleration acc);

    // CUDA management functions
    void setGPUID(int new_gpuID);
    int getGPUID() const;
    static int getNumGPUs();
    static int getGPUMemory(int gpu); // returns the MB of a given GPU
    static char* getGPUName(int gpu); // returns a null-terminated char array

    // functions assuming vectors are pre-allocated on the GPU (stable for the implemented cases)
    void evaluateHierarchicalFunctionsGPU(const double gpu_x[], int cpu_num_x, double gpu_y[]) const; // does not work for Global or LocalPolynomial with order > 2
    void evaluateSparseHierarchicalFunctionsGPU(const double gpu_x[], int cpu_num_x, int* &gpu_pntr, int* &gpu_indx, double* &gpu_vals, int &num_nz) const; // LocalPolynomial only with order 0, 1, or 2


    // Do not use these functions within C++, not as efficient as evaluateSparseHierarchicalFunctions()
    // these functions are needed for interfaces with other languages
    int evaluateSparseHierarchicalFunctionsGetNZ(const double x[], int num_x) const;
    void evaluateSparseHierarchicalFunctionsStatic(const double x[], int num_x, int pntr[], int indx[], double vals[]) const;


    // EXPERIMENTAL: works only for LocalPolynomial grids (will be moved to sequences and others)
    void removePointsByHierarchicalCoefficient(double tolerance, int output = -1, const double *scale_correction = 0); // have python error tests, but needs math consistency test


    // WARNING: the functions below are mostly for debugging and research purposes
    // do not modify the returned pointers, unless you really know what you are doing
    const int* getPointsIndexes() const;
    const int* getNeededIndexes() const;

    // make these protected for a release
    inline GridGlobal*          getGridGlobal(){          return (GridGlobal*)          base.get(); }
    inline GridSequence*        getGridSequence(){        return (GridSequence*)        base.get(); }
    inline GridLocalPolynomial* getGridLocalPolynomial(){ return (GridLocalPolynomial*) base.get(); }
    inline GridFourier*         getGridFourier(){         return (GridFourier*)         base.get(); }
    inline GridWavelet*         getGridWavelet(){         return (GridWavelet*)         base.get(); }
    inline const GridGlobal*          getGridGlobal() const{          return (const GridGlobal*)          base.get(); }
    inline const GridSequence*        getGridSequence() const{        return (const GridSequence*)        base.get(); }
    inline const GridLocalPolynomial* getGridLocalPolynomial() const{ return (const GridLocalPolynomial*) base.get(); }
    inline const GridFourier*         getGridFourier() const{         return (const GridFourier*)         base.get(); }
    inline const GridWavelet*         getGridWavelet() const{         return (const GridWavelet*)         base.get(); }

protected:
    void clear();

    void mapCanonicalToTransformed(int num_dimensions, int num_points, TypeOneDRule rule, double x[]) const;
    void mapTransformedToCanonical(int num_dimensions, int num_points, TypeOneDRule rule, double x[]) const;
    double getQuadratureScale(int num_dimensions, TypeOneDRule rule) const;

    void mapConformalCanonicalToTransformed(int num_dimensions, int num_points, double x[]) const;
    void mapConformalTransformedToCanonical(int num_dimensions, int num_points, Data2D<double> &x) const;
    void mapConformalWeights(int num_dimensions, int num_points, double weights[]) const;

    const double* formCanonicalPoints(const double *x, Data2D<double> &x_temp, int num_x) const;
    #ifdef Tasmanian_ENABLE_CUDA
    const double* formCanonicalPointsGPU(const double *gpu_x, int num_x, cudaDoubles &gpu_x_temp) const;
    #endif
    void formTransformedPoints(int num_points, double x[]) const; // when calling get***Points()

    void writeAscii(std::ofstream &ofs) const;
    void readAscii(std::ifstream &ifs);

    void writeBinary(std::ofstream &ofs) const;
    void readBinary(std::ifstream &ifs);

private:
    std::unique_ptr<BaseCanonicalGrid> base;

    std::vector<double> domain_transform_a, domain_transform_b;
    std::vector<int> conformal_asin_power;
    std::vector<int> llimits;

    TypeAcceleration acceleration;
    int gpuID;

    bool usingDynamicConstruction;

    #ifdef Tasmanian_ENABLE_CUDA
    mutable AccelerationDomainTransform acc_domain;
    #endif
};

}

#endif
