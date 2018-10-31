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

#ifndef __TASMANIAN_SPARSE_GRID_GLOBAL_NESTED_HPP
#define __TASMANIAN_SPARSE_GRID_GLOBAL_NESTED_HPP

#include <cstdlib>

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgIndexManipulator.hpp"
#include "tsgLinearSolvers.hpp"
#include "tsgCacheLagrange.hpp"
#include "tsgOneDimensionalWrapper.hpp"
#include "tsgGridCore.hpp"

#include "tsgAcceleratedDataStructures.hpp"

namespace TasGrid{

class GridSequence : public BaseCanonicalGrid{
public:
    GridSequence();
    GridSequence(const GridSequence &seq);
    ~GridSequence();

    bool isSequence() const{ return true; }

    void write(std::ofstream &ofs) const;
    void read(std::ifstream &ifs);

    void writeBinary(std::ofstream &ofs) const;
    void readBinary(std::ifstream &ifs);

    void makeGrid(int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, TypeOneDRule crule, const std::vector<int> &anisotropic_weights, const std::vector<int> &level_limits);
    void copyGrid(const GridSequence *seq);
    void setPoints(MultiIndexSet &pset, int cnum_outputs, TypeOneDRule crule);

    void updateGrid(int depth, TypeDepth type, const std::vector<int> &anisotropic_weights, const std::vector<int> &level_limits);
    void updateGrid(MultiIndexSet &update);

    int getNumDimensions() const;
    int getNumOutputs() const;
    TypeOneDRule getRule() const;

    int getNumLoaded() const;
    int getNumNeeded() const;
    int getNumPoints() const; // returns the number of loaded points unless no points are loaded, then returns the number of needed points

    void getLoadedPoints(double *x) const;
    void getNeededPoints(double *x) const;
    void getPoints(double *x) const; // returns the loaded points unless no points are loaded, then returns the needed points

    void getQuadratureWeights(double weights[]) const;
    void getInterpolationWeights(const double x[], double weights[]) const;

    void loadNeededPoints(const double *vals, TypeAcceleration acc = accel_none);

    void evaluate(const double x[], double y[]) const;
    void integrate(double q[], double *conformal_correction) const;

    void evaluateFastCPUblas(const double x[], double y[]) const;
    void evaluateFastGPUcublas(const double x[], double y[]) const;
    void evaluateFastGPUcuda(const double x[], double y[]) const;
    void evaluateFastGPUmagma(int gpuID, const double x[], double y[]) const;

    void evaluateBatch(const double x[], int num_x, double y[]) const;
    void evaluateBatchCPUblas(const double x[], int num_x, double y[]) const;
    void evaluateBatchGPUcublas(const double x[], int num_x, double y[]) const;
    void evaluateBatchGPUcuda(const double x[], int num_x, double y[]) const;
    void evaluateBatchGPUmagma(int gpuID, const double x[], int num_x, double y[]) const;

    void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const;
    #ifdef Tasmanian_ENABLE_CUDA
    void evaluateHierarchicalFunctionsGPU(const double x[], int num_x, double y[]) const;
    #endif

    void estimateAnisotropicCoefficients(TypeDepth type, int output, std::vector<int> &weights) const;
    void setAnisotropicRefinement(TypeDepth type, int min_growth, int output, const std::vector<int> &level_limits);
    void setSurplusRefinement(double tolerance, int output, const std::vector<int> &level_limits);
    void clearRefinement();
    void mergeRefinement();

    void setHierarchicalCoefficients(const double c[], TypeAcceleration acc);

    void getPolynomialSpace(bool interpolation, int &n, int* &poly) const;

    const double* getSurpluses() const;
    const int* getPointIndexes() const;

    void clearAccelerationData();

protected:
    void reset();

    void evalHierarchicalFunctions(const double x[], double fvalues[]) const;

    void prepareSequence();
    void cacheBasisIntegrals(std::vector<double> &integ) const;

    template<typename T>
    void cacheBasisValues(const T x[], std::vector<std::vector<T>> &cache) const{
        cache.resize(num_dimensions);
        for(int j=0; j<num_dimensions; j++){
            cache[j].resize(max_levels[j] + 1);
            T b = 1.0;
            T this_x = x[j];
            cache[j][0] = b;
            for(int i=0; i<max_levels[j]; i++){
                b *= (this_x - nodes[i]);
                cache[j][i+1] = b;
            }
            for(int i=1; i<=max_levels[j]; i++){
                cache[j][i] /= coeff[i];
            }
        }
    }

    void recomputeSurpluses();
    void applyTransformationTransposed(double weights[]) const;

    double evalBasis(const int f[], const int p[]) const; // evaluate function corresponding to f at p

    #ifdef Tasmanian_ENABLE_CUDA
    void loadCudaNodes() const{
        if (cuda_num_nodes.size() != 0) return;
        int maxl = max_levels[0];
        for(auto l : max_levels) if (maxl < l) maxl = l;
        cuda_nodes.load(nodes);
        cuda_coeffs.load(coeff);
        std::vector<int> num_nodes = max_levels;
        for(auto &n : num_nodes) n++;
        cuda_num_nodes.load(num_nodes);
        const MultiIndexSet *work = (points.empty()) ? &needed : &points;
        int num_points = work->getNumIndexes();
        Data2D<int> transpoints; transpoints.resize(work->getNumIndexes(), num_dimensions);
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++){
                transpoints.getStrip(j)[i] = work->getIndex(i)[j];
            }
        }
        cuda_points.load(*(transpoints.getVector()));
    }
    void clearCudaNodes(){
        cuda_nodes.clear();
        cuda_coeffs.clear();
        cuda_num_nodes.clear();
        cuda_points.clear();
    }
    #endif

private:
    int num_dimensions, num_outputs;
    TypeOneDRule rule;

    MultiIndexSet points;
    MultiIndexSet needed;

    std::vector<double> surpluses;
    std::vector<double> nodes;
    std::vector<double> coeff;

    StorageSet values;

    std::vector<int> max_levels;

    #ifdef Tasmanian_ENABLE_CUDA
    mutable LinearAlgebraEngineGPU cuda_engine;
    mutable cudaDoubles cuda_surpluses;
    mutable cudaInts cuda_num_nodes, cuda_points;
    mutable cudaDoubles cuda_nodes, cuda_coeffs;
    #endif
};

}

#endif
