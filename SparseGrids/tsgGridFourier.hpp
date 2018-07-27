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

#ifndef __TASMANIAN_SPARSE_GRID_FOURIER_HPP
#define __TASMANIAN_SPARSE_GRID_FOURIER_HPP

#include <cstdlib>
#include <math.h>
#include <complex>

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgIndexManipulator.hpp"
#include "tsgLinearSolvers.hpp"
#include "tsgOneDimensionalWrapper.hpp"
#include "tsgGridCore.hpp"

#include "tsgAcceleratedDataStructures.hpp"

namespace TasGrid{

class GridFourier : public BaseCanonicalGrid {
public:
    GridFourier();
    GridFourier(const GridFourier &fourier);
    ~GridFourier();

    void write(std::ofstream &ofs) const;
    void read(std::ifstream &ifs, std::ostream *logstream = 0);

    void writeBinary(std::ofstream &ofs) const;
    void readBinary(std::ifstream &ifs, std::ostream *logstream = 0);

    void makeGrid(int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, const int* anisotropic_weights = 0, const int* level_limits = 0);
    void copyGrid(const GridFourier *fourier);

    void setTensors(IndexSet* &tset, int cnum_outputs);
    int* referenceExponents(const int levels[], const IndexSet *list);

    int getNumDimensions() const;
    int getNumOutputs() const;
    TypeOneDRule getRule() const;

    int getNumLoaded() const;
    int getNumNeeded() const;
    int getNumPoints() const; // returns the number of loaded points unless no points are loaded, then returns the number of needed points

    void loadNeededPoints(const double *vals, TypeAcceleration acc = accel_none);

    double* getLoadedPoints() const;
    void getLoadedPoints(double *x) const;
    double* getNeededPoints() const;
    void getNeededPoints(double *x) const;
    double* getPoints() const;
    void getPoints(double *x) const; // returns the loaded points unless no points are loaded, then returns the needed points

    double* getInterpolationWeights(const double x[]) const;
    void getInterpolationWeights(const double x[], double weights[]) const;

    double* getQuadratureWeights() const;
    void getQuadratureWeights(double weights[]) const;

    void evaluate(const double x[], double y[]) const;
    void evaluateBatch(const double x[], int num_x, double y[]) const;

    void evaluateFastCPUblas(const double x[], double y[]) const;
    void evaluateFastGPUcublas(const double x[], double y[], std::ostream *os) const;
    void evaluateFastGPUcuda(const double x[], double y[], std::ostream *os) const;
    void evaluateFastGPUmagma(int gpuID, const double x[], double y[], std::ostream *os) const;

    void evaluateBatchCPUblas(const double x[], int num_x, double y[]) const;
    void evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *os) const;
    void evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream *os) const;
    void evaluateBatchGPUmagma(int gpuID, const double x[], int num_x, double y[], std::ostream *os) const;

    void integrate(double q[], double *conformal_correction) const;

    void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const;
    void setHierarchicalCoefficients(const double c[], TypeAcceleration acc, std::ostream *os);

    void clearAccelerationData();
    void clearRefinement();
    void mergeRefinement();

    const int* getPointIndexes() const;
    const IndexSet* getExponents() const;
    const double* getFourierCoefs() const;

protected:
    void reset();
    void calculateFourierCoefficients();

    std::complex<double>* getBasisFunctions(const double x[]) const;
    void getBasisFunctions(const double x[], std::complex<double> weights[]) const;
    void getBasisFunctions(const double x[], double weights[]) const;   // for evaluateHierarchicalFunctions (parallelized)

    int convertIndexes(const int i, const int levels[]) const;

private:
    int num_dimensions, num_outputs;

    OneDimensionalWrapper *wrapper;

    IndexSet *tensors;
    IndexSet *active_tensors;
    int *active_w;
    int *max_levels;
    IndexSet *points;
    IndexSet *needed;
    IndexSet *exponents;

    double *fourier_coefs;
    int **exponent_refs;
    int **tensor_refs;

    StorageSet *values;

    mutable BaseAccelerationData *accel;

};

}

#endif
