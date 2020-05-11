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

#ifndef __TSG_BASE_CLASS_HPP
#define __TSG_BASE_CLASS_HPP

/*!
 * \internal
 * \file tsgGridCore.hpp
 * \brief Definition for the base grid class.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSets
 *
 * Definition of the BaseCanonicalGrid class used by all grids.
 * \endinternal
 */

#include "tsgDConstructGridGlobal.hpp"
#include "tsgCudaLoadStructures.hpp"
#include "tsgHierarchyManipulator.hpp"

#ifndef __TASMANIAN_DOXYGEN_SKIP
namespace TasGrid{

class BaseCanonicalGrid{
public:
    BaseCanonicalGrid(AccelerationContext const *acc) : acceleration(acc), num_dimensions(0), num_outputs(0){}
    BaseCanonicalGrid(AccelerationContext const *acc, int cnum_dimensions, int cnum_outputs, MultiIndexSet &&cpoints, MultiIndexSet &&cneeded, StorageSet &&cvalues)
        : acceleration(acc), num_dimensions(cnum_dimensions), num_outputs(cnum_outputs),
          points(std::forward<MultiIndexSet>(cpoints)), needed(std::forward<MultiIndexSet>(cneeded)),
          values(std::forward<StorageSet>(cvalues)){}
    BaseCanonicalGrid(AccelerationContext const *acc, BaseCanonicalGrid const &other, int ibegin, int iend)
        : acceleration(acc), num_dimensions(other.num_dimensions), num_outputs(iend - ibegin),
          points(other.points), needed(other.needed),
          values((iend - ibegin == other.num_outputs) ? other.values : other.values.splitValues(ibegin, iend))
          {}
    virtual ~BaseCanonicalGrid() = default;

    virtual bool isGlobal() const{ return false; }
    virtual bool isSequence() const{ return false; }
    virtual bool isLocalPolynomial() const{ return false; }
    virtual bool isWavelet() const{ return false; }
    virtual bool isFourier() const{ return false; }

    int getNumDimensions() const{ return num_dimensions; }
    int getNumOutputs() const{ return num_outputs; }
    virtual TypeOneDRule getRule() const = 0;

    int getNumLoaded() const{ return (num_outputs == 0) ? 0 : points.getNumIndexes(); }
    int getNumNeeded() const{ return needed.getNumIndexes(); }
    int getNumPoints() const{ return ((points.empty()) ? needed.getNumIndexes() : points.getNumIndexes()); }
    const int* getPointIndexes() const{ return ((points.empty()) ? needed.getIndex(0) : points.getIndex(0)); }

    virtual void write(std::ostream&, bool) const = 0;

    virtual void getLoadedPoints(double *x) const = 0;
    virtual void getNeededPoints(double *x) const = 0;
    virtual void getPoints(double *x) const = 0;

    virtual void getQuadratureWeights(double weights[]) const = 0;
    virtual void getInterpolationWeights(const double x[], double weights[]) const = 0;

    virtual void loadNeededPoints(const double *vals) = 0;
    const double* getLoadedValues() const{ return (points.empty()) ? nullptr : values.getValues(0); }

    virtual void evaluate(const double x[], double y[]) const = 0;
    virtual void integrate(double q[], double *conformal_correction) const = 0;

    virtual void evaluateBatch(const double x[], int num_x, double y[]) const = 0;

    virtual void evaluateBatchGPU(const double[], int, double[]) const = 0;
    virtual void evaluateBatchGPU(const float[], int, float[]) const = 0;
    virtual void evaluateHierarchicalFunctionsGPU(const double[], int, double[]) const = 0;
    virtual void evaluateHierarchicalFunctionsGPU(const float[], int, float[]) const = 0;

    virtual void clearRefinement() = 0;
    virtual void mergeRefinement() = 0;

    virtual void beginConstruction() = 0;
    virtual void writeConstructionData(std::ostream&, bool) const = 0;
    virtual void readConstructionData(std::istream&, bool) = 0;
    virtual void loadConstructedPoint(const double[], const std::vector<double> &) = 0;
    virtual void loadConstructedPoint(const double[], int, const double[]) = 0;
    virtual void finishConstruction() = 0;

    virtual void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const = 0; // add acceleration here
    virtual std::vector<double> getSupport() const{ return std::vector<double>(Utils::size_mult(getNumPoints(), getNumDimensions()), 2.0); }
    virtual void setHierarchicalCoefficients(const double c[]) = 0;
    virtual void integrateHierarchicalFunctions(double integrals[]) const = 0;

    virtual void updateAccelerationData(AccelerationContext::ChangeType change) const = 0;

protected:
    AccelerationContext const *acceleration;
    int num_dimensions, num_outputs;
    MultiIndexSet points;
    MultiIndexSet needed;
    StorageSet values;
};

// For purposes of reading grids with older file formats
// Specialize the reader class for each grid type
// Befriend the reader class and each grid type
template<class GridType>
struct GridReaderVersion5{ // 5 refers to the file format version, not the Tasmanian version
    template<typename iomode> static std::unique_ptr<BaseCanonicalGrid> read(AccelerationContext const*, std::istream &){
        return std::unique_ptr<BaseCanonicalGrid>();
    }
};

// Factory reader method that instantiates the GridReaderVersion5 class
template<class GridType, typename iomode> std::unique_ptr<GridType> readGridVersion5(AccelerationContext const *acc, std::istream &is, iomode){
    return GridReaderVersion5<GridType>::template read<iomode>(acc, is);
}

}
#endif

#endif
