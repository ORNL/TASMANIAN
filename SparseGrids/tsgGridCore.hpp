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

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"

namespace TasGrid{

class BaseCanonicalGrid{
public:
    BaseCanonicalGrid();
    virtual ~BaseCanonicalGrid();

    virtual int getNumDimensions() const = 0;
    virtual int getNumOutputs() const = 0;
    virtual TypeOneDRule getRule() const = 0;

    virtual int getNumLoaded() const = 0;
    virtual int getNumNeeded() const = 0;
    virtual int getNumPoints() const = 0;

    virtual double* getLoadedPoints() const = 0;
    virtual void getLoadedPoints(double *x) const = 0;
    virtual double* getNeededPoints() const = 0;
    virtual void getNeededPoints(double *x) const = 0;
    virtual double* getPoints() const = 0;
    virtual void getPoints(double *x) const = 0;

    virtual double* getQuadratureWeights() const = 0;
    virtual void getQuadratureWeights(double weights[]) const = 0;
    virtual double* getInterpolationWeights(const double x[]) const = 0;
    virtual void getInterpolationWeights(const double x[], double weights[]) const = 0;

    virtual void loadNeededPoints(const double *vals, TypeAcceleration acc) = 0;

    virtual void evaluate(const double x[], double y[]) const = 0;
    virtual void integrate(double q[], double *conformal_correction) const = 0;

    virtual void evaluateFastCPUblas(const double x[], double y[]) const = 0;

    virtual void evaluateFastGPUcublas(const double x[], double y[], std::ostream *os) const = 0;
    virtual void evaluateFastGPUcuda(const double x[], double y[], std::ostream *os) const = 0;
    virtual void evaluateFastGPUmagma(const double x[], double y[], std::ostream *os) const = 0;

    virtual void evaluateBatch(const double x[], int num_x, double y[]) const = 0;
    virtual void evaluateBatchCPUblas(const double x[], int num_x, double y[]) const = 0;

    virtual void evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *os) const = 0;
    virtual void evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream *os) const = 0;
    virtual void evaluateBatchGPUmagma(const double x[], int num_x, double y[], std::ostream *os) const = 0;

    virtual void clearRefinement() = 0;
    virtual void mergeRefinement() = 0;

    virtual void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const = 0; // add acceleration here
    virtual void setHierarchicalCoefficients(const double c[], TypeAcceleration acc, std::ostream *os) = 0;

    virtual void clearAccelerationData() = 0;
};

class SplitDirections{
public:
    SplitDirections(const IndexSet *points);
    ~SplitDirections();

    int getNumJobs() const;
    int getJobDirection(int job) const;
    int getJobNumPoints(int job) const;
    const int* getJobPoints(int job) const;

protected:
    bool doesBelongSameLine(const int a[], const int b[], int direction) const;

private:
    int num_dimensions, num_allocated_jobs, num_jobs;
    int *job_directions, *num_job_pnts, **job_pnts;
};

}

#endif
