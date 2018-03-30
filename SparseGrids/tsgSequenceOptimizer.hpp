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

#ifndef __TASMANIAN_SPARSE_GRID_SEQUENCE_OPTIMIZER_HPP
#define __TASMANIAN_SPARSE_GRID_SEQUENCE_OPTIMIZER_HPP

#include "tsgEnumerates.hpp"

namespace TasGrid{

class GreedySequences{
public:
    GreedySequences();
    ~GreedySequences();

    double* getLejaNodes(int n) const;
    double* getMaxLebesgueNodes(int n) const;
    double* getMinLebesgueNodes(int n) const;
    double* getMinDeltaNodes(int n) const;

    double findLebesgueConstant(int n, const double nodes[]) const;

    int getNumMinLebesgueStored() const;
    double getMinLebesgueStored(int i) const;
    int getNumMinDeltaStored() const;
    double getMinDeltaStored(int i) const;
};

class Functional{
public:
    Functional();
    ~Functional();

    virtual double getValue(double x) const = 0;

    virtual bool hasDerivative() const = 0;
    virtual double getDiff(double x) const;

    virtual int getNumIntervals() const = 0;
    virtual double* getIntervals() const = 0;

    static double* makeCoeff(int num_nodes, const double nodes[]); // may not be the best place to put this
    static double* evalLag(int num_nodes, const double nodes[], const double coeff[], double x);
    static double basisDx(int num_nodes, const double nodes[], const double coeff[], int inode, double x);
    static double* sortIntervals(int num_nodes, const double nodes[]);
    static int nodeCompar(const void * a, const void * b);
};

class Residual : public Functional{
public:
    Residual(int cnum_nodes, const double cnodes[]);
    ~Residual();

    double getValue(double x) const;

    bool hasDerivative() const;
    double getDiff(double x) const;

    int getNumIntervals() const;
    double* getIntervals() const;

private:
    int num_nodes;
    double *nodes;
};

class MaxLebesgue : public Functional{
public:
    MaxLebesgue(int cnum_nodes, const double cnodes[], double new_node);
    MaxLebesgue(int cnum_nodes, const double cnodes[]);
    ~MaxLebesgue();

    double getValue(double x) const;

    bool hasDerivative() const;
    double getDiff(double x) const;

    int getNumIntervals() const;
    double* getIntervals() const;

private:
    int num_nodes;
    double *nodes;
    double *coeff;
};

class MinLebesgue : public Functional{
public:
    MinLebesgue(int cnum_nodes, const double cnodes[]);
    ~MinLebesgue();

    double getValue(double x) const;
    bool hasDerivative() const;

    int getNumIntervals() const;
    double* getIntervals() const;

private:
    int num_nodes;
    double *nodes;
};

class MaxDelta : public Functional{
public:
    MaxDelta(int cnum_nodes, const double cnodes[], double new_node);
    MaxDelta(int cnum_nodes, const double cnodes[]);
    ~MaxDelta();

    double getValue(double x) const;

    bool hasDerivative() const;
    double getDiff(double x) const;

    int getNumIntervals() const;
    double* getIntervals() const;

private:
    int num_nodes;
    double *nodes;
    double *coeff_plus, *coeff_minus;
};

class MinDelta : public Functional{
public:
    MinDelta(int cnum_nodes, const double cnodes[]);
    ~MinDelta();

    double getValue(double x) const;
    bool hasDerivative() const;

    int getNumIntervals() const;
    double* getIntervals() const;

private:
    int num_nodes;
    double *nodes;
};

struct OptimizerResult{
    double xmax, fmax;
};

class Optimizer{
public:
    Optimizer();
    ~Optimizer();

    static OptimizerResult argMaxGlobal(const Functional *F); // assumes range is [-1,1] and both points are included
    static OptimizerResult argMaxLocalPattern(const Functional *F, double left, double right);
    static double argMaxLocalSecant(const Functional *F, double left, double right);
};

}

#endif
