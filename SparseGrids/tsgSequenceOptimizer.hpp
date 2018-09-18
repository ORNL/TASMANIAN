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

#include <vector>

#include "tsgEnumerates.hpp"

namespace TasGrid{

class GreedySequences{
public:
    GreedySequences();
    ~GreedySequences();

    void getLejaNodes(int n, std::vector<double> &nodes) const;
    void getMaxLebesgueNodes(int n, std::vector<double> &nodes) const;
    void getMinLebesgueNodes(int n, std::vector<double> &nodes) const;
    void getMinDeltaNodes(int n, std::vector<double> &nodes) const;

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

class VectorFunctional{
public:
    VectorFunctional();
    virtual ~VectorFunctional();

    virtual double getValue(double x) const = 0;

    virtual bool hasDerivative() const = 0;
    virtual double getDiff(double x) const = 0;

    virtual void getIntervals(std::vector<double> &intervals) const = 0;

    static void makeCoeff(const std::vector<double> &nodes, std::vector<double> &coeffs);
    static void evalLag(const std::vector<double> &nodes, const std::vector<double> &coeffs, double x, std::vector<double> &lag);
    static double basisDx(const std::vector<double> &nodes, const std::vector<double> &coeffs, int inode, double x);
    static void sortIntervals(const std::vector<double> &nodes, std::vector<double> &intervals);
};

namespace Optimizer{
    OptimizerResult argMaxGlobal(const Functional *F); // assumes range is [-1,1] and both points are included
    OptimizerResult argMaxLocalPattern(const Functional *F, double left, double right);
    double argMaxLocalSecant(const Functional *F, double left, double right);

    OptimizerResult argMaxGlobal(const VectorFunctional &F); // assumes range is [-1,1] and both points are included
    OptimizerResult argMaxLocalPattern(const VectorFunctional &F, double left, double right);
    double argMaxLocalSecant(const VectorFunctional &F, double left, double right);
};

template<TypeOneDRule rule>
class tempFunctional : public VectorFunctional{
public:
    tempFunctional(const std::vector<double> &cnodes) : nodes(cnodes){
        if (rule != rule_leja) makeCoeff(nodes, coeff);
    }
    tempFunctional(const std::vector<double> &cnodes, double new_node){
        nodes.resize(cnodes.size() + 1);
        std::copy(cnodes.begin(), cnodes.end(), nodes.data());
        nodes[cnodes.size()] = new_node;
        makeCoeff(nodes, coeff);
        if (rule == rule_mindeltaodd){ // hack, rule_mindeltaodd indicates the companion functional to rule_mindelta (for the min-max problem)
            nodes_less1 = cnodes;
            makeCoeff(nodes_less1, coeff_less1);
        }
    }
    ~tempFunctional(){}

    double getValue(double x) const{
        if (rule == rule_leja){
            double p = 1.0;
            for(auto n : nodes) p *= (x - n);
            return fabs(p);
        }else if (rule == rule_maxlebesgue){
            std::vector<double> lag;
            evalLag(nodes, coeff, x, lag);
            double sum = 0.0;
            for(auto l : lag) sum += fabs(l);
            return sum;
        }else if (rule == rule_minlebesgue){
            // points almost overlap, will cause an explosion in the Lebesgue constant
            for(auto n : nodes) if (fabs(x - n) < 10*TSG_NUM_TOL) return -1.E+100;
            tempFunctional<rule_maxlebesgue> M(nodes, x);
            OptimizerResult R = Optimizer::argMaxGlobal(M);
            return -R.fmax;
        }else if (rule == rule_mindeltaodd){
            std::vector<double> lag, lag_less1;
            evalLag(nodes, coeff, x, lag);
            evalLag(nodes_less1, coeff_less1, x, lag_less1);
            auto il = lag.begin();
            double sum = 0.0;
            for(auto l : lag_less1) sum += fabs(l - *il++);
            return sum + fabs(*il);
        }else if (rule == rule_mindelta){
            // points almost overlap, will cause an explosion in the Lebesgue constant
            for(auto n : nodes) if (fabs(x - n) < 10*TSG_NUM_TOL) return -1.E+100;
            tempFunctional<rule_mindeltaodd> M(nodes, x);
            OptimizerResult R = Optimizer::argMaxGlobal(M);
            return -R.fmax;
        }
    }

    bool hasDerivative() const{
        return ((rule != rule_minlebesgue) && (rule != rule_mindelta));
    }

    double getDiff(double x) const{
        if (rule == rule_leja){
            double s = 1.0, p = 1.0, n = (x - nodes[0]);
            for(size_t j=1; j<nodes.size(); j++){
                p *= n;
                n = (x - nodes[j]);
                s *= n;
                s += p;
            }
            return s;
        }else if (rule == rule_maxlebesgue){
            std::vector<double> lag;
            evalLag(nodes, coeff, x, lag);
            double sum = 0.0;
            int i = 0;
            for(auto l : lag){
                double diff = basisDx(nodes, coeff, i++, x);
                sum += (l > 0.0) ? diff : -diff;
            }
            return sum;
        }else if (rule == rule_mindeltaodd){
            std::vector<double> lag, lag_less1;
            evalLag(nodes, coeff, x, lag);
            evalLag(nodes_less1, coeff_less1, x, lag_less1);
            auto il = lag.begin();
            int i = 0;
            double sum = 0.0;
            for(auto l : lag_less1){
                double diff = basisDx(nodes, coeff, i, x) - basisDx(nodes_less1, coeff_less1, i, x);
                sum += ((*il++ - l) > 0.0) ? diff : -diff;
                i++;
            }
            sum += basisDx(nodes, coeff, (int) (nodes.size() - 1), x) * ((*il > 0.0) ? 1.0 : -1.0);
            return sum;
        }
        return 0.0;
    }

    void getIntervals(std::vector<double> &intervals) const{
        sortIntervals(nodes, intervals);
    }

private:
    std::vector<double> nodes;
    std::vector<double> coeff;
    std::vector<double> nodes_less1;
    std::vector<double> coeff_less1;
};

}

#endif
