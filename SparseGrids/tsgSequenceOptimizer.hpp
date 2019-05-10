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

namespace Optimizer{

class VectorFunctional{
public:
    VectorFunctional();
    virtual ~VectorFunctional();

    virtual double getValue(double x) const = 0;

    virtual bool hasDerivative() const = 0;
    virtual double getDiff(double x) const = 0;

    virtual void getIntervals(std::vector<double> &intervals) const = 0;
};

struct OptimizerResult{
    double xmax, fmax;
};

OptimizerResult argMaxGlobal(const VectorFunctional &F); // assumes range is [-1,1] and both points are included
OptimizerResult argMaxLocalPattern(const VectorFunctional &F, double left, double right);
double argMaxLocalSecant(const VectorFunctional &F, double left, double right);


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
            return std::abs(p);
        }else if (rule == rule_maxlebesgue){
            std::vector<double> lag;
            evalLag(nodes, coeff, x, lag);
            double sum = 0.0;
            for(auto l : lag) sum += std::abs(l);
            return sum;
        }else if (rule == rule_minlebesgue){
            // points almost overlap, will cause an explosion in the Lebesgue constant
            for(auto n : nodes) if (std::abs(x - n) < 10 * Maths::num_tol) return -1.E+100;
            tempFunctional<rule_maxlebesgue> M(nodes, x);
            OptimizerResult R = Optimizer::argMaxGlobal(M);
            return -R.fmax;
        }else if (rule == rule_mindeltaodd){
            std::vector<double> lag, lag_less1;
            evalLag(nodes, coeff, x, lag);
            evalLag(nodes_less1, coeff_less1, x, lag_less1);
            auto il = lag.begin();
            double sum = 0.0;
            for(auto l : lag_less1) sum += std::abs(l - *il++);
            return sum + std::abs(*il);
        }else{ // rule == rule_mindelta
            // points almost overlap, will cause an explosion in the Lebesgue constant
            for(auto n : nodes) if (std::abs(x - n) < 10 * Maths::num_tol) return -1.E+100;
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
        }else{ // rule == rule_mindeltaodd
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
    }

    void getIntervals(std::vector<double> &intervals) const{
        sortIntervals(nodes, intervals);
    }

    static void makeCoeff(const std::vector<double> &nodes, std::vector<double> &coeffs){
        size_t num_nodes = nodes.size();
        coeffs.resize(num_nodes);
        for(size_t i=0; i<num_nodes; i++){
            double c = 1.0;
            for(size_t j=0; j<i; j++){
                c *= (nodes[i] - nodes[j]);
            }
            for(size_t j=i+1; j<num_nodes; j++){
                c *= (nodes[i] - nodes[j]);
            }
            coeffs[i] = c;
        }
    }
    static void evalLag(const std::vector<double> &nodes, const std::vector<double> &coeffs, double x, std::vector<double> &lag){
        int num_nodes = (int) nodes.size();
        lag.resize(num_nodes);
        lag[0] = 1.0;
        for(int i=0; i<num_nodes-1; i++){
            lag[i+1] = (x - nodes[i]) * lag[i];
        }
        double w = 1.0;
        lag[num_nodes-1] /= coeffs[num_nodes-1];
        for(int i= num_nodes-2; i>=0; i--){
            w *= (x - nodes[i+1]);
            lag[i] *= w / coeffs[i];
        }
    }
    static double basisDx(const std::vector<double> &nodes, const std::vector<double> &coeffs, int inode, double x){
        size_t num_nodes = nodes.size();
        double s = 1.0;
        double p = 1.0;
        double n = (inode != 0) ? (x - nodes[0]) : (x - nodes[1]);

        for(int j=1; j<inode; j++){
            p *= n;
            n = (x - nodes[j]);
            s *= n;
            s += p;
        }
        for(size_t j = (size_t) ((inode == 0) ? 2 : inode+1); j<num_nodes; j++){
            p *= n;
            n = (x - nodes[j]);
            s *= n;
            s += p;
        }
        return s / coeffs[inode];
    }
    static void sortIntervals(const std::vector<double> &nodes, std::vector<double> &intervals){
        size_t num_nodes = nodes.size();

        std::vector<double> v1(num_nodes), v2(num_nodes);
        size_t loop_end = (num_nodes % 2 == 1) ? num_nodes - 1 : num_nodes;
        for(size_t i=0; i<loop_end; i+=2){
            if (nodes[i] < nodes[i+1]){
                v1[i]   = nodes[i];
                v1[i+1] = nodes[i+1];
            }else{
                v1[i]   = nodes[i+1];
                v1[i+1] = nodes[i];
            }
        }
        if (loop_end != num_nodes){
            v1[num_nodes-1] = nodes[num_nodes-1];
        }

        size_t stride = 2;
        while(stride < num_nodes){
            auto ic = v2.begin();
            auto i2 = v1.begin();
            while(ic < v2.end()){
                auto ia = i2;
                auto ib = i2;
                std::advance(ib, stride);
                auto iaend = ia;
                std::advance(iaend, stride);
                if (iaend > v1.end()) iaend = v1.end();
                auto ibend = ib;
                std::advance(ibend, stride);
                if (ibend > v1.end()) ibend = v1.end();
                bool aless = (ia < iaend);
                bool bless = (ib < ibend);
                while(aless || bless){
                    bool picka;
                    if (aless && bless){
                        picka = *ia < *ib;
                    }else{
                        picka = aless;
                    }

                    if (picka){
                        *ic++ = *ia++;
                    }else{
                        *ic++ = *ib++;
                    }

                    aless = (ia < iaend);
                    bless = (ib < ibend);
                }

                std::advance(i2, 2 * stride);
            }
            stride *= 2;
            std::swap(v1, v2);
        }

        intervals = std::move(v1);
    }

private:
    std::vector<double> nodes;
    std::vector<double> coeff;
    std::vector<double> nodes_less1;
    std::vector<double> coeff_less1;
};

template<TypeOneDRule rule> double getNextNode(std::vector<double> const &nodes);

std::vector<double> getPrecomputedMinLebesgueNodes();
std::vector<double> getPrecomputedMinDeltaNodes();

template<TypeOneDRule rule> std::vector<double> getGreedyNodes(int n);

}

}

#endif
