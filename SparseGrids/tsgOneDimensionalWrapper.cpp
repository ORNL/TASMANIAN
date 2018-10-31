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

#ifndef __TSG_ONE_DIMENSIONAL_WRAPPER_CPP
#define __TSG_ONE_DIMENSIONAL_WRAPPER_CPP

#include <stdexcept>
#include <string>

#include "tsgOneDimensionalWrapper.hpp"

namespace TasGrid{

OneDimensionalWrapper::OneDimensionalWrapper() : num_levels(0), rule(rule_none){}

void OneDimensionalWrapper::load(const CustomTabulated &custom, int max_level, TypeOneDRule crule, double alpha, double beta){
    num_levels = max_level + 1;
    rule = crule;

    // find the points per level and the cumulative pointers
    isNonNested = OneDimensionalMeta::isNonNested(rule);
    num_points.resize(num_levels);
    pntr.resize(num_levels+1); pntr[0] = 0;
    if (rule != rule_customtabulated){
        for(int l=0; l<num_levels; l++){
            num_points[l] = OneDimensionalMeta::getNumPoints(l, rule);
            pntr[l+1] = pntr[l] + num_points[l];
        }
    }else{
        for(int l=0; l<num_levels; l++){
            num_points[l] = custom.getNumPoints(l);
            pntr[l+1] = pntr[l] + num_points[l];
        }
    }
    int num_total = pntr[num_levels];

    weights.resize(num_levels);
    nodes.resize(num_levels);
    coeff.resize(num_levels);

    if (isNonNested){
        indx.reserve(num_total);
        unique.reserve(num_total);

        if (rule == rule_customtabulated){
            if (num_levels > custom.getNumLevels()){
                std::string message = "ERROR: custom-tabulated rule needed with levels ";
                message += std::to_string(num_levels);
                message += ", but only ";
                message += std::to_string(custom.getNumLevels());
                message += " are provided.";
                throw std::runtime_error(message);
            }
        }

        for(int l=0; l<num_levels; l++){
            int n = num_points[l];
            if ((rule == rule_chebyshev) || (rule == rule_chebyshevodd)){
                OneDimensionalNodes::getChebyshev(n, weights[l], nodes[l]);
            }else if ((rule == rule_gausslegendre) || (rule == rule_gausslegendreodd)){
                OneDimensionalNodes::getGaussLegendre(n, weights[l], nodes[l]);
            }else if ((rule == rule_gausschebyshev1) || (rule == rule_gausschebyshev1odd)){
                OneDimensionalNodes::getGaussChebyshev1(n, weights[l], nodes[l]);
            }else if ((rule == rule_gausschebyshev2) || (rule == rule_gausschebyshev2odd)){
                OneDimensionalNodes::getGaussChebyshev2(n, weights[l], nodes[l]);
            }else if ((rule == rule_gaussgegenbauer) || (rule == rule_gaussgegenbauerodd)){
                OneDimensionalNodes::getGaussJacobi(n, weights[l], nodes[l], alpha, alpha);
            }else if ((rule == rule_gausshermite) || (rule == rule_gausshermiteodd)){
                OneDimensionalNodes::getGaussHermite(n, weights[l], nodes[l], alpha);
            }else if ((rule == rule_gaussjacobi) || (rule == rule_gaussjacobiodd)){
                OneDimensionalNodes::getGaussJacobi(n, weights[l], nodes[l], alpha, beta);
            }else if ((rule == rule_gausslaguerre) || (rule == rule_gausslaguerreodd)){
                OneDimensionalNodes::getGaussLaguerre(n, weights[l], nodes[l], alpha);
            }else if (rule == rule_customtabulated){
                custom.getWeightsNodes(l, weights[l], nodes[l]);
            }

            for(auto x : nodes[l]){
                int point = -1;
                for(int j=0; j<(int) unique.size(); j++){
                    if (fabs(x - unique[j]) < TSG_NUM_TOL){
                        point = j;
                        break;
                    }
                }
                if (point == - 1){ // new point found
                    unique.push_back(x);
                    point = (int) (unique.size() - 1);
                }
                indx.push_back(point);
            }
        }

        unique.shrink_to_fit();

        // make coefficients
        for(int l=0; l<num_levels; l++){
            int n = num_points[l];
            coeff[l].resize(n);
            double *x = nodes[l].data();
            for(int i=0; i<n; i++){
                double c = 1.0;
                for(int j=0; j<i; j++){
                    c /= (x[i] - x[j]);
                }
                for(int j=i+1; j<n; j++){
                    c /= (x[i] - x[j]);
                }
                coeff[l][i] = c;
            }
        }
    }else{
        if (rule == rule_clenshawcurtis){
            OneDimensionalNodes::getClenshawCurtisNodes(max_level, unique);
        }else if (rule == rule_clenshawcurtis0){
            OneDimensionalNodes::getClenshawCurtisNodesZero(max_level, unique);
        }else if (rule == rule_fejer2){
            OneDimensionalNodes::getFejer2Nodes(max_level, unique);
        }else if (rule == rule_gausspatterson){
            TableGaussPatterson gp;
            if (num_levels > gp.getNumLevels()){
                std::string message = "ERROR: gauss-patterson rule needed with level ";
                message += std::to_string(max_level);
                message += ", but only ";
                message += std::to_string(gp.getNumLevels());
                message += " are hardcoded.";
                throw std::runtime_error(message);
            }
            gp.getNodes(max_level, unique);

            // move here to avoid a warning
            for(int l=0; l<num_levels; l++){
                weights[l].resize(num_points[l]);
                for(int i=0; i<num_points[l]; i++){
                    weights[l][i] = gp.getWeight(l, i);
                }
            }
        }else if (rule == rule_rleja){
            OneDimensionalNodes::getRLeja(OneDimensionalMeta::getNumPoints(max_level,rule), unique);
        }else if ((rule == rule_rlejaodd) || (rule == rule_rlejadouble2) || (rule == rule_rlejadouble4)){
            OneDimensionalNodes::getRLejaCentered(OneDimensionalMeta::getNumPoints(max_level,rule), unique);
        }else if ((rule == rule_rlejashifted) || (rule == rule_rlejashiftedeven) || (rule == rule_rlejashifteddouble)){
            OneDimensionalNodes::getRLejaShifted(OneDimensionalMeta::getNumPoints(max_level,rule), unique);
        }else if ((rule == rule_leja) || (rule == rule_lejaodd)){
            Optimizer::getGreedyNodes<rule_leja>(OneDimensionalMeta::getNumPoints(max_level, rule), unique);
        }else if ((rule == rule_maxlebesgue) || (rule == rule_maxlebesgueodd)){
            Optimizer::getGreedyNodes<rule_maxlebesgue>(OneDimensionalMeta::getNumPoints(max_level, rule), unique);
        }else if ((rule == rule_minlebesgue) || (rule == rule_minlebesgueodd)){
            Optimizer::getGreedyNodes<rule_minlebesgue>(OneDimensionalMeta::getNumPoints(max_level, rule), unique);
        }else if ((rule == rule_mindelta) || (rule == rule_mindeltaodd)){
            Optimizer::getGreedyNodes<rule_mindelta>(OneDimensionalMeta::getNumPoints(max_level, rule), unique);
        }else{ // if (rule==rule_fourier)
            OneDimensionalNodes::getFourierNodes(max_level, unique);
        }

        for(int l=0; l<num_levels; l++){
            int n = num_points[l];
            weights[l].resize(n);
            coeff[l].resize(n);
            for(int i=0; i<n; i++){
                double c = 1.0;
                for(int j=0; j<i; j++){
                    c /= (unique[i] - unique[j]);
                }
                for(int j=i+1; j<n; j++){
                    c /= (unique[i] - unique[j]);
                }
                if (rule == rule_clenshawcurtis0){
                    c /= (unique[i] - 1.0) * (unique[i] + 1.0);
                }
                coeff[l][i] = c;
            }
        }

        if (rule == rule_clenshawcurtis){
            for(int l=0; l<num_levels; l++){
                for(int i=0; i<num_points[l]; i++){
                    weights[l][i] = OneDimensionalNodes::getClenshawCurtisWeight(l, i);
                }
            }
        }else if (rule == rule_clenshawcurtis0){
            for(int l=0; l<num_levels; l++){
                for(int i=0; i<num_points[l]; i++){
                    weights[l][i] = OneDimensionalNodes::getClenshawCurtisWeightZero(l, i);
                }
            }
        }else if (rule == rule_fejer2){
            for(int l=0; l<num_levels; l++){
                for(int i=0; i<num_points[l]; i++){
                    weights[l][i] = OneDimensionalNodes::getFejer2Weight(l, i);
                }
            }
        }else /*if (rule == rule_gausspatterson){
            // weights are normally computed here, but I put the GP weights above, since they are only read once from a table
        }else*/ if ((rule == rule_rleja) || (rule == rule_rlejaodd) || (rule == rule_rlejadouble2) || (rule == rule_rlejadouble4) || (rule == rule_leja) || (rule == rule_lejaodd) ||
                (rule == rule_rlejashifted) || (rule == rule_rlejashiftedeven) || (rule == rule_rlejashifteddouble) ||
               (rule == rule_maxlebesgue) || (rule == rule_maxlebesgueodd) || (rule == rule_minlebesgue) || (rule == rule_minlebesgueodd) || (rule == rule_mindelta) || (rule == rule_mindeltaodd)){
            int n = 1 + num_points[max_level] / 2; // number of Gauss-Legendre points needed to integrate the basis functions
            std::vector<double> lag_x, lag_w;
            OneDimensionalNodes::getGaussLegendre(n, lag_w, lag_x);
            for(int l=0; l<num_levels; l++){
                std::fill(weights[l].begin(), weights[l].end(), 0.0);
                int npl = num_points[l];
                std::vector<double> v(npl);
                for(int i=0; i<n; i++){
                    v[0] = 1.0;
                    for(int j=0; j<npl-1; j++){
                        v[j+1] = (lag_x[i] - unique[j]) * v[j];
                    }
                    v[npl-1] *= coeff[l][npl-1];
                    double s = 1.0;
                    for(int j=npl-2; j>=0; j--){
                        s *= (lag_x[i] - unique[j+1]);
                        v[j] *= s * coeff[l][j];
                    }
                    for(int j=0; j<npl; j++){
                        weights[l][j] += lag_w[i] * v[j];
                    }
                }
            }
        }
    }
}

OneDimensionalWrapper::~OneDimensionalWrapper(){}

int OneDimensionalWrapper::getNumPoints(int level) const{ return num_points[level]; }
int OneDimensionalWrapper::getPointIndex(int level, int j) const{ return indx[pntr[level] + j]; }

double OneDimensionalWrapper::getNode(int j) const{ return unique[j]; }
double OneDimensionalWrapper::getWeight(int level, int j) const{ return weights[level][j];  }

const double* OneDimensionalWrapper::getNodes(int level) const{ return (isNonNested) ? nodes[level].data() : unique.data(); }
const double* OneDimensionalWrapper::getCoefficients(int level) const{ return coeff[level].data(); }

const std::vector<int>* OneDimensionalWrapper::getPointsCount() const{ return &pntr; }

TypeOneDRule OneDimensionalWrapper::getType() const{ return rule; }
int OneDimensionalWrapper::getNumLevels() const{ return num_levels; }

}

#endif
