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

#include "tsgOneDimensionalWrapper.hpp"

namespace TasGrid{

OneDimensionalWrapper::OneDimensionalWrapper(const OneDimensionalMeta *meta, int max_level, TypeOneDRule crule, double alpha, double beta, std::ostream *logstream) :
    num_levels(max_level+1), rule(crule), indx(0), nodes(0)
{
    // find the points per level and the cumulative pointers
    isNonNested = OneDimensionalMeta::isNonNested(rule);
    num_points = new int[num_levels];
    pntr = new int[num_levels+1]; pntr[0] = 0;
    for(int l=0; l<num_levels; l++){
        num_points[l] = meta->getNumPoints(l, rule);
        pntr[l+1] = pntr[l] + num_points[l];
    }
    int num_total = pntr[num_levels];

    weights = new double[num_total];
    coeff = new double[num_total];

    if (isNonNested){
        indx    = new int[num_total];
        nodes   = new double[num_total];

        int num_unique = 0;
        unique = new double[num_total];

        double *x = 0, *w = 0;

        if (rule == rule_customtabulated){
            if (num_levels > meta->getCustom()->getNumLevels()){
                if (logstream != 0){ (*logstream) << "ERROR: custom-tabulated rule needed with levels " << num_levels << ", but only " << meta->getCustom()->getNumLevels() << " are provided." << endl; }
                #ifndef TASMANIAN_XSDK
                exit(1); // disabled by TASMANIAN_XSDK
                #else
                max_level = meta->getCustom()->getNumLevels();
                num_levels = max_levels;
                #endif // TASMANIAN_XSDK
            }
        }

        for(int l=0; l<num_levels; l++){
            int n = num_points[l];
            if ((rule == rule_chebyshev) || (rule == rule_chebyshevodd)){
                OneDimensionalNodes::getChebyshev(n, w, x);
            }else if ((rule == rule_gausslegendre) || (rule == rule_gausslegendreodd)){
                OneDimensionalNodes::getGaussLegendre(n, w, x);
            }else if ((rule == rule_gausschebyshev1) || (rule == rule_gausschebyshev1odd)){
                OneDimensionalNodes::getGaussChebyshev1(n, w, x);
            }else if ((rule == rule_gausschebyshev2) || (rule == rule_gausschebyshev2odd)){
                OneDimensionalNodes::getGaussChebyshev2(n, w, x);
            }else if ((rule == rule_gaussgegenbauer) || (rule == rule_gaussgegenbauerodd)){
                OneDimensionalNodes::getGaussJacobi(n, w, x, alpha, alpha);
            }else if ((rule == rule_gausshermite) || (rule == rule_gausshermiteodd)){
                OneDimensionalNodes::getGaussHermite(n, w, x, alpha);
            }else if ((rule == rule_gaussjacobi) || (rule == rule_gaussjacobiodd)){
                OneDimensionalNodes::getGaussJacobi(n, w, x, alpha, beta);
            }else if ((rule == rule_gausslaguerre) || (rule == rule_gausslaguerreodd)){
                OneDimensionalNodes::getGaussLaguerre(n, w, x, alpha);
            }else if (rule == rule_customtabulated){
                meta->getCustom()->getWeightsNodes(l, w, x);
            }

            for(int i=0; i<n; i++){
                weights[pntr[l] + i] = w[i];
                nodes[pntr[l] + i] = x[i];

                int point = -1;
                for(int j=0; j<num_unique; j++){
                    if (fabs(x[i] - unique[j]) < TSG_NUM_TOL){
                        point = j;
                        break;
                    }
                }
                if (point == - 1){ // new point found
                    unique[num_unique] = x[i];
                    point = num_unique++;
                }
                indx[pntr[l] + i] = point;
            }

        }

        double *t = unique;
        unique = new double[num_unique];
        std::copy(t, t + num_unique, unique);

        delete[] t;
        delete[] x;
        delete[] w;

        // make coefficients
        for(int l=0; l<num_levels; l++){
            int n = num_points[l];
            x = &(nodes[pntr[l]]);
            w = &(coeff[pntr[l]]);
            for(int i=0; i<n; i++){
                w[i] = 1.0;
                for(int j=0; j<i; j++){
                    w[i] /= (x[i] - x[j]);
                }
                for(int j=i+1; j<n; j++){
                    w[i] /= (x[i] - x[j]);
                }
            }
        }
    }else{
        TableGaussPatterson *gp;
        if (rule == rule_clenshawcurtis){
            unique = OneDimensionalNodes::getClenshawCurtisNodes(max_level);
        }else if (rule == rule_clenshawcurtis0){
            unique = OneDimensionalNodes::getClenshawCurtisNodesZero(max_level);
        }else if (rule == rule_fejer2){
            unique = OneDimensionalNodes::getFejer2Nodes(max_level);
        }else if (rule == rule_gausspatterson){
            gp = new TableGaussPatterson();
            if (num_levels > gp->getNumLevels()){
                if (logstream != 0){ (*logstream) << "ERROR: gauss-patterson rule needed with level " << max_level << ", but only " << gp->getNumLevels() << " are hardcoded." << endl; }
                #ifndef TASMANIAN_XSDK
                exit(1); // disabled by TASMANIAN_XSDK
                #else
                max_level = gp->getNumLevels()-1;
                num_levels = max_level+1;
                #endif // TASMANIAN_XSDK
            }
            unique = gp->getNodes(max_level);

            // move here to avoid a warning
            for(int l=0; l<num_levels; l++){
                for(int i=0; i<num_points[l]; i++){
                    weights[pntr[l] + i] = gp->getWeight(l, i);
                }
            }
            delete gp;
        }else if (rule == rule_rleja){
            unique = OneDimensionalNodes::getRLeja(meta->getNumPoints(max_level,rule));
        }else if ((rule == rule_rlejaodd) || (rule == rule_rlejadouble2) || (rule == rule_rlejadouble4)){
            unique = OneDimensionalNodes::getRLejaCentered(meta->getNumPoints(max_level,rule));
        }else if ((rule == rule_rlejashifted) || (rule == rule_rlejashiftedeven) || (rule == rule_rlejashifteddouble)){
            unique = OneDimensionalNodes::getRLejaShifted(meta->getNumPoints(max_level,rule));
        }else if ((rule == rule_leja) || (rule == rule_lejaodd)){
            GreedySequences greedy;
            unique = greedy.getLejaNodes(meta->getNumPoints(max_level, rule));
        }else if ((rule == rule_maxlebesgue) || (rule == rule_maxlebesgueodd)){
            GreedySequences greedy;
            unique = greedy.getMaxLebesgueNodes(meta->getNumPoints(max_level, rule));
        }else if ((rule == rule_minlebesgue) || (rule == rule_minlebesgueodd)){
            GreedySequences greedy;
            unique = greedy.getMinLebesgueNodes(meta->getNumPoints(max_level, rule));
        }else if ((rule == rule_mindelta) || (rule == rule_mindeltaodd)){
            GreedySequences greedy;
            unique = greedy.getMinDeltaNodes(meta->getNumPoints(max_level, rule));
        }

        for(int l=0; l<num_levels; l++){
            int n = num_points[l];
            double *c = &(coeff[pntr[l]]);
            for(int i=0; i<n; i++){
                c[i] = 1.0;
                for(int j=0; j<i; j++){
                    c[i] /= (unique[i] - unique[j]);
                }
                for(int j=i+1; j<n; j++){
                    c[i] /= (unique[i] - unique[j]);
                }
                if (rule == rule_clenshawcurtis0){
                    c[i] /= (unique[i] - 1.0) * (unique[i] + 1.0);
                }
            }
        }

        if (rule == rule_clenshawcurtis){
            for(int l=0; l<num_levels; l++){
                for(int i=0; i<num_points[l]; i++){
                    weights[pntr[l] + i] = OneDimensionalNodes::getClenshawCurtisWeight(l, i);
                }
            }
        }else if (rule == rule_clenshawcurtis0){
            for(int l=0; l<num_levels; l++){
                for(int i=0; i<num_points[l]; i++){
                    weights[pntr[l] + i] = OneDimensionalNodes::getClenshawCurtisWeightZero(l, i);
                }
            }
        }else if (rule == rule_fejer2){
            for(int l=0; l<num_levels; l++){
                for(int i=0; i<num_points[l]; i++){
                    weights[pntr[l] + i] = OneDimensionalNodes::getFejer2Weight(l, i);
                }
            }
        }else /*if (rule == rule_gausspatterson){
            // weights are normally computed here, but I put the GP weights above, since they are only read once from a table
        }else*/ if ((rule == rule_rleja) || (rule == rule_rlejaodd) || (rule == rule_rlejadouble2) || (rule == rule_rlejadouble4) || (rule == rule_leja) || (rule == rule_lejaodd) ||
                (rule == rule_rlejashifted) || (rule == rule_rlejashiftedeven) || (rule == rule_rlejashifteddouble) ||
               (rule == rule_maxlebesgue) || (rule == rule_maxlebesgueodd) || (rule == rule_minlebesgue) || (rule == rule_minlebesgueodd) || (rule == rule_mindelta) || (rule == rule_mindeltaodd)){
            int n = 1 + num_points[max_level] / 2; // number of Gauss-Legendre points needed to integrate the basis functions
            double *lag_x = 0, *lag_w = 0;
            OneDimensionalNodes::getGaussLegendre(n, lag_w, lag_x);
            std::fill(weights, weights+num_total, 0.0);
            for(int l=0; l<num_levels; l++){
                int npl = num_points[l];
                double *v = new double[npl];
                double *c = &(coeff[pntr[l]]);
                double *w = &(weights[pntr[l]]);
                for(int i=0; i<n; i++){
                    v[0] = 1.0;
                    for(int j=0; j<npl-1; j++){
                        v[j+1] = (lag_x[i] - unique[j]) * v[j];
                    }
                    v[npl-1] *= c[npl-1];
                    double s = 1.0;
                    for(int j=npl-2; j>=0; j--){
                        s *= (lag_x[i] - unique[j+1]);
                        v[j] *= s * c[j];
                    }
                    for(int j=0; j<npl; j++){
                        w[j] += lag_w[i] * v[j];
                    }
                }
                delete[] v;
            }
            delete[] lag_w;
            delete[] lag_x;
        }
    }
}

OneDimensionalWrapper::~OneDimensionalWrapper(){
    delete[] num_points;
    delete[] pntr;
    delete[] indx;
    delete[] weights;
    delete[] nodes;
    delete[] unique;
    delete[] coeff;
}

int OneDimensionalWrapper::getNumPoints(int level) const{  return num_points[level];  }
int OneDimensionalWrapper::getPointIndex(int level, int j) const{  return indx[pntr[level] + j];  }

double OneDimensionalWrapper::getNode(int j) const{  return unique[j];  }
double OneDimensionalWrapper::getWeight(int level, int j) const{  return weights[pntr[level] + j];  }

const double* OneDimensionalWrapper::getNodes(int level) const{  return (isNonNested) ? &(nodes[pntr[level]]) : unique;  }
const double* OneDimensionalWrapper::getCoefficients(int level) const{  return &(coeff[pntr[level]]);  }

int OneDimensionalWrapper::getPointsCount(int level) const{  return pntr[level]; }

TypeOneDRule OneDimensionalWrapper::getType() const{ return rule; }
int OneDimensionalWrapper::getNumLevels() const{ return num_levels; }

}

#endif
