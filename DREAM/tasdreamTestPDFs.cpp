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

#ifndef __TASMANIAN_TASDREAM_TEST_PDFS_CPP
#define __TASMANIAN_TASDREAM_TEST_PDFS_CPP

#include "tasdreamTestPDFs.hpp"

UnscaledUniform1D::UnscaledUniform1D(){}
UnscaledUniform1D::~UnscaledUniform1D(){}
int UnscaledUniform1D::getNumDimensions() const{ return 1; }
void UnscaledUniform1D::evaluate(const std::vector<double> &x, std::vector<double> &y, bool){
    if (y.size() < x.size()) y.resize(x.size());
    for(size_t i=0; i<x.size(); i++) y[i] = 1.0;
}
void UnscaledUniform1D::getDomainBounds(std::vector<bool> &lower, std::vector<bool> &upper){
    if (lower.size() < 1) lower.resize(1);
    if (upper.size() < 1) upper.resize(1);
    lower[0] = true; // Gaussian is unbounded
    upper[0] = true;
}
void UnscaledUniform1D::getDomainBounds(std::vector<double> &lower, std::vector<double> &upper){
    if (lower.size() < 1) lower.resize(1);
    if (upper.size() < 1) upper.resize(1);
    lower[0] = -1.0; // since bounds above give false,
    upper[0] =  1.0;
}
void UnscaledUniform1D::getInitialSample(double y[]){ y[0] = -1.0 + 2.0 * u.getSample01(); }

Beta1D::Beta1D(){
    b = new TasDREAM::BetaPDF(-1.0, 1.0, 2.0, 5.0);
}
Beta1D::~Beta1D(){ delete b; }
int Beta1D::getAPIversion() const{ return 5; }
int Beta1D::getNumDimensions() const{ return 1; }
void Beta1D::evaluate(int num_points, const double x[], double y[], bool useLogForm){
    if (useLogForm){
        for(int i=0; i<num_points; i++) y[i] = b->getDensityLog(x[i]);
    }else{
        for(int i=0; i<num_points; i++) y[i] = b->getDensity(x[i]);
    }
}
void Beta1D::getDomainBounds(bool* lower_bound, bool* upper_bound){ lower_bound[0] = true; upper_bound[0] = true; }
void Beta1D::getDomainBounds(double* lower_bound, double* upper_bound){ lower_bound[0] = -1.0; upper_bound[0] = 1.0; }
void Beta1D::getInitialSample(double y[]){ y[0] = -1.0 + 2.0 * u.getSample01(); }

Gamma1D::Gamma1D(){
    g = new TasDREAM::GammaPDF(-2.0, 9.0, 2.0);
}
Gamma1D::~Gamma1D(){ delete g; }
int Gamma1D::getNumDimensions() const{ return 1; }
void Gamma1D::evaluate(const std::vector<double> &x, std::vector<double> &y, bool useLogForm){
    size_t num_points = x.size();
    if (y.size() < x.size()) y.resize(x.size());
    if (useLogForm){
        for(size_t i=0; i<num_points; i++) y[i] = g->getDensityLog(x[i]);
    }else{
        for(size_t i=0; i<num_points; i++) y[i] = g->getDensity(x[i]);
    }
}
void Gamma1D::getDomainBounds(std::vector<bool> &lower, std::vector<bool> &upper){
    if (lower.size() < 1) lower.resize(1);
    if (upper.size() < 1) upper.resize(1);
    lower[0] = true;
    upper[0] = false;
}
void Gamma1D::getDomainBounds(std::vector<double> &lower, std::vector<double> &upper){
    if (lower.size() < 1) lower.resize(1);
    if (upper.size() < 1) upper.resize(1);
    lower[0] = -2.0;
    upper[0] =  0.0;
}
void Gamma1D::getInitialSample(double y[]){ y[0] = -1.0 + 10.0 * u.getSample01(); }

Gaussian2D::Gaussian2D(){
    g1 = new TasDREAM::TruncatedGaussianPDF(-0.5, 0.1, -1.0, 1.0);
    g2 = new TasDREAM::TruncatedGaussianPDF(-0.5, 0.1, -1.0, 1.0);
}
Gaussian2D::~Gaussian2D(){ delete g1; delete g2; }
int Gaussian2D::getNumDimensions() const{ return 2; }
void Gaussian2D::evaluate(const std::vector<double> &x, std::vector<double> &y, bool useLogForm){
    size_t num_points = x.size() / 2;
    if (y.size() < x.size()) y.resize(num_points);
    if (useLogForm){
        for(size_t i=0; i<num_points; i++) y[i] = g1->getDensityLog(x[2*i]) + g2->getDensityLog(x[2*i+1]);
    }else{
        for(size_t i=0; i<num_points; i++) y[i] = g1->getDensity(x[2*i]) * g2->getDensity(x[2*i+1]);
    }
}
void Gaussian2D::getDomainBounds(std::vector<bool> &lower, std::vector<bool> &upper){
    if (lower.size() < 2) lower.resize(2);
    if (upper.size() < 2) upper.resize(2);
    lower[0] = true; lower[1] = true;
    upper[0] = true; upper[1] = true;
}
void Gaussian2D::getDomainBounds(std::vector<double> &lower, std::vector<double> &upper){
    if (lower.size() < 2) lower.resize(2);
    if (upper.size() < 2) upper.resize(2);
    lower[0] = -1.0; lower[1] = -1.0;
    upper[0] =  1.0; upper[1] =  1.0;
}
void Gaussian2D::getInitialSample(double x[]){
    x[0] = -1.0 + 2.0 * u.getSample01();
    x[1] = -1.0 + 2.0 * u.getSample01();
}

#endif
