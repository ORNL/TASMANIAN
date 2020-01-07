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

#ifndef __TASMANIAN_DREAM_STATE_CPP
#define __TASMANIAN_DREAM_STATE_CPP

#include "tsgDreamState.hpp"

namespace TasDREAM{

TasmanianDREAM::TasmanianDREAM() : num_chains(0), num_dimensions(0), init_state(false), init_values(false), accepted(0){}

TasmanianDREAM::TasmanianDREAM(int cnum_chains, int cnum_dimensions) :
num_chains(cnum_chains), num_dimensions(cnum_dimensions), init_state(false), init_values(false), accepted(0){
    if (cnum_chains < 1) throw std::invalid_argument("ERROR: num_chains must be positive");
    if (cnum_dimensions < 1) throw std::invalid_argument("ERROR: num_dimensions must be positive");
}
TasmanianDREAM::TasmanianDREAM(int cnum_chains, const TasGrid::TasmanianSparseGrid &grid) :
num_chains(cnum_chains), num_dimensions(grid.getNumDimensions()), init_state(false), init_values(false), accepted(0){
    if (cnum_chains < 1) throw std::invalid_argument("ERROR: num_chains must be positive");
    if (grid.getNumDimensions() < 1) throw std::invalid_argument("ERROR: num_dimensions must be positive");
}
TasmanianDREAM::~TasmanianDREAM(){}

void TasmanianDREAM::setState(const std::vector<double> &new_state){
    if (new_state.size() != num_chains * num_dimensions) throw std::runtime_error("ERROR: new state has incorrect dimension, must be num_chains times num_dimensions.");
    state = new_state;
    init_state = true;
    init_values = false;
}
void TasmanianDREAM::setState(std::function<void(double *)> update_state){
    state.resize(num_chains * num_dimensions);
    auto istate = state.begin();
    for(size_t i=0; i<num_chains; i++){
        update_state(&*istate);
        std::advance(istate, num_dimensions);
    }
    init_state = true;
    init_values = false;
}

void TasmanianDREAM::setPDFvalues(const std::vector<double> &new_values){
    if (new_values.size() != num_chains) throw std::runtime_error("ERROR: new_values has incorrect dimension, must match num_chains.");
    pdf_values = new_values;
    init_values = true;
}
void TasmanianDREAM::setPDFvalues(std::function<void(const std::vector<double> &state, std::vector<double> &values)> probability_distribution){
    if (!init_state) throw std::runtime_error("ERROR: calling setPDFvalues() with a function requires that the state is set first.");
    pdf_values.resize(num_chains);
    probability_distribution(state, pdf_values);
    init_values = true;
}

void TasmanianDREAM::clearPDFvalues(){
    pdf_values = std::vector<double>();
    init_values = false;
}

void TasmanianDREAM::getIJKdelta(size_t i, size_t j, size_t k, double w, std::vector<double> &x) const{
    std::copy_n(state.begin() + i * num_dimensions, num_dimensions, x.data());
    if (w != 0.0){
        auto ik = state.begin() + k * num_dimensions;
        auto ij = state.begin() + j * num_dimensions;
        for(auto &xv : x) xv += w * (*ik++ - *ij++);
    }
}

void TasmanianDREAM::expandHistory(int num_snapshots){
    history.reserve(history.size() + num_snapshots * num_dimensions * num_chains);
    pdf_history.reserve(pdf_history.size() + num_snapshots * num_chains);
}

void TasmanianDREAM::saveStateHistory(size_t num_accepted){
    history.insert(history.end(), state.begin(), state.end());
    pdf_history.insert(pdf_history.end(), pdf_values.begin(), pdf_values.end());
    accepted += num_accepted;
}

void TasmanianDREAM::getHistoryMeanVariance(std::vector<double> &mean, std::vector<double> &var) const{
    mean.resize(num_dimensions);
    var.resize(num_dimensions);
    std::fill(mean.begin(), mean.end(), 0.0);
    std::fill(var.begin(), var.end(), 0.0);
    if (num_dimensions == 0) return;

    auto ih = history.begin();
    while(ih != history.end()){
        auto im = mean.begin();
        auto iv = var.begin();
        for(size_t i=0; i<num_dimensions; i++){
            *iv++ += *ih * *ih;
            *im++ += *ih++;
        }
    }

    double n = (double) (history.size() / num_dimensions);
    for(auto &m : mean) m /= n;

    auto im = mean.begin();
    for(auto &v : var){
        v /= n;
        v -= *im * *im;
        im++;
    }
}

void TasmanianDREAM::getApproximateMode(std::vector<double> &mode) const{
    auto imax = std::max_element(pdf_history.begin(), pdf_history.end());
    mode.resize(num_dimensions);
    std::copy_n(history.begin() + std::distance(pdf_history.begin(), imax) * num_dimensions, num_dimensions, mode.data());
}

void TasmanianDREAM::clearHistory(){
    history = std::vector<double>();
    pdf_history = std::vector<double>();
    accepted = 0;
}

extern "C"{ // for python purposes
void* tsgMakeDreamState(int num_chains, int num_dimensions){
    return (void*) new TasmanianDREAM(num_chains, num_dimensions);
}
void tsgDeleteDreamState(void* state){
    delete reinterpret_cast<TasmanianDREAM*>(state);
}
int tsgDreamStateGetDims(void *state){
    return reinterpret_cast<TasmanianDREAM*>(state)->getNumDimensions();
}
int tsgDreamStateGetChains(void *state){
    return reinterpret_cast<TasmanianDREAM*>(state)->getNumChains();
}
int tsgDreamStateGetNumHistory(void *state){
    return (int) reinterpret_cast<TasmanianDREAM*>(state)->getNumHistory();
}
void tsgDreamStateSet(void *state, double const x[]){
    std::vector<double> vx(x, x + Utils::size_mult(reinterpret_cast<TasmanianDREAM*>(state)->getNumChains(),
                                                   reinterpret_cast<TasmanianDREAM*>(state)->getNumDimensions()));
    reinterpret_cast<TasmanianDREAM*>(state)->setState(vx);
}

void tsgDreamStateGetHistory(void *state, double hist[]){
    auto h = reinterpret_cast<TasmanianDREAM*>(state)->getHistory();
    std::copy(h.begin(), h.end(), hist);
}
void tsgDreamStateGetHistoryPDF(void *state, double histpdf[]){
    auto h = reinterpret_cast<TasmanianDREAM*>(state)->getHistoryPDF();
    std::copy(h.begin(), h.end(), histpdf);
}
void tsgDreamStateGetMeanVar(void *state, double mean[], double variance[]){
    std::vector<double> mn, var;
    reinterpret_cast<TasmanianDREAM*>(state)->getHistoryMeanVariance(mn, var);
    std::copy(mn.begin(), mn.end(), mean);
    std::copy(var.begin(), var.end(), variance);
}
void tsgDreamStateGetMode(void *state, double mode[]){
    auto m = reinterpret_cast<TasmanianDREAM*>(state)->getApproximateMode();
    std::copy(m.begin(), m.end(), mode);
}
double tsgDreamStateGetRate(void *state){
    return reinterpret_cast<TasmanianDREAM*>(state)->getAcceptanceRate();
}

}

}

#endif
