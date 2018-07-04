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

#ifndef __TASMANIAN_SPARSE_GRID_GLOBAL_CPP
#define __TASMANIAN_SPARSE_GRID_GLOBAL_CPP

#include "tsgGridGlobal.hpp"

#include "tsgHiddenExternals.hpp"

#ifdef Tasmanian_ENABLE_CUBLAS
#include <cuda_runtime_api.h>
#include <cuda.h>
#endif // Tasmanian_ENABLE_CUBLAS

#include "tsgCudaMacros.hpp"

namespace TasGrid{

GridGlobal::GridGlobal() : num_dimensions(0), num_outputs(0), alpha(0.0), beta(0.0), wrapper(0), tensors(0), active_tensors(0), active_w(0),
    points(0), needed(0), tensor_refs(0), max_levels(0), values(0),
    updated_tensors(0), updated_active_tensors(0), updated_active_w(0), custom(0), accel(0)
{}
GridGlobal::GridGlobal(const GridGlobal &global) : num_dimensions(0), num_outputs(0), alpha(0.0), beta(0.0), wrapper(0), tensors(0), active_tensors(0), active_w(0),
    points(0), needed(0), tensor_refs(0), max_levels(0), values(0),
    updated_tensors(0), updated_active_tensors(0), updated_active_w(0), custom(0), accel(0)
{
    copyGrid(&global);
}
GridGlobal::~GridGlobal(){ reset(true); }

void GridGlobal::write(std::ofstream &ofs) const{
    ofs << std::scientific; ofs.precision(17);
    ofs << num_dimensions << " " << num_outputs << " " << alpha << " " << beta << endl;
    if (num_dimensions > 0){
        ofs << OneDimensionalMeta::getIORuleString(rule) << endl;
        if (rule == rule_customtabulated){
            custom->write(ofs);
        }
        tensors->write(ofs);
        active_tensors->write(ofs);
        ofs << active_w[0];
        for(int i=1; i<active_tensors->getNumIndexes(); i++){
            ofs << " " << active_w[i];
        }
        ofs << endl;
        if (points == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            points->write(ofs);
        }
        if (needed == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            needed->write(ofs);
        }
        ofs << max_levels[0];
        for(int j=1; j<num_dimensions; j++){
            ofs << " " << max_levels[j];
        }
        ofs << endl;
        if (num_outputs > 0) values->write(ofs);
        if (updated_tensors != 0){
            ofs << "1" << endl;
            updated_tensors->write(ofs);
            updated_active_tensors->write(ofs);
            ofs << updated_active_w[0];
            for(int i=1; i<updated_active_tensors->getNumIndexes(); i++){
                ofs << " " << updated_active_w[i];
            }
        }else{
            ofs << "0";
        }
        ofs << endl;
    }
}
void GridGlobal::writeBinary(std::ofstream &ofs) const{
    int num_dim_out[2];
    num_dim_out[0] = num_dimensions;
    num_dim_out[1] = num_outputs;
    ofs.write((char*) num_dim_out, 2*sizeof(int));
    double alpha_beta[2];
    alpha_beta[0] = alpha;
    alpha_beta[1] = beta;
    ofs.write((char*) alpha_beta, 2*sizeof(double));
    if (num_dimensions > 0){
        num_dim_out[0] = OneDimensionalMeta::getIORuleInt(rule);
        ofs.write((char*) num_dim_out, sizeof(int));
        if (rule == rule_customtabulated){
            custom->writeBinary(ofs);
        }
        tensors->writeBinary(ofs);
        active_tensors->writeBinary(ofs);
        ofs.write((char*) active_w, active_tensors->getNumIndexes() * sizeof(int));
        char flag;
        if (points == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            points->writeBinary(ofs);
        }
        if (needed == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            needed->writeBinary(ofs);
        }
        ofs.write((char*) max_levels, num_dimensions * sizeof(int));

        if (num_outputs > 0) values->writeBinary(ofs);
        if (updated_tensors != 0){
            flag = 'y'; ofs.write(&flag, sizeof(char));
            updated_tensors->writeBinary(ofs);
            updated_active_tensors->writeBinary(ofs);
            ofs.write((char*) updated_active_w, updated_active_tensors->getNumIndexes() * sizeof(int));
        }else{
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }
    }
}
void GridGlobal::read(std::ifstream &ifs, std::ostream *logstream){
    reset(true); // true deletes any custom rule
    ifs >> num_dimensions >> num_outputs >> alpha >> beta;
    if (num_dimensions > 0){
        int flag;
        std::string T;
        ifs >> T;
        rule = OneDimensionalMeta::getIORuleString(T.c_str());
        if (rule == rule_customtabulated){
            custom = new CustomTabulated(logstream);
            custom->read(ifs);
        }
        tensors = new IndexSet(num_dimensions);  tensors->read(ifs);
        active_tensors = new IndexSet(num_dimensions);  active_tensors->read(ifs);
        active_w = new int[active_tensors->getNumIndexes()];  for(int i=0; i<active_tensors->getNumIndexes(); i++){ ifs >> active_w[i]; }
        ifs >> flag; if (flag == 1){ points = new IndexSet(num_dimensions); points->read(ifs); }
        ifs >> flag; if (flag == 1){ needed = new IndexSet(num_dimensions); needed->read(ifs); }
        max_levels = new int[num_dimensions];  for(int j=0; j<num_dimensions; j++){ ifs >> max_levels[j]; }
        if (num_outputs > 0){ values = new StorageSet(0, 0); values->read(ifs); }
        ifs >> flag;
        IndexManipulator IM(num_dimensions);
        int oned_max_level;
        if (flag == 1){
            updated_tensors = new IndexSet(num_dimensions);  updated_tensors->read(ifs);
            updated_active_tensors = new IndexSet(num_dimensions);  updated_active_tensors->read(ifs);
            updated_active_w = new int[updated_active_tensors->getNumIndexes()];  for(int i=0; i<updated_active_tensors->getNumIndexes(); i++){ ifs >> updated_active_w[i]; }
            IM.getMaxLevels(updated_tensors, 0, oned_max_level);
        }else{
            oned_max_level = max_levels[0];
            for(int j=1; j<num_dimensions; j++) if (oned_max_level < max_levels[j]) oned_max_level = max_levels[j];
        }
        OneDimensionalMeta meta(custom);
        wrapper = new OneDimensionalWrapper(&meta, oned_max_level, rule, alpha, beta, logstream);

        int nz_weights = active_tensors->getNumIndexes();
        IndexSet *work = (points != 0) ? points : needed;
        if (OneDimensionalMeta::isNonNested(rule)){
            tensor_refs = new int*[nz_weights];
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<nz_weights; i++){
                tensor_refs[i] = IM.referenceGenericPoints(active_tensors->getIndex(i), wrapper, work);
            }
        }else{
            tensor_refs = new int*[nz_weights];
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<nz_weights; i++){
                tensor_refs[i] = IM.referenceNestedPoints(active_tensors->getIndex(i), wrapper, work);
            }
        }
        work = 0;
    }
}

void GridGlobal::readBinary(std::ifstream &ifs, std::ostream *logstream){
    reset(true); // true deletes any custom rule
    int num_dim_out[2];
    ifs.read((char*) num_dim_out, 2*sizeof(int));
    num_dimensions = num_dim_out[0];
    num_outputs = num_dim_out[1];
    double alpha_beta[2];
    ifs.read((char*) alpha_beta, 2*sizeof(double));
    alpha = alpha_beta[0];
    beta = alpha_beta[1];
    if (num_dimensions > 0){
        ifs.read((char*) num_dim_out, sizeof(int));
        rule = OneDimensionalMeta::getIORuleInt(num_dim_out[0]);
        if (rule == rule_customtabulated){
            custom = new CustomTabulated(logstream);
            custom->readBinary(ifs);
        }

        tensors = new IndexSet(num_dimensions);  tensors->readBinary(ifs);
        active_tensors = new IndexSet(num_dimensions);  active_tensors->readBinary(ifs);
        active_w = new int[active_tensors->getNumIndexes()];
        ifs.read((char*) active_w, active_tensors->getNumIndexes() * sizeof(int));

        char flag;
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ points = new IndexSet(num_dimensions); points->readBinary(ifs); }
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ needed = new IndexSet(num_dimensions); needed->readBinary(ifs); }

        max_levels = new int[num_dimensions];
        ifs.read((char*) max_levels, num_dimensions * sizeof(int));

        if (num_outputs > 0){ values = new StorageSet(0, 0); values->readBinary(ifs); }

        ifs.read((char*) &flag, sizeof(char));
        IndexManipulator IM(num_dimensions);
        int oned_max_level;
        if (flag == 'y'){
            updated_tensors = new IndexSet(num_dimensions);  updated_tensors->readBinary(ifs);
            updated_active_tensors = new IndexSet(num_dimensions);  updated_active_tensors->readBinary(ifs);
            updated_active_w = new int[updated_active_tensors->getNumIndexes()];
            ifs.read((char*) updated_active_w, updated_active_tensors->getNumIndexes() * sizeof(int));
            IM.getMaxLevels(updated_tensors, 0, oned_max_level);
        }else{
            oned_max_level = max_levels[0];
            for(int j=1; j<num_dimensions; j++) if (oned_max_level < max_levels[j]) oned_max_level = max_levels[j];
        }
        OneDimensionalMeta meta(custom);
        wrapper = new OneDimensionalWrapper(&meta, oned_max_level, rule, alpha, beta, logstream);

        int nz_weights = active_tensors->getNumIndexes();
        IndexSet *work = (points != 0) ? points : needed;
        if (OneDimensionalMeta::isNonNested(rule)){
            tensor_refs = new int*[nz_weights];
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<nz_weights; i++){
                tensor_refs[i] = IM.referenceGenericPoints(active_tensors->getIndex(i), wrapper, work);
            }
        }else{
            tensor_refs = new int*[nz_weights];
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<nz_weights; i++){
                tensor_refs[i] = IM.referenceNestedPoints(active_tensors->getIndex(i), wrapper, work);
            }
        }
        work = 0;
    }
}

void GridGlobal::reset(bool includeCustom){
    clearAccelerationData();
    if (tensor_refs != 0){ for(int i=0; i<active_tensors->getNumIndexes(); i++){ delete[] tensor_refs[i]; tensor_refs[i] = 0; } delete[] tensor_refs; tensor_refs = 0; }
    if (wrapper != 0){ delete wrapper; wrapper = 0; }
    if (tensors != 0){ delete tensors; tensors = 0; }
    if (active_tensors != 0){ delete active_tensors; active_tensors = 0; }
    if (active_w != 0){ delete[] active_w; active_w = 0; }
    if (points != 0){ delete points; points = 0; }
    if (needed != 0){ delete needed; needed = 0; }
    if (max_levels != 0){ delete[] max_levels; max_levels = 0; }
    if (values != 0){ delete values;  values = 0; }
    if (updated_tensors != 0){ delete updated_tensors; updated_tensors = 0; }
    if (updated_active_tensors != 0){ delete updated_active_tensors; updated_active_tensors = 0; }
    if (updated_active_w != 0){ delete[] updated_active_w; updated_active_w = 0; }
    if (includeCustom && (custom != 0)){ delete custom; custom = 0; }
    num_dimensions = num_outputs = 0;
}

void GridGlobal::clearRefinement(){
    if (needed != 0){ delete needed; needed = 0; }
    if (updated_tensors != 0){ delete updated_tensors; updated_tensors = 0; }
    if (updated_active_tensors != 0){ delete updated_active_tensors; updated_active_tensors = 0; }
    if (updated_active_w != 0){ delete[] updated_active_w; updated_active_w = 0; }
}

void GridGlobal::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, TypeOneDRule crule, const int *anisotropic_weights, double calpha, double cbeta, const char* custom_filename, const int *level_limits){
    if ((crule == rule_customtabulated) && (custom == 0)){
        custom = new CustomTabulated(custom_filename);
    }
    IndexManipulator IM(cnum_dimensions, custom);

    IndexSet *tset = IM.selectTensors(depth, type, anisotropic_weights, crule);
    if (level_limits != 0){
        IndexSet *limited = IM.removeIndexesByLimit(tset, level_limits);
        if (limited != 0){
            delete tset;
            tset = limited;
        }
    }
    setTensors(tset, cnum_outputs, crule, calpha, cbeta);
}
void GridGlobal::copyGrid(const GridGlobal *global){
    if (custom != 0){
        delete custom;  custom = 0;
    }
    if (global->rule == rule_customtabulated) custom = new CustomTabulated(*(global->custom));

    IndexSet *tset = new IndexSet(global->tensors);

    setTensors(tset, global->num_outputs, global->rule, global->alpha, global->beta);

    if ((num_outputs > 0) && (global->points != 0)){ // if there are values inside the source object
        loadNeededPoints(global->values->getValues(0));
    }

    // if there is an active refinement awaiting values
    if (global->updated_tensors != 0){
        updated_tensors = new IndexSet(global->updated_tensors);
        updated_active_tensors = new IndexSet(global->updated_active_tensors);
        updated_active_w = new int[updated_active_tensors->getNumIndexes()];
        std::copy(global->updated_active_w, global->updated_active_w + updated_active_tensors->getNumIndexes(), updated_active_w);

        needed = new IndexSet(global->needed);

        OneDimensionalMeta meta(custom);
        delete wrapper;
        wrapper = new OneDimensionalWrapper(&meta, global->wrapper->getNumLevels(), rule, alpha, beta);
    }
}

void GridGlobal::setTensors(IndexSet* &tset, int cnum_outputs, TypeOneDRule crule, double calpha, double cbeta){
    reset(false);
    num_dimensions = tset->getNumDimensions();
    num_outputs = cnum_outputs;
    rule = crule;
    alpha = calpha;  beta = cbeta;

    tensors = tset;
    tset = 0;

    IndexManipulator IM(num_dimensions, custom);

    OneDimensionalMeta meta(custom);
    max_levels = new int[num_dimensions];
    int max_level; IM.getMaxLevels(tensors, max_levels, max_level);
    wrapper = new OneDimensionalWrapper(&meta, max_level, rule, alpha, beta);


    int* tensors_w = IM.makeTensorWeights(tensors);
    active_tensors = IM.nonzeroSubset(tensors, tensors_w);

    int nz_weights = active_tensors->getNumIndexes();

    active_w = new int[nz_weights];
    int count = 0;
    for(int i=0; i<tensors->getNumIndexes(); i++){ if (tensors_w[i] != 0) active_w[count++] = tensors_w[i]; }

    delete[] tensors_w;

    if (OneDimensionalMeta::isNonNested(rule)){
        needed = IM.generateGenericPoints(active_tensors, wrapper);

        tensor_refs = new int*[nz_weights];
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<nz_weights; i++){
            tensor_refs[i] = IM.referenceGenericPoints(active_tensors->getIndex(i), wrapper, needed);
        }
    }else{
        needed = IM.generateNestedPoints(tensors, wrapper); // nested grids exploit nesting

        tensor_refs = new int*[nz_weights];
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<nz_weights; i++){
            tensor_refs[i] = IM.referenceNestedPoints(active_tensors->getIndex(i), wrapper, needed);
        }
    }

    if (num_outputs == 0){
        points = needed;
        needed = 0;
    }else{
        values = new StorageSet(num_outputs, needed->getNumIndexes());
    }
}

void GridGlobal::updateGrid(int depth, TypeDepth type, const int *anisotropic_weights, const int *level_limits){
    if ((num_outputs == 0) || (points == 0)){
        makeGrid(num_dimensions, num_outputs, depth, type, rule, anisotropic_weights, alpha, beta, 0, level_limits);
    }else{
        clearRefinement();

        IndexManipulator IM(num_dimensions, custom);

        updated_tensors = IM.selectTensors(depth, type, anisotropic_weights, rule);
        if (level_limits != 0){
            IndexSet *limited = IM.removeIndexesByLimit(updated_tensors, level_limits);
            if (limited != 0){
                delete updated_tensors;
                updated_tensors = limited;
            }
        }

        needed = updated_tensors->diffSets(tensors);
        if ((needed != 0) && (needed->getNumIndexes() > 0)){
            delete needed;

            updated_tensors->addIndexSet(tensors); // avoids the case where existing points in tensor are not included in the update

            OneDimensionalMeta meta(custom);
            int max_level; IM.getMaxLevels(updated_tensors, 0, max_level);
            delete wrapper;
            wrapper = new OneDimensionalWrapper(&meta, max_level, rule, alpha, beta);

            int* updates_tensor_w = IM.makeTensorWeights(updated_tensors);
            updated_active_tensors = IM.nonzeroSubset(updated_tensors, updates_tensor_w);

            int nz_weights = updated_active_tensors->getNumIndexes();

            updated_active_w = new int[nz_weights];
            int count = 0;
            for(int i=0; i<updated_tensors->getNumIndexes(); i++){ if (updates_tensor_w[i] != 0) updated_active_w[count++] = updates_tensor_w[i];  }

            IndexSet *new_points = 0;
            if (OneDimensionalMeta::isNonNested(rule)){
                new_points = IM.generateGenericPoints(updated_active_tensors, wrapper);
            }else{
                new_points = IM.generateNestedPoints(updated_tensors, wrapper);
            }

            needed = new_points->diffSets(points);

            delete new_points;
            delete[] updates_tensor_w;
        }else{
            clearRefinement();
        }
    }
}

int GridGlobal::getNumDimensions() const{ return num_dimensions; }
int GridGlobal::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridGlobal::getRule() const{ return rule; }
const char* GridGlobal::getCustomRuleDescription() const{ return (custom != 0) ? custom->getDescription() : "";  }

double GridGlobal::getAlpha() const{ return alpha; }
double GridGlobal::getBeta() const{ return beta; }

int GridGlobal::getNumLoaded() const{ return (((points == 0) || (num_outputs == 0)) ? 0 : points->getNumIndexes()); }
int GridGlobal::getNumNeeded() const{ return ((needed == 0) ? 0 : needed->getNumIndexes()); }
int GridGlobal::getNumPoints() const{ return ((points == 0) ? getNumNeeded() : points->getNumIndexes()); }

double* GridGlobal::getLoadedPoints() const{
    if (points == 0) return 0;
    int num_points = points->getNumIndexes();
    double *x = new double[num_dimensions * num_points];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = wrapper->getNode(p[j]);
        }
    }
    return x;
}
void GridGlobal::getLoadedPoints(double *x) const{
    int num_points = points->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = wrapper->getNode(p[j]);
        }
    }
}
double* GridGlobal::getNeededPoints() const{
    if (needed == 0) return 0;
    int num_points = needed->getNumIndexes();
    if (num_points == 0) return 0;
    double *x = new double[num_dimensions * num_points];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = wrapper->getNode(p[j]);
        }
    }
    return x;
}
void GridGlobal::getNeededPoints(double *x) const{
    int num_points = needed->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = wrapper->getNode(p[j]);
        }
    }
}

double* GridGlobal::getPoints() const{
    return ((points == 0) ? getNeededPoints() : getLoadedPoints());
}
void GridGlobal::getPoints(double *x) const{
    if (points == 0){ getNeededPoints(x); }else{ getLoadedPoints(x); };
}

double* GridGlobal::getQuadratureWeights() const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
    double *weights = new double[num_points];
    getQuadratureWeights(weights);
    return weights;
}
void GridGlobal::getQuadratureWeights(double weights[]) const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
    std::fill(weights, weights + num_points, 0.0);

    int *num_oned_points = new int[num_dimensions];
    for(int n=0; n<active_tensors->getNumIndexes(); n++){
        const int* levels = active_tensors->getIndex(n);
        num_oned_points[0] = wrapper->getNumPoints(levels[0]);
        int num_tensor_points = num_oned_points[0];
        for(int j=1; j<num_dimensions; j++){
            num_oned_points[j] = wrapper->getNumPoints(levels[j]);
            num_tensor_points *= num_oned_points[j];
        }
        double tensor_weight = (double) active_w[n];

        #pragma omp parallel for schedule(static)
        for(int i=0; i<num_tensor_points; i++){
            int t = i;
            double w = 1.0;
            for(int j=num_dimensions-1; j>=0; j--){
                w *= wrapper->getWeight(levels[j], t % num_oned_points[j]);
                t /= num_oned_points[j];
            }
            weights[tensor_refs[n][i]] += tensor_weight * w;
        }
    }
    delete[] num_oned_points;
    work = 0;
}

double* GridGlobal::getInterpolationWeights(const double x[]) const{
    IndexSet *work = (points == 0) ? needed : points;

    int num_points = work->getNumIndexes();
    double *weights = new double[num_points];

    getInterpolationWeights(x, weights);

    return weights;
}
void GridGlobal::getInterpolationWeights(const double x[], double weights[]) const{
    IndexSet *work = (points == 0) ? needed : points;

    CacheLagrange<double> *lcache = new CacheLagrange<double>(num_dimensions, max_levels, wrapper, x);

    int num_points = work->getNumIndexes();
    std::fill(weights, weights + num_points, 0.0);

    int *num_oned_points = new int[num_dimensions];
    for(int n=0; n<active_tensors->getNumIndexes(); n++){
        const int* levels = active_tensors->getIndex(n);
        num_oned_points[0] = wrapper->getNumPoints(levels[0]);
        int num_tensor_points = num_oned_points[0];
        for(int j=1; j<num_dimensions; j++){
            num_oned_points[j] = wrapper->getNumPoints(levels[j]);
            num_tensor_points *= num_oned_points[j];
        }
        double tensor_weight = (double) active_w[n];
        for(int i=0; i<num_tensor_points; i++){
            int t = i;
            double w = 1.0;
            for(int j=num_dimensions-1; j>=0; j--){
                w *= lcache->getLagrange(j, levels[j], t % num_oned_points[j]);
                t /= num_oned_points[j];
            }
            weights[tensor_refs[n][i]] += tensor_weight * w;
        }
    }
    delete[] num_oned_points;
    work = 0;

    delete lcache;
}

void GridGlobal::loadNeededPoints(const double *vals, TypeAcceleration){
    if (accel != 0) accel->resetGPULoadedData();
    if (points == 0){
        values->setValues(vals);
        points = needed;
        needed = 0;
    }else if (needed == 0){ // resetting the points
        values->setValues(vals);
    }else{
        values->addValues(points, needed, vals);
        points->addIndexSet(needed);
        delete needed; needed = 0;
        for(int i=0; i<active_tensors->getNumIndexes(); i++){ delete[] tensor_refs[i]; } delete[] tensor_refs;

        delete tensors; tensors = updated_tensors; updated_tensors = 0;
        delete active_tensors; active_tensors = updated_active_tensors; updated_active_tensors = 0;
        delete[] active_w; active_w = updated_active_w; updated_active_w = 0;

        int nz_weights = active_tensors->getNumIndexes();
        IndexManipulator IM(num_dimensions, custom);
        int m; IM.getMaxLevels(tensors, max_levels, m);
        if (OneDimensionalMeta::isNonNested(rule)){
            tensor_refs = new int*[nz_weights];
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<nz_weights; i++){
                tensor_refs[i] = IM.referenceGenericPoints(active_tensors->getIndex(i), wrapper, points);
            }
        }else{
            tensor_refs = new int*[nz_weights];
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<nz_weights; i++){
                tensor_refs[i] = IM.referenceNestedPoints(active_tensors->getIndex(i), wrapper, points);
            }
        }
    }
}
void GridGlobal::mergeRefinement(){
    if (needed == 0) return; // nothing to do
    int num_all_points = getNumLoaded() + getNumNeeded();
    double *vals = new double[num_all_points * num_outputs];
    std::fill(vals, vals + num_all_points * num_outputs, 0.0);
    values->setValuesPointer(vals, num_all_points);
    if (points == 0){
        points = needed;
        needed = 0;
    }else{
        points->addIndexSet(needed);
        delete needed; needed = 0;
        for(int i=0; i<active_tensors->getNumIndexes(); i++){ delete[] tensor_refs[i]; } delete[] tensor_refs;

        delete tensors; tensors = updated_tensors; updated_tensors = 0;
        delete active_tensors; active_tensors = updated_active_tensors; updated_active_tensors = 0;
        delete[] active_w; active_w = updated_active_w; updated_active_w = 0;

        int nz_weights = active_tensors->getNumIndexes();
        IndexManipulator IM(num_dimensions, custom);
        int m; IM.getMaxLevels(tensors, max_levels, m);
        if (OneDimensionalMeta::isNonNested(rule)){
            tensor_refs = new int*[nz_weights];
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<nz_weights; i++){
                tensor_refs[i] = IM.referenceGenericPoints(active_tensors->getIndex(i), wrapper, points);
            }
        }else{
            tensor_refs = new int*[nz_weights];
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<nz_weights; i++){
                tensor_refs[i] = IM.referenceNestedPoints(active_tensors->getIndex(i), wrapper, points);
            }
        }
    }
}
const double* GridGlobal::getLoadedValues() const{
    if (getNumLoaded() == 0) return 0;
    return values->getValues(0);
}

void GridGlobal::evaluate(const double x[], double y[]) const{
    double *w = getInterpolationWeights(x);
    TasBLAS::setzero(num_outputs, y);
    for(int k=0; k<num_outputs; k++){
        for(int i=0; i<points->getNumIndexes(); i++){
            const double *v = values->getValues(i);
            y[k] += w[i] * v[k];
        }
    }
    delete[] w;
}

void GridGlobal::evaluateFastCPUblas(const double x[], double y[]) const{
    #ifdef Tasmanian_ENABLE_BLAS
    double *w = getInterpolationWeights(x);
    TasBLAS::dgemv(num_outputs, points->getNumIndexes(), values->getValues(0), w, y);
    delete[] w;
    #else
    evaluate(x, y);
    #endif // Tasmanian_ENABLE_BLAS
}
#ifdef Tasmanian_ENABLE_CUBLAS
void GridGlobal::evaluateFastGPUcublas(const double x[], double y[], std::ostream *os) const{
    makeCheckAccelerationData(accel_gpu_cublas, os);

    AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
    double *weights = getInterpolationWeights(x);

    gpu->cublasDGEMV(num_outputs, points->getNumIndexes(), weights, y);

    delete[] weights;
}
#else
void GridGlobal::evaluateFastGPUcublas(const double x[], double y[], std::ostream *) const{ evaluateFastCPUblas(x, y); }
#endif // Tasmanian_ENABLE_CUBLAS
void GridGlobal::evaluateFastGPUcuda(const double x[], double y[], std::ostream *os) const{
    evaluateFastGPUcublas(x, y, os);
}
void GridGlobal::evaluateFastGPUmagma(const double x[], double y[], std::ostream *os) const{
    evaluateFastGPUcublas(x, y, os);
}

void GridGlobal::evaluateBatch(const double x[], int num_x, double y[]) const{
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        evaluate(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(y[((size_t) i) * ((size_t) num_outputs)]));
    }
}

void GridGlobal::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    #ifdef Tasmanian_ENABLE_BLAS
    int num_points = points->getNumIndexes();
    double *weights = new double[num_points * num_x];
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        getInterpolationWeights(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(weights[((size_t) i) * ((size_t) num_points)]));
    }

    TasBLAS::dgemm(num_outputs, num_x, num_points, 1.0, values->getValues(0), weights, 0.0, y);

    delete[] weights;
    #else
    evaluateBatch(x, num_x, y);
    #endif // Tasmanian_ENABLE_BLAS
}
#ifdef Tasmanian_ENABLE_CUBLAS
void GridGlobal::evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *os) const{
    int num_points = points->getNumIndexes();
    makeCheckAccelerationData(accel_gpu_cublas, os);

    AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
    double *weights = new double[((size_t) num_points) * ((size_t) num_x)];
    evaluateHierarchicalFunctions(x, num_x, weights);

    double *gpu_weights = TasCUDA::cudaSend(((size_t) num_points) * ((size_t) num_x), weights, os);
    double *gpu_result = TasCUDA::cudaNew<double>(((size_t) num_outputs) * ((size_t) num_x), os);

    gpu->cublasDGEMM(num_outputs, num_x, num_points, gpu_weights, gpu_result);

    TasCUDA::cudaRecv<double>(num_outputs * num_x, gpu_result, y, os);

    TasCUDA::cudaDel<double>(gpu_result, os);
    TasCUDA::cudaDel<double>(gpu_weights, os);

    delete[] weights;
}
#else
void GridGlobal::evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *) const{ evaluateBatchCPUblas(x, num_x, y); }
#endif // Tasmanian_ENABLE_CUBLAS
void GridGlobal::evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream *os) const{
    evaluateBatchGPUcublas(x, num_x, y, os);
}
void GridGlobal::evaluateBatchGPUmagma(const double x[], int num_x, double y[], std::ostream *os) const{
    evaluateBatchGPUcublas(x, num_x, y, os);
}
#if defined(Tasmanian_ENABLE_CUBLAS) || defined(Tasmanian_ENABLE_CUDA)
void GridGlobal::makeCheckAccelerationData(TypeAcceleration acc, std::ostream *os) const{
    if (AccelerationMeta::isAccTypeFullMemoryGPU(acc)){
        if ((accel != 0) && (!accel->isCompatible(acc))){
            delete accel;
            accel = 0;
        }
        if (accel == 0){ accel = (BaseAccelerationData*) (new AccelerationDataGPUFull()); }
        AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
        gpu->setLogStream(os);
        double *gpu_values = gpu->getGPUValues();
        if (gpu_values == 0) gpu->loadGPUValues(((size_t) points->getNumIndexes()) * ((size_t) values->getNumOutputs()), values->getValues(0));
    }
}
#else
void GridGlobal::makeCheckAccelerationData(TypeAcceleration, std::ostream *) const{}
#endif // Tasmanian_ENABLE_CUBLAS || Tasmanian_ENABLE_CUDA

void GridGlobal::integrate(double q[], double *conformal_correction) const{
    double *w = getQuadratureWeights();
    if (conformal_correction != 0) for(int i=0; i<points->getNumIndexes(); i++) w[i] *= conformal_correction[i];
    std::fill(q, q+num_outputs, 0.0);
    #pragma omp parallel for schedule(static)
    for(int k=0; k<num_outputs; k++){
        for(int i=0; i<points->getNumIndexes(); i++){
            const double *v = values->getValues(i);
            q[k] += w[i] * v[k];
        }
    }
    delete[] w;
}

void GridGlobal::evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const{
    int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        getInterpolationWeights(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(y[((size_t) i) * ((size_t) num_points)]));
    }
}

double* GridGlobal::computeSurpluses(int output, bool normalize) const{
    int num_points = points->getNumIndexes();
    double *surp = new double[num_points];

    if (OneDimensionalMeta::isSequence(rule)){
        double max_surp = 0.0;
        for(int i=0; i<num_points; i++){
            const double* v = values->getValues(i);
            surp[i] = v[output];
            if (fabs(surp[i]) > max_surp) max_surp = fabs(surp[i]);
        }

        IndexManipulator IM(num_dimensions);
        int *level = IM.computeLevels(points);
        int top_level = level[0];  for(int i=1; i<num_points; i++){ if (top_level < level[i]) top_level = level[i];  }
        int top_1d = 0; const int *id = points->getIndex(0); for(int i=0; i<num_points*num_dimensions; i++) if (top_1d < id[i]) top_1d = id[i];

        int *parents = IM.computeDAGup(points);

        const double* nodes = wrapper->getNodes(0);
        double *coeff = new double[top_1d+1];
        coeff[0] = 1.0;
        for(int i=1; i<=top_1d; i++){
            coeff[i] = 1.0;
            for(int j=0; j<i; j++) coeff[i] *= (nodes[i] - nodes[j]);
        }

        for(int l=1; l<=top_level; l++){
            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<num_points; i++){
                if (level[i] == l){
                    const int* p = points->getIndex(i);

                    int *monkey_count = new int[top_level + 1];
                    int *monkey_tail = new int[top_level + 1];
                    bool *used = new bool[num_points];
                    std::fill(used, used + num_points, false);

                    int current = 0;

                    monkey_count[0] = 0;
                    monkey_tail[0] = i;

                    while(monkey_count[0] < num_dimensions){
                        if (monkey_count[current] < num_dimensions){
                            int branch = parents[monkey_tail[current] * num_dimensions + monkey_count[current]];
                            if ((branch == -1) || (used[branch])){
                                monkey_count[current]++;
                            }else{
                                double basis_value = 1.0;
                                const int* f = points->getIndex(branch);
                                for(int j=0; j<num_dimensions; j++){
                                    double x = nodes[p[j]];
                                    double w = 1.0;
                                    for(int k=0; k<f[j]; k++){
                                        w *= (x - nodes[k]);
                                    }
                                    basis_value *= w / coeff[f[j]];
                                }
                                surp[i] -= basis_value * surp[branch];

                                used[branch] = true;

                                monkey_count[++current] = 0;
                                monkey_tail[current] = branch;
                            }
                        }else{
                            monkey_count[--current]++;
                        }
                    }

                    delete[] used;
                    delete[] monkey_tail;
                    delete[] monkey_count;
                }
            }
        }

        delete[] coeff;
        delete[] parents;
        delete[] level;

        if (normalize){
            #pragma omp parallel for schedule(static)
            for(int i=0; i<num_points; i++) surp[i] /= max_surp;
        }
    }else{
        IndexManipulator IM(num_dimensions, custom);
        IndexSet* polynomial_set = IM.getPolynomialSpace(active_tensors, rule, true);

        GridGlobal *gg = new GridGlobal();
        IndexSet *quadrature_tensors = IM.selectTensors(polynomial_set, true, rule_gausspatterson);

        int *max_lvl = new int[num_dimensions];
        int ml; IM.getMaxLevels(quadrature_tensors, max_lvl, ml);
        delete[] max_lvl;

        if (ml < TableGaussPatterson::getNumLevels()-1){
            gg->setTensors(quadrature_tensors, 0, rule_gausspatterson, 0.0, 0.0);
        }else{
            delete quadrature_tensors;
            quadrature_tensors = IM.selectTensors(polynomial_set, true, rule_clenshawcurtis);
            gg->setTensors(quadrature_tensors, 0, rule_clenshawcurtis, 0.0, 0.0);
        }
        delete polynomial_set;

        int qn = gg->getNumPoints();
        double *w = gg->getQuadratureWeights();
        double *x = gg->getPoints();
        double *I = new double[qn];
        delete gg;
        #pragma omp parallel for schedule(static)
        for(int i=0; i<qn; i++){
            double *y = new double[num_outputs];
            evaluate(&(x[i*num_dimensions]), y);
            I[i] = w[i] * y[output];
            delete[] y;
        }

        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<num_points; i++){
            // for each surp, do the quadrature
            const int* p = points->getIndex(i);
            double c = 0.0;
            for(int k=0; k<qn; k++){
                double v = legendre(p[0], x[k*num_dimensions]);
                for(int j=1; j<num_dimensions; j++){
                    v *= legendre(p[j], x[k*num_dimensions+j]);
                }
                c += v * I[k];
            }
            double nrm = sqrt((double) p[0] + 0.5);
            for(int j=1; j<num_dimensions; j++){
                nrm *= sqrt((double) p[j] + 0.5);
            }
            surp[i] = c * nrm;
        }

        delete[] I;
        delete[] x;
        delete[] w;
    }

    return surp;
}

int* GridGlobal::estimateAnisotropicCoefficients(TypeDepth type, int output) const{
    double tol = 1000.0 * TSG_NUM_TOL;
    double *surp = computeSurpluses(output, false);

    int num_points = points->getNumIndexes();

    int n = 0, m;
    for(int j=0; j<num_points; j++){
        surp[j] = fabs(surp[j]);
        if (surp[j] > tol) n++;
    }

    double *A, *b;

    if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
        m = 2*num_dimensions + 1;
        A = new double[n * m];
        b = new double[n];

        int count = 0;
        for(int c=0; c<num_points; c++){
            const int *indx = points->getIndex(c);
            if (surp[c] > tol){
                for(int j=0; j<num_dimensions; j++){
                    A[j*n + count] = ((double) indx[j]);
                }
                for(int j=0; j<num_dimensions; j++){
                    A[(num_dimensions + j)*n + count] = log((double) (indx[j] + 1));
                }
                A[2*num_dimensions*n + count] = 1.0;
                b[count++] = -log(surp[c]);
            }
        }
    }else{
        m = num_dimensions + 1;
        A = new double[n * m];
        b = new double[n];

        int count = 0;
        for(int c=0; c<num_points; c++){
            const int *indx = points->getIndex(c);
            if (surp[c] > tol){
                for(int j=0; j<num_dimensions; j++){
                    A[j*n + count] = ((double) indx[j]);
                }
                A[num_dimensions*n + count] = 1.0;
                b[count++] = - log(surp[c]);
            }
        }
    }
    delete[] surp;

    double *x = new double[m];
    TasmanianDenseSolver::solveLeastSquares(n, m, A, b, 1.E-5, x);

    int *weights = new int[m--];
    for(int j=0; j<m; j++){
        weights[j] = (int)(x[j] * 1000.0 + 0.5);
    }

    int min_weight = weights[0];  for(int j=1; j<num_dimensions; j++){ if (min_weight < weights[j]) min_weight = weights[j];  } // start by finding the largest weight
    if (min_weight < 0){ // all directions are diverging, default to isotropic total degree
        for(int j=0; j<num_dimensions; j++)  weights[j] = 1;
        if (m == 2*num_dimensions) for(int j=num_dimensions; j<2*num_dimensions; j++)  weights[j] = 0;
    }else{
        for(int j=0; j<num_dimensions; j++)  if((weights[j]>0) && (weights[j]<min_weight)) min_weight = weights[j];  // find the smallest positive weight
        for(int j=0; j<num_dimensions; j++){
            if (weights[j] <= 0){
                weights[j] = min_weight;
                if (m == 2*num_dimensions){
                    if (std::abs(weights[num_dimensions + j]) > weights[j]){
                        weights[num_dimensions + j] = (weights[num_dimensions + j] > 0.0) ? weights[j] : -weights[j];
                    }
                }
            }
        }
    }

    delete[] A;
    delete[] b;
    delete[] x;

    return weights;
}

void GridGlobal::setAnisotropicRefinement(TypeDepth type, int min_growth, int output, const int *level_limits){
    clearRefinement();
    int *weights = estimateAnisotropicCoefficients(type, output); //for(int i=0; i<2*num_dimensions; i++) cout << weights[i] << "  "; cout << endl;

    IndexManipulator IM(num_dimensions);
    int level = IM.getMinChildLevel(tensors, type, weights, rule);

    updated_tensors = IM.selectTensors(level, type, weights, rule);
    needed = updated_tensors->diffSets(tensors); // this exploits the 1-1 correspondence between points and tensors for sequence rules (see the correction below)

    while((needed == 0) || (needed->getNumIndexes() < min_growth)){ // for CC min_growth is lots of points
        delete updated_tensors;
        if (needed != 0) delete needed;
        updated_tensors = IM.selectTensors(++level, type, weights, rule);
        needed = updated_tensors->diffSets(tensors);
    }
    delete[] weights;

    if (level_limits != 0){
        IndexSet *limited = IM.removeIndexesByLimit(updated_tensors, level_limits);
        if (limited != 0){
            delete updated_tensors;
            updated_tensors = limited;
        }
    }
    updated_tensors->addIndexSet(tensors); // avoids the case where existing points in tensor are not included in the update

    OneDimensionalMeta meta(custom);
    int max_level; IM.getMaxLevels(updated_tensors, 0, max_level);
    delete wrapper;
    wrapper = new OneDimensionalWrapper(&meta, max_level, rule, alpha, beta);

    int* updates_tensor_w = IM.makeTensorWeights(updated_tensors);
    updated_active_tensors = IM.nonzeroSubset(updated_tensors, updates_tensor_w);

    int nz_weights = updated_active_tensors->getNumIndexes();

    updated_active_w = new int[nz_weights];
    int count = 0;
    for(int i=0; i<updated_tensors->getNumIndexes(); i++){ if (updates_tensor_w[i] != 0) updated_active_w[count++] = updates_tensor_w[i];  }

    if (!(OneDimensionalMeta::isSequence(rule))){ // non-sequence rules need to recompute needed
        delete needed;

        IndexSet *new_points = IM.generateNestedPoints(updated_tensors, wrapper);

        needed = new_points->diffSets(points);

        delete new_points;
    }
    delete[] updates_tensor_w;
}

void GridGlobal::setSurplusRefinement(double tolerance, int output, const int *level_limits){
    clearRefinement();
    double *surp = computeSurpluses(output, true);

    int n = points->getNumIndexes();
    bool *flagged = new bool[n];

    for(int i=0; i<n; i++){
        flagged[i] = (fabs(surp[i]) > tolerance);
    }

    IndexManipulator IM(num_dimensions);
    IndexSet *kids = IM.selectFlaggedChildren(points, flagged, level_limits);

    if ((kids != 0) && (kids->getNumIndexes() > 0)){
        kids->addIndexSet(points);

        updated_tensors = IM.getLowerCompletion(kids);
        if (updated_tensors == 0){
            updated_tensors = kids;
            kids = 0;
        }else{
            updated_tensors->addIndexSet(kids);
            delete kids;
        }

        OneDimensionalMeta meta(custom);
        int max_level; IM.getMaxLevels(updated_tensors, 0, max_level);
        delete wrapper;
        wrapper = new OneDimensionalWrapper(&meta, max_level, rule, alpha, beta);

        int* updates_tensor_w = IM.makeTensorWeights(updated_tensors);
        updated_active_tensors = IM.nonzeroSubset(updated_tensors, updates_tensor_w);

        int nz_weights = updated_active_tensors->getNumIndexes();

        updated_active_w = new int[nz_weights];
        int count = 0;
        for(int i=0; i<updated_tensors->getNumIndexes(); i++){ if (updates_tensor_w[i] != 0) updated_active_w[count++] = updates_tensor_w[i];  }

        delete[] updates_tensor_w;

        needed = updated_tensors->diffSets(tensors);
    }

    delete[] flagged;
    delete[] surp;
}
void GridGlobal::setHierarchicalCoefficients(const double c[], TypeAcceleration acc, std::ostream*){
    if (accel != 0) accel->resetGPULoadedData();
    if (points != 0) clearRefinement();
    loadNeededPoints(c, acc);
}

double GridGlobal::legendre(int n, double x){
    if (n == 0) return 1.0;
    if (n == 1) return x;

    double lm = 1, l = x;
    for(int i=2; i<=n; i++){
        double lp = ((2*i - 1) * x * l) / ((double) i) - ((i - 1) * lm) / ((double) i);
        lm = l; l = lp;
    }
    return l;
}

void GridGlobal::clearAccelerationData(){
    if (accel != 0){
        delete accel;
        accel = 0;
    }
}

void GridGlobal::getPolynomialSpace(bool interpolation, int &n, int* &poly) const{
    IndexManipulator IM(num_dimensions, custom);
    IndexSet* set = IM.getPolynomialSpace(active_tensors, rule, interpolation);

    n = set->getNumIndexes();

    poly = new int[n * num_dimensions];
    const int* p = set->getIndex(0);

    std::copy(p, p + n * num_dimensions, poly);

    delete set;
}
const int* GridGlobal::getPointIndexes() const{
    return ((points == 0) ? needed->getIndex(0) : points->getIndex(0));
}

}

#endif
