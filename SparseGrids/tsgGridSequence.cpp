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

#ifndef __TASMANIAN_SPARSE_GRID_GLOBAL_NESTED_CPP
#define __TASMANIAN_SPARSE_GRID_GLOBAL_NESTED_CPP

#include "tsgGridSequence.hpp"

#include "tsgHiddenExternals.hpp"

#include "tsgCudaMacros.hpp"

namespace TasGrid{

GridSequence::GridSequence() : num_dimensions(0), num_outputs(0){}
GridSequence::~GridSequence(){}

void GridSequence::write(std::ofstream &ofs) const{
    using std::endl;

    ofs << std::scientific; ofs.precision(17);
    ofs << num_dimensions << " " << num_outputs << endl;
    if (num_dimensions > 0){
        ofs << OneDimensionalMeta::getIORuleString(rule) << endl;
        if (points.empty()){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            points.write(ofs);
        }
        if (needed.empty()){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            needed.write(ofs);
        }
        if (surpluses.empty()){
            ofs << "0";
        }else{
            ofs << "1";
            for(auto s : surpluses) ofs << " " << s;
        }
        ofs << endl;
        if (num_outputs > 0) values.write(ofs);
    }
}
void GridSequence::writeBinary(std::ofstream &ofs) const{
    int num_dim_out[2];
    num_dim_out[0] = num_dimensions;
    num_dim_out[1] = num_outputs;
    ofs.write((char*) num_dim_out, 2*sizeof(int));
    if (num_dimensions > 0){
        num_dim_out[0] = OneDimensionalMeta::getIORuleInt(rule);
        ofs.write((char*) num_dim_out, sizeof(int));

        char flag;
        if (points.empty()){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            points.writeBinary(ofs);
        }
        if (needed.empty()){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            needed.writeBinary(ofs);
        }
        if (surpluses.empty()){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            ofs.write((char*) surpluses.data(), surpluses.size() * sizeof(double));
        }
        if (num_outputs > 0) values.writeBinary(ofs);
    }
}
void GridSequence::read(std::ifstream &ifs){
    reset();
    ifs >> num_dimensions >> num_outputs;
    if (num_dimensions > 0){
        int flag;
        std::string T;
        ifs >> T;
        rule = OneDimensionalMeta::getIORuleString(T.c_str());
        ifs >> flag;  if (flag == 1) points.read(ifs);
        ifs >> flag;  if (flag == 1) needed.read(ifs);
        ifs >> flag;
        if (flag == 1){
            surpluses.resize(((size_t) num_outputs) * ((size_t) points.getNumIndexes()));
            for(auto &s : surpluses) ifs >> s;
        }
        if (num_outputs > 0) values.read(ifs);

        prepareSequence();
    }
}
void GridSequence::readBinary(std::ifstream &ifs){
    reset();
    int num_dim_out[2];
    ifs.read((char*) num_dim_out, 2*sizeof(int));
    num_dimensions = num_dim_out[0];
    num_outputs = num_dim_out[1];
    if (num_dimensions > 0){
        ifs.read((char*) num_dim_out, sizeof(int));
        rule = OneDimensionalMeta::getIORuleInt(num_dim_out[0]);

        char flag;
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y') points.readBinary(ifs);
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y') needed.readBinary(ifs);

        ifs.read((char*) &flag, sizeof(char));
        if (flag == 'y'){
            surpluses.resize(((size_t) num_outputs) * ((size_t) points.getNumIndexes()));
            ifs.read((char*) surpluses.data(), surpluses.size() * sizeof(double));
        }

        if (num_outputs > 0) values.readBinary(ifs);

        prepareSequence();
    }
}

void GridSequence::reset(){
    clearAccelerationData();
    points = MultiIndexSet();
    needed = MultiIndexSet();
    values = StorageSet();
    nodes.clear();
    coeff.clear();
    surpluses.resize(0);
    surpluses.shrink_to_fit();
}
void GridSequence::clearRefinement(){ needed = MultiIndexSet(); }

void GridSequence::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, TypeOneDRule crule, const std::vector<int> &anisotropic_weights, const std::vector<int> &level_limits){

    MultiIndexSet pset(cnum_dimensions);
    if ((type == type_qptotal) || (type == type_qptensor) || (type == type_qpcurved) || (type == type_qphyperbolic)){
        MultiIndexManipulations::selectTensors(depth, type, [&](int i) -> long long{ return OneDimensionalMeta::getQExact(i, crule); },
                                               anisotropic_weights, pset);
    }else{
        MultiIndexManipulations::selectTensors(depth, type, [&](int i) -> long long{ return i; }, anisotropic_weights, pset);
    }

    if (!level_limits.empty()) MultiIndexManipulations::removeIndexesByLimit(level_limits, pset);

    setPoints(pset, cnum_outputs, crule);
}
void GridSequence::copyGrid(const GridSequence *seq){
    MultiIndexSet pset;
    if (seq->points.empty()){
        pset = seq->needed;
    }else{
        pset = seq->points;
    }
    setPoints(pset, seq->num_outputs, seq->rule);

    if ((num_outputs > 0) && (!seq->points.empty())){ // if there are values inside the source object
        loadNeededPoints(seq->values.getValues(0));
    }

    if ((!seq->points.empty()) && (!seq->needed.empty())){ // there is a refinement
        needed = seq->needed;
        prepareSequence();
    }
    surpluses = seq->surpluses;
}

void GridSequence::setPoints(MultiIndexSet &pset, int cnum_outputs, TypeOneDRule crule){
    reset();
    num_dimensions = pset.getNumDimensions();
    num_outputs = cnum_outputs;
    rule = crule;

    if (num_outputs == 0){
        points = std::move(pset);
    }else{
        needed = std::move(pset);
        values.resize(num_outputs, needed.getNumIndexes());
    }
    prepareSequence();
}

void GridSequence::updateGrid(int depth, TypeDepth type, const std::vector<int> &anisotropic_weights, const std::vector<int> &level_limits){
    MultiIndexSet pset(num_dimensions);
    if ((type == type_qptotal) || (type == type_qptensor) || (type == type_qpcurved) || (type == type_qphyperbolic)){
        MultiIndexManipulations::selectTensors(depth, type, [&](int i) -> long long{ return OneDimensionalMeta::getQExact(i, rule); },
                                               anisotropic_weights, pset);
    }else{
        MultiIndexManipulations::selectTensors(depth, type, [&](int i) -> long long{ return i; }, anisotropic_weights, pset);
    }

    if (!level_limits.empty()) MultiIndexManipulations::removeIndexesByLimit(level_limits, pset);

    updateGrid(pset);
}

void GridSequence::updateGrid(MultiIndexSet &update){
    clearRefinement();
    if ((num_outputs == 0) || (points.empty())){
        setPoints(update, num_outputs, rule);
    }else{
        update.addSortedInsexes(*points.getVector());
        update.diffSets(points, needed);

        if (!needed.empty()) prepareSequence();
    }
}

int GridSequence::getNumDimensions() const{ return num_dimensions; }
int GridSequence::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridSequence::getRule() const{ return rule; }

int GridSequence::getNumLoaded() const{ return (((points.empty()) || (num_outputs == 0)) ? 0 : points.getNumIndexes()); }
int GridSequence::getNumNeeded() const{ return needed.getNumIndexes(); }
int GridSequence::getNumPoints() const{ return ((points.empty()) ? needed.getNumIndexes() : points.getNumIndexes()); }

void GridSequence::getLoadedPoints(double *x) const{
    int num_points = points.getNumIndexes();
    Data2D<double> split;
    split.load(num_dimensions, num_points, x);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points.getIndex(i);
        double *xx = split.getStrip(i);
        for(int j=0; j<num_dimensions; j++){
            xx[j] = nodes[p[j]];
        }
    }
}
void GridSequence::getNeededPoints(double *x) const{
    int num_points = needed.getNumIndexes();
    Data2D<double> split;
    split.load(num_dimensions, num_points, x);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed.getIndex(i);
        double *xx = split.getStrip(i);
        for(int j=0; j<num_dimensions; j++){
            xx[j] = nodes[p[j]];
        }
    }
}
void GridSequence::getPoints(double *x) const{
    if (points.empty()){ getNeededPoints(x); }else{ getLoadedPoints(x); }
}

void GridSequence::getQuadratureWeights(double *weights) const{
    const MultiIndexSet& work = (points.empty()) ? needed : points;
    std::vector<double> integ;
    cacheBasisIntegrals(integ);
    int n = work.getNumIndexes();
    for(int i=0; i<n; i++){
        const int* p = work.getIndex(i);
        weights[i] = integ[p[0]];
        for(int j=1; j<num_dimensions; j++){
            weights[i] *= integ[p[j]];
        }
    }

    applyTransformationTransposed(weights);
}

void GridSequence::getInterpolationWeights(const double x[], double *weights) const{
    std::vector<std::vector<double>> cache;
    cacheBasisValues<double>(x, cache);
    const MultiIndexSet& work = (points.empty()) ? needed : points;
    int n = work.getNumIndexes();
    weights[0] = 1.0;
    for(int i=1; i<n; i++){
        const int* p = work.getIndex(i);
        weights[i] = cache[0][p[0]];
        for(int j=1; j<num_dimensions; j++){
            weights[i] *= cache[j][p[j]];
        }
    }
    applyTransformationTransposed(weights);
}

void GridSequence::loadNeededPoints(const double *vals, TypeAcceleration){
    #ifdef Tasmanian_ENABLE_CUDA
    cuda_surpluses.clear();
    clearCudaNodes();
    #endif
    if (points.empty()){
        values.setValues(vals);
        points = std::move(needed);
        needed = MultiIndexSet();
    }else if (needed.empty()){
        values.setValues(vals);
    }else{
        values.addValues(points, needed, vals);
        points.addSortedInsexes(*needed.getVector());
        needed = MultiIndexSet();
        prepareSequence();
    }
    recomputeSurpluses();
}
void GridSequence::mergeRefinement(){
    if (needed.empty()) return; // nothing to do
    int num_all_points = getNumLoaded() + getNumNeeded();
    size_t num_vals = ((size_t) num_all_points) * ((size_t) num_outputs);
    std::vector<double> vals(num_vals, 0.0);
    values.setValues(vals);
    if (points.empty()){
        points = std::move(needed);
        needed = MultiIndexSet();
    }else{
        points.addMultiIndexSet(needed);
        needed = MultiIndexSet();
        prepareSequence();
    }
    surpluses.resize(num_vals);
    surpluses.shrink_to_fit();
    std::fill(surpluses.begin(), surpluses.end(), 0.0);
}

void GridSequence::beginConstruction(){
    dynamic_values = std::unique_ptr<SequenceConstructData>(new SequenceConstructData);
    if (points.empty()){
        dynamic_values->initial_points = std::move(needed);
        needed = MultiIndexSet();
    }
}
void GridSequence::writeConstructionDataBinary(std::ofstream &ofs) const{}
void GridSequence::writeConstructionData(std::ofstream &ofs) const{}
void GridSequence::readConstructionDataBinary(std::ifstream &ifs){}
void GridSequence::readConstructionData(std::ifstream &ifs){}
void GridSequence::getCandidateConstructionPoints(TypeDepth type, const std::vector<int> &weights, std::vector<double> &x, const std::vector<int> &level_limits){}
void GridSequence::getCandidateConstructionPoints(TypeDepth type, int output, std::vector<double> &x, const std::vector<int> &level_limits){}
void GridSequence::getCandidateConstructionPoints(std::function<double(const int *)> getTensorWeight, std::vector<double> &x, const std::vector<int> &level_limits){}
void GridSequence::loadConstructedPoint(const double x[], const std::vector<double> &y){}
void GridSequence::finishConstruction(){
    dynamic_values = std::unique_ptr<SequenceConstructData>();
}

void GridSequence::evaluate(const double x[], double y[]) const{
    std::vector<std::vector<double>> cache;
    cacheBasisValues<double>(x, cache);

    std::fill(y, y + num_outputs, 0.0);

    int num_points = points.getNumIndexes();

    Data2D<double> surps;
    surps.cload(num_outputs, num_points, surpluses.data());

    for(int i=0; i<num_points; i++){
        const int* p = points.getIndex(i);
        const double *s = surps.getCStrip(i);
        double basis_value = cache[0][p[0]];
        for(int j=1; j<num_dimensions; j++){
            basis_value *= cache[j][p[j]];
        }

        for(int k=0; k<num_outputs; k++){
            y[k] += basis_value * s[k];
        }
    }
}
void GridSequence::evaluateBatch(const double x[], int num_x, double y[]) const{
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    Data2D<double> yy; yy.load(num_outputs, num_x, y);
    #pragma omp parallel for
    for(int i=0; i<num_x; i++)
        evaluate(xx.getCStrip(i), yy.getStrip(i));
}

#ifdef Tasmanian_ENABLE_BLAS
void GridSequence::evaluateFastCPUblas(const double x[], double y[]) const{
    std::vector<double> fvalues(getNumPoints());
    evalHierarchicalFunctions(x, fvalues.data());
    TasBLAS::dgemv(num_outputs, points.getNumIndexes(), surpluses.data(), fvalues.data(), y);
}
void GridSequence::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    int num_points = points.getNumIndexes();
    Data2D<double> weights; weights.resize(num_points, num_x);
    evaluateHierarchicalFunctions(x, num_x, weights.getStrip(0));

    TasBLAS::dgemm(num_outputs, num_x, num_points, 1.0, surpluses.data(), weights.getStrip(0), 0.0, y);
}
#endif // Tasmanian_ENABLE_BLAS

#ifdef Tasmanian_ENABLE_CUDA
void GridSequence::evaluateFastGPUcublas(const double x[], double y[]) const{
    if (cuda_surpluses.size() == 0) cuda_surpluses.load(surpluses);

    std::vector<double> hweights(points.getNumIndexes());
    evalHierarchicalFunctions(x, hweights.data());

    cuda_engine.cublasDGEMM(num_outputs, 1, points.getNumIndexes(), 1.0, cuda_surpluses, hweights, 0.0, y);
}
void GridSequence::evaluateFastGPUcuda(const double x[], double y[]) const{ evaluateFastGPUcublas(x, y); }
void GridSequence::evaluateBatchGPUcublas(const double x[], int num_x, double y[]) const{
    if (cuda_surpluses.size() == 0) cuda_surpluses.load(surpluses);

    Data2D<double> hweights; hweights.resize(points.getNumIndexes(), num_x);
    evaluateHierarchicalFunctions(x, num_x, hweights.getStrip(0));

    cuda_engine.cublasDGEMM(num_outputs, num_x, points.getNumIndexes(), 1.0, cuda_surpluses, *(hweights.getVector()), 0.0, y);
}
void GridSequence::evaluateBatchGPUcuda(const double x[], int num_x, double y[]) const{
    if (cuda_surpluses.size() == 0) cuda_surpluses.load(surpluses);
    loadCudaNodes();

    int num_points = points.getNumIndexes();
    cudaDoubles gpu_x(num_dimensions, num_x, x);
    cudaDoubles gpu_basis(num_x, num_points);
    cudaDoubles gpu_result(num_x, num_outputs);

    evaluateHierarchicalFunctionsGPU(gpu_x.data(), num_x, gpu_basis.data());

    cuda_engine.cublasDGEMM(num_outputs, num_x, points.getNumIndexes(), 1.0, cuda_surpluses, gpu_basis, 0.0, gpu_result);
    gpu_result.unload(y);
}
#endif // Tasmanian_ENABLE_CUDA

#ifdef Tasmanian_ENABLE_MAGMA
void GridSequence::evaluateFastGPUmagma(int gpuID, const double x[], double y[]) const{
    if (cuda_surpluses.size() == 0) cuda_surpluses.load(surpluses);

    std::vector<double> hweights((size_t) points.getNumIndexes());
    evalHierarchicalFunctions(x, hweights.data());

    cuda_engine.magmaCudaDGEMM(gpuID, num_outputs, 1, points.getNumIndexes(), 1.0, cuda_surpluses, hweights, 0.0, y);
}
void GridSequence::evaluateBatchGPUmagma(int gpuID, const double x[], int num_x, double y[]) const{
    if (cuda_surpluses.size() == 0) cuda_surpluses.load(surpluses);
    loadCudaNodes();

    int num_points = points.getNumIndexes();
    cudaDoubles gpu_x(num_dimensions, num_x, x);
    cudaDoubles gpu_basis(num_x, num_points);
    cudaDoubles gpu_result(num_x, num_outputs);

    evaluateHierarchicalFunctionsGPU(gpu_x.data(), num_x, gpu_basis.data());

    cuda_engine.magmaCudaDGEMM(gpuID, num_outputs, num_x, points.getNumIndexes(), 1.0, cuda_surpluses, gpu_basis, 0.0, gpu_result);
    gpu_result.unload(y);
}
#endif

void GridSequence::integrate(double q[], double *conformal_correction) const{
    int num_points = points.getNumIndexes();
    std::fill(q, q + num_outputs, 0.0);

    Data2D<double> surp;
    surp.cload(num_outputs, num_points, surpluses.data());

    // for sequence grids, quadrature weights are expensive,
    // if using simple integration use the basis integral + surpluses, which is fast
    // if using conformal map, then we have to compute the expensive weights
    if (conformal_correction == 0){
        std::vector<double> integ;
        cacheBasisIntegrals(integ);
        for(int i=0; i<num_points; i++){
            const int* p = points.getIndex(i);
            double w = integ[p[0]];
            const double *s = surp.getCStrip(i);
            for(int j=1; j<num_dimensions; j++){
                w *= integ[p[j]];
            }
            for(int k=0; k<num_outputs; k++){
                q[k] += w * s[k];
            }
        }
    }else{
        std::vector<double> w(num_points);
        getQuadratureWeights(w.data());
        for(int i=0; i<num_points; i++){
            w[i] *= conformal_correction[i];
            const double *vals = values.getValues(i);
            for(int k=0; k<num_outputs; k++){
                q[k] += w[i] * vals[k];
            }
        }
    }
}

void GridSequence::evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const{
    int num_points = (points.empty()) ? needed.getNumIndexes() : points.getNumIndexes();
    Data2D<double> yy; yy.load(num_points, num_x, y);
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        evalHierarchicalFunctions(xx.getCStrip(i), yy.getStrip(i));
    }
}
void GridSequence::evalHierarchicalFunctions(const double x[], double fvalues[]) const{
    const MultiIndexSet& work = (points.empty()) ? needed : points;
    int num_points = work.getNumIndexes();

    std::vector<std::vector<double>> cache;
    cacheBasisValues<double>(x, cache);

    for(int i=0; i<num_points; i++){
        const int* p = work.getIndex(i);
        fvalues[i] = cache[0][p[0]];
        for(int j=1; j<num_dimensions; j++){
            fvalues[i] *= cache[j][p[j]];
        }
    }
}
#ifdef Tasmanian_ENABLE_CUDA
void GridSequence::evaluateHierarchicalFunctionsGPU(const double gpu_x[], int num_x, double gpu_y[]) const{
    loadCudaNodes();
    TasCUDA::devalseq(num_dimensions, num_x, max_levels, gpu_x, cuda_num_nodes, cuda_points, cuda_nodes, cuda_coeffs, gpu_y);
}
#endif
void GridSequence::setHierarchicalCoefficients(const double c[], TypeAcceleration acc){
    #ifdef Tasmanian_ENABLE_CUDA
    cuda_surpluses.clear();
    clearCudaNodes();
    #endif
    std::vector<double> *vals = 0;
    size_t num_ponits = (size_t) getNumPoints();
    size_t num_vals = num_ponits * ((size_t) num_outputs);
    if (!points.empty()){
        clearRefinement();
    }else{
        points = std::move(needed);
        needed = MultiIndexSet();
    }
    vals = values.aliasValues();
    vals->resize(num_vals);
    surpluses.resize(num_vals);
    std::copy_n(c, num_vals, surpluses.data());
    std::vector<double> x(((size_t) getNumPoints()) * ((size_t) num_dimensions));
    getPoints(x.data());
    switch(acc){
        #ifdef Tasmanian_ENABLE_BLAS
        case accel_cpu_blas: evaluateBatchCPUblas(x.data(), points.getNumIndexes(), vals->data()); break;
        #endif
        #ifdef Tasmanian_ENABLE_CUDA
        case accel_gpu_cublas: evaluateBatchGPUcublas(x.data(), points.getNumIndexes(), vals->data()); break;
        case accel_gpu_cuda:   evaluateBatchGPUcuda(x.data(), points.getNumIndexes(), vals->data()); break;
        #endif
        default:
            evaluateBatch(x.data(), points.getNumIndexes(), vals->data());
    }
}

void GridSequence::estimateAnisotropicCoefficients(TypeDepth type, int output, std::vector<int> &weights) const{
    double tol = 1000 * TSG_NUM_TOL;
    int num_points = points.getNumIndexes();
    std::vector<double> max_surp(num_points);

    if (output == -1){
        std::vector<double> nrm(num_outputs, 0.0);
        for(int i=0; i<num_points; i++){
            const double *val = values.getValues(i);
            int k=0;
            for(auto &n : nrm){
                double v = fabs(val[k++]);
                if (n < v) n = v;
            }
        }
        Data2D<double> surp;
        surp.cload(num_outputs, num_points, surpluses.data());
        #pragma omp parallel for
        for(int i=0; i<num_points; i++){
            const double *s = surp.getCStrip(i);
            double smax = 0.0;
            for(int k=0; k<num_outputs; k++){
                double v = fabs(s[k]) / nrm[k];
                if (smax < v) smax = v;
            }
            max_surp[i] = smax;
        }
    }else{
        Data2D<double> surp;
        surp.cload(num_outputs, num_points, surpluses.data());
        int i = 0;
        for(auto &m : max_surp) m = surp.getCStrip(i++)[output];
    }

    int n = 0, m;
    for(int i=0; i<num_points; i++){
        n += (max_surp[i] > tol) ? 1 : 0;
    }

    Data2D<double> A;
    std::vector<double> b(n);

    if ((type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved)){
        m = 2*num_dimensions + 1;
        A.resize(n, m);

        int count = 0;
        for(int c=0; c<num_points; c++){
            const int *indx = points.getIndex(c);
            if (max_surp[c] > tol){
                for(int j=0; j<num_dimensions; j++){
                    A.getStrip(j)[count] = ((double) indx[j]);
                }
                for(int j=0; j<num_dimensions; j++){
                    A.getStrip(j + num_dimensions)[count] = log((double) (indx[j] + 1));
                }
                A.getStrip(2*num_dimensions)[count] = 1.0;
                b[count++] = -log(max_surp[c]);
            }
        }
    }else{
        m = num_dimensions + 1;
        A.resize(n, m);

        int count = 0;
        for(int c=0; c<num_points; c++){
            const int *indx = points.getIndex(c);
            if (max_surp[c] > tol){
                for(int j=0; j<num_dimensions; j++){
                    A.getStrip(j)[count] = - ((double) indx[j]);
                }
                A.getStrip(num_dimensions)[count] = 1.0;
                b[count++] = log(max_surp[c]);
            }
        }
    }

    std::vector<double> x(m);
    TasmanianDenseSolver::solveLeastSquares(n, m, A.getStrip(0), b.data(), 1.E-5, x.data());

    weights.resize(--m);
    for(int j=0; j<m; j++){
        weights[j] = (int)(x[j] * 1000.0 + 0.5);
    }

    int min_weight = weights[0];  for(int j=1; j<num_dimensions; j++){  if (min_weight < weights[j]) min_weight = weights[j];  } // start by finding the largest weight
    if (min_weight < 0){ // all directions are diverging, default to isotropic total degree
        for(int j=0; j<num_dimensions; j++)  weights[j] = 1;
        if (m == 2*num_dimensions) for(int j=num_dimensions; j<2*num_dimensions; j++)  weights[j] = 0;
    }else{
        for(int j=0; j<num_dimensions; j++)  if ((weights[j]>0) && (weights[j]<min_weight)) min_weight = weights[j];  // find the smallest positive weight
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
}

void GridSequence::setAnisotropicRefinement(TypeDepth type, int min_growth, int output, const std::vector<int> &level_limits){
    clearRefinement();

    std::vector<int> weights;
    estimateAnisotropicCoefficients(type, output, weights);

    int level = 0; // find a better way to guess the offset, will tie to new refinement procedures

    do{
        MultiIndexSet total(num_dimensions);

        if ((type == type_qptotal) || (type == type_qptensor) || (type == type_qpcurved) || (type == type_qphyperbolic)){
            MultiIndexManipulations::selectTensors(++level, type, [&](int i) -> long long{ return OneDimensionalMeta::getQExact(i, rule); },
                                                   weights, total);
        }else{
            MultiIndexManipulations::selectTensors(++level, type, [&](int i) -> long long{ return i; }, weights, total);
        }

        total.diffSets(points, needed);

        if (!level_limits.empty()) MultiIndexManipulations::removeIndexesByLimit(level_limits, needed);

    }while(needed.empty() || (needed.getNumIndexes() < min_growth));

    prepareSequence();
}
void GridSequence::setSurplusRefinement(double tolerance, int output, const std::vector<int> &level_limits){
    clearRefinement();

    int num_points = points.getNumIndexes();
    std::vector<bool> flagged(num_points);

    std::vector<double> norm(num_outputs, 0.0);
    for(int i=0; i<num_points; i++){
        const double *val = values.getValues(i);
        for(int k=0; k<num_outputs; k++){
            double v = fabs(val[k]);
            if (norm[k] < v) norm[k] = v;
        }
    }

    Data2D<double> surp;
    surp.cload(num_outputs, num_points, surpluses.data());

    if (output == -1){
        for(int i=0; i<num_points; i++){
            const double *s = surp.getCStrip(i);
            double smax = fabs(s[0]) / norm[0];
            for(int k=1; k<num_outputs; k++){
                double v = fabs(s[k]) / norm[k];
                if (smax < v) smax = v;
            }
            flagged[i] = (smax > tolerance);
        }
    }else{
        for(int i=0; i<num_points; i++){
            flagged[i] = ((fabs(surp.getCStrip(i)[output]) / norm[output]) > tolerance);
        }
    }

    MultiIndexSet kids(num_dimensions);
    MultiIndexManipulations::selectFlaggedChildren(points, flagged, level_limits, kids);
    if (kids.getNumIndexes() > 0){
        kids.addMultiIndexSet(points);
        MultiIndexManipulations::completeSetToLower<int>(kids);

        kids.diffSets(points, needed);
        if (!needed.empty()) prepareSequence();
    }
}

void GridSequence::getPolynomialSpace(bool interpolation, int &n, int* &poly) const{
    MultiIndexSet space; // used only when interpolation is false
    const MultiIndexSet &work = (points.empty()) ? needed : points;
    if (!interpolation){ // when using interpolation, the polynomial space coincides with points/needed
        MultiIndexManipulations::createPolynomialSpace(work, [&](int l) -> int{ return OneDimensionalMeta::getQExact(l, rule); }, space);
    }
    const MultiIndexSet &result = (interpolation) ? work : space;

    n = result.getNumIndexes();
    poly = new int[result.getVector()->size()];
    std::copy(result.getVector()->begin(), result.getVector()->end(), poly);
}
const double* GridSequence::getSurpluses() const{
    return surpluses.data();
}
const int* GridSequence::getPointIndexes() const{
    return ((points.empty()) ? needed.getIndex(0) : points.getIndex(0));
}

void GridSequence::prepareSequence(){
    int mp = 0, mn = 0, max_level;
    if (needed.empty()){ // points must be non-empty
        MultiIndexManipulations::getMaxIndex(points, max_levels, mp);
    }else if (points.empty()){ // only needed, no points (right after creation)
        MultiIndexManipulations::getMaxIndex(needed, max_levels, mn);
    }else{ // both points and needed are set
        MultiIndexManipulations::getMaxIndex(points, max_levels, mp);
        mn = *std::max_element(needed.getVector()->begin(), needed.getVector()->end());
    }
    max_level = (mp > mn) ? mp : mn;
    max_level++;

    if ((size_t) max_level > nodes.size()){
        if (rule == rule_leja){
            Optimizer::getGreedyNodes<rule_leja>(max_level, nodes);
        }else if (rule == rule_maxlebesgue){
            Optimizer::getGreedyNodes<rule_maxlebesgue>(max_level, nodes);
        }else if (rule == rule_minlebesgue){
            Optimizer::getGreedyNodes<rule_minlebesgue>(max_level, nodes);
        }else if (rule == rule_mindelta){
            Optimizer::getGreedyNodes<rule_mindelta>(max_level, nodes);
        }else if (rule == rule_rleja){
            OneDimensionalNodes::getRLeja(max_level, nodes);
        }else if (rule == rule_rlejashifted){
            OneDimensionalNodes::getRLejaShifted(max_level, nodes);
        }
    }
    coeff.resize((size_t) max_level);
    coeff[0] = 1.0;
    for(int i=1; i<max_level; i++){
        coeff[i] = 1.0;
        for(int j=0; j<i; j++) coeff[i] *= (nodes[i] - nodes[j]);
    }
}

void GridSequence::cacheBasisIntegrals(std::vector<double> &integ) const{
    int max_level = max_levels[0];

    for(auto l: max_levels) if (max_level < l) max_level = l;

    integ.resize(++max_level, 0.0); // integrals of basis functions

    int n = 1 + max_level / 2; // number of Gauss-Legendre points needed to integrate the basis functions
    std::vector<double> lag_x, lag_w;
    OneDimensionalNodes::getGaussLegendre(n, lag_w, lag_x);

    for(int i=0; i<n; i++){
        double v = 1.0;
        for(int j=1; j<max_level; j++){
            v *= (lag_x[i] - nodes[j-1]);
            integ[j] += lag_w[i] * v / coeff[j];
        }
    }
    integ[0] = 2.0;
}

double GridSequence::evalBasis(const int f[], const int p[]) const{
    double v = 1.0;
    for(int j=0; j<num_dimensions; j++){
        double x = nodes[p[j]];
        double w = 1.0;
        for(int i=0; i<f[j]; i++){
            w *= (x - nodes[i]);
        }
        v *= w / coeff[f[j]];
    }
    return v;
}

void GridSequence::recomputeSurpluses(){
    int num_points = points.getNumIndexes();
    surpluses = *(values.aliasValues());

    Data2D<double> surp;
    surp.load(num_outputs, num_points, surpluses.data());

    std::vector<int> level;
    MultiIndexManipulations::computeLevels(points, level);
    int top_level = *std::max_element(level.begin(), level.end());

    Data2D<int> parents;
    MultiIndexManipulations::computeDAGup(points, parents);

    for(int l=1; l<=top_level; l++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<num_points; i++){
            if (level[i] == l){
                const int* p = points.getIndex(i);
                double *surpi = surp.getStrip(i);

                std::vector<int> monkey_count(top_level + 1);
                std::vector<int> monkey_tail(top_level + 1);
                std::vector<bool> used(num_points, false);

                int current = 0;

                monkey_count[0] = 0;
                monkey_tail[0] = i;

                while(monkey_count[0] < num_dimensions){
                    if (monkey_count[current] < num_dimensions){
                        int branch = parents.getStrip(monkey_tail[current])[monkey_count[current]];
                        if ((branch == -1) || (used[branch])){
                            monkey_count[current]++;
                        }else{
                            const double *branch_surp = surp.getCStrip(branch);;
                            double basis_value = evalBasis(points.getIndex(branch), p);
                            for(int k=0; k<num_outputs; k++){
                                 surpi[k] -= basis_value * branch_surp[k];
                            }
                            used[branch] = true;

                            monkey_count[++current] = 0;
                            monkey_tail[current] = branch;
                        }
                    }else{
                        monkey_count[--current]++;
                    }
                }
            }
        }
    }
}

void GridSequence::applyTransformationTransposed(double weights[]) const{
    const MultiIndexSet& work = (points.empty()) ? needed : points;
    int num_points = work.getNumIndexes();

    std::vector<int> level;
    MultiIndexManipulations::computeLevels(work, level);
    int top_level = *std::max_element(level.begin(), level.end());

    Data2D<int> parents;
    MultiIndexManipulations::computeDAGup(work, parents);

    std::vector<int> monkey_count(top_level + 1);
    std::vector<int> monkey_tail(top_level + 1);
    std::vector<bool> used(num_points);

    for(int l=top_level; l>0; l--){
        for(int i=0; i<num_points; i++){
            if (level[i] == l){
                const int* p = work.getIndex(i);
                int current = 0;

                monkey_count[0] = 0;
                monkey_tail[0] = i;
                std::fill(used.begin(), used.end(), false);

                while(monkey_count[0] < num_dimensions){
                    if (monkey_count[current] < num_dimensions){
                        int branch = parents.getStrip(monkey_tail[current])[monkey_count[current]];
                        if ((branch == -1) || used[branch]){
                            monkey_count[current]++;
                        }else{
                            weights[branch] -= weights[i] * evalBasis(work.getIndex(branch), p);
                            used[branch] = true;

                            monkey_count[++current] = 0;
                            monkey_tail[current] = branch;
                        }
                    }else{
                        monkey_count[--current]++;
                    }
                }
            }
        }
    }
}

void GridSequence::clearAccelerationData(){
    #ifdef Tasmanian_ENABLE_CUDA
    cuda_engine.reset();
    cuda_surpluses.clear();
    clearCudaNodes();
    #endif
}


}

#endif
