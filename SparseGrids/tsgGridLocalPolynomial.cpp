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

#ifndef __TASMANIAN_SPARSE_GRID_LPOLY_CPP
#define __TASMANIAN_SPARSE_GRID_LPOLY_CPP

#include "tsgGridLocalPolynomial.hpp"

#include "tsgHiddenExternals.hpp"
#include "tsgCudaMacros.hpp" // remove later

namespace TasGrid{

GridLocalPolynomial::GridLocalPolynomial() : num_dimensions(0), num_outputs(0), order(1), top_level(0),
                         points(0), needed(0), values(0), rule(0),
                         accel(0), sparse_affinity(0)  {}
GridLocalPolynomial::GridLocalPolynomial(const GridLocalPolynomial &pwpoly) : num_dimensions(0), num_outputs(0), order(1), top_level(0),
                         points(0), needed(0), values(0), rule(0),
                         accel(0), sparse_affinity(0)
{
    copyGrid(&pwpoly);
}

GridLocalPolynomial::~GridLocalPolynomial(){
    reset();
}
void GridLocalPolynomial::reset(bool clear_rule){
    clearAccelerationData();
    num_dimensions = num_outputs = top_level = 0;
    if (points != 0){ delete points; points = 0; }
    if (needed != 0){ delete needed; needed = 0; }
    if (values != 0){ delete values; values = 0; }
    if (clear_rule){ rule = 0; order = 1; }
    parents.load(0, 0, 0);
    sparse_affinity = 0;
    surpluses.clear();
}

void GridLocalPolynomial::write(std::ofstream &ofs) const{
    using std::endl;

    ofs << std::scientific; ofs.precision(17);
    ofs << num_dimensions << " " << num_outputs << " " << order << " " << top_level << endl;
    if (num_dimensions > 0){
        ofs << OneDimensionalMeta::getIORuleString(rule->getType()) << endl;
        if (points == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            points->write(ofs);
        }
        if (surpluses.getNumStrips() == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            for(auto s: *surpluses.getVector()) ofs << " " << s;
            ofs << endl;
        }
        if (needed == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            needed->write(ofs);
        }
        if (parents.getNumStrips() == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1";
            size_t n = parents.getTotalEntries();
            const int *p = parents.getCStrip(0);
            for(size_t i=0; i<n; i++){ ofs << " " << p[i]; } ofs << endl;
        }
        int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();
        ofs << roots.size(); for(auto const &r : roots){ ofs << " " << r; } ofs << endl;
        ofs << pntr[0]; for(int i=1; i<=num_points; i++){  ofs << " " << pntr[i];  }  ofs << endl;
        ofs << indx[0]; for(int i=1; i<pntr[num_points]; i++){  ofs << " " << indx[i];  }  ofs << endl;
        if (num_outputs > 0) values->write(ofs);
    }
}
void GridLocalPolynomial::writeBinary(std::ofstream &ofs) const{
    int dims[4];
    dims[0] = num_dimensions;
    dims[1] = num_outputs;
    dims[2] = order;
    dims[3] = top_level;
    ofs.write((char*) dims, 4 * sizeof(int));
    if (num_dimensions > 0){
        dims[0] = OneDimensionalMeta::getIORuleInt(rule->getType());
        ofs.write((char*) dims, sizeof(int));

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
        if (surpluses.getNumStrips() == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            ofs.write((char*) surpluses.getCStrip(0), surpluses.getTotalEntries() * sizeof(double));
        }
        if (parents.getTotalEntries() == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            ofs.write((char*) parents.getCStrip(0), parents.getTotalEntries() * sizeof(int));
        }

        int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();
        int num_roots = (int) roots.size();
        ofs.write((char*) &num_roots, sizeof(int));
        ofs.write((char*) roots.data(), num_roots * sizeof(int));
        ofs.write((char*) pntr.data(), (num_points+1) * sizeof(int));
        ofs.write((char*) indx.data(), pntr[num_points] * sizeof(int));

        if (num_outputs > 0) values->writeBinary(ofs);
    }
}
void GridLocalPolynomial::read(std::ifstream &ifs){
    reset();
    ifs >> num_dimensions >> num_outputs >> order >> top_level;
    if (num_dimensions > 0){
        int flag;
        std::string T;
        ifs >> T;
        TypeOneDRule crule = OneDimensionalMeta::getIORuleString(T.c_str());
        if (crule == rule_localp){
            rule = &rpoly;
        }else if (crule == rule_semilocalp){
            rule = &rsemipoly;
        }else if (crule == rule_localp0){
            rule = &rpoly0;
        }else if (crule == rule_localpb){
            rule = &rpolyb;
        }
        if (order == 0){
            rule = &rpolyc;
        }
        rule->setMaxOrder(order);

        ifs >> flag; if (flag == 1){ points = new IndexSet(num_dimensions); points->read(ifs); }
        ifs >> flag;
        if (flag == 1){
            surpluses.resize(num_outputs, points->getNumIndexes());
            for(auto &s: *surpluses.getVector()) ifs >> s;
        }
        ifs >> flag;
        if (flag == 1){
            needed = new IndexSet(num_dimensions);
            needed->read(ifs);
        }
        ifs >> flag;
        if (flag == 1){
            parents.resize(rule->getMaxNumParents() * num_dimensions, points->getNumIndexes());
            size_t n = parents.getTotalEntries();
            int *p = parents.getStrip(0);
            for(size_t i=0; i<n; i++){ ifs >> p[i]; }
        }

        int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();

        int num_roots;
        ifs >> num_roots;
        roots.resize(num_roots);
        for(int i=0; i<num_roots; i++) ifs >> roots[i];
        pntr.resize(num_points + 1);
        for(int i=0; i<=num_points; i++) ifs >> pntr[i];
        if (pntr[num_points] > 0){
            indx.resize(pntr[num_points]);
            for(int i=0; i<pntr[num_points]; i++) ifs >> indx[i];
        }else{
            indx.resize(1);  ifs >> indx[0]; // there is a special case when the grid has only one point without any children
        }

        if (num_outputs > 0){ values = new StorageSet(0, 0); values->read(ifs); }
    }
}
void GridLocalPolynomial::readBinary(std::ifstream &ifs){
    reset();
    int dims[4];
    ifs.read((char*) dims, 4 * sizeof(int));
    num_dimensions = dims[0];
    num_outputs = dims[1];
    order = dims[2];
    top_level = dims[3];
    if (num_dimensions > 0){
        ifs.read((char*) dims, sizeof(int));
        TypeOneDRule crule = OneDimensionalMeta::getIORuleInt(dims[0]);
        if (crule == rule_localp){
            rule = &rpoly;
        }else if (crule == rule_semilocalp){
            rule = &rsemipoly;
        }else if (crule == rule_localp0){
            rule = &rpoly0;
        }else if (crule == rule_localpb){
            rule = &rpolyb;
        }
        if (order == 0){
            rule = &rpolyc;
        }
        rule->setMaxOrder(order);

        char flag;
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ points = new IndexSet(num_dimensions); points->readBinary(ifs); }
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ needed = new IndexSet(num_dimensions); needed->readBinary(ifs); }

        ifs.read((char*) &flag, sizeof(char));
        if (flag == 'y'){
            surpluses.resize(num_outputs, points->getNumIndexes());
            ifs.read((char*) surpluses.getStrip(0), surpluses.getTotalEntries() * sizeof(double));
        }

        ifs.read((char*) &flag, sizeof(char));
        if (flag == 'y'){
            parents.resize(rule->getMaxNumParents() * num_dimensions, points->getNumIndexes());
            ifs.read((char*) parents.getStrip(0), parents.getTotalEntries() * sizeof(int));
        }

        int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();
        int num_roots;
        ifs.read((char*) &num_roots, sizeof(int));
        roots.resize(num_roots);
        ifs.read((char*) roots.data(), num_roots * sizeof(int));

        pntr.resize(num_points + 1);
        ifs.read((char*) pntr.data(), (num_points+1) * sizeof(int));

        if (pntr[num_points] > 0){
            indx.resize(pntr[num_points]);
            ifs.read((char*) indx.data(), pntr[num_points] * sizeof(int));
        }else{
            indx.resize(1);
            ifs.read((char*) indx.data(), sizeof(int));
        }

        if (num_outputs > 0){ values = new StorageSet(0, 0); values->readBinary(ifs); }
    }
}

void GridLocalPolynomial::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, int corder, TypeOneDRule crule, const std::vector<int> &level_limits){
    reset();
    num_dimensions = cnum_dimensions;
    num_outputs = cnum_outputs;
    order = corder;

    TypeOneDRule effective_rule = ((crule == rule_semilocalp) && ((order == -1) || (order > 1))) ? rule_semilocalp : rule_localp; // semi-localp of order 1 is identical to localp
    if (crule == rule_localp0) effective_rule = rule_localp0;
    if (crule == rule_localpb) effective_rule = rule_localpb;

    if (effective_rule == rule_localp){
        rule = &rpoly;
    }else if (effective_rule == rule_semilocalp){
        rule = &rsemipoly;
    }else if (effective_rule == rule_localp0){
        rule = &rpoly0;
    }else if (effective_rule == rule_localpb){
        rule = &rpolyb;
    }
    if (order == 0){ // if the rule is zero-order
        rule = &rpolyc;
    }
    rule->setMaxOrder(order);

    IndexManipulator IM(num_dimensions);
    IndexSet* deltas = IM.selectTensors(depth, type_level, std::vector<int>(), rule_leja);

    // Limits come here
    if (!level_limits.empty()){
        IndexSet *limited = IM.removeIndexesByLimit(deltas, level_limits);
        if (limited != 0){
            delete deltas;
            deltas = limited;
        }
    }

    needed = IM.generatePointsFromDeltas(deltas, [&](int level) -> int { return rule->getNumPoints(level); });
    delete deltas;

    buildTree();

    if (num_outputs == 0){
        points = needed;
        needed = 0;
        IM.computeDAGupLocal(points, rule, parents);
    }else{
        values = new StorageSet(num_outputs, needed->getNumIndexes());
    }
}

void GridLocalPolynomial::copyGrid(const GridLocalPolynomial *pwpoly){
    reset();
    num_dimensions = pwpoly->num_dimensions;
    num_outputs = pwpoly->num_outputs;
    order = pwpoly->order;

    if (pwpoly->rule->getType() == rule_localp){
        rule = &rpoly;
    }else if (pwpoly->rule->getType() == rule_semilocalp){
        rule = &rsemipoly;
    }else if (pwpoly->rule->getType() == rule_localp0){
        rule = &rpoly0;
    }else if (pwpoly->rule->getType() == rule_localpb){
        rule = &rpolyb;
    }
    if (pwpoly->rule->getMaxOrder() == 0){
        rule = &rpolyc;
    }
    rule->setMaxOrder(pwpoly->rule->getMaxOrder());

    if (pwpoly->points != 0) points = new IndexSet(pwpoly->points);
    if (pwpoly->needed != 0) needed = new IndexSet(pwpoly->needed);
    parents = pwpoly->parents;

    buildTree();

    if (pwpoly->values != 0) values = new StorageSet(pwpoly->values);

    if ((points != 0) && (num_outputs > 0)){ // points are loaded
        surpluses.resize(num_outputs, points->getNumIndexes());
        *surpluses.getVector() = *pwpoly->surpluses.getVector(); // copy assignment
    }
}

int GridLocalPolynomial::getNumDimensions() const{ return num_dimensions; }
int GridLocalPolynomial::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridLocalPolynomial::getRule() const{ return rule->getType(); }
int GridLocalPolynomial::getOrder() const{ return order; }

int GridLocalPolynomial::getNumLoaded() const{ return (((points == 0) || (num_outputs == 0)) ? 0 : points->getNumIndexes()); }
int GridLocalPolynomial::getNumNeeded() const{ return ((needed == 0) ? 0 : needed->getNumIndexes()); }
int GridLocalPolynomial::getNumPoints() const{ return ((points == 0) ? getNumNeeded() : points->getNumIndexes()); }

void GridLocalPolynomial::getLoadedPoints(double *x) const{
    int num_points = points->getNumIndexes();
    Data2D<double> split;
    split.load(num_dimensions, num_points, x);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        double *xx = split.getStrip(i);
        for(int j=0; j<num_dimensions; j++){
            xx[j] = rule->getNode(p[j]);
        }
    }
}
void GridLocalPolynomial::getNeededPoints(double *x) const{
    int num_points = needed->getNumIndexes();
    Data2D<double> split;
    split.load(num_dimensions, num_points, x);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed->getIndex(i);
        double *xx = split.getStrip(i);
        for(int j=0; j<num_dimensions; j++){
            xx[j] = rule->getNode(p[j]);
        }
    }
}
void GridLocalPolynomial::getPoints(double *x) const{
    if (points == 0){ getNeededPoints(x); }else{ getLoadedPoints(x); }
}

void GridLocalPolynomial::evaluate(const double x[], double y[]) const{
    std::vector<int> monkey_count(top_level+1);
    std::vector<int> monkey_tail(top_level+1);

    bool isSupported;
    int offset;

    std::fill(y, y + num_outputs, 0.0);

    for(auto const &r : roots){
        double basis_value = evalBasisSupported(points->getIndex(r), x, isSupported);

        if (isSupported){
            const double *s = surpluses.getCStrip(r);
            for(int k=0; k<num_outputs; k++) y[k] += basis_value * s[k];

            int current = 0;
            monkey_tail[0] = r;
            monkey_count[0] = pntr[r];

            while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                    offset = indx[monkey_count[current]];
                    basis_value = evalBasisSupported(points->getIndex(offset), x, isSupported);
                    if (isSupported){
                        s = surpluses.getCStrip(offset);
                        for(int k=0; k<num_outputs; k++) y[k] += basis_value * s[k];

                        monkey_tail[++current] = offset;
                        monkey_count[current] = pntr[offset];
                    }else{
                        monkey_count[current]++;
                    }
                }else{
                    monkey_count[--current]++;
                }
            }
        }
    }
}
void GridLocalPolynomial::evaluateFastCPUblas(const double x[], double y[]) const{ evaluate(x, y); }
// standard BLAS cannot accelerate dense matrix times a sparse vector, fallback to regular evaluate()

#ifdef Tasmanian_ENABLE_CUDA
void GridLocalPolynomial::evaluateFastGPUcublas(const double x[], double y[]) const{
    int num_points = points->getNumIndexes();
    if (cuda_surpluses.size() == 0) loadCudaData();

    std::vector<int> sindx;
    std::vector<double> svals;
    int num_nz;
    buildSparseVector<0>(points, x, num_nz, sindx, svals);
    buildSparseVector<1>(points, x, num_nz, sindx, svals);
    cuda_engine.cusparseMatveci(num_outputs, num_points, 1.0, cuda_surpluses, sindx, svals, 0.0, y);
}
#else
void GridLocalPolynomial::evaluateFastGPUcublas(const double[], double[]) const{}
#endif
// evaluation of a single x cannot be accelerated with a gpu (not parallelizable), do that on the CPU and use the GPU only for the case of many outputs
void GridLocalPolynomial::evaluateFastGPUcuda(const double x[], double y[]) const{ evaluateFastGPUcublas(x, y); }
void GridLocalPolynomial::evaluateFastGPUmagma(int, const double x[], double y[]) const{ evaluate(x, y); }

void GridLocalPolynomial::evaluateBatch(const double x[], int num_x, double y[]) const{
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    Data2D<double> yy; yy.load(num_outputs, num_x, y);
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        evaluate(xx.getCStrip(i), yy.getStrip(i));
    }
}
void GridLocalPolynomial::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    if ((sparse_affinity == 1) || ((sparse_affinity == 0) && (num_outputs <= TSG_LOCALP_BLAS_NUM_OUTPUTS))){
        evaluateBatch(x, num_x, y);
        return;
    }

    int *sindx, *spntr;
    double *svals;
    buildSpareBasisMatrix(x, num_x, 32, spntr, sindx, svals); // build sparse matrix corresponding to x

    int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();
    double nnz = (double) spntr[num_x];
    double total_size = ((double) num_x) * ((double) num_points);

    if ((sparse_affinity == -1) || ((sparse_affinity == 0) && (nnz / total_size > 0.1))){
        // potentially wastes a lot of memory
        Data2D<double> A;
        A.resize(num_points, num_x, 0.0);
        for(int i=0; i<num_x; i++){
            double *row = A.getStrip(i);
            for(int j=spntr[i]; j<spntr[i+1]; j++) row[sindx[j]] = svals[j];
        }
        TasBLAS::dgemm(num_outputs, num_x, num_points, 1.0, surpluses.getCStrip(0), A.getCStrip(0), 0.0, y);
    }else{
        Data2D<double> yy; yy.load(num_outputs, num_x, y);
        #pragma omp parallel for
        for(int i=0; i<num_x; i++){
            double *this_y = yy.getStrip(i);
            std::fill(this_y, this_y + num_outputs, 0.0);
            for(int j=spntr[i]; j<spntr[i+1]; j++){
                double v = svals[j];
                const double *s = surpluses.getCStrip(sindx[j]);
                for(int k=0; k<num_outputs; k++) this_y[k] += v * s[k];
            }
        }
    }

    delete[] sindx;
    delete[] spntr;
    delete[] svals;
}

#ifdef Tasmanian_ENABLE_CUDA
void GridLocalPolynomial::evaluateBatchGPUcublas(const double x[], int num_x, double y[]) const{
    int num_points = points->getNumIndexes();
    makeCheckAccelerationData(accel_gpu_cublas);
    checkAccelerationGPUValues();
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;

    int *sindx, *spntr;
    double *svals;
    buildSpareBasisMatrix(x, num_x, 32, spntr, sindx, svals); // build sparse matrix corresponding to x

    gpu_acc->cusparseMatmul(true, num_points, num_outputs, num_x, spntr, sindx, svals, 0, y); // true for pointers residing on the cpu

    delete[] svals;
    delete[] sindx;
    delete[] spntr;
}
void GridLocalPolynomial::evaluateBatchGPUcuda(const double x[], int num_x, double y[]) const{
    if ((order == -1) || (order > 2)){ // GPU evaluations are availabe only for order 0, 1, and 2. Cubic will come later, but higher order will not be supported
        evaluateBatchGPUcublas(x, num_x, y);
        return;
    }
    int num_points = points->getNumIndexes();
    makeCheckAccelerationData(accel_gpu_cuda);
    checkAccelerationGPUValues();
    checkAccelerationGPUNodes();
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;

    bool useDense = (sparse_affinity == -1) || ((sparse_affinity == 0) && (num_dimensions > 6)); // dimension is the real criteria here

    if (useDense){
        double *gpu_x = TasCUDA::cudaSend<double>(num_x * num_dimensions, x);
        double *gpu_weights = TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) points->getNumIndexes()));
        double *gpu_result = TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) values->getNumOutputs()));

        buildDenseBasisMatrixGPU(gpu_x, num_x, gpu_weights);
        gpu_acc->cublasDGEMM(false, values->getNumOutputs(), num_x, num_points, gpu_weights, gpu_result);
        //TasCUDA::cudaDgemm(values->getNumOutputs(), num_x, num_points, gpu_acc->getGPUValues(), gpu_weights, gpu_result);

        TasCUDA::cudaRecv<double>(num_x * values->getNumOutputs(), gpu_result, y);

        TasCUDA::cudaDel<double>(gpu_result);
        TasCUDA::cudaDel<double>(gpu_weights);
        TasCUDA::cudaDel<double>(gpu_x);
    }else{
        double *gpu_x = TasCUDA::cudaSend<double>(num_x * num_dimensions, x);
        double *gpu_y = TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) num_outputs));

        int *gpu_spntr, *gpu_sindx, num_nz = 0;
        double *gpu_svals;
        buildSparseBasisMatrixGPU(gpu_x, num_x, gpu_spntr, gpu_sindx, gpu_svals, num_nz);

        if (num_outputs == 1){
            gpu_acc->cusparseMatvec(num_points, num_x, gpu_spntr, gpu_sindx, gpu_svals, num_nz, gpu_y);
        }else{
            gpu_acc->cusparseMatmul(false, num_points, num_outputs, num_x, gpu_spntr, gpu_sindx, gpu_svals, num_nz, gpu_y); // false for pointers already living on the gpu
        }
        //TasCUDA::cudaSparseMatmul(num_x, num_outputs, num_nz, gpu_spntr, gpu_sindx, gpu_svals, gpu_acc->getGPUValues(), gpu_y);

        TasCUDA::cudaRecv<double>(((size_t) num_x) * ((size_t) num_outputs), gpu_y, y);

        TasCUDA::cudaDel<int>(gpu_spntr);
        TasCUDA::cudaDel<int>(gpu_sindx);
        TasCUDA::cudaDel<double>(gpu_svals);
        TasCUDA::cudaDel<double>(gpu_y);
        TasCUDA::cudaDel<double>(gpu_x);
    }
}
#else
void GridLocalPolynomial::evaluateBatchGPUcublas(const double[], int, double[]) const{}
void GridLocalPolynomial::evaluateBatchGPUcuda(const double[], int, double[]) const{}
#endif // Tasmanian_ENABLE_CUDA

void GridLocalPolynomial::evaluateBatchGPUmagma(int, const double x[], int num_x, double y[]) const{
    evaluateBatchGPUcublas(x, num_x, y);
}

void GridLocalPolynomial::loadNeededPoints(const double *vals, TypeAcceleration){
    #ifdef Tasmanian_ENABLE_CUDA
    clearCudaLoadedData();
    #endif
    if (points == 0){
        values->setValues(vals);
        points = needed;
        needed = 0;
    }else if (needed == 0){
        values->setValues(vals);
    }else{
        values->addValues(points, needed, vals);
        points->addIndexSet(needed);
        delete needed; needed = 0;
        buildTree();
    }
    if (accel != 0) accel->resetGPULoadedData();
    recomputeSurpluses();
}
void GridLocalPolynomial::mergeRefinement(){
    if (needed == 0) return; // nothing to do
    #ifdef Tasmanian_ENABLE_CUDA
    clearCudaLoadedData();
    #endif
    int num_all_points = getNumLoaded() + getNumNeeded();
    size_t num_vals = ((size_t) num_all_points) * ((size_t) num_outputs);
    std::vector<double> vals(num_vals, 0.0);
    values->setValues(vals);
    if (points == 0){
        points = needed;
        needed = 0;
    }else{
        points->addIndexSet(needed);
        delete needed; needed = 0;
        buildTree();
    }
    surpluses.resize(num_outputs, num_all_points);
    std::fill(surpluses.getVector()->begin(), surpluses.getVector()->end(), 0.0);
}

void GridLocalPolynomial::getInterpolationWeights(const double x[], double *weights) const{
    IndexSet *work = (points == 0) ? needed : points;

    std::vector<int> active_points(0);
    std::fill(weights, weights + work->getNumIndexes(), 0.0);

    std::vector<int> monkey_count(top_level+1);
    std::vector<int> monkey_tail(top_level+1);

    double basis_value;
    bool isSupported;
    int offset;

    for(unsigned r=0; r<roots.size(); r++){

        basis_value = evalBasisSupported(work->getIndex(roots[r]), x, isSupported);

        if (isSupported){
            active_points.push_back(roots[r]);
            weights[roots[r]] = basis_value;

            int current = 0;
            monkey_tail[0] = roots[r];
            monkey_count[0] = pntr[roots[r]];

            while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                    offset = indx[monkey_count[current]];

                    basis_value = evalBasisSupported(work->getIndex(offset), x, isSupported);

                    if (isSupported){
                        active_points.push_back(offset);
                        weights[offset] = basis_value;

                        monkey_tail[++current] = offset;
                        monkey_count[current] = pntr[offset];
                    }else{
                        monkey_count[current]++;
                    }
                }else{
                    monkey_count[--current]++;
                }
            }
        }
    }

    // apply the transpose of the surplus transformation
    const Data2D<int> *dagUp;
    Data2D<int> lparents;
    if (parents.getNumStrips() == work->getNumIndexes()){
        dagUp = &parents;
    }else{
        IndexManipulator IM(num_dimensions);
        IM.computeDAGupLocal(work, rule, lparents);
        dagUp = &lparents;
    }

    std::vector<int> level(active_points.size());
    int active_top_level = 0;
    for(size_t i=0; i<active_points.size(); i++){
        const int *p = work->getIndex(active_points[i]);
        int current_level = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            current_level += rule->getLevel(p[j]);
        }
        if (active_top_level < current_level) active_top_level = current_level;
        level[i] = current_level;
    }

    std::vector<bool> used(work->getNumIndexes());
    std::vector<double> node(num_dimensions);
    int max_parents = rule->getMaxNumParents() * num_dimensions;

    for(int l=active_top_level; l>0; l--){
        for(size_t i=0; i<active_points.size(); i++){
            if (level[i] == l){
                const int* p = work->getIndex(active_points[i]);
                for(int j=0; j<num_dimensions; j++) node[j] = rule->getNode(p[j]);

                std::fill(used.begin(), used.end(), false);

                monkey_count[0] = 0;
                monkey_tail[0] = active_points[i];
                int current = 0;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp->getCStrip(monkey_tail[current])[monkey_count[current]];
                        if ((branch == -1) || used[branch]){
                            monkey_count[current]++;
                        }else{
                            const int *func = work->getIndex(branch);
                            basis_value = rule->evalRaw(func[0], node[0]);
                            for(int j=1; j<num_dimensions; j++) basis_value *= rule->evalRaw(func[j], node[j]);
                            weights[branch] -= weights[active_points[i]] * basis_value;
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

void GridLocalPolynomial::evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
    Data2D<double> yy; yy.load(num_points, num_x, y);
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        const double *this_x = xx.getCStrip(i);
        double *this_y = yy.getStrip(i);
        bool dummy;
        for(int j=0; j<num_points; j++){
            this_y[j] = evalBasisSupported(work->getIndex(j), this_x, dummy);
        }
    }
}

void GridLocalPolynomial::recomputeSurpluses(){
    int num_points = points->getNumIndexes();

    surpluses.resize(num_outputs, num_points);
    *surpluses.getVector() = *values->aliasValues(); // copy assignment

    Data2D<int> dagUp;
    IndexManipulator IM(num_dimensions);
    IM.computeDAGupLocal(points, rule, dagUp);

    int max_parents = rule->getMaxNumParents() * num_dimensions;

    std::vector<int> level(num_points);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        int current_level = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            current_level += rule->getLevel(p[j]);
        }
        level[i] = current_level;
    }

    for(int l=1; l<=top_level; l++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<num_points; i++){
            if (level[i] == l){
                const int* p = points->getIndex(i);
                std::vector<double> x(num_dimensions);
                double *surpi = surpluses.getStrip(i);
                for(int j=0; j<num_dimensions; j++) x[j] = rule->getNode(p[j]);

                std::vector<int> monkey_count(top_level + 1);
                std::vector<int> monkey_tail(top_level + 1);
                std::vector<bool> used(num_points);
                std::fill(used.begin(), used.end(), false);

                int current = 0;

                monkey_count[0] = 0;
                monkey_tail[0] = i;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp.getStrip(monkey_tail[current])[monkey_count[current]];
                        if ((branch == -1) || (used[branch])){
                            monkey_count[current]++;
                        }else{
                            const double *branch_surp = surpluses.getCStrip(branch);
                            double basis_value = evalBasisRaw(points->getIndex(branch), x.data());
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

double GridLocalPolynomial::evalBasisRaw(const int point[], const double x[]) const{
    double f = rule->evalRaw(point[0], x[0]);
    for(int j=1; j<num_dimensions; j++) f *= rule->evalRaw(point[j], x[j]);
    return f;
}
double GridLocalPolynomial::evalBasisSupported(const int point[], const double x[], bool &isSupported) const{
    double f = rule->evalSupport(point[0], x[0], isSupported);
    if (!isSupported) return 0.0;
    for(int j=1; j<num_dimensions; j++){
        f *= rule->evalSupport(point[j], x[j], isSupported);
        if (!isSupported) return 0.0;
    }
    return f;
}

void GridLocalPolynomial::buildSpareBasisMatrix(const double x[], int num_x, int num_chunk, int* &spntr, int* &sindx, double* &svals) const{
    std::vector<std::vector<int>> tindx;
    std::vector<std::vector<double>> tvals;
    std::vector<int> numnz;
    buildSparseMatrixBlockForm(x, num_x, num_chunk, numnz, tindx, tvals);

    spntr = new int[num_x + 1];

    int nz = 0;
    for(int i=0; i<num_x; i++){
        spntr[i] = nz;
        nz += numnz[i];
    }
    spntr[num_x] = nz;

    sindx = new int[nz];
    svals = new double[nz];

    int c = 0;
    for(auto &idx : tindx) for(auto i: idx) sindx[c++] = i;
    c = 0;
    for(auto &vls : tvals) for(auto v: vls) svals[c++] = v;
}
void GridLocalPolynomial::buildSpareBasisMatrix(const double x[], int num_x, int num_chunk, std::vector<int> &spntr, std::vector<int> &sindx, std::vector<double> &svals) const{
    std::vector<std::vector<int>> tindx;
    std::vector<std::vector<double>> tvals;
    std::vector<int> numnz;
    buildSparseMatrixBlockForm(x, num_x, num_chunk, numnz, tindx, tvals);

    spntr.resize(num_x + 1);

    int nz = 0;
    auto inz = numnz.begin();
    for(auto &s: spntr){
        s = nz;
        nz += *inz++;
    }

    sindx.resize(nz);
    svals.resize(nz);

    auto ii = sindx.begin();
    for(auto &idx : tindx) for(auto i: idx) *ii++ = i;
    auto iv = svals.begin();
    for(auto &vls : tvals) for(auto v: vls) *iv++ = v;
}
void GridLocalPolynomial::buildSpareBasisMatrixStatic(const double x[], int num_x, int num_chunk, int *spntr, int *sindx, double *svals) const{
    std::vector<std::vector<int>> tindx;
    std::vector<std::vector<double>> tvals;
    std::vector<int> numnz;
    buildSparseMatrixBlockForm(x, num_x, num_chunk, numnz, tindx, tvals);

    int nz = 0;
    for(int i=0; i<num_x; i++){
        spntr[i] = nz;
        nz += numnz[i];
    }
    spntr[num_x] = nz;
    int c = 0;
    for(auto &idx : tindx) for(auto i: idx) sindx[c++] = i;
    c = 0;
    for(auto &vls : tvals) for(auto v: vls) svals[c++] = v;
}
int GridLocalPolynomial::getSpareBasisMatrixNZ(const double x[], int num_x) const{
    IndexSet *work = (points == 0) ? needed : points;

    std::vector<int> sindx; // dummy vectors, never referenced
    std::vector<double> svals;

    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    std::vector<int> num_nz(num_x);

    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        buildSparseVector<0>(work, xx.getCStrip(i), num_nz[i], sindx, svals);
    }

    int total_nz = 0;
    for(auto n: num_nz) total_nz += n;

    return total_nz;
}

void GridLocalPolynomial::buildSparseMatrixBlockForm(const double x[], int num_x, int num_chunk, std::vector<int> &numnz, std::vector<std::vector<int>> &tindx, std::vector<std::vector<double>> &tvals) const{
    // numnz will be resized to (num_x + 1) with the last entry set to a dummy zero
    numnz.resize(num_x + 1);
    int num_blocks = num_x / num_chunk + ((num_x % num_chunk != 0) ? 1 : 0);

    tindx.resize(num_blocks);
    tvals.resize(num_blocks);

    const IndexSet *work = (points != 0) ? points : needed;
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);

    #pragma omp parallel for
    for(int b=0; b<num_blocks; b++){
        tindx[b].resize(0);
        tvals[b].resize(0);
        int chunk_size = (b < num_blocks - 1) ? num_chunk : (num_x - (num_blocks - 1) * num_chunk);
        for(int i = b * num_chunk; i < b * num_chunk + chunk_size; i++){
            buildSparseVector<2>(work, xx.getCStrip(i), numnz[i], tindx[b], tvals[b]);
        }
    }
    numnz[num_x] = 0;
}

#ifdef Tasmanian_ENABLE_CUDA
void GridLocalPolynomial::buildDenseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, double gpu_y[]) const{
    int num_points = getNumPoints();
    makeCheckAccelerationData(accel_gpu_cuda);
    checkAccelerationGPUNodes();
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;
    TasCUDA::devalpwpoly(order, rule->getType(), num_dimensions, cpu_num_x, num_points, gpu_x, gpu_acc->getGPUNodes(), gpu_acc->getGPUSupport(), gpu_y);
}
void GridLocalPolynomial::buildSparseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nz) const{
    int num_points = getNumPoints();
    makeCheckAccelerationData(accel_gpu_cuda);
    checkAccelerationGPUNodes();
    checkAccelerationGPUHierarchy();
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;
    TasCUDA::devalpwpoly_sparse(order, rule->getType(), num_dimensions, cpu_num_x, num_points, gpu_x, gpu_acc->getGPUNodes(), gpu_acc->getGPUSupport(),
                                gpu_acc->getGPUpntr(), gpu_acc->getGPUindx(), (int) roots.size(), gpu_acc->getGPUroots(),
                                gpu_spntr, gpu_sindx, gpu_svals, num_nz);
}

#else
void GridLocalPolynomial::buildDenseBasisMatrixGPU(const double*, int, double*) const{}
void GridLocalPolynomial::buildSparseBasisMatrixGPU(const double*, int, int*&, int*&, double*&, int&) const{}
#endif // Tasmanian_ENABLE_CUDA

void GridLocalPolynomial::buildTree(){
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();

    std::vector<int> level(num_points);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = work->getIndex(i);
        int current_level =rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            current_level += rule->getLevel(p[j]);
        }
        level[i] = current_level;
    }

    top_level = 0;
    for(auto l: level) if (top_level < l) top_level = l;

    int max_1d_kids = rule->getMaxNumKids();
    int max_kids = max_1d_kids*num_dimensions;
    std::vector<int> monkey_count(top_level + 1);
    std::vector<int> monkey_tail(top_level + 1);

    std::vector<int>  tree(max_kids * num_points, -1);
    std::vector<bool> free(num_points, true);
    std::vector<int>  kid(num_dimensions);

    int next_root = 0;
    roots.resize(0);

    while(next_root != -1){
        roots.push_back(next_root);
        free[next_root] = false;

        monkey_tail[0] = next_root;
        monkey_count[0] = 0;
        int current = 0;

        while(monkey_count[0] < max_kids){
            if (monkey_count[current] < max_kids){
                const int *p = work->getIndex(monkey_tail[current]);

                int dir = monkey_count[current] / max_1d_kids;
                int ikid = rule->getKid(p[dir], monkey_count[current] % max_1d_kids);

                if (ikid == -1){
                    monkey_count[current]++;
                }else{
                    std::copy(p, p + num_dimensions, kid.data());
                    kid[dir] = ikid;
                    int t = work->getSlot(kid.data());
                    if ((t == -1) || (!free[t])){
                        monkey_count[current]++;
                    }else{
                        tree[monkey_tail[current] * max_kids + monkey_count[current]] =  t;
                        monkey_count[++current] = 0;
                        monkey_tail[current] = t;
                        free[t] = false;
                    }
                }
            }else{
                monkey_count[--current]++;
            }
        }

        next_root = -1;
        int next_level = top_level+1;
        for(int i=0; i<num_points; i++){
            if (free[i] && (level[i] < next_level)){
                next_root = i;
                next_level = level[i];
            }
        }
    }

    pntr.resize(num_points + 1);
    pntr[0] = 0;
    for(int i=0; i<num_points; i++){
        pntr[i+1] = pntr[i];
        for(int j=0; j<max_kids; j++) if (tree[i * max_kids + j] > -1) pntr[i+1]++;
    }

    indx.resize((pntr[num_points] > 0) ? pntr[num_points] : 1);
    indx[0] = 0;
    int count = 0;
    for(int i=0; i<num_points; i++){
        for(int j=0; j<max_kids; j++){
            int t = tree[i * max_kids + j];
            if (t > -1) indx[count++] = t;
        }
    }
}

void GridLocalPolynomial::getBasisIntegrals(double *integrals) const{
    IndexSet *work = (points == 0) ? needed : points;

    int n = 0;
    double *w = 0, *x = 0;

    if ((rule->getMaxOrder() == -1) || (rule->getMaxOrder() > 3) ){
        n = top_level / 2 + 1;
        OneDimensionalNodes::getGaussLegendre(n, w, x);
    }

    for(int i=0; i<work->getNumIndexes(); i++){
        const int* p = work->getIndex(i);
        integrals[i] = rule->getArea(p[0], n, w, x);
        for(int j=1; j<num_dimensions; j++){
            integrals[i] *= rule->getArea(p[j], n, w, x);
        }
    }

    if ((rule->getMaxOrder() == -1) || (rule->getMaxOrder() > 3)){
        delete[] x;
        delete[] w;
    }
}
void GridLocalPolynomial::getQuadratureWeights(double *weights) const{
    IndexSet *work = (points == 0) ? needed : points;
    getBasisIntegrals(weights);

    std::vector<int> monkey_count(top_level+1);
    std::vector<int> monkey_tail(top_level+1);

    double basis_value;

    const Data2D<int> *dagUp;
    Data2D<int> lparents;
    if (parents.getNumStrips() == work->getNumIndexes()){
        dagUp = &parents;
    }else{
        IndexManipulator IM(num_dimensions);
        IM.computeDAGupLocal(work, rule, lparents);
        dagUp = &lparents;
    }

    int num_points = work->getNumIndexes();
    std::vector<int> level(num_points);
    for(int i=0; i<num_points; i++){
        const int *p = work->getIndex(i);
        level[i] = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            level[i] += rule->getLevel(p[j]);
        }
    }

    std::vector<double> node(num_dimensions);
    int max_parents = rule->getMaxNumParents() * num_dimensions;

    for(int l=top_level; l>0; l--){
        for(int i=0; i<num_points; i++){
            if (level[i] == l){
                const int* p = work->getIndex(i);
                for(int j=0; j<num_dimensions; j++) node[j] = rule->getNode(p[j]);

                std::vector<bool> used(work->getNumIndexes(), false);

                monkey_count[0] = 0;
                monkey_tail[0] = i;
                int current = 0;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp->getCStrip(monkey_tail[current])[monkey_count[current]];
                        if ((branch == -1) || used[branch]){
                            monkey_count[current]++;
                        }else{
                            const int *func = work->getIndex(branch);
                            basis_value = rule->evalRaw(func[0], node[0]);
                            for(int j=1; j<num_dimensions; j++) basis_value *= rule->evalRaw(func[j], node[j]);
                            weights[branch] -= weights[i] * basis_value;
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

void GridLocalPolynomial::integrate(double q[], double *conformal_correction) const{
    int num_points = points->getNumIndexes();
    std::fill(q, q + num_outputs, 0.0);

    if (conformal_correction == 0){
        std::vector<double> integrals(num_points);
        getBasisIntegrals(integrals.data());
        for(int i=0; i<num_points; i++){
            const double *s = surpluses.getCStrip(i);
            double wi = integrals[i];
            for(int k=0; k<num_outputs; k++) q[k] += wi * s[k];
        }
    }else{
        std::vector<double> w(num_points);
        getQuadratureWeights(w.data());
        for(int i=0; i<num_points; i++){
            double wi = w[i] * conformal_correction[i];
            const double *vals = values->getValues(i);
            for(int k=0; k<num_outputs; k++) q[k] += wi * vals[k];
        }
    }
}

void GridLocalPolynomial::getNormalization(std::vector<double> &norms) const{
    norms.resize(num_outputs);
    std::fill(norms.begin(), norms.end(), 0.0);
    for(int i=0; i<points->getNumIndexes(); i++){
        const double *v = values->getValues(i);
        for(int j=0; j<num_outputs; j++){
            if (norms[j] < fabs(v[j])) norms[j] = fabs(v[j]);
        }
    }
}

void GridLocalPolynomial::buildUpdateMap(double tolerance, TypeRefinement criteria, int output, const double *scale_correction, Data2D<int> &map2) const{
    int num_points = points->getNumIndexes();
    map2.resize(num_dimensions, num_points);
    if (tolerance == 0.0){
        std::fill(map2.getVector()->begin(), map2.getVector()->end(), 1); // if tolerance is 0, refine everything
        return;
    }else{
        std::fill(map2.getVector()->begin(), map2.getVector()->end(), 0);
    }

    std::vector<double> norm;
    getNormalization(norm);

    int active_outputs = (output == -1) ? num_outputs : 1;
    Data2D<double> scale;
    if (scale_correction == 0){
        scale.resize(active_outputs, num_points, 1.0);
    }else{
        scale.cload(active_outputs, num_points, scale_correction);
    }

    if ((criteria == refine_classic) || (criteria == refine_parents_first)){
        #pragma omp parallel for
        for(int i=0; i<num_points; i++){
            bool small = true;
            const double *s = surpluses.getCStrip(i);
            const double *c = scale.getCStrip(i);
            if (output == -1){
                for(int k=0; k<num_outputs; k++) small = small && ((c[k] * fabs(s[k]) / norm[k]) <= tolerance);
            }else{
                small = ((c[0] * fabs(s[output]) / norm[output]) <= tolerance);
            }
            if (!small){
                int *m = map2.getStrip(i);
                std::fill(m, m + num_dimensions, 1);
            }
        }
    }else{
        // construct a series of 1D interpolants and use a refinement criteria that is a combination of the two hierarchical coefficients
        IndexManipulator IM(num_dimensions);
        Data2D<int> dagUp;
        IM.computeDAGupLocal(points, rule, dagUp);

        int max_1D_parents = rule->getMaxNumParents();

        SplitDirections split(points);

        #pragma omp parallel for
        for(int j=0; j<split.getNumJobs(); j++){ // split.getNumJobs() gives the number of 1D interpolants to construct
            int d = split.getJobDirection(j);
            int nump = split.getJobNumPoints(j);
            const int *pnts = split.getJobPoints(j);

            std::vector<int> global_to_pnts(num_points);
            std::vector<int> levels(nump);

            int max_level = 0;

            Data2D<double> vals;
            vals.resize(active_outputs, nump);

            for(int i=0; i<nump; i++){
                const double* v = values->getValues(pnts[i]);
                const int *p = points->getIndex(pnts[i]);
                if (output == -1){
                    std::copy(v, v + num_outputs, vals.getStrip(i));
                }else{
                    vals.getStrip(i)[0] = v[output];
                }
                global_to_pnts[pnts[i]] = i;
                levels[i] = rule->getLevel(p[d]);
                if (max_level < levels[i]) max_level = levels[i];
            }

            std::vector<int> monkey_count(max_level + 1);
            std::vector<int> monkey_tail(max_level + 1);

            for(int l=1; l<=max_level; l++){
                for(int i=0; i<nump; i++){
                    if (levels[i] == l){
                        const int *p = points->getIndex(pnts[i]);
                        double x = rule->getNode(p[d]);
                        double *valsi = vals.getStrip(i);

                        int current = 0;
                        monkey_count[0] = d * max_1D_parents;
                        monkey_tail[0] = pnts[i]; // uses the global indexes
                        std::vector<bool> used(nump, false);

                        while(monkey_count[0] < (d+1) * max_1D_parents){
                            if (monkey_count[current] < (d+1) * max_1D_parents){
                                int branch = dagUp.getStrip(monkey_tail[current])[monkey_count[current]];
                                if ((branch == -1) || (used[global_to_pnts[branch]])){
                                    monkey_count[current]++;
                                }else{
                                    const int *branch_point = points->getIndex(branch);
                                    double basis_value = rule->evalRaw(branch_point[d], x);
                                    const double *branch_vals = vals.getCStrip(global_to_pnts[branch]);
                                    for(int k=0; k<active_outputs; k++) valsi[k] -= basis_value * branch_vals[k];

                                    used[global_to_pnts[branch]] = true;
                                    monkey_count[++current] = d * max_1D_parents;
                                    monkey_tail[current] = branch;
                                }
                            }else{
                                monkey_count[--current]++;
                            }
                        }
                    }
                }
            }

            // at this point, vals contains the one directional surpluses
            for(int i=0; i<nump; i++){
                const double *s = surpluses.getCStrip(pnts[i]);
                const double *c = scale.getCStrip(pnts[i]);
                const double *v = vals.getCStrip(i);
                bool small = true;
                if (output == -1){
                    for(int k=0; k<num_outputs; k++){
                        small = small && (((c[k] * fabs(s[k]) / norm[k]) <= tolerance) || ((c[k] * fabs(v[k]) / norm[k]) <= tolerance));
                    }
                }else{
                    small = ((c[0] * fabs(s[output]) / norm[output]) <= tolerance) || ((c[0] * fabs(v[0]) / norm[output]) <= tolerance);
                }
                map2.getStrip(pnts[i])[d] = (small) ? 0 : 1;;
            }
        }
    }
}

bool GridLocalPolynomial::addParent(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude) const{
    std::vector<int> dad(num_dimensions);  std::copy(point, point + num_dimensions, dad.data());
    bool added = false;
    dad[direction] = rule->getParent(point[direction]);
    if ((dad[direction] != -1) && (exclude->getSlot(dad) == -1)){
        destination->addIndex(dad);
        added = true;
    }
    dad[direction] = rule->getStepParent(point[direction]);
    if ((dad[direction] != -1) && (exclude->getSlot(dad) == -1)){
        destination->addIndex(dad);
        added = true;
    }
    return added;
}
void GridLocalPolynomial::addChild(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude)const{
    std::vector<int> kid(num_dimensions);  std::copy(point, point + num_dimensions, kid.data());
    int max_1d_kids = rule->getMaxNumKids();
    for(int i=0; i<max_1d_kids; i++){
        kid[direction] = rule->getKid(point[direction], i);
        if ((kid[direction] != -1) && (exclude->getSlot(kid) == -1)){
            destination->addIndex(kid);
        }
    }
}
void GridLocalPolynomial::addChildLimited(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude, const std::vector<int> &level_limits) const{
    std::vector<int> kid(num_dimensions);  std::copy(point, point + num_dimensions, kid.data());
    int max_1d_kids = rule->getMaxNumKids();
    for(int i=0; i<max_1d_kids; i++){
        kid[direction] = rule->getKid(point[direction], i);
        if ((kid[direction] != -1)
            && ((level_limits[direction] == -1) || (rule->getLevel(kid[direction]) <= level_limits[direction]))
            && (exclude->getSlot(kid) == -1)){
            destination->addIndex(kid);
        }
    }
}

void GridLocalPolynomial::clearRefinement(){
    if (needed != 0){ delete needed; needed = 0; }
}
const double* GridLocalPolynomial::getSurpluses() const{
    return surpluses.getVector()->data();
}
const int* GridLocalPolynomial::getPointIndexes() const{
    return ((points == 0) ? needed->getIndex(0) : points->getIndex(0));
}
const int* GridLocalPolynomial::getNeededIndexes() const{
    return ((needed != 0) ? needed->getIndex(0) : 0);
}
void GridLocalPolynomial::setSurplusRefinement(double tolerance, TypeRefinement criteria, int output, const std::vector<int> &level_limits, const double *scale_correction){
    clearRefinement();

    Data2D<int> pmap;
    buildUpdateMap(tolerance, criteria, output, scale_correction, pmap);

    bool useParents = (criteria == refine_fds) || (criteria == refine_parents_first);

    GranulatedIndexSet *refined = new GranulatedIndexSet(num_dimensions);

    int num_points = points->getNumIndexes();

    if (level_limits.empty()){
        for(int i=0; i<num_points; i++){
            const int *map = pmap.getCStrip(i);
            for(int j=0; j<num_dimensions; j++){
                if (map[j] == 1){ // if this dimension needs to be refined
                    if (!(useParents && addParent(points->getIndex(i), j, refined, points))){
                        addChild(points->getIndex(i), j, refined, points);
                    }
                }
            }
        }
    }else{
        for(int i=0; i<num_points; i++){
            const int *map = pmap.getCStrip(i);
            for(int j=0; j<num_dimensions; j++){
                if (map[j] == 1){ // if this dimension needs to be refined
                    if (!(useParents && addParent(points->getIndex(i), j, refined, points))){
                        addChildLimited(points->getIndex(i), j, refined, points, level_limits);
                    }
                }
            }
        }
    }

    //////////// TEST CODE /////////////////////
    // complete the set of points to ensure no missing parents
    // often times it results in real excess of points,
    // but I had to test
//    if (refined->getNumIndexes() > 0){
//        IndexManipulator IM(num_dimensions);
//        IndexSet* total = new IndexSet(refined);
//        IndexSet* completion = IM.getLowerCompletion(total);
//        cout << "total = " << total->getNumIndexes() << "  " << completion->getNumIndexes() << endl;
//        if ((completion != 0) && (completion->getNumIndexes() > 0)){
//            total->addIndexSet(completion);
//            delete completion;
//        }
//        needed = total->diffSets(points);
//        cout << " needed = " << needed->getNumIndexes() << endl;
//    }
    //////// END TEST CODE /////////////////
    if (refined->getNumIndexes() > 0){
        needed = new IndexSet(refined);
    }
    delete refined;
}
int GridLocalPolynomial::removePointsByHierarchicalCoefficient(double tolerance, int output, const double *scale_correction){
    clearRefinement();
    int num_points = points->getNumIndexes();
    std::vector<bool> pmap(num_points); // point map, set to true if the point is to be kept, false otherwise

    std::vector<double> norm;
    getNormalization(norm);

    Data2D<double> scale;
    int active_outputs = (output == -1) ? num_outputs : 1;
    if (scale_correction == 0){
        scale.resize(active_outputs, num_points, 1.0);
    }else{
        scale.cload(active_outputs, num_points, scale_correction);
    }

    for(int i=0; i<num_points; i++){
        bool small = true;
        const double *s = surpluses.getCStrip(i);
        const double *c = scale.getCStrip(i);
        if (output == -1){
            for(int k=0; k<num_outputs; k++) small = small && ((c[k] * fabs(s[k]) / norm[k]) <= tolerance);
        }else{
            small = ((c[0] * fabs(s[output]) / norm[output]) <= tolerance);
        }
        pmap[i] = !small;
    }

    int num_kept = 0; for(int i=0; i<num_points; i++) num_kept += (pmap[i]) ? 1 : 0;

    if (num_kept == 0) return 0; // trivial case, remove all
    if (num_kept == num_points) return num_points; // trivial case, remove nothing

    // save a copy of the points and the values
    std::vector<int> point_kept(((size_t) num_kept) * ((size_t) num_dimensions));
    Data2D<int> pp;
    pp.load(num_dimensions, num_kept, point_kept.data()); // to handle double indexing

    StorageSet *values_kept = new StorageSet(num_outputs, num_kept);
    values_kept->aliasValues()->resize(((size_t) num_kept) * ((size_t) num_outputs));

    num_kept = 0;
    for(int i=0; i<num_points; i++){
        if (pmap[i]){
            std::copy(points->getIndex(i), points->getIndex(i) + num_dimensions, pp.getStrip(num_kept));
            std::copy(values->getValues(i), values->getValues(i) + num_outputs, values_kept->getValues(num_kept));
            num_kept++;
        }
    }

    int dims = num_dimensions, outs = num_outputs;

    reset(false);
    num_dimensions = dims; num_outputs = outs;

    points = new IndexSet(num_dimensions, point_kept);
    values = values_kept;

    buildTree();
    recomputeSurpluses();

    return points->getNumIndexes();
}

void GridLocalPolynomial::setHierarchicalCoefficients(const double c[], TypeAcceleration acc){
    #ifdef Tasmanian_ENABLE_CUDA
    clearCudaLoadedData();
    #endif
    if (accel != 0) accel->resetGPULoadedData();
    if (points != 0){
        clearRefinement();
    }else{
        points = needed;
        needed = 0;
    }
    surpluses.resize(num_outputs, getNumPoints());
    std::copy(c, c + surpluses.getTotalEntries(), surpluses.getVector()->data());

    std::vector<double> *vals = values->aliasValues();
    vals->resize(surpluses.getTotalEntries());

    std::vector<double> x(((size_t) getNumPoints()) * ((size_t) num_dimensions));
    getPoints(x.data());
    if (acc == accel_cpu_blas){
        evaluateBatchCPUblas(x.data(), points->getNumIndexes(), vals->data());
    }else if (acc == accel_gpu_cublas){
        evaluateBatchGPUcublas(x.data(), points->getNumIndexes(), vals->data());
    }else if (acc == accel_gpu_cuda){
        evaluateBatchGPUcuda(x.data(), points->getNumIndexes(), vals->data());
    }else{
        evaluateBatch(x.data(), points->getNumIndexes(), vals->data());
    }
}

#ifdef Tasmanian_ENABLE_CUDA
void GridLocalPolynomial::makeCheckAccelerationData(TypeAcceleration acc) const{
    if (AccelerationMeta::isAccTypeFullMemoryGPU(acc)){
        if ((accel != 0) && (!accel->isCompatible(acc))){
            delete accel;
            accel = 0;
        }
        if (accel == 0){ accel = (BaseAccelerationData*) (new AccelerationDataGPUFull()); }
    }
}
void GridLocalPolynomial::checkAccelerationGPUValues() const{
    AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
    double *gpu_values = gpu->getGPUValues();
    if (gpu_values == 0) gpu->loadGPUValues(points->getNumIndexes() * values->getNumOutputs(), surpluses.getCStrip(0));
}
void GridLocalPolynomial::checkAccelerationGPUNodes() const{
    AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
    if (gpu->getGPUNodes() == 0){
        IndexSet *work = (points == 0) ? needed : points;
        int num_entries = num_dimensions * work->getNumIndexes();
        double *cpu_nodes = new double[getNumPoints() * num_dimensions];
        getPoints(cpu_nodes);
        double *cpu_support = new double[num_entries];
        if (rule->getType() == rule_localp){
            switch(order){
            case 0: encodeSupportForGPU<0, rule_localp>(work, cpu_support); break;
            case 2: encodeSupportForGPU<2, rule_localp>(work, cpu_support); break;
            default:
                encodeSupportForGPU<1, rule_localp>(work, cpu_support);
            }
        }else if (rule->getType() == rule_semilocalp){
            encodeSupportForGPU<2, rule_semilocalp>(work, cpu_support);
        }else if (rule->getType() == rule_localpb){
            switch(order){
            case 2: encodeSupportForGPU<2, rule_localpb>(work, cpu_support); break;
            default:
                encodeSupportForGPU<1, rule_localpb>(work, cpu_support);
            }
        }else{
            switch(order){
            case 2: encodeSupportForGPU<2, rule_localp0>(work, cpu_support); break;
            default:
                encodeSupportForGPU<1, rule_localp0>(work, cpu_support);
            }
        }
        gpu->loadGPUNodesSupport(num_entries, cpu_nodes, cpu_support);
        delete[] cpu_nodes;
        delete[] cpu_support;
    }
}
void GridLocalPolynomial::checkAccelerationGPUHierarchy() const{
    AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
    int num_points = getNumPoints();
    gpu->loadGPUHierarchy(num_points, pntr.data(), indx.data(), (int) roots.size(), roots.data());
}
#else
void GridLocalPolynomial::makeCheckAccelerationData(TypeAcceleration) const{}
void GridLocalPolynomial::checkAccelerationGPUValues() const{}
void GridLocalPolynomial::checkAccelerationGPUNodes() const{}
void GridLocalPolynomial::checkAccelerationGPUHierarchy() const{}
#endif // Tasmanian_ENABLE_CUDA

void GridLocalPolynomial::clearAccelerationData(){
    #ifdef Tasmanian_ENABLE_CUDA
    cuda_engine.reset();
    clearCudaLoadedData();
    #endif
    if (accel != 0){
        delete accel;
        accel = 0;
    }
}

void GridLocalPolynomial::setFavorSparse(bool favor){
    // sparse_affinity == -1: use dense algorithms
    // sparse_affinity ==  1: use sparse algorithms
    // sparse_affinity ==  0: let Tasmanian decide
    // favor true/false bumps you upper or lower on the scale
    if (favor && (sparse_affinity < 1)) sparse_affinity++;
    if (!favor && (sparse_affinity > -1)) sparse_affinity--;
}

}

#endif
